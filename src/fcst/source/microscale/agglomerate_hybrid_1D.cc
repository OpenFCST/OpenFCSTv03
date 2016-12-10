//---------------------------------------------------------------------------
// C++ Interface: agglomerate_hybrid_1D.h
//
// Description: Used to solve a system of equations representing a
// 			spherical hybrid agglomerate.
//
// Author: Peter Dobson <pdobson@ualberta.ca>, (C) 2011
// 			University of Alberta
//
// Copyright: See COPYING file that comes with this distribution
//
//---------------------------------------------------------------------------

#include <microscale/agglomerate_hybrid_1D.h>
#ifdef _OPENMP
	#include <omp.h>
#else
	#define omp_get_thread_num() 0
#endif


namespace NAME = FuelCellShop::MicroScale;

const std::string NAME::HybridAgglomerate::concrete_name("HybridAgglomerateNumerical");
NAME::HybridAgglomerate const* NAME::HybridAgglomerate::PROTOTYPE = new NAME::HybridAgglomerate(concrete_name);


NAME::HybridAgglomerate::HybridAgglomerate(std::string concrete_name){
    //FcstUtilities::log <<" Register " + concrete_name +  " to FactoryMap"<<std::endl;
       this->get_mapFactory()->insert(
               std::pair<std::string, NAME::HybridAgglomerate*>(
                       concrete_name, this));
}


//---------------------------------------------------------------------------
NAME::HybridAgglomerate::HybridAgglomerate(int verbose)
:
WaterAgglomerate(verbose)
{
	//Same as water filled constructor
}

//---------------------------------------------------------------------------
void
NAME::HybridAgglomerate::set_structure()
{

    electrolyte = this->layer->get_resource<FuelCellShop::Material::PolymerElectrolyteBase>();
    kinetics = this->layer->get_resource<FuelCellShop::Kinetics::BaseKinetics>();


    //We need to recalculate the film thickness or porosity again
    //Since it was called before we had initia

    if (fixed_agg_variable.compare("Thickness") == 0)
    {

        epsilon_agg = compute_epsilon_agg();
        // n_agg depends on epsilon_agg, therefore the former needs to be computed first.
        n_agg = compute_n(this->layer);
    }
    else if (fixed_agg_variable.compare("Porosity") == 0)
    {
        n_agg = compute_n(this->layer);
        // thickenss_agg depends on n_agg, therefore it is computed last:
        delta_agg = compute_thickness_agg();
    }
    else
    {
        FcstUtilities::log<<"Undefined execution path in HybridAgglomerate::initialize "<< std::endl;
        abort();
    }



    r_agg *=1e-7;// convert to cm
    delta_agg *=1e-7; // convert to cm
    interface = r_agg/(r_agg+delta_agg);
    hybrid_interface = std::pow((1-hybrid_core_volume_fraction),1.0/3.0)*interface;

}



//---------------------------------------------------------------------------
void
NAME::HybridAgglomerate::fsub(double &x, double z[], double [], double f[])
{

	double D_O2_EFF = D_O2_N;
	double D_H_EFF = D_H_NAF; //9.2e-5; //SIGMA_EFF*R*T/(c_h[0]*pow(F,2.0));
	double rel_permittivity = rel_permittivity_Naf;//electrolyte->get_permittivity();
	
	std::vector<double> J(1,0.0);
	double c_fixed = c_H;

	
	if (x < interface)
	{

	    //if (z[0] > R_tol)
	    //{
        std::vector<SolutionVariable> c_reactants;
        c_reactants.push_back(SolutionVariable (z[0],1, oxygen_concentration)); //TODO: this is not generic
        c_reactants.push_back(SolutionVariable (z[2],1, proton_concentration)); //TODO: how efficient is this data copying?
        SolutionVariable v_membrane(z[4], 1, protonic_electrical_potential);
        SolutionVariable v_solid(phi_S, 1, electronic_electrical_potential);

        this->kinetics->set_reactant_concentrations(c_reactants);
        this->kinetics->set_electrolyte_potential(v_membrane);
        this->kinetics->set_solid_potential(v_solid);

        this->kinetics->current_density(J);

        J[0] = getAV(x)*J[0];
	    //}
		D_O2_EFF = pow(epsilon_agg,1.5) * DO2_Water;
		D_H_EFF = D_H_Water;

		rel_permittivity = rel_permittivity_water; //water.get_permittivity();
		c_fixed = 0.0;
	}


	double T = this->solutions[temperature_of_REV][this->sol_index];

	f[0] = z[1]/D_O2_EFF;
	f[1] = (((pow((r_agg+delta_agg),2.0))*J[0])/(4.0*F))-((2.0*z[1])/x);

	if ((x < interface) and (x > hybrid_interface)){
	    D_H_EFF = D_H_NAF;
	    rel_permittivity = rel_permittivity_Naf;
	    c_fixed = c_H;
	    D_O2_EFF = pow(epsilon_agg,1.5) *D_O2_N;

	}

	f[2] = z[3]/D_H_EFF - (F/(R*T))*z[2]*z[5];

	f[3] = (pow((r_agg+delta_agg),2.0)*J[0]/(1.0*F))-((2.0*z[3])/x);
	f[4] = -((F/(permittivity_0*rel_permittivity))*(z[2]-c_fixed)*pow((r_agg+delta_agg),2.0)) - ((2.0*z[5])/x);

	//6 equations
	//f[4] = z[5]/(permittivity_0*rel_permittivity);
	//f[5] = (-1.0*F*(z[2]-c_fixed)*pow((r_agg+delta_agg),2.0)) - ((2.0*z[5])/x);
}

//---------------------------------------------------------------------------
void
NAME::HybridAgglomerate::dfsub(double &x, double z[], double [], double df[])
{

    double D_O2_EFF = D_O2_N;
    double D_H_EFF = D_H_NAF; //SIGMA_EFF*R*T/(c_h[0]*pow(F,2.0));
    double rel_permittivity = rel_permittivity_Naf;//electrolyte->get_permittivity();

    int n_eq = n_comp + n_y;
    int n_var = m_star + n_y;
    //For the two-dimensional array (or matrix), declare an array of pointers.
    double** dfc;
    dfc = new double*[n_eq];
    //For each ptr in the array, assign it an array of 4 doubles.
    for (int i=0; i < n_eq; ++i)
        dfc[i] = new double[n_var];
    // Initialize the matrix with zeros.
    for (int i=0; i < n_eq; i++)
        for (int j=0; j <n_var; j++)
            dfc[i][j]= 0.0;


	double c_fixed = c_H;
	double PARTIALJ0(0.0);
	double PARTIALJ1(0.0);
	double PARTIALJ2(0.0);
	if (x < interface)
	{
		D_O2_EFF = pow(epsilon_agg,1.5) * DO2_Water;
		D_H_EFF = pow(epsilon_agg,1.5) * D_H_Water;
		rel_permittivity = rel_permittivity_water;//water.get_permittivity();
		c_fixed = 0.0;

		std::vector<SolutionVariable> c_reactants;
		c_reactants.push_back(SolutionVariable (z[0], 1,oxygen_concentration)); //TODO: make this generic
		c_reactants.push_back(SolutionVariable (z[2], 1,proton_concentration)); //TODO: make this generic
		SolutionVariable v_membrane(z[4],1, protonic_electrical_potential);
		SolutionVariable v_solid(phi_S, 1, electronic_electrical_potential);

		this->kinetics->set_reactant_concentrations(c_reactants);
		this->kinetics->set_electrolyte_potential(v_membrane);
		this->kinetics->set_solid_potential(v_solid);
		std::map< VariableNames, std::vector<double> > derivatives;

		this->kinetics->derivative_current(derivatives);

		PARTIALJ0 = getAV(x) * (derivatives[oxygen_molar_fraction][0] *(HO2N/P) );
		PARTIALJ1 = getAV(x) * derivatives[proton_concentration][0];
		PARTIALJ2 = getAV(x) * derivatives[protonic_electrical_potential][0];

	}


	if ((x < interface) and (x > hybrid_interface)){
	        D_H_EFF = D_H_NAF;
	        rel_permittivity = rel_permittivity_Naf;
	        c_fixed = c_H;
	        D_O2_EFF = pow(epsilon_agg,1.5) *D_O2_N;
	    }

	double T = this->solutions[temperature_of_REV][this->sol_index];
    dfc[0][1] = 1.0/D_O2_EFF;
    // All others = 0

    dfc[1][0] = ((pow((r_agg+delta_agg),2.))/(4.0*F)) * PARTIALJ0;
    dfc[1][1] = -((2.0)/(x));
    dfc[1][2] = ((pow((r_agg+delta_agg),2.))/(4.0*F)) * PARTIALJ1;
    dfc[1][4] = ((pow((r_agg+delta_agg),2.))/(4.0*F)) * PARTIALJ2;
    // All others = 0

    dfc[2][2] = -(F/(R*T))*z[5];
    dfc[2][3] = 1.0/D_H_EFF;
    dfc[2][5] = -(F/(R*T))*z[2];
    // All others = 0

    dfc[3][0] = (pow((r_agg+delta_agg),2.)/(1.0*F)) * PARTIALJ0;
    dfc[3][2] = (pow((r_agg+delta_agg),2.)/(1.0*F))  * PARTIALJ1;
    dfc[3][3] = -((2.0)/(x));
    dfc[3][4] = (pow((r_agg+delta_agg),2.)/(1.0*F))  * PARTIALJ2;
    // All others = 0

    // 5 equations
    dfc[4][2] = -(F/(permittivity_0*rel_permittivity)) * pow((r_agg+delta_agg),2.0);
    dfc[4][5] = -((2.0)/x);
    // All others = 0

    // 6 equations
    //dfc[4][5] = 1.0/(permittivity_0*rel_permittivity);
    //dfc[5][2] = -F * pow((r_agg+delta_agg),2.0);
    //dfc[5][5] = -2.0/x;

    // All others = 0

    //Convert the matrix dfc to a one dimensional array that Fortran will understand.
    FuelCell::ApplicationCore::c_to_for_matrix(n_eq,n_var,dfc,df);

    // Free the memory
    for (int i=0; i < n_eq; ++i)
        delete [] dfc[i];
    delete [] dfc;
}



//---------------------------------------------------------------------------
void
NAME::HybridAgglomerate::fsub_wrapper (double &x, double z[], double y[], double f[])
{
	HybridAgglomerate *ptr_Agglomerate = ((HybridAgglomerate*) (ptr_DAE_object[omp_get_thread_num()]));
	ptr_Agglomerate->fsub(x, z, y, f);
}

//---------------------------------------------------------------------------
void
NAME::HybridAgglomerate::dfsub_wrapper (double &x, double z[], double y[], double f[])
{
	HybridAgglomerate *ptr_Agglomerate = ((HybridAgglomerate*) (ptr_DAE_object[omp_get_thread_num()]));
	ptr_Agglomerate->dfsub(x, z, y, f);
}

//---------------------------------------------------------------------------



double NAME::HybridAgglomerate::compute_thickness_agg()
{
    //set default and starting values
       double thickness_0 = r_agg / 10.;
       const double tol = 1e-6;
       const double eps_tol = 1e-6;
       const double iterations = 1000;

       CL_Properties props = this->layer->get_properties();

       // Agglomerate is partially ionomer filed within the agglomerate.
       double porosity =  epsilon_agg * hybrid_core_volume_fraction;

       double thickness_i = thickness_0;

       for (int i = 0; i < iterations; ++i)
       {
           thickness_i = thickness_0
                   - ((compute_epsilon_N(thickness_0, porosity) - props[CLPropNames::ionomer_fraction])
                           / compute_depsilonN_dthickness(thickness_0));

           if (fabs((thickness_i - thickness_0) / thickness_0) < tol
                   || fabs(compute_epsilon_N(thickness_i, porosity) - compute_epsilon_N(thickness_0, porosity))
                           < eps_tol)
           {
               if (thickness_i < 0)
               {
                   FcstUtilities::log<<"Warning:: Thickness is negative -- Results to not make sense"<<std::endl;
                   return -1;
               }
               else
               {
                   //FcstUtilities::log<<"True eN:"<<compute_epsilon_N(thickness_i)<<std::endl;
                   //FcstUtilities::log<<"eN:"<<this->epsilon_N<<std::endl;
                   return thickness_i;
               }
           }
           else
               thickness_0 = thickness_i;

       }
       //if value does not converge, return value that will drive epsilon_N_cat above 1
       //error handling will set porosity to a minimum
       //BEWARE: program will run as normal - however, FcstUtilities::log output will output this negative value
       return -1;
}


double NAME::HybridAgglomerate::compute_epsilon_agg()
{
    AssertThrow(false, ExcMessage("Hybrid agglomerate cannot have a constant Thickness"));

}


void
NAME::HybridAgglomerate::print_properties(){
    FcstUtilities::log << "=========== CL MICROSTRUCTURE =========" << std::endl;
    FcstUtilities::log <<  "Agglomerate Type: " << get_name() << std::endl;
    FcstUtilities::log <<  "Agglomerate Outer Core Radius: " << std::setw(6) <<   get_radius() << " [nm]" << std::endl;
    FcstUtilities::log <<  "Agglomerate Inner Core Radius: " << std::setw(6) <<   get_radius()*(hybrid_interface/interface) << " [nm]" << std::endl;
    FcstUtilities::log <<  "Agglomerate Porosity: " << std::setw(6) << epsilon_agg << std::endl;
    FcstUtilities::log <<  "Agglomerate Thin Film: " << std::setw(6) << get_film_thickness() << " [nm]" << std::endl;
    FcstUtilities::log << "=======================================" << std::endl;
}


