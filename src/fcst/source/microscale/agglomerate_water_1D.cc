//---------------------------------------------------------------------------
// C++ Interface: agglomerate_water_1D.h
//
// Description: Used to solve a system of equations representing a
// 			spherical water-filled agglomerate.
//
// Author: Peter Dobson <pdobson@ualberta.ca>, (C) 2011
// 			University of Alberta
//
// Copyright: See COPYING file that comes with this distribution
//
//---------------------------------------------------------------------------

#include <microscale/agglomerate_water_1D.h>
#ifdef _OPENMP
	#include <omp.h>
#else
	#define omp_get_thread_num() 0
#endif


namespace NAME = FuelCellShop::MicroScale;

const std::string NAME::WaterAgglomerate::concrete_name("WaterAgglomerateNumerical");
NAME::WaterAgglomerate const* NAME::WaterAgglomerate::PROTOTYPE = new NAME::WaterAgglomerate(concrete_name);


NAME::WaterAgglomerate::WaterAgglomerate(std::string concrete_name){
    //FcstUtilities::log <<" Register " + concrete_name +  " to FactoryMap"<<std::endl;
       this->get_mapFactory()->insert(
               std::pair<std::string, NAME::WaterAgglomerate*>(
                       concrete_name, this));
}


//---------------------------------------------------------------------------
NAME::WaterAgglomerate::WaterAgglomerate(int verbose)
:

H_tol(1e-12),
lambda(1.0),
lambda2(1.0)
{
    R_tol =1e-12;
	verbosity(verbose);
	this->has_derivatives_ = false;
	sol_names.resize(4);

	column_names.push_back("x"); column_names.push_back("z0");
    column_names.push_back("z1");  column_names.push_back("z2");
    column_names.push_back("z3"); column_names.push_back("z4");
    column_names.push_back("z5"); column_names.push_back("i");

    //concrete_name = "numericalwater";

}

//---------------------------------------------------------------------------
void
NAME::WaterAgglomerate::set_structure()
{

    electrolyte = this->layer->get_resource<FuelCellShop::Material::PolymerElectrolyteBase>();
    kinetics = this->layer->get_resource<FuelCellShop::Kinetics::BaseKinetics>();

    r_agg *=1e-7;// convert to cm
    delta_agg *=1e-7; // convert to cm
    interface = r_agg/(r_agg+delta_agg);

}

//---------------------------------------------------------------------------
FuelCellShop::SolutionMap
NAME::WaterAgglomerate::compute_current ()
{

    //In water filled case will have changed because the
    //layer's porosity is changed by the water filled assumption.
    //At set structure time, the layer porosity was not
    //yet corrected.

    CL_Properties props = this->layer->get_properties();
    AV = props[CLPropNames::active_area_scaled];
    AV /= pow(interface,3.0); //new scaling

    //True pressure is now available from the layer.
    //Was 0 at set structure time.
    P= props[CLPropNames::pressure];


    double E_r = 0.0;
    if(this->solutions[this->reactant][this->sol_index] <= 0.0){
        SolutionMap sols;
        sols.push_back(SolutionVariable(0, 1, current_density));
        sols.push_back(SolutionVariable(0,1, CL_effectiveness));

        return sols; //To stabalize FEM solution
    }


    if(this->solutions[this->reactant].get_variablename() == hydrogen_molar_fraction)
        AssertThrow(false, ExcMessage("Hydrogen reaction not yet implemented in numerical ionomer agglomerate."));

    // Initialize variables
    //electrolyte->proton_conductivity(sigma_p);
    D_H_Water = water.get_DH(); // FOR VALIDATION: 9.2e-5
    DO2_Water = water.get_DO2();  // FOR VALIDATION: 9.19e-5

    //Necessary call to base initial solution functions
    update_initial_solution();


    electrolyte->set_T(this->solutions[temperature_of_REV][this->sol_index]);
    electrolyte->set_lambda(this->solutions[membrane_water_content][this->sol_index]);


    electrolyte->proton_diffusivity(D_H_NAF);
    electrolyte->oxygen_diffusivity(D_O2_N);

    rel_permittivity_water = water.get_Relative_Permittivity(); // FOR VALIDATION: 80
    rel_permittivity_Naf = electrolyte->get_permittivity();
    HO2N = electrolyte->get_H_O2();


    sol_names.push_back(proton_concentration); //Add an extra solution name for the water filled model
    kinetics->set_derivative_flags(sol_names);

    // Set the oxygen concentration and membrane potential at the boundary

    c_R = (this->solutions[this->reactant][this->sol_index]*P)/HO2N;// FOR VALIDATION: c_O2 = 1.2e-6;
    c_H = electrolyte->get_density() / electrolyte->get_EW(); // FOR VALIDATION: c_H = 1.2e-4;

    // FOR VALIDATION: Over potential should be 0.4, change phi_S and phi_M accordingly from application
    // Set the solid phase potential for the domain
    phi_S = this->solutions[electronic_electrical_potential][this->sol_index];
    phi_M = this->solutions[protonic_electrical_potential][this->sol_index];



    SolutionVariable temperature(this->solutions[temperature_of_REV][this->sol_index],2, VariableNames::temperature_of_REV); //TODO change this 2 hard coded variable
    kinetics->set_temperature(temperature);
    kinetics->set_p_t(P);
    /*
    Solve the problem.  An indicator of success or failure
    is returned.  See IFLAG in COLDAE.f for meaning of
    the return value.
    Once the problem is solved, a
    continuous numerical solution is determined.
    */
    bool obtained_solution = false;
    int flag = 0;

    //Setup the solver
    lambda = 1.0; // try full solution first
    setup_DAE_solver();

    //Solve the problem
    flag = prob->DAE_solve();
    if (flag <= 0)
    {
        //Solver failed
        obtained_solution = false;
    }
    else obtained_solution = true;

    if (!obtained_solution)
    {
        #ifdef DEBUG
            FcstUtilities::log << get_name () + "attempting Continuation." << std::endl;
        #endif

        clear_memory();
        // a simple continuation strategy based on tolerance
        if (cont_tolerance(1e-3,1e-9)) //cont_tolerance(0.1,1e-4,1e-9)
        {
            obtained_solution = true;
        }
        // use simple continuation on curent density.
        else if (cont_cd())//cont_cd()
        {
            obtained_solution = true;
        }
        else
        {
            //all continuation strategies failed
            obtained_solution = false;
        }
    }


    double unit_size = boundary_1 - boundary_0;
    double agg_size = interface - boundary_0;
    double volume = (4.0/3.0)*pi*pow(unit_size,3.0);
    double agg_volume = (4.0/3.0)*pi*pow(agg_size,3.0);
    double I = 0.0; // Integral of current over the mesh
    double I_max = 0.0;
    double C_OH = 0.0;
    double C_O = 0.0;


    std::vector<double> J(1,0.0);

    //Get the maximum theoretical current density
    std::vector<SolutionVariable> c_reactants;
    c_reactants.push_back(SolutionVariable (c_R,1,oxygen_concentration)); //TODO: make this generic
    c_reactants.push_back(SolutionVariable (c_H,1,proton_concentration));
    SolutionVariable v_membrane(phi_M, 1, protonic_electrical_potential);
    SolutionVariable v_solid(phi_S, 1, electronic_electrical_potential);
    //SolutionVariable temperature(this->solutions[temperature_of_REV][this->sol_index],2, VariableNames::temperature_of_REV); //TODO change this 2 hard coded variable

    this->kinetics->set_reactant_concentrations(c_reactants);
    this->kinetics->set_electrolyte_potential(v_membrane);
    this->kinetics->set_solid_potential(v_solid);
    this->kinetics->current_density(J);

    I_max = (AV) * J[0] * agg_volume/volume;


    if (obtained_solution)
    {
        int n_mesh_final = prob->get_size_final_mesh();
        mesh = new double[n_mesh_final];
        prob->get_copy_final_mesh(mesh);
        //output_solution_to_file("agglomerate_sol.txt");

        // loop over each subinterval, defined by the mesh points
        // i.e. Composite integration.
        // NOTE: Must be implemented in each class... values are very problem specific.
        final_results.clear();
        std::vector<std::vector<double>> new_lineS;
        std::vector<double> new_line;

        //Temp Variables
        std::vector<double> c_oxygen;
        std::vector<double> c_proton;
        std::vector<double> V_membrane;
        std::vector<double> V_solid;


        //For each cell
        for (int j=0;j<n_mesh_final;++j)
        {
            double left, right;
            left = mesh[j];
            right = mesh[j+1];

            // Get the quadrature points based on the order of the solution on the subinterval
            std::vector<double> x_quad, w_quad;
            get_quadrature_points(left,right,x_quad,w_quad,prob);

            //get the solution at the quadrature points
            J.clear(); c_oxygen.clear(); c_proton.clear(); V_membrane.clear(); V_solid.clear();
            J.resize(x_quad.size());

            c_oxygen.resize(x_quad.size());
            c_proton.resize(x_quad.size());
            V_membrane.resize(x_quad.size());
            V_solid.resize(x_quad.size());

            new_lineS.clear();
            //Get data for each quad point per cell
            for (unsigned int i=0; i<x_quad.size(); ++i)
            {
                new_line.clear();
                double *z = new double[m_star];
                double *y = new double[n_y];
                prob->DAE_solution(x_quad[i],z,y); //Solution is returned in z.
                c_oxygen[i] = z[0];
                c_proton[i] = z[2];
                V_membrane[i] = z[4];//phi_M;
                V_solid[i] = phi_S;


                new_line.push_back(x_quad[i]); new_line.push_back(z[0]); new_line.push_back(z[1]);
                new_line.push_back(z[2]); new_line.push_back(z[3]); new_line.push_back(z[4]);
                new_line.push_back(z[5]);
                new_lineS.push_back(new_line);

                delete [] z;
                delete [] y;
            }

            //Compute the current
            c_reactants[0] = SolutionVariable(c_oxygen, oxygen_concentration); //TODO: make this generic
            c_reactants[1] = SolutionVariable(c_proton, proton_concentration);
            v_membrane = SolutionVariable(V_membrane, protonic_electrical_potential);
            v_solid = SolutionVariable(V_solid, electronic_electrical_potential);

            this->kinetics->set_reactant_concentrations(c_reactants);
            this->kinetics->set_electrolyte_potential(v_membrane);
            this->kinetics->set_solid_potential(v_solid);
            this->kinetics->set_temperature(temperature);
            this->kinetics->current_density(J);
            std::vector<double> c_OH(2, 0.0), c_O(2, 0.0);
            if(this->kinetics->has_coverage(OH_coverage))
                this->kinetics->OH_coverage(c_OH);
            if(this->kinetics->has_coverage(O_coverage))
                this->kinetics->O_coverage(c_O);

            // scale by Active area and CL volume, setup for integral over a sphere
            for (unsigned int i=0; i<x_quad.size(); ++i)
            {
                if ((c_oxygen[i]<R_tol) or (c_proton[i]<H_tol) or ((right > interface) or ((j+1) == (n_mesh_final))) )
                    J[i] = 0.0;
                else
                    J[i] = getAV(x_quad[i]) * J[i] * (4.0*pi*pow(x_quad[i],2.0));
            }

            for(unsigned int z =0; z < new_lineS.size(); z++){
                new_lineS.at(z).push_back(J[z]);
                final_results.push_back(new_lineS.at(z));
            }

            // Add the contribution of the spherical shell
            I += integrate(left,right,w_quad,J);
            C_O += integrate(left,right,w_quad,c_O);
            C_OH += integrate(left,right,w_quad,c_OH);

        }

        //Necessary call to base save solution functions
        save_initial_solution();


        delete [] mesh;
        clear_memory();

        double I_avg = I/volume;
        E_r = I_avg/I_max;
        SolutionMap sols;
        sols.push_back(SolutionVariable(I_avg, 1, current_density));
        sols.push_back(SolutionVariable(E_r,1, CL_effectiveness));
        sols.push_back(SolutionVariable(C_OH/volume,1, OH_coverage));
        sols.push_back(SolutionVariable(C_O/volume,1, O_coverage));
        return sols;

    }
    else
    {
        SolutionMap sols;
        sols.push_back(SolutionVariable(I_max, 1, current_density));
        sols.push_back(SolutionVariable(E_r,1, CL_effectiveness));

        DAE_Error(flag);
        return sols;
    }

}

//---------------------------------------------------------------------------
int 
NAME::WaterAgglomerate::cont_tolerance (double start_tol, double end_tol)
{
	//FcstUtilities::log << "Atempting continuation on tolerance " << std::endl;
	
	//Attempt to solve problem using continuation on tolerances
	int flag=0;
	double new_tol=0;
	int solve_attempts=0;
	bool obtained_sol = false;
	int save_attempts=0;
	//Set new meshs sizes and reset work arrays.
	do 
	{
		setup_DAE_solver();
		new_tol = end_tol*pow(1.20,++solve_attempts);
		//FcstUtilities::log << "Start tol: " << new_tol << std::endl;
		for (int i=0; i<n_comp;i++) tol[i] = new_tol;
		tol[3] = 1e-6;
		if(new_tol >=start_tol)
		{
			//failed to find a good starting point
			return 0;
		}
		flag=prob->DAE_solve();
		if (flag == 1) obtained_sol = true;
		else clear_memory();
	}while(!obtained_sol);
	
	clear_memory();
	setup_DAE_solver();
	// 	prob->set_integer_space(1000000);
	// 	prob->set_float_space(10000000);
	//Set initial tolerances
	for (int i=0; i<n_comp;i++) tol[i] = new_tol;
	tol[3] = 1e-6;
	solve_attempts=0;
	do 
	{
		//Apply continuation	
		flag = prob->DAE_solve();
		if (flag <= 0)
		{
			//this continuation strategy has failed
			//FcstUtilities::log << "CONT ENDED TOO SOON AT: " << tol[1] << std::endl;
			new_tol = new_tol*pow(1.05,++save_attempts);	
			if (save_attempts > 20)
			{
				clear_memory();
				//FcstUtilities::log << "tolerance continuation strategu failed." << std::endl;
				return 0;
			}
			prob->use_simple_cont();
		} 
		else
		{
			save_attempts=0;
			new_tol = new_tol/pow(1.10,++solve_attempts); //pow(1.05,++solve_attempts);
			if(solve_attempts > 40)
			{
				clear_memory();
				//FcstUtilities::log << "tolerance continuation strategy failed." << std::endl;
				return 0;
			}
			prob->use_simple_cont();
			
		}
		
		for (int i=0; i<n_comp;i++) tol[i] = new_tol;
		tol[3] = 1e-6;
		//FcstUtilities::log << "New tol: " << tol[1] << std::endl;
		
		
	}while(new_tol >= end_tol);
	
	return 1;
}

//---------------------------------------------------------------------------
int
NAME::WaterAgglomerate::cont_cd()
{
	double E0 = 1.1; // Approximately - add set OCV later
	double eta = E0 - (phi_S - phi_M);
	double eta_max = 0.6;
	
	bool obtained_solution = false;
	int flag = 0;
	int solve_attempts = 0;
	
	do 
	{
		lambda = 1.0/pow(10.,solve_attempts++); // try full solution first
		lambda2 = 1.0;
		if (eta > eta_max) lambda2 = lambda;
		//		FcstUtilities::log << "Parameter Continuation - lambda = " << lambda << " - " << lambda2;
		setup_DAE_solver();
		if(solve_attempts >= 10)
		{
			//failed to find a good starting point
			lambda = 1.0;
			return 0;
		}
		flag = prob->DAE_solve();
		if (flag == 1) obtained_solution = true;
		else clear_memory();
	}while(!obtained_solution );

	
	while (lambda < 1.0)
	{
		lambda = pow(lambda,0.7);
		if (lambda > 0.99) 
			lambda = 1.0;
		if (eta > eta_max) lambda2 = lambda;
		
		//				FcstUtilities::log << "\tParameter Continuation - lambda =  " << lambda;
		prob->use_simple_cont();
		flag = prob->DAE_solve();
		if(flag<=0)
		{
			// strategy failed
			lambda = 1.0;
			lambda2 = 1.0;
			clear_memory();
			return 0;
		}
	}

	lambda = 1.0;
	lambda2 = 1.0;
	return 1;

}

//---------------------------------------------------------------------------
void
NAME::WaterAgglomerate::setup_DAE_solver ()
{
//Declare an array that holds the number and order of each of the ODEs
n_comp = 5;
n_y = 0;
mm = new int[n_comp];
mm[0] = 1;
mm[1] = 1;
mm[2] = 1;
mm[3] = 1;
mm[4] = 2;
//mm[5] = 1;
m_star = 0;
for (int i=0;i<n_comp;++i)
	m_star += mm[i];

// Set the mesh parameters
n_mesh = 100;
n_colloc = 3;
boundary_0 = 0.0;
boundary_1 = 1.0;

bool set_fixpoint = false;
if (interface < 1.0) set_fixpoint = true;
fixpnt = new double[1]; // Use a mesh with a fixpoint 
fixpnt[0] = interface;

//declare the location of each boundary point
zeta = new double[m_star];
zeta[0] = boundary_0;
zeta[1] = boundary_0;
zeta[2] = boundary_0;
zeta[3] = boundary_1;
zeta[4] = boundary_1;
zeta[5] = boundary_1;

ltol = new int[n_comp + n_y];
tol = new double[n_comp + n_y];
for (int i=0; i<n_comp;++i)
{
	ltol[i] = i+1;
	tol[i] = 1e-9;
}
tol[4] = 1e-6;

// WARNING!!  This must be done in other implementations using the DAEWrapper class
//	Failing to do so will cause runtime errors that are difficult to locate.
// Assign the global object pointer to the WaterAgglomerate object!!!
ptr_DAE_object[omp_get_thread_num()] = (void*) this;

//Create an instance of the DAE solver class.  
prob = new FuelCell::ApplicationCore::DAESolver
(n_comp,  //Number of ODES
n_y, // Number of algebraic constraints (0 for a BVP)
mm, // Array that holds the order of each of the ODEs
boundary_0, // Leftmost boundary point
boundary_1, // Rightmost boundary point
zeta, // Location of boundary points
&fsub_wrapper,  // ptr to ODE function
&dfsub_wrapper, //ptr to Jacobian of ODE function
&gsub_wrapper, //ptr to boundary-condition function
&dgsub_wrapper, //ptr to derivatives of boundary-condition function
&guess_wrapper //ptr to optional initial guess function
);

/* 
These are methods used to set certain parameters of COLDAE.
They must be called before the problem is actually solved. 
These methods are optional and are not required to be used 
for every problem. 
*/

prob->set_collocation_points(n_colloc); //Set the number of collocation points 
prob->set_initial_mesh_size(n_mesh); //Set the initial mesh size
if (set_fixpoint) prob->set_fixpnts(1, fixpnt);
prob->set_tolerance(n_comp,ltol,tol);
prob->set_output(n_output); //Set the output level 
				//-1 is full output
				//0 is selected output
				//1 is no output
}

//---------------------------------------------------------------------------
void
NAME::WaterAgglomerate::fsub(double &x, double z[], double [], double f[])
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
        c_reactants.push_back(SolutionVariable (z[0],1, oxygen_concentration)); //TODO: make this generic
        c_reactants.push_back(SolutionVariable (z[2],1, proton_concentration));
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
		c_fixed = core_charge_factor*c_H;
	}


	double T = this->solutions[temperature_of_REV][this->sol_index];

	f[0] = z[1]/D_O2_EFF;
	f[1] = (((pow((r_agg+delta_agg),2.0))*J[0])/(4.0*F))-((2.0*z[1])/x);

	f[2] = z[3]/D_H_EFF - (F/(R*T))*z[2]*z[5];
	f[3] = (pow((r_agg+delta_agg),2.0)*J[0]/(1.0*F))-((2.0*z[3])/x);
	f[4] = -((F/(permittivity_0*rel_permittivity))*(z[2]-c_fixed)*pow((r_agg+delta_agg),2.0)) - ((2.0*z[5])/x);

	//6 equations
	//f[4] = z[5]/(permittivity_0*rel_permittivity);
	//f[5] = (-1.0*F*(z[2]-c_fixed)*pow((r_agg+delta_agg),2.0)) - ((2.0*z[5])/x);
}

//---------------------------------------------------------------------------
void
NAME::WaterAgglomerate::dfsub(double &x, double z[], double [], double df[])
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
		c_fixed =  core_charge_factor*c_H;

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
NAME::WaterAgglomerate::gsub (int &i, double z[], double &g)
{
// 	std::vector<double> c_h(1);
// 	c_h[0] = z[2];
    double T = this->solutions[temperature_of_REV][this->sol_index];
	double D_O2_EFF = pow(epsilon_agg,1.5) * DO2_Water;
	//double SIGMA_EFF = pow(epsilon_agg,1.5) * sigma_p;
	double D_H_EFF = pow(epsilon_agg,1.5) * D_H_Water; //SIGMA_EFF*R*T/(c_h[0]*pow(F,2.0));
	//double rel_permittivity = rel_permittivity_water;

	if (i == 1) g = z[1]/D_O2_EFF;
	else if (i == 2) g = z[3]/D_H_EFF - (F/(R*T))*z[2]*z[5];
	else if (i == 3) g = z[5];//(permittivity_0*rel_permittivity);
	else if (i == 4) g = z[0] - c_R;
	else if (i == 5) g = z[2] - c_H;
	else if (i == 6) g = z[4] - phi_M;
}

//---------------------------------------------------------------------------
void
NAME::WaterAgglomerate::dgsub (int &i, double z[], double dg[])
{
// 	std::vector<double> c_h(1);
// 	c_h[0] = z[2];
    double T = this->solutions[temperature_of_REV][this->sol_index];
	double D_O2_EFF = pow(epsilon_agg,1.5) * DO2_Water;
	//double SIGMA_EFF = pow(epsilon_agg,1.5) * sigma_p;
	double D_H_EFF = pow(epsilon_agg,1.5) * D_H_Water; //SIGMA_EFF*R*T/(c_h[0]*pow(F,2.0));
	//double rel_permittivity = rel_permittivity_water;

	for (int j=0; j < 6; j++) dg[j]=0.;

	if (i == 1) dg[1] = 1.0/D_O2_EFF;
	else if (i == 2) 
	{
		dg[2] = -(F/(R*T))*z[5];
		dg[3] = 1.0/D_H_EFF;
		dg[5] = -(F/(R*T))*z[2];
	}
	else if (i == 3) dg[5] = 1.0;//(permittivity_0*rel_permittivity);
	else if (i == 4) dg[0] = 1.0;
	else if (i == 5) dg[2] = 1.0;
	else if (i == 6) dg[4] = 1.0;
}

//---------------------------------------------------------------------------
void
NAME::WaterAgglomerate::guess (double &x, double z[], double [], double dz[])
{

    if(!use_initial_data(z, x))
    {
        //This initial guess strategy will only work for very high cell OCV
        double D_O2_EFF = D_O2_N;
        double E0 = 1.1;
        double T = this->solutions[temperature_of_REV][this->sol_index];
        double d_film = 1.0;
        double d_agg = 0.0;
        double poly = 5.0;

        for (int j=0; j < 6; j++)
        {
            z[j] = 0;
            dz[j] = 0;
        }

        if (x > interface)
        {
            double slope = (1.0-d_film)/(boundary_1-interface);
            z[0] = c_R * (slope*x + (1.0 - slope));
            dz[0] = slope;

            z[1] = 0.0;
            dz[1] = 0.0;
            z[2] = c_H;
            dz[2] = 0.0;
            z[3] = 0.0;
            dz[3] = 0.0;
            z[4] = phi_M;
            dz[4] = 0.0;
            z[5] = 0.0;
            dz[5] = 0.0;
        }

        else
        {
            double D_H_EFF = pow(epsilon_agg,1.5) * D_H_Water;
            D_O2_EFF =  pow(epsilon_agg,1.5) * DO2_Water;
            double slope = (d_film - d_agg)/pow(interface,poly);

            z[0] = c_R; //* (d_agg + slope*pow(x,poly));
            dz[0] = 0.0; //c_O2 * poly * slope *pow(x,poly-1);
            z[1] = dz[0]*D_O2_EFF;
            dz[1] = 0.0;

            z[2] = c_H * (d_agg + slope*pow(x,poly));
            dz[2] = c_H * poly * slope *pow(x,poly-1);
            z[3] = D_H_EFF * (dz[2] + (F/(R*T)*z[2]*dz[4]));
            dz[3] = 0.0;

            double step = 0.5 *(E0 - (phi_S - phi_M));
            poly = 3.0;
            double slope2 = -step/pow(interface,poly);
            z[4] = (phi_M+step) + slope2 * pow(x,poly);
            dz[4] = poly*slope2*pow(x,poly-1);
            z[5] = dz[4];
            dz[5] = poly*(poly-1)*slope2*pow(x,poly-2);
        }

        // if (lambda == 0.001)
        // {
        // 	std::string filename = "guess.txt";
        // 	std::stringstream inputss;
        // 	inputss << "x\t" << "O2\t" << "dO2\t" << "H\t" << "dH\t" << "V\t" << "dV\n";
        // 	inputss << x << "\t";
        // 	for (int i=0; i<6; ++i)
        // 		inputss << z[i] << "\t";
        //
        // 	std::fstream fout;
        // 	fout.open(filename.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
        // 	fout << inputss.str().c_str() << std::endl;
        // 	fout.close();
        // }
   }
}

//---------------------------------------------------------------------------
void
NAME::WaterAgglomerate::fsub_wrapper (double &x, double z[], double y[], double f[])
{
	WaterAgglomerate *ptr_Agglomerate = ((WaterAgglomerate*) (ptr_DAE_object[omp_get_thread_num()]));
	ptr_Agglomerate->fsub(x, z, y, f);
}

//---------------------------------------------------------------------------
void
NAME::WaterAgglomerate::dfsub_wrapper (double &x, double z[], double y[], double f[])
{
	WaterAgglomerate *ptr_Agglomerate = ((WaterAgglomerate*) (ptr_DAE_object[omp_get_thread_num()]));
	ptr_Agglomerate->dfsub(x, z, y, f);
}

//---------------------------------------------------------------------------
void
NAME::WaterAgglomerate::gsub_wrapper (int &i, double z[], double &g)
{
	WaterAgglomerate *ptr_Agglomerate = ((WaterAgglomerate*) (ptr_DAE_object[omp_get_thread_num()]));
	ptr_Agglomerate->gsub(i, z, g);
}

//---------------------------------------------------------------------------
void
NAME::WaterAgglomerate::dgsub_wrapper (int &i, double z[], double dg[])
{
	WaterAgglomerate *ptr_Agglomerate = ((WaterAgglomerate*) (ptr_DAE_object[omp_get_thread_num()]));
	ptr_Agglomerate->dgsub(i, z, dg);
}

//---------------------------------------------------------------------------
void
NAME::WaterAgglomerate::guess_wrapper (double &x, double z[], double y[], double df[])
{
	WaterAgglomerate *ptr_Agglomerate = ((WaterAgglomerate*) (ptr_DAE_object[omp_get_thread_num()]));
	ptr_Agglomerate->guess(x, z, y, df);
}



double NAME::WaterAgglomerate::compute_thickness_agg()
{
    //set default and starting values
       double thickness_0 = r_agg / 10.;
       const double tol = 1e-6;
       const double eps_tol = 1e-6;
       const double iterations = 1000;

       CL_Properties props = this->layer->get_properties();

       // Agglomerate is water-filled, there is no ionomer within the agglomerate.
       double porosity =  0.0;

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
                   throw std::runtime_error("Water filled agglomerate film thickness has converged to a negative value.");
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


double NAME::WaterAgglomerate::compute_epsilon_agg()
{
    AssertThrow(false, ExcMessage("Water filled agglomerate cannot have a constant Thickness"));

}



