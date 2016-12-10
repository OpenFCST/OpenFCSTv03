//---------------------------------------------------------------------------
// C++ Interface: agglomerate_ionomer_1D.cc
//
// Description: Used to solve a system of equations representing a
// 			spherical ionomer-filled agglomerate.
//
// Author: Jason Boisvert <jjb701@mail.usask.ca>, (C) 2011
//			University of Saskatchewan
//
//		Adapted by Peter Dobson <pdobson@ualberta.ca>,
//		University of Alberta for implementation in 
//		the Fuel Cell Simulation Toolbox (FCST).
//
// Copyright: See COPYING file that comes with this distribution
//
//---------------------------------------------------------------------------

#include <microscale/agglomerate_ionomer_1D.h>
#include <typeinfo>
#include <chrono>

#ifdef _OPENMP
	#include <omp.h>
#else
	#define omp_get_thread_num() 0
#endif
// #define _1_EQ_
namespace NAME = FuelCellShop::MicroScale;


const std::string NAME::IonomerAgglomerate::concrete_name("IonomerAgglomerateNumerical");
NAME::IonomerAgglomerate const* NAME::IonomerAgglomerate::PROTOTYPE = new NAME::IonomerAgglomerate(concrete_name);


NAME::IonomerAgglomerate::IonomerAgglomerate(std::string name){
    //FcstUtilities::log <<" Register " + concrete_name +  " to FactoryMap"<<std::endl;
       this->get_mapFactory()->insert(
               std::pair<std::string, NAME::IonomerAgglomerate*>(
                       concrete_name, this));
}

//---------------------------------------------------------------------------
NAME::IonomerAgglomerate::IonomerAgglomerate(int verbose)
:
endTol (1.e-5),
uFD(true)
{
    R_tol = 1.e-20;
	verbosity(1);
	this->has_derivatives_ = false;

	column_names.push_back("x"); column_names.push_back("z0");
	column_names.push_back("z1");  column_names.push_back("z2");
	column_names.push_back("z3"); column_names.push_back("i");
	prev_effectiveness = 0.5;
}

//---------------------------------------------------------------------------
void
NAME::IonomerAgglomerate::set_structure()
{

    CL_Properties props = this->layer->get_properties();

    electrolyte = this->layer->get_resource<FuelCellShop::Material::PolymerElectrolyteBase>();
    kinetics = this->layer->get_resource<FuelCellShop::Kinetics::BaseKinetics>();
    AV = props[CLPropNames::active_area_scaled]/cm_to_m; // convert to 1/m


    r_agg *=1e-9;//convert to m
    delta_agg *=1e-9; // convert to m
    interface = r_agg/(r_agg+delta_agg);// set to r_agg to remove
    AV /= pow(interface,3.0); //new scaling

}

//---------------------------------------------------------------------------
FuelCellShop::SolutionMap
NAME::IonomerAgglomerate::compute_current ( )
{


    double E_r = 0.0;

    //True pressure is now available from the layer.
    //Was 0 at set structure time.
    P= this->layer->get_properties()[CLPropNames::pressure];


    if(this->solutions[this->reactant][this->sol_index] <= 0.0){
        SolutionMap sols;
        sols.push_back(SolutionVariable(0, 1, current_density));
        sols.push_back(SolutionVariable(0,1, CL_effectiveness));

        return sols; //To stabalize FEM solution
    }


    electrolyte->set_T(this->solutions[temperature_of_REV][this->sol_index]);
    electrolyte->set_lambda(this->solutions[membrane_water_content][this->sol_index]);
    electrolyte->proton_conductivity(sigma_p);

    //Get Reactant diffusivity and henry's constant
    if (this->solutions[this->reactant].get_variablename() == oxygen_molar_fraction)
    {
        electrolyte->oxygen_diffusivity(D_R_N);
        H_R_N = electrolyte->get_H_O2();
        tempReactantName = oxygen_concentration;
    }
    else if(this->solutions[this->reactant].get_variablename() == hydrogen_molar_fraction)
    {
        AssertThrow(false, ExcMessage("Hydrogen reaction not yet implemented in numerical ionomer agglomerate."));
        //electrolyte->hydrogen_diffusivity(D_R_N);
        //H_R_N = electrolyte->get_H_H2();
        //tempReactantName = hydrogen_concentration;
    }
    SolutionVariable temperature(this->solutions[temperature_of_REV][this->sol_index],2, VariableNames::temperature_of_REV); //TODO Address hard coded '2'
    kinetics->set_temperature(temperature);
    kinetics->set_derivative_flags(sol_names);
    kinetics->set_p_t(P);
    sigma_p = sigma_p/cm_to_m; // convert to S/m
    D_R_N = D_R_N*cm2_to_m2; // convert to m^2/s
    H_R_N = H_R_N*cm3_to_m3; // convert to Pa m^3/mol

    sigma_p *= this->cond_factor;

    // Set the oxygen concentration and membrane potential at the boundary
    // Set the solid phase potential for the domain
    c_R = (this->solutions[this->reactant][this->sol_index]*P)/H_R_N; //c_O2 will now be in mol/m^3
    phi_S = this->solutions[electronic_electrical_potential][this->sol_index];
    phi_M = this->solutions[protonic_electrical_potential][this->sol_index];



    update_initial_solution();




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
    prob->set_integer_space(1e5);//was at 100000 1e5
    prob->set_float_space(1e6);//was at 1000000 1e6
//     if (this->non_equil_bc)
//         uFD = false;
//     else
        uFD = true;

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
            FcstUtilities::log << "First attempt failed, using a continuation strategy.\n";
        #endif

        clear_memory();
        if (cont_tolerance(1.e-2,endTol))
            obtained_solution = true;
        else
            obtained_solution = false; //all continuation strategies failed
    }

    //Deal with consequences of solver

    double unit_size = boundary_1 - boundary_0;// + delta_agg/(r_agg); //####(no thin film)
    double agg_size = interface - boundary_0;
    double volume = (4.0/3.0)*pi*pow(unit_size,3.0);
    double agg_volume = (4.0/3.0)*pi*pow(agg_size,3.0);


    double I = 0.0; // Integral of current over the mesh
    double C_OH = 0.0;
    double C_O = 0.0;

    double I_max = 0.0;

    std::vector<double> J(1,0.0);


    std::vector<SolutionVariable> c_reactants;
    c_reactants.push_back(SolutionVariable (c_R*cm3_to_m3,1,tempReactantName));
    SolutionVariable v_membrane(phi_M, 1, protonic_electrical_potential);
    SolutionVariable v_solid(phi_S, 1, electronic_electrical_potential);

    this->kinetics->set_reactant_concentrations(c_reactants);
    this->kinetics->set_electrolyte_potential(v_membrane);
    this->kinetics->set_solid_potential(v_solid);
    this->kinetics->current_density(J);

    I_max = (AV*cm_to_m) * J[0] * agg_volume/volume; //convert the current obtained from kinetics to A/m^2


    if (obtained_solution)
    {
        int n_mesh_final = prob->get_size_final_mesh();
        mesh = new double[n_mesh_final];
        prob->get_copy_final_mesh(mesh);


        //|| Get the total current produced by the agglomerate ||//
        // loop over each subinterval, defined by the mesh points
        // i.e. Composite integration.
        // NOTE: Must be implemented in each class... values are very problem specific.

        final_results.clear();
        std::vector<std::vector<double>> new_lineS;
        std::vector<double> new_line;

        //Temp Variables
        std::vector<double> temp_reactant;
        std::vector<double> temp_v_membrane;
        std::vector<double> temp_v_solid;

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
            J.clear();
            temp_reactant.clear(); temp_v_membrane.clear(); temp_v_solid.clear();
            J.resize(x_quad.size());

            temp_reactant.resize(x_quad.size());
            temp_v_membrane.resize(x_quad.size());
            temp_v_solid.resize(x_quad.size());

            new_lineS.clear();
            //Get data for each quad point per cell
            for (unsigned int i=0; i<x_quad.size(); ++i)
            {
                new_line.clear();
                double *z = new double[m_star];
                double *y = new double[n_y];
                prob->DAE_solution(x_quad[i],z,y); //Solution is returned in z.
                temp_reactant[i] = z[0]*cm3_to_m3; //convert back to cm before passing to kinetics class
                temp_v_membrane[i] = z[2];
                temp_v_solid[i] = phi_S;


                new_line.push_back(x_quad[i]); new_line.push_back(z[0]); new_line.push_back(z[1]);
                new_line.push_back(z[2]); new_line.push_back(z[3]);
                new_lineS.push_back(new_line);

                delete [] z;
                delete [] y;

            }


            //Compute the current
            c_reactants[0] = SolutionVariable(temp_reactant, tempReactantName);
            v_membrane = SolutionVariable(temp_v_membrane, protonic_electrical_potential);
            v_solid = SolutionVariable(temp_v_solid, electronic_electrical_potential);

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

            //Scale by Active area and CL volume, setup for integral over a sphere
            for (unsigned int i=0; i<x_quad.size(); ++i)
            {
                if (temp_reactant[i] < R_tol or ((right > interface) or ((j+1) == (n_mesh_final))))
                    J[i] = 0.0;
                else
                    J[i] = getAV(x_quad[i]) * (J[i]/cm2_to_m2 ) * (4.0*pi*pow(x_quad[i],2.0));//convert the current obtained from kinetics to A/m^2
            }

            //Save the results
            for(unsigned int z =0; z < new_lineS.size(); z++){
                new_lineS.at(z).push_back(J[z]);
                final_results.push_back(new_lineS.at(z));
            }

            // Add the contribution of the spherical shell
            I += integrate(left,right,w_quad,J);
            C_O += integrate(left,right,w_quad,c_O);
            C_OH += integrate(left,right,w_quad,c_OH);
        }




        save_initial_solution();

        double volume_cm = volume/cm3_to_m3;
        double I_avg = I/volume_cm;
        delete [] mesh;
        clear_memory();

        E_r = I_avg/I_max;
		prev_effectiveness = E_r;
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
        sols.push_back(SolutionVariable(I_max*prev_effectiveness, 1, current_density));
        sols.push_back(SolutionVariable(E_r,1, CL_effectiveness));

        DAE_Error(flag);
        return sols;
//
    }
}

//---------------------------------------------------------------------------
int 
NAME::IonomerAgglomerate::cont_tolerance (double start_tol, double end_tol)
{
	
	
	//Attempt to solve problem using continuation on tolerances
	int flag=0;
	double new_tol=0;
	int solve_attempts=0;
	bool obtained_sol = false;
	int save_attempts=0;
// 	if (this->non_equil_bc)
// 	    uFD = false;
// 	else
	    uFD = true;
	setup_DAE_solver();
 	prob->set_integer_space(1e6);// was at 10000000, 1e7
 	prob->set_float_space(1e7);// was at 100000000, 1e8
	//Set initial tolerances
	new_tol=1e-3;
	for (int i=0; i<n_comp;i++) tol[i] = new_tol;
	flag = prob->DAE_solve();
	solve_attempts=0;
	do 
	{
		//Apply continuation	
		
		if (flag <= 0)
		{
			//this continuation strategy has failed
			//FcstUtilities::log << "CONT ENDED TOO SOON AT: " << tol[1] << std::endl;
                        //FcstUtilities::log << "save attempt tol: " << new_tol << std::endl;
			new_tol = new_tol*2.;
			save_attempts++;	
			if (save_attempts > 80)
			{
				clear_memory();

                #ifdef DEBUG
				    FcstUtilities::log << "tolerance continuation strategy failed (saved)." << std::endl;
                #endif

				return 0;
			}
			prob->use_simple_cont();
		} 
		else
		{
			save_attempts=0;
			new_tol = new_tol/4.; //pow(1.05,++solve_attempts);
			solve_attempts++;
			//FcstUtilities::log << "global attempt tol: " << new_tol << std::endl;
			if(solve_attempts > 40)
			{
				clear_memory();
                #ifdef DEBUG
                   FcstUtilities::log  << "tolerance continuation strategy failed (global)." << std::endl;
                #endif
				return 0;
			}
			prob->use_simple_cont();
		}
		
		for (int i=0; i<n_comp;i++) tol[i] = new_tol;
		if (tol[0] < 1e-5 /*|| this->non_equil_bc*/) uFD = false;
		else uFD = true;
		flag = prob->DAE_solve();
		
	}while(new_tol >= end_tol ||flag <=0 );

// 	FcstUtilities::log << "Tol: " << new_tol << std::endl;
    #ifdef DEBUG
	    FcstUtilities::log << "Continuation worked " << std::endl;
	#endif

	return 1;
}

//---------------------------------------------------------------------------
void
NAME::IonomerAgglomerate::setup_DAE_solver ()
{
	//Declare an array that holds the number and order of each of the ODEs
	#ifdef _1_EQ_
	n_comp = 2;
	#else
	n_comp = 4;
	#endif 

	n_y = 0;
	mm = new int[n_comp];
	mm[0] = 1;
	mm[1] = 1;
	mm[2] = 1;
	mm[3] = 1;
	m_star = 0;
	for (int i=0; i<n_comp;++i)
		m_star += mm[i];

	// Set the mesh parameters
	n_mesh = 100;
	n_colloc = 4;
	boundary_0 = 0.0; 
	boundary_1 = 1.0;
	
// 	FcstUtilities::log << "interface " << interface << std::endl;
	
	bool set_fixpoint = false;
	if (interface < 1.0) set_fixpoint = true;
	fixpnt = new double[1]; // Use a mesh with a fixpoint
	fixpnt[0] =  interface;
	
	//declare the location of each boundary point
	zeta = new double[4];
	zeta[0] = boundary_0;
	#ifdef _1_EQ_
	zeta[1] = boundary_1;
	#else
	zeta[1] = boundary_0;
	zeta[2] = boundary_1;
	zeta[3] = boundary_1;
	#endif
	
	ltol = new int[n_comp + n_y];
	tol = new double[n_comp + n_y];
	for (int i=0; i<n_comp;++i)
	{
		ltol[i] = i+1;
		tol[i] = endTol;
	}

	
	
	// WARNING!!  This must be done in other implementations using the DAEWrapper class
	//	Failing to do so will cause runtime errors that are difficult to locate.
	// Assign the global object pointer to the IonomerAgglomerate object!!!
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
	These are methods used to set *certain parameters of COLDAE.
	They must be called before the problem is actually solved. 
	These methods are optional and are not required to be used 
	for every problem. 
	*/
	prob->set_solver_control(-1); //Newton solver control
	prob->set_collocation_points(n_colloc); //Set the number of collocation points 
	prob->set_initial_mesh_size(n_mesh); //Set the initial mesh size
	if (set_fixpoint) prob->set_fixpnts(1, fixpnt);
	prob->set_tolerance(n_comp,ltol,tol);
        //n_output
	prob->set_output(n_output); //Set the output level 
	//-1 is full output
	//0 is selected output
	//1 is no output
	
	//Limit our meshsize for the first attempt
	//A larger mesh will have to be used for continuation
	prob->set_integer_space(1e6);//was at 100000000 1e8
	prob->set_float_space(1e7);//was at 100000000  1e8
}

//---------------------------------------------------------------------------
void
NAME::IonomerAgglomerate::fsub(double &x, double z[], double y[], double f[])
{


	double epsilon = 1.0; //******
	std::vector<double> J(1,0.0);

	if (x <= interface)
	{
	    if ((z[0]*cm3_to_m3)>R_tol)
	    {
	        std::vector<SolutionVariable> c_reactants;
	        c_reactants.push_back(SolutionVariable (z[0]*cm3_to_m3, 1, tempReactantName));
	        SolutionVariable v_membrane(z[2], 1, protonic_electrical_potential);
	        SolutionVariable v_solid(phi_S, 1, electronic_electrical_potential);

	        this->kinetics->set_reactant_concentrations(c_reactants);
	        this->kinetics->set_electrolyte_potential(v_membrane);
	        this->kinetics->set_solid_potential(v_solid);
	        this->kinetics->current_density(J);

	        J[0] = getAV(x)*(J[0] * lambda)/cm2_to_m2;//convert to m before passing to COLDAE
	    }
		epsilon = pow(epsilon_agg,1.5);//******
	}
	
	double molarNumerator;

	if(this->reactant == oxygen_molar_fraction)
	{
	    molarNumerator = 4.0;
	}
	else if (this->reactant == hydrogen_molar_fraction)
	{
	    molarNumerator = 2.0;
	}
	else
	    AssertThrow(false, ExcMessage("Ionomer agglomerate cannot solve for this type of reactant!"));
                
	f[0] = z[1]/epsilon;
	f[1] = (( ( pow((r_agg +delta_agg),2.0)) *J[0])/(molarNumerator*F*D_R_N))-((2.0*z[1])/(x)); // consumption must be proportional to z[0]

    #ifndef _1_EQ_
        f[2] = (z[3])/(epsilon);
        f[3] = (pow((r_agg +delta_agg),2.0)*J[0])/(sigma_p)-((2.0*z[3])/x); // ####thin film
    #endif

}

//---------------------------------------------------------------------------
void
NAME::IonomerAgglomerate::dfsub(double &x, double z[], double [], double df[])
{
	//size of Jacobian
	int n_eq = n_comp + n_y;
	int n_var = m_star + n_y;

	//variables for finite differences 
	int i = 0; 
	int j = 0;
	double dummy=0;
	double *w1 = new double[n_eq];
	double *w2 = new double[n_eq];

	double esp=1.e-7; // finite stepsize
   
    //create space for Jacobian
    double** dfc;
	dfc = new double*[n_eq];
	//For each ptr in the array, assign it an array of 4 doubles. 
	for ( i=0; i < n_eq; ++i)
		dfc[i] = new double[n_var];
	// Initialize the matrix with zeros.
	for (i=0; i < n_eq; i++)
		for (int j=0; j <n_var; j++)
			dfc[i][j]= 0.0;


    //Construct Jacobian
	for (i=0;i<n_var;i++)
	{
		double temp;
		temp = z[i];
		//std << z[i] << "\n";
	 	z[i] = z[i] + esp;
		fsub(x,z,&dummy,w1);
		z[i] = z[i]-2.0*esp;
		fsub(x,z,&dummy,w2);
		z[i]=temp;
		for (j=0;j < n_eq; j++)
		{
			dfc[j][i] = (w1[j]-w2[j])/(2.0*esp);
		}
	}

	//create space for Jacobian
	double** dfc1;
	dfc1 = new double*[n_eq];
	//For each ptr in the array, assign it an array of 4 doubles. 
	for ( i=0; i < n_eq; ++i)
	    dfc1[i] = new double[n_var];
	// Initialize the matrix with zeros.0
	for (i=0; i < n_eq; i++)
	    for (int j=0; j <n_var; j++)
	        dfc1[i][j]= 0.0;
		
	double epsilon = 1.0; //******
	//std::vector<double> c_o2(1);
	//std::vector<double> phi_m(1);
	//std::vector<double> phi_s(1);
	//double HO2N = electrolyte->get_H_O2();

	//phi_m[0] = z[2];//phi_M;
	//phi_s[0] = phi_S;

	double PARTIALJ0(0.0);
	double PARTIALJ2(0.0);
	
	if (x <= interface)
	{
		epsilon = pow(epsilon_agg,1.5);
		if ((z[0]*cm3_to_m3)>R_tol)
		{
		    std::vector<SolutionVariable> c_reactants;
		    c_reactants.push_back(SolutionVariable (z[0]*cm3_to_m3, 1, tempReactantName));
		    SolutionVariable v_membrane(z[2], 1, protonic_electrical_potential);
		    SolutionVariable v_solid(phi_S, 1, electronic_electrical_potential);

		    this->kinetics->set_reactant_concentrations(c_reactants);
		    this->kinetics->set_electrolyte_potential(v_membrane);
		    this->kinetics->set_solid_potential(v_solid);
		    std::map< VariableNames, std::vector<double> > derivatives;
			this->kinetics->derivative_current(derivatives);

			PARTIALJ0 = (getAV(x) * (derivatives[this->reactant][0] *(H_R_N/P) *cm_to_m) * lambda);//convert to concentration in m before passing to COLDAE (Tafel returns wrt mole fraction)
			PARTIALJ2 = (getAV(x) * (derivatives[protonic_electrical_potential][0]/cm2_to_m2) * lambda);//convert to m before passing to COLDAE
		}
	}
	
	double molarNumerator;
    if(this->reactant == oxygen_molar_fraction)
    {
        molarNumerator = 4.0;
    }
    else if (this->reactant == hydrogen_molar_fraction)
    {
        molarNumerator = 2.0;
    }
    else
        AssertThrow(false, ExcMessage("Ionomer agglomerate cannot solve for this type of reactant!"));


	dfc1[0][0] = 0.0;
	dfc1[0][1] = 1.0/epsilon;
    #ifndef _1_EQ_
	    dfc1[0][2] = 0.0;
	    dfc1[0][3] = 0.0;
    #endif

	dfc1[1][0] = ((pow((r_agg+delta_agg),2.))/(molarNumerator*F*D_R_N)) * PARTIALJ0;
	dfc1[1][1] = -((2.0)/(x));


    #ifndef _1_EQ_
        dfc1[1][2] = ((pow((r_agg+delta_agg),2.))/(molarNumerator*F*D_R_N)) * PARTIALJ2;
        dfc1[1][3] = 0.0;

        dfc1[2][0] = 0.0;
        dfc1[2][1] = 0.0;
        dfc1[2][2] = 0.0;
        dfc1[2][3] = 1.0/(epsilon);

        dfc1[3][0] = (pow((r_agg+delta_agg),2.) * PARTIALJ0)/(sigma_p);
        dfc1[3][1] = 0.0;
        dfc1[3][2] = (pow((r_agg+delta_agg),2.)  * PARTIALJ2)/(sigma_p);
        dfc1[3][3] = -((2.0)/(x));

    #endif

    // 	Convert the matrix dfc to a one dimensional array that Fortran will understand.
	if (uFD)
	    FuelCell::ApplicationCore::c_to_for_matrix(n_eq,n_var,dfc,df);
	else
	    FuelCell::ApplicationCore::c_to_for_matrix(n_eq,n_var,dfc1,df);

	// Free the memory
	for (int i=0; i < n_eq; ++i)
		delete [] dfc[i];
	delete [] dfc;
	for (int i=0; i < n_eq; ++i)
		delete [] dfc1[i];
	delete [] dfc1;	
	
	delete [] w1;
	delete [] w2;

}

//---------------------------------------------------------------------------
void
NAME::IonomerAgglomerate::gsub (int &i, double z[], double &g)
{
	double DIFF = pow(epsilon_agg,1.5);// * D_O2_N;
	double SIGMA_EFF = pow(epsilon_agg,1.5);// * sigma_p;
	
	if (i == 1) g = z[1];
	else if (i == 2) g = z[3];
	else if (i == 3)
	{
		if (this->non_equil_bc)
			g = z[0] - c_R + (z[1]*D_R_N)/(r_agg*this->k_O2); //reverse dimensionalization from N_O2^bar to N_O2
		else
			g = z[0] - c_R ;// +  (z[1]*delta_agg)/(r_agg) ;//####(no thin film)
	}
	else g = z[2] - phi_M;// +  (z[3]*delta_agg)/r_agg;//####(no thin film)
}

//---------------------------------------------------------------------------
void
NAME::IonomerAgglomerate::dgsub (int &i, double [], double dg[])
{
	double DIFF = pow(epsilon_agg,1.5);// * D_O2_N;
	double SIGMA_EFF = pow(epsilon_agg,1.5);// * sigma_p;
	
	for (int j=0; j < 4; j++) dg[j]=0.;
	
	if (i == 1) dg[1] = 1.0;
	else if (i == 2) dg[3] = 1.0;
	else if (i == 3)
	{
		dg[0] = 1.0;
		if (this->non_equil_bc) 
		{
			dg[1] = (1.0*D_R_N)/(r_agg*this->k_O2);
		}
		else
		{
			if (uFD)
				dg[1] =  0.0; //Jason:: pertubation for convergence
			else
				dg[1] = 0.0;//
		}
	}
	else 
	{
		dg[2] = 1.0;
// 				if (uFD) dg[3] = 0.0; // Jason :: pertubation for convergence
// 					else dg[3] = 0.0;
	}	
}


//---------------------------------------------------------------------------
void
NAME::IonomerAgglomerate::guess (double &x, double z[], double [], double df[])
{
    if(!use_initial_data(z, x))
    {

        double SIGMA_EFF = sigma_p;
        double DIFF = D_R_N;

        double E0 = 1.1; // Approximately - add set OCV later
        double eta = E0 - (phi_S - phi_M);
        double d_film;
        double d_agg;
        double poly;

        if(eta <= 0.3)
        {
            d_film = 1.0;
            d_agg = 1.0;
            poly = 1.0;
        }
        else if(eta <= 0.6)
        {
            d_film = 0.9;
            d_agg = 0.5;
            poly = 3.0;
        }
        else
        {
            d_film = 0.3;
            d_agg = 0.0;
            poly = 9.0;
        }


        if (x > interface)
        {
            DIFF =  pow(epsilon_agg,1.5) * D_R_N;
            SIGMA_EFF = pow(epsilon_agg,1.5) * sigma_p;
            double slope = (1.0-d_film)/(boundary_1-interface);
            z[0] = c_R * (slope*x + (1.0 - slope)) ;
            df[0] = slope;
        }
        else
        {
            double slope = (d_film - d_agg)/pow(interface,poly);
            z[0] = c_R * (d_agg + slope*pow(x,poly));
            df[0] = c_R * poly * slope * pow(x,poly-1);
        }

        z[1] = 4.*df[0]*pow(epsilon_agg,1.5);// 0.1;
        df[1] = 0.0;
        z[2] = phi_M;
        df[2] = 0.0;
        z[3] = 0.1;
        df[3] = 0.0;
    }



}


//---------------------------------------------------------------------------
void
NAME::IonomerAgglomerate::fsub_wrapper (double &x, double z[], double y[], double f[])
{
	IonomerAgglomerate *ptr_Agglomerate = ((IonomerAgglomerate*)  (ptr_DAE_object[omp_get_thread_num()]));
	ptr_Agglomerate->fsub(x, z, y, f);
}


//---------------------------------------------------------------------------
void
NAME::IonomerAgglomerate::dfsub_wrapper (double &x, double z[], double y[], double f[])
{
	IonomerAgglomerate *ptr_Agglomerate = ((IonomerAgglomerate*) (ptr_DAE_object[omp_get_thread_num()]));;
	ptr_Agglomerate->dfsub(x, z, y, f);
}


//---------------------------------------------------------------------------
void
NAME::IonomerAgglomerate::gsub_wrapper (int &i, double z[], double &g)
{
	IonomerAgglomerate *ptr_Agglomerate = ((IonomerAgglomerate*) (ptr_DAE_object[omp_get_thread_num()]));;
	ptr_Agglomerate->gsub(i, z, g);
}


//---------------------------------------------------------------------------
void
NAME::IonomerAgglomerate::dgsub_wrapper (int &i, double z[], double dg[])
{
	IonomerAgglomerate *ptr_Agglomerate = ((IonomerAgglomerate*)(ptr_DAE_object[omp_get_thread_num()]));;
	ptr_Agglomerate->dgsub(i, z, dg);
}


//---------------------------------------------------------------------------
void
NAME::IonomerAgglomerate::guess_wrapper (double &x, double z[], double y[], double df[])
{
	IonomerAgglomerate *ptr_Agglomerate = ((IonomerAgglomerate*) (ptr_DAE_object[omp_get_thread_num()]));;
	ptr_Agglomerate->guess(x, z, y, df);
}
