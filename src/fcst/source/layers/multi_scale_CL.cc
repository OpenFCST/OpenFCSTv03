//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: multi_scale_CL.cc
//    - Description: Class characterizing the catalyst layer and defining effective properties
//    - Developers: M. Secanell, Peter Dobson, Philip Wardlaw, Michael Moore, Madhur Bhaiya
//    - $Id: multi_scale_CL.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include <layers/multi_scale_CL.h>

#ifdef _OPENMP
#include <omp.h>
#define PARALLEL 1
#define agg_threads() 9
#else
#define omp_get_thread_num() 0
#define agg_threads() 1
#endif

namespace NAME = FuelCellShop::Layer;

template<int dim>
const std::string NAME::MultiScaleCL<dim>::concrete_name("MultiScaleCL");

template<int dim>
NAME::MultiScaleCL<dim> const* NAME::MultiScaleCL<dim>::PROTOTYPE =
        new NAME::MultiScaleCL<dim>();

//---------------------------------------------------------------------------
template<int dim>
NAME::MultiScaleCL<dim>::MultiScaleCL() :
NAME::ConventionalCL<dim>() 
{

    this->get_mapFactory()->insert(
            std::pair<std::string, FuelCellShop::Layer::CatalystLayer<dim>*>(
                    concrete_name, this));
}

//---------------------------------------------------------------------------
template<int dim>
NAME::MultiScaleCL<dim>::MultiScaleCL(std::string cl_section_name) :
NAME::ConventionalCL<dim>(cl_section_name) 
{


}

//---------------------------------------------------------------------------
template<int dim>
NAME::MultiScaleCL<dim>::~MultiScaleCL() {
    micro.clear();
}

//---------------------------------------------------------------------------
template<int dim>
void NAME::MultiScaleCL<dim>::declare_parameters(

    const std::string& cl_section_name, ParameterHandler &param) const {
    // Declare a material class for liquid water:

    NAME::ConventionalCL<dim>::declare_parameters(cl_section_name, param);

    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection(cl_section_name);
        {
            param.enter_subsection(this->concrete_name);
            {
                FuelCellShop::MicroScale::MicroScaleBase::declare_MicroScale_parameters(param);
                //-- Agglomerate Structure
                param.declare_entry("Average current in cell", "false",
                        Patterns::Bool(),
                        "Decide whether to take the average current density in the cell");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
template<int dim>
void NAME::MultiScaleCL<dim>::initialize(ParameterHandler &param)
{
    NAME::ConventionalCL<dim>::initialize(param);

    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection(this->name);
        {
            param.enter_subsection(this->concrete_name);
            {
                average_cell_current = param.get_bool("Average current in cell");
                initialize_micro_scale(param);
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }

    param.leave_subsection();

    compute_void_fraction();
}


//---------------------------------------------------------------------------
template<int dim>
void NAME::MultiScaleCL<dim>::initialize_micro_scale(ParameterHandler &param) {

    micro.clear();

    //for(unsigned int i =0; i< agg_threads(); i++)
    //    micro.push_back(FuelCellShop::MicroScale::MicroScaleBase::create_MicroStructure(param,this));

    #ifdef _OPENMP
        FcstUtilities::log << "MultiScaleCL running in OpenMP mode." << std::endl;
    #endif

    for(unsigned int i = 0; i < this->material_ids.size(); i++){
        this->set_local_material_id(this->material_ids.at(i));
    	for(unsigned int j =0 ; j< agg_threads(); j++){
    	    micro[this->local_material_id()].push_back(FuelCellShop::MicroScale::MicroScaleBase::create_MicroStructure(param,this));
            #ifdef _OPENMP
    	    param.leave_subsection();
    	    micro[this->local_material_id()].at(j)->make_thread_safe(param,j);
    	    param.enter_subsection(this->concrete_name);
            #endif


    	}
    }

    //un-set id to invalid to discourage coders relying on some default unknown id...
    this->unset_local_material_id();

}

//---------------------------------------------------------------------------
template<int dim>
void NAME::MultiScaleCL<dim>::current_density(std::vector<double>& current) {
    std::vector<double> effectiveness(current.size());

    //call the current density function for the class
    current_density(current, effectiveness);
}

//---------------------------------------------------------------------------
template<int dim>
void NAME::MultiScaleCL<dim>::current_density(std::vector<double>& current,
        std::vector<double>& Er) {

    //update_agglomerate_guess();

    current.clear();
    current.resize(this->solutions[this->reactant].size());
    Er.resize(this->solutions[this->reactant].size());
    coverage_map.clear();
    std::vector<double> coverage_OH(this->solutions[this->reactant].size(), 0.0);
    std::vector<double> coverage_O(this->solutions[this->reactant].size(), 0.0);

    if (average_cell_current)
    {
        // ##### average the solution variables over the all quadrature points, use average to compute single value for current density, apply to all points

        std::map<VariableNames, SolutionVariable > averagedSol;
        double cell_current;

        double x_R = 0.;
        double Vm = 0.;
        double Vs = 0.;
        double lambda_cell = 0.;
        double Eff = 0.;
        double t_cell = 0.;
        //Take a simple average of the values in the cell
        for (unsigned int j = 0; j < this->solutions[this->reactant].size(); ++j) {
            x_R += this->solutions[this->reactant][j] / this->solutions[this->reactant].size();

            Vm +=  this->solutions[protonic_electrical_potential][j]
                                                                  / this->solutions[protonic_electrical_potential].size();

            Vs +=  this->solutions[electronic_electrical_potential][j]
                                                                    / this->solutions[electronic_electrical_potential].size();

            lambda_cell += this->solutions[membrane_water_content][j]
                                                                   / this->solutions[membrane_water_content].size();

            t_cell += this->solutions[temperature_of_REV][j]
                                                          / this->solutions[temperature_of_REV].size();
        }

        averagedSol[this->reactant] = SolutionVariable(x_R, 1, this->reactant);
        averagedSol[protonic_electrical_potential] = SolutionVariable(Vm, 1, protonic_electrical_potential);
        averagedSol[electronic_electrical_potential] = SolutionVariable(Vs, 1, electronic_electrical_potential);
        averagedSol[membrane_water_content] = SolutionVariable(lambda_cell, 1, membrane_water_content);
        averagedSol[temperature_of_REV]           = SolutionVariable(t_cell, 1, temperature_of_REV);
        SolutionMap s = micro_scale_current(averagedSol, 0, 0);
        cell_current = s.at(VariableNames::current_density)[0];
        Eff = s.at(VariableNames::CL_effectiveness)[0];

        double C_OH = 0.0; double C_O =0.0;
        if (s.has(VariableNames::OH_coverage))
            C_OH = s.at(VariableNames::OH_coverage)[0];

        if (s.has(VariableNames::O_coverage))
            C_O  = s.at(VariableNames::O_coverage)[0];

        for (unsigned int i = 0; i < current.size(); ++i) {
            current[i] = cell_current;
            Er[i] = Eff;

            coverage_OH[i] = C_OH;
            coverage_O[i] = C_O;

            if (std::isnan(current[i]))
                current[i] = 0.0;
        }

    }
    else
    {
        #pragma omp parallel for  shared(current, Er) num_threads(agg_threads())
        for (unsigned int i = 0; i < current.size(); ++i) {
            SolutionMap s = micro_scale_current(this->solutions, i, i);
            current[i] = s.at(VariableNames::current_density)[0];
            Er[i] = s.at(VariableNames::CL_effectiveness)[0];
            if (s.has(VariableNames::OH_coverage))
                coverage_OH[i] = s.at(VariableNames::OH_coverage)[0];

            if (s.has(VariableNames::O_coverage))
                coverage_O[i] = s.at(VariableNames::O_coverage)[0];

            if (std::isnan(current[i]))
                current[i] = 0.0;
        }


    }

    coverage_map.push_back(SolutionVariable(coverage_OH, OH_coverage));
    coverage_map.push_back(SolutionVariable(coverage_O, O_coverage));
}


//---------------------------------------------------------------------------
template<int dim>
FuelCellShop::SolutionMap
NAME::MultiScaleCL<dim>::micro_scale_current(std::map<VariableNames ,SolutionVariable>& solutionMap,
        const unsigned int& sol_index, const unsigned int& thread_index){

    #ifdef _OPENMP
        unsigned int idx = thread_index;
    #else
        unsigned int idx = 0;
    #endif

    //Generate solutions from micro scale
    micro.at(this->local_material_id()).at(idx)->set_solution(solutionMap, this->reactant, sol_index);
    SolutionMap answer = micro.at(this->local_material_id()).at(idx)->compute_current();


    //Make some additional checks when compiling in debug
    Assert(answer.has(VariableNames::current_density),
              ExcMessage("Micro scale object does not supply necessary solution for VariableNames::current_density."));
    Assert(answer.has(VariableNames::CL_effectiveness),
              ExcMessage("Micro scale object does not supply necessary solution for VariableNames::CL_effectiveness."));


    //Scale Current density
    SolutionVariable curr = answer.pop(VariableNames::current_density);
    answer.push_back(SolutionVariable(curr[0]*(1.0 - this->epsilon_V.at(this->local_material_id())),1,VariableNames::current_density));

    return answer;
}


//---------------------------------------------------------------------------
template<int dim>
void NAME::MultiScaleCL<dim>::solve_current_derivatives_average(
        std::map< VariableNames, std::vector<double> >& Dcurrent) {


    std::vector<double> dcurrent_cell(3, 0.0);
    double cell_current = 0.;
    double Er_dummy;
    double x_R = 0.;
    double Vm = 0.;
    double Vs = 0.;
    double lambda_cell = 0.;
    double t_cell = 0;

    //Take a simple average of the values in the cell
    for (unsigned int j = 0; j < this->solutions[this->reactant].size(); ++j) {
        x_R += this->solutions[this->reactant][j] / this->solutions[this->reactant].size(); //this->x_O2->at(j)/this->x_O2->size();

        Vm += this->solutions[protonic_electrical_potential][j]
                / this->solutions[protonic_electrical_potential].size(); //this->phi_m->at(j)/this->phi_m->size();

        Vs += this->solutions[electronic_electrical_potential][j]
                        / this->solutions[electronic_electrical_potential].size(); //this->phi_s->at(j)/this->phi_s->size();

        lambda_cell += this->solutions[membrane_water_content][j]
                / this->solutions[membrane_water_content].size(); //this->lambda->at(j)/this->lambda->size();

        t_cell += this->solutions[temperature_of_REV][j]
                / this->solutions[temperature_of_REV].size();

    }

    std::vector<std::map<VariableNames, SolutionVariable >> averagedSol;
    averagedSol.clear();

    std::map<VariableNames, SolutionVariable > averagedSolMap;
    averagedSolMap[this->reactant]                                 = SolutionVariable(x_R, 1, this->reactant);
    averagedSolMap[protonic_electrical_potential]    = SolutionVariable(Vm, 1, protonic_electrical_potential);
    averagedSolMap[electronic_electrical_potential]  = SolutionVariable(Vs, 1, electronic_electrical_potential);
    averagedSolMap[membrane_water_content]           = SolutionVariable(lambda_cell, 1, membrane_water_content);
    averagedSolMap[temperature_of_REV]           = SolutionVariable(t_cell, 1, temperature_of_REV);

    //Make copies if averagedSol for each thread
    for(unsigned int i =0; i< agg_threads(); i++){
        averagedSol.push_back(averagedSolMap);
    }



    //Derivative Pertubations
    double h = 1.0e-4; //TODO: This should not be hard coded
    double x_R_h = x_R + h;
    double Vm_h = Vm + h;
    double Vs_h = Vs + h;
    double x_R_h2 = x_R - h;
    double Vm_h2 = Vm - h;
    double Vs_h2 = Vs - h;

    //Used to store the location of the protonic and electronic derivatives.
    unsigned int location_phi_m = 0;
    unsigned int location_phi_s = 0;

    if (micro.at(this->local_material_id()).at(0)->has_derivatives()) {
        micro.at(this->local_material_id()).at(0)->set_solution(averagedSol.at(0), this->reactant, 0);
        dcurrent_cell = micro.at(this->local_material_id()).at(0)->compute_derivative_current(); //Derivatives are in the order x02, phi_s, phi_m
        VariableNames name;
        for (unsigned int i = 0; i < this->derivative_flags.size(); ++i) {
           if(i==0)
               name = this->reactant;
           else if(i==1)
               name = protonic_electrical_potential;
           else if(i==2)
               name = electronic_electrical_potential;

           #pragma omp parallel for  shared(Dcurrent) num_threads(agg_threads())
           for (unsigned int j = 0; j < this->solutions[this->reactant].size(); ++j) {
               Dcurrent[name][j] = dcurrent_cell[i] * (1.0 - this->epsilon_V.at(this->local_material_id())); //*(x_O2[j]/x);//*(kinetics_dcurrent[i][j]/ave_dcurrent[i]);// / Dcurrent[i].size();
               if (std::isnan(Dcurrent[name][j]))
                   Dcurrent[name][j] = 0.0;
           }
        }
    }
    else {

        #pragma omp parallel for  shared(Dcurrent)  num_threads(agg_threads())
        for (unsigned int i = 0; i < this->derivative_flags.size(); ++i) {

            if (this->derivative_flags[i] == this->reactant) {

                //Forward pertubation
                averagedSol.at(i%agg_threads())[this->reactant] =
                        SolutionVariable(x_R_h,1,this->reactant);
                SolutionMap s = micro_scale_current(averagedSol.at(i%agg_threads()), 0, i);
                double Dcurrent_node = s.at(VariableNames::current_density)[0];

                //Backward pertubation
                averagedSol.at(i%agg_threads())[this->reactant] =
                        SolutionVariable(x_R_h2,1,this->reactant);

                s = micro_scale_current(averagedSol.at(i%agg_threads()), 0 ,i);
                cell_current = s.at(VariableNames::current_density)[0];

                //Set value back to default averaged value
                averagedSol.at(i%agg_threads())[this->reactant]=
                        SolutionVariable(x_R,1,this->reactant);

                for (unsigned int j = 0; j < this->solutions[this->reactant].size(); ++j) {

                    Dcurrent[this->reactant][j] = ((Dcurrent_node - cell_current)
                            / (x_R_h - x_R_h2));
                    if (std::isnan(Dcurrent[this->reactant][j]))
                        Dcurrent[this->reactant][j] = 0.0;
                }
            }
            else if (this->derivative_flags[i] == electronic_electrical_potential) {


                //Forward pertubation
                averagedSol.at(i%agg_threads())[electronic_electrical_potential] =
                        SolutionVariable(Vs_h, 1,electronic_electrical_potential);

                SolutionMap s = micro_scale_current(averagedSol.at(i%agg_threads()), 0, i );
                double Dcurrent_node = s.at(VariableNames::current_density)[0];

                //Backward pertubation
                averagedSol.at(i%agg_threads())[electronic_electrical_potential] =
                        SolutionVariable(Vs_h2, 1,electronic_electrical_potential);

                s = micro_scale_current(averagedSol.at(i%agg_threads()), 0, i);
                cell_current = s.at(VariableNames::current_density)[0];
                //Set value back to default averaged value
                averagedSol.at(i%agg_threads())[electronic_electrical_potential] =
                        SolutionVariable(Vs, 1,electronic_electrical_potential);


                for (unsigned int j = 0; j < this->solutions[electronic_electrical_potential].size(); ++j) {

                    Dcurrent[electronic_electrical_potential][j] = ((Dcurrent_node - cell_current)
                            / (Vs_h - Vs_h2));
                    if (std::isnan(Dcurrent[electronic_electrical_potential][j]))
                        Dcurrent[electronic_electrical_potential][j] = 0.0;
                }

            } else //undefined case
            {
                for (unsigned int j = 0; j < this->solutions[this->reactant].size(); ++j)
                    Dcurrent[this->derivative_flags[i]][j] = 0.0;
            }
        }

        //Give the proton derivative the negative value of the electronic derivative
        for (unsigned int j = 0; j < this->solutions[electronic_electrical_potential].size(); ++j)
            Dcurrent[protonic_electrical_potential][j] = -Dcurrent[electronic_electrical_potential][j];

    }
}


//---------------------------------------------------------------------------
template<int dim>
void NAME::MultiScaleCL<dim>::solve_current_derivatives_at_each_node(
        std::map< VariableNames, std::vector<double> >& Dcurrent) {

    //update_agglomerate_guess();
    std::vector<double> cell_current(this->solutions[this->reactant].size(), 0.0);
    std::vector<std::vector<double> > d_current(this->derivative_flags.size(),
            std::vector<double>(this->solutions[this->reactant].size(), 0.0));


    std::vector<double> x_R_h(this->solutions[this->reactant].size(), 0.0);
    std::vector<double> phi_m_h(this->solutions[protonic_electrical_potential].size(), 0.0);
    std::vector<double> phi_s_h(this->solutions[electronic_electrical_potential].size(), 0.0);

    double Er_dummy;

    if (micro.at(this->local_material_id()).at(0)->has_derivatives()) {
        #pragma omp parallel for  shared(Dcurrent) num_threads(agg_threads())
        for (unsigned int j = 0; j < this->solutions[this->reactant].size(); ++j) {
            unsigned int idx = 0;


            #ifdef _OPENMP
            idx = j;
            #endif

            micro.at(this->local_material_id()).at(idx)->set_solution(this->solutions, this->reactant, j);
            std::vector<double> dcurrent_q_point = micro.at(this->local_material_id()).at(idx)->compute_derivative_current();
            for (unsigned int i = 0; i < this->derivative_flags.size(); ++i) {
                if (this->derivative_flags[i] == this->reactant)
                    Dcurrent[this->reactant][j] = dcurrent_q_point[0]
                            * (1.0 - this->epsilon_V.at(this->local_material_id()));
                else if (this->derivative_flags[i] == protonic_electrical_potential)
                    Dcurrent[protonic_electrical_potential][j] = dcurrent_q_point[1]
                            * (1.0 - this->epsilon_V.at(this->local_material_id()));
                else if (this->derivative_flags[i] == electronic_electrical_potential)
                    Dcurrent[electronic_electrical_potential][j] = dcurrent_q_point[2]
                            * (1.0 - this->epsilon_V.at(this->local_material_id()));
                else
                    Dcurrent[this->derivative_flags[i]][j] = 0.0;

                if (std::isnan(Dcurrent[this->derivative_flags[i]][j]))
                    Dcurrent[this->derivative_flags[i]][j] = 0.0;
            }
        }
    }
    else //Compute Numerical Derivatives using forward differencing
    {

        double h = 1.0e-4; //TODO: This should not be hard coded
        #pragma omp parallel for  shared(cell_current,x_R_h,phi_m_h,phi_s_h)  num_threads(agg_threads())
        for (unsigned int j = 0; j < this->solutions[this->reactant].size(); ++j) {
            //Compute mid points
            SolutionMap s = micro_scale_current(this->solutions, j, j);
            cell_current[j] = s.at(VariableNames::current_density)[0];

            //Fill Forward Pertubations
            x_R_h[j] = this->solutions[this->reactant][j] + h;
            phi_m_h[j] = this->solutions[protonic_electrical_potential][j] + h;
            phi_s_h[j] = this->solutions[electronic_electrical_potential][j] + h;

        }

        for (unsigned int i = 0; i < this->derivative_flags.size(); ++i) {
            std::map<VariableNames, SolutionVariable > perturbedSol = this->solutions; //Set the perturbed Solution map to the solution

            if (this->derivative_flags[i] == this->reactant) {

                //Set Pertubations for Oxygen_molar_fraction
                perturbedSol[this->reactant] = SolutionVariable(&x_R_h, this->reactant);
                #pragma omp parallel for  shared(Dcurrent) num_threads(agg_threads())
                for (unsigned int j = 0; j < this->solutions[this->reactant].size(); ++j) {
                    //Compute forward point for Oxygen_molar_fraction

                    SolutionMap s = micro_scale_current(perturbedSol, j, j);
                    double Dcurrent_node = s.at(VariableNames::current_density)[0];

                    //Compute Derivative
                    Dcurrent[this->reactant][j] = ((Dcurrent_node - cell_current[j])
                            / (x_R_h[j] - this->solutions[this->reactant][j]));
                    if (std::isnan(Dcurrent[this->reactant][j]))
                        Dcurrent[this->reactant][j] = 0.0;
                }

            }
            else if (this->derivative_flags[i] == protonic_electrical_potential) {

                //Set Pertubations for Protonic_electrical_potential
                perturbedSol[protonic_electrical_potential] = SolutionVariable(&phi_m_h, protonic_electrical_potential);

                #pragma omp parallel for  shared(Dcurrent) num_threads(agg_threads())
                for (unsigned int j = 0; j < this->solutions[protonic_electrical_potential].size(); ++j) {
                    //Compute forward point for Protonic_electrical_potential
                    SolutionMap s = micro_scale_current(perturbedSol, j, j);
                    double Dcurrent_node = s.at(VariableNames::current_density)[0];

                    //Compute Derivative
                    Dcurrent[protonic_electrical_potential][j] = ((Dcurrent_node - cell_current[j])
                            / (phi_m_h[j] - this->solutions[protonic_electrical_potential][j]));
                    if (std::isnan(Dcurrent[protonic_electrical_potential][j]))
                        Dcurrent[protonic_electrical_potential][j] = 0.0;
                }
            }
            else if (this->derivative_flags[i] == electronic_electrical_potential) {

                //Set Pertubations for Electronic_electrical_potential
                perturbedSol[electronic_electrical_potential] = SolutionVariable(&phi_s_h, electronic_electrical_potential);

                #pragma omp parallel for  shared(Dcurrent) num_threads(agg_threads())
                for (unsigned int j = 0; j < this->solutions[electronic_electrical_potential].size(); ++j) {
                    //Compute forward point for Electronic_electrical_potential
                    SolutionMap s = micro_scale_current(perturbedSol, j, j);
                    double Dcurrent_node = s.at(VariableNames::current_density)[0];

                    //Compute Derivative
                    Dcurrent[electronic_electrical_potential][j] = ((Dcurrent_node - cell_current[j])
                            / (phi_s_h[j] - this->solutions[electronic_electrical_potential][j]));
                    if (std::isnan(Dcurrent[electronic_electrical_potential][j]))
                        Dcurrent[electronic_electrical_potential][j] = 0.0;
                } //For all quadrature points
            } //If this variable flag
        } //For all derivative flags
    }// If numerical derivatives
}


//---------------------------------------------------------------------------
template<int dim>
void NAME::MultiScaleCL<dim>::derivative_current_density(
        std::map< VariableNames, std::vector<double> >& Dcurrent) {



    //update_agglomerate_guess();

    for( unsigned int i =0; i < this->derivative_flags.size(); i++){
        Dcurrent[this->derivative_flags[i]] = std::vector<double>(this->solutions[this->reactant].size(), 0.0);
    }

    if (average_cell_current)
        solve_current_derivatives_average(Dcurrent);
    else
        solve_current_derivatives_at_each_node(Dcurrent);

}


//---------------------------------------------------------------------------
template<int dim>
void NAME::MultiScaleCL<dim>::compute_void_fraction() {
    // Porosity
	 for(unsigned int i = 0; i < this->material_ids.size(); i++){
	     this->set_local_material_id(this->material_ids.at(i));
	     this->epsilon_V[this->local_material_id()] = 1. - this->epsilon_S.at(this->local_material_id()) -
	             this->epsilon_N.at(this->local_material_id())-
	             micro.at(this->local_material_id()).at(0)->aux_volume_fraction();
	 }
	 //unset id to invalid to discourage coders rellying on some default unknown id...
	 this->unset_local_material_id();
}


//---------------------------------------------------------------------------
template <int dim>
FuelCellShop::SolutionMap
NAME::MultiScaleCL<dim>::get_coverages(){
    return coverage_map;
}


//---------------------------------------------------------------------------
template<int dim>
void NAME::MultiScaleCL<dim>::print_layer_properties() const {
    FuelCellShop::Layer::ConventionalCL<dim>::print_layer_properties();

    for(unsigned i = 0; i< this->material_ids.size(); i++){

    	FcstUtilities::log<<"-----Agglomerate for sub layer (" << i +1 << ")-----" <<std::endl;
    	micro.at(this->material_ids.at(i)).at(0)->print_properties();
    	FcstUtilities::log<<"---------------------------------------"<<std::endl;
    }


}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations.
template class NAME::MultiScaleCL<deal_II_dimension>;
