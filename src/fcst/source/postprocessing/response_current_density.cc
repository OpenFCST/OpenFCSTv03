// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: responses.h
// - Description: This is a file containing the declaration of several FCST response evaluators
// - Developers: Marc Secanell Gallart,    University of Alberta
// - $Id: response_current_density.cc 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#include <postprocessing/response_current_density.h>

namespace NAME = FuelCellShop::PostProcessing;

//---------------------------------------------------------------------------

template <int dim>
void 
NAME::ORRCurrentDensityResponse<dim>::declare_parameters(ParameterHandler& param) const
{
    param.enter_subsection("PostProcessing"); 
    {
        param.enter_subsection("ORRCurrentDensityResponse"); 
        {
            param.declare_entry ("Catalyst layer surface area [cm^2]",  // volume fraction of void space in the CL
                                 "1.0", // [-]
                                 Patterns::Double(),
                                 "Active area of the CL. Used to scale the volumetric current density");
            param.declare_entry ("Catalyst layer volume [cm^3]",  // volume fraction of void space in the CL
                                 "1.0", // [-]
                                 Patterns::Double(),
                                 "Volume of the CL. Used to average the coverages");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::ORRCurrentDensityResponse<dim>::initialize(ParameterHandler&param )
{
    
    param.enter_subsection("PostProcessing"); 
    {
        param.enter_subsection("ORRCurrentDensityResponse"); 
        {
            S_CL = param.get_double("Catalyst layer surface area [cm^2]");
            V_CL = param.get_double("Catalyst layer volume [cm^3]");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    if ( this->system_management->solution_in_userlist("oxygen_molar_fraction") )
    {
        xi.solution_index = this->system_management->solution_name_to_index("oxygen_molar_fraction"); 
        xi.fetype_index = this->system_management->block_info->base_element[xi.solution_index];
        xi.indices_exist = true;
    }
    /*
    // Not necessarily the case if multi-component used.
    else
        throw std::runtime_error("oxygen_molar_fraction variable required for ORRCurrentDensityResponse");
    */
    
    if ( this->system_management->solution_in_userlist("protonic_electrical_potential") )
    {
        phiM.solution_index = this->system_management->solution_name_to_index("protonic_electrical_potential");
        phiM.fetype_index = this->system_management->block_info->base_element[phiM.solution_index];
        phiM.indices_exist = true;
    }
    else
        throw std::runtime_error("protonic_electrical_potential variable required for ORRCurrentDensityResponse");
    
    if ( this->system_management->solution_in_userlist("electronic_electrical_potential") )
    {
        phiS.solution_index = this->system_management->solution_name_to_index("electronic_electrical_potential");
        phiS.fetype_index = this->system_management->block_info->base_element[phiS.solution_index];
        phiS.indices_exist = true;
    }
    else
        throw std::runtime_error("electronic_electrical_potential variable required for ORRCurrentDensityResponse");
    
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::ORRCurrentDensityResponse<dim>::compute_responses(std::vector< FuelCellShop::SolutionVariable > solution_variables,
                                                        const typename DoFApplication<dim>::CellInfo& info,
                                                        FuelCellShop::Layer::BaseLayer<dim>* const layer, 
                                                        std::map<FuelCellShop::PostProcessing::ResponsesNames, double>& resp) const
{
    // Make sure you are in a CL       
    const std::type_info& base_layer = layer->get_base_type();
    const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
    AssertThrow(base_layer == CatalystLayer,
                ExcMessage("CurrentDensityResponse can only be used with a CatalystLayer object"));
    
    // Clear resp
    resp.clear();
    
    // Create a vector where the values of coefficients at quadrature points are stored:
    unsigned int n_q_points_cell = (info.fe(phiM.fetype_index)).n_quadrature_points;
    std::vector<double> values(n_q_points_cell);
    
    // Create CL:
    FuelCellShop::Layer::CatalystLayer<dim>* catalyst_layer = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);

    catalyst_layer->set_solution(solution_variables);
    catalyst_layer->set_cell_id(info.dof_active_cell->index());
    catalyst_layer->set_local_material_id(info.dof_active_cell->material_id());
    catalyst_layer->current_density(values);
    
    for (unsigned int q=0; q<n_q_points_cell; ++q)
    {
        double JxW = info.fe(phiM.fetype_index).JxW(q);
        //integrate RR = \nabla*i at the quadrature points. Since all components of f_values have
        //the value RR stored, I just use the first one.
        resp[FuelCellShop::PostProcessing::ResponsesNames::ORR_current] += -(1.0/S_CL)*(values[q] * JxW);
    }
    
    
    // NEXT: Check if kinetic model is DT and if so, output coverages
    if (catalyst_layer->get_kinetics_type() == "DoubleTrapKinetics")
    {
       
        //std::vector<double> O_coverage(n_q_points_cell, 0.0);
        //std::vector<double> OH_coverage(n_q_points_cell, 0.0);
        //catalyst_layer->get_kinetics()->O_coverage(O_coverage);
        //catalyst_layer->get_kinetics()->OH_coverage(OH_coverage);
        SolutionMap coverages = catalyst_layer->get_coverages();
        
        for (unsigned int q=0; q<n_q_points_cell; ++q)
        {
            double JxW = info.fe(phiM.fetype_index).JxW(q);
            double volume = info.dof_active_cell->measure();
            if(coverages.has(VariableNames::O_coverage))
                resp[FuelCellShop::PostProcessing::ResponsesNames::O_coverage] += (1/V_CL)*coverages.at(VariableNames::O_coverage)[q]*JxW;
            if(coverages.has(VariableNames::OH_coverage))
                resp[FuelCellShop::PostProcessing::ResponsesNames::OH_coverage] += (1/V_CL)*coverages.at(VariableNames::OH_coverage)[q]*JxW;
        }
    }
    else{
        resp[FuelCellShop::PostProcessing::ResponsesNames::O_coverage] = 0.0;
        resp[FuelCellShop::PostProcessing::ResponsesNames::OH_coverage] = 0.0;
    }
    
}
//---------------------------------------------------------------------------
template <int dim>
void 
NAME::ORRCurrentDensityResponse<dim>::compute_responses(const typename DoFApplication<dim>::CellInfo& info,
                                                        FuelCellShop::Layer::BaseLayer<dim>* const layer, 
                                                        std::map<FuelCellShop::PostProcessing::ResponsesNames, double>& resp) const
{ 
    // Clear resp
    resp.clear();
    
    // Make sure you are in a CL       
    const std::type_info& base_layer = layer->get_base_type();
    const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
    AssertThrow(base_layer == CatalystLayer,
                ExcMessage("CurrentDensityResponse can only be used with a CatalystLayer object"));
                   
    // Creating solution variable vector to passed to catalyst layer classes
    unsigned int solution_cell = info.global_data->find_vector("Solution");
    std::vector< FuelCellShop::SolutionVariable > solution_variables;
    if ( this->system_management->solution_in_userlist("oxygen_molar_fraction") )
        solution_variables.push_back( FuelCellShop::SolutionVariable(&info.values[solution_cell][xi.solution_index], oxygen_molar_fraction) );
    else //if ( this->system_management->solution_in_userlist("oxygen_molar_fraction") )
        AssertThrow( false , ExcMessage("Error in ORRCurrentDensityResponse<dim>::compute_responses. oxygen_molar_fraction needed in current implementation") );
    
    solution_variables.push_back( FuelCellShop::SolutionVariable(&info.values[solution_cell][phiM.solution_index], protonic_electrical_potential) );
    solution_variables.push_back( FuelCellShop::SolutionVariable(&info.values[solution_cell][phiS.solution_index], electronic_electrical_potential) );
    if ( this->system_management->solution_in_userlist("temperature_of_REV") )
        solution_variables.push_back( FuelCellShop::SolutionVariable(&info.values[solution_cell][this->system_management->solution_name_to_index("temperature_of_REV")], temperature_of_REV) );
    
    compute_responses(solution_variables, info, layer, resp);                                                        
    
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

template <int dim>
void 
NAME::HORCurrentDensityResponse<dim>::declare_parameters(ParameterHandler& param) const
{
    param.enter_subsection("PostProcessing"); 
    {
        param.enter_subsection("HORCurrentDensityResponse"); 
        {
            param.declare_entry ("Catalyst layer surface area [cm^2]",  // volume fraction of void space in the CL
                                 "1.0", // [-]
                                 Patterns::Double(),
                                 "Active area of the CL. Used to scale the volumetric current density");
            param.declare_entry ("Catalyst layer volume [cm^3]",  // volume fraction of void space in the CL
                                 "1.0", // [-]
                                 Patterns::Double(),
                                 "Volume of the CL. Used to average the coverages");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::HORCurrentDensityResponse<dim>::initialize(ParameterHandler&param )
{
    param.enter_subsection("PostProcessing"); 
    {
        param.enter_subsection("HORCurrentDensityResponse"); 
        {
            S_CL = param.get_double("Catalyst layer surface area [cm^2]");
            V_CL = param.get_double("Catalyst layer volume [cm^3]");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    phiM.solution_index = this->system_management->solution_name_to_index("protonic_electrical_potential");    
    phiS.solution_index = this->system_management->solution_name_to_index("electronic_electrical_potential");
    phiS.fetype_index = this->system_management->block_info->base_element[phiS.solution_index];
    phiS.indices_exist = true;
    
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::HORCurrentDensityResponse<dim>::compute_responses(std::vector< FuelCellShop::SolutionVariable > solution_variables,
                                                        const typename DoFApplication<dim>::CellInfo& info,
                                                        FuelCellShop::Layer::BaseLayer<dim>* const layer, 
                                                        std::map<FuelCellShop::PostProcessing::ResponsesNames, double>& resp) const
{
    // Make sure you are in a CL       
    const std::type_info& base_layer = layer->get_base_type();
    const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
    AssertThrow(base_layer == CatalystLayer,
                ExcMessage("CurrentDensityResponse can only be used with a CatalystLayer object"));
    
    // Clear resp
    resp.clear();
    
    // Create a vector where the values of coefficients at quadrature points are stored:
    unsigned int n_q_points_cell = (info.fe(phiS.fetype_index)).n_quadrature_points;
    std::vector<double> values(n_q_points_cell);
    
    // Create CL:
    FuelCellShop::Layer::CatalystLayer<dim>* catalyst_layer = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);

    catalyst_layer->set_solution(solution_variables);
    catalyst_layer->set_cell_id(info.dof_active_cell->index());
    catalyst_layer->current_density(values);

    for (unsigned int q=0; q<n_q_points_cell; ++q)
    {
        double JxW = info.fe(phiS.fetype_index).JxW(q);
        //integrate RR = \nabla*i at the quadrature points. Since all components of f_values have
        //the value RR stored, I just use the first one.
        resp[FuelCellShop::PostProcessing::ResponsesNames::HOR_current] += (1.0/S_CL)*(values[q] * JxW);
    }
}
//---------------------------------------------------------------------------
template <int dim>
void 
NAME::HORCurrentDensityResponse<dim>::compute_responses(const typename DoFApplication<dim>::CellInfo& info,
                                                        FuelCellShop::Layer::BaseLayer<dim>* const layer, 
                                                        std::map<FuelCellShop::PostProcessing::ResponsesNames, double>& resp) const
{ 
    // Clear resp
    resp.clear();
    
    // Make sure you are in a CL       
    const std::type_info& base_layer = layer->get_base_type();
    const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
    AssertThrow(base_layer == CatalystLayer,
                ExcMessage("CurrentDensityResponse can only be used with a CatalystLayer object"));
     
    // Creating solution variable vector to passed to catalyst layer classes:
    unsigned int solution_cell = info.global_data->find_vector("Solution");
    std::vector< FuelCellShop::SolutionVariable > solution_variables;

    unsigned int n_q_points_cell = (info.fe(phiS.fetype_index)).n_quadrature_points;
    std::vector< double > x_H2(n_q_points_cell);
    
    if ( this->system_management->solution_in_userlist("hydrogen_molar_fraction") )
    {
        x_H2 = info.values[solution_cell][this->system_management->solution_name_to_index("hydrogen_molar_fraction")];
    }
    else if ( this->system_management->solution_in_userlist("water_molar_fraction") )        
    {
        for(unsigned int q = 0; q < n_q_points_cell; ++q)
            x_H2[q] = 1.0 - info.values[solution_cell][this->system_management->solution_name_to_index("water_molar_fraction")][q];
    }
    else
        AssertThrow(false, ExcMessage("Possiblity not implemented either water or hydrogen mole fraction needed") );
    
    solution_variables.push_back( FuelCellShop::SolutionVariable(x_H2, hydrogen_molar_fraction) );
    solution_variables.push_back( FuelCellShop::SolutionVariable(&info.values[solution_cell][phiM.solution_index], protonic_electrical_potential) );
    solution_variables.push_back( FuelCellShop::SolutionVariable(&info.values[solution_cell][phiS.solution_index], electronic_electrical_potential) );
    if ( this->system_management->solution_in_userlist("temperature_of_REV") )
        solution_variables.push_back( FuelCellShop::SolutionVariable(&info.values[solution_cell][this->system_management->solution_name_to_index("temperature_of_REV")], temperature_of_REV) );
    
    compute_responses(solution_variables, info, layer, resp);                                                        
    
}



//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations.
template class NAME::ORRCurrentDensityResponse<deal_II_dimension>;
template class NAME::HORCurrentDensityResponse<deal_II_dimension>;
