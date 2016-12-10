//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2014 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: data_out.cc
//    - Description: 
//    - Developers: M. Secanell
//    - $Id: data_out.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include <postprocessing/data_out.h>

namespace NAME = FuelCellShop::PostProcessing;

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

template <int dim>
NAME::ORRCurrentDensityDataOut<dim>::ORRCurrentDensityDataOut(FuelCell::SystemManagement* sm,
                                                              boost::shared_ptr< FuelCellShop::Layer::CatalystLayer<dim> > cl,
                                                              FuelCell::OperatingConditions* operating)
{
    this->system_management = sm;
    catalyst_layer = cl;
    opCond = operating;
    oxygen_density_name = "";
}

//---------------------------------------------------------------------------
template <int dim>
std::vector<std::string>
NAME::ORRCurrentDensityDataOut<dim>::get_names() const
{
    std::vector<std::string> solution_names;
    solution_names.push_back ("ORR_volumetric_current_density");
    solution_names.push_back ("ORR_overpotential");
    solution_names.push_back ("ORR_effectiveness");
    solution_names.push_back ("O_coverage");
    solution_names.push_back ("OH_coverage");
    
    return solution_names;
}

//---------------------------------------------------------------------------
template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
NAME::ORRCurrentDensityDataOut<dim>::get_data_component_interpretation () const
{
    std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation;
    interpretation.push_back (DataComponentInterpretation::component_is_scalar);
    interpretation.push_back (DataComponentInterpretation::component_is_scalar);
    interpretation.push_back (DataComponentInterpretation::component_is_scalar);
    interpretation.push_back (DataComponentInterpretation::component_is_scalar);
    interpretation.push_back (DataComponentInterpretation::component_is_scalar);
    
    return interpretation;
}

//---------------------------------------------------------------------------
template <int dim>
UpdateFlags
NAME::ORRCurrentDensityDataOut<dim>::get_needed_update_flags() const
{
return update_values | update_gradients | update_q_points;
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::ORRCurrentDensityDataOut<dim>::compute_derived_quantities_vector (const std::vector< Vector< double > > &solution,
                                                                        const std::vector< std::vector< Tensor< 1, dim > > > & /*duh*/,
                                                                        const std::vector< std::vector< Tensor< 2, dim > > > & /*dduh*/,
                                                                        const std::vector< Point< dim > > & /*normals*/,
                                                                        const std::vector< Point<dim> > & /*evaluation_points*/,
                                                                        const types::material_id & mat_id,
                                                                        std::vector< Vector< double > > &computed_quantities) const
{
    Assert(this->system_management->solution_in_userlist ("protonic_electrical_potential"), ExcInvalidState());
    Assert(this->system_management->solution_in_userlist ("electronic_electrical_potential"), ExcInvalidState());
 
    unsigned int n_quad_points(solution.size());
    
    std::vector<double> current(n_quad_points, 0.0);
    std::vector<double> overpotential(n_quad_points, 0.0);
    std::vector<double> effectiveness(n_quad_points, 0.0);
    std::vector<double> O_coverage(n_quad_points, 0.0);
    std::vector<double> OH_coverage(n_quad_points, 0.0);

    if( catalyst_layer->belongs_to_material(mat_id) )
    {
        Assert( catalyst_layer->get_kinetics()->get_reaction_name() == ORR, ExcInvalidState() );
        
        std::vector<double> c_O2;
        std::vector<double> phi_m;
        std::vector<double> phi_s;
        std::vector<double> t_rev;

        // Obtain the necessary solution variables at each quadrature point:
        for (unsigned int i = 0; i<n_quad_points; ++i)
        {           
            // -- Get solution variables
            // oxygen concentration:
            if (this->system_management->solution_in_userlist ("oxygen_molar_fraction"))
            {
                c_O2.push_back(opCond->get_c_c() * solution[i](this->system_management->solution_name_to_index("oxygen_molar_fraction")) );            
            }
            else if (this->system_management->solution_in_userlist (this->oxygen_density_name) )
            {
                FuelCellShop::Material::Oxygen oxygen;
                double molar_mass = Units::convert(oxygen.get_molar_mass(), Units::UNIT, Units::K_UNIT);
                double density = solution[i](this->system_management->solution_name_to_index(this->oxygen_density_name));  //go from density to concentration
                double c_o2_i = density/molar_mass;
                c_O2.push_back( c_o2_i );   
            }         
            // potentials
            phi_m.push_back(solution[i](this->system_management->solution_name_to_index("protonic_electrical_potential")) );
            phi_s.push_back(solution[i](this->system_management->solution_name_to_index("electronic_electrical_potential")) );
            // temp
            double T(0.0);
            if ( this->system_management->solution_in_userlist("temperature_of_REV") )
            {
                T = solution[i](this->system_management->solution_name_to_index("temperature_of_REV"));
                t_rev.push_back(T);
            }
            else
                T = opCond->get_T();
               
            overpotential[i] = phi_s[i] - phi_m[i] - catalyst_layer->get_kinetics()->get_cat()->voltage_cell_th(T);
        }

        // -- Create solution variable vector to compute the current_density
        std::vector< FuelCellShop::SolutionVariable > solution_variables;
        solution_variables.push_back( FuelCellShop::SolutionVariable(&c_O2, VariableNames::oxygen_concentration) );
        solution_variables.push_back( FuelCellShop::SolutionVariable(&phi_m, protonic_electrical_potential) );
        solution_variables.push_back( FuelCellShop::SolutionVariable(&phi_s, electronic_electrical_potential ) );
                    
        if ( this->system_management->solution_in_userlist("temperature_of_REV") )
        {
            solution_variables.push_back( FuelCellShop::SolutionVariable(&t_rev, temperature_of_REV ) );
        }

        // -- Compute current:
        // It cannot happen in parallel, so lock object:
        catalyst_layer_mutex.lock();
        {        
            catalyst_layer->set_solution(solution_variables);
            catalyst_layer->current_density(current, effectiveness);
            
            if (catalyst_layer->get_kinetics_type() == "DoubleTrapKinetics")
            {
                SolutionMap coverages = catalyst_layer->get_coverages();
                
                if(coverages.has(VariableNames::O_coverage))
                    O_coverage = coverages.at(VariableNames::O_coverage).get_default_data();
                if(coverages.has(VariableNames::OH_coverage))
                    OH_coverage = coverages.at(VariableNames::OH_coverage).get_default_data();
            }
        }
        catalyst_layer_mutex.unlock();
    }

    for (unsigned int i = 0; i<n_quad_points; ++i)
    {
        computed_quantities[i](0) = current[i];
        computed_quantities[i](1) = overpotential[i];
        computed_quantities[i](2) = effectiveness[i];
        computed_quantities[i](3) = O_coverage[i];
        computed_quantities[i](4) = OH_coverage[i];
    }
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

template <int dim>
NAME::HORCurrentDensityDataOut<dim>::HORCurrentDensityDataOut(FuelCell::SystemManagement* sm,
                                                              boost::shared_ptr< FuelCellShop::Layer::CatalystLayer<dim> > cl,
                                                              FuelCell::OperatingConditions* operating)
{
    this->system_management = sm;
    catalyst_layer = cl;
    opCond = operating;
}

//---------------------------------------------------------------------------
template <int dim>
std::vector<std::string>
NAME::HORCurrentDensityDataOut<dim>::get_names() const
{
    std::vector<std::string> solution_names;
    solution_names.push_back ("HOR_volumetric_current_density");
    solution_names.push_back ("HOR_overpotential");
    solution_names.push_back ("HOR_effectiveness");
    return solution_names;
}

//---------------------------------------------------------------------------
template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
NAME::HORCurrentDensityDataOut<dim>::get_data_component_interpretation () const
{
    std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation;
    interpretation.push_back (DataComponentInterpretation::component_is_scalar);
    interpretation.push_back (DataComponentInterpretation::component_is_scalar);
    interpretation.push_back (DataComponentInterpretation::component_is_scalar);
    
    return interpretation;
}

//---------------------------------------------------------------------------
template <int dim>
UpdateFlags
NAME::HORCurrentDensityDataOut<dim>::get_needed_update_flags() const
{
return update_values | update_gradients | update_q_points;
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::HORCurrentDensityDataOut<dim>::compute_derived_quantities_vector (const std::vector< Vector< double > > &solution,
                                                                        const std::vector< std::vector< Tensor< 1, dim > > > & /*duh*/,
                                                                        const std::vector< std::vector< Tensor< 2, dim > > > & /*dduh*/,
                                                                        const std::vector< Point< dim > > & /*normals*/,
                                                                        const std::vector< Point<dim> > & /*evaluation_points*/,
                                                                        const types::material_id & mat_id,
                                                                        std::vector< Vector< double > > &computed_quantities) const
{
    Assert(this->system_management->solution_in_userlist ("protonic_electrical_potential"), ExcInvalidState());
    Assert(this->system_management->solution_in_userlist ("electronic_electrical_potential"), ExcInvalidState());
    
    unsigned int n_quad_points(solution.size());
    
    std::vector<double> current(n_quad_points,0.0);
    std::vector<double> overpotential(n_quad_points, 0.0);
    std::vector<double> effectiveness(n_quad_points,0.0);
    
    if(catalyst_layer->belongs_to_material(mat_id))
    {
        Assert( catalyst_layer->get_kinetics()->get_reaction_name() == HOR, ExcInvalidState() );
        
        std::vector<double> c_H2;
        std::vector<double> phi_m;
        std::vector<double> phi_s;
        std::vector<double> t_rev;

        for (unsigned int i = 0; i<n_quad_points; ++i)
        {
            if (this->system_management->solution_in_userlist ("hydrogen_molar_fraction"))
            {
                c_H2.push_back(opCond->get_c_a() * solution[i](this->system_management->solution_name_to_index("hydrogen_molar_fraction")) );            
            }
            else if (!this->system_management->solution_in_userlist ("hydrogen_molar_fraction") && this->system_management->solution_in_userlist ("water_molar_fraction"))
            {
                c_H2.push_back(opCond->get_c_a() * (1 - solution[i](this->system_management->solution_name_to_index("water_molar_fraction")) ) ); 
            }
            else if (this->system_management->solution_in_userlist (this->hydrogen_density_name) )
            {
                FuelCellShop::Material::Hydrogen hydrogen;
                double molar_mass = Units::convert(hydrogen.get_molar_mass(), Units::UNIT, Units::K_UNIT);
                double density = solution[i](this->system_management->solution_name_to_index(this->hydrogen_density_name));  //go from density to concentration
                double c_H2_i = density/molar_mass;
                c_H2.push_back( c_H2_i );   
            }              
            
            phi_m.push_back(solution[i](this->system_management->solution_name_to_index("protonic_electrical_potential")) );
            phi_s.push_back(solution[i](this->system_management->solution_name_to_index("electronic_electrical_potential")) );
            
            double T(0.0);
            if ( this->system_management->solution_in_userlist("temperature_of_REV") )
            {
                T = solution[i](this->system_management->solution_name_to_index("temperature_of_REV"));
                t_rev.push_back(T);
            }
            else
                T = opCond->get_T();
               
            overpotential[i] = phi_s[i] - phi_m[i] - catalyst_layer->get_kinetics()->get_cat()->voltage_cell_th(T);
        }

        std::vector< FuelCellShop::SolutionVariable > solution_variables;
        //solution_variables.push_back( FuelCellShop::SolutionVariable(&x_H2, hydrogen_molar_fraction)) ;
        solution_variables.push_back( FuelCellShop::SolutionVariable(&c_H2, hydrogen_concentration)) ;
        solution_variables.push_back( FuelCellShop::SolutionVariable(&phi_m, protonic_electrical_potential) );
        solution_variables.push_back( FuelCellShop::SolutionVariable(&phi_s, electronic_electrical_potential ) );
        if ( this->system_management->solution_in_userlist("temperature_of_REV") )
                solution_variables.push_back( FuelCellShop::SolutionVariable(&t_rev, temperature_of_REV ) );

        catalyst_layer_mutex.lock();
        catalyst_layer->set_solution(solution_variables);
        catalyst_layer->current_density(current, effectiveness);
        catalyst_layer_mutex.unlock();
    }

    for (unsigned int i = 0; i<n_quad_points; ++i)
    {
        computed_quantities[i](0) = current[i];
        computed_quantities[i](1) = overpotential[i];
        computed_quantities[i](2) = effectiveness[i];
    }
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

template <int dim>
NAME::RelativeHumidityDataOut<dim>::RelativeHumidityDataOut(FuelCell::SystemManagement* sm,
                                                            std::vector< boost::shared_ptr< FuelCellShop::Layer::PorousLayer<dim> > > pls,
                                                            FuelCell::OperatingConditions* operating)
:
DataPostprocessorScalar<dim> ("Relative_humidity", update_values | update_q_points )
{
    this->system_management = sm;
    porous_layers = pls;
    opCond = operating;
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::RelativeHumidityDataOut<dim>::compute_derived_quantities_vector (const std::vector< Vector< double > > &solution,
                                                                       const std::vector< std::vector< Tensor< 1, dim > > > & /*duh*/,
                                                                       const std::vector< std::vector< Tensor< 2, dim > > > & /*dduh*/,
                                                                       const std::vector< Point< dim > > & /*normals*/,
                                                                       const std::vector< Point<dim> > & /*evaluation_points*/,
                                                                       const types::material_id & mat_id,
                                                                       std::vector< Vector< double > > &computed_quantities) const
{    
    std::vector< double > RH( computed_quantities.size(), 0.0);
    
    for (unsigned int pl = 0; pl < porous_layers.size(); ++pl)
    {
        if( porous_layers.at(pl)->belongs_to_material(mat_id) )
        {
            // Iterate over every point:
            for (unsigned int i=0; i<computed_quantities.size(); i++)
            {
                Assert(computed_quantities[i].size() == 1, ExcDimensionMismatch (computed_quantities[i].size(), 1));
                
                double x_w(0.0), p_t(0.0), saturation_pressure(0.0);
                
                if (this->system_management->solution_in_userlist ("water_molar_fraction"))
                    x_w = solution[i](this->system_management->solution_name_to_index("water_molar_fraction"));
                else
                    x_w = opCond->get_x_wv();
                    
                if (this->system_management->solution_in_userlist ("temperature_of_REV"))
                {
                    double T = solution[i](this->system_management->solution_name_to_index("temperature_of_REV"));
                    saturation_pressure = water.get_water_vapor_saturation_pressure(T);
                }
                else
                    saturation_pressure = opCond->saturation_pressure();
                    
                porous_layers.at(pl)->get_p(p_t);
                
                RH[i] = ( (p_t * Units::convert(1.,Units::ATM_to_PA) * x_w) / saturation_pressure ) ;
            }
            
            goto OUTPUT; // Computation for this cell is done.
        }
    }

    OUTPUT:
    for (unsigned int i = 0; i<computed_quantities.size(); ++i)
        computed_quantities[i](0) = RH[i];
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


template <int dim>
NAME::CapillaryPressureDataOut<dim>::CapillaryPressureDataOut(FuelCell::SystemManagement* sm,
                                                              std::vector< boost::shared_ptr< FuelCellShop::Layer::PorousLayer<dim> > > pls,
                                                              FuelCell::OperatingConditions* operating)
:
DataPostprocessorScalar<dim> ("capillary_pressure", update_values | update_q_points )
{
    this->system_management = sm;
    opCond = operating;
    porous_layers = pls;
    
}

//---------------------------------------------------------------------------
template <int dim>
UpdateFlags
NAME::CapillaryPressureDataOut<dim>::get_needed_update_flags() const
{
return update_values;
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::CapillaryPressureDataOut<dim>::compute_derived_quantities_vector (const std::vector< Vector< double > > &solution,
                                                                       const std::vector< std::vector< Tensor< 1, dim > > > & /*duh*/,
                                                                       const std::vector< std::vector< Tensor< 2, dim > > > & /*dduh*/,
                                                                       const std::vector< Point< dim > > & /*normals*/,
                                                                       const std::vector< Point<dim> > & /*evaluation_points*/,
                                                                       const types::material_id & mat_id,
                                                                       std::vector< Vector< double > > &computed_quantities) const
{   
    Assert(this->system_management->solution_in_userlist ("liquid_water_saturation"), ExcInvalidState());
    
    unsigned int n_quad_points(solution.size());
    
    std::vector<double> capillary_pressure(n_quad_points,0.0);  
    
    std::vector<double> l_saturation;
    std::vector<double> t_rev;
    
    for (unsigned int pl = 0; pl < porous_layers.size(); ++pl)
    {
        if( porous_layers.at(pl)->belongs_to_material(mat_id) )
        {
            // Iterate over every point:
            for (unsigned int i=0; i<computed_quantities.size(); i++)
            {
                Assert(computed_quantities[i].size() == 1, ExcDimensionMismatch (computed_quantities[i].size(), 1));
                
                
                l_saturation.push_back(solution[i](this->system_management->solution_name_to_index("liquid_water_saturation")) );
                t_rev.push_back(solution[i](this->system_management->solution_name_to_index("temperature_of_REV")) );
                
                porous_layers.at(pl)->set_saturation(FuelCellShop::SolutionVariable(&l_saturation, liquid_water_saturation));
                porous_layers.at(pl)->set_temperature(FuelCellShop::SolutionVariable(&t_rev, temperature_of_REV));
                porous_layers.at(pl)->pcapillary(capillary_pressure);
            }
            
            goto OUTPUT; // Computation for this cell is done.
        }
    }
    
    OUTPUT:
    for (unsigned int i = 0; i<computed_quantities.size(); ++i)
        computed_quantities[i](0) = capillary_pressure[i];
    
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


template <int dim>
NAME::SaturationDataOut<dim>::SaturationDataOut(FuelCell::SystemManagement* sm,
                                                              std::vector< boost::shared_ptr< FuelCellShop::Layer::PorousLayer<dim> > > pls,
                                                              FuelCell::OperatingConditions* operating)
:
DataPostprocessorScalar<dim> ("Saturation", update_values | update_q_points )
{
    this->system_management = sm;
    opCond = operating;
    porous_layers = pls;
    
}

//---------------------------------------------------------------------------
template <int dim>
std::vector<std::string>
NAME::SaturationDataOut<dim>::get_names() const
{
    std::vector<std::string> solution_names;
    solution_names.push_back ("velocity_vector_x");
    solution_names.push_back ("velocity_vector_y");
    solution_names.push_back ("Saturation");
    solution_names.push_back ("Derivative_Saturation");
    solution_names.push_back ("Derivative_Permeability");
    solution_names.push_back ("Derivative_Interfacial_surface_area");
    
    return solution_names;
}

//---------------------------------------------------------------------------
template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
NAME::SaturationDataOut<dim>::get_data_component_interpretation () const
{
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation (dim,
                    DataComponentInterpretation::component_is_part_of_vector);
    interpretation.push_back (DataComponentInterpretation::component_is_scalar);
    interpretation.push_back (DataComponentInterpretation::component_is_scalar);
    interpretation.push_back (DataComponentInterpretation::component_is_scalar);
    interpretation.push_back (DataComponentInterpretation::component_is_scalar);
    
    return interpretation;
}

//---------------------------------------------------------------------------
template <int dim>
UpdateFlags
NAME::SaturationDataOut<dim>::get_needed_update_flags() const
{
    return update_values | update_gradients | update_q_points;
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::SaturationDataOut<dim>::compute_derived_quantities_vector (const std::vector< Vector< double > > &solution,
                                                                       const std::vector< std::vector< Tensor< 1, dim > > > & duh,
                                                                       const std::vector< std::vector< Tensor< 2, dim > > > & /*dduh*/,
                                                                       const std::vector< Point< dim > > & /*normals*/,
                                                                       const std::vector< Point<dim> > & /*evaluation_points*/,
                                                                       const types::material_id & mat_id,
                                                                       std::vector< Vector< double > > &computed_quantities) const
{   
    Assert(this->system_management->solution_in_userlist ("capillary_pressure"), ExcInvalidState());
    
    unsigned int n_quad_points(solution.size());
    
    std::vector<double> saturation(n_quad_points,0.0);
    
    std::vector<double> derivative_saturation(n_quad_points,0.0);
    
    std::vector<double> derivative_Permeability(n_quad_points,0.0);
    
    std::vector<double> derivative_interfafical_surface_area_PSD(n_quad_points,0.0);
    
    std::vector<double> l_capillary;
    
    for (unsigned int pl = 0; pl < porous_layers.size(); ++pl)
    {
        if( porous_layers.at(pl)->belongs_to_material(mat_id) )
        {
            
            Assert(computed_quantities.size() == solution.size(), ExcDimensionMismatch (computed_quantities.size(), solution.size()));
            
            // Iterate over every point:
            for (unsigned int i=0; i<computed_quantities.size(); i++)
            {
                l_capillary.push_back(solution[i](this->system_management->solution_name_to_index("capillary_pressure")) );
                
                porous_layers.at(pl)->set_capillary_pressure(FuelCellShop::SolutionVariable(&l_capillary, capillary_pressure));
                
                porous_layers.at(pl)->saturation_from_capillary_equation(saturation);
                
                porous_layers.at(pl)->derivative_saturation_from_capillary_equation_PSD(derivative_saturation);
                
                porous_layers.at(pl)->derivative_relative_liquid_permeablity_PSD(derivative_Permeability);
                
                porous_layers.at(pl)->derivative_interfacial_surface_area_PSD(derivative_interfafical_surface_area_PSD);

            }
            

            
            goto OUTPUT; // Computation for this cell is done.
        }
    }
    
    OUTPUT:
    for (unsigned int i = 0; i<n_quad_points; ++i)
    {
        for (unsigned int d=0; d<dim; ++d)
            computed_quantities[i](d) = duh[i][this->system_management->solution_name_to_index("capillary_pressure")][d] ;
        
        computed_quantities[i](dim) = saturation[i];
        
        computed_quantities[i](dim+1) = derivative_saturation[i];
        
        computed_quantities[i](dim+2) = derivative_Permeability[i];
        
        computed_quantities[i](dim+3) = derivative_interfafical_surface_area_PSD[i];
    }
    
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

// Explicit instantiations.
template class NAME::ORRCurrentDensityDataOut<deal_II_dimension>;
template class NAME::HORCurrentDensityDataOut<deal_II_dimension>;
template class NAME::RelativeHumidityDataOut<deal_II_dimension>;
template class NAME::CapillaryPressureDataOut<deal_II_dimension>;
template class NAME::SaturationDataOut<deal_II_dimension>;
