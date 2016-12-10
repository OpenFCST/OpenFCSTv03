//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2015 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: data_out_mass_transport.cc
//    - Description: Data postprocessing routine to output total density, concentration and molar fractions
//    - Developers: M. Secanell
//
//---------------------------------------------------------------------------

#include <postprocessing/data_out_mass_transport.h>

namespace NAME = FuelCellShop::PostProcessing;

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template <int dim>
NAME::MolarFractionDataOut<dim>::MolarFractionDataOut(const FuelCell::SystemManagement& sm,
                                                      FuelCellShop::Material::GasMixture& gas_in)
:
DataPostprocessorScalar<dim> ("Total concentration", update_values),
system_management(&sm),
fluid(&gas_in)
{ }

//---------------------------------------------------------------------------
template <int dim>
NAME::MolarFractionDataOut<dim>::~MolarFractionDataOut() 
{}

//---------------------------------------------------------------------------
template <int dim>
std::vector<std::string>
NAME::MolarFractionDataOut<dim>::get_names() const
{    
    std::vector<std::string> solution_names;
    solution_names.push_back ("Total_density");
    solution_names.push_back ("Total_concentration");
    solution_names.push_back ("Total_pressure");
    
    for (unsigned int i=0; i < fluid->get_gases().size(); i++)
    {
        std::string name = fluid->get_gases()[i]->name_material()+"_molar_fraction";
        solution_names.push_back (name);
    }
    
    return solution_names;                
}


//---------------------------------------------------------------------------
template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
NAME::MolarFractionDataOut<dim>::get_data_component_interpretation () const
{
    std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation;
    interpretation.push_back (DataComponentInterpretation::component_is_scalar);
    interpretation.push_back (DataComponentInterpretation::component_is_scalar);
    interpretation.push_back (DataComponentInterpretation::component_is_scalar);

    for (unsigned int i=0; i < fluid->get_gases().size(); i++)
    {
        interpretation.push_back (DataComponentInterpretation::component_is_scalar);
    }
    
    return interpretation;
}

//---------------------------------------------------------------------------
template <int dim>
UpdateFlags
NAME::MolarFractionDataOut<dim>::get_needed_update_flags() const
{
    return update_values;
}


//---------------------------------------------------------------------------
template <int dim>
void
NAME::MolarFractionDataOut<dim>::compute_derived_quantities_vector (const std::vector< Vector< double > > &solution,
                                                                    const std::vector< std::vector< Tensor< 1, dim > > > & /*duh*/,
                                                                    const std::vector< std::vector< Tensor< 2, dim > > > & /*dduh*/,
                                                                    const std::vector< Point< dim > > & /*normals*/,
                                                                    const std::vector< Point<dim> > & /*evaluation_points*/,
                                                                    const types::material_id & mat_id,
                                                                    std::vector< Vector< double > > &computed_quantities) const
{
    // Check that you are solving a problem with fluid flow, i.e. densities should exist:
    for (unsigned int s=0; s < fluid->get_gases().size(); s++)
    {
        std::string name = "density_species_"+std::to_string(s+1);
        AssertThrow(system_management->solution_in_userlist (name), 
                    ExcMessage("MolarFractionDataOut should be used with applications that solve for species densities"));
    }
    
    unsigned int n_quad_points(solution.size());
    
    for (unsigned int q = 0; q<n_quad_points; ++q)
    {
        unsigned int n_species = fluid->get_gases().size();
        
        std::vector<double> density;
        density.resize(n_species);
        
        std::vector<double> concentration;
        concentration.resize(n_species);
        
        double density_t = 0.0;
        double concentration_t = 0.0;
        
        for (unsigned int s=0; s < n_species; s++)
        {
            std::string name = "density_species_"+std::to_string(s+1);
            density[s] = solution[q](system_management->solution_name_to_index(name));
            density_t += density[s];
            concentration[s] = density[s]/(fluid->get_gases()[s]->get_molar_mass()*1e3); //Note: 1e3 to go from kg to g as density is reported in g/cm3
            concentration_t += concentration[s];                        
        }
        
        // Total_density:
        computed_quantities[q](0) = density_t;
        // Total_concentration
        computed_quantities[q](1) = concentration_t;    
        // Total_concentration
        computed_quantities[q](2) = concentration_t*Constants::R()*fluid->get_temperature()*Units::convert(1.,Units::C_UNIT3, Units::UNIT3);    
        
        // Species molar fractions:
        for (unsigned int s=0; s < n_species; s++)
        {
            computed_quantities[q](s+3) = concentration[s]/concentration_t;
        }
    }
}    

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

// Explicit instantiations.
template class NAME::MolarFractionDataOut<deal_II_dimension>;