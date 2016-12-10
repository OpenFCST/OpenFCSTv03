// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: response_mass_flux.cc
// - Description: It contains definitions of mass flux response evaluator classes.
// - Developers: Chad Balen
//
// ----------------------------------------------------------------------------

#include <postprocessing/response_mass_flux.h>

namespace NAME = FuelCellShop::PostProcessing;

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::MassFluxResponse<dim>::compute_responses(const typename DoFApplication<dim>::CellInfo& info,
                                               FuelCellShop::Layer::BaseLayer<dim>* const layer, 
                                               std::map<FuelCellShop::PostProcessing::ResponsesNames, double>& respMap) const
{
    FcstUtilities::log << "MassFluxResponse<dim>::compute_responses() has not been implemented yet." << std::endl;
    Assert(false, ExcInternalError());
}

template <int dim>
void 
NAME::MassFluxResponse<dim>::compute_responses(std::vector< FuelCellShop::SolutionVariable > solution_variables,
                                               const typename DoFApplication<dim>::CellInfo& info,
                                               FuelCellShop::Layer::BaseLayer<dim>* const layer, 
                                               std::map<FuelCellShop::PostProcessing::ResponsesNames, double>& respMap) const
{
    FcstUtilities::log << "MassFluxResponse<dim>::compute_responses() has not been implemented yet." << std::endl;
    Assert(false, ExcInternalError());
}

template <int dim>
void 
NAME::MassFluxResponse<dim>::bdry_responses(std::vector<double>&                                                     dst,
                                            const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                            std::vector<std::string>&                                                responseNames,
                                            std::vector<unsigned int>&                                               bdryIDs) const
{
    const unsigned int current_bdry_id = bdry_info.dof_face->boundary_indicator(); // Get boundary ID of current location
    
    for(unsigned int i = 0; i < bdryIDs.size(); ++i) //Loop through vector of boundaries to apply response and see if mass flux needs to be calculated here
    {
        if (current_bdry_id == bdryIDs[i]) //If along boundary user wishes to apply bdry response calc mass flux
        {
            unsigned int n_q_points_bdry = bdry_info.get_fe_val_unsplit().n_quadrature_points; //Set number of quadrature points
            std::vector< Point<dim> > normal_vectors = bdry_info.get_fe_val_unsplit().get_normal_vectors();
            
            for(unsigned int j = 0; j < responseNames.size(); ++j) //Loop through responseNames and see which indecies are for mass flux
            {
                /**Looking for string of form: total_mass_flux_species_* at boundary *
                 * To do this there is two checks:
                 * 1) Is the response name of type total_mass_flux_species_*
                 * 2) get last position of string and confirm that it's boundary is same as current boundary
                 * NOTE: atoi will return 0 if string cannot be converted to an integer, but will be fine because first part would fail then.
                 */
                if ( (responseNames[j].compare(0,24,"total_mass_flux_species_") == 0) && (atoi(responseNames[j].substr(38).c_str()) == current_bdry_id) )
                {
                    std::string species     = responseNames[j].substr(24); //get substring from position 24 to end which should be the full species number
                    unsigned int speciesKey = atoi(species.c_str()) - 1; //convert species string number to int to be used to index position of velocity and density components
                    
                    unsigned int density_index        = speciesKey*(dim+1) + this->system_management->solution_name_to_index("density_species_1"); //Use speciesKey to get correct density and velocity
                    unsigned int velocity_first_index = density_index + 1;
                    
                    FEValuesExtractors::Scalar densityExtractor  = FEValuesExtractors::Scalar(density_index);
                    FEValuesExtractors::Vector velocityExtractor = FEValuesExtractors::Vector(velocity_first_index);
                    
                    std::vector<double> density_bdry(n_q_points_bdry);
                    std::vector< Tensor<1,dim> > velocity_bdry(n_q_points_bdry);
                    
                    bdry_info.get_fe_val_unsplit()[ densityExtractor ].get_function_values(bdry_info.global_data->vector(bdry_info.global_data->n_vectors()-1 ), density_bdry);
                    bdry_info.get_fe_val_unsplit()[ velocityExtractor ].get_function_values(bdry_info.global_data->vector(bdry_info.global_data->n_vectors()-1 ), velocity_bdry);
                    
                    
                    for(unsigned int q = 0; q < n_q_points_bdry; ++q)
                    {
                        dst[j]+= velocity_bdry[q] * normal_vectors[q]
                                 *
                                 density_bdry[q]
                                 *
                                 bdry_info.get_fe_val_unsplit().JxW(q);

                    }
                }
            }
        }
    }
}

/////////////////////////////
// EXPLICIT INSTANTIATIONS //
/////////////////////////////

// --- MassFluxResponse ---
template class NAME::MassFluxResponse<deal_II_dimension>;