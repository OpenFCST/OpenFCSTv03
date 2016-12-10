//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2015 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: data_out_mass_transport.h
//    - Description: Data postprocessing routine to output total density, concentration and molar fractions
//    - Developers: M. Secanell
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__POSTPROCESSING__DATAOUT_MASSTRANSPORT__H
#define _FUELCELLSHOP__POSTPROCESSING__DATAOUT_MASSTRANSPORT__H

// Include deal.II classes
#include <deal.II/numerics/data_out.h>

//Include STL
#include <boost/shared_ptr.hpp>
#include <mutex>
#include <typeinfo>

// Include OpenFCST routines:
#include <layers/porous_layer.h>
#include <materials/PureGas.h>
#include <utils/operating_conditions.h>

using namespace dealii;

namespace FuelCellShop
{
    namespace PostProcessing
    {
        /**
         * Class used to output molar fractions based on the density of each species.
         * 
         * This class is necessary in mass transport applications.
         */
        template <int dim>
        class MolarFractionDataOut 
        : 
        public dealii::DataPostprocessorScalar<dim>
        {
        public:
            
            ///@name Constructors, destructor, and initalization
            //@{
            /**
             * Constructor of our class
             */
            MolarFractionDataOut(const FuelCell::SystemManagement& sm,
                                 FuelCellShop::Material::GasMixture& gas_in);
            /**
             * Destructor
             */
            virtual ~MolarFractionDataOut();
            
            /**
             * Function that provides the names of the output variables.
             * In this case, current density, overpotential, effectiveness, and oxygen coverages (if applicable). The latter
             * is set to one if the model is macro-homogeneous.
             */
            virtual std::vector<std::string> get_names() const;
            
            /**
             * Function that states if the output functions are a scalar or a vector.
             */
            virtual std::vector<DataComponentInterpretation::DataComponentInterpretation> get_data_component_interpretation () const;
            
            /**
             * Flags to be updated in each cell when computing the solution.
             */
            virtual UpdateFlags get_needed_update_flags() const;
            
            /**
             * 
             */                                                            
            virtual void compute_derived_quantities_vector (const std::vector< Vector< double > > &solution,
                                                            const std::vector< std::vector< Tensor< 1, dim > > > & /*duh*/,
                                                            const std::vector< std::vector< Tensor< 2, dim > > > & /*dduh*/,
                                                            const std::vector< Point< dim > > & /*normals*/,
                                                            const std::vector< Point<dim> > & /*evaluation_points*/,
                                                            const types::material_id & mat_id,
                                                            std::vector< Vector< double > > &computed_quantities) const;
            
        private:
            /**
             * Pointer to system management
             */
            const FuelCell::SystemManagement* system_management;
            
            /**
             * 
             */
            FuelCellShop::Material::GasMixture* fluid;
        };
    }
}

#endif