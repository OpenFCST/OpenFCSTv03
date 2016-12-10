//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2014 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: data_out.h
//    - Description: 
//    - Developers: M. Secanell
//    - $Id: data_out.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__POSTPROCESSING__DATAOUT_H
#define _FUELCELLSHOP__POSTPROCESSING__DATAOUT_H

// Include deal.II classes
#include <deal.II/numerics/data_out.h>

//Include STL
#include <boost/shared_ptr.hpp>
#include <mutex>
#include <typeinfo>

// Include OpenFCST routines:

#include <layers/porous_layer.h>
#include <layers/catalyst_layer.h>
#include <layers/gas_diffusion_layer.h>
#include <layers/micro_porous_layer.h>

#include <materials/PureGas.h>

#include <utils/operating_conditions.h>

using namespace dealii;

namespace FuelCellShop
{
    /**
     * Namespace used for all classes use for post-processing either for evaluating
     * a new quantity at a quadrature point for outputting with the solution,
     * for the evaluation of functions or for the evaluation of functionals
     * such as the current density.
     * 
     */
    namespace PostProcessing
    {
        
        /**
         * Class used to evaluate the ORR current density, overpotential, effectiveness and the oxygen coverages (when applicable)
         * at catalyst layer DoF point in the finite element mesh.
         * 
         * Child of DataPostprocessor (deal.II post-processing class used in DataOut class).
         * 
         * See step-29, step-32 and step-33 in the deal.II website for more details about this class.
         * 
         * 
         * <h3>Usage</h3>
         * 
         * This class is used in the data_out section of any application. It requires a pointer to SystemManagement, a
         * pointer to the catalyst layer object for which you would like to compute the current density, and a pointer to OperatingConditions object.
         *
         * @code
         * AppCathode<dim>::data_out(const std::string& filename,
         *                           const FuelCell::ApplicationCore::FEVectors& src)
         * 
         * // --- Find solution and assign solution interpretations ---
         * (...)
         * 
         * // --- Create vector of PostProcessing objects ---
         * std::vector< DataPostprocessor<dim>* > PostProcessing;
         * 
         * // --- cathode_current --- For instance, on CCL object.
         * // Create the post processing object and push it to the vector.
         * FuelCellShop::PostProcessing::ORRCurrentDensityDataOut<dim> cathode_current(&system_management, CCL, &OC);
         * PostProcessing.push_back(&cathode_current);
         *
         * // Outputting
         * DoFApplication<dim>::data_out(filename,
         *                               solution, // This is initialized in the initial part of the data_out function
         *                               system_management.get_solution_names(),
         *                               PostProcessing);
         * @endcode
         * 
         * \author M. Secanell, 2014
         */
        template <int dim>
        class ORRCurrentDensityDataOut 
        : 
        public dealii::DataPostprocessor<dim>
        {
        public:
            
            ///@name Constructors, destructor, and initalization
            //@{
            /**
             * Constructor of our class
             */
            ORRCurrentDensityDataOut(FuelCell::SystemManagement* ,
                                     boost::shared_ptr< FuelCellShop::Layer::CatalystLayer<dim> > ,
                                     FuelCell::OperatingConditions* );
            /**
             * Destructor
             */
            virtual ~ORRCurrentDensityDataOut() {}
            
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
             * For the case of multi-component solvers, we need to specify which one of the variables contains the oxygen concentration,
             * e.g., density_species_1, density_species_2 and so on. Provide this value here.
             */
            void set_oxygen_density_name(std::string name)
            {
                oxygen_density_name = name;
            }
            
            /**
             * Member function used to calculate the current density. 
             * Its inputs are 
             * - a vector representing values of the function (which is here vector-valued) representing the data 
             * vector given to DataOut::add_data_vector, evaluated at all evaluation points where we generate output
             * - a tensor objects representing the first derivative (not used and not given properly since in the constructor
             * we specified that only solution values need to be updated)
             * 
             * The derived quantities are returned in the computed_quantities vector. 
             * 
             * Remember that this function may only use data for which the respective update flag is specified by either the
             * constructor or get_needed_update_flags.
             */
            virtual void compute_derived_quantities_vector (const std::vector< Vector< double > > &uh,
                                                    const std::vector< std::vector< Tensor< 1, dim > > > & /*duh*/,
                                                    const std::vector< std::vector< Tensor< 2, dim > > > & /*dduh*/,
                                                    const std::vector< Point< dim > > & /*normals*/,
                                                    const std::vector< Point<dim> > & /*evaluation_points*/,
                                                    const types::material_id & mat_id,
                                                    std::vector< Vector< double > > &computed_quantities) const;

                                                            
        private:
            /**
             * Pointer to catalyst layer object
             */
            boost::shared_ptr< FuelCellShop::Layer::CatalystLayer<dim> > catalyst_layer;
            
            /**
             * Pointer to system management
             */
            FuelCell::SystemManagement* system_management;

            /**
             * Pointer to operating conditions class
             */
            FuelCell::OperatingConditions* opCond;
            
            /**
             * The mutex class used to prevent simultaneous calls to the layer objects from multiple threads. 
             * 
             * @note The mutable keywork is used since this data member needs to be modified in a const function
             * (the function needs to be const because of deal.ii)
             */
            mutable std::mutex catalyst_layer_mutex;
            
            /**
             * Set oxygen concentration name
             */
            std::string oxygen_density_name;
        };
        
        //-------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------
        
        /**
         * Class used to evaluate the HOR current density, overpotential and effectiveness
         * at catalyst layer DoF point in the finite element mesh.
         * 
         * Child of DataPostprocessor (deal.II post-processing class used in DataOut class).
         * 
         * See step-29, step-32 and step-33 in the deal.II website for more details about this class.
         * 
         * 
         * <h3>Usage</h3>
         * 
         * This class is used in the data_out section of any application. It requires a pointer to SystemManagement, a
         * pointer to the catalyst layer object for which you would like to compute
         * the current density, and a pointer to OperatingConditions object.
         *
         * @code
         * AppCathode<dim>::data_out(const std::string& filename,
         *                           const FuelCell::ApplicationCore::FEVectors& src)
         * 
         * // --- Find solution and assign solution interpretations ---
         * (...)
         * 
         * // --- Create vector of PostProcessing objects ---
         * std::vector< DataPostprocessor<dim>* > PostProcessing;
         * 
         * // --- anode_current --- For instance, on ACL object.
         * // Create the post processing object and push it to the vector.
         * FuelCellShop::PostProcessing::HORCurrentDensityDataOut<dim> anode_current(&system_management, ACL, &OC);
         * PostProcessing.push_back(&anode_current);
         *
         * // Outputting
         * DoFApplication<dim>::data_out(filename,
         *                               solution, // This is initialized in the initial part of the data_out function
         *                               system_management.get_solution_names(),
         *                               PostProcessing);
         * @endcode
         * 
         * \author M. Secanell, 2014
         */
        template <int dim>
        class HORCurrentDensityDataOut 
        : 
        public dealii::DataPostprocessor<dim>
        {
        public:
            
            ///@name Constructors, destructor, and initalization
            //@{
            /**
             * Constructor of our class
             */
            HORCurrentDensityDataOut(FuelCell::SystemManagement* ,
                                     boost::shared_ptr< FuelCellShop::Layer::CatalystLayer<dim> > ,
                                     FuelCell::OperatingConditions* );
            /**
             * Destructor
             */
            virtual ~HORCurrentDensityDataOut() {}
            
            /**
             * Function that provides the names of the output variables.
             * In this case, current density, overpotential and effectiveness. The latter is set to 
             * one if the model is macro-homogeneous.
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
             * For the case of multi-component solvers, we need to specify which one of the variables contains the oxygen concentration,
             * e.g., density_species_1, density_species_2 and so on. Provide this value here.
             */
            void set_hydrogen_density_name(std::string name)
            {
                hydrogen_density_name = name;
            }
            /**
             * Member function used to calculate the current density. 
             * Its inputs are 
             * - a vector representing values of the function (which is here vector-valued) representing the data 
             * vector given to DataOut::add_data_vector, evaluated at all evaluation points where we generate output
             * - a tensor objects representing the first derivative (not used and not given properly since in the constructor
             * we specified that only solution values need to be updated)
             * 
             * The derived quantities are returned in the computed_quantities vector. 
             * 
             * Remember that this function may only use data for which the respective update flag is specified by either the
             * constructor or get_needed_update_flags.
             */
            virtual void compute_derived_quantities_vector (const std::vector< Vector< double > > &uh,
                                                    const std::vector< std::vector< Tensor< 1, dim > > > & /*duh*/,
                                                    const std::vector< std::vector< Tensor< 2, dim > > > & /*dduh*/,
                                                    const std::vector< Point< dim > > & /*normals*/,
                                                    const std::vector< Point<dim> > & /*evaluation_points*/,
                                                    const types::material_id & mat_id,
                                                    std::vector< Vector< double > > &computed_quantities) const;

                                                            
        private:
            /**
             * Pointer to catalyst layer object
             */
            boost::shared_ptr< FuelCellShop::Layer::CatalystLayer<dim> > catalyst_layer;
            
            /**
             * Pointer to system management
             */
            FuelCell::SystemManagement* system_management;

            /**
             * Pointer to operating conditions class
             */
            FuelCell::OperatingConditions* opCond;
            
            /**
             * The mutex class used to prevent simultaneous calls to the layer objects from multiple threads. 
             * 
             * @note The mutable keywork is used since this data member needs to be modified in a const function
             * (the function needs to be const because of deal.ii)
             */
            mutable std::mutex catalyst_layer_mutex;
            
            /**
             * Set hydrogen density name
             */
            std::string hydrogen_density_name;
        };
       
        //-------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------
        
       /**
         * Class used to evaluate the relative humidity at PorousLayer (\em viz., GDL, MPL and CL) DoF point in the
         * finite element mesh.
         * 
         * Child of DataPostprocessorScalar (deal.II post-processing class used in DataOut class).
         * 
         * See step-29, step-32 and step-33 in the deal.II website for more details about this class.
         * 
         * 
         * <h3>Usage</h3>
         * 
         * This class is used in the data_out section of any application. It requires a pointer to SystemManagement 
         * object, a vector of pointer to the PorousLayer objects for which you would like to compute
         * the current density, and pointer to OperatingConditions object.
         *
         * @code
         * AppCathode<dim>::data_out(const std::string& filename,
         *                           const FuelCell::ApplicationCore::FEVectors& src)
         * 
         * // --- Find solution and assign solution interpretations ---
         * (...)
         * 
         * // --- Create vector of PostProcessing objects ---
         * std::vector< DataPostprocessor<dim>* > PostProcessing;
         * 
         * // --- relative_humidity ---
         * // First create a vector of PorousLayer objects. We have CGDL, CMPL and CCL in the cathode application.
         * std::vector< boost::shared_ptr< FuelCellShop::Layer::PorousLayer<dim> > > porous_layers;
         * porous_layers.push_back(CGDL);
         * porous_layers.push_back(CMPL);
         * porous_layers.push_back(CCL);
         *
         * // Then create the post processing object and push it to the vector.
         * FuelCellShop::PostProcessing::RelativeHumidityDataOut<dim> relative_humidity(&system_management, porous_layers, &OC);
         * PostProcessing.push_back(&relative_humidity);
         *
         * // Outputting
         * DoFApplication<dim>::data_out(filename,
         *                               solution, // This is initialized in the initial part of the data_out function
         *                               system_management.get_solution_names(),
         *                               PostProcessing);
         * @endcode
         * 
         * \author M. Secanell, 2014
         */
        template <int dim>
        class RelativeHumidityDataOut 
        : 
        public dealii::DataPostprocessorScalar<dim>
        {
        public:
            
            ///@name Constructors, destructor, and initalization
            //@{
            /**
             * Constructor of our class
             */
            RelativeHumidityDataOut(FuelCell::SystemManagement* ,
                                    std::vector< boost::shared_ptr< FuelCellShop::Layer::PorousLayer<dim> > > ,
                                    FuelCell::OperatingConditions* );
            /**
             * Destructor
             */
            virtual ~RelativeHumidityDataOut() {}
           
                                                            
            virtual void compute_derived_quantities_vector (const std::vector< Vector< double > > &uh,
                                                            const std::vector< std::vector< Tensor< 1, dim > > > & /*duh*/,
                                                            const std::vector< std::vector< Tensor< 2, dim > > > & /*dduh*/,
                                                            const std::vector< Point< dim > > & /*normals*/,
                                                            const std::vector< Point<dim> > & /*evaluation_points*/,
                                                            const types::material_id & mat_id,
                                                            std::vector< Vector< double > > &computed_quantities) const;
            
                                                            
        private:
            
            /**
             * WaterVapor object, used to compute saturation pressure as a function of temperature.
             */
            FuelCellShop::Material::WaterVapor water;
            
            /**
             * Pointer to system management
             */
            FuelCell::SystemManagement* system_management;
            
            /**
             * Vector of pointer to PorousLayer objects.
             */
            std::vector< boost::shared_ptr< FuelCellShop::Layer::PorousLayer<dim> > > porous_layers;
            
            /**
             * Pointer to operating conditions class
             */
            FuelCell::OperatingConditions* opCond;
            
        };
	
	//-------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------
	
       /**
         * Class used to evaluate the capillary pressure at PorousLayer (\em viz., GDL, MPL and CL) DoF point in the
         * finite element mesh.
         * 
         * Child of DataPostprocessorScalar (deal.II post-processing class used in DataOut class).
         * 
         * See step-29, step-32 and step-33 in the deal.II website for more details about this class.
         * 
         * 
         * <h3>Usage</h3>
         * 
         * This class is used in the data_out section of two phase application. It requires a pointer to SystemManagement 
         * object, a vector of pointer to the PorousLayer objects for which you would like to compute
         * the current density, and pointer to OperatingConditions object.
         *
         * @code
         * AppCathode<dim>::data_out(const std::string& filename,
         *                           const FuelCell::ApplicationCore::FEVectors& src)
         * 
         * // --- Find solution and assign solution interpretations ---
         * (...)
         * 
         * // --- Create vector of PostProcessing objects ---
         * std::vector< DataPostprocessor<dim>* > PostProcessing;
         * 
         * // --- relative_humidity ---
         * // First create a vector of PorousLayer objects. We have CGDL, CMPL and CCL in the cathode application.
         * std::vector< boost::shared_ptr< FuelCellShop::Layer::PorousLayer<dim> > > porous_layers;
         * porous_layers.push_back(CGDL);
         * porous_layers.push_back(CMPL);
         * porous_layers.push_back(CCL);
         *
         * // Then create the post processing object and push it to the vector.
         * FuelCellShop::PostProcessing::CapillaryPressureDataOut<dim> capillary_pressure(&system_management, porous_layers, &OC);
         * PostProcessing.push_back(&capillary_pressure);
         *
         * // Outputting
         * DoFApplication<dim>::data_out(filename,
         *                               solution, // This is initialized in the initial part of the data_out function
         *                               system_management.get_solution_names(),
         *                               PostProcessing);
         * @endcode
         * 
         * \author J. Zhou
         * M. Secanell, 2014
         */
        template <int dim>
        class CapillaryPressureDataOut 
        : 
        public dealii::DataPostprocessorScalar<dim>
        {
        public:
            
            ///@name Constructors, destructor, and initalization
            //@{
            /**
             * Constructor of our class
             */
            CapillaryPressureDataOut(FuelCell::SystemManagement* ,
                                     std::vector< boost::shared_ptr< FuelCellShop::Layer::PorousLayer<dim> > > pls,
                                     FuelCell::OperatingConditions* );
            /**
             * Destructor
             */
            virtual ~CapillaryPressureDataOut() {}
            
            /**
             * Flags to be updated in each cell when computing the solution.
             */
            virtual UpdateFlags get_needed_update_flags() const;
            
            /**
             * 
             */                                                            
            virtual void compute_derived_quantities_vector (const std::vector< Vector< double > > &uh,
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
            FuelCell::SystemManagement* system_management;
                                    
            /**
             * Vector of pointer to PorousLayer objects.
             */
            std::vector< boost::shared_ptr< FuelCellShop::Layer::PorousLayer<dim> > > porous_layers;
            
            /**
             * Pointer to operating conditions class
             */
            FuelCell::OperatingConditions* opCond;

        };
        
        /**
         * Class used to output saturation inside the layer.
         */
        template <int dim>
        class SaturationDataOut 
        : 
        public dealii::DataPostprocessorScalar<dim>
        {
        public:
            
            ///@name Constructors, destructor, and initalization
            //@{
            /**
             * Constructor of our class
             */
            SaturationDataOut(FuelCell::SystemManagement* ,
                                     std::vector< boost::shared_ptr< FuelCellShop::Layer::PorousLayer<dim> > > pls,
                                     FuelCell::OperatingConditions* );
            /**
             * Destructor
             */
            virtual ~SaturationDataOut() {}
            
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
            virtual void compute_derived_quantities_vector (const std::vector< Vector< double > > &uh,
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
            FuelCell::SystemManagement* system_management;
                                    
            /**
             * Vector of pointer to PorousLayer objects.
             */
            std::vector< boost::shared_ptr< FuelCellShop::Layer::PorousLayer<dim> > > porous_layers;
            
            /**
             * Pointer to operating conditions class
             */
            FuelCell::OperatingConditions* opCond;

        };        
    }
    
}
    


#endif
