//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2014 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: multi_scale_CL.h
//    - Description: Class characterizing the catalyst layer and defining effective properties
//    - Developers: M. Secanell, Peter Dobson, Philip Wardlaw, Michael Moore, Madhur Bhaiya
//    - $Id: multi_scale_CL.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__LAYER__MULTISCALE_CL__LAYER_H
#define _FUELCELLSHOP__LAYER__MULTISCALE_CL__LAYER_H

//Include Boost classes
#include <boost/smart_ptr.hpp>

// Include FCST classes
#include <utils/fcst_utilities.h>
#include <layers/conventional_CL.h>
#include <microscale/micro_scale_base.h>

//Include STL
#include <stdexcept>
#include <map>


namespace FuelCellShop
{
    namespace Layer
    {
        /**
         * This class characterizes a catalyst layer and uses this information
         * to compute effective transport properties and interfacial areas for phase
         * change or electrochemical reactions.
         * 
         * @author M. Secanell, 2009-13
         * @author P. Dobson, 2009-11
         * @author M. Moore, 2010-12
         * @author M. Bhaiya, 2011-13
         * @author P. Wardlaw, 2012-14
         * 
         */
        template <int dim>
        class MultiScaleCL :
        public ConventionalCL<dim>
        {
        public:
            /**
             * Concrete name used for objects of this class. This name is used when
             * setting up the subsection where the data is stored in the input file.
             * 
             * The data will be store under
             * \code
             * subsection name_specified_in_constructor
             *    set Material id = 2
             *    set Catalyst layer type = DummyCL # <-here I select the type of object of type CatalystLayer
             *    subsection DummyCL # <- this is the concrete_name for this class
             *       set all info relevant to this object
             *    end
             * end
             * \endcode
             */
            static const std::string concrete_name;
            
	    
            /**
             * Function for setting current cell_id from applications.
             *
             * - \param unsigned int id is the id of the current cell from the application's perspective
             *
             * @note Reimplemented here as cell_id is needed for poly disperse agglomerates
             */
            virtual void set_cell_id(const unsigned int& id){
                cell_id_ = id;
            }

            /**
             * Destructor
             */
            ~MultiScaleCL();
            
            ///@name Accessors and info
            //@{

            /**
             * This member function will use a FuelCellShop::Kinetics class in order to compute the current density
             * production in the CL
             * Returns the current density at each quadrature point.  If the current is averaged over the cell,
             * the average current is assigned to each quadrature point in the cell.
             *
             * @param current is an empty vector to be filled with current density values.  Will be resized to the size of
             *              the solution given in set_solution()
             *
             */
            virtual void current_density ( std::vector<double>& current );
            
            /**
             * This member function will use a FuelCellShop::Kinetics class in order to compute the current density
             * production in the CL
             * Returns the current density at each quadrature point.  If the current is averaged over the cell,
             * the average current is assigned to each quadrature point in the cell.
             *
             * @param current is an empty vector to be filled with current density values.  Will be resized to the size of
             * 				the solution given in set_solution()
             * 
             * @param effectiveness is an empty vector to be filled with micro scale effectiveness values.
             *
             */
            virtual void current_density (std::vector<double>& current, std::vector<double>& effectiveness );           
            
            /**
             * This member function will use a FuelCellShop::Kinetics class in order to compute the derivative of the current density
             * with respect to the variables setup in \ref set_derivative_flags (std::vector< std::string > &flags)
             *
             */
            virtual void derivative_current_density ( std::map< VariableNames, std::vector<double> >& );
            
            /**
             * Print out composition and micro-structural properties of the catalyst layer
             */
            virtual void print_layer_properties() const;
            

            /*
             * Public member function which returns a boost::shared pointer of template type T.
             *
             *
             * \param T is the type of object you are requestiong, e.g., FuelCellShop::Material::CatalystSupportBase
             *
             * <h3>Usage</h3>
             * To be used by MicroScale objects in order to get resources such as materials or
             * kinetics. Currently the MicroScaleCL provides resources of the following types:
             *
             * -FuelCellShop::Kinetics::BaseKinetics
             * -FuelCellShop::Material::CatalystBase
             * -FuelCellShop::Material::PolymerElectrolyteBase
             * -FuelCellShop::Material::CatalystSupportBase
             *
             * The above types may be used as template arguments for this function.
             *
             * @code
             *
             * //We wish to link the following kinetics pointer to the
             * //kinetics object pointed to by the MultiScaleCL
             *
             * boost::shared_ptr< FuelCellShop::Kinetics::BaseKinetics > kinetics;
             *
             * //To get a kinetics pointer from an instance of MultiScaleCL we do
             * //the following:
             *  kinetics = layer->get_resource<FuelCellShop::Kinetics::BaseKinetics>();
             *
             *
             * @endcode
             *
             */
            template<typename T>
            inline boost::shared_ptr<T>
            get_resource(){

                boost::shared_ptr<void> ptr;


                if(typeid(T) == typeid(FuelCellShop::Kinetics::BaseKinetics))
                {
                    ptr = this->kinetics;
                }
                else if(typeid(T) == typeid(FuelCellShop::Material::CatalystBase))
                {
                    ptr = this->catalyst;
                }
                else if(typeid(T) == typeid(FuelCellShop::Material::PolymerElectrolyteBase))
                {
                    ptr = this->electrolyte;
                }
                else if(typeid(T) == typeid(FuelCellShop::Material::CatalystSupportBase))
                {
                    ptr = this->catalyst_support;
                }
                else
                {
                    std::string msg = std::string(typeid(*this).name()) + " does not have object type " + std::string(typeid(T).name());
                    throw std::runtime_error(msg);
                }

                return boost::static_pointer_cast<T>(ptr);

            }


            /*
             * A list of properties that the MultiScaleCl shares publicly using
             * the get_properties interface.
             *
             * Used by MicroScale objects
             */
            enum Properties{
                void_fraction=0,
                solid_fraction,
                ionomer_fraction,
                active_area_scaled,
                pressure,
                cell_id
            };


            /*
             * Public member function for getting properties of the MultiScaleCL
             * needed by MicroScale object.
             */
            inline std::map<Properties, double> get_properties(){

                std::map<Properties, double> properties;

                properties[solid_fraction] = this->epsilon_S.at(this->local_material_id());
                properties[void_fraction] = this->epsilon_V.at(this->local_material_id());
                properties[ionomer_fraction] = this->epsilon_N.at(this->local_material_id());
                //The following active area scaling is part of the scaling described by equation 2.90 of Philip Wardlaw's MSc. thesis, and previously by Peter Dobson
                properties[active_area_scaled] = this->Av.at(this->local_material_id()) / (1.0 - this->epsilon_V.at(this->local_material_id()));
                properties[pressure] = this->constant_solutions.at(total_pressure); // <- Not always the case that it is constant...
                properties[cell_id] = double(cell_id_);
                return properties;
            }

            /**Method for getting coverages from micro scale objects
             * Overloaded here since kinetics cannot be used directly
             */
            virtual SolutionMap get_coverages();


           
        protected:
	  
	  
	   ///@name Constructors, destructor, and initalization
            //@{    
            /**
             * Prototype Constructor
             * 
             */
            MultiScaleCL ( );
            
            /**
             * Constructor
             * 

             */
            MultiScaleCL ( std::string name );
            
            
            
            /**
             * Declare all necessary parameters in order to compute the coefficients
             * 
             */            
            void declare_parameters ( ParameterHandler &param ) const
            {
               declare_parameters(this->name, param);                
            }

            /**
             * Member function used to read in data and initialize the necessary data
             * to compute the coefficients.
             */
            void initialize ( ParameterHandler &param );

            //@}

            ///@name Constructors and declarations
            //@{
            /**
            * Constructor
            */
            //MultiScaleCL(const std::string& cl_section_name);

            /**
            * Constructor that also intialises the pointers to the catalyst, catalyst
            * support and electrolyte objects.
            */
            /*MultiScaleCL(const std::string& name,
                          FuelCellShop::Material::PolymerElectrolyteBase*,
                          FuelCellShop::Material::CatalystSupportBase*,
                          FuelCellShop::Material::CatalystBase*);   
            */
            /**
             * Declare parameters for a parameter file.
             * 
             * Parameters that can be declared are defined in:
             * 
             * @code 
             * subsection Fuel cell data
             *  (...)
             *   subsection Cathode catalyst layer          <- This is the name of the CL subsection specified in cl_section_name
             *     (...)
             *     NOTE: Here you will need info on ConventionalCL since information from ConventionalCL is used for MultiScaleCL
             *     (...)
             *     subsection MultiScaleCL                 <- This is the subsection specified by concrete_name
             *      set Average current in cell = false     # Decide whether to take the average current density in the cell
             *     end
             *   end
             * end
             * @endcode
             */
            void declare_parameters (const std::string& cl_section_name, 
                                     ParameterHandler &param) const;
                                     
            //@}
            
            ///@name Instance Delivery (Replica creator)
            //@{
            /**
             * This member function is used to create an object of type gas diffusion layer
             * 
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::Layer::CatalystLayer<dim> > create_replica (const std::string &cl_section_name)
            {
                return boost::shared_ptr<FuelCellShop::Layer::CatalystLayer<dim> > (new FuelCellShop::Layer::MultiScaleCL<dim> (cl_section_name));
            }   
            //@}           
            ///@name Instance Delivery (Prototype)
            //@{
            /**
             *
             */            
            static MultiScaleCL<dim> const* PROTOTYPE;
            //@}
            
            ///@name Internal member functions
            //@{
            /**
             * Compute porosity and volume fraction of solid and ionomer in the catalyst layer
             */
            void compute_void_fraction();
            
            /**
             * Creates the microscale objects populating the microscale object
             */
            void initialize_micro_scale(ParameterHandler &param);
            
            /**
             * Private member functions for solving current density given an microscale.
             */
            SolutionMap micro_scale_current(std::map<VariableNames ,SolutionVariable>& solutionMap, const unsigned int& sol_index, const unsigned int& thread_index);
            
            
            /**
             * Private member functions for solving for current derivatives in an averaging approach.
             */
            void solve_current_derivatives_average(std::map< VariableNames, std::vector<double> >& Dcurrent);
            
            /**
             * Private member functions for solving for current derivatives in a per node approach.
             */
            void solve_current_derivatives_at_each_node(std::map< VariableNames, std::vector<double> >& Dcurrent);
            //@}

            /** Boolean value to choose whether to average the current over the cell */
            bool average_cell_current;
            /** 
             * Vector of shared_ptr microscale objects used for calculating current density and current density derivatives
             */

            std::map<unsigned int, std::vector<boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase>>> micro;
            
            /*Solution map used for storing coverages */
            SolutionMap coverage_map;

            unsigned int cell_id_;

        };
        
    } // Layer
    
}  // FuelCellShop

#endif
    
