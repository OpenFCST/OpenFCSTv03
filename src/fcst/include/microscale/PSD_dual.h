//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: PSD_base.h
//    - Description: Base class for pore size distribution model.
//    - Developers: 2009-13 by Marc Secanell, University of Alberta
//                  2013-14 by Jie Zhou, University of Alberta
//
//---------------------------------------------------------------------------
#ifndef _FUELCELLSHOP__DUAL__PSD_H
#define _FUELCELLSHOP__DUAL__PSD_H

// Include OpenFCST routines:
#include <microscale/PSD_HI.h>
#include <microscale/PSD_HO.h>



using namespace dealii; 

namespace FuelCellShop
{
    
    
    namespace MicroScale
    {
    /**
     * 
     * @brief Dual Pore Size Distribution
     * 
     * 
     * Based on the results of the mecury intrusion experiment, this class calculates the coefficients such as
     * relative liquid permeability, permeability and knudsen radius...
     * 
     * This class is the child class of the all the BasePSD class which reimplements all the functions that are decleared in the BasePSD class.
     * 
     * <h3> Input parameters </h3> 
     * The input parameters for this class are:Mode probability global,
     * Mode characteristic radius global,
     * Mode width global,
     * Volume fraction Hydrophobic,
     * Volume fraction Hydrophilic...
     * All data below is from the Hydrophilic unit test
     * @code
     * subsection PSD parameters
     *   subsection BasePSD
     *     set psd type = DualPSD
     *     set Mode probability global = 0.72, 0.28
     *     set Mode characteristic radius global = 34.0, 14.2
     *     set Mode width global = 0.35, 1.0
     *     set Gamma = 0.24 
     *     set contact_angle = 1.396
     *     set lamda = 1.0
     *     set probability P_b = 1
     *     subsection HIPSD 
     *          set Hydrophilic Mode probability global = 0.72, 0.28
     *          set Hydrophilic Mode characteristic radius global = 34.0, 14.2
     *          set Hydrophilic Mode width global = 0.35, 1.0
     *     end
     *     subsection HOPSD
     *          set Hydrophobic Mode probability global = 0.72, 0.28
     *          set Hydrophobic Mode characteristic radius global = 34.0, 14.2
     *          set Hydrophobic Mode width global = 0.35, 1.0
     *     end
     *   end
     * end
     * @endcode
     * 
     * <h3> Usage details</h3>
     * If you want to use the object of the DualPSD, use the following code,
     * but you should not create an object of the DualPSD out of the PSD scope.
     * If you want to use the PSD in the layer level, use the code in BasePSD.
     * 
     * @code
     * //Create an object of TemplateClass 
     * FuelCellShop::MicroScale::DualPSD<dim> PSD_object; 
     * @endcode 
     * 
     * <h3> References </h3>
     *
     * [1]  Pedro Abdiel Mateo Villanueva, 
     * A MIXED WETTABILITY PORE SIZE DISTRIBUTION MODEL FOR 
     * THE ANALYSIS OF WATER TRANSPORT IN PEMFC MATERIALS, M. Sc. thesis, University of Alberta, 2013
     * 
     * @author J. Zhou
     * 
     *   Marc Secanell
     * @date 2014
     */
        template <int dim>
        class DualPSD : public BasePSD<dim>
        {
        public:
            
            ///@name Constructors, declarations and Initalization
            //@{
            /** 
             * Consturctor
             */
            DualPSD();
            
            /**
             * Constructor.
             */
            DualPSD (std::string name);
            
            /**
             * Destructor.
             */
            virtual ~DualPSD();
            
            /**
             * Declare parameters for a parameter file.
             */ 
            void declare_parameters (ParameterHandler &param) const;
            /**
             * Member function used to read in data and initialize the necessary data
             * to compute the coefficients.
             */
            void initialize (ParameterHandler &param );
            
            /**
             * Member function used to set the temperature [\p Kelvin] at every quadrature point
             * inside the cell. This function should particulary be used in the case of non-isothermal application.
             */
            inline void set_temperature (const SolutionVariable& T_in)
            {
                Assert( T_in.get_variablename() == temperature_of_REV, ExcMessage("Wrong solution variable passed in PSD::set_temperature.") );
                this->T_vector = T_in;
            }
            
            /**
             * Member function used to set the capillary pressure [\p psi] at every quadrature point
             * inside the cell. This function should particulary be used in the case of non-isothermal application.
             */
            inline void set_capillary_pressure (const SolutionVariable& C_in)
            {
                Assert( C_in.get_variablename() == capillary_pressure, ExcMessage("Wrong solution variable passed in PSD::capillary_pressure.") );
                this->Capillary_pressure_vector = C_in;
                
                psd_hi.set_capillary_pressure(C_in);
                psd_ho.set_capillary_pressure(C_in);
            }
            
            /**
             * Member function used to set the critical radius [\p nm] at every quadrature point
             * inside the cell. This function should particulary be used in the case of non-isothermal application.
             * It needs to be implemented at each iteration to reduce the computational time.
             */
            
            inline void set_critical_radius()
            {
                psd_hi.set_critical_radius();
                psd_ho.set_critical_radius();
                
            }
            
            /**
             * Member function used to set the saturation at every quadrature point
             * inside the cell. This function should particulary be used in the case of non-isothermal application.
             * It needs to be implemented at each iteration to reduce the computational time.
             */
            
            inline void set_saturation()
            {
                psd_hi.set_saturation();
                psd_ho.set_saturation();
            }
            

            //@}
             
             ///@name concrete_name
             //@{
             /**
              * Concrete name used for objects of this class. This name is used when
              * setting up the subsection where the data is stored in the input file.
              * 
              * The data will be store under
              * \code
              * subsection name_specified_in_constructor
              *    set psd type = DualPSD # <-here I select the type of object of type psd
              *    subsection DualPSD # <- this is the concrete_name for this class
              *       set all info relevant to this object
              *    end
              * end
              * \endcode
              */
             static const std::string concrete_name;
             
             //@}
             
            ///@name Accessors and info
             //@{
             /**
              * This member function returns a type_info object with the name of the 
              * base layer type the inherited class belongs to, i.e.
              * HIPSD
              * HOPSD
              * DualPSD
              * 
              * Note that this is necessary if we want to find out not the name of the actual class which can be obtain using
              * @code const std::type_info& name = typeid(*this) @endcode
              * but the name of the parent class.
              * 
              * @note Do not re-implement this class in children classes
              */
             const std::type_info& get_base_type() const
             {
                 return typeid(DualPSD<dim>);
             }
             
            /**
             * This function is used to compute the saturation by using PSD
             * The saturation iof a porous medium can be obtained by integrating 
             * the change of cumulative pore volume fraction as a function of 
             * effective pore radius over the whole pore size domain
             * 
             * \f$ S_{total} = S_{HI} + S_{HO}  \f$
             */

            virtual  void get_saturation(std::vector<double>& S) const ;
            virtual  void get_derivative_saturation(std::vector<double>& S) const ;
            /**
             * This function is used to compute the saturated_permeability by using PSD
             * 
             * \f$ k_{sat} = \frac{1}{8}\left[\frac{\varepsilon_{o}}{\lambda}\right]^{2} \sum_{k} \,\exp{(-2\,s_{k}^{2})}\,r_{k}^{2}f_{k} \f$ 
             */
            
            virtual  void get_global_saturated_permeability(double& saturated_permeability) const ;
            
            /**
             * This function is used to compute the liquid_permeability by using PSD
             * 
             * \f$ k_{r,L} = \frac{k_{L,HI} + k_{L,HO}}{k_{sat}}  \f$
             */

            virtual  void get_relative_liquid_permeability(std::vector<double>& liquid_permeability) const ;
            virtual  void get_derivative_relative_liquid_permeability(std::vector<double>& liquid_permeability) const ;
            /**
             * This function is used to compute the gas_permeability by using PSD
             * 
             * \f$  k_{r,G} = \frac{k_{G,HI} + k_{G,HO}}{k_{sat}}  \f$
             */
            
            virtual  void get_relative_gas_permeability(std::vector<double>& gas_permeability) const ;
            
            /**
             * This function is used to compute the liquid_gas_interfacial_surface by using PSD
             * 
             * \f$  a_{LV} = \left[ \frac{a(r)_{HI}}{V_{T}} + \frac{a(r)_{HO}}{V_{T}} \right]\frac{1}{\lambda_{2}}   \f$
             */
            
            virtual  void get_liquid_gas_interfacial_surface(std::vector<double>& liquid_gas_interfacial_surface) const ;
            
            /**
             * This function is used to compute the liquid_gas_interfacial_surface by using PSD
             * 
             * \f$  a_{LV} = \left[ \frac{a(r)_{HI}}{V_{T}} + \frac{a(r)_{HO}}{V_{T}} \right]\frac{1}{\lambda_{2}}   \f$
             */
            
            virtual  void get_derivative_liquid_gas_interfacial_surface(std::vector<double>&) const ;
            
            /**
             * This function is used to compute the pore_wetted_wall by using PSD
             * 
             * \f$  a_{wall} = a_{wall,HI} + a_{wall,HO}    \f$
             */
            
            virtual  void get_wetted_wall_surface_area(std::vector<double>& wetted_wall_surface_area) const ;
            
            /**
             * This function is used to compute the knudsen_radius by using PSD
             * 
             * \f$  r_{Kn} = \frac{C_{1} + C_{2}}{C_{3} + C_{4}}   \f$ 
             */
            
            virtual  void get_knudsen_radius(std::vector<double>& knudsen_radius) const ;
            
            /**
             * This function is used to compute the diffusivity by using PSD
             */
            virtual  void get_diffusivity() const ;
            
            virtual void get_maximum_cross_sectional_areas(double&) const;
            
            /**
             * This function is used to plot PSD configuration
             */
            virtual void get_PSD_plot(const std::vector<double> ,std::vector<double>&,std::vector<double>&,std::vector<double>&) const;
            
            
            //@}
            
        protected:
            ///@name Instance Delivery (Replica creator)
            //@{
            /**
             * This member function is used to create an object of type psd
             * 
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::MicroScale::BasePSD<dim> > create_replica (const std::string &psd_section_name)
            {
                return boost::shared_ptr<FuelCellShop::MicroScale::BasePSD<dim> > (new FuelCellShop::MicroScale::DualPSD<dim> (psd_section_name));
            }   
            //@}           
            
            ///@name Instance Delivery (Prototype)
            //@{
            /**
             * PROTOTYPE is the pointer is the dynamic pointer pointing to the DualPSD class itself.
             */            
            static DualPSD<dim> const* PROTOTYPE;
            //@}
            
            //////////
            // DATA //
            //////////
            
            ///@name PSD properties
            //@{
            
            /** Creating a pointer of PSD hydrophilic class. */
            
            HIPSD<dim> psd_hi;
            
            /** Creating a pointer of PSD hydrophobic class. */
           
            HOPSD<dim> psd_ho;
            
            double pressure_c;
            
            /** Temperature at every quadrature point inside the cell. */
            SolutionVariable T_vector;
            
            /** Capillary pressure at every quadrature point inside the cell. */
            SolutionVariable Capillary_pressure_vector;

            //@}
        };
        
    } // PSD
    
}  // FuelCellShop

#endif