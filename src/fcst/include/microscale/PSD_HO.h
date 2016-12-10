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
//    - $ $
//
//---------------------------------------------------------------------------
#ifndef _FUELCELLSHOP__HO__PSD_H
#define _FUELCELLSHOP__HO__PSD_H

// Include OpenFCST routines:
#include <microscale/PSD_base.h>

using namespace dealii;


namespace FuelCellShop
{
    
    
    namespace MicroScale
    {
    /**
     * 
     * @brief Hydrophobic Pore Size Distribution
     * 
     * 
     * Based on the results of the mecury intrusion experiment, this class calculates the coefficients such as
     * relative liquid permeability, permeability and knudsen radius...
     * 
     * This class is the child class of the all the BasePSD class which reimplements all the functions that are decleared in the BasePSD class.
     * 
     * 
     * <h3> Input parameters </h3> 
     * The input parameters for this class are:Mode probability global,
     * Mode characteristic radius global,
     * Mode width global,
     * Volume fraction Hydrophobic,
     * @code
     * subsection PSD parameters
     *   subsection BasePSD
     *     set psd type = HOPSD
     *     set Mode probability global = 0.72, 0.28
     *     set Mode characteristic radius global = 34.0, 14.2
     *     set Mode width global = 0.35, 1.0
     *     set Gamma = 0.24 
     *     set contact_angle = 1.396
     *     set lamda = 1.0
     *     set probability P_b = 1
     *     subsection HIPSD 
     *          set Hydrophobic Mode probability global = 0.72, 0.28
     *          set Hydrophobic Mode characteristic radius global = 34.0, 14.2
     *          set Hydrophobic Mode width global = 0.35, 1.0
     *     end
     *   end
     * end
     * @endcode
     * 
     * <h3> Usage details</h3>
     * 
     * If you want to use the object of the HOPSD, use the following code,
     * but you should not create an object of the HOPSD out of the PSD scope.
     * If you want to use the PSD in the layer level, use the code in BasePSD.
     * 
     * @code
     * //Create an object of TemplateClass
     * FuelCellShop::MicroScale::HOPSD<dim> PSD_object; 
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
     *  Marc Secanell
     * @date 2014
     */
        template <int dim>
        class HOPSD : public BasePSD<dim>
        {
        public:            
            ///@name Constructor,Initalization, Declarations
            
            //@{
            /**
             * Constructor
             * 
             * \warning For internal use only.
             * 
             * Constructor used only to create a prototype. Do not use
             * in general since this will not include the name of the section
             * in the parameter file you need. 
             */
            HOPSD();    
            
            /**
             * Constructor.
             */
            HOPSD (std::string name);
            
            /**
             * Destructor.
             */
            virtual ~HOPSD();
             
             /**
              * Member function used to read in data and initialize the necessary data
              * to compute the coefficients.
              */
             void initialize ( ParameterHandler &param );
             
             /**
              * Declare all necessary parameters in order to compute the coefficients
              * 
              */    
             void declare_parameters (ParameterHandler &param) const;
             
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
                 
                 critical_radius_is_initialized = false;
                 saturation_is_initialized = false;
                 critical_radius_computed.clear();
                 saturation_computed.clear();
             }
             
             /**
              * Member function used to set the critical radius [\p nm] at every quadrature point
              * inside the cell. This function should particulary be used in the case of non-isothermal application.
              * It needs to be implemented at each iteration time to reduce the computational time.
              */
             
             inline void set_critical_radius()
             {
                 get_critical_radius(critical_radius_computed);
                 
                 critical_radius_is_initialized = true;
             }
             
             /**
              * Member function used to set the saturation at every quadrature point
              * inside the cell. This function should particulary be used in the case of non-isothermal application.
              * It needs to be implemented at each iteration time to reduce the computational time.
              */
             
             inline void set_saturation()
             {
                 get_saturation(saturation_computed);
                 
                 saturation_is_initialized = true;
             }
             /**
              * Concrete name used for objects of this class. This name is used when
              * setting up the subsection where the data is stored in the input file.
              * 
              * The data will be store under
              * \code
              * subsection name_specified_in_constructor
              *    set psd type = HOPSD # <-here I select the type of object of type psd
              *    subsection HOPSD # <- this is the concrete_name for this class
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
                 return typeid(HOPSD<dim>);
             }
            
            /**
             * This function is used to compute the saturation by using PSD
             * The saturation iof a porous medium can be obtained by integrating 
             * the change of cumulative pore volume fraction as a function of 
             * effective pore radius over the whole pore size domain
             * 
             * \f$ S_{HO} = F_{HO}\sum_{k} \frac{f_{HO,k}}{2} \left[1 - \text{erf}\left(\frac{\text{ln}(r_{c,HO}) - \text{ln}(r_{HO,k})}{s_{HO,k}\sqrt{2}} \right)  \right]  \f$
             */
            
            virtual  void get_saturation(std::vector<double>& S) const ;
            virtual  void get_derivative_saturation(std::vector<double>& ) const ;
            
            /**
             * This function is used to compute the Hydrophobic saturated_permeability by using PSD
             * 
             * \f$ k_{sat} = \frac{1}{8}\left[\frac{\varepsilon_{o}}{\lambda}\right]^{2} \sum_{k} \,\exp{(-2\,s_{k}^{2})}\,r_{k}^{2}f_{k} \f$ 
             */
            
            virtual  void get_global_saturated_permeability(double& saturated_permeability) const ;
            virtual  void get_global_saturated_permeability(const double,double& ) const ;
            /**
             * This function is used to compute the Hydrophobic liquid_permeability by using PSD
             * 
             *  \f$ k_{L,HO} =  \frac{F_{HO}}{16}\left[\frac{\varepsilon_{o}\,S_{e}}{\lambda}\right]^{2}\sum_{k} \,\exp{(-2\,s_{HO,k}^{2})}\,r_{HO,k}^{2}f_{HO,k}  \left[-\text{erf}\left(\frac{\ln{(r_{c,HO})} - \ln{(r_{HO,k})}}{s_{HO,k}\sqrt{2}} - s_{HO,k}\sqrt{2} \right) + 1 \right]  \f$
             * 
             */
            
            virtual  void get_pore_HO_liquid_saturated_permeability(std::vector<double>& saturated_HO_permeability) const ;
            
            virtual  void get_derivative_pore_HO_liquid_saturated_permeability(std::vector<double>& ) const ;
            
            virtual  void get_derivative_pore_HO_liquid_saturated_permeability(const double porosity, const std::vector<double> S,
                                                                               const std::vector<double> ds_dp,
                                                                               std::vector<double>& ) const ;
            
            virtual  void get_pore_HO_liquid_saturated_permeability(const double porosity, const std::vector<double> S,
                                                                    std::vector<double>& saturated_HO_permeability) const ;
            
            /**
             * This function is used to compute the Hydrophobic liquid_permeability by using PSD
             */
            
            virtual  void get_relative_liquid_permeability(std::vector<double>& liquid_permeability) const ;
            
            virtual  void get_derivative_relative_liquid_permeability(std::vector<double>& ) const ;
            
            /**
             * This function is used to compute the Hydrophobic gas_permeability by using PSD 
             * 
             *  \f$ k_{G,HO} =  \frac{F_{HO}}{16}\left[\frac{\varepsilon_{o}\,\left(1 - S_{e}\right)}{\lambda}\right]^{2}\sum_{k} \,\exp{(-2\,s_{HO,k}^{2})}\,r_{HO,k}^{2}f_{HO,k}  \left[\text{erf}\left(\frac{\ln{(r_{c,HO})} - \ln{(r_{HO,k})}}{s_{HO,k}\sqrt{2}} - s_{HO,k}\sqrt{2} \right) + 1 \right]  \f$ 
             */
            
            virtual  void get_pore_HO_gas_saturated_permeability(std::vector<double>& saturated_HO_permeability) const ;
            
            virtual void  get_pore_HO_gas_saturated_permeability(const double porosity, const std::vector<double> S,
                                                                 std::vector<double>& saturated_HO_permeability) const ;
            
            /**
             * This function is used to compute the Hydrophobic gas_permeability by using PSD 
             */
            
            virtual  void get_relative_gas_permeability(std::vector<double>& gas_permeability) const ;
            
            /**
             * This function is used to compute the Hydrophobic liquid_gas_interfacial_surface by using PSD
             * 
             * \f$ \frac{a(r)_{HO}}{V_{T}} =  \frac{a(r)_{c}}{a_{T}} \left( 1 - \frac{a(r)_{c}}{a_{T}}\right)\, F_{HO}\sum_{k}\frac{f_{k,HO}\exp{\left(\frac{s_{k,HO}^{2}}{2}  \right)}}{8\,r_{k,HO}}   \left[1 - \text{erf}\left(\frac{\text{ln}(r_{cr}) - \text{ln}(r_{k,HO})}{s_{k,HO}\sqrt{2}} + \frac{s_{k,HO}\sqrt{2}}{2}    \right)   \right]   \f$
             */
            
            virtual  void get_liquid_gas_interfacial_surface(std::vector<double>& HO_liquid_gas_interfacial_surface_a) const ;
            
            virtual  void get_liquid_gas_interfacial_surface_withoutPb(std::vector<double>& HO_liquid_gas_interfacial_surface) const ;
            
            /**
             * This function is used to compute the Hydrophobic liquid_gas_interfacial_surface by using PSD
             * 
             * \f$ \frac{a(r)_{HO}}{V_{T}} =  \frac{a(r)_{c}}{a_{T}} \left( 1 - \frac{a(r)_{c}}{a_{T}}\right)\, F_{HO}\sum_{k}\frac{f_{k,HO}\exp{\left(\frac{s_{k,HO}^{2}}{2}  \right)}}{8\,r_{k,HO}}   \left[1 - \text{erf}\left(\frac{\text{ln}(r_{cr}) - \text{ln}(r_{k,HO})}{s_{k,HO}\sqrt{2}} + \frac{s_{k,HO}\sqrt{2}}{2}    \right)   \right]   \f$
             */
            
            virtual  void get_derivative_liquid_gas_interfacial_surface_increment(std::vector<double>& ) const ;
            virtual  void get_derivative_liquid_gas_interfacial_surface(std::vector<double>& ) const ;
            
            /**
             * This function is used to compute the Hydrophobic pore_wetted_wall by using PSD
             * 
             * \f$ a_{wall,HO} = \sum_{k} \frac{F_{HO}\,f_{k,HO}}{r_{k,HO}}\text{exp}\left(\frac{s_{k,HO}^{2}}{2} \right)     \left[ 1 - \text{erf} \left( \frac{\text{ln}(r_{c,HO}) - \text{ln}(r_{k,HO})}{s_{k,HO}\sqrt{2}} + \frac{s_{k,HO}}{\sqrt{2}} \right)   \right]  \f$
             */
            
            virtual  void get_pore_HO_wetted_wall_surface_area(std::vector<double>& HO_wetted_wall_surface_area) const ;
            
            /**
             * This function is used to compute the Hydrophobic pore_wetted_wall by using PSD
             */
            
            virtual  void get_wetted_wall_surface_area(std::vector<double>& wetted_wall_surface_area) const ;
            
            /**
             * This function is used to compute the Hydrophobic knudsen_radius C2 by using PSD
             * 
             * \f$ C_{2} = F_{HO} \sum_{k} \frac{f_{HO,k}}{2} \left[1 + \text{erf}\left(\frac{\text{ln}(r_{c,HO}) - \text{ln}(r_{HO,k})}{s_{HO,k}\sqrt{2}} \right)  \right]  \f$
             */
            
            virtual  void get_pore_knudsen_radius_C2(std::vector<double>& knudsen_radius_C2) const ;
            
            /**
             * This function is used to compute the Hydrophobic knudsen_radius C4 by using PSD
             * 
             * \f$ C_{4} = \sum_{k} \frac{F_{HO}\,f_{k,HO}}{r_{k,HO}}\text{exp}\left(\frac{s_{k,HO}^{2}}{2} \right) \left[ 1 + \text{erf} \left( \frac{\text{ln}(r_{c,HO}) - \text{ln}(r_{k,HO})}{s_{k,HO}\sqrt{2}} + \frac{s_{k,HO}}{\sqrt{2}} \right)   \right]    \f$
             */
            
            virtual  void get_pore_knudsen_radius_C4(std::vector<double>& knudsen_radius_C4) const ;
            
            /**
             * This function is used to compute the Hydrophobic knudsen_radius by using PSD
             */
            
            virtual  void get_knudsen_radius(std::vector<double>& knudsen_radius) const ;
            
            /**
             * This function is used to compute the diffusivity by using PSD
             */
            virtual  void get_diffusivity() const ;
            
            
            /**
             * This function is used to compute the diffusivity by using PSD
             * 
             * @note the p_c is the capillary pressure which is one of the 
             * solution variable in the equation liquid water
             * therefore we need to use set_capillary function to compute 
             * the p_c at each quatrature points
             * Units of the critical_radius is in  \f$  r_c \quad \left[ um \right] \f$
             * and the capillary pressure 
             * imposed by the solution is \f$ P_c \quad \left[ Pascal \right] \f$
             */
            
            virtual  void get_critical_radius(std::vector<double>& dst) const;
            virtual void get_derivate_critical_radius(std::vector<double>& dst) const;
            virtual  void get_maximum_cross_sectional_areas(double&) const;
            
            /**
             * This function is used to plot PSD configuration
             */
            virtual void get_PSD_plot(const std::vector<double> ,std::vector<double>&) const;
            
            //@}
            
        protected:
            
            ///@name Instance Delivery (Replica creator)
            //@{
            /**
             * This member function is used to create an object of type psd
             * 
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::MicroScale::BasePSD <dim>> create_replica (const std::string &psd_section_name)
            {
                return boost::shared_ptr<FuelCellShop::MicroScale::BasePSD <dim>> (new FuelCellShop::MicroScale::HOPSD<dim> (psd_section_name));
            }   
            //@}         
            
            ///@name Instance Delivery (Prototype)
            //@{
            /**
             *  PROTOTYPE is the pointer is the dynamic pointer pointing to the HOPSD class itself.
             */            
            static HOPSD<dim> const* PROTOTYPE;
            //@}
            
            //////////
            // DATA //
            //////////
            
            ///@name Psd properties
            //@{
            
            /**
             * The f_k is the contribution of the log-normal distribution k to the total PSD
             */            
            std::vector<double> fHO_k;
            
            /**
             * The r_k is the characteristic pore size of the distribution k
             */            
            std::vector<double> rHO_k;            
            
            /**
             * The s_k is the spread of the distribution k
             */
            std::vector<double> sHO_k;
            
            double contact_angle_HO;
            
            /** Temperature at every quadrature point inside the cell. */
            
            SolutionVariable T_vector;
            
            /** Capillary pressure at every quadrature point inside the cell. */
            
            SolutionVariable Capillary_pressure_vector;
            
            /** Constant capillary pressure only for unit_test use */
            
            double pressure_c;
            
            /** Critical_radius_computed by the get_critical_radius function */
            
            std::vector<double> critical_radius_computed;
            
            /** 
             * Check if the critical radius has already been computed by set_critical_radius function.
             * It has to be reset each time when the solution variable capillary pressure is updated
             */
            
            bool critical_radius_is_initialized;
            
            /** Saturation_computed by the get_critical_radius function */
            std::vector<double> saturation_computed;
            
            /** 
             * Check if the saturation has already been computed by set_saturation function.
             * It has to be reset each time when the solution variable capillary pressure is updated
             */
            bool saturation_is_initialized;

            //@}
        };
        
    } // PSD
    
}  // FuelCellShop

#endif