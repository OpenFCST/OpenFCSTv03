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
#ifndef _FUELCELLSHOP__BASE__PSD_H
#define _FUELCELLSHOP__BASE__PSD_H

// Include deal.II classes
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/fe/fe_values.h>

//Include STL
#include <cmath>
#include <iostream>
#include <math.h> 

// Include OpenFCST routines:
#include <application_core/fcst_variables.h>
#include <application_core/system_management.h>
#include <utils/fcst_utilities.h>
#include <utils/fcst_constants.h>

using namespace dealii;

namespace FuelCellShop
{
    
    
    namespace MicroScale
    {
    /**
     * 
     * @brief Pore Size Distribution
     * 
     * Based on the results of the mecury intrusion experiment, this class calculates effective transport properties
     * such as relative liquid permeability, permeability and knudsen radius, using a bundle of capilaries analogy as
     * described in reference [1].
     * 
     * @note This class is a virtual class, i.e. you should not create an object of this class. This class only
     * implements the declare_parameters and initialization of data sections for the children. All other
     * functions are virtual, i.e. this class only declares the interface for all child classes. Therefore, 
     * the PSD children classes need to reimplemented all the functions declared here. 
     * 
     * <h3> Input parameters </h3> 
     * 
     * The input parameter section declared in the PSD base class is shown below.
     * 
     * An exampe of a PSD section is provided below:
     * @code
     * subsection PSD parameters
     *   subsection BasePSD
     *     set psd type = HIPSD
     *     set Gamma = 0.24 
     *     set lamda = 1.0
     *     set probability P_b = 1
     *     set Contact angle = 90
     *     set Volume fraction Hydrophilic = 0.5
     *     set Volume fraction Hydrophobic
     *     set Mode probability global = 0.72, 0.28
     *     set Mode characteristic radius global = 34.0, 14.2
     *     set Mode width global = 0.35, 1.0
     *   end
     * end
     * @endcode
     * 
     * <h3> Usage details</h3>
     * 
     * To create a PSD object in the layer class, you need to call the static function declare_PSD_parameters
     * first to declare all the parameters in its children classes. Then using create_PSD with the concrete_name
     * to generate the shared pointer which you can use. The PSD_type is the concrete_name for the particular 
     * PSD you want to create.
     *         
     * @code
     * //--------- IN DECLARE_PARAMETERS ------------------------------------------------------
     * template <int dim>
     * void 
     * NAME::GasDiffusionLayer<dim>::declare_parameters(ParameterHandler& param)
     * {
     *   (...)
     *   // Declare section on the input file where all info will be stored.
     *   FuelCellShop::MicroScale::BasePSD<dim>::declare_PSD_parameters(param); 
     *   (...)
     * }
     * 
     * //--------- IN INITIALIZE ------------------------------------------------------
     * template <int dim>
     * void
     * NAME::GasDiffusionLayer<dim>::_initialize(ParameterHandler& param)
     * {   
     *      param.enter_subsection("PSD parameters");
     *          param.enter_subsection("BasePSD");
     *              PSD_type = param.get("psd type");
     *          param.leave_subsection();
     *      param.leave_subsection();
     * 
     *      PSD = FuelCellShop::MicroScale::BasePSD<dim>::create_PSD(PSD_type,param);
     * 
     *   (...)
     * }
     * @endcode
     * 
     * <h3> References </h3>
     *
     * [1]  Pedro Abdiel Mateo Villanueva, A MIXED WETTABILITY PORE SIZE DISTRIBUTION MODEL FOR 
     * THE ANALYSIS OF WATER TRANSPORT IN PEMFC MATERIALS, M. Sc. thesis, University of Alberta, 2013
     * 
     * @author 
     *         J. Zhou
     *  
     *  Marc Secanell 
     * @date 2014
     */
    template <int dim>
        class BasePSD : public Subscriptor
        {
        public:
            
            ///@name Destructor
            //@{
            /**
             * Destructor
             */
            virtual ~BasePSD() {}
            //@}
            
            ///@name Instance Delivery (Functions)
            //@{    
             
            /**
             * Function used to declare all the data necessary in the parameter files for
             * all BasePSD children.
             * 
             * This member function should be used instead of declare_parameters() when we want
             * to use a BasePSD pointer that selects the type of Psd to run at runtime.
             * 
             * \param psd_section_name Name of the section that will encapuslate all the information about the PSD
             * \param param ParameterHandler object used to store all information about the simulation. Used
             * to read the parameter file.
             * 
             * The parameter file would look as follows:
             * 
             * subsection psd_section_name
             *   set psd type = PSD_dual   # Options: PSD_HI | PSD_HO | PSD_dual | NonePSD
             * 
             */
            
            static void declare_PSD_parameters (ParameterHandler &param)
            {
                for (typename FuelCellShop::MicroScale::BasePSD<dim>::_mapFactory::iterator iterator = FuelCellShop::MicroScale::BasePSD<dim>::get_mapFactory()->begin(); 
                     iterator != FuelCellShop::MicroScale::BasePSD<dim>::get_mapFactory()->end(); 
                iterator++)
                     {
                         iterator->second->declare_parameters(param);
                     }
            }
            
            /**
             * 
             * Function used to select the appropriate CatalystLayer type as specified in the ParameterHandler under
             * line 
             * @code 
             * set psd type = PSD_HI 
             * @endcode
             * current options are [ PSD_HI | PSD_HO | PSD_dual | NonePSD ]
             * 
             * The class will read the appropriate section in the parameter file, i.e. the one with name \param psd_section_name ,
             * create an object of the desired type and return it.
             * 
             */
            static boost::shared_ptr<FuelCellShop::MicroScale::BasePSD<dim> > create_PSD (const std::string& psd_section_name, 
                                                                                                     ParameterHandler &param)
            {
                boost::shared_ptr<FuelCellShop::MicroScale::BasePSD<dim> > pointer;
                
                std::string concrete_name;
                
                param.enter_subsection("PSD parameters");
                {
                    param.enter_subsection("BasePSD");
                    {
                        concrete_name = param.get("psd type");
                    }
                    param.leave_subsection();
                }
                param.leave_subsection();
                
                typename FuelCellShop::MicroScale::BasePSD<dim>::_mapFactory::iterator iterator = FuelCellShop::MicroScale::BasePSD<dim>::get_mapFactory()->find(concrete_name);
                
                if (iterator != FuelCellShop::MicroScale::BasePSD<dim>::get_mapFactory()->end())
                {
                    if (iterator->second)
                    {
                        pointer = iterator->second->create_replica(psd_section_name);
                    }
                    else 
                    {
                        FcstUtilities::log<<"Pointer not initialized"<<std::endl;
                        abort();    
                    }
                }
                else
                {
                    FcstUtilities::log<<"Concrete name in FuelCellShop::MicroScale::BasePSD::create_psd does not exist"<<std::endl;
                    abort();
                }
                
                pointer->initialize(param);
                
                return pointer;
            }
            //@}
            
            ///@name Initalization
            //@{
            
            /**
             * This member function return the name of the type of layer, i.e. 
             * HIPSD
             * HOPSD
             * DualPSD
             * sss
             * 
             * Note that this is necessary if we want to find out not the name of the actual class which can be obtain using
             * @code const std::type_info& name = typeid(*this) @endcode
             * but the name of the parent class.
             * 
             * @note Do not implement this class anywhere other than the following "base" classes:
             * HIPSD
             * HOPSD
             * DualPSD
             */
            
            virtual const std::type_info& get_base_type() const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            }    
            
            
            /**
             * - \p por is the porosity for the given layer that contains this PSD 
             */
            inline void set_porosity(double porosity )
            {
                this->por = porosity;
            }
            
            /**
             * - \p por is the porosity for the given layer that contains this PSD 
             */
            inline double get_porosity() const
            {
                return this->por;
            }
            
            /**
             * Set the names of FCST solution variables
             * with respect to which you would like
             * to compute the derivatives of material properties.
             */
            inline void set_derivative_flags(const std::vector<VariableNames>& flags)
            {
                this->derivative_flags = flags;
            }
            
            /**
             * Set those solution variables which are constant in the particular application. If the effective properties in the psd depend on other variables that are usually
             * part of the solution vector but are assumed to be constant in this simulation, the const solution value should be passed to the class using this member function. This 
             * method should be called in the initialization section of the application. This function takes value to be set as the first argument and the #VariableNames as 
             * second argument. For instance, it's required to store constant temperature value for an isothermal application, in that case this method can be used. \em e.g., in 
             * order to set temperature as \p 353.0 [\p Kelvin] in the psd, you can use the following code:
             * @code
             * // In the initialization section of the application.
             * psd.set_constant_solution(353.0, VariableNames::temperature_of_REV);
             * @endcode
             */
            virtual void set_constant_solution(const double& value, const VariableNames& name)
            {
                constant_solutions[name] = value;
            }
            
            /**
             * If the effective properties in the psd depend on the solution, the solution for a given cell should be passed to
             * the class using this member function. It is used to set SolutionVariable structure inside the psd. This structure stores the solution variable values
             * at all quadrature points in the cell. For sample usage details, please see documentation of FuelCellShop::SolutionVariable structure.
             * 
             * Note, this function in the base psd sets the interface. It has to be reimplemented in respective child psd classes for respective uses.
             * 
             * \note Use only for solution variables.
             */
            virtual void set_solution(const std::vector< SolutionVariable >&)
            {
                const std::type_info& info = typeid(*this);
                deallog << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            }
            
            /**
             * Member function used to set the temperature [\p Kelvin] at every quadrature point
             * inside the cell. This function should particulary be used in the case of non-isothermal application.
             */
            virtual inline void set_temperature (const SolutionVariable& T_in)
            {
                const std::type_info& info = typeid(*this);
                deallog << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            }
            
            /**
             * Member function used to set the capillary pressure [\p psi] at every quadrature point
             * inside the cell. This function should particulary be used in the case of non-isothermal application.
             */
            virtual inline void set_capillary_pressure (const SolutionVariable& C_in)
            {
                const std::type_info& info = typeid(*this);
                deallog << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            }
            
            /**
             * Member function used to set the critical radius [\p nm] at every quadrature point
             * inside the cell. This function should particulary be used in the case of non-isothermal application.
             * It needs to be implemented at each iteration to reduce the computational time.
             */
            
            virtual inline void set_critical_radius()
            {             
                const std::type_info& info = typeid(*this);
                deallog << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
                
            }
            /**
             * This function is used to create PSD configuration plot by outputing all the numbers.
             */
            
            virtual inline void get_PSD_plot()
            {               
                const std::type_info& info = typeid(*this);
                deallog << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            } 
            
            /**
             * Member function used to set the saturation at every quadrature point
             * inside the cell. This function should particulary be used in the case of non-isothermal application.
             * It needs to be implemented at each iteration to reduce the computational time.
             */
            
            virtual inline void set_saturation()
            {               
                const std::type_info& info = typeid(*this);
                deallog << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            } 
            
            /**
             * Member function used to set the saturation at every quadrature point
             * inside the cell. This function should particulary be used in the case of non-isothermal application.
             * It needs to be implemented at each iteration to reduce the computational time.
             */
            
            virtual inline void initialize(ParameterHandler &param)
            {               
                const std::type_info& info = typeid(*this);
                deallog << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            } 
            //@}
            
            ///@name Accessors and info
            //@{
            /**
             * Return the name of the PSD.
             */
            inline const std::string& name_psd() const
            {
                return name;
            }
            
            /**
             * This function prints out
             * the psd properties.
             */
            virtual void print_psd_properties() const
            {
                const std::type_info& info = typeid(*this);
                deallog << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            } 
            
            /**
             * This function is used to compute the saturation by using PSD, the saturation is represented by \f$  \frac {V_water} { V_total}  \quad \f$  
             * 
             * The saturation of a porous medium can be obtained by integrating 
             * the change of cumulative pore volume fraction as a function of 
             * effective pore radius over the whole pore size domain
             * The saturation is represented by \f$  \frac {V_w} { V_t}  \quad \f$  
             * 
             * @warning This function is a virtual function and it must be reimplemented in the child class
             */
            
            
            virtual  void get_saturation(std::vector<double>& ) const = 0;
            virtual  void get_derivative_saturation(std::vector<double>& ) const = 0;
            /**
             * This function is used to compute the saturated_permeability by using PSD.
             * 
             * @warning This function is a virtual function and it must be reimplemented in the child class
             */
            virtual  void get_global_saturated_permeability(double& ) const = 0;
            
            /**
             * This function is used to compute the liquid_permeability by using PSD.
             * 
             * @warning This function is a virtual function and it must be reimplemented in the child class
             */
            virtual  void get_relative_liquid_permeability(std::vector<double>& ) const = 0;
            
                        /**
             * This function is used to compute the liquid_permeability by using PSD.
             * 
             * @warning This function is a virtual function and it must be reimplemented in the child class
             */
            virtual  void get_derivative_relative_liquid_permeability(std::vector<double>& ) const = 0;
            
            /**
             * This function is used to compute the gas_permeability by using PSD.
             * 
             * @warning This function is a virtual function and it must be reimplemented in the child class
             */
            virtual  void get_relative_gas_permeability(std::vector<double>& ) const = 0;
            
            /**
             * This function is used to compute the liquid_gas_interfacial_surface by using PSD.
             * 
             * @warning This function is a virtual function and it must be reimplemented in the child class
             */
            virtual  void get_liquid_gas_interfacial_surface(std::vector<double>& ) const = 0;
            
            /**
             * This function is used to compute the liquid_gas_interfacial_surface by using PSD.
             * 
             * @warning This function is a virtual function and it must be reimplemented in the child class
             */
            virtual  void get_derivative_liquid_gas_interfacial_surface(std::vector<double>& ) const = 0;
            
            /**
             * This function is used to compute the pore_wetted_wall by using PSD.
             * 
             * @warning This function is a virtual function and it must be reimplemented in the child class
             */
            virtual  void get_wetted_wall_surface_area(std::vector<double>& ) const = 0;
            
            /**
             * This function is used to compute the knudsen_radius by using PSD.
             * 
             * @warning This function is a virtual function and it must be reimplemented in the child class
             */
            virtual  void get_knudsen_radius(std::vector<double>& ) const = 0;
            
            //@}
            
        protected:
            ///@name Constructors, destructor, and parameter initalization
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
            BasePSD()
            {};

            /**
             * Constructor
             */
            BasePSD(const std::string& name);

            /**
             * Declare parameters for a parameter file. 
             * 
             * \note This member function must be called in all children using
             * FuelCellShop::MicroScale::BasePSD<dim>::declare_PSD_parameters(param).
             * 
             */            
            void virtual declare_parameters (ParameterHandler &param) const;
            
            /**
             * Member function used to read in data and initialize the necessary data
             * to compute the coefficients.
             */
            void _initialize (ParameterHandler &param) ;
            
            //@}
            ///@name Instance Delivery (Types)
            //@{         
            /** 
             * This object is used to store all objects of type psd. 
             */
            typedef std::map< std::string, BasePSD<dim>* > _mapFactory;      
            //@}
            
            ///@name Instance Delivery (Function)
            //@{         
            /**
             * Return the map library that stores all childrens of this class. The declare_parameters of
             * each one of the children that are in the map are called in declare_all_psd.
             * 
             * @warning In order for children of this class to appear in the map the following four things are 
             * necessary
             * - a static PROTOTYPE object has to be created. For example, 
             * in the .h file:
             * @code
             * static PSD_dual<dim> const* PROTOTYPE;
             * @endcode
             * in the .cc file:
             * @code
             * template <int dim>
             * NAME::PSD_dual<dim> const* NAME::PSD_dual<dim>::PROTOTYPE = new NAME::PSD_dual<dim>();
             * @endcode
             * - a default constructor which creates the PROTOTYPE is needed, e.g.
             * @code
             * template <int dim>
             * NAME::PSD_dual<dim>::PSD_dual()
             * : FuelCellShop::psd::BasePSD<dim> ()
             * {
             *    deallog<<" Register PSD_dual to FactoryMap"<<std::endl;
             *    this->get_mapFactory()->insert(std::pair<std::string, FuelCellShop::psd::BasePSD<dim>* >(concrete_name, this));
             * }
             * @endcode
             * - a static std::string concrete_name must be declared and initialized to the desired name of
             * the children class. This name is used as the name in the map and as the name of the subsection 
             * where its parameters are declared. For example, 
             * in the .h file:
             * @code
             * static const std::string concrete_name;
             * @endcode
             * in the .cc file:
             * @code
             * template <int dim>
             * const std::string NAME::PSD_dual<dim>::concrete_name ("PSD_dual");
             * @endcode
             * - virtual boost::shared_ptr<FuelCellShop::psd::BasePSD<dim> > create_replica () needs to
             * be re-implemented in the child. For example, in the .h file
             * @code
             * virtual boost::shared_ptr<FuelCellShop::psd::BasePSD<dim> > create_replica (std::string &name)
             * {
             *   return boost::shared_ptr<FuelCellShop::psd::PSD_dual<dim> > (new FuelCellShop::psd::PSD_dual<dim> (name));
             * }   
             * @endcode
             */
            static _mapFactory * get_mapFactory()
            {
                static _mapFactory mapFactory;
                return &mapFactory;
            }                               
            ///@name Instance Delivery (Private functions)
            //@{                     
            /**
             * This member function is used to create an object of type gas diffusion layer
             * 
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::MicroScale::BasePSD<dim> > create_replica (const std::string &name)
            {
                const std::type_info& info = typeid(*this);
                deallog << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            //@}
            
            ///@name Internal variables
            //@{            
            /** Name of the psd. This value is used as a header in a subsection for the parameter file
             * where all the information concerning this psd is held.*/
            
            const std::string name;
            
            /** Flags for derivatives: These flags are used to request derivatives of material properties.*/
            std::vector<VariableNames> derivative_flags;
            
            /**
             * Map storing values of solution variables constant in a particular application. 
             */
            std::map< VariableNames, double > constant_solutions;

            /**
             * The gamma is the Water-air interface surface tension
             */
            double gamma;
            
            /**
             * The contact_angle is the contact angle between air and water
             */
            double contact_angle;
            
            /**
             * The lamda is measured from the mercury intrusion experiment
             */
            double lamda;
            
            /**
             * The F_HI is the fraction of the total volume corresponding to the hydrophilic pores
             */
            double F_HI;
            
            /**
             * The F_HO is the fraction of the total volume corresponding to the hydrophobic pores
             */
            double F_HO;
            
            /**
             * The f_k is the contribution of the log-normal distribution k to the total PSD
             */            
            std::vector<double> f_k;
            
            /**
             * The r_k is the characteristic pore size of the distribution k
             */            
            std::vector<double> r_k;            
            
            /**
             * The s_k is the spread of the distribution k
             */
            std::vector<double> s_k;
            
            /**
             * The por is the porosity inherited from the layer class
             */
            double por;
            //@}
        };
        
    } // PSD
    
}  // FuelCellShop

#endif