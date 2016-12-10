//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: polymer_electrolyte_material_base.h
//    - Description: base class for polymer electrolyte materials e.g. Nafion
//    - Developers: M. Secanell and Madhur Bhaiya
//    - Id: $Id: polymer_electrolyte_material_base.h 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP_POLYMER_ELECTROLYTE_MATERIAL_BASE__H
#define _FUELCELLSHOP_POLYMER_ELECTROLYTE_MATERIAL_BASE__H

// Include deal.II classes
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/fe/fe_values.h>

// Include FCST classes
#include <materials/base_material.h>
#include <materials/PureGas.h>
#include <utils/fcst_utilities.h>



namespace FuelCellShop
{
    namespace Material
    {
        /**
         * 
         * This class implements the interface to compute the properties of a "standard" polymer electrolyte membrane material.
         * 
         * Note that there are many types of polymer electrolye materials such as Nafion, Gore , Flemion membranes etc.. This class serves a common interface for each one of these types. 
         * It implements member functions to declare the parameter file for each catalyst, read the parameter file and select the appropriate
         * polymer electrolyte materials based on the user input. It basically has accessor methods for various physical
         * properties and return bulk transport propeties of the materials. These methods are implemented in the children classes. Since, almost all of 
         * the transport properties are dependent on solution variables, \em e.g. #membrane_water_content, \f$ \lambda \f$, 
         * #temperature_of_REV, \f$ T \f$ etc. Hence it's important to set these solution variables values inside the class, prior to accessing
         * the bulk transport properties, \em e.g. proton conductivity, \f$ \sigma_{H^+} \f$ [\p S/cm], water diffusivity, \f$ D_{\lambda} \f$ [\p cm^2/s], etc.
         * Besides, the derivatives of these transport properties w.r.t various solution variables can also be returned from this class.
         * 
         * <h3>Usage Details:</h3>  
         * 
         * 
         * 
         * 
         * \author M. Bhaiya, 2012-13
         * \author M. Secanell, 2013
         */
        class PolymerElectrolyteBase 
        :
        public BaseMaterial
        {
        public:
            ///@name Instance Delivery (Public functions)
            //@{
            /**
             * Function used to declare all the data necessary in the parameter files for
             * all PolymerElectrolyteBase children.
             * 
             */
            static void declare_PolymerElectrolyte_parameters (ParameterHandler &param)
            {
                
                for (typename FuelCellShop::Material::PolymerElectrolyteBase::_mapFactory::iterator iterator = FuelCellShop::Material::PolymerElectrolyteBase::get_mapFactory()->begin(); 
                     iterator != FuelCellShop::Material::PolymerElectrolyteBase::get_mapFactory()->end(); 
                iterator++)
                     {
                         iterator->second->declare_parameters(param);
                     }        
            }
                        
            /**
              * 
              * Function called in create_CatalystLayer and used to select the appropriate PolymerElectrolyteBase children that will be used 
              * in the layer. The name of the PolymerElectrolyteBase children object to be used is provided in polymer_electrolyte_name.
              *
              * The name of the PolymerElectrolyteBase children object is provided in the ParameterHandler in the CatalystLayer subsection as follows:
              * 
              * @code 
              * subsection Catalyst Layer Properties <- This name is the name of the catalyst layer subsection where the kinetics are taking place.
              * (...)
              *   set Polymer Electrolyte type = Nafion
              * (...)
              * end
              * @endcode
              * current options are [ Nafion ]
              * 
              * 
              */
             static boost::shared_ptr<FuelCellShop::Material::PolymerElectrolyteBase > create_PolymerElectrolyte (ParameterHandler &param,
                                                                                                                  std::string polymer_electrolyte_name)
             {       
                 boost::shared_ptr<FuelCellShop::Material::PolymerElectrolyteBase > pointer;

                 typename FuelCellShop::Material::PolymerElectrolyteBase::_mapFactory::iterator iterator = FuelCellShop::Material::PolymerElectrolyteBase::get_mapFactory()->find(polymer_electrolyte_name);
                 
                 if (iterator != FuelCellShop::Material::PolymerElectrolyteBase::get_mapFactory()->end())
                 {
                     if (iterator->second)
                     {
                         pointer = iterator->second->create_replica();
                     }
                     else 
                     {
                         FcstUtilities::log<<"Pointer not initialized"<<std::endl;
                         abort();
                     }
                 }
                 else
                 {
                     FcstUtilities::log<<"Concrete name in FuelCellShop::Material::PolymerElectrolyteBase::create_PolymerElectrolyte does not exist"<<std::endl;
                     abort();
                 }
                 
                 pointer->initialize(param);
                 
                 return pointer;
             }
            //@}
  
            ///@name Transport properties and derivatives accessor methods
            //@{
            
            /**
             * Compute the equilibrium water content, \f$ \lambda_{eq} \f$, inside the polymer electrolyte for vapor-equilibriated case, at
             * every quadrature point in the cell. It takes a vector as an input argument and values are passed by reference.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void sorption_isotherm(std::vector<double>&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            }
            /**
             * Compute the derivatives for water sorption source terms, \f$ \frac{\partial \lambda_{eq}}{\partial u} \f$, at every quadrature 
             * point in the cell. The derivatives are computed based on the flags set by the #set_derivative_flags method. It takes map as an input argument by reference, in 
             * which \p Key corresponds to the variable about which derivative is being computed, \em i.e., #VariableNames and \p Value corresponds to
             * the vector storing the derivative values.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void sorption_isotherm_derivative(std::map < VariableNames, std::vector<double> >&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            
            /**
             * Compute the proton conductivity, \f$ \sigma_{H^+} \f$ [\p S/cm], inside the polymer electrolyte for constant case. It takes a double as an 
             * input argument and value is passed by reference.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void proton_conductivity(double&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            /**
             * Compute the proton conductivity, \f$ \sigma_{H^+} \f$ [\p S/cm], inside the polymer electrolyte, at every quadrature point in the cell. It takes vector as an 
             * input argument and values are passed by reference.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void proton_conductivity(std::vector<double>&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            /**
             * Compute the derivatives of proton conductivity at every quadrature point in the cell. The derivatives are computed based on the flags
             * set by the #set_derivative_flags method. It takes map as an input argument by reference, in  which \p Key corresponds to the variable about 
             * which derivative is being computed, \em i.e., #VariableNames and \p Value corresponds to the vector storing the derivative values.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void proton_conductivity_derivative(std::map< VariableNames, std::vector<double> >&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            
            /**
             * Compute the water diffusivity, \f$ D_{\lambda} \f$ [\p cm^2/s], inside the polymer electrolyte for constant case. It takes a double as an 
             * input argument and value is passed by reference.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void water_diffusivity(double&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            /**
             * Compute the water diffusivity, \f$ D_{\lambda} \f$ [\p cm^2/s], inside the polymer electrolyte at every quadrature point in the cell. It takes vector as an 
             * input argument and values are passed by reference.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void water_diffusivity(std::vector<double>&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            /**
             * Compute the derivatives of water diffusivity at every quadrature point in the cell. The derivatives are computed based on the flags
             * set by the #set_derivative_flags method. It takes map as an input argument by reference, in  which \p Key corresponds to the variable about 
             * which derivative is being computed, \em i.e., #VariableNames and \p Value corresponds to the vector storing the derivative values.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void water_diffusivity_derivative(std::map< VariableNames, std::vector<double> >&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            
            /**
             * Compute the electro-osmotic drag coefficient inside the polymer electrolyte at every quadrature point in
             * the cell. It takes vector as an input argument and values are passed by reference.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void electroosmotic_drag(std::vector<double>&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            /**
             * Compute the derivatives of electro-osmotic drag coefficient, at every quadrature point in the cell. The derivatives are computed based on the flags
             * set by the #set_derivative_flags method. It takes map as an input argument by reference, in  which \p Key corresponds to the variable about 
             * which derivative is being computed, \em i.e., #VariableNames and \p Value corresponds to the vector storing the derivative values.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void electroosmotic_drag_derivative(std::map< VariableNames, std::vector<double> >&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            
            /**
             * Compute the thermo-osmotic diffusion coefficient, [\p gm/ \p (cm-s-K \p )],  inside the polymer electrolyte at every quadrature 
             * point in the cell. It takes vector as an input argument and values are passed by reference.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void thermoosmotic_coeff(std::vector<double>&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            /**
             * Compute the derivatives of thermo-osmotic diffusion coefficient, at every quadrature point in the cell. The derivatives are computed based on the flags
             * set by the #set_derivative_flags method. It takes map as an input argument by reference, in  which \p Key corresponds to the variable about 
             * which derivative is being computed, \em i.e., #VariableNames and \p Value corresponds to the vector storing the derivative values.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void thermoosmotic_coeff_derivative(std::map< VariableNames, std::vector<double> >&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;  
            }
            
            /**
             * Compute the oxygen diffusivity, [\p cm^2/s], inside the polymer electrolyte for constant case. It takes a double as an 
             * input argument and value is passed by reference.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void oxygen_diffusivity(double&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            /**
             * Compute the oxygen diffusivity, [\p cm^2/s], inside the polymer electrolyte for constant case. It takes a double as an
             * input argument and value is passed by reference.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void hydrogen_diffusivity(double&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            /**
             * Compute the oxygen diffusivity, [\p cm^2/s], inside the polymer electrolyte at every quadrature point in the cell. It takes
             * vector as an input argument and values are passed by reference.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void oxygen_diffusivity(std::vector<double>&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            /**
             * Compute the derivatives of oxygen diffusivity, at every quadrature point in the cell. The derivatives are computed based on the flags
             * set by the #set_derivative_flags method. It takes map as an input argument by reference, in  which \p Key corresponds to the variable about 
             * which derivative is being computed, \em i.e., #VariableNames and \p Value corresponds to the vector storing the derivative values.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void oxygen_diffusivity_derivative(std::map< VariableNames, std::vector<double> >&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            
            /**
             * Compute the proton diffusivity, [\p cm^2/s], inside the polymer electrolyte for constant case. It takes a double as an 
             * input argument and value is passed by reference.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void proton_diffusivity(double&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }

            /**
             * Compute the enthalpy of sorption [\p J/mol] of water, at all quadrature points in the cell. It takes
             * vector as an input argument and values are passed by reference.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void sorption_enthalpy(std::vector<double>& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            /**
             * Compute the derivatives of enthalpy of sorption of water, at every quadrature point in the cell. The derivatives are computed based on the flags
             * set by the #set_derivative_flags method. It takes map as an input argument by reference, in  which \p Key corresponds to the variable about 
             * which derivative is being computed, \em i.e., #VariableNames and \p Value corresponds to the vector storing the derivative values.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void sorption_enthalpy_derivative(std::map< VariableNames, std::vector<double> >&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }

            //@}
            
            ///@name Accessor methods for molar enthalpy of sorbed water and its derivatives
            //@{
            
            /**
             * Compute the molar enthalpy, \f$ \bar{H}_{\lambda} ~\f$ [\p J/mol] of sorbed water in the polymer electrolyte as a function of \f$ T \f$. It takes \b Temperature, \f$ T \f$ as 
             * input argument by reference and returns the double value corresponding to molar enthalpy.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual double get_Hlambda(const double&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            
            /**
             * Compute \f$ \frac{\partial  \bar{H}_{\lambda}}{\partial T} \f$ of sorbed water in the polymer electrolyte as a function of \f$ T \f$. It takes \b Temperature, \f$ T \f$ as
             * input argument by reference and returns the double value corresponding to the derivative.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual double get_dHlambda_dT(const double&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            
            /**
             * Compute \f$ \frac{\partial^2  \bar{H}_{\lambda}}{\partial T^2} \f$ of sorbed water in the polymer electrolyte as a function of \f$ T \f$. It takes \b Temperature, \f$ T \f$ as 
             * input argument by reference and returns the double value corresponding to the second derivative.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual double get_d2Hlambda_dT2(const double&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            //@}
            
            ///@name Physical properties accessor methods
            //@{
            
            /** Get the density [\p gm/cm^3] of the dry polymer electrolyte material. */
            inline double get_density() const
            {
                return rho_M;
            }
            
            /** Get Henry's constant [\p Pa-cm^3/mol] for oxygen dissolution in the polymer electrolyte. */
            inline double get_H_O2() const
            {
                return H_O2;
            }
            
            /** Get Henry's constant [\p Pa-cm^3/mol] for hydrogen dissolution in the polymer electrolyte. */
            inline double get_H_H2() const
            {
                return H_H2;
            }
            
            /** Get Equivalent Weight (grams of dry polymer electrolyte per moles of \f$ SO_3^- \f$) of the polymer electrolyte material.*/
            inline double get_EW() const
            {
                return EW;
            }
            
            /** Get permittivity of the polymer electrolyte material. */
            inline double get_permittivity() const
            {
                return permittivity;
            }
            
            //@}
            
            ///@name Solution setting methods
            //@{
            
            /**
             * Specify the total pressure in Pascals.
             */
            inline void set_p_t(const double& p_t)
            {
                p_total = p_t;
            }
            
            /**
             * Specify the temperature [\p Kelvin] in the polymer electrolyte for isothermal case.
             */
            inline void set_T (const double& Temp)
            {
                T = Temp;
            }
            
            /**
             * Specify the water content in the polymer electrolyte for constant lambda case.
             */
            inline void set_lambda(const double& L)
            {
                lambda = L;
            }
            
            /**
             * Set the solution variable, membrane water content \f$ \lambda \f$. It also initializes the temperature solution variable, if it
             * isn't initialized yet. This will help in the case of isothermal applications. #set_T method should be used in the initialization of the
             * application to ensure that temperature value is set for isothermal case, otherwise default value \p 0.0 will be used leading to wrong results.
             */
            inline void set_membrane_water_content(const FuelCellShop::SolutionVariable& l_in)
            {
                Assert( ((l_in.is_initialized()) && (l_in.get_variablename() == membrane_water_content)),
                        ExcMessage("Input solution variable not initialized properly in PolymerElectrolyteBase::set_membrane_water_content method.") );
                
                lambda_var = l_in;
                
                if ( !T_var.is_initialized() )
                    T_var = FuelCellShop::SolutionVariable(T, lambda_var.size(), temperature_of_REV);
            }
            
            /**
             * Set the solution variable, temperature \f$ T \f$ [\p Kelvin]. It also initializes the lambda solution variable, it it isn't 
             * initialized yet. This will help in the case of applications where \f$ \lambda \f$ is not being solved for. Membrane water content value is 
             * defaulted to \p 12.0, in case if different value is required, use #set_lambda method.
             */
            inline void set_temperature(const FuelCellShop::SolutionVariable& T_in)
            {
                Assert( ((T_in.is_initialized()) && (T_in.get_variablename() == temperature_of_REV)),
                        ExcMessage("Input solution variable not initialized properly in PolymerElectrolyteBase::set_temperature method.") );
                
                T_var = T_in;
                
                if ( !lambda_var.is_initialized() )
                    lambda_var = FuelCellShop::SolutionVariable(lambda, T_var.size(), membrane_water_content);
            }
            
            /**
             * Set the solution variable, water vapor molar fraction \f$ x_{H_2O} \f$. It also initializes the temperature solution variable, if it
             * isn't initialized yet. This will help in the case of isothermal applications. #set_T method should be used in the initialization of the
             * application to ensure that temperature value is set for isothermal case, otherwise default value \p 0.0 will be used leading to 
             * wrong results. Also, total pressure should be set using #set_p_t method.
             */
            inline void set_water_molar_fraction(const FuelCellShop::SolutionVariable& x_in)
            {
                Assert( ((x_in.is_initialized()) && (x_in.get_variablename() == water_molar_fraction)),
                        ExcMessage("Input solution variable not initialized properly in PolymerElectrolyteBase::set_water_molar_fraction method.") );
                
                xwater_var = x_in;
                
                if ( !T_var.is_initialized() )
                    T_var = FuelCellShop::SolutionVariable(T, xwater_var.size(), temperature_of_REV);
            }
            //@}
            
        protected:
            ///@name Constructors, destructor and parameter initialization
            //@{
            /**
             * Constructor
             */
            PolymerElectrolyteBase(std::string name)
            : BaseMaterial(name)
            {
                lambda = 12.0;
                T = 0.0;        // Used for assertion to throw error, if temperature values are not set.
                p_total = 0.0;  // Used for assertion to throw error, if pressure values are not set.
                
                water_mat = new FuelCellShop::Material::WaterVapor;
            };
            
            /**
             * Constructor
             */
            PolymerElectrolyteBase()
            : FuelCellShop::Material::BaseMaterial()
            {
                lambda = 12.0;
                T = 0.0;        // Used for assertion to throw error, if temperature values are not set.
                p_total = 0.0;  // Used for assertion to throw error, if pressure values are not set.
                
                water_mat = new FuelCellShop::Material::WaterVapor;
            };
            
            /**
             * Destructor
             */
            virtual ~PolymerElectrolyteBase()
            {
                delete water_mat;
            };
            
            /**
             * Declare parameters for a parameter file
             * \warning This is a PureFunction and it does not declare anything, so please do not call this function in the children.
             */
            virtual void declare_parameters(ParameterHandler &param) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            };
                       
            /** Initialize parameters.
             * \warning This is a PureFunction and it does not declare anything, so please do not call this function in the children.
             */
            virtual void initialize (ParameterHandler &param)
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            };
            //@}

             ///@name Instance Delivery (Types)
            //@{         
            /** 
             * This object is used to store all objects of type PolymerElectrolyteBase. 
             */
            typedef std::map< std::string, FuelCellShop::Material::PolymerElectrolyteBase* > _mapFactory;      
            //@}       
            
            ///@name Instance Delivery (Private and static)
            //@{         
            /**
             * 
             */
            static _mapFactory * get_mapFactory()
            {
                static _mapFactory mapFactory;
                return &mapFactory;
            }  
            //@}
            
            ///@name Instance Delivery (Private functions)
            //@{         
            /**
             * This member function is used to create an object of type PolymerElectrolyteBase.
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::Material::PolymerElectrolyteBase > create_replica ()
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            //@}
            
            ///@name Numerical constants
            //@{
            /** Equivalent weight.*/
            double EW;
            
            /** Dry polymer electrolyte density [\p gm/cm^3].*/
            double rho_M;
            
            /**
             * Henry's Constant [\p Pa-cm^3/mol] for dissolution of oxygen.
             */
            double H_O2;
            /**
             * Henry's Constant [\p Pa-cm^3/mol] for dissolution of hydroge.
             */
            double H_H2;
            /**
             * Permittivity of the material.
             */
            double permittivity;
            
            
            /** Temperature [\p Kelvins] for isothermal case. In order to update this value, use #set_T method. */
            double T;
            
            /** Membrane water content, \f$ \lambda \f$ for constant lambda case. In order to update
             * this value, use #set_lambda method.
             */
            double lambda;
            
            /** Total pressure (in Pascals). In order to update this value, use #set_p_t method. */
            double p_total;
            
            
            /**
             * Proton conductivity value [\p S/cm], given in the parameter file; used with "Constant" method.
             */
            double sigma_p;
            
            /**
             * Diffusion of water in polymer electrolyte (if constant option is used), [\p cm2/s].
             */
            double diffusion_w;
            
            /**
             * Given electroosmotic drag value in the parameter file. This
             * value will be used when method_electroosmotic_drag is set to "Constant".
             */
            double given_n_drag;
            
            /**
             * Given thermo-osmotic diffusion coefficient [\p gm/\p(cm-s-K \p)] in the parameter file. This
             * value will be used when method_thermoosmosis is set to "Constant".
             */
            double given_thermoosmotic_coeff;
            
            /**
             * Given enthalpy of sorption of water [\p J/mol] in the parameter file. This
             * value will be used when method_enthalpy_sorption is set to "Constant".
             */
            double given_enthalpy_sorption;
            
            /**
             * Diffusion coefficient of oxygen [\p cm^2/s], given in the parameter file.
             */
            double D_O2;
            
            /**
             * Effective diffusion coefficient of hydrogen [\p cm^2/s], given in the parameter file.
             */
            double D_H2;

            /**
             * Diffusion coefficient of protons [\p cm^2/s], given in the parameter file.
             */
            double D_Protons;
            //@}
            
            ///@name Empirical method name strings
            //@{
            
            /**
             * Method/Semi-empirical relation to compute protonic conductivity.
             */
            std::string method_conductivity;

            /**
             * Method to compute equilibrium water content value from sorption isotherm.
             */
            std::string method_sorption;
            
            /**
             * Method/semi-empirical relation to compute water diffusivity
             */
            std::string method_diffusivity;
            
            /**
             * Method/semi-empirical relation to compute Electro-osmotic drag
             */
            std::string method_electroosmotic_drag;
            
            /**
             * Method/semi-empirical relation to compute thermo-osmotic diffusion coefficient.
             */
            std::string method_thermoosmosis;
            
            /**
             * Method/semi-empirical relation to compute enthalpy of sorption of water.
             */
            std::string method_enthalpy_sorption;
            //@}
            
            ///@name Solution variables
            //@{
            
            /**
             * Solution variable, membrane water content \f$ \lambda \f$.
             */
            FuelCellShop::SolutionVariable lambda_var;
            
            /**
             * Solution variable, temperature \f$ T \f$.
             */
            FuelCellShop::SolutionVariable T_var;
            
            /**
             * Solution variable, water vapor molar fraction \f$ x_{H_2O} \f$.
             */
            FuelCellShop::SolutionVariable xwater_var;
            //@}


            /**
             * Pointer to FuelCellShop::Material::WaterVapor object, dynamically allocated when the class is constructed.
             */
            FuelCellShop::Material::WaterVapor* water_mat;    
        };
    }
}

#endif
