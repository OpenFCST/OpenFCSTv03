//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2011-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: catalyst_base.h
//    - Description: Base catalyst material class
//    - Developers: M. Secanell and Madhur Bhaiya
//    - $Id: catalyst_base.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP_CATALYST_BASE__H
#define _FUELCELLSHOP_CATALYST_BASE__H

// Include deal.II classes
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/fe/fe_values.h>

// Include FCST classes
#include <materials/base_material.h>
#include <utils/fcst_constants.h>


namespace FuelCellShop
{
    namespace Material
    {
        /**
         * 
         * This class implements the interface to compute the properties of a "standard" catalyst.
         * 
         * Note that there are many types of catalyst such as platinum, ruthinium, core-shell Pt alloys and 
         * non-precious metal alloys. This class serves a common interface for each one of these types. 
         * It implements member functions to declare the parameter file for each catalyst, read the parameter
         * file and select the appropriate catalysts based on the user input.
         * 
         * The properties of each catalyst are implemented in the children classes. The most important properties 
         * that define our catalyst are the activity and reaction order parameters. Therefore, before the layer 
         * can be used, the reaction that it is catalyzing should be specified.
         * 
         * After the reaction is known, the class return the exchange current density, transfer coefficient, reaction 
         * order etc. for the type of catalysts. If the reaction is a multi-step reaction, it might return additional
         * parameters such as the free-energies of activation for each elementary reaction.
         * 
         * This class should never be used in an application. Instead, the class is used inside FuelCellShop::Kinetics::BaseKinetics objects. To use
         * the class, follow the instructions below.
         * 
         * 
         * <h3>Usage Details:</h3>  
         * 
         * 
         * 
         * 
         * \author M. Bhaiya, 2012-13
         * \author M. Secanell, 2013 
         */
        class CatalystBase 
        :
        public BaseMaterial
        {
        public:
            ///@name Instance Delivery (Public functions)
            //@{
            /**
             * Function used to declare all the data necessary in the parameter files for
             * all BaseKinetics children.
             * 
             */
            static void declare_Catalyst_parameters (ParameterHandler &param)
            {
                
                for (typename FuelCellShop::Material::CatalystBase::_mapFactory::iterator iterator = FuelCellShop::Material::CatalystBase::get_mapFactory()->begin(); 
                     iterator != FuelCellShop::Material::CatalystBase::get_mapFactory()->end(); 
                iterator++)
                     {
                         iterator->second->declare_parameters(param);
                     }        
            }

             /**
              * 
              * Function called in create_CatalystLayer and used to select the appropriate CatalystBase children that will be used 
              * in the layer. The name of the CatalystBase children object to be used is provided in catalyst_name.
              *
              * The name of the CatalystBase children object is provided in the ParameterHandler in the CatalystLayer subsection as follows:
              * 
              * @code 
              * subsection Catalyst Layer Properties <- This name is the name of the catalyst layer subsection where the kinetics are taking place.
              * (...)
              *   set Catalyst type = Platinum
              * (...)
              * end
              * @endcode
              * current options are [ Platinum ]
              * 
              * 
              */
            static boost::shared_ptr<FuelCellShop::Material::CatalystBase > create_Catalyst (ParameterHandler &param,
                                                                                              std::string catalyst_name)
             {       
                 boost::shared_ptr<FuelCellShop::Material::CatalystBase > pointer;

                 typename FuelCellShop::Material::CatalystBase::_mapFactory::iterator iterator = FuelCellShop::Material::CatalystBase::get_mapFactory()->find(catalyst_name);
                 
                 if (iterator != FuelCellShop::Material::CatalystBase::get_mapFactory()->end())
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
                     FcstUtilities::log<<"Concrete name in FuelCellShop::Material::CatalystBase::create_Catalyst does not exist"<<std::endl;
                     abort();
                 }
                 
                 pointer->initialize(param);
                 
                 return pointer;
             }
            //@}
            ///@name Initializaton
            //@{            
            /**
             * Member function used to specify the reaction for which the kinetic parameters are needed.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void set_reaction_kinetics(const ReactionNames)
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            };
            //@}            
            ///@name Accessors and Info
            //@{
            /**
             * Return anodic transfer coefficient for the reaction specified using #set_reaction_kinetics method.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void alpha_anodic(double&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            };
            
            /**
             * Return derivative of anodic transfer coefficient for the reaction specified using #set_reaction_kinetics method
             * with respect to the solution and design parameters specified using #set_derivative_flags method.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void derivative_alpha_anodic(std::vector<double>& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            };
            
            /**
             * Return cathodic transfer coefficient for the reaction specified using #set_reaction_kinetics method.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void alpha_cathodic(double&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            };
            
            /**
             * Return derivative of cathodic transfer coefficient for the reaction specified using #set_reaction_kinetics method
             * with respect to the solution and design parameters specified using #set_derivative_flags method.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void derivative_alpha_cathodic(std::vector<double>& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            };
            
            /**
             * Compute the exchange current density [\p A/cm^2] for the reaction specified using #set_reaction_kinetics method. It takes temperature [\p Kelvin]
             * as an input argument by reference and returns the exchange current density [\p A/cm^2].
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual double exchange_current_density(const double&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            };
            /**
             * Compute the derivative of exchange current density [\p A/cm^2] w.r.t temperature [\p Kelvin] for the reaction specified using #set_reaction_kinetics method. It takes
             * temperature [\p Kelvin] as an input argument by reference and returns the derivative.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual double derivative_exchange_current_density(const double&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            };
                  
            /**
             * Compute the reference concentration for the reaction specified using #set_reaction_kinetics method.
             * The reaction might depend on more than one species, therefore this class returns a std::map of reference concentrations referenced using #VariableNames 
             * as \p Key. This map is returned by reference (second argument). It takes a vector of #VariableNames as first argument, corresponding 
             * to those solution variables for whose reference concentrations are required.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void reference_concentration(const std::vector<VariableNames>&,
                                                 std::map<VariableNames, double>& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            };
            
            /** Compute the theroretical cell voltage [\p Volts] as a function of temperature, for the reaction specified using #set_reaction_kinetics
             * method. It takes temperature [\p Kelvin] as an input argument and returns the voltage.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual double voltage_cell_th(const double&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            };
            
            /** Compute the derivative of theoretical cell voltage [\p Volts] w.r.t temperature [\p Kelvin], for the reaction specified using 
             * #set_reaction_kinetics method. It takes temperature [\p Kelvin] as an input argument and returns the derivative.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual double dvoltage_cell_th_dT(const double&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            };
            
            /**
             * Compute the reaction order of the electrochemical reaction with respect to each species involved in the reaction specified
             * using #set_reaction_kinetics method. The reaction might depend on more than one species, therefore this class returns a std::map of reaction orders referenced using #VariableNames 
             * as \p Key. This map is returned by reference (second argument). It takes a vector of #VariableNames as first argument, corresponding 
             * to those solution variables for whose reaction orders are required.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual void reaction_order(const std::vector<VariableNames>&,
                                        std::map<VariableNames, double>&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            };
                       
            /**
             * Member function to return the name of the reaction set in the catalyst material class.
             */
            inline ReactionNames get_reaction_name() const
            {
                Assert(name_reaction_kinetics != noReaction, ExcMessage("Reaction name not yet set in the CatalystBase object."));
                return name_reaction_kinetics;
            }
            
            /**
             * Member function to return method for kinetics parameters (ORR).
             */
            inline std::string get_kinetic_parameter_method() const
            {
                Assert(method_kinetics_ORR.size() != 0, ExcMessage("Kinetic parameter method not yet set in the CatalystBase object."));
                return method_kinetics_ORR;
            }


            /** Obtain the density [\p gm/cm^3].*/
            inline double get_density() const
            { 
                return density; 
            };
            //@}
              
        protected:
            ///@name Constructors, destructor, and parameter initalization
            //@{
            /**
             * Constructor
             * 
             * \deprecated
             */
            CatalystBase(std::string name)
            : FuelCellShop::Material::BaseMaterial(name)
            {};

            /**
             * Constructor
             */
            CatalystBase()
            : FuelCellShop::Material::BaseMaterial()
            {};

            /**
             * Destructor
             */
            virtual ~CatalystBase()
            {};  
            
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
             * This object is used to store all objects of type CatalystBase. 
             */
            typedef std::map< std::string, FuelCellShop::Material::CatalystBase* > _mapFactory;      
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
             * This member function is used to create an object of type CatalystBase
             * 
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::Material::CatalystBase > create_replica ()
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            //@}
            ///@name Protected member functions:
            //@{                        
            /**
             * Check that whether a particular reaction is implemented in the class or not. It takes a std::string as an input argument, against
             * which the check is done. 
             * <b>For Developers:</b> As more reactions are implemented in the class, this method should be extended. Moreover, it is recommended to use this
             * method in #set_reaction_kinetics method as an assertation check.
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual bool check_reaction_implementation(const std::string) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            };
            //@}
            
            
            ///@name Internal variables
            //@{            
            /**
             * Reaction name for which the class returns kinetic parameters.
             */
            ReactionNames name_reaction_kinetics;
            
            /** Method for kinetics parameters (ORR), given in the parameter file. */
            std::string method_kinetics_ORR;
            
            /** Density of catalyst particles [\p gm/cm^3].*/
            double density;
            //@}
        };
        
    }
}


#endif
