//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2011-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: base_kinetics.h
//    - Description: Implements the properties of a "standard" carbon black support
//    - Developers: M. Secanell and Madhur Bhaiya
//    - $Id: catalyst_support_base.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__CATALYST_SUPPORT_BASE_H
#define _FUELCELLSHOP__CATALYST_SUPPORT_BASE_H

// Include FCST classes
#include <materials/base_material.h>


namespace FuelCellShop
{
    namespace Material
    {
        /**
         * 
         * This class implements the interface to compute the properties of a "standard" catalyst support.
         * 
         * Note that there are many types of catalyst supports such as carbon black, carbon nanotubes and 3M organic
         * nano-wiskers. This class serves a common interface for each one of these types. Children should be implemented
         * for each type with appropriate propertes. The properties of each catalyst support materials are implemented in the children classes. This class is used
         * to return physical properties such as density, electrical conductivity and thermal conductivity. This class is normally not required to be used inside
         * the application. Instead, it's used inside the catalyst layer classes to access these above-mentioned physical properties.
         *
         * <h3>Usage Details:</h3>
         * Three accessors methods are available, \em viz., #get_density, #get_electrical_conductivity, #get_thermal_conductivity, which return \p double values
         * corresponding to density [\p gm/cm^3], electrical conductivity [\p S/cm] and thermal conductivity [\p W/\p (cm-K \p )] respectively, of the catalyst support material.
         * 
         * 
         * \author M. Secanell, 2011-13
         * \author Madhur Bhaiya, 2013
         * 
         */
        class CatalystSupportBase
        :
        public BaseMaterial
        {
        public:
            ///@name Instance Delivery (Public functions)
            /**
             * Function used to declare all the data necessary in the parameter files for
             * all CatalystSupportBase children.
             * 
             */
            static void declare_CatalystSupport_parameters (ParameterHandler &param)
            {
                
                for (typename FuelCellShop::Material::CatalystSupportBase::_mapFactory::iterator iterator = FuelCellShop::Material::CatalystSupportBase::get_mapFactory()->begin(); 
                     iterator != FuelCellShop::Material::CatalystSupportBase::get_mapFactory()->end();
                    iterator++)
                     {
                         iterator->second->declare_parameters(param);
                     }        
            }

             /**
              * 
              * Function called in create_CatalystLayer and used to select the appropriate CatalystSupportBase children that will be used 
              * in the layer. The name of the CatalystSupportBase children object to be used is provided in support_name.
              *
              * The name of the CatalystSupportBase children object is provided in the ParameterHandler in the CatalystLayer subsection as follows:
              * 
              * @code 
              * subsection Catalyst Layer Properties <- This name is the name of the catalyst layer subsection where the kinetics are taking place.
              * (...)
              *   set Catalyst support type = CarbonBlack
              * (...)
              * end
              * @endcode
              * current options are [ CarbonBlack ]
              * 
              * 
              */
             static boost::shared_ptr<FuelCellShop::Material::CatalystSupportBase > create_CatalystSupport (ParameterHandler &param,
                                                                                                            std::string support_name)
             {       
                 
                 boost::shared_ptr<FuelCellShop::Material::CatalystSupportBase > pointer;
                 
                 typename FuelCellShop::Material::CatalystSupportBase::_mapFactory::iterator iterator = FuelCellShop::Material::CatalystSupportBase::get_mapFactory()->find(support_name);
                 
                 if (iterator != FuelCellShop::Material::CatalystSupportBase::get_mapFactory()->end())
                 {
                     if (iterator->second)
                     {
                         pointer =  iterator->second->create_replica();
                     }
                     else 
                     {
                         FcstUtilities::log<<"Pointer not initialized"<<std::endl;
                         abort();
                     }
                 }
                 else
                 {
                     FcstUtilities::log<<"Concrete name in FuelCellShop::Material::CatalystSupportBase::create_CatalystSupport does not exist"<<std::endl;
                     abort();
                 }
                 
                 pointer->initialize(param);
                 
                 return pointer;
             }
            //@}


            ///@name Information and accessors
            //@{             
            
            /** Obtain the electrical conductivity [\p S/cm].
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual double get_electrical_conductivity() const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
                return 0;
            };
            
            /** Obtain the thermal conductivity [\p W/\p (cm-K \p )].
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual double get_thermal_conductivity() const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
                return 0;
            };
            
            /** Obtain the density [\p gm/cm^3].
             * \note This is an abstract virtual function. Please look at the child classes for its reimplementations (if any).
             */
            virtual double get_density() const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
                return 0;
            };   
            //@}
            
        protected:
            ///@name Constructors, destructor, and parameter initalization
            //@{
            /** Constructor 
             */
            CatalystSupportBase()
            :
            BaseMaterial()
            {};
            
            /** Constructor 
             * The constructor initialize parameters using the default values. This is
             * so that if I do not want to call declare_parameters and initialize, I can
             * still use the routine with the hard coded values.
             * 
             * \todo1 What are we going to do with the Carbon, Nafion, Platinum, etc. classes? i.e. puresolids
             */
            CatalystSupportBase(std::string name)
            :
            BaseMaterial(name)
            {};
            
            /**  Destructor  */
            ~CatalystSupportBase()
            {};
            
            /**
             * Declare parameters for a parameter file
             * 
             * \warning This is a PureFunction and it does not declare anything, so please do not call this function in the children.
             * \deprecated
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
             * This object is used to store all objects of type CatalystSupportBase. 
             */
            typedef std::map< std::string, FuelCellShop::Material::CatalystSupportBase* > _mapFactory;      
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
             * This member function is used to create an object of type CatalystSupportBase
             * 
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::Material::CatalystSupportBase > create_replica ()
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            //@}
            ///@name Internal variables
            //@{
            
            /** Electrical conductivity [\p S/cm] of catalyst support extrapolated to 100% solid phase */
            double electrical_conductivity;
            
            /** Thermal conductivity [\p W/\p (cm-K \p )] of catalyst support extrapolated to 100% solid phase */
            double thermal_conductivity;
            
            /** Density of catalyst support [\p gm/cm^3] */
            double density;      
            
            //@}
        };
    }
}

#endif
