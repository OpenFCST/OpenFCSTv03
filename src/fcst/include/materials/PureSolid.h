//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: PureSolid.h
//    - Description: Material class for pure liquids.
//    - Developers: Marc Secanell, Jie Zhou
//    - $Id: PureSolid.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__PURESOLID__H
#define _FUELCELLSHOP__PURESOLID__H

//Include STL
#include <cstdlib>
#include <vector>
#include <fstream>
#include <string>
#include <map>

#include <deal.II/base/parameter_handler.h>

#include <materials/base_material.h>
#include <application_core/system_management.h>

enum SolidMaterialTypes
{
    nosolidmaterial = 0,
    Graphite,
    DummySolid
};

namespace FuelCellShop

{
    namespace Material
    {
        
        /**
         * This class is a base class for all pure solid materials
         * used in FCST.
         *
         * This class contains data provided by its children
         * and implements the methods by means of which
         * different properties of pure solid materials
         * can be computed.
         * 
         * \warning Do not create an object of this class. Use its children instead.
         *  *
         * \author Marc Secanell, Jie Zhou, 2014
         */
        
        class PureSolid : public BaseMaterial
        {
        public:
            
            ///@name Constructors, destructor and initialization
            //@{
            
            /** 
             *Consturctor
             */
            PureSolid()
            :
            BaseMaterial()
            {};
            
            /** 
             *Destructor 
             */
            virtual ~PureSolid();
            
            /**
             * This routine is used to create a PureSolid with the desired properties. \param concrete_name should
             * contain a string with the concrete_name of the PureSolid child that you would like to create. \param param
             * should contain the ParameterHandler handler from which you would harvest any necessary data under subsection             * 
             * @code
             * subsection Material Database
             *   subsection Pure Solid
             *     subsection concrete_name
             * 
             *     end 
             *   end
             * end
             * @endcode
             * 
             * \note Concrete_name should be obtained form the layer class.             */
            static void declare_PureSolid_parameters (ParameterHandler &param)
            {
                
                for (typename FuelCellShop::Material::PureSolid::_mapFactory::iterator iterator = FuelCellShop::Material::PureSolid::get_mapFactory()->begin(); 
                     iterator != FuelCellShop::Material::PureSolid::get_mapFactory()->end(); 
                     iterator++)
                     {
                         param.enter_subsection("Material Database");
                         {
                             param.enter_subsection("PureSolid");
                             {
                                 iterator->second->declare_parameters(param);
                             }
                             param.leave_subsection();                             
                         }
                         param.leave_subsection();
                     }    
                  
            }
            
            /**
             * This function returns a boost shared ptr of a certain material with the name of concrete_name
             * and it initializes the PureSolid material with the desired properties
             * 
             */
            static boost::shared_ptr<FuelCellShop::Material::PureSolid > create_PureSolid (std::string concrete_name,
                                                                                           ParameterHandler &param)
             {       
                 
                 boost::shared_ptr<FuelCellShop::Material::PureSolid> pointer;
                 
                 typename FuelCellShop::Material::PureSolid::_mapFactory::iterator iterator = FuelCellShop::Material::PureSolid::get_mapFactory()->find(concrete_name);
                 
                 if (iterator != FuelCellShop::Material::PureSolid::get_mapFactory()->end())
                 {
                     if (iterator->second)
                     {
                         pointer = iterator->second->create_replica(concrete_name);
                     }
                     else 
                     {
                         FcstUtilities::log<<"Pointer not initialized"<<std::endl;
                         abort();
                     }
                 }
                 else
                 {
                     AssertThrow(false, ExcMessage("Concrete name in FuelCellShop::Material::PureSolid::create_PureSolid does not exist"));
                 }
                 
                 param.enter_subsection("Material Database");
                 {
                     param.enter_subsection("PureSolid");
                     {
                         pointer->initialize(param);
                     }
                     param.leave_subsection();                             
                 }
                 param.leave_subsection();                
                 
                 return pointer;
                 
             }
            //@} 
            
            ///@name Service functions. 
            //@{
            
            /** Obtain the density */
            virtual  double get_density() const = 0;
            
            /** Obtain the electrical conductivity units (S/M)*/
            virtual  double get_electrical_conductivity(double temperature) const = 0 ;
            
            /** 
             * Obtain the electrical conductivity (S/M). 
             */
            virtual  void get_electrical_conductivity(std::vector<double> temperature,std::vector<double>& dst) const = 0;
            
            /** Obtain the derivative of the electrical conductivity */
            virtual  double get_Delectrical_conductivity_Dtemperature(double temperature) const = 0;
            
            /** Obtain the derivative of the electrical conductivity */
            virtual  void get_Delectrical_conductivity_Dtemperature(std::vector<double> temperature,std::vector<double> &dst) const = 0;
            
            /** Obtain the thermal conductivity (watts/m K)*/
            virtual   double get_thermal_conductivity(double temperature) const = 0;
            
            /** Obtain the thermal conductivity (watts/m K)*/
            virtual  void get_thermal_conductivity(std::vector<double> temperature,std::vector<double>&dst) const = 0;
            
            /** Obtain the derivative of the thermal conductivity */
            virtual  double get_Dthermal_conductivity_Dtemperature(double temperature) const = 0;
            
            /** Obtain the derivative of the thermal conductivity */
            virtual  void get_Dthermal_conductivity_Dtemperature(std::vector<double> temperature,std::vector<double> &dst) const = 0;
            
            /** Obtain the coefficient_thermal_expansion (microns/m °C)*/
            virtual   double get_coefficient_thermal_expansion(double temperature) const = 0;
            
            /** Obtain the compressive_strength (unit N/mm2)*/
            virtual   double get_compressive_strength(double temperature) const = 0;
            
            /** Obtain the H2_permeability (unit cm3*cm-2*s-1)*/
            virtual   double get_H2_permeability(double temperature) const = 0;
            
            /** Obtain the Poissons ratio (unitless)*/
            virtual   double get_Poissons_ratio(double temperature) const = 0;

            //@}
            
            ///@name Instance Delivery (Types)
            //@{         
            /** 
             * This object is used to store all objects of type PureSolid. 
             * The mapping is from string which is basic called concrete_name to a pointer pointing to the children of PureSolid
             */
            typedef std::map< std::string, PureSolid* > _mapFactory;      
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
            


        protected:

            ///@name Instance Delivery (Public variables)
            //@{         
            /**
             * Concrete name used for objects of this class. This name is used when
             * setting up the subsection where the data is stored in the input file.
             * 
             * The data will be store under
             * @code
             * subsection Material Database
             *   subsection Pure Solid
             *     subsection concrete_name
             * 
             *     end 
             *   end
             * end
             * @endcode
             */
            static const std::string concrete_name;
            //@} 
            
            ///@name Instance Delivery (Public variables)
            //@{         
            /**
             * PROTOTYPE used for pointing to this class. This pointer is used when
             * setting up the subsection where the data is stored in the input file.
             * 
             * The data will be store under
             * @code
             * subsection Material Database
             *   subsection Pure Solid
             *     subsection concrete_name
             * 
             *     end 
             *   end
             * end
             * @endcode
             */
            static PureSolid const* PROTOTYPE;
            //@} 
            
            ///@name Instance Delivery (Private functions)
            //@{         
            /**
             * This member function is used to create an object of type PureSolid Material
             * 
             * \warning This class MUST be redeclared in every child.
             */                      
            virtual boost::shared_ptr<FuelCellShop::Material::PureSolid > create_replica (const std::string &name)
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            
            //@}
            
            ///@name Constructors, destructor and initialization
            //@{
            
            /** 
             *Constructor 
             */
            PureSolid(const std::string& name);
           
            /**
             * Declare parameters for a parameter file.
             * 
             */
            virtual void declare_parameters (ParameterHandler &param) const {};
            
            /**
             * Member function used to read in data and initialize the necessary data
             * to compute the coefficients.
             */
            virtual void initialize (ParameterHandler &param) {};
            
            //@}
            
            //////////
            // DATA //
            //////////
            
            ///@name Accessors and info
            //@{
            /** Density */
            double density;
            
            /** electrical_conductivity */
            double electrical_conductivity;
            
            /** thermal_conductivity */
            double thermal_conductivity;
            
            /** Coefficient of Thermal Expansion */
            double coefficient_thermal_expansion;
            
            /** Compressive strength */
            double compressive_strength;        
            
            /** H2 permeability */
            double H2_permeability; 
            
            /** Poissons_ratio */
            double Poissons_ratio; 
            
//             /** Solutions at every quadrature point inside the cell. 
//             * for this mapping key is the Enumerator of VariableName_of_REV
//             * the mapped value is of type SolutionVariable
//             */
//             std::map <VariableNames, SolutionVariable> Solutions;
//             
//             SolutionVariable T_vector;
            //@}
            
        };

        /**
         * This class describes properties of pure Poco \p Graphite.
         * Use this name in the parameters file if needed.
         * 
         * All default material properties are from: 
         * http://www.poco.com/MaterialsandServices/Graphite/SemiconductorGrades/HPD.aspx   
         *
         * \author Marc Secanell, Jie Zhou
         */
        
        class Graphite : public PureSolid
        {
        public:
            
            ///@name Constructors, destructor and initialization
            //@{
            /**
             * Constructor.
             */
            Graphite();
            /**
             * Destructor.
             */
            virtual ~Graphite();
            //@}
            
            ///@name Service functions.
            //@{
                        
            /** Obtain the density units (g/cm3)**/
            virtual inline double get_density() const;
            
            /** Obtain the electrical conductivity units (S/M)*/
            virtual inline double get_electrical_conductivity(double temperature) const;
            
            /** Obtain the electrical conductivity */
            virtual inline void get_electrical_conductivity(std::vector<double> temperature,std::vector<double>& dst) const;
            
            /** Obtain the derivative of the electrical conductivity */
            virtual inline double get_Delectrical_conductivity_Dtemperature(double temperature) const 
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
                return 0;
            }
            
            /** Obtain the derivative of the electrical conductivity */
            virtual inline void get_Delectrical_conductivity_Dtemperature(std::vector<double>,std::vector<double>&) const {}
            
            /** Obtain the thermal conductivity units (watts/m K)*/
            virtual  inline double get_thermal_conductivity(double temperature) const;
            
            /** Obtain the thermal conductivity */
            virtual inline void get_thermal_conductivity(std::vector<double>,std::vector<double>&) const;
            
            /** Obtain the derivative of the thermal conductivity */
            virtual inline double get_Dthermal_conductivity_Dtemperature(double temperature) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
                return 0;
            }
            
            /** Obtain the derivative of the thermal conductivity */
            virtual inline void get_Dthermal_conductivity_Dtemperature(std::vector<double>,std::vector<double>&) const {}
            
            /** Obtain the coefficient_thermal_expansion units (microns/m °C)*/
            virtual  inline double get_coefficient_thermal_expansion(double temperature) const;
            
            /** Obtain the coefficient_thermal_expansion units (N/mm2)*/
            virtual  inline double get_compressive_strength(double temperature) const;
            
            /** Obtain the H2_permeability units (cm3*cm-2*s-1)*/
            virtual  inline double get_H2_permeability(double temperature) const;
            
            /** Obtain the Poissons_ratio units ()*/
            virtual  inline double get_Poissons_ratio(double temperature) const;
            //@}
            
        protected:
            
            ///@name Instance Delivery (Private functions)
            //@{         
            /**
             * This member function is used to create an object of type PureSolid Material
             * 
             * \warning This class MUST be redeclared in every child.
             */    
            virtual boost::shared_ptr<FuelCellShop::Material::PureSolid > create_replica (const std::string &name)
            {
                return boost::shared_ptr<FuelCellShop::Material::PureSolid > (new FuelCellShop::Material::Graphite ());
            }
            //@}
            
            ///@name Instance Delivery (Public variables)
            //@{         
            /**
             * Concrete name used for objects of this class. This name is used when
             * setting up the subsection where the data is stored in the input file.
             */
            static const std::string concrete_name;
            /** 
             * PROTOTYPE used for pointing to this class. This pointer is used when
             * setting up the subsection where the data is stored in the input file.
             */
            static Graphite const* PROTOTYPE;
            //@}
            
            ///@name Constructors, destructor, and initialization
            //@{            
            /**
             * Constructor.
             */
            Graphite(const std::string& name);
            
            /**
             * Declare parameters for a parameter file.
             * 
             */
            virtual void declare_parameters (ParameterHandler &param) const;
            
            /**
             * Member function used to read in data and initialize the necessary data
             * to compute the coefficients.
             */
            virtual void initialize (ParameterHandler &param);
            //@}
            
        };
        
        /**
         * This class describes properties of pure \p Dummy.
         * Use this name in the parameters file if needed.
         *
         * \author Marc Secanell, Jie Zhou
         */
        
        class DummySolid : public PureSolid
        {
        public: 
            ///@name Constructors, destructor and initialization
            //@{
            /**
             * Constructor.
             */
            DummySolid();
            /**
             * Destructor.
             */
            virtual ~DummySolid();
            //@}
            
            ///@name Service functions.
            //@{
            /** Obtain the density units (g/cm3)**/
            virtual inline double get_density() const;
            
            /** Obtain the electrical conductivity units (S/M)*/
            virtual inline double get_electrical_conductivity(double temperature) const ;
            
            /** Obtain the electrical conductivity */
            virtual inline void get_electrical_conductivity(std::vector<double> temperature,std::vector<double>& dst) const;
            
            /** Obtain the derivative of the electrical conductivity */
            virtual inline double get_Delectrical_conductivity_Dtemperature(double temperature) const;
            
            /** Obtain the derivative of the electrical conductivity */
            virtual inline void get_Delectrical_conductivity_Dtemperature(std::vector<double> temperature ,std::vector<double>& dst) const;
            
            /** Obtain the thermal conductivity units (watts/m K)*/
            virtual  inline double get_thermal_conductivity(double temperature) const;
            
            /** Obtain the thermal conductivity */
            virtual inline void get_thermal_conductivity(std::vector<double> temperature,std::vector<double>& dst) const;
            
            /** Obtain the derivative of the thermal conductivity */
            virtual inline double get_Dthermal_conductivity_Dtemperature(double temperature) const;
            
            /** Obtain the derivative of the thermal conductivity */
            virtual inline void get_Dthermal_conductivity_Dtemperature(std::vector<double> temperature,std::vector<double>& dst) const;
            
            /** Obtain the coefficient_thermal_expansion units (microns/m °C)*/
            virtual  inline double get_coefficient_thermal_expansion(double temperature) const;
            
            /** Obtain the coefficient_thermal_expansion units (N/mm2)*/
            virtual  inline double get_compressive_strength(double temperature) const;
            
            /** Obtain the H2_permeability units (cm3*cm-2*s-1)*/
            virtual  inline double get_H2_permeability(double temperature) const;
            
            /** Obtain the Poissons_ratio units ()*/
            virtual  inline double get_Poissons_ratio(double temperature) const;
            //@}
            
        protected:
            
            ///@name Instance Delivery (Private functions)
            //@{         
            /**
             * This member function is used to create an object of type PureSolid Material
             * 
             * \warning This class MUST be redeclared in every child.
             */  
            virtual boost::shared_ptr<FuelCellShop::Material::PureSolid > create_replica (const std::string &name)
            {
                return boost::shared_ptr<FuelCellShop::Material::PureSolid > (new FuelCellShop::Material::DummySolid ());
            } 
            //@}
            
            /**
             * Concrete name used for objects of this class. This name is used when
             * setting up the subsection where the data is stored in the input file.
             */
            static const std::string concrete_name;
            /** 
             * PROTOTYPE used for pointing to this class. This pointer is used when
             * setting up the subsection where the data is stored in the input file.
             */
            static DummySolid const* PROTOTYPE;
            //@}
            
            
            ///@name Constructors, destructor, and initialization
            //@{
            
            /**
             * Constructor.
             */
            DummySolid(const std::string& name);
            /**
             * Declare parameters for a parameter file.
             * 
             */
            virtual void declare_parameters (ParameterHandler &param) const;
            
            /**
             * Member function used to read in data and initialize the necessary data
             * to compute the coefficients.
             */
            virtual void initialize (ParameterHandler &param);
            //@}
        };
        
    }
}

#endif