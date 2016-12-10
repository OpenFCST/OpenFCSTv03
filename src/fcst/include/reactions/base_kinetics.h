//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: base_kinetics.h
//    - Description: Base Kinetics. It implements the interface for other kinetics classes
//        and some common methods.
//    - Developers: M. Moore, M. Bhaiya, M. Secanell and V. Zingan
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__BASE__KINETICS_H
#define _FUELCELLSHOP__BASE__KINETICS_H

// Include deal.II classes
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/fe/fe_values.h>

// Include STL
#include <cmath>
#include <iostream>
#include <map>
#include <algorithm>

// Include OpenFCST routines:
#include <layers/base_layer.h>
#include <materials/catalyst_base.h>
#include <materials/polymer_electrolyte_material_base.h>
#include <utils/fcst_constants.h>
#include <application_core/system_management.h>

using namespace dealii;

namespace FuelCellShop
{
    namespace Kinetics
    {
        /**
         * Virtual class used to provide the interface for all kinetic/reaction children.
         *
         * The main purpose of kinetic objects is to calculate the rate of the reaction. In the case
         * of fuel cells, instead of the rate we return the current that is produced per unit area of
         * catalyst in the electrode.
         *
         * To use this class, first the object needs to be initialized by doing two
         * very important things:
         * - Assign the materials to be used by the class using the member functions: #set_catalyst and #set_electrolyte
         * - Initialize the solution variables using the member functions: #set_electrolyte_potential,
         * #set_solid_potential and so on.
         *
         * No object of type BaseKinetics should ever be created, instead this class
         * is used to initialize pointers of type BaseKinetics and to provide a common
         * interface for all Kinetics classes.
         *
         * The class has a database of children such that it will declare all
         * necessary parameters for all children in the
         * input file, read the input file, create the appropriate children and return a pointer
         * to BaseKinetics with the children selected.
         *
         *
         * <h3>Usage Details:</h3>
         *
         * In order to create a kinetics object within an reactive layer, the following steps need to be taken.
         *
         * First, in the CatalystLayer object .h file that needs a kinetics object, create a pointer to a BaseKinetics object, i.e.
         * @code
         * boost::shared_ptr<FuelCellShop::Kinetics::BaseKinetics<dim> >ORR;
         * @endcode
         *
         * This pointer object will be available anywhere inside the reactive layer. Because we do not want to
         * worry about deleting the pointer afterwards, we use a Boost pointer which has its own memory management
         * algorithms. See the <a href="http://www.boost.org/doc/libs/1_54_0/libs/smart_ptr/smart_ptr.htm">Boost website</a> for more information
         *
         * Once the pointer is available, we need to do three things in the reactive layer
         *
         * - Call declare_Kinetics_parameters in the reactive layer declare_parameters() member function. This
         * member function is used to define all parameters that can be read from the input file
         *
         * - Call create_Kinetics and initialize. The former member function will fill the pointer created above with the appropriate
         * gas diffusion layer you have selected in the input file. Then, initialize reads from the input file all the data from the file in order
         * to setup the object.
         *
         * The object is ready for use now.
         *
         * For a code example see unit_tests/source/DT_test.cc.
         *
         * @endcode
         *
         * \author M. Moore, 2011-12
         * \author M. Bhaiya, 2013
         * \author M. Secanell, 2013
         *
         */
        class BaseKinetics
        {
        public:

            /**
             * Destructor
             */
            virtual ~BaseKinetics()
            {};
            
            ///@name Instance Delivery (Public functions)
            /**
             * Function used to declare all the data necessary in the parameter files for
             * all BaseKinetics children.
             *
             */
            static void declare_Kinetics_parameters (ParameterHandler &param)
            {

                for (typename FuelCellShop::Kinetics::BaseKinetics::_mapFactory::iterator iterator = FuelCellShop::Kinetics::BaseKinetics::get_mapFactory()->begin();
                     iterator != FuelCellShop::Kinetics::BaseKinetics::get_mapFactory()->end();
                iterator++)
                     {
                         iterator->second->declare_parameters(param);
                     }
            }

             /**
              *
              * Function called in create_CatalystLayer and used to select the appropriate BaseKinetics type that will be used
              * in the layer. The name of the kinetics object to be used is provided in kinetics_name.
              *
              * The name of the kinetics object is provided in the ParameterHandler in the CatalystLayer subsection as follows:
              *
              * @code
              * subsection Catalyst Layer Properties <- This name is the name of the catalyst layer subsection where the kinetics are taking place.
              * (...)
              *   set Kinetics type = TafelKinetics
              * (...)
              * end
              * @endcode
              * current options are [ TafelKinetics | DoubleTrapKinetics | DualPathKinetics | BVKinetics ]
              *
              *
              */
             static boost::shared_ptr<FuelCellShop::Kinetics::BaseKinetics > create_Kinetics (ParameterHandler &param,
                                                                                              std::string kinetics_name)
             {

                 boost::shared_ptr<FuelCellShop::Kinetics::BaseKinetics > pointer;

                 typename FuelCellShop::Kinetics::BaseKinetics::_mapFactory::iterator iterator = FuelCellShop::Kinetics::BaseKinetics::get_mapFactory()->find(kinetics_name);

                 if (iterator != FuelCellShop::Kinetics::BaseKinetics::get_mapFactory()->end())
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
                     FcstUtilities::log<<"Concrete name does not exist"<<std::endl;
                     abort();
                 }

                 pointer->initialize(param);

                 return pointer;
             }
            //@}

            ///@name Initalization
            //@{

            /**
             * Set the electrolyte phase potential. The potential should be in Volts.
             *
             */
            void set_electrolyte_potential(const SolutionVariable& phi)
            {
                Assert( phi.get_variablename() == protonic_electrical_potential, ExcMessage("Wrong solution variable passed in BaseKinetics::set_electrolyte_potential.") );
                phi_m = phi;
            }
            /**
             * Set the solid phase potential. The potential should be in Volts.
             */
            void set_solid_potential(const SolutionVariable& phi)
            {
                Assert( phi.get_variablename() == electronic_electrical_potential, ExcMessage("Wrong solution variable passed in BaseKinetics::set_solid_potential.") );
                phi_s = phi;
            }
            /**
             * Set temperature. The temperature should be in Kelvin.
             */
            void set_temperature(const SolutionVariable& temperature)
            {
                Assert( temperature.get_variablename() == temperature_of_REV, ExcMessage("Wrong solution variable passed in BaseKinetics::set_temperature") );
                T = temperature;
            }
            /**
             * Set reactant concentrations at the catalyst/electrolyte interface. This should be the concentration in Nafion.
             * 
             * The concentrations should be in \f$ mol/cm^3 \f$. 
             * 
             * \warning It should be ensured that no molar fractions are set using this method.
             */
            void set_reactant_concentrations(const std::vector<SolutionVariable>& conc_vec)
            {
                Assert( conc_vec.size() > 0, ExcMessage("Atleast one reactant concentration should be passed inside BaseKinetics::set_reactant_concentrations. !!") );

                for ( unsigned int i=0; i<conc_vec.size(); ++i )
                {
                    Assert( (conc_vec[i].get_variablename() != hydrogen_molar_fraction &&
                            conc_vec[i].get_variablename() != oxygen_molar_fraction   &&
                            conc_vec[i].get_variablename() != water_molar_fraction), ExcMessage("Molar fractions can't be setup using BaseKinetics::set_reactant_concentrations.") );

                    reactants_map[ conc_vec[i].get_variablename() ] = conc_vec[i];
                }
            }

            /**
             * Set the variables for which you would like to compute the derivaitives
             */
            void set_derivative_flags(const std::vector<VariableNames>& flags)
            {
                derivative_flags = flags;
            };

            /**
             * Set a pointer to the catalyst that will be used. The catalyst is specified
             * at the application level and passed to the kinetics class using this function.
             */
            void set_catalyst(FuelCellShop::Material::CatalystBase* cat_in)
            {
                catalyst = cat_in;
            }

            /**
             * Member function used to set the electrolyte pointer to that used by the application.
             * Needed so that conversion from molar fraction to concentration can be done. (Henry's
             * constant in particular is needed).
             */
            void set_electrolyte(FuelCellShop::Material::PolymerElectrolyteBase* electrolyte_in)
            {
                electrolyte = electrolyte_in;
            }

            /**
             * Member function used to specify the reaction for which the kinetic parameters are needed, for example
             * for a Platinum catalyst, we can specify that we need the kinetic parameters for either the oxygen reduction reaction (ORR)
             * or the hydrogen oxidation reaction (HOR)
             */
            virtual void set_reaction_kinetics(const ReactionNames  name)
            {
                name_reaction_kinetics = name;
            }

            /** Set the total gas pressure [\p Pascals] in the cell*/
            void set_p_t(const double& P_Tot)
            {
                p_total = P_Tot;
            }
            //@}
            ///@name Information and accessors
            //@{
            /**
             * Returns the name of the reaction that this kinetics class is using.
             */
            inline ReactionNames get_reaction_name() const
            {
                return name_reaction_kinetics;
            }

            /** Function to get pointer to catalyst class. This function is basically needed by thermal source class. */
            FuelCellShop::Material::CatalystBase* get_cat() const
            {
                return catalyst;

            }

            /**
             * Function to return the value of the current [\p A/cm^2].
             */
            virtual void current_density (std::vector<double>&)
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;

            };

            /**
             * Function to return the derivative of the current density w.r.t solution variables. It returns a map of vectors containing
             * derivative w.r.t solution variables / design variables set using #set_derivative_flags method. Each vector can be accessed by using \p Key
             * of the map, which correpsonds to the #VariableNames (solution/design variable). This method takes input map by reference, hence the map is needed
             * to be created at application/equation level with default arguments and passed inside this method.
             */
            virtual void derivative_current (std::map< VariableNames, std::vector<double> >&)
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            }

            /**
             * Used to return the coverage of the intermediate species if they are computed.
             */
            virtual void compute_coverages(const std::string& name_species,
                                           std::vector<double>& coverage) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            };
            /**
             * Used to return the coverage of the intermediate species if they are computed.
             *
             * \deprecated
             */
            virtual void OH_coverage(std::vector<double>& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            };

            /**
             * Used to return the coverage of the intermediate species if they are computed.
             *
             * \deprecated
             */
            virtual void O_coverage(std::vector<double>& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            }

            //@}

            virtual bool has_coverage(const VariableNames& type){
                return false;
            }


        protected:

            ///@name Instance Delivery (Types)
            //@{
            /**
             * This object is used to store all objects of type BaseKinetics.
             */
            typedef std::map< std::string, BaseKinetics* > _mapFactory;
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
             * This member function is used to create an object of type BaseKinetics
             *
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::Kinetics::BaseKinetics > create_replica ()
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            //@}

            ///@name Constructors, destructor, and initalization
            //@{
            /**
             * Constructor
             */
            BaseKinetics()
            {
                //FcstUtilities::log << "Creating kinetics object of type ";

                catalyst = NULL;
                electrolyte = NULL;

                // Universal constants
                R = Constants::R();
                F = Constants::F();
                K = Constants::K();

                kin_param_initialized = false;
            };

            /**
             * Declare parameters for a parameter file
             */
            virtual void declare_parameters (ParameterHandler&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;

            };

            /**
             * Member function used to read in data and initialize the necessary data
             * to compute the coefficients.
             */
            virtual void initialize (ParameterHandler&)
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            };
            //@}

            ///@name Material objects and other internal data/methods
            //@{

            /**
             * Pure abstract method used to initialize certain kinetics parameters which are normally constant, \em e.g. reference concentrations, reaction order
             * etc. This function is necessarily required to be reimplemented in all of the child kinetics classes.
             */
            virtual void init_kin_param() = 0;

            /**
             * Boolean variable to determine whether #init_kin_param has been already called or not.
             */
            bool kin_param_initialized;

            /**
             * Enumeration with the reaction name for which the class returns kinetic parameters.
             */
            ReactionNames name_reaction_kinetics;

            /**
             * Flags for derivatives: These flags are used to request derivatives which are
             * computed using the derivative_current function.
             */
            std::vector<VariableNames> derivative_flags;

            /**
             * Number of quadrature points in the cell.
             */
            unsigned int n_quad;

            /**
             * Pointer to the catalyst object that is created at the application level and
             * passed to the kinetics class using the set_catalyst function.
             */
            FuelCellShop::Material::CatalystBase* catalyst;
            /**
             * Pointer to the electrolyte object created in the application that is used
             * to calculate the properties of the electrolyte in the catalyst layer.
             */
            FuelCellShop::Material::PolymerElectrolyteBase* electrolyte;
            //@}

            ///@name Universal constants
            //@{
            /** Universal gas constant */
            double R;
            /** Universal Farday's constant */
            double F;
            //Boltzmann constant
            /** Boltzmann constant */
            double K;
            //@}

            ///@name Solution variables
            //@{
            /** Total gas pressure [\p Pascals] in the cell for isobaric case. */
            double p_total;

            /**
             * Map of SolutionVariables storing a pointer to the solution vector
             * storing the concentration of each one of the reactants implemented
             */
            std::map< VariableNames, SolutionVariable > reactants_map;

            /** Struct storing a pointer to the solution vector for the electrolyte potential */
            SolutionVariable phi_m;

            /** Struct stroing a pointer to the solution vector for the electronic/solid potential */
            SolutionVariable phi_s
            ;
            /** Struct stroing a pointer to the solution vector for the temperature */
            SolutionVariable T;
            //@}

        }; // BaseKinetics
    } // Kinetics
} // FuelCellShop

#endif
