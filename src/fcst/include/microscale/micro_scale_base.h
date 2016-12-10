// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: MicroScaleBase
// - Description: Base class for micro scale objects
// - Developers: Philip Wardlaw, University of Alberta
// - $Id: micro_scale_base.h 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__LAYER__MICROSCALE_BASE_H
#define _FUELCELLSHOP__LAYER__MICROSCALE_BASE_H


// The following definitions are necessary for compilation
// MultiScaleCL must be included in this header
// MultiScaleCL makes references to this class
// This class must be declared before we read MultiScaleCL file
namespace FuelCellShop
{
    namespace MicroScale
    {
        class MicroScaleBase;
    }
    namespace Layer
    {
        template <int dim>
        class MultiScaleCL;
    }
}




//FCST classes and utils
#include <utils/fcst_constants.h>
#include <utils/fcst_utilities.h>
#include <layers/multi_scale_CL.h>
#include <utils/logging.h>

//deal.ii classes and utils
#include  <deal.II/base/utilities.h>
#include  <deal.II/base/parameter_handler.h>

//STL containers, functions, streams
#include <vector>
#include <map>
#include <string.h>
#include <cmath>
#include <iostream>

//Boost shared pointer
#include <boost/shared_ptr.hpp>



namespace FuelCellShop
{
    namespace MicroScale
    {
        /**
        *
        * @brief The base class for micro scale objects in OpenFCST.
        *
        *
        * This polymorphic class specifies all the necessary interfaces that must be implemented by a micro
        * scale object in order for it to work with OpenFCST's MultiScaleCl class.
        *
        *
        * <h3> Theory </h3>
        * Micro scale objects primary purpose are to determine reaction rates within the Cl
        * whilst accounting for micro scale mass transport and other losses.
        *
        * Used to calculate reaction rate of CL micro structure [A/cm^3]
        * and micro scale effectiveness (a ratio of CL micro structure
        * reaction rate normalized by  reaction rate calculated using a
        * macro scale homogeneous CL expression).
        *
        *
        * <h3> Input parameters </h3>
        * LIST OF INPUT PARAMETERS FOR THE CLASS.
        * @code
        * ...
        *   subsection MicroScale
        *     set Microscale type = IonomerAgglomerateAnalytical # Chose micro scale type,
        *     #Options as of October 2014 are [ICCP|HybridAgglomerateNumerical|IonomerAgglomerateAnalytical|
        *     #IonomerAgglomerateNumerical|IonomerAgglomerateSun|PolyAgglomerate|WaterAgglomerateNumerical]
        *     #Further subsections for micro scale child classes
        *   end
        * ...
        * @endcode
        *
        * <h3> Usage details</h3>
        * The base class should be used to create an instance of itself, using a deal.ii parameter interface.
        * First <b>declare_MicroScale_parameters</b> must be called. This function will declare all parameters for all children of
        * the microscale base class (Ionomer filled analytical agglomerate, etc.). <b>N.b</b> The parameters will be entered
        * within the current subsection of the parameter object.
        *   Once you have read/initialized your param object use the <b>create_MicroStructure</b> function to create your
        * micro scale object. The object will be automatically initialized.
        *   Use the <b>compute_current</b> function to obtain values for current density and micro scale effectiveness.
        * Function <b>set_solution</b> should be called before calling <b>compute_current</b> or <b>compute_derivative_current</b>, in order
        * for the results to be valid for your present solution.
        *
        *
        * @code
        *
        * //Set up micro scale parameters within catalyst layer context
        * param.enter_subsection(multi_scale_catalyst_layer);
        * FuelCellShop::MicroScale::MicroScaleBase::declare_MicroScale_parameters(param);
        *
        * //Create an object of MicroScaleBase. In this context <b>this</b> is a pointer to type MultiScaleCL
        * boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase> micro = FuelCellShop::MicroScale::MicroScaleBase::create_MicroStructure(param,this);
        * param.leave_subsection();
        *
        * //Set solution to micro scale
        * micro->set_solution(solution_map, reactant_name, solution_index);
        *
        * //Solve
        * double effectiveness;
        * double current_density = micro->compute_current(effectiveness);
        *
        * @endcode
        *
        *
        * @author Philip Wardlaw
        * @date 2014
        */
        class MicroScaleBase
        {

            public:

            /** Destructor */
            virtual ~MicroScaleBase() {};



            ///@name Instance Delivery (Functions)
            //@{

            /**
            * Function used to declare all the data necessary in the parameter files for
            * all MicroScale children.
            *
            * This member function should be used instead of declare_parameters().
            *
            * \param param ParameterHandler object used to store all information about the simulation. Used
            * to read the parameter file. The parameter object should already be set to the desired subsection.
            *
            * The parameter subsection would look as follows:
            *
            * @code
            * subsection MicroScale
            *   set Microscale type = IonomerAgglomerateAnalytical
            *   subsection AgglomerateBase
            *     set Radius of the agglomerate [nm] = 200
            *     set Constant agglomerate parameter [Thickness | Porosity] =  Porosity #Variable used to select if thickness or porosity should be constant
            *     set Thickness of the agglomerate film [nm] = 15
            *     set Agglomerate porosity = 0.25
            *   end
            *
            *   subsection NumericalAgglomerateBase
            *     set Initial condition tolerance factor = 0.05 #Factor determining how often initial solutions are store/updated
            *     set Database name = main_db #Database existing in FCST/Databases. Setting an unused name will create a new database.
            *     set Agglomerate Loading Profile = 1 #A list of doubles of the form "1,2.2,3...", used to weight the Pt. loading profile
            *   end
            *
            *   subsection IonomerAgglomerateNumerical
            *      set Proton Conductivity factor = 1.0
            *      set Use non equilibrium BC = false
            *      set Non Equilibrium BC Rate constant = 1.0
            *   end
            *
            *   #No parameters so far for other microscale types
            *
            * end
            *
            * @endcode
            */
            static void declare_MicroScale_parameters (ParameterHandler &param)
            {
                param.enter_subsection("MicroScale");
                std::string micro_type_list;

                for (typename FuelCellShop::MicroScale::MicroScaleBase::_mapFactory::iterator iterator = FuelCellShop::MicroScale::MicroScaleBase::get_mapFactory()->begin();
                        iterator != FuelCellShop::MicroScale::MicroScaleBase::get_mapFactory()->end(); iterator++)
                {
                    iterator->second->declare_parameters(param);
                    micro_type_list += iterator->first +"|";
                }

                micro_type_list =  micro_type_list.substr(0,  micro_type_list.size()-1);//Trim string 1 '|' char
                param.declare_entry("Microscale type", "IonomerAgglomerateAnalytical", Patterns::Selection(micro_type_list));
                param.leave_subsection();

            }

            /**
            * Function used to select the appropriate MicroScale type as specified in the ParameterHandler under
            * line
            * @code
            * set Microscale type = IonomerAgglomerateAnalytical
            * @endcode
            * current options are [ IonomerAgglomerateAnalytical | IonomerAgglomerateNumerical | IonomerAgglomerateSun
            * | WaterAgglomerateNumerical | WaterConicalPoreAgglomerate]
            *
            * The parameter object must already be set to the appropriate subsection.
            *
            * The MultiScaleCL that calls this function must pass a pointer of itself so that the micro structure
            * class can call MultiScaleCl->get_resources() and ->get_properties().
            *
            */
            static boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase > create_MicroStructure (ParameterHandler &param,
                    FuelCellShop::Layer::MultiScaleCL<deal_II_dimension>* layer)
            {
                boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase > pointer;

                param.enter_subsection("MicroScale");
                std::string concrete_name = param.get("Microscale type");


                typename FuelCellShop::MicroScale::MicroScaleBase::_mapFactory::iterator iterator = FuelCellShop::MicroScale::MicroScaleBase::get_mapFactory()->find(concrete_name);

                if (iterator != FuelCellShop::MicroScale::MicroScaleBase::get_mapFactory()->end())
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
                    FcstUtilities::log<<"Concrete name in FuelCellShop::MicroScale::MicroScaleBase::create_MicroStructure does not exist"<<std::endl;
                    abort();
                }

                pointer->set_layer(layer);
                pointer->initialize(param);
                pointer->set_structure();
                param.leave_subsection();
                return pointer;
            }
            //@}

            ///@name Initalization
            //@{
            /**
            * Function for setting the solution map(reactant concentration, phi_s, phi_m, etc.).
            * First argument provide the map of SolutionVariable. The second argument provide
            * the name of the primary reactant. The final argument is the index of the solution
            * map that the micro scale object should solve for.
            *
            * This function should be called immediatly before <b>compute_current</b> or
            * <b>compute_derivative_current</b>.
            *
            */
            virtual void set_solution(const std::map<VariableNames,SolutionVariable>&,const VariableNames&, const int&) = 0;
            //@}

            /**
            * Function used to compute the current density produced by the micro structure.
            * Several other functionals may be returned (effectiveness, coverages)
            * Returns current density in [A/cm^3 of sold volume]
            *
            * Solves for solution variables set by the last call to <b>set_solution</b>.
            */
            virtual SolutionMap compute_current () = 0;

            /**
            * Function to compute the derivative of the current density at the local operating conditions.
            * Returns current density derivatives in [A/cm^3 of solid/liquid phase micro structure]
            *

            * Solves for solution variables set by the last call to <b>set_solution</b>.
            *
            * <h3> Usage details</h3>
            * Call <b>has_derivatives</b> to check if it is OK to call this function.
            *
            */
            virtual std::vector<double> compute_derivative_current ()
            {
                Assert(false, ExcPureFunctionCalled());
                return std::vector<double>(3, 0.0);
            }

            /**
            * Returns true if the class instance can calculate
            * current density derivatives.
            *
            * <h3> Usage details</h3>
            * Call this function if unsure whether or not to call <b>compute_derivative_current</b>.
            *
            */
            virtual bool has_derivatives() = 0;

            /**
            * Return name of class instance, i.e. concrete name.
            */
            virtual std::string get_name() = 0;

            /**
            * Print out key micro-structural dimensions, defined by child.
            */
            virtual void print_properties() = 0;


            /**
            * MicroScale object may have extra contribution to volume of layer, e.g. water.
            *
            * Returns a fraction.
            */
            virtual double aux_volume_fraction() = 0;


            /*
            * MicroScale will deep create new instances of pointers (kinetics, materials), to prevent race condition during parallel execution.
            */
            virtual void make_thread_safe(ParameterHandler &param, unsigned int thread_index) = 0;


        protected:

            /*
            * Protected member function, for storing address to catalyst layer which
            * constructed this MicroScale object.
            */
            void set_layer(FuelCellShop::Layer::MultiScaleCL<deal_II_dimension>* _layer){

                if(_layer!= NULL){
                    layer = _layer;
                }
                else
                {
                    std::string msg = "MultiScaleCl pointer passed to MicroScaleBase is NULL!";
                    throw std::runtime_error(msg);
                }
            }

            /*
            * Protected member pointer, for storing address to catalyst layer which
            * constructed this MicroScale object.
            */
            FuelCellShop::Layer::MultiScaleCL<deal_II_dimension>* layer;

            /*
            * Protected pure virtual member function for declaring parameters, must be implemented
            * by all children. Parameter structure is hierarchical, therefore children
            * should call their parent's declare parameter function from their own declare
            * parameter function.
            */
            virtual void declare_parameters (ParameterHandler &param) const =0;

            /*
            * Protected pure virtual member function for initializing parameters, must be implemented
            * by all children. Parameter structure is hierarchical, therefore children
            * should call their parent's initialize function from their own initialize function.
            */
            virtual void initialize (ParameterHandler &param) =0;


            /*
            *  Constructor
            */
            MicroScaleBase(){}


            /*
            * Protected member function for setting up MicroScale structure after object
            * initialization.
            *
            * Responsible for pulling structural properties from catalyst layer.
            *
            * Responsible for setting up pointer to materials and kinetics.
            *
            */
            virtual void set_structure () = 0;

            ///@name Instance Delivery (Types)
            //@{
            /**
            * This object is used to store all objects of type MicroScaleBase.
            */
            typedef std::map< std::string, MicroScaleBase* > _mapFactory;
            //@}


            ///@name Instance Delivery (Protected functions)
            //@{
            /**
            * This member function is used to create an object of type MicroScaleBase
            *
            * \warning This class MUST be redeclared in every child.
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
            * This member function is used to create an object of MicroScaleBase
            *
            * \warning This class MUST be redeclared in every child.
            */
            virtual boost::shared_ptr<MicroScaleBase> create_replica () = 0;
            //@}

        };
    }
}

#endif
