// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: initial_and_boundary_data.h
// - Description: This namespace contains data and methods
//                that handle initial and boundary data of
//                a problem at hand
// - Developers: Valentin N. Zingan, University of Alberta
//               M. Secanell (documentation), UofA
// - $Id: initial_and_boundary_data.h 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#ifndef _FCST_FUELCELL_INITIAL_AND_BOUNDARY_DATA_H_
#define _FCST_FUELCELL_INITIAL_AND_BOUNDARY_DATA_H_

#define _IS_NOT_CONSTANT_ 1.e300
//-- dealII
#include <deal.II/base/point.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/dofs/function_map.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

//-- OpenFCST
#include <application_core/application_data.h>
#include <application_core/system_management.h>
#include <solvers/solver_utils.h>

using namespace dealii;
using namespace FuelCell::ApplicationCore;

/**
 * The typedef for the std::map that reflects the following
 * structure:
 *
 * - \p first  \p argument : name of the solution component,
 * - \p second \p argument : material id,
 * - \p third  \p argument : value of the solution component.
 * 
 * This map is usually used to initialize the solution of the problem. If you would like an initial
 * solution that is different for each matrial id, you can specify this here. For example, imagine
 * I want to specify for material_id 0, an oxygen molar fraction of 0.1 and for material_id 2,
 * an oxygen molar fraction of 0.2. Then, I would write
 * @code
 * std::map<types::material_id, double> values;
 * values[0] = 0.1;
 * values[1] = 0.2;
 * component_materialID_value_map['oxygen_molar_fraction'] = values;
 * @endcode
 * 
 * In general, this object is initialized directly from the input file in the Equations section of the data file.
 * For example, in the input file if I write:
 * @code
 * subsection Equations
 *   subsection Ficks Transport Equation - oxygen
 *     subsection Initial data
 *       set oxygen_molar_fraction = 0: 0.1, 1: 0.2
 *     end
 *   end
 * end
 * @endcode
 * This would be equivalent to the code above.
 */
typedef std::map<   std::string , std::map<types::material_id, double>   > component_materialID_value_map;

/**
 * The typedef for the std::map
 * that reflects the following
 * structure:
 *
 * - \p first  \p argument : name of the solution component,
 * - \p second \p argument : boundary id,
 * - \p third  \p argument : value of the solution component.
 */
typedef std::map<   std::string , std::map<types::boundary_id, double>   > component_boundaryID_value_map;

/**
 * This namespace contains data and methods
 * that handle initial and boundary data of
 * a problem at hand.
 *
 * \author Valentin N. Zingan, 2013
 */

namespace FuelCell
{
    namespace InitialAndBoundaryData
    {
        
        /**
         * This function does some checkings on
         * its argument \p maps.
         */
        template<typename COMPONENT_xxxID_VALUE_MAP>
        const bool check(const std::vector< COMPONENT_xxxID_VALUE_MAP >& maps)
        {
            // --- types info ---
            
            const std::type_info& material = typeid(component_materialID_value_map);
            const std::type_info& boundary = typeid(component_boundaryID_value_map);
            
            const std::type_info& info = typeid(COMPONENT_xxxID_VALUE_MAP);
            
            // --- checkings ---
            
            // --- 1 ---
            
            if( maps.size() == 0 )
            {
                FcstUtilities::log << "The argument is empty" << std::endl;
                return false;
            }
            
            // --- 2 ---
            
            for(unsigned int i = 0; i < maps.size(); ++i)
                if( maps[i].empty() )
                {
                    FcstUtilities::log << "Outer map/maps of the argument is/are empty" << std::endl;
                    return false;
                }
                
                // --- 3 ---
                
                for(unsigned int i = 0; i < maps.size(); ++i)
                    for(typename COMPONENT_xxxID_VALUE_MAP::const_iterator iter  = maps[i].begin();
                        iter != maps[i].end();
                    ++iter)
                        if( iter->second.empty() )
                        {
                            FcstUtilities::log << "Inner map/maps of the argument is/are empty" << std::endl;
                            return false;
                        }
                        
                        // --- 4 ---
                        
                        if( info == material )
                            for(unsigned int i = 0; i < maps.size(); ++i)
                                for(typename COMPONENT_xxxID_VALUE_MAP::const_iterator iter  = maps[i].begin();
                                    iter != maps[i].end();
                                ++iter)
                                    if( iter->second.find(numbers::invalid_material_id) != iter->second.end() )
                                    {
                                        FcstUtilities::log << "Invalid material id/ids of the argument" << std::endl;
                                        return false;
                                    }
                                    
                                    if( info == boundary )
                                        for(unsigned int i = 0; i < maps.size(); ++i)
                                            for(typename COMPONENT_xxxID_VALUE_MAP::const_iterator iter  = maps[i].begin();
                                                iter != maps[i].end();
                                            ++iter)
                                                if( iter->second.find(numbers::invalid_boundary_id) != iter->second.end() )
                                                {
                                                    FcstUtilities::log << "Invalid boundary id/ids of the argument" << std::endl;
                                                    return false;
                                                }
                                                
                                                return true;
        }
        
        /**
         * This function makes piece wise constant initial data.
         * 
         * @param dst : Solution vector of size DOF. This is usually an FEVector. 
         * @param mapping : This parameter contains the geometrical map. It is stored in DoFApplication under *this->mapping
         * @param dof: DoFHandler object. It is stored in DoFApplication under *this->dof
         * @param system_management: This object stores the solution variable names and number as well as equation couplings. It should be in your application.
         * @param maps: This object is a vector of component_materialID_value_map. The component_materialID_value_map object is generated in the linear
         * applications and should be initialized with the EquationBase objects using a statement as follows:
         * @code
         * component_materialID_value_maps.push_back( ficks_transport_equation.get_component_materialID_value()    );
         *  component_materialID_value_maps.push_back( electron_transport_equation.get_component_materialID_value() );
         *  component_materialID_value_maps.push_back( proton_transport_equation.get_component_materialID_value()   );
         * @endcode
         * where electron_transport_equation is an object of #FuelCellShop::Equation::ElectronTransportEquation for example.
         * 
         */
        template<typename VECTOR, typename DH>
        void make_piece_wise_constant_initial_data(VECTOR&                                              dst,
                                                   const Mapping<DH::dimension, DH::space_dimension>&   mapping,
                                                   const DH&                                            dof,
                                                   const FuelCell::SystemManagement&                    system_management,
                                                   const std::vector< component_materialID_value_map >& maps)
        {
            // --- checking ---
            
            AssertThrow( FuelCell::InitialAndBoundaryData::check< component_materialID_value_map >(maps),
                         ExcMessage("The last argument in the FuelCell::InitialAndBoundaryData::make_piece_wise_constant_initial_data function "
                         "is wrong : check the previous message") );
            
            // --- main ---
            
            const unsigned int n_components = system_management.get_number_of_solution_names();
            AssertThrow( n_components == dof.get_fe().n_components(),
                         ExcDimensionMismatch( n_components , dof.get_fe().n_components() ) );
            
            for(unsigned int i = 0; i < maps.size(); ++i)
                for(component_materialID_value_map::const_iterator iter  = maps[i].begin(); iter != maps[i].end(); ++iter)
                    {
                        std::map<types::material_id, double> tmp = iter->second;
                        
                        for(std::map<types::material_id, double>::const_iterator iter2  = tmp.begin(); iter2 != tmp.end(); ++iter2)
                            {
                                if( iter2->second != _IS_NOT_CONSTANT_ )
                                {
                                    std::map< types::material_id, const Function<DH::space_dimension>* > function_map;
                                    
                                    const ConstantFunction<DH::space_dimension> constant_function(iter2->second, n_components);
                                    function_map[iter2->first] = &constant_function;
                                    
                                    std::vector<bool> component_mask(n_components, false);
                                    component_mask[system_management.solution_name_to_index(iter->first)] = true;
                                    
                                    VectorTools::interpolate_based_on_material_id( mapping,
                                                                                   dof,
                                                                                   function_map,
                                                                                   dst,
                                                                                   component_mask );
                                }
                            }
                    }
        }
        /**
         * @author: Mayank Sabharwal
         * This function makes piece wise constant Dirichlet BCs.
         * It is currently not being used in any appplication but has been templated based on the above functions.
         * 
         * @param dst : Solution vector of size DOF. This is usually an FEVector. 
         * @param mapping : This parameter contains the geometrical map. It is stored in DoFApplication under *this->mapping
         * @param dof: DoFHandler object. It is stored in DoFApplication under *this->dof
         * @param system_management: This object stores the solution variable names and number as well as equation couplings. It should be in your application.
         * @param maps: This object is a vector of component_materialID_value_map. The component_materialID_value_map object is generated in the linear
         * applications and should be initialized with the EquationBase objects using a statement as follows:
         */
        template<typename DH>
        void make_constant_DirichletBC_values(std::map<unsigned int, double>&                      dst,
                                              const Mapping<DH::dimension, DH::space_dimension>&   mapping,
                                              const DH&                                            dof,
                                              const FuelCell::SystemManagement&                    system_management,
                                              const std::vector< component_boundaryID_value_map >& maps)
        {
            // --- checking ---
            
            AssertThrow( FuelCell::InitialAndBoundaryData::check< component_boundaryID_value_map >(maps),
                         ExcMessage("The last argument in the FuelCell::InitialAndBoundaryData::make_constant_DirichletBC_values  function "
                         "is wrong : check the previous message") );
            
            // --- main ---
            
            const unsigned int n_components = system_management.get_number_of_solution_names();
            AssertThrow( n_components == dof.get_fe().n_components(),
                         ExcDimensionMismatch( n_components , dof.get_fe().n_components() ) );
            
            for(unsigned int i = 0; i < maps.size(); ++i)
                for(component_boundaryID_value_map::const_iterator iter  = maps[i].begin(); iter != maps[i].end(); ++iter)
                    {
                        std::map<types::boundary_id, double> tmp = iter->second;
                        
                        for(std::map<types::boundary_id, double>::const_iterator iter2  = tmp.begin(); iter2 != tmp.end(); ++iter2)
                            {
                                typename FunctionMap<DH::space_dimension>::type function_map;
                                
                                const ConstantFunction<DH::space_dimension> constant_function(iter2->second, n_components);
                                function_map[iter2->first] = &constant_function;

                                
                                std::vector<bool> component_mask(n_components, false);
                                component_mask[system_management.solution_name_to_index(iter->first)] = true;
                                
                                VectorTools::interpolate_boundary_values( mapping,
                                                                          dof,
                                                                          function_map,
                                                                          dst,
                                                                          component_mask );
                            }
                    }
        }
        
        /**
         * This function applies piece wise constant Dirichlet BCs.
         * 
         * @param dst : Solution vector of size DOF. This is usually an FEVector. 
         * @param mapping : This parameter contains the geometrical map. It is stored in DoFApplication under *this->mapping
         * @param dof: DoFHandler object. It is stored in DoFApplication under *this->dof
         * @param system_management: This object stores the solution variable names and number as well as equation couplings. It should be in your application.
         * @param maps: This object is a vector of component_materialID_value_map. The component_materialID_value_map object is generated in the linear
         * applications and should be initialized with the EquationBase objects using a statement as follows:
         * @code
         *  component_boundaryID_value_map.push_back( ficks_transport_equation.get_component_boundaryID_value()    );
         *  component_boundaryID_value_map.push_back( electron_transport_equation.get_component_boundaryID_value() );
         *  component_boundaryID_value_map.push_back( proton_transport_equation.get_component_boundaryID_value()   );
         * @endcode
         * where electron_transport_equation is an object of #FuelCellShop::Equation::ElectronTransportEquation for example.
         * 
         */
        template<typename VECTOR, typename DH>
        void apply_piece_wise_constant_DirichletBCs(VECTOR&                                              dst,
                                                    const Mapping<DH::dimension, DH::space_dimension>&   mapping,
                                                    const DH&                                            dof,
                                                    const FuelCell::SystemManagement&                    system_management,
                                                    const std::vector< component_boundaryID_value_map >& maps)
        {
            // --- checking ---
            /*
            AssertThrow( FuelCell::InitialAndBoundaryData::check< component_boundaryID_value_map >(maps),
                         ExcMessage("The last argument in the FuelCell::InitialAndBoundaryData::apply_piece_wise_constant_DirichletBCs function "
                         "is wrong : check the previous message") );
            */
            // --- main ---
            
            const unsigned int n_components = system_management.get_number_of_solution_names();
            AssertThrow( n_components == dof.get_fe().n_components(),
                         ExcDimensionMismatch( n_components , dof.get_fe().n_components() ) );
            
            for(unsigned int i = 0; i < maps.size(); ++i)
                for(component_boundaryID_value_map::const_iterator iter  = maps[i].begin(); iter != maps[i].end(); ++iter)
                    {
                        std::map<types::boundary_id, double> tmp = iter->second;
                        
                        for(std::map<types::boundary_id, double>::const_iterator iter2  = tmp.begin(); iter2 != tmp.end(); ++iter2)
                            {
                                if( iter2->second != _IS_NOT_CONSTANT_ )
                                {
                                    typename FunctionMap<DH::space_dimension>::type function_map;
                                    
                                    const ConstantFunction<DH::space_dimension> constant_function(iter2->second, n_components);
                                    function_map[iter2->first] = &constant_function;
                                    
                                    std::vector<bool> component_mask(n_components, false);
                                    component_mask[system_management.solution_name_to_index(iter->first)] = true;
                                    
                                    std::map<unsigned int, double> boundary_values;
                                    
                                    VectorTools::interpolate_boundary_values( mapping,
                                                                              dof,
                                                                              function_map,
                                                                              boundary_values,
                                                                              component_mask );
                                    
                                    for(std::map<unsigned int, double>::const_iterator iter3  = boundary_values.begin();
                                        iter3 != boundary_values.end(); ++iter3)
                                        {
                                            dst(iter3->first) = iter3->second;
                                        }
                                }
                            }
                    }
        }
        
        /**
         * This function makes zero boundary values.
         */
        template<typename DH>
        void make_zero_boundary_values(std::map<unsigned int, double>&                      dst,
                                       const Mapping<DH::dimension, DH::space_dimension>&   mapping,
                                       const DH&                                            dof,
                                       const FuelCell::SystemManagement&                    system_management,
                                       const std::vector< component_boundaryID_value_map >& maps)
        {
            // --- checking ---
            /*
            AssertThrow( FuelCell::InitialAndBoundaryData::check< component_boundaryID_value_map >(maps),
                         ExcMessage("The last argument in the FuelCell::InitialAndBoundaryData::make_zero_boundary_values function "
                         "is wrong : check the previous message") );
            */
            // --- main ---
            
            const unsigned int n_components = system_management.get_number_of_solution_names();
            AssertThrow( n_components == dof.get_fe().n_components(),
                         ExcDimensionMismatch( n_components , dof.get_fe().n_components() ) );
            
            for(unsigned int i = 0; i < maps.size(); ++i)
                for(component_boundaryID_value_map::const_iterator iter  = maps[i].begin(); iter != maps[i].end(); ++iter)
                    {
                        std::map<types::boundary_id, double> tmp = iter->second;
                        
                        for(std::map<types::boundary_id, double>::const_iterator iter2  = tmp.begin(); iter2 != tmp.end(); ++iter2)
                            {
                                typename FunctionMap<DH::space_dimension>::type function_map;
                                
                                const ZeroFunction<DH::space_dimension> zero_function(n_components);
                                function_map[iter2->first] = &zero_function;
                                
                                std::vector<bool> component_mask(n_components, false);
                                component_mask[system_management.solution_name_to_index(iter->first)] = true;
                                
                                VectorTools::interpolate_boundary_values( mapping,
                                                                          dof,
                                                                          function_map,
                                                                          dst,
                                                                          component_mask );
                            }
                    }
        }
        
        /**
         * This function applies
         * zero boundary values to
         * the linear system of equations.
         *
         * \warning The function name does not accurately describe the function
         * The function is no longer used anywhere. Suggest deprecation
         *
         */
        template<typename MATRIX, typename VECTOR, typename DH>
        void apply_zero_boundary_values_to_linear_system(MATRIX&                                              matrix,
                                                         VECTOR&                                              solution,
                                                         VECTOR&                                              rhs,
                                                         const Mapping<DH::dimension, DH::space_dimension>&   mapping,
                                                         const DH&                                            dof,
                                                         const FuelCell::SystemManagement&                    system_management,
                                                         const std::vector< component_boundaryID_value_map >& maps,
                                                         const bool&                                          repair_diagonal = false)
        {
            if( repair_diagonal )
            {
                SolverUtils::repair_diagonal( matrix,
                                              solution,
                                              rhs );
            }
            
            std::map<unsigned int, double> boundary_values;
            FuelCell::InitialAndBoundaryData::make_zero_boundary_values<DH>( boundary_values,
                                                                             mapping,
                                                                             dof,
                                                                             system_management,
                                                                             maps );
            MatrixTools::apply_boundary_values( boundary_values,
                                                matrix,
                                                solution,
                                                rhs );
        }
        
        /**
         * This class is a means to make variable
         * initial or boundary data.
         *
         * \author Valentin N. Zingan, 2013
         */
        
        template<int dim>
        class InitialOrBoundaryDataBase : public Function<dim>
        {
        public:
            
            ///@name Service functions
            //@{
            
            /**
             * This function calls the \p math_expression() function of
             * this class and has NOT to be overriden
             * in the derived classes.
             */
            virtual double value(const Point<dim>&  point,
                                 const unsigned int no_component = 0) const;
                                 
                                 //@}
                                 
        protected:
            
            ///@name Constructors, destructor, and initialization
            //@{
            
            /**
             * Constructor.
             *
             * \param data         - data of the application,
             * \param n_components - the total number of
             *                       the solution components.
             */
            InitialOrBoundaryDataBase(boost::shared_ptr<ApplicationData> data,
                                      const unsigned int                 n_components = 1);
            
            /**
             * Destructor.
             */
            ~InitialOrBoundaryDataBase();
            
            //@}
            
            ///@name Service functions
            //@{
            
            /**
             * This function implements the main functionality of
             * this class and has to be overriden
             * in the derived classes.
             */
            virtual double math_expression(const Point<dim>&  point,
                                           const unsigned int no_component = 0) const;
                                           
            //@}
                                           
            ///@name Minor functions
            //@{
            /**
             * This function is used to print out the name of another function
             * that has been declared in the scope of this class,
             * but not yet been implemented.
             */
            void print_caller_name(const std::string& caller_name) const;
            
            //@}
            
            //////////
            // DATA //
            //////////
            
            ///@name GENERIC DATA
            //@{
            
            /**
             * Data of the application serves to exchange
             * the information between YourApplication<dim>
             * and children of this class.
             */
            boost::shared_ptr<ApplicationData> data;
            
            //@}
            
        };
        
    } // InitialAndBoundaryData

} // FuelCell

#endif
