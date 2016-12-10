// ------------------------------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: fcst_variables.h
// - Description: This structure keeps a name of an FCST solution variable
//                and
//                values of this variable in the quadrature points of a mesh entity
// - Developers: Phil Wardlaw,       University of Alberta
//               Madhur Bhaiya,      University of Alberta
//               Valentin N. Zingan, University of Alberta
//
// ------------------------------------------------------------------------------------------------

#ifndef _FCST_FUELCELLSHOP_SOLUTION_VARIABLE_H_
#define _FCST_FUELCELLSHOP_SOLUTION_VARIABLE_H_
//-- C++ Standard libraries
#include <cmath>
#include <iostream>

//-- dealII
#include <deal.II/base/point.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/fe/fe_values.h>

//--OpenFCST
#include <application_core/system_management.h>

using namespace dealii;

namespace FuelCellShop
{
    /**
     * This structure is used to encapsulate data from constant values and variable solution data that is used
     * in Layer classes to compute effective transport properties. This structure is usually generated inside
     * Equation classes and then passed to the layer class as appropriate. The layer then harvests the appropriate
     * information.
     * 
     * This structure is used to store values for a particular solution variable, at all quadrature points in the cell. It is
     * best utilized while setting previous Newton iteration values from application to the layer classes etc. It is recommended that
     * those classes which require certain old solution values to compute effective properties etc. should contain this structure as their 
     * data member. This structure has four constructors. Default constructor doesn't set any value. It also sets the Boolean member #initialized to \p \b false. This
     * can be checked by using #is_initialized member function and is helpful to determine whether any values are stored or not inside the object. The other constructor
     * takes std::vector<double>* corresponding to the solution variable values and enumeration #VariableNames representing the name of the
     * solution variable. The #VariableNames of the solution variable can be accessed by get_variablename() method. In order
     * to access value at a particular quadrature point, <b>operator [ ]</b> is provided. Other two constructors are useful in the case when we set some default values
     * for a particular variable (which is not being solved in the current application).
     * 
     * <h3>Sample Usage:</h3>
     * @code
     *  // Inside the application, let's pass temperature vector from old solution to the layer class.
     *  // Old solution is accessed from FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo cell_info object.
     * 
     *  SolutionVariable temperature_old( &cell_info.values[last_iter][temp_index], VariableNames::temperature_of_REV );
     *  layer.set_solution( temperature_old );
     * 
     *  // Let us compute the temperature in Celsius at first quadrature point in the layer/kinetics etc. class:
     * 
     *  SolutionVariable temperature_cell;
     *  // Inside the layer/kinetics etc. class, temperature values are set into the temperature_cell structure using the set_solution method from the application.
     * 
     *  if (temperature_cell.is_initialized() && (temperature_cell.get_variablename() == VariableNames::temperature_of_REV)
     *  {
     *     double temp_celsius = temperature_cell[0] - 273.15;
     *  } 
     * @endcode
     * 
     * \TODO 
     * Can we create a "super" SolutionVariable object that stores all constant and variable data needed at each 
     * quadrature point? The amount of data needed would be provided by SystemManagement in combination with
     * Equation classes.
     * 
     * 
     * \author Philip Wardlaw
     * \author Madhur Bhaiya
     * \author Marc Secanell
     * \author Valentin N. Zingan
     */
    
    struct SolutionVariable
    {
        
        ///@name Initialization
        //@{
        
        /** Default Constructor */
        SolutionVariable()
        {
            data = NULL;
            initialized = false;
            initialized_default_data = false;
            initialized_data = false;
        }
        
        /**
         * Constructor for setting up the pointer to solution variable values and name of the solution variable.
         */
        SolutionVariable(const std::vector<double>* data_in, const VariableNames& name_in)
        {
            data = data_in;
            name = name_in;
            initialized = true;
            initialized_default_data = false;
            initialized_data = true;
        }
        
        /**
         * Constructor to initialize the solution variable values, taking a default value and size of the vector as an input argument.
         * \note This constructor is recommended when the solution variable values are not being solved for in the application, while
         * it is a required variable to compute certain properties. Then, it can be used to set some default values.
         */
        SolutionVariable(const double& value, const unsigned int& length, const VariableNames& name_in)
        {
            data = NULL;
            default_data = std::vector<double>(length, value);
            name = name_in;
            initialized = true;
            initialized_default_data = true;
            initialized_data = false;
        }
        
        /**
         * Constructor to initialize the solution variable values, taking values as an input vector argument.
         * 
         * \note This constructor is recommended when the solution variable values are not being solved for in the application, while
         * it is a required variable to compute certain properties. Then, it can be used to set some default values.
         */
        SolutionVariable(const std::vector<double>& data_in, const VariableNames& name_in)
        {
            data = NULL;
            default_data = data_in;
            name = name_in;
            initialized = true;
            initialized_default_data = true;
            initialized_data = false;
        }
        
        //@}
        
        ///@name Accessors
        //@{
        /**
         * Return a reference to the default_data stored in the class. The default_data is the data that was passed to the class
         * in the constructor.
         */
        const std::vector<double>& get_default_data() const
        {
            Assert( initialized_default_data, ExcMessage("default_data is not initialized") );
            return default_data;
        }
        
        /**
         * 
         */
        const std::vector<double>* get_data() const
        {
            Assert( initialized_data, ExcMessage("data is not initialized") );
            return data;
        }
        
        /**
         * Function to get the #VariableNames enumeration corresponding to this struct.
         */
        VariableNames get_variablename() const
        {
            Assert( initialized, ExcMessage("SolutionVariable not initialized !!!") );
            return name;
        }
        
        /**
         * Function to determine whether the structure is initialized or not.
         */
        const bool is_initialized() const
        {
            return initialized;
        }
        
        const bool is_default_data_initialized() const
        {
            return initialized_default_data;
        }
        
        const bool is_data_initialized() const
        {
            return initialized_data;
        }
        
        /**
         * Function to the length of the internal data element.
         */
        unsigned int size() const
        {
            unsigned int answer = 0;
            //If initialized and pointer is not NULL
            if (initialized && (data != NULL))
                answer = data->size();
            else if (initialized && (data == NULL))
                answer = default_data.size();
            
            return answer;
        }
        
        //@}
        
        ///@name Operators
        //@{
        
        /**
         * Operator to access the value at a particular quadrature point in the cell.
         * \note Solution values vector should be initialized before using this operator. Also, the input index should not be out of range.
         */
        const double& operator[](const unsigned int& i) const
        {
            Assert( initialized, ExcMessage("SolutionVariables struct is not initialized !!!") );
            
            if (data != NULL)
            {
                Assert( i < data->size(), ExcMessage("Index is out of range in operator[] for SolutionVariables struct.") );
                return data->at(i);
            }
            else if (data == NULL)
            {
                Assert( i < default_data.size(), ExcMessage("Index is out of range in operator[] for SolutionVariables struct.") );
                return default_data.at(i);
            }
        }
        
        /**
         * Operator \p SolutionVariable * \p double.
         * All \p default_data at quadrature points is multilpied by \p right operand.
         * The multiplication is performed on \p default_data only.
         * Therefore, to use this operator, \p default_data is initialized
         * and \p data is NOT.
         */
        friend SolutionVariable operator* (const SolutionVariable& left,
                                           const double&           right)
        {
            Assert( left.is_default_data_initialized(), ExcMessage("SolutionVariable left operand default_data is not initialized") );
            
            std::vector<double> tmp(left.size());
            
            const VariableNames name = left.get_variablename();
            
            if( left.is_default_data_initialized() && !left.is_data_initialized() )
            {
                for(unsigned int q = 0; q < tmp.size(); ++q)
                    tmp[q] = left.get_default_data().at(q);
                
                for(unsigned int q = 0; q < tmp.size(); ++q)
                    tmp[q] *= right;
                
                return SolutionVariable(tmp,
                                        name);
            }
            else if( !left.is_default_data_initialized() && left.is_data_initialized() )
            {
                AssertThrow( false, ExcNotImplemented() );
            }
            else
            {
                AssertThrow( false, ExcInternalError() );
            }
        }
        
        /**
         * Operator \p double * \p SolutionVariable.
         * All \p default_data at quadrature points is multilpied by \p left operand.
         * The multiplication is performed on \p default_data only.
         * Therefore, to use this operator, \p default_data is initialized
         * and \p data is NOT.
         */
        friend SolutionVariable operator* (const double&           left,
                                           const SolutionVariable& right)
        {
            Assert( right.is_default_data_initialized(), ExcMessage("SolutionVariable right operand default_data is not initialized") );
            
            std::vector<double> tmp(right.size());
            
            const VariableNames name = right.get_variablename();
            
            if( right.is_default_data_initialized() && !right.is_data_initialized() )
            {
                for(unsigned int q = 0; q < tmp.size(); ++q)
                    tmp[q] = right.get_default_data().at(q);
                
                for(unsigned int q = 0; q < tmp.size(); ++q)
                    tmp[q] *= left;
                
                return SolutionVariable(tmp,
                                        name);
            }
            else if( !right.is_default_data_initialized() && right.is_data_initialized() )
            {
                AssertThrow( false, ExcNotImplemented() );
            }
            else
            {
                AssertThrow( false, ExcInternalError() );
            }
        }
        
        /**
         * Operator \p SolutionVariable / \p double.
         * All \p default_data at quadrature points is divided by \p right operand.
         * The division is performed on \p default_data only.
         * Therefore, to use this operator, \p default_data is initialized
         * and \p data is NOT.
         */
        friend SolutionVariable operator/ (const SolutionVariable& left,
                                           const double&           right)
        {
            Assert( left.is_default_data_initialized(), ExcMessage("SolutionVariable left operand default_data is not initialized") );
            
            std::vector<double> tmp(left.size());
            
            const VariableNames name = left.get_variablename();
            
            if( left.is_default_data_initialized() && !left.is_data_initialized() )
            {
                for(unsigned int q = 0; q < tmp.size(); ++q)
                    tmp[q] = left.get_default_data().at(q);
                
                for(unsigned int q = 0; q < tmp.size(); ++q)
                    tmp[q] /= right;
                
                return SolutionVariable(tmp,
                                        name);
            }
            else if( !left.is_default_data_initialized() && left.is_data_initialized() )
            {
                AssertThrow( false, ExcNotImplemented() );
            }
            else
            {
                AssertThrow( false, ExcInternalError() );
            }
        }
        
        
        //@}
        
    private:
        
        ///@name DATA
        //@{
        
        /**
         * Constant data. 
         * 
         * \warning Only one of either default_data or data is initialize. The selection occurs depending on the constructor used.
         */
        std::vector<double> default_data;
        
        /**
         * Data in quadrature points of a mesh entity.
         * 
         * \warning Only one of either default_data or data is initialize. The selection occurs depending on the constructor used.
         * 
         */
        const std::vector<double>* data;
        
        /**
         * FCST variable name stored in VariableNames enumeration.
         */
        VariableNames name;
        
        /**
         * \p true if either \p default_data
         * or \p data is initialized.
         */
        bool initialized;
        
        /**
         * \p true if \p default_data
         * is initialized.
         */
        bool initialized_default_data;
        
        /**
         * \p true if \p data
         * is initialized.
         */
        bool initialized_data;
        
        //@}
        
    };
    
    ///@name Unary Predicate Functions for SolutionVariable structure
    //@{
    
    /**
     * Unary Predicate to return true if a SolutionVariable object belongs to #protonic_electrical_potential.
     */
    static bool is_phiM(const SolutionVariable& sol_var)
    {
        return (sol_var.get_variablename() == protonic_electrical_potential);
    }
    
    /**
     * Unary Predicate to return true if a SolutionVariable object belongs to #electronic_electrical_potential.
     */
    static bool is_phiS(const SolutionVariable& sol_var)
    {
        return (sol_var.get_variablename() == electronic_electrical_potential);
    }
    
    //@}
    
    
    /**
     * 
     * @brief Convenient storage object for SolutionVariables
     *
     * Used for storing and passing sets of solutions.
     *
     * Inherits of std::map< VariableNames, SolutionVariable> privately,
     * as to no expose operator[] since we wish to check all map inserts
     *
     * <h3> Usage details</h3>
     * Simply add Solutions using push_back(), and access them as normal using the
     * at() function. To determine if a SolutionMap has a certain SolutionVariable
     * corresponding to a VariableNames use the has() function.
     *
     * @code
     *
     * //Create a SolutionVariable
     * FuelCellShop::SolutionVariable phi_s(1.5, 8, electronic_electrical_potential);
     *
     * //Create an object of SolutionMap
     * FuelCellShop::SolutionMap map;
     *
     * //Add the solution variable to the map
     * map.push_back(phi_s);
     *
     * //Adding the same solution again will throw an exception
     * //map.push_back(phi_s) here will throw std::runtime_error
     *
     * //Access the solution variable from the map
     * double first_value = map.at(electronic_electrical_potential)[0];
     *
     * //Query the solution map to see if it has a solution for a given variable type
     * if(map.has(oxygen_molar_fraction) //false
     *    ...
     *
     * @endcode
     *
     * @author Philip Wardlaw
     * @date 2014
     */
    class SolutionMap : private std::map< VariableNames, SolutionVariable>
    {
    public:
        
        /**
         * Public function for adding SolutionVariable, uses the VariableNames stored within the SolutionVariable as a key.
         *
         * @warning
         * SolutionMap will not accept multiple definitions of a solution and will throw std::runtime_error if the user tries
         * to store two solutions of corresponding to the same VariableNames type.
         *
         * To circumvent this either clear the map using clear(), or to replace a single element first using erase(key) to
         * clear existing SolutionVariable corresponding to the same VariableNames type or pop(key) to get the SolutionVariable
         * whilst simulataneously erasing it.
         * @endwarning
         *
         */
        void push_back(const SolutionVariable& a){
            
            
            //Check that the key hasn't already been entered
            if (this->count(a.get_variablename()) > 0)
                throw std::runtime_error("You have already added a SolutionVariable of corresponding VariableNames type to SolutionMap");
            
            
            //add the variable
            this->operator [](a.get_variablename()) = a;
        }
        
        /**
         * Expose std::map<VariableNames, SolutionVariable>::at() interface publicly
         */
        SolutionVariable& at(VariableNames key){
            return std::map<VariableNames, SolutionVariable>::at(key);
        }
        
        /**
         * Expose std::map<VariableNames, SolutionVariable>::clear() interface publicly
         */
        void clear(){
            std::map<VariableNames, SolutionVariable>::clear();
        }
        
        /**
         * Expose std::map<VariableNames, SolutionVariable>::erase() interface publicly
         */
        void erase(const VariableNames& v){
            std::map<VariableNames, SolutionVariable>::erase(v);
        }
        
        
        /**
         * Find if a solution corresponding VariableNames type exist inside map.
         */
        bool has(const VariableNames& v) const{
            bool answer = false;
            
            if(this->count(v) >0)
                answer = true;
            
            return answer;
        }
        
        /**
         * Returns and entry whilst removing it from the list
         */        
        SolutionVariable pop(const VariableNames& v){
            SolutionVariable temp = this->at(v);
            std::map<VariableNames, SolutionVariable>::erase(v);
            return temp;
        }
    };
}
#endif