// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: system_management.h
// - Description: This class manages systems of equations
// - Developers: Valentin N. Zingan,    University of Alberta
//               Marc Secanell Gallart, University of Alberta
//
// ----------------------------------------------------------------------------

#ifndef _FCST_FUELCELL_SYSTEM_MANAGEMENT_H_
#define _FCST_FUELCELL_SYSTEM_MANAGEMENT_H_
//-- dealII 
#include <deal.II/base/parameter_handler.h>

//-- OpenFCST
#include <application_core/matrix_block.h>
#include <application_core/mesh_loop_info_objects.h>
#include <utils/logging.h>

//-- boost libraries
#include <boost/algorithm/string/join.hpp>

using namespace dealii;
using namespace FuelCell::ApplicationCore;

/**
 * The typedef for the map of cell or flux (DG FEM only)
 * couplings stored in the actual equation classes.
 *
 * The first std::string means the name of available FCST equation.
 *
 * The second std::string means the name of available FCST solution variable.
 */
typedef std::map<   std::string , std::map<std::string, DoFTools::Coupling>   > couplings_map;

/**
 * The enumeration containing the names of some of the available FCST solution variables
 * and their derivatives to be used in kinetics, materials, and other related classes.
 *
 * The enumeration is usually used in conjunction with FuelCellShop::SolutionVariable class.
 * For example, if the only available FCST solution variable is \p electronic_electrical_potential
 * then the respective piece of code will look like:
 *
 * @code
 * FuelCellShop::SolutionVariable solution;
 * (...)
 * if( solution.get_variablename() == VariableNames::electronic_electrical_potential )
 * {
 * (...)
 * }
 * @endcode
 *
 * In this case, if \p solution is initialized such that it contains \p electronic_electrical_potential
 * then the \p if statement will return \p TRUE.
 */
enum VariableNames
{
    nothing = 0,
    
    oxygen_molar_fraction,
    water_molar_fraction,
    hydrogen_molar_fraction,
    nitrogen_molar_fraction,
    protonic_electrical_potential,
    electronic_electrical_potential,
    membrane_water_content,
    temperature_of_solid_phase,
    temperature_of_fluid_phase,
    temperature_of_REV,
    liquid_water_saturation,
    capillary_pressure,
    
    current_density,
    CL_effectiveness,
    
    OH_coverage,
    O_coverage,
    
    oxygen_concentration,
    water_concentration,
    hydrogen_concentration,
    nitrogen_concentration,
    proton_concentration,
    
    total_pressure,
    
    density_species_1,
    velocity_species_1_X,
    velocity_species_1_Y,
    velocity_species_1_Z,
    
    density_species_2,
    velocity_species_2_X,
    velocity_species_2_Y,
    velocity_species_2_Z,
    
    density_species_3,
    velocity_species_3_X,
    velocity_species_3_Y,
    velocity_species_3_Z,
    
    density_species_4,
    velocity_species_4_X,
    velocity_species_4_Y,
    velocity_species_4_Z,
    
    density_species_5,
    velocity_species_5_X,
    velocity_species_5_Y,
    velocity_species_5_Z,
    
    single_fluid_pressure,
    single_fluid_velocity_X,
    single_fluid_velocity_Y,
    single_fluid_velocity_Z
};


/**
 * Enumeartion used to specify the desired solver to use.
 * 
 * This enumeration is used in FuelCellShop::Equation::EquationBase::set_solver(SolverName )
 * 
 * <h3> Usage </h3>
 * Simply create an instance of an enumeration, then use it in logical expressions.
 * @code
 * FuelCell::ApplicationCore::SolverName name_solver = FuelCell::ApplicationCore::SolverName::LINEAR 
 * 
 * if (name_solver == FuelCell::ApplicationCore::SolverName::LINEAR )
 * @endcode
 * 
 * \author M. Secanell, 2015
 */
namespace FuelCell
{
    namespace ApplicationCore
    {
        enum SolverName 
        {
            noSolver = 0,
            LINEAR, 
            NEWTON, 
            PICARD
        };
    }
}

/**
 * <h3> Theory </h3>
 * An enumeration used for expressing a type of reaction. Used by layers,
 * equations and catalyst materials for specifying a type of reaction
 * in a unified manner.
 *
 * <h3> Usage </h3>
 * Simply create an instance of an enumeration, then use it in logical expressions.
 *
 * @code
 * ReactionNames my_reaction = HOR;
 *
 * ...
 *
 * if (my_reaction == HOR)...
 * @endcode
 *
 * By Default the enumeration will be created as noReaction, therefore you must explicitly
 * set the reaction name. See unit tests in enumeration_test class.
 *
 * \author P. Wardlaw, 2014
 */
enum ReactionNames{
    noReaction = 0,
    
    ORR,
    HOR
};


namespace FuelCell
{
    
/**
 * IMPORTANT: Add all new solution variables
 * and equations here !
 *
 * This class contains
 * all the information on available FCST solution variables
 * and equations and feeds the actual equation classes with
 * different types of information on the linear system
 * structure.
 *
 * Available FCST solution variables:
 *
 * -  GROUP 1
 *
 * - "oxygen_molar_fraction"           \f$ \left( x_{O_2}          \right) \f$
 * - "water_molar_fraction"            \f$ \left( x_{H_2O}         \right) \f$
 * - "hydrogen_molar_fraction"         \f$ \left( x_{H_2}          \right) \f$
 * - "nitrogen_molar_fraction"         \f$ \left( x_{N_2}          \right) \f$
 * - "protonic_electrical_potential"   \f$ \left( \Phi_m           \right) \f$
 * - "electronic_electrical_potential" \f$ \left( \Phi_s           \right) \f$
 * - "membrane_water_content"          \f$ \left( \lambda          \right) \f$
 * - "temperature_of_solid_phase"      \f$ \left( T_{\text{solid}} \right) \f$
 * - "temperature_of_fluid_phase"      \f$ \left( T_{\text{fluid}} \right) \f$
 * - "temperature_of_REV"              \f$ \left( T_{\text{REV}}   \right) \f$
 * - "liquid_water_saturation"         \f$ \left( s                \right) \f$
 *
 * -  GROUP 2
 *
 * - "density_species_1"    \f$ \left( \rho_1  \right) \f$
 * - "velocity_species_1_X" \f$ \left( u_{1,x} \right) \f$
 * - "velocity_species_1_Y" \f$ \left( u_{1,y} \right) \f$
 * - "velocity_species_1_Z" \f$ \left( u_{1,z} \right) \f$
 *
 * - "density_species_2"    \f$ \left( \rho_2  \right) \f$
 * - "velocity_species_2_X" \f$ \left( u_{2,x} \right) \f$
 * - "velocity_species_2_Y" \f$ \left( u_{2,y} \right) \f$
 * - "velocity_species_2_Z" \f$ \left( u_{2,z} \right) \f$
 *
 * - "density_species_3"    \f$ \left( \rho_3  \right) \f$
 * - "velocity_species_3_X" \f$ \left( u_{3,x} \right) \f$
 * - "velocity_species_3_Y" \f$ \left( u_{3,y} \right) \f$
 * - "velocity_species_3_Z" \f$ \left( u_{3,z} \right) \f$
 *
 * - "density_species_4"    \f$ \left( \rho_4  \right) \f$
 * - "velocity_species_4_X" \f$ \left( u_{4,x} \right) \f$
 * - "velocity_species_4_Y" \f$ \left( u_{4,y} \right) \f$
 * - "velocity_species_4_Z" \f$ \left( u_{4,z} \right) \f$
 *
 * - "density_species_5"    \f$ \left( \rho_5  \right) \f$
 * - "velocity_species_5_X" \f$ \left( u_{5,x} \right) \f$
 * - "velocity_species_5_Y" \f$ \left( u_{5,y} \right) \f$
 * - "velocity_species_5_Z" \f$ \left( u_{5,z} \right) \f$
 *
 * -  GROUP 3
 *
 * - "single_fluid_pressure"   \f$ \left( p   \right) \f$
 * - "single_fluid_velocity_X" \f$ \left( u_x \right) \f$
 * - "single_fluid_velocity_Y" \f$ \left( u_y \right) \f$
 * - "single_fluid_velocity_Z" \f$ \left( u_z \right) \f$
 *
 * Available FCST equations:
 *
 * -  GROUP 1
 *
 * - "Ficks Transport Equation - oxygen"
 * - "Ficks Transport Equation - water"
 * - "Ficks Transport Equation - hydrogen"
 * - "Ficks Transport Equation - nitrogen"
 * - "Proton Transport Equation"
 * - "Electron Transport Equation"
 * - "Membrane Water Content Transport Equation"
 * - "Thermal Transport Equation"
 * - "Liquid Water Saturation Transport Equation"
 *
 * -  GROUP 2
 *
 * - "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - mass conservation - species 1"       \f$ \left( \rho_1, u_{1,x}, u_{1,y}, u_{1,z}               \right) \f$
 * - "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation X - species 1" \f$ \left( \{ \rho_s, u_{s,x}, u_{s,y}, u_{s,z} \}_{s=1}^5 \right) \f$
 * - "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation Y - species 1" \f$ \left( \{ \rho_s, u_{s,x}, u_{s,y}, u_{s,z} \}_{s=1}^5 \right) \f$
 * - "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation Z - species 1" \f$ \left( \{ \rho_s, u_{s,x}, u_{s,y}, u_{s,z} \}_{s=1}^5 \right) \f$
 *
 * - "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - mass conservation - species 2"       \f$ \left( \rho_2, u_{2,x}, u_{2,y}, u_{2,z}               \right) \f$
 * - "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation X - species 2" \f$ \left( \{ \rho_s, u_{s,x}, u_{s,y}, u_{s,z} \}_{s=1}^5 \right) \f$
 * - "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation Y - species 2" \f$ \left( \{ \rho_s, u_{s,x}, u_{s,y}, u_{s,z} \}_{s=1}^5 \right) \f$
 * - "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation Z - species 2" \f$ \left( \{ \rho_s, u_{s,x}, u_{s,y}, u_{s,z} \}_{s=1}^5 \right) \f$
 *
 * - "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - mass conservation - species 3"       \f$ \left( \rho_3, u_{3,x}, u_{3,y}, u_{3,z}               \right) \f$
 * - "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation X - species 3" \f$ \left( \{ \rho_s, u_{s,x}, u_{s,y}, u_{s,z} \}_{s=1}^5 \right) \f$
 * - "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation Y - species 3" \f$ \left( \{ \rho_s, u_{s,x}, u_{s,y}, u_{s,z} \}_{s=1}^5 \right) \f$
 * - "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation Z - species 3" \f$ \left( \{ \rho_s, u_{s,x}, u_{s,y}, u_{s,z} \}_{s=1}^5 \right) \f$
 *
 * - "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - mass conservation - species 4"       \f$ \left( \rho_4, u_{4,x}, u_{4,y}, u_{4,z}               \right) \f$
 * - "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation X - species 4" \f$ \left( \{ \rho_s, u_{s,x}, u_{s,y}, u_{s,z} \}_{s=1}^5 \right) \f$
 * - "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation Y - species 4" \f$ \left( \{ \rho_s, u_{s,x}, u_{s,y}, u_{s,z} \}_{s=1}^5 \right) \f$
 * - "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation Z - species 4" \f$ \left( \{ \rho_s, u_{s,x}, u_{s,y}, u_{s,z} \}_{s=1}^5 \right) \f$
 *
 * - "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - mass conservation - species 5"       \f$ \left( \rho_5, u_{5,x}, u_{5,y}, u_{5,z}               \right) \f$
 * - "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation X - species 5" \f$ \left( \{ \rho_s, u_{s,x}, u_{s,y}, u_{s,z} \}_{s=1}^5 \right) \f$
 * - "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation Y - species 5" \f$ \left( \{ \rho_s, u_{s,x}, u_{s,y}, u_{s,z} \}_{s=1}^5 \right) \f$
 * - "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation Z - species 5" \f$ \left( \{ \rho_s, u_{s,x}, u_{s,y}, u_{s,z} \}_{s=1}^5 \right) \f$
 *
 * -  GROUP 3
 *
 * - "Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - mass conservation"       \f$ \left(    u_x, u_y, u_z \right) \f$
 * - "Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation X" \f$ \left( p, u_x, u_y, u_z \right) \f$
 * - "Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation Y" \f$ \left( p, u_x, u_y, u_z \right) \f$
 * - "Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation Z" \f$ \left( p, u_x, u_y, u_z \right) \f$
 *
 * \author Valentin N. Zingan, 2012
 * \author Marc Secanell Gallart, 2012
 */

class SystemManagement : public Subscriptor
{
public:

///@name Constructors, destructor, and initialization
//@{

  /**
   * Constructor assembling
   *
   * - \p all_solution_names,
   * - \p all_equation_names.
   */
  SystemManagement();

  /**
   * Constructor assembling
   *
   * - \p all_solution_names,
   * - \p all_equation_names,
   * - \p block_info     (AppFrame structure) will be initialized by the AppFrame itself,
   * - \p cell_couplings (AppFrame structure) will be initialized by the \p make_cell_couplings function of this class called by a user defined application,
   * - \p flux_couplings (AppFrame structure) will be initialized by the \p make_flux_couplings function of this class called by a user defined application (DG FEM only).
   */
  SystemManagement(FuelCell::ApplicationCore::BlockInfo& block_info,
                   Table<2, DoFTools::Coupling>&     cell_couplings,
                   Table<2, DoFTools::Coupling>&     flux_couplings);

  /**
   * Destructor.
   */
 ~SystemManagement();

  /**
   * This function assembles
   *
   * - \p block_info     (AppFrame structure) will be initialized by the AppFrame itself,
   * - \p cell_couplings (AppFrame structure) will be initialized by the \p make_cell_couplings function of this class called by a user defined application,
   * - \p flux_couplings (AppFrame structure) will be initialized by the \p make_flux_couplings function of this class called by a user defined application (DG FEM only)
   *
   * if the empty constructor was used.
   */
  void initialize(FuelCell::ApplicationCore::BlockInfo& rblock_info,
                  Table<2, DoFTools::Coupling>&     rcell_couplings,
                  Table<2, DoFTools::Coupling>&     rflux_couplings)
  {
    block_info     = &rblock_info;
    cell_couplings = &rcell_couplings;
    flux_couplings = &rflux_couplings;
  }

  /**
   * Declare parameters.
   */
  void declare_parameters(ParameterHandler& param) const;

  /**
   * Initialize parameters.
   */
  void initialize(ParameterHandler& param);

//@}

///@name Service functions
//@{
// The # and ## preprocessor operators are available in C++ and ANSI/ISO C.
// The # operator causes a replacement-text token to be converted to a string surrounded by quotes.
const std::string& get_VariableNames(const VariableNames value)
{
    static std::map<VariableNames, std::string> VariableNames_strings;
    if (VariableNames_strings.size() == 0){
#define INSERT_ELEMENT(p) VariableNames_strings[p] = #p
        INSERT_ELEMENT(oxygen_molar_fraction);
        INSERT_ELEMENT(water_molar_fraction);
        INSERT_ELEMENT(hydrogen_molar_fraction);
#undef INSERT_ELEMENT
    }

    return VariableNames_strings[value];
}
  /**
   * This function returns
   * a boolean stating
   * whether an available FCST solution variable
   * with name \p name stored in \p solution_names.
   */
  const bool solution_in_userlist(const std::string& name) const;

  /**
   * This function returns
   * the index of
   * a user defined solution variable
   * with name \p name stored in \p solution_names.
   */
  const unsigned int solution_name_to_index(const std::string& name) const;

  /**
   * This function returns
   * the index of
   * a user defined solution variable
   * with name \p name stored in \p solution_names.
   */
  //const unsigned int solution_name_to_index(const VariableNames) const;

  /**
   * This function returns
   * the index of
   * a user defined equation
   * with name \p name stored in \p equation_names.
   */
  const unsigned int equation_name_to_index(const std::string& name) const;

  /**
   * This function returns
   * the index of
   * a system matrix block
   * based on
   *
   * - \p equation_name (line),
   * - \p solution_name (column),
   * - \p cell_couplings.
   *
   * \note A system matrix block
   * with
   * \p cell_couplings(\p equation_name, \p solution_name) = DoFTools::none
   * is missed.
   */
  const unsigned int matrix_block_index(const std::string& equation_name,
                                        const std::string& solution_name) const;
  
  /**
   * Check if, based on cell_couplings, this block must be assembled
   */
  inline bool matrix_block_index_exists(const std::string& equation_name,
                                        const std::string& solution_name) const
  {     
      const unsigned int line   = this->equation_name_to_index(equation_name);
      const unsigned int column = this->solution_name_to_index(solution_name);
      
      bool value = false;
      if ( (*cell_couplings)(line,column) != 0 )
          value = true;
      
      return value;
   }

  /**
   * This function fills out the* \p cell_couplings object
   * based on the information drawn from the actual equation classes in use.  
   */
  void make_cell_couplings(const std::vector<couplings_map>& src);

  /**
   * This function fills out the
   * \p flux_couplings (DG FEM only) object
   * based on the information
   * drawn from the actual equation classes
   * in use.
   */
  void make_flux_couplings(const std::vector<couplings_map>& src);

//@}

///@name Accessors and info
//@{

  /**
   * This function returns
   * \p all_solution_names.
   */
  const std::vector<std::string>& get_all_solution_names() const
  {
    return all_solution_names;
  }

  /**
   * This function returns
   * \p all_equation_names.
   */
  const std::vector<std::string>& get_all_equation_names() const
  {
    return all_equation_names;
  }

  /**
   * This function returns
   * \p solution_names.
   */
  const std::vector<std::string>& get_solution_names() const
  {
    return solution_names;
  }

  /**
   * This function returns
   * \p equation_names.
   */
  const std::vector<std::string>& get_equation_names() const
  {
    return equation_names;
  }

  /**
   * This function returns
   * \p n_solution_names.
   */
  const unsigned int& get_number_of_solution_names() const
  {
    return n_solution_names;
  }

  /**
   * This function prints out
   * the information on
   * available FCST solution variables
   * and equations
   * and what you have picked out.
   *
   * It also displays
   * the basic linear
   * system structure.
   */
  void print_system_info() const;

//@}

///@name Exceptions
//@{

  /**
   * Exception thrown when
   * a user defined solution variable
   * is not found among available FCST solution variables.
   */
  DeclException2(VariableNotFoundInFCSTVariables,
                 std::string,
                 std::string,
                 << "A "  << arg1 << " with name \"" << arg2 << "\" is not stored in available FCST solution variables");

  /**
   * Exception thrown when
   * a user defined equation
   * is not found among available FCST equations.
   */
  DeclException2(EquationNotFoundInFCSTEquations,
                 std::string,
                 std::string,
                 << "An " << arg1 << " with name \"" << arg2 << "\" is not stored in available FCST equations");

  /**
   * Exception thrown when
   * a \p name in \p solution_name_to_index function
   * is not found among user defined solution variables.
   */
  DeclException2(VariableNotFoundInUserVariables,
                 std::string,
                 std::string,
                 << "A "  << arg1 << " with name \"" << arg2 << "\" is not stored in user defined solution variables");

  /**
   * Exception thrown when
   * a \p name in \p equation_name_to_index function
   * is not found among user defined equations.
   */
  DeclException2(EquationNotFoundInUserEquations,
                 std::string,
                 std::string,
                 << "An " << arg1 << " with name \"" << arg2 << "\" is not stored in user defined equations");

  /**
   * Exception thrown when
   * a requested system matrix block
   * does not exist.
   *
   * \note Such a block has
   * \p cell_couplings(\p equation_name, \p solution_name) = DoFTools::none.
   */
  DeclException5(SystemMatrixBlockDoesNotExist,
                 std::string,
                 unsigned int,
                 unsigned int,
                 std::string,
                 std::string,
                 << "A "  << arg1 << " (" << arg2 << "," << arg3 << ") " << "does not exist, because the variable \"" << arg4 << "\" is not coupled with the equation \"" << arg5 << "\"");

//@}

  //////////
  // DATA //
  //////////

///@name External data pointers
//@{

  /**
   * Pointer to the external
   * YourApplication<dim>::block_info object.
   * Direct access is provided.
   */
  FuelCell::ApplicationCore::BlockInfo* block_info;

  /**
   * Pointer to the external
   * YourApplication<dim>::cell_couplings object.
   * Direct access is provided.
   */
  Table<2, DoFTools::Coupling>* cell_couplings;

  /**
   * Pointer to the external
   * YourApplication<dim>::flux_couplings (DG FEM only) object.
   * Direct access is provided.
   */
  Table<2, DoFTools::Coupling>* flux_couplings;

//@}

protected:

///@name Other functions
//@{

  /**
   * This function fills out \p all_solution_names.
   *
   * If you have a name of a new FCST solution variable, add it here.
   */
  void set_all_solution_names();

  /**
   * This function fills out \p all_equation_names.
   *
   * If you have a name of a new FCST equation, add it here.
   */
  void set_all_equation_names();

  /**
   * This function throws an exception if
   * a user defined solution variable
   * is not found among available FCST solution variables.
   */
  void check_solution_names() const;

  /**
   * This function throws an exception if
   * a user defined equation
   * is not found among available FCST equations.
   */
  void check_equation_names() const;

//@}

  //////////
  // DATA //
  //////////

///@name System properties
//@{
   /**
    *
    */
   //std::vector<VariableNames > variable_names;
  /**
   * Vector storing the names of available FCST solution variables.
   */
  std::vector<std::string> all_solution_names;

  /**
   * Vector storing the names of available FCST equations.
   */
  std::vector<std::string> all_equation_names;

  /**
   * Vector storing the names of user defined solution variables.
   */
  std::vector<std::string> solution_names;

  /**
   * Vector storing the names of user defined equations.
   */
  std::vector<std::string> equation_names;

  /**
   * The number of user defined solution variables.
   * The number of user defined equations
   * is implicitly assumed to be the same.
   */
  unsigned int n_solution_names;

  std::map<VariableNames, std::string> VariableNames_strings;
//@}

};

} // FuelCell

#endif