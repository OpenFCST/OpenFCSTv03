// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: system_management.cc
// - Description: This class manages systems of equations
// - Developers: Valentin N. Zingan,    University of Alberta
//               Marc Secanell Gallart, University of Alberta
//
// ----------------------------------------------------------------------------

#include <application_core/system_management.h>

       //////////////////////////////////////////////////
       //////////////////////////////////////////////////
       // CONSTRUCTORS, DESTRUCTOR, AND INITIALIZATION //
       //////////////////////////////////////////////////
       //////////////////////////////////////////////////

// ---             ---
// --- Constructor ---
// ---             ---

FuelCell::SystemManagement::SystemManagement()
{
  this->set_all_solution_names();
  this->set_all_equation_names();


    if (VariableNames_strings.size() == 0){
        #define INSERT_ELEMENT(p) VariableNames_strings[p] = #p
        INSERT_ELEMENT(nothing);
        INSERT_ELEMENT(oxygen_molar_fraction);
        INSERT_ELEMENT(water_molar_fraction);
        INSERT_ELEMENT(hydrogen_molar_fraction);
        INSERT_ELEMENT(nitrogen_molar_fraction);
        INSERT_ELEMENT(protonic_electrical_potential);
        INSERT_ELEMENT(electronic_electrical_potential);
        INSERT_ELEMENT(membrane_water_content);
        INSERT_ELEMENT(temperature_of_solid_phase);
        INSERT_ELEMENT(temperature_of_fluid_phase);
        INSERT_ELEMENT(temperature_of_REV);
        INSERT_ELEMENT(liquid_water_saturation);

        INSERT_ELEMENT(current_density);
        INSERT_ELEMENT(CL_effectiveness);

        INSERT_ELEMENT(OH_coverage);
        INSERT_ELEMENT(O_coverage);

        INSERT_ELEMENT(oxygen_concentration);
        INSERT_ELEMENT(water_concentration);
        INSERT_ELEMENT(hydrogen_concentration);
        INSERT_ELEMENT(nitrogen_concentration);
        INSERT_ELEMENT(proton_concentration);

        INSERT_ELEMENT(total_pressure);

        INSERT_ELEMENT(density_species_1);
        INSERT_ELEMENT(velocity_species_1_X);
        INSERT_ELEMENT(velocity_species_1_Y);
        INSERT_ELEMENT(velocity_species_1_Z);

        INSERT_ELEMENT(density_species_2);
        INSERT_ELEMENT(velocity_species_2_X);
        INSERT_ELEMENT(velocity_species_2_Y);
        INSERT_ELEMENT(velocity_species_2_Z);

        INSERT_ELEMENT(density_species_3);
        INSERT_ELEMENT(velocity_species_3_X);
        INSERT_ELEMENT(velocity_species_3_Y);
        INSERT_ELEMENT(velocity_species_3_Z);

        INSERT_ELEMENT(density_species_4);
        INSERT_ELEMENT(velocity_species_4_X);
        INSERT_ELEMENT(velocity_species_4_Y);
        INSERT_ELEMENT(velocity_species_4_Z);

        INSERT_ELEMENT(density_species_5);
        INSERT_ELEMENT(velocity_species_5_X);
        INSERT_ELEMENT(velocity_species_5_Y);
        INSERT_ELEMENT(velocity_species_5_Z);

        INSERT_ELEMENT(single_fluid_pressure);
        INSERT_ELEMENT(single_fluid_velocity_X);
        INSERT_ELEMENT(single_fluid_velocity_Y);
        INSERT_ELEMENT(single_fluid_velocity_Z);
        #undef INSERT_ELEMENT
    }
}

// ---             ---
// --- Constructor ---
// ---             ---

FuelCell::SystemManagement::SystemManagement(FuelCell::ApplicationCore::BlockInfo& block_info,
                                             Table<2, DoFTools::Coupling>&     cell_couplings,
                                             Table<2, DoFTools::Coupling>&     flux_couplings)
:
block_info(&block_info),
cell_couplings(&cell_couplings),
flux_couplings(&flux_couplings)
{
  this->set_all_solution_names();
  this->set_all_equation_names();
}

// ---            ---
// --- Destructor ---
// ---            ---

FuelCell::SystemManagement::~SystemManagement()
{ }

// ---                    ---
// --- declare_parameters ---
// ---                    ---

void
FuelCell::SystemManagement::declare_parameters(ParameterHandler& param) const
{
  param.enter_subsection("System management");
  {
    param.declare_entry("Number of solution variables",
                        "1",
                         Patterns::Integer(),
                        "Number of variables that you would like to solve for. You need to specify the name of the Solution variables and Equations in the sections below.");

    std::string joinedString = "no_name|"+boost::algorithm::join(all_solution_names, "|");
    
    param.enter_subsection("Solution variables");
    {
      const unsigned int n_sol_vars_max = 11;
      for(unsigned int index = 1; index <= n_sol_vars_max; ++index)
      {
             std::ostringstream streamOut;
             streamOut << index;
             const std::string name = "Solution variable " + streamOut.str();

             param.declare_entry(name.c_str(),
                                "no_name",
                                 Patterns::Selection(joinedString),
                                "Name of the solution variables we are solving for. These names depend on the application.");
      }
    }
    param.leave_subsection();

    joinedString = "no_name|"+boost::algorithm::join(all_equation_names, "|");
    
    param.enter_subsection("Equations");
    {
      const unsigned int n_eqs_max = 11;
      for(unsigned int index = 1; index <= n_eqs_max; ++index)
      {
             std::ostringstream streamOut;
             streamOut << index;
             const std::string name = "Equation " + streamOut.str();

             param.declare_entry(name.c_str(),
                                "no_name",
                                 Patterns::Selection(joinedString),
                                "Name of the equations you would like to solve for. The names depend on the application. Some equations are necessary");
      }
    }
    param.leave_subsection();
  }
  param.leave_subsection();
}

// ---            ---
// --- initialize ---
// ---            ---

void
FuelCell::SystemManagement::initialize(ParameterHandler& param)
{
    // clear data first (in case initialize is called twice)
    solution_names.clear();
    equation_names.clear();

    param.enter_subsection("System management");
    {
        n_solution_names = param.get_integer("Number of solution variables");

        param.enter_subsection("Solution variables");
        {
            for(unsigned int index = 1; index <= n_solution_names; ++index)
            {
                std::ostringstream streamOut;
                streamOut << index;
                const std::string name = "Solution variable " + streamOut.str();

                solution_names.push_back( param.get(name.c_str()) );
            }
        }
        param.leave_subsection();

        param.enter_subsection("Equations");
        {
            for(unsigned int index = 1; index <= n_solution_names; ++index)
            {
                std::ostringstream streamOut;
                streamOut << index;
                const std::string name = "Equation " + streamOut.str();

                equation_names.push_back( param.get(name.c_str()) );
            }
        }
        param.leave_subsection();
    }
    param.leave_subsection();

    this->check_solution_names();
    this->check_equation_names();
}

       ///////////////////////
       ///////////////////////
       // SERVICE FUNCTIONS //
       ///////////////////////
       ///////////////////////

// ---                      ---
// --- solution_in_userlist ---
// ---                      ---

const bool
FuelCell::SystemManagement::solution_in_userlist(const std::string& name) const
{
  const std::vector<std::string>::const_iterator iter = std::find(solution_names.begin(),
                                                                  solution_names.end(),
                                                                  name);
  if( iter != solution_names.end() )
    return true;
  else
    return false;
}

// ---                        ---
// --- solution_name_to_index ---
// ---                        ---

const unsigned int
FuelCell::SystemManagement::solution_name_to_index(const std::string& name) const
{
  const std::vector<std::string>::const_iterator iter = std::find(solution_names.begin(),
                                                                  solution_names.end(),
                                                                  name);

  AssertThrow( iter != solution_names.end(), VariableNotFoundInUserVariables("variable", name) );

  return iter - solution_names.begin();
}

/*
const unsigned int
FuelCell::SystemManagement::solution_name_to_index(const VariableNames) const
{

}
*/

// ---                        ---
// --- equation_name_to_index ---
// ---                        ---

const unsigned int
FuelCell::SystemManagement::equation_name_to_index(const std::string& name) const
{
  const std::vector<std::string>::const_iterator iter = std::find(equation_names.begin(),
                                                                  equation_names.end(),
                                                                  name);

  AssertThrow( iter != equation_names.end(), EquationNotFoundInUserEquations("equation", name) );

  return iter - equation_names.begin();
}

// ---                    ---
// --- matrix_block_index ---
// ---                    ---

const unsigned int
FuelCell::SystemManagement::matrix_block_index(const std::string& equation_name,
                                               const std::string& solution_name) const
{
  const unsigned int line   = this->equation_name_to_index(equation_name);
  const unsigned int column = this->solution_name_to_index(solution_name);

  AssertThrow( (*cell_couplings)(line,column) != 0, SystemMatrixBlockDoesNotExist("requested system matrix block",
                                                                                   line,
                                                                                   column,
                                                                                   solution_name,
                                                                                   equation_name) );

  unsigned int number = 0;

  for(unsigned int eq = 0; eq < equation_names.size(); ++eq)
    for(unsigned int var = 0; var < solution_names.size(); ++var)
    {
      if( eq == line && var == column )
        return number;

      if( (*cell_couplings)(eq,var) != 0 )
        number++;
    }
}

// ---                     ---
// --- make_cell_couplings ---
// ---                     ---

void
FuelCell::SystemManagement::make_cell_couplings(const std::vector<couplings_map>& src)
{
    if( src.size() == 0 )
    {
        FcstUtilities::log << "FuelCell::SystemManagement::make_cell_couplings function: The argument is empty" << std::endl;
        return;
    }
    
    for(unsigned int i = 0; i < src.size(); ++i)
        if( src[i].empty() )
        {
            FcstUtilities::log << "FuelCell::SystemManagement::make_cell_couplings function: Outer map/maps of the argument is/are empty" << std::endl;
            return;
        }
        
    for(unsigned int i = 0; i < src.size(); ++i)
        for(typename couplings_map::const_iterator iter  = src[i].begin(); iter != src[i].end(); ++iter)
        {
            std::map<std::string, DoFTools::Coupling> tmp = iter->second;
            if( tmp.empty() )
            {
                FcstUtilities::log << "FuelCell::SystemManagement::make_cell_couplings function: Inner map/maps of the argument is/are empty" << std::endl;
                return;
            }
        }
        
    (*cell_couplings).reinit(equation_names.size(), solution_names.size());
            
    for(unsigned int index = 0; index < src.size(); ++index)
    {
        couplings_map cell_couplings_i = src[index];
        couplings_map::const_iterator iter;
        
        for( iter = cell_couplings_i.begin(); iter != cell_couplings_i.end(); ++iter )
        {
            std::map<std::string, DoFTools::Coupling> int_map = iter->second;
            std::map<std::string, DoFTools::Coupling>::const_iterator int_iter;
            
            for( int_iter = int_map.begin(); int_iter != int_map.end(); ++int_iter )
            {
                const unsigned int line   = this->equation_name_to_index(iter->first);
                const unsigned int column = this->solution_name_to_index(int_iter->first);
                
                (*cell_couplings)(line,column) = int_iter->second;
            }
        }
    }
}

// ---                     ---
// --- make_flux_couplings ---
// ---                     ---

void
FuelCell::SystemManagement::make_flux_couplings(const std::vector<couplings_map>& src)
{
  if( src.size() == 0 )
  {
         FcstUtilities::log << "FuelCell::SystemManagement::make_flux_couplings function: The argument is empty" << std::endl;
         return;
  }

  for(unsigned int i = 0; i < src.size(); ++i)
         if( src[i].empty() )
         {
                FcstUtilities::log << "FuelCell::SystemManagement::make_flux_couplings function: Outer map/maps of the argument is/are empty" << std::endl;
                return;
         }

  for(unsigned int i = 0; i < src.size(); ++i)
         for(typename couplings_map::const_iterator iter  = src[i].begin();
                                                    iter != src[i].end();
                                                  ++iter)
         {
                std::map<std::string, DoFTools::Coupling> tmp = iter->second;
                if( tmp.empty() )
                {
                       FcstUtilities::log << "FuelCell::SystemManagement::make_flux_couplings function: Inner map/maps of the argument is/are empty" << std::endl;
                       return;
                }
         }

  (*flux_couplings).reinit(equation_names.size(), solution_names.size());

  for(unsigned int index = 0; index < src.size(); ++index)
  {
    couplings_map cell_couplings_i = src[index];
    couplings_map::const_iterator iter;

    for( iter = cell_couplings_i.begin(); iter != cell_couplings_i.end(); ++iter )
    {
      std::map<std::string, DoFTools::Coupling> int_map = iter->second;
      std::map<std::string, DoFTools::Coupling>::const_iterator int_iter;

      for( int_iter = int_map.begin(); int_iter != int_map.end(); ++int_iter )
      {
               const unsigned int line   = this->equation_name_to_index(iter->first);
               const unsigned int column = this->solution_name_to_index(int_iter->first);

               (*flux_couplings)(line,column) = int_iter->second;
      }
    }
  }
}

       /////////////////////
       /////////////////////
       // OTHER FUNCTIONS //
       /////////////////////
       /////////////////////

// ---                        ---
// --- set_all_solution_names ---
// ---                        ---

void
FuelCell::SystemManagement::set_all_solution_names()
{
    all_solution_names.push_back("test_var");
    // GROUP 1
    
    all_solution_names.push_back("oxygen_molar_fraction");
    all_solution_names.push_back("water_molar_fraction");
    all_solution_names.push_back("hydrogen_molar_fraction");
    all_solution_names.push_back("nitrogen_molar_fraction");
    all_solution_names.push_back("protonic_electrical_potential");
    all_solution_names.push_back("electronic_electrical_potential");
    all_solution_names.push_back("membrane_water_content");
    all_solution_names.push_back("temperature_of_solid_phase");
    all_solution_names.push_back("temperature_of_fluid_phase");
    all_solution_names.push_back("temperature_of_REV");
    all_solution_names.push_back("liquid_water_saturation");
    all_solution_names.push_back("capillary_pressure");
    
    // GROUP 2
    
    all_solution_names.push_back("density_species_1");
    all_solution_names.push_back("velocity_species_1_X");
    all_solution_names.push_back("velocity_species_1_Y");
    all_solution_names.push_back("velocity_species_1_Z");
    
    all_solution_names.push_back("density_species_2");
    all_solution_names.push_back("velocity_species_2_X");
    all_solution_names.push_back("velocity_species_2_Y");
    all_solution_names.push_back("velocity_species_2_Z");
    
    all_solution_names.push_back("density_species_3");
    all_solution_names.push_back("velocity_species_3_X");
    all_solution_names.push_back("velocity_species_3_Y");
    all_solution_names.push_back("velocity_species_3_Z");
    
    all_solution_names.push_back("density_species_4");
    all_solution_names.push_back("velocity_species_4_X");
    all_solution_names.push_back("velocity_species_4_Y");
    all_solution_names.push_back("velocity_species_4_Z");
    
    all_solution_names.push_back("density_species_5");
    all_solution_names.push_back("velocity_species_5_X");
    all_solution_names.push_back("velocity_species_5_Y");
    all_solution_names.push_back("velocity_species_5_Z");
    
    // GROUP 3
    
    all_solution_names.push_back("single_fluid_pressure");
    all_solution_names.push_back("single_fluid_velocity_X");
    all_solution_names.push_back("single_fluid_velocity_Y");
    all_solution_names.push_back("single_fluid_velocity_Z");
    
    
}

// ---                        ---
// --- set_all_equation_names ---
// ---                        ---

void
FuelCell::SystemManagement::set_all_equation_names()
{
    // GROUP 1
    all_equation_names.push_back("Test Equation");
    all_equation_names.push_back("Ficks Transport Equation - oxygen");
    all_equation_names.push_back("Ficks Transport Equation - water");
    all_equation_names.push_back("Ficks Transport Equation - hydrogen");
    all_equation_names.push_back("Ficks Transport Equation - nitrogen");
    all_equation_names.push_back("Proton Transport Equation");
    all_equation_names.push_back("Electron Transport Equation");
    all_equation_names.push_back("Membrane Water Content Transport Equation");
    all_equation_names.push_back("Thermal Transport Equation");
    all_equation_names.push_back("Liquid Water Saturation Transport Equation");
    all_equation_names.push_back("Liquid Water Capillary Transport Equation");
    all_equation_names.push_back("Liquid Water Cathode Capillary Transport Equation");
    all_equation_names.push_back("CapillaryTesting");
    
    // GROUP 2
    
    all_equation_names.push_back("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - mass conservation - species 1");
    all_equation_names.push_back("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation X - species 1");
    all_equation_names.push_back("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation Y - species 1");
    all_equation_names.push_back("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation Z - species 1");
    
    all_equation_names.push_back("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - mass conservation - species 2");
    all_equation_names.push_back("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation X - species 2");
    all_equation_names.push_back("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation Y - species 2");
    all_equation_names.push_back("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation Z - species 2");
    
    all_equation_names.push_back("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - mass conservation - species 3");
    all_equation_names.push_back("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation X - species 3");
    all_equation_names.push_back("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation Y - species 3");
    all_equation_names.push_back("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation Z - species 3");
    
    all_equation_names.push_back("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - mass conservation - species 4");
    all_equation_names.push_back("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation X - species 4");
    all_equation_names.push_back("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation Y - species 4");
    all_equation_names.push_back("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation Z - species 4");
    
    all_equation_names.push_back("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - mass conservation - species 5");
    all_equation_names.push_back("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation X - species 5");
    all_equation_names.push_back("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation Y - species 5");
    all_equation_names.push_back("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - momentum conservation Z - species 5");
    
    // GROUP 3
    
    all_equation_names.push_back("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - mass conservation");
    all_equation_names.push_back("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation X");
    all_equation_names.push_back("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation Y");
    all_equation_names.push_back("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation Z");
}

// ---                      ---
// --- check_solution_names ---
// ---                      ---

void
FuelCell::SystemManagement::check_solution_names() const
{
    for(unsigned int index = 0; index < solution_names.size(); ++index)
    {
        const std::string name = solution_names[index];
        
        const std::vector<std::string>::const_iterator iter = std::find(all_solution_names.begin(),
                                                                        all_solution_names.end(),
                                                                        name);
        
        AssertThrow( iter != all_solution_names.end(), VariableNotFoundInFCSTVariables("variable", name) );
    }
}

// ---                      ---
// --- check_equation_names ---
// ---                      ---

void
FuelCell::SystemManagement::check_equation_names() const
{
    for(unsigned int index = 0; index < equation_names.size(); ++index)
    {
        const std::string name = equation_names[index];
        
        const std::vector<std::string>::const_iterator iter = std::find(all_equation_names.begin(),
                                                                        all_equation_names.end(),
                                                                        name);
        
        AssertThrow( iter != all_equation_names.end(), EquationNotFoundInFCSTEquations("equation", name) );
    }
}

////////////////////////
////////////////////////
// ACCESSORS AND INFO //
////////////////////////
////////////////////////

// ---                   ---
// --- print_system_info ---
// ---                   ---

void
FuelCell::SystemManagement::print_system_info() const
{
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "------------------------------";
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "Available FCST solution variables:" << std::endl;
    FcstUtilities::log << std::endl;
    
    for(unsigned int index = 0; index < all_solution_names.size(); ++index)
        FcstUtilities::log << all_solution_names[index] << std::endl;
    
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "Available FCST equations:" << std::endl;
    FcstUtilities::log << std::endl;
    
    for(unsigned int index = 0; index < all_equation_names.size(); ++index)
        FcstUtilities::log << all_equation_names[index] << std::endl;
    
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "User defined solution variables:" << std::endl;
    FcstUtilities::log << std::endl;
    
    for(unsigned int index = 0; index < solution_names.size(); ++index)
        FcstUtilities::log << solution_names[index] << std::endl;
    
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "User defined equations:" << std::endl;
    FcstUtilities::log << std::endl;
    
    for(unsigned int index = 0; index < equation_names.size(); ++index)
        FcstUtilities::log << equation_names[index] << std::endl;
    
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "Number of blocks: " << block_info->global.size();
    FcstUtilities::log << std::endl;
    
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "Base element: " << std::endl;
    FcstUtilities::log << std::endl;
    
    for(unsigned int index = 0; index < block_info->base_element.size(); ++index)
        FcstUtilities::log << "No block: " << index << ", Name of block: " << solution_names[index] << ", No fe: " << block_info->base_element[index] << std::endl;
    
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "Cell couplings table is represented by (line - equation, column - variable):" << std::endl;
    FcstUtilities::log << std::endl;
    
    for(unsigned int eq = 0; eq < equation_names.size(); ++eq)
    {
        for(unsigned int var = 0; var < solution_names.size(); ++var)
            FcstUtilities::log << (*cell_couplings)(eq,var) << "   ";
        
        FcstUtilities::log << std::endl;
    }
    
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "Flux couplings table is represented by (line - equation, column - variable):" << std::endl;
    FcstUtilities::log << std::endl;
    
    for(unsigned int eq = 0; eq < equation_names.size(); ++eq)
    {
        for(unsigned int var = 0; var < solution_names.size(); ++var)
            FcstUtilities::log << (*flux_couplings)(eq,var) << "   ";
        
        FcstUtilities::log << std::endl;
    }
    
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "FCST is using the mode: ";
    if(block_info->local_renumbering.size() == 0)
        FcstUtilities::log << "standard deal.II mode" << std::endl;
    else
        FcstUtilities::log << "block-wise renumbering mode" << std::endl;
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "------------------------------";
    FcstUtilities::log << std::endl;
}