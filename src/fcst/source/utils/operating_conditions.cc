//---------------------------------------------------------------------------
//    $Id: operating_conditions.cc 2605 2014-08-15 03:36:44Z secanell $
//
//    Copyright (C) 2009 by Marc Secanell
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include "utils/operating_conditions.h"

//---------------------------------------------------------------------------

FuelCell::OperatingConditions::OperatingConditions()
{
    R = Constants::R();
}

//---------------------------------------------------------------------------
FuelCell::OperatingConditions::~OperatingConditions()
{}

//---------------------------------------------------------------------------
void
FuelCell::OperatingConditions::declare_parameters (ParameterHandler &param) const
{
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Operating conditions");
        {
            param.declare_entry("Adjust initial solution and boundary conditions",
                                "false",
                                Patterns::Bool(),
                                "Use the parameters in Operating conditions to create an initial solution"
                                "and overwrite the boundary conditions for the problem specified in Equations>>Initial Data"
                                "using the parameters in this section");
            param.declare_entry ("Temperature cell [K]",
                                "353", // K, or 80 Celsius
                                Patterns::Double());
            param.declare_entry ("Cathode pressure [Pa]",
                                 "101325", // Pa, or 1 atm
                                 Patterns::Double(),
                                 "Pressure at the cathode channel in Pa");
            param.declare_entry ("Cathode initial oxygen mole fraction (prior to humidification)",
                                 "0.21", // 
                                 Patterns::Double());
            param.declare_entry ("Cathode relative humidity",
                                 "0.7",
                                 Patterns::Double(),
                                 "Relative humidity (fraction) in the cathode channel");
            param.declare_entry ("Anode pressure [Pa]",
                                 "101325", // Pa, or 1 atm
                                 Patterns::Double(),
                                 "Pressure at the anode channel in Pa");
            param.declare_entry ("Anode relative humidity",
                                 "0.7",
                                 Patterns::Double(),
                                 "Relative humidity (fraction) in the anode channel");
            param.declare_entry ("Anode dry hydrogen mole fraction",
                                 "0.9", // 
                                 Patterns::Double());
            param.declare_entry ("Voltage cell [V]",
                                 "0.6",
                                 Patterns::Double());
            param.declare_entry("Voltage drop in the anode [V]",
                                "0.015",
                                Patterns::Double());
            param.declare_entry ("Open circuit voltage [V]",
                                 "1.23",
                                 Patterns::Double());
        }
        param.leave_subsection();
    }
    param.leave_subsection();

    // Declare application specific parameters:
    param.enter_subsection("Application");
    {
        param.enter_subsection("Electrode KG");
        {
            param.declare_entry("Anode simulation",
                                "false",
                                Patterns::Bool(),
                                "Set to true if you want to simulate an anode. Otherwise, it is assumed that the user would like to run a cathode case.");
        }
        param.leave_subsection();
    }
    param.leave_subsection();    
}

//---------------------------------------------------------------------------
void
FuelCell::OperatingConditions::initialize (ParameterHandler& param)
{   
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Operating conditions");
        {
            adjust_BC = param.get_bool("Adjust initial solution and boundary conditions");
            T_cell = param.get_double("Temperature cell [K]");
            // Cathode
            p_c = param.get_double("Cathode pressure [Pa]");
            c_c = p_c/(R*T_cell)*1E-6; //mol/cm^3
            RH_c = param.get_double("Cathode relative humidity");
            channel_oxygen_mole_fraction = param.get_double("Cathode initial oxygen mole fraction (prior to humidification)");
            // Anode
            p_a = param.get_double("Anode pressure [Pa]");
            c_a = p_a/(R*T_cell)*1E-6; //mol/cm^3
            RH_a = param.get_double("Anode relative humidity");
            channel_dry_hydrogen_mole_fraction = param.get_double("Anode dry hydrogen mole fraction");
            // Cell voltage
            V_cell = param.get_double("Voltage cell [V]");
            dV_a = param.get_double("Voltage drop in the anode [V]");
            OCV = param.get_double("Open circuit voltage [V]");
            
            if (OCV > voltage_cell_th()) // NOTE: This also initializes E_th
                OCV = E_th;
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    param.enter_subsection("Application");
    {
        param.enter_subsection("Electrode KG");
        {
            anode = param.get_bool("Anode simulation");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    if ( anode )
        compartment = ANODE;
    else
        compartment = CATHODE;
}
//---------------------------------------------------------------------------
void
FuelCell::OperatingConditions::adjust_initial_solution(std::vector< component_materialID_value_map >& maps,
                                                       const boost::shared_ptr< FuelCellShop::Geometry::GridBase<dim> > grid) const
{
    for(unsigned int i = 0; i < maps.size(); ++i)
    {
        for(component_materialID_value_map::iterator iter  = maps[i].begin(); iter != maps[i].end(); ++iter)
        {      
            if ((iter->first.compare("oxygen_molar_fraction") == 0) && adjust_BC)
            {
                std::vector<unsigned int> material_ID = grid->get_material_id("Cathode GDL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_x_o2();
                
                material_ID = grid->get_material_id("Cathode MPL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_x_o2();
                
                material_ID = grid->get_material_id("Cathode CL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_x_o2();
                
            }
            else if ((iter->first.compare("hydrogen_molar_fraction") == 0) && adjust_BC)
            {
                std::string mesh_type = grid->get_mesh_type();
                std::string name = "Anode";
                if ( mesh_type.compare("CathodeMPL") == 0 )        // Sometimes, a cathode mesh is used for an anode. In this case, set mole fractions in anode.
                    name = "Cathode";
                    
                std::vector<unsigned int> material_ID = grid->get_material_id(name+" GDL");                
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_x_h2();
                
                material_ID = grid->get_material_id(name+" MPL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_x_h2();
                
                material_ID = grid->get_material_id(name+" CL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_x_h2();     
            }
            else if ((iter->first.compare("nitrogen_molar_fraction") == 0) && adjust_BC)
            {
                std::string mesh_type = grid->get_mesh_type();
                std::string name = "Anode";
                if ( mesh_type.compare("CathodeMPL") == 0 )        // Sometimes, a cathode mesh is used for an anode. In this case, set mole fractions in anode.
                    name = "Cathode";
                    
                std::vector<unsigned int> material_ID = grid->get_material_id(name+" GDL");                
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_x_n2(compartment);
                
                material_ID = grid->get_material_id(name+" MPL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_x_n2(compartment);
                
                material_ID = grid->get_material_id(name+" CL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_x_n2(compartment);     
            }
            else if ((iter->first.compare("water_molar_fraction") == 0) && adjust_BC)
            {
                std::vector<unsigned int> material_ID = grid->get_material_id("Cathode GDL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_x_wv(compartment);
                
                material_ID = grid->get_material_id("Cathode MPL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_x_wv(compartment);

                material_ID = grid->get_material_id("Cathode CL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)            
                    iter->second[material_ID[ind]] = this->get_x_wv(compartment);
                
                material_ID = grid->get_material_id("Anode CL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_x_wv(compartment);
                
                material_ID = grid->get_material_id("Anode MPL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_x_wv(compartment);
                
                material_ID = grid->get_material_id("Anode GDL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_x_wv(compartment);
            }
            //
            else if ((iter->first.compare("electronic_electrical_potential") == 0) && adjust_BC)
            {
                std::vector<unsigned int> material_ID = grid->get_material_id("Cathode gas channel");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                   iter->second[material_ID[ind]] = this->get_V();
                
                material_ID = grid->get_material_id("Cathode GDL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                   iter->second[material_ID[ind]] = this->get_V();      
                
                material_ID = grid->get_material_id("Cathode MPL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_V(); 
                
                material_ID = grid->get_material_id("Cathode CL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_V(); 

            }
            //
            else if ((iter->first.compare("protonic_electrical_potential") == 0) && adjust_BC)
            {
                std::vector<unsigned int> material_ID = grid->get_material_id("Cathode gas channel");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = (this->get_V() - OCV)/4.0;
                
                material_ID = grid->get_material_id("Membrane");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = (this->get_V() - OCV)/4.0;
                
                material_ID = grid->get_material_id("Cathode CL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = (this->get_V() - OCV)/2.0; 

                material_ID = grid->get_material_id("Anode CL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = -0.0001;

            }
            //
            else if ((iter->first.compare("membrane_water_content") == 0) && adjust_BC)
            {
                std::vector<unsigned int> material_ID = grid->get_material_id("Membrane");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = 5.0;
                
                material_ID = grid->get_material_id("Cathode CL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = 4.0; 
                
                material_ID = grid->get_material_id("Anode CL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = 7.0;
            }
            //
            else if ((iter->first.compare("temperature_of_REV") == 0) && adjust_BC)
            {
                std::vector<unsigned int> material_ID = grid->get_material_id("Cathode GDL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                   iter->second[material_ID[ind]] = this->get_T();

                material_ID = grid->get_material_id("Cathode MPL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_T();

                material_ID = grid->get_material_id("Cathode CL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)            
                    iter->second[material_ID[ind]] = this->get_T();
                
                material_ID = grid->get_material_id("Membrane");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_T();
                
                material_ID = grid->get_material_id("Anode CL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_T();
                
                material_ID = grid->get_material_id("Anode MPL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_T();
                
                material_ID = grid->get_material_id("Anode GDL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_T();

            }
            else if ((iter->first.compare(0,16,"density_species_") == 0) && adjust_BC)
            {
                //-- Get index:
                // Chad: std::string species = iter->first.substr(16); //get substring from position 16 to end which should be the full species number
                // Chad: double density = get_density_from_map(species); //Take species number and return the respective species density
                unsigned int species_index = get_species_index(iter->first);                
                
                //-- Check index within range of cathode and anode gas mix.
                unsigned int total_n_gases = cathode_mix->n_gases() + anode_mix->n_gases();
                AssertThrow(species_index < total_n_gases, ExcMessage("OperatingConditions::There are more species than gases in anode and cathode"));
                
                // I assume species index corresponds with cathode gas mix index or if larger, it belongs to anode:
                // Cathode:
                if (species_index < cathode_mix->n_gases())
                {
                    double density = get_cathode_gas_density(species_index);
                    
                    std::vector<unsigned int> material_ID = grid->get_material_id("Cathode gas channel");
                    for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                        iter->second[material_ID[ind]] = density;
                    material_ID = grid->get_material_id("Cathode GDL");
                    for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                        iter->second[material_ID[ind]] = density;
                    material_ID = grid->get_material_id("Cathode MPL");
                    for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                        iter->second[material_ID[ind]] = density;
                    material_ID = grid->get_material_id("Cathode CL");
                    for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                            iter->second[material_ID[ind]] = density;     
                }
                //-- Anode:   
                else
                {                 
                    double density = get_anode_gas_density(species_index);
                                
                    std::vector<unsigned int> material_ID = grid->get_material_id("Anode CL");
                    for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                        iter->second[material_ID[ind]] = density;
                    material_ID = grid->get_material_id("Anode MPL");
                    for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                        iter->second[material_ID[ind]] = density;
                    material_ID = grid->get_material_id("Anode GDL");
                    for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                        iter->second[material_ID[ind]] = density;
                    material_ID = grid->get_material_id("Anode gas channel");
                    for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                        iter->second[material_ID[ind]] = density;
                }
            }
        }
    }
}
//---------------------------------------------------------------------------
void
FuelCell::OperatingConditions::adjust_boundary_conditions(std::vector< component_boundaryID_value_map >& maps,
                                                          const boost::shared_ptr< FuelCellShop::Geometry::GridBase<dim> > grid) const
{
    for(unsigned int i = 0; i < maps.size(); ++i)
        for(component_boundaryID_value_map::iterator iter  = maps[i].begin(); iter != maps[i].end(); ++iter)
        {                       
            if ((iter->first.compare("oxygen_molar_fraction") == 0) && adjust_BC)
            {
                unsigned int boundary_ID = grid->get_boundary_id("c_Ch/GDL");
                iter->second[boundary_ID] = this->get_x_o2();
                boundary_ID = grid->get_boundary_id("a_GDL/Ch");
                iter->second[boundary_ID] = 0.0;
            }
            else if ((iter->first.compare("hydrogen_molar_fraction") == 0) && adjust_BC)
            {
                std::string mesh_type = grid->get_mesh_type();
                if ( mesh_type.compare("CathodeMPL") == 0 )        // Sometimes, a cathode mesh is used for an anode. In this case, set mole fractions in cathode as if it were an anode..
                {
                    unsigned int boundary_ID = grid->get_boundary_id("c_Ch/GDL");
                    iter->second[boundary_ID] = this->get_x_h2();
                }
                else // Otherwise, set it in the anode and set it to zero in the cathode:
                {                
                    unsigned int boundary_ID = grid->get_boundary_id("a_Ch/GDL");
                    iter->second[boundary_ID] = this->get_x_h2();
                    boundary_ID = grid->get_boundary_id("c_GDL/Ch");
                    iter->second[boundary_ID] = 0.0;
                }
            }
            else if ((iter->first.compare("water_molar_fraction") == 0) && adjust_BC)
            {
                unsigned int boundary_ID = grid->get_boundary_id("c_Ch/GDL");
                iter->second[boundary_ID] = this->get_x_wv(compartment);
                boundary_ID = grid->get_boundary_id("a_GDL/Ch");
                iter->second[boundary_ID] = this->get_x_wv(compartment);
            }
            else if ((iter->first.compare("electronic_electrical_potential") == 0) && adjust_BC)
            {
                unsigned int boundary_ID = grid->get_boundary_id("c_BPP/GDL");
                iter->second[boundary_ID] = this->get_V();
                boundary_ID = grid->get_boundary_id("c_Ch_Roof");
                iter->second[boundary_ID] = this->get_V();
                boundary_ID = grid->get_boundary_id("a_GDL/BPP");
                iter->second[boundary_ID] = 0.0;
            }
            else if ((iter->first.compare("temperature_of_REV") == 0) && adjust_BC)
            {
                unsigned int boundary_ID = grid->get_boundary_id("c_BPP/GDL");
                iter->second[boundary_ID] = this->get_T();  
                boundary_ID = grid->get_boundary_id("a_GDL/BPP");
                iter->second[boundary_ID] = this->get_T();
            }
            //Check substring and see if it is one of the possible density species
            else if ((iter->first.compare(0, 16, "density_species_") == 0) && adjust_BC)
            {
                //-- Get index:
                unsigned int species_index = get_species_index(iter->first);
                // Chad: std::string species = iter->first.substr(16); //get substring from position 16 to end which should be the full species number
                // Chad: double density = get_density_from_map(species); //Take species number and return the respective species density
                
                //-- Check index within range of cathode and anode gas mix.
                unsigned int total_n_gases = cathode_mix->n_gases() + anode_mix->n_gases();
                AssertThrow(species_index < total_n_gases, ExcMessage("OperatingConditions::There are more species than gases in anode and cathode"));
                
                // I assume species index corresponds with cathode gas mix index or if larger, it belongs to anode:
                // Cathode:
                if (species_index < cathode_mix->n_gases())
                {
                    double density = get_cathode_gas_density(species_index);
                    
                    unsigned int boundary_ID = grid->get_boundary_id("c_Ch_Inlet"); //If with channel simulation                    
                    iter->second[boundary_ID] = density;
                    boundary_ID = grid->get_boundary_id("c_Ch/GDL"); //If without channel simulation
                    iter->second[boundary_ID] = density;               
                
                }
                //-- Anode:   
                else
                {                 
                    double density = get_anode_gas_density(species_index);
                                
                    unsigned int boundary_ID = grid->get_boundary_id("a_Ch_Inlet"); //If with channel simulation                    
                    iter->second[boundary_ID] = density;
                    boundary_ID = grid->get_boundary_id("a_GDL/Ch"); //If without channel simulation
                    iter->second[boundary_ID] = density;  
                }
            }
        }
}

//---------------------------------------------------------------------------
double
FuelCell::OperatingConditions::saturation_pressure() const
{
    double T_celsius = T_cell - 273;
    return pow(10,(-2.1794+0.02953*T_celsius-0.000091837*pow(T_celsius,2)+0.00000014454*pow(T_celsius,3)));
}

//---------------------------------------------------------------------------
double 
FuelCell::OperatingConditions::get_x_wv(int compartment) const 
{
    double x_wv = 0.0;
    if (compartment == ANODE)
        x_wv = this->get_x_wv_anode();
    else if (compartment == CATHODE)
        x_wv = this->get_x_wv();
    else
    {
        FcstUtilities::log << "Feature not implemented yet." << std::endl;
        Assert(false, ExcNotImplemented());
    }
    
    return x_wv; 
}

//---------------------------------------------------------------------------
double 
FuelCell::OperatingConditions::get_x_wv() const 
{
    return RH_c*(saturation_pressure()*Units::convert(1.,Units::ATM_to_PA))/p_c; 
}

//---------------------------------------------------------------------------
double 
FuelCell::OperatingConditions::get_x_wv_anode() const 
{
    return RH_a*(saturation_pressure()*Units::convert(1.,Units::ATM_to_PA))/p_a; 
}

//---------------------------------------------------------------------------
double 
FuelCell::OperatingConditions::get_x_o2() const
{
    return channel_oxygen_mole_fraction*(p_c-(saturation_pressure()*Units::convert(1.,Units::ATM_to_PA))*RH_c)/p_c; 
}

//---------------------------------------------------------------------------
double 
FuelCell::OperatingConditions::get_x_h2() const
{
    return channel_dry_hydrogen_mole_fraction*(p_a-(saturation_pressure()*Units::convert(1.,Units::ATM_to_PA))*RH_a)/p_a; 
}

//---------------------------------------------------------------------------
double 
FuelCell::OperatingConditions::get_x_n2(int compartment) const
{
    double x_n2 = 0.0;
    if (compartment == ANODE)
        x_n2 = 1 - this->get_x_h2() - this->get_x_wv(compartment);
    else if (compartment == CATHODE)
        x_n2 = 1 - this->get_x_o2() - this->get_x_wv(compartment);
    else
    {
        FcstUtilities::log << "Feature not implemented yet." << std::endl;
        Assert(false, ExcNotImplemented());
    }
    
    return x_n2; 
}

//---------------------------------------------------------------------------
double 
FuelCell::OperatingConditions::get_rho_wv(int compartment) const 
{
    double rho_wv = 0.0;
    if (compartment == ANODE)
        rho_wv = this->get_rho_wv_anode();
    else if (compartment == CATHODE)
        rho_wv = this->get_rho_wv();
    else
    {
        FcstUtilities::log << "Feature not implemented yet." << std::endl;
        Assert(false, ExcNotImplemented());
    }        
    
    return rho_wv;
}

//---------------------------------------------------------------------------
double 
FuelCell::OperatingConditions::get_rho_wv() const 
{
    double totalConcentration = get_c_c();                      // [mol/cm^3]
    double wvMoleFraction     = get_x_wv(compartment);          // [mol]
    
    return totalConcentration * wvMoleFraction * M_water;       // [g/cm^3]
}

//---------------------------------------------------------------------------
double 
FuelCell::OperatingConditions::get_rho_wv_anode() const 
{
    double totalConcentration = get_c_a();                      // [mol/cm^3]
    double wvMoleFraction     = get_x_wv(compartment);          // [mol]
    
    return totalConcentration * wvMoleFraction * M_water;       // [g/cm^3]
}

//---------------------------------------------------------------------------
double 
FuelCell::OperatingConditions::get_rho_o2() const 
{
    double totalConcentration = get_c_c();                      // [mol/cm^3]
    double o2MoleFraction     = get_x_o2();                     // [mol]
    
    return totalConcentration * o2MoleFraction * M_oxygen;      // [g/cm^3]
}

//---------------------------------------------------------------------------
double 
FuelCell::OperatingConditions::get_rho_h2() const 
{
  double totalConcentration = get_c_a();                        // [mol/cm^3]
  double h2MoleFraction     = get_x_h2();                       // [mol]
  
  return totalConcentration * h2MoleFraction * M_hydrogen;      // [g/cm^3]
}

//---------------------------------------------------------------------------
double 
FuelCell::OperatingConditions::get_rho_n2(int compartment) const
{
    double rho = 0.0;
    if (compartment == ANODE)
        rho = this->get_rho_n2_anode();
    else if (compartment == CATHODE)
        rho = this->get_rho_n2();
    else
    {
        FcstUtilities::log << "Feature not implemented yet." << std::endl;
        Assert(false, ExcNotImplemented());
    }
    return rho;
}

//---------------------------------------------------------------------------
double 
FuelCell::OperatingConditions::get_rho_n2() const 
{
    double totalConcentration = get_c_c();                      // [mol/cm^3]
    double n2MoleFraction     = this->get_x_n2(compartment);     // [mol]
  
    return totalConcentration * n2MoleFraction * M_nitrogen;    // [g/cm^3]
}

//---------------------------------------------------------------------------
double 
FuelCell::OperatingConditions::get_rho_n2_anode() const 
{
    double totalConcentration = get_c_a();                      // [mol/cm^3]
    double n2MoleFraction     = this->get_x_n2(compartment);    // [mol]
 
    return totalConcentration * n2MoleFraction * M_nitrogen;    // [g/cm^3]
}

//---------------------------------------------------------------------------
void 
FuelCell::OperatingConditions::set_gas_map(std::map<unsigned int, std::string> tmp)
{
    gasSpeciesMap = tmp;
}

//---------------------------------------------------------------------------
std::map<unsigned int, std::string>
FuelCell::OperatingConditions::get_gas_map() const 
{
    return gasSpeciesMap;
}

//---------------------------------------------------------------------------
double
FuelCell::OperatingConditions::get_density_from_map(std::string speciesStr) const 
{    
    unsigned int speciesKey = atoi(speciesStr.c_str()); //convert species string number to int to use in map
    std::map<unsigned int, std::string> tmp = this->get_gas_map(); //get the gas species map
    double rho; //Initialize var to hold the density
    
    //Takes speciesKey and gas species map and determine what material density we need
    if(tmp[speciesKey] == "oxygen")
        rho = this->get_rho_o2();
    else if(tmp[speciesKey] == "water")
        rho = this->get_rho_wv(compartment);
    else if(tmp[speciesKey] == "nitrogen")
        rho = this->get_rho_n2(compartment);
    else
    {
        //TODO: add some break or output later
    }
    
    return rho; //return the density once found
}

//---------------------------------------------------------------------------
double 
FuelCell::OperatingConditions::voltage_cell_th()
{
    double p_h2 = (p_a/Units::convert(1.,Units::ATM_to_PA))*get_x_h2();
    double p_o2 = (p_c/Units::convert(1.,Units::ATM_to_PA))*get_x_o2();
    E_th  = 1.229 - (T_cell-298.15)*8.456E-4 + T_cell*4.31E-5*(log(p_h2)+0.5*log(p_o2));
    return E_th;
}


//---------------------------------------------------------------------------
void 
FuelCell::OperatingConditions::print_operating_conditions() const
{
    FcstUtilities::log<<"========= OPERATING CONDITIONS ========"<<std::endl;
    FcstUtilities::log<<"Temperature: "<< T_cell <<std::endl;
    FcstUtilities::log<<"Cathode Pressure: "<< p_c <<std::endl;
    FcstUtilities::log<<"Cathode RH: "<< RH_c <<std::endl;
    FcstUtilities::log<<"Anode Pressure: "<< p_a <<std::endl;
    FcstUtilities::log<<"Anode RH: "<< RH_a <<std::endl;
    /*
    FcstUtilities::log<<"Open Circuit Voltage: "<< OCV <<std::endl;
    FcstUtilities::log<<"Cell Voltage: "<< V_cell <<std::endl;
    FcstUtilities::log<<"O2 Density [g/cm^3]: "<< get_rho_o2() <<std::endl;
    FcstUtilities::log<<"H2O Density [g/cm^3]: "<< get_rho_wv() <<std::endl;
    FcstUtilities::log<<"N2 Density [g/cm^3]: "<< get_rho_n2() <<std::endl;
    */
    FcstUtilities::log<<"======================================="<<std::endl;
}