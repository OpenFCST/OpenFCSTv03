// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: PureGas.cc
// - Description: This class describes properties of pure gases
// - Developers: Valentin N. Zingan, University of Alberta
// - Id: $Id: PureGas.cc 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#include <materials/PureGas.h>

namespace NAME = FuelCellShop::Material;

//////////////////////////////////////////////////
//////////////////////////////////////////////////
// CONSTRUCTORS, DESTRUCTOR, AND INITIALIZATION //
//////////////////////////////////////////////////
//////////////////////////////////////////////////

// ---             ---
// --- Constructor ---
// ---             ---

NAME::PureGas::PureGas(const std::string& name)
:
NAME::BaseMaterial(name)
{ }

// ---            ---
// --- Destructor ---
// ---            ---

NAME::PureGas::~PureGas()
{ }

// ---                    ---
// --- declare_parameters ---
// ---                    ---

void
NAME::PureGas::declare_parameters(ParameterHandler& param) const
{
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Materials");
        {
            param.enter_subsection(this->name);
            {
                param.declare_entry("Dynamic viscosity mode",
                                    "Sutherland",
                                    Patterns::Selection("Sutherland | ChapmanEnskog"),
                                    " ");

                param.declare_entry("Bulk viscosity mode",
                                    "zero",
                                    Patterns::Selection("zero | Stokes"),
                                    " ");

                param.declare_entry("Thermal conductivity mode",
                                    "Sutherland",
                                    Patterns::Selection("Sutherland | ChapmanEnskog"),
                                    " ");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

// ---            ---
// --- initialize ---
// ---            ---

void
NAME::PureGas::initialize(ParameterHandler& param)
{
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Materials");
        {
            param.enter_subsection(this->name);
            {
                dynamic_viscosity_mode    = param.get("Dynamic viscosity mode");
                bulk_viscosity_mode       = param.get("Bulk viscosity mode");
                thermal_conductivity_mode = param.get("Thermal conductivity mode");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

///////////////////////
///////////////////////
// SERVICE FUNCTIONS //
///////////////////////
///////////////////////

/////////
// EoS //
/////////

// ---              ---
// --- get_pressure ---
// ---              ---

const double
NAME::PureGas::get_pressure(const double& density,
                            const double& temperature) const
{
    return density * ( Constants::R() / molar_mass ) * temperature;
}

// ---              ---
// --- get_pressure ---
// ---              ---

void
NAME::PureGas::get_pressure(const std::vector<double>& density,
                            const double&              temperature,
                            std::vector<double>&       pressure) const
{
    AssertThrow( density.size() == pressure.size() , ExcDimensionMismatch(density.size(), pressure.size()) );

    for(unsigned int q = 0; q < density.size(); ++q)
        pressure[q] = density[q] * ( Constants::R() / molar_mass ) * temperature;
}

// ---              ---
// --- get_pressure ---
// ---              ---

void
NAME::PureGas::get_pressure(const std::vector<double>& density,
                            const std::vector<double>& temperature,
                            std::vector<double>&       pressure) const
{
    AssertThrow( density.size()     == pressure.size() , ExcDimensionMismatch(density.size(),     pressure.size()) );
    AssertThrow( temperature.size() == pressure.size() , ExcDimensionMismatch(temperature.size(), pressure.size()) );

    for(unsigned int q = 0; q < density.size(); ++q)
        pressure[q] = density[q] * ( Constants::R() / molar_mass ) * temperature[q];
}

// ---                        ---
// --- get_Dpressure_Ddensity ---
// ---                        ---

const double
NAME::PureGas::get_Dpressure_Ddensity(const double& temperature) const
{
    return ( Constants::R() / molar_mass ) * temperature;
}

// ---                        ---
// --- get_Dpressure_Ddensity ---
// ---                        ---

void
NAME::PureGas::get_Dpressure_Ddensity(const std::vector<double>& temperature,
                                      std::vector<double>&       dst) const
{
    AssertThrow( temperature.size() == dst.size() , ExcDimensionMismatch(temperature.size(), dst.size()) );

    for(unsigned int q = 0; q < temperature.size(); ++q)
        dst[q] = ( Constants::R() / molar_mass ) * temperature[q];
}

// ---                            ---
// --- get_Dpressure_Dtemperature ---
// ---                            ---

const double
NAME::PureGas::get_Dpressure_Dtemperature(const double& density) const
{
    return density * ( Constants::R() / molar_mass );
}

// ---                            ---
// --- get_Dpressure_Dtemperature ---
// ---                            ---

void
NAME::PureGas::get_Dpressure_Dtemperature(const std::vector<double>& density,
                                          std::vector<double>&       dst) const
{
    AssertThrow( density.size() == dst.size() , ExcDimensionMismatch(density.size(), dst.size()) );

    for(unsigned int q = 0; q < density.size(); ++q)
        dst[q] = density[q] * ( Constants::R() / molar_mass );
}

// ---             ---
// --- get_density ---
// ---             ---

const double
NAME::PureGas::get_density(const double& temperature, // K
                           const double& pressure, //Pa
                           const double& molarMass // kg/mol
                          ) const
{
    double density = molarMass * pressure / (Constants::R() * temperature); // kg/m^3
    return density;
}

//////////////////////////////////
// Sutherland dynamic viscosity //
//////////////////////////////////

// ---                                  ---
// --- get_Sutherland_dynamic_viscosity ---
// ---                                  ---

const double
NAME::PureGas::get_Sutherland_dynamic_viscosity(const double& temperature) const
{
    return A_Sutherland * std::pow(temperature, 1.5) / (temperature + B_Sutherland);
}

// ---                                  ---
// --- get_Sutherland_dynamic_viscosity ---
// ---                                  ---

void
NAME::PureGas::get_Sutherland_dynamic_viscosity(const std::vector<double>& temperature,
                                                std::vector<double>&       dynamic_viscosity) const
{
    AssertThrow( temperature.size() == dynamic_viscosity.size() , ExcDimensionMismatch(temperature.size(), dynamic_viscosity.size()) );

    for(unsigned int q = 0; q < temperature.size(); ++q)
        dynamic_viscosity[q] = A_Sutherland * std::pow(temperature[q], 1.5) / (temperature[q] + B_Sutherland);
}

// ---                                                ---
// --- get_DSutherland_dynamic_viscosity_Dtemperature ---
// ---                                                ---

const double
NAME::PureGas::get_DSutherland_dynamic_viscosity_Dtemperature(const double& temperature) const
{
    return 1.5 * A_Sutherland * std::pow(temperature, 0.5) / (temperature + B_Sutherland)
            -
            A_Sutherland * std::pow(temperature, 1.5) / std::pow((temperature + B_Sutherland), 2.0);
}

// ---                                                ---
// --- get_DSutherland_dynamic_viscosity_Dtemperature ---
// ---                                                ---

void
NAME::PureGas::get_DSutherland_dynamic_viscosity_Dtemperature(const std::vector<double>& temperature,
                                                              std::vector<double>&       dst) const
{
    AssertThrow( temperature.size() == dst.size() , ExcDimensionMismatch(temperature.size(), dst.size()) );

    for(unsigned int q = 0; q < temperature.size(); ++q)
        dst[q] = 1.5 * A_Sutherland * std::pow(temperature[q], 0.5) / (temperature[q] + B_Sutherland)
                -
                A_Sutherland * std::pow(temperature[q], 1.5) / std::pow((temperature[q] + B_Sutherland), 2.0);
}

//////////////////////////////////////
// Chapman Enskog dynamic viscosity //
//////////////////////////////////////

// ---                                     ---
// --- get_ChapmanEnskog_dynamic_viscosity ---
// ---                                     ---

const double
NAME::PureGas::get_ChapmanEnskog_dynamic_viscosity(const double& temperature) const
{
    const double omega = get_collision_integral(temperature);

    return 2.6693e-6 * std::sqrt(molar_mass*1.0e3*temperature) / (collision_diameter*collision_diameter*omega);
}

// ---                                     ---
// --- get_ChapmanEnskog_dynamic_viscosity ---
// ---                                     ---

void
NAME::PureGas::get_ChapmanEnskog_dynamic_viscosity(const std::vector<double>& temperature,
                                                   std::vector<double>&       dynamic_viscosity) const
{
    AssertThrow( temperature.size() == dynamic_viscosity.size() , ExcDimensionMismatch(temperature.size(), dynamic_viscosity.size()) );

    std::vector<double> omega;
    omega.resize(temperature.size());

    get_collision_integral(temperature, omega);

    for(unsigned int q = 0; q < temperature.size(); ++q)
        dynamic_viscosity[q] = 2.6693e-6 * std::sqrt(molar_mass*1.0e3*temperature[q]) / (collision_diameter*collision_diameter*omega[q]);
}

// ---                                                   ---
// --- get_DChapmanEnskog_dynamic_viscosity_Dtemperature ---
// ---                                                   ---

const double
NAME::PureGas::get_DChapmanEnskog_dynamic_viscosity_Dtemperature(const double& temperature) const
{
    const double omega     = get_collision_integral(temperature);
    const double Domega_DT = get_Dcollision_integral_Dtemperature(temperature);

    return ( 2.6693e-6 * std::sqrt(molar_mass*1.0e3*temperature) / (collision_diameter*collision_diameter*omega) ) * ( 0.5/temperature - Domega_DT/omega );
}

// ---                                                   ---
// --- get_DChapmanEnskog_dynamic_viscosity_Dtemperature ---
// ---                                                   ---

void
NAME::PureGas::get_DChapmanEnskog_dynamic_viscosity_Dtemperature(const std::vector<double>& temperature,
                                                                 std::vector<double>&       dst) const
{
    AssertThrow( temperature.size() == dst.size() , ExcDimensionMismatch(temperature.size(), dst.size()) );

    std::vector<double> omega;
    omega.resize(temperature.size());

    get_collision_integral(temperature, omega);

    std::vector<double> Domega_DT;
    Domega_DT.resize(temperature.size());

    get_Dcollision_integral_Dtemperature(temperature, Domega_DT);

    for(unsigned int q = 0; q < temperature.size(); ++q)
        dst[q] = ( 2.6693e-6 * std::sqrt(molar_mass*1.0e3*temperature[q]) / (collision_diameter*collision_diameter*omega[q]) ) * ( 0.5/temperature[q] - Domega_DT[q]/omega[q] );
}

///////////////////////
// dynamic viscosity //
///////////////////////

// ---                       ---
// --- get_dynamic_viscosity ---
// ---                       ---

const double
NAME::PureGas::get_dynamic_viscosity(const double& temperature) const
{
    if( dynamic_viscosity_mode == "Sutherland" )
        return get_Sutherland_dynamic_viscosity(temperature);

    else if( dynamic_viscosity_mode == "ChapmanEnskog" )
        return get_ChapmanEnskog_dynamic_viscosity(temperature);

    else
    {
        FcstUtilities::log << "Dynamic viscosity mode you specified does not exist" << std::endl;
        Assert(false, ExcInternalError());
    }
}

// ---                       ---
// --- get_dynamic_viscosity ---
// ---                       ---

void
NAME::PureGas::get_dynamic_viscosity(const std::vector<double>& temperature,
                                     std::vector<double>&       dynamic_viscosity) const
{
    AssertThrow( temperature.size() == dynamic_viscosity.size() , ExcDimensionMismatch(temperature.size(), dynamic_viscosity.size()) );

    if( dynamic_viscosity_mode == "Sutherland" )
        get_Sutherland_dynamic_viscosity(temperature, dynamic_viscosity);

    else if( dynamic_viscosity_mode == "ChapmanEnskog" )
        get_ChapmanEnskog_dynamic_viscosity(temperature, dynamic_viscosity);

    else
    {
        FcstUtilities::log << "Dynamic viscosity mode you specified does not exist" << std::endl;
        Assert(false, ExcInternalError());
    }
}

// ---                                     ---
// --- get_Ddynamic_viscosity_Dtemperature ---
// ---                                     ---

const double
NAME::PureGas::get_Ddynamic_viscosity_Dtemperature(const double& temperature) const
{
    if( dynamic_viscosity_mode == "Sutherland" )
        return get_DSutherland_dynamic_viscosity_Dtemperature(temperature);

    else if( dynamic_viscosity_mode == "ChapmanEnskog" )
        return get_DChapmanEnskog_dynamic_viscosity_Dtemperature(temperature);

    else
    {
        FcstUtilities::log << "Dynamic viscosity mode you specified does not exist" << std::endl;
        Assert(false, ExcInternalError());
    }
}

// ---                                     ---
// --- get_Ddynamic_viscosity_Dtemperature ---
// ---                                     ---

void
NAME::PureGas::get_Ddynamic_viscosity_Dtemperature(const std::vector<double>& temperature,
                                                   std::vector<double>&       dst) const
{
    AssertThrow( temperature.size() == dst.size() , ExcDimensionMismatch(temperature.size(), dst.size()) );

    if( dynamic_viscosity_mode == "Sutherland" )
        get_DSutherland_dynamic_viscosity_Dtemperature(temperature, dst);

    else if( dynamic_viscosity_mode == "ChapmanEnskog" )
        get_DChapmanEnskog_dynamic_viscosity_Dtemperature(temperature, dst);

    else
    {
        FcstUtilities::log << "Dynamic viscosity mode you specified does not exist" << std::endl;
        Assert(false, ExcInternalError());
    }
}

////////////////////
// bulk viscosity //
////////////////////

// ---                    ---
// --- get_bulk_viscosity ---
// ---                    ---

const double
NAME::PureGas::get_bulk_viscosity(const double& dynamic_viscosity) const
{
    if( bulk_viscosity_mode == "zero" )
        return 0.0;

    else if( bulk_viscosity_mode == "Stokes" )
        return -(2.0/3.0)*dynamic_viscosity;

    else
    {
        FcstUtilities::log << "Bulk viscosity mode you specified does not exist" << std::endl;
        Assert(false, ExcInternalError());
    }
}

// ---                    ---
// --- get_bulk_viscosity ---
// ---                    ---

void
NAME::PureGas::get_bulk_viscosity(const std::vector<double>& dynamic_viscosity,
                                  std::vector<double>&       bulk_viscosity) const
{
    AssertThrow( dynamic_viscosity.size() == bulk_viscosity.size() , ExcDimensionMismatch(dynamic_viscosity.size(), bulk_viscosity.size()) );

    for(unsigned int q = 0; q < dynamic_viscosity.size(); ++q)
    {
        if( bulk_viscosity_mode == "zero" )
            bulk_viscosity[q] = 0.0;

        else if( bulk_viscosity_mode == "Stokes" )
            bulk_viscosity[q] = -(2.0/3.0)*dynamic_viscosity[q];

        else
        {
            FcstUtilities::log << "Bulk viscosity mode you specified does not exist" << std::endl;
            Assert(false, ExcInternalError());
        }
    }
}

// ---                                  ---
// --- get_Dbulk_viscosity_Dtemperature ---
// ---                                  ---

const double
NAME::PureGas::get_Dbulk_viscosity_Dtemperature(const double& src) const
{
    if( bulk_viscosity_mode == "zero" )
        return 0.0;

    else if( bulk_viscosity_mode == "Stokes" )
        return -(2.0/3.0)*src;

    else
    {
        FcstUtilities::log << "Bulk viscosity mode you specified does not exist" << std::endl;
        Assert(false, ExcInternalError());
    }
}

// ---                                  ---
// --- get_Dbulk_viscosity_Dtemperature ---
// ---                                  ---

void
NAME::PureGas::get_Dbulk_viscosity_Dtemperature(const std::vector<double>& src,
                                                std::vector<double>&       dst) const
{
    AssertThrow( src.size() == dst.size() , ExcDimensionMismatch(src.size(), dst.size()) );

    for(unsigned int q = 0; q < src.size(); ++q)
    {
        if( bulk_viscosity_mode == "zero" )
            dst[q] = 0.0;

        else if( bulk_viscosity_mode == "Stokes" )
            dst[q] = -(2.0/3.0)*src[q];

        else
        {
            FcstUtilities::log << "Bulk viscosity mode you specified does not exist" << std::endl;
            Assert(false, ExcInternalError());
        }
    }
}

/////////////////////////////////////
// Sutherland thermal conductivity //
/////////////////////////////////////

// ---                                     ---
// --- get_Sutherland_thermal_conductivity ---
// ---                                     ---

const double
NAME::PureGas::get_Sutherland_thermal_conductivity(const double& temperature) const
{
    return get_specific_heat_capacity(temperature) * get_Sutherland_dynamic_viscosity(temperature) / Prandtl;
}

// ---                                     ---
// --- get_Sutherland_thermal_conductivity ---
// ---                                     ---

void
NAME::PureGas::get_Sutherland_thermal_conductivity(const std::vector<double>& temperature,
                                                   std::vector<double>&       thermal_conductivity) const
{
    AssertThrow( temperature.size() == thermal_conductivity.size() , ExcDimensionMismatch(temperature.size(), thermal_conductivity.size()) );

    std::vector<double> specific_heat_capacity;
    specific_heat_capacity.resize(temperature.size());

    get_specific_heat_capacity(temperature, specific_heat_capacity);

    std::vector<double> dynamic_viscosity;
    dynamic_viscosity.resize(temperature.size());

    get_Sutherland_dynamic_viscosity(temperature, dynamic_viscosity);

    for(unsigned int q = 0; q < temperature.size(); ++q)
        thermal_conductivity[q] = specific_heat_capacity[q] * dynamic_viscosity[q] / Prandtl;
}

// ---                                                   ---
// --- get_DSutherland_thermal_conductivity_Dtemperature ---
// ---                                                   ---

const double
NAME::PureGas::get_DSutherland_thermal_conductivity_Dtemperature(const double& temperature) const
{
    return ( get_specific_heat_capacity(temperature)*get_DSutherland_dynamic_viscosity_Dtemperature(temperature)
             +
             get_Dspecific_heat_capacity_Dtemperature(temperature)*get_Sutherland_dynamic_viscosity(temperature) ) / Prandtl;
}

// ---                                                   ---
// --- get_DSutherland_thermal_conductivity_Dtemperature ---
// ---                                                   ---

void
NAME::PureGas::get_DSutherland_thermal_conductivity_Dtemperature(const std::vector<double>& temperature,
                                                                 std::vector<double>&       dst) const
{
    AssertThrow( temperature.size() == dst.size() , ExcDimensionMismatch(temperature.size(), dst.size()) );

    std::vector<double> specific_heat_capacity;
    specific_heat_capacity.resize(temperature.size());

    get_specific_heat_capacity(temperature, specific_heat_capacity);

    std::vector<double> dynamic_viscosity;
    dynamic_viscosity.resize(temperature.size());

    get_Sutherland_dynamic_viscosity(temperature, dynamic_viscosity);

    std::vector<double> Dspecific_heat_capacity_DT;
    Dspecific_heat_capacity_DT.resize(temperature.size());

    get_Dspecific_heat_capacity_Dtemperature(temperature, Dspecific_heat_capacity_DT);

    std::vector<double> Ddynamic_viscosity_DT;
    Ddynamic_viscosity_DT.resize(temperature.size());

    get_DSutherland_dynamic_viscosity_Dtemperature(temperature, Ddynamic_viscosity_DT);

    for(unsigned int q = 0; q < temperature.size(); ++q)
        dst[q] = ( specific_heat_capacity[q]*Ddynamic_viscosity_DT[q] + dynamic_viscosity[q]*Dspecific_heat_capacity_DT[q] ) / Prandtl;
}

/////////////////////////////////////////
// Chapman Enskog thermal conductivity //
/////////////////////////////////////////

// ---                                        ---
// --- get_ChapmanEnskog_thermal_conductivity ---
// ---                                        ---

const double
NAME::PureGas::get_ChapmanEnskog_thermal_conductivity(const double& temperature) const
{
    const double omega = get_collision_integral(temperature);

    return 8.3224e-3 * std::sqrt(temperature/(molar_mass*1.0e3)) / (collision_diameter*collision_diameter*omega);
}

// ---                                        ---
// --- get_ChapmanEnskog_thermal_conductivity ---
// ---                                        ---

void
NAME::PureGas::get_ChapmanEnskog_thermal_conductivity(const std::vector<double>& temperature,
                                                      std::vector<double>&       thermal_conductivity) const
{
    AssertThrow( temperature.size() == thermal_conductivity.size() , ExcDimensionMismatch(temperature.size(), thermal_conductivity.size()) );

    std::vector<double> omega;
    omega.resize(temperature.size());

    get_collision_integral(temperature, omega);

    for(unsigned int q = 0; q < temperature.size(); ++q)
        thermal_conductivity[q] = 8.3224e-3 * std::sqrt(temperature[q]/(molar_mass*1.0e3)) / (collision_diameter*collision_diameter*omega[q]);
}

// ---                                                      ---
// --- get_DChapmanEnskog_thermal_conductivity_Dtemperature ---
// ---                                                      ---

const double
NAME::PureGas::get_DChapmanEnskog_thermal_conductivity_Dtemperature(const double& temperature) const
{
    const double omega     = get_collision_integral(temperature);
    const double Domega_DT = get_Dcollision_integral_Dtemperature(temperature);

    return ( 8.3224e-3 * std::sqrt(temperature/(molar_mass*1.0e3)) / (collision_diameter*collision_diameter*omega) ) * ( 0.5/temperature - Domega_DT/omega );
}

// ---                                                      ---
// --- get_DChapmanEnskog_thermal_conductivity_Dtemperature ---
// ---                                                      ---

void
NAME::PureGas::get_DChapmanEnskog_thermal_conductivity_Dtemperature(const std::vector<double>& temperature,
                                                                    std::vector<double>&       dst) const
{
    AssertThrow( temperature.size() == dst.size() , ExcDimensionMismatch(temperature.size(), dst.size()) );

    std::vector<double> omega;
    omega.resize(temperature.size());

    get_collision_integral(temperature, omega);

    std::vector<double> Domega_DT;
    Domega_DT.resize(temperature.size());

    get_Dcollision_integral_Dtemperature(temperature, Domega_DT);

    for(unsigned int q = 0; q < temperature.size(); ++q)
        dst[q] = ( 8.3224e-3 * std::sqrt(temperature[q]/(molar_mass*1.0e3)) / (collision_diameter*collision_diameter*omega[q]) ) * ( 0.5/temperature[q] - Domega_DT[q]/omega[q] );
}

//////////////////////////
// thermal conductivity //
//////////////////////////

// ---                          ---
// --- get_thermal_conductivity ---
// ---                          ---

const double
NAME::PureGas::get_thermal_conductivity(const double& temperature) const
{
    if( thermal_conductivity_mode == "Sutherland" )
        return get_Sutherland_thermal_conductivity(temperature);

    else if( thermal_conductivity_mode == "ChapmanEnskog" )
        return get_ChapmanEnskog_thermal_conductivity(temperature);

    else
    {
        FcstUtilities::log << "Thermal conductivity mode you specified does not exist" << std::endl;
        Assert(false, ExcInternalError());
    }
}

// ---                          ---
// --- get_thermal_conductivity ---
// ---                          ---

void
NAME::PureGas::get_thermal_conductivity(const std::vector<double>& temperature,
                                        std::vector<double>&       thermal_conductivity) const
{
    AssertThrow( temperature.size() == thermal_conductivity.size() , ExcDimensionMismatch(temperature.size(), thermal_conductivity.size()) );

    if( thermal_conductivity_mode == "Sutherland" )
        get_Sutherland_thermal_conductivity(temperature, thermal_conductivity);

    else if( thermal_conductivity_mode == "ChapmanEnskog" )
        get_ChapmanEnskog_thermal_conductivity(temperature, thermal_conductivity);

    else
    {
        FcstUtilities::log << "Thermal conductivity mode you specified does not exist" << std::endl;
        Assert(false, ExcInternalError());
    }
}

// ---                                        ---
// --- get_Dthermal_conductivity_Dtemperature ---
// ---                                        ---

const double
NAME::PureGas::get_Dthermal_conductivity_Dtemperature(const double& temperature) const
{
    if( thermal_conductivity_mode == "Sutherland" )
        return get_DSutherland_thermal_conductivity_Dtemperature(temperature);

    else if( thermal_conductivity_mode == "ChapmanEnskog" )
        return get_DChapmanEnskog_thermal_conductivity_Dtemperature(temperature);

    else
    {
        FcstUtilities::log << "Thermal conductivity mode you specified does not exist" << std::endl;
        Assert(false, ExcInternalError());
    }
}

// ---                                        ---
// --- get_Dthermal_conductivity_Dtemperature ---
// ---                                        ---

void
NAME::PureGas::get_Dthermal_conductivity_Dtemperature(const std::vector<double>& temperature,
                                                      std::vector<double>&       dst) const
{
    AssertThrow( temperature.size() == dst.size() , ExcDimensionMismatch(temperature.size(), dst.size()) );

    if( thermal_conductivity_mode == "Sutherland" )
        get_DSutherland_thermal_conductivity_Dtemperature(temperature, dst);

    else if( thermal_conductivity_mode == "ChapmanEnskog" )
        get_DChapmanEnskog_thermal_conductivity_Dtemperature(temperature, dst);

    else
    {
        FcstUtilities::log << "Thermal conductivity mode you specified does not exist" << std::endl;
        Assert(false, ExcInternalError());
    }
}

////////////////////
// molar enthalpy //
////////////////////

// ---                    ---
// --- get_molar_enthalpy ---
// ---                    ---

const double
NAME::PureGas::get_molar_enthalpy(const double& temperature) const
{
    return H_ref + get_specific_heat_capacity(temperature)*molar_mass*( temperature - T_ref );
}

// ---                    ---
// --- get_molar_enthalpy ---
// ---                    ---

void
NAME::PureGas::get_molar_enthalpy(const std::vector<double>& temperature,
                                  std::vector<double>&       molar_enthalpy) const
{
    AssertThrow( temperature.size() == molar_enthalpy.size() , ExcDimensionMismatch(temperature.size(), molar_enthalpy.size()) );

    std::vector<double> specific_heat_capacity;
    specific_heat_capacity.resize(temperature.size());

    get_specific_heat_capacity(temperature, specific_heat_capacity);

    for(unsigned int q = 0; q < temperature.size(); ++q)
        molar_enthalpy[q] = H_ref + specific_heat_capacity[q]*molar_mass*( temperature[q] - T_ref );
}

// ---                                  ---
// --- get_Dmolar_enthalpy_Dtemperature ---
// ---                                  ---

const double
NAME::PureGas::get_Dmolar_enthalpy_Dtemperature(const double& temperature) const
{
    return molar_mass*( get_specific_heat_capacity(temperature) + get_Dspecific_heat_capacity_Dtemperature(temperature)*(temperature - T_ref) );
}

// ---                                  ---
// --- get_Dmolar_enthalpy_Dtemperature ---
// ---                                  ---

void
NAME::PureGas::get_Dmolar_enthalpy_Dtemperature(const std::vector<double>& temperature,
                                                std::vector<double>&       dst) const
{
    AssertThrow( temperature.size() == dst.size() , ExcDimensionMismatch(temperature.size(), dst.size()) );

    std::vector<double> specific_heat_capacity;
    specific_heat_capacity.resize(temperature.size());

    get_specific_heat_capacity(temperature, specific_heat_capacity);

    std::vector<double> Dspecific_heat_capacity_DT;
    Dspecific_heat_capacity_DT.resize(temperature.size());

    get_Dspecific_heat_capacity_Dtemperature(temperature, Dspecific_heat_capacity_DT);

    for(unsigned int q = 0; q < temperature.size(); ++q)
        dst[q] = molar_mass*( specific_heat_capacity[q] + Dspecific_heat_capacity_DT[q]*(temperature[q] - T_ref) );
}

// ---                                    ---
// --- get_D2molar_enthalpy_Dtemperature2 ---
// ---                                    ---

const double
NAME::PureGas::get_D2molar_enthalpy_Dtemperature2(const double& temperature) const
{
    return molar_mass*( 2.0*get_Dspecific_heat_capacity_Dtemperature(temperature) + get_D2specific_heat_capacity_Dtemperature2(temperature)*(temperature - T_ref) );
}

// ---                                    ---
// --- get_D2molar_enthalpy_Dtemperature2 ---
// ---                                    ---

void
NAME::PureGas::get_D2molar_enthalpy_Dtemperature2(const std::vector<double>& temperature,
                                                  std::vector<double>&       dst) const
{
    AssertThrow( temperature.size() == dst.size() , ExcDimensionMismatch(temperature.size(), dst.size()) );

    std::vector<double> Dspecific_heat_capacity_DT;
    Dspecific_heat_capacity_DT.resize(temperature.size());

    get_Dspecific_heat_capacity_Dtemperature(temperature, Dspecific_heat_capacity_DT);

    std::vector<double> D2specific_heat_capacity_DT2;
    D2specific_heat_capacity_DT2.resize(temperature.size());

    get_D2specific_heat_capacity_Dtemperature2(temperature, D2specific_heat_capacity_DT2);

    for(unsigned int q = 0; q < temperature.size(); ++q)
        dst[q] = molar_mass*( 2.0*Dspecific_heat_capacity_DT[q] + D2specific_heat_capacity_DT2[q]*(temperature[q] - T_ref) );
}

/////////////////////////////////////
// water vapor saturation pressure //
/////////////////////////////////////

// ---                                     ---
// --- get_water_vapor_saturation_pressure ---
// ---                                     ---

const double
NAME::PureGas::get_water_vapor_saturation_pressure(const double& temperature) const
{
    AssertThrow( get_chemical_formula() == "H2O Vapor",
                ExcMessage("Saturation pressure is computed only for water vapor") );
    
    return 101325.0 * std::pow( 10.0,
                                Constants::b_0()
                                +
                                Constants::b_1()*(temperature-273.15)
                                +
                                Constants::b_2()*(temperature-273.15)*(temperature-273.15)
                                +
                                Constants::b_3()*(temperature-273.15)*(temperature-273.15)*(temperature-273.15) );
}

// ---                                     ---
// --- get_water_vapor_saturation_pressure ---
// ---                                     ---

void
NAME::PureGas::get_water_vapor_saturation_pressure(const std::vector<double>& temperature,
                                                   std::vector<double>&       water_vapor_saturation_pressure) const
{
    AssertThrow( get_chemical_formula() == "H2O Vapor",
                ExcMessage("Saturation pressure is computed only for water vapor") );

    AssertThrow( temperature.size() == water_vapor_saturation_pressure.size() , ExcDimensionMismatch(temperature.size(), water_vapor_saturation_pressure.size()) );

    for(unsigned int q = 0; q < temperature.size(); ++q)
        water_vapor_saturation_pressure[q] = 101325.0 * std::pow( 10.0,
                                                                  Constants::b_0()
                                                                  +
                                                                  Constants::b_1()*(temperature[q]-273.15)
                                                                  +
                                                                  Constants::b_2()*(temperature[q]-273.15)*(temperature[q]-273.15)
                                                                  +
                                                                  Constants::b_3()*(temperature[q]-273.15)*(temperature[q]-273.15)*(temperature[q]-273.15) );
}

// ---                                                   ---
// --- get_Dwater_vapor_saturation_pressure_Dtemperature ---
// ---                                                   ---

const double
NAME::PureGas::get_Dwater_vapor_saturation_pressure_Dtemperature(const double& temperature) const
{
    AssertThrow( get_chemical_formula() == "H2O Vapor",
                ExcMessage("Derivative of saturation pressure is computed only for water vapor") );

    const double water_vapor_saturation_pressure = get_water_vapor_saturation_pressure(temperature);

    return water_vapor_saturation_pressure * std::log(10.0) * ( 1.0*Constants::b_1()
                                                                +
                                                                2.0*Constants::b_2()*(temperature-273.15)
                                                                +
                                                                3.0*Constants::b_3()*(temperature-273.15)*(temperature-273.15) );
}

// ---                                                   ---
// --- get_Dwater_vapor_saturation_pressure_Dtemperature ---
// ---                                                   ---

void
NAME::PureGas::get_Dwater_vapor_saturation_pressure_Dtemperature(const std::vector<double>& temperature,
                                                                 std::vector<double>&       dst) const
{
    AssertThrow( get_chemical_formula() == "H2O Vapor",
                ExcMessage("Derivative of saturation pressure is computed only for water vapor") );

    AssertThrow( temperature.size() == dst.size() , ExcDimensionMismatch(temperature.size(), dst.size()) );

    std::vector<double> water_vapor_saturation_pressure;
    water_vapor_saturation_pressure.resize(temperature.size());

    get_water_vapor_saturation_pressure(temperature, water_vapor_saturation_pressure);

    for(unsigned int q = 0; q < temperature.size(); ++q)
        dst[q] = water_vapor_saturation_pressure[q] * std::log(10.0) * ( 1.0*Constants::b_1()
                                                                         +
                                                                         2.0*Constants::b_2()*(temperature[q]-273.15)
                                                                         +
                                                                         3.0*Constants::b_3()*(temperature[q]-273.15)*(temperature[q]-273.15) );
}

////////////////////////
// collision integral //
////////////////////////

// ---                        ---
// --- get_collision_integral ---
// ---                        ---

const double
NAME::PureGas::get_collision_integral(const double& temperature) const
{
    const double tmp = temperature/eps_BY_k;

    return Constants::A_vk() / std::pow(tmp, Constants::B_vk())
           +
           Constants::C_vk() / std::exp(Constants::D_vk()*tmp)
           +
           Constants::E_vk() / std::exp(Constants::F_vk()*tmp);
}

// ---                        ---
// --- get_collision_integral ---
// ---                        ---

void
NAME::PureGas::get_collision_integral(const std::vector<double>& temperature,
                                      std::vector<double>&       collision_integral) const
{
    AssertThrow( temperature.size() == collision_integral.size() , ExcDimensionMismatch(temperature.size(), collision_integral.size()) );

    for(unsigned int q = 0; q < temperature.size(); ++q)
    {
        const double tmp = temperature[q]/eps_BY_k;

        collision_integral[q] = Constants::A_vk() / std::pow(tmp, Constants::B_vk())
                                +
                                Constants::C_vk() / std::exp(Constants::D_vk()*tmp)
                                +
                                Constants::E_vk() / std::exp(Constants::F_vk()*tmp);
    }
}

// ---                                      ---
// --- get_Dcollision_integral_Dtemperature ---
// ---                                      ---

const double
NAME::PureGas::get_Dcollision_integral_Dtemperature(const double& temperature) const
{
    const double tmp = temperature/eps_BY_k;

    return - (1.0/eps_BY_k)*( Constants::A_vk()*Constants::B_vk()*std::pow( tmp, (-Constants::B_vk()-1.0) )
                              +
                              Constants::C_vk()*Constants::D_vk()*std::exp( -Constants::D_vk()*tmp )
                              +
                              Constants::E_vk()*Constants::F_vk()*std::exp( -Constants::F_vk()*tmp )   );
}

// ---                                      ---
// --- get_Dcollision_integral_Dtemperature ---
// ---                                      ---

void
NAME::PureGas::get_Dcollision_integral_Dtemperature(const std::vector<double>& temperature,
                                                    std::vector<double>&       dst) const
{
    AssertThrow( temperature.size() == dst.size() , ExcDimensionMismatch(temperature.size(), dst.size()) );

    for(unsigned int q = 0; q < temperature.size(); ++q)
    {
        const double tmp = temperature[q]/eps_BY_k;

        dst[q] = - (1.0/eps_BY_k)*( Constants::A_vk()*Constants::B_vk()*std::pow( tmp, (-Constants::B_vk()-1.0) )
                                    +
                                    Constants::C_vk()*Constants::D_vk()*std::exp( -Constants::D_vk()*tmp )
                                    +
                                    Constants::E_vk()*Constants::F_vk()*std::exp( -Constants::F_vk()*tmp )   );
    }
}

/////////////////////////////////////////////////
// specific heat capacity at constant pressure //
/////////////////////////////////////////////////

// ---                            ---
// --- get_specific_heat_capacity ---
// ---                            ---

const double
NAME::PureGas::get_specific_heat_capacity(const double& temperature) const
{
    return c_0
           +
           c_1 * (temperature/1000.0)
           +
           c_2 * (temperature/1000.0) * (temperature/1000.0)
           +
           c_3 * (temperature/1000.0) * (temperature/1000.0) * (temperature/1000.0);
}

// ---                            ---
// --- get_specific_heat_capacity ---
// ---                            ---

void
NAME::PureGas::get_specific_heat_capacity(const std::vector<double>& temperature,
                                          std::vector<double>&       specific_heat_capacity) const
{
    AssertThrow( temperature.size() == specific_heat_capacity.size() , ExcDimensionMismatch(temperature.size(), specific_heat_capacity.size()) );

    for(unsigned int q = 0; q < temperature.size(); ++q)
        specific_heat_capacity[q] = c_0
                                    +
                                    c_1 * (temperature[q]/1000.0)
                                    +
                                    c_2 * (temperature[q]/1000.0) * (temperature[q]/1000.0)
                                    +
                                    c_3 * (temperature[q]/1000.0) * (temperature[q]/1000.0) * (temperature[q]/1000.0);
}

// ---                                          ---
// --- get_Dspecific_heat_capacity_Dtemperature ---
// ---                                          ---

const double
NAME::PureGas::get_Dspecific_heat_capacity_Dtemperature(const double& temperature) const
{
    return ( 1.0*c_1
             +
             2.0*c_2 * (temperature/1000.0)
             +
             3.0*c_3 * (temperature/1000.0) * (temperature/1000.0) ) / 1000.0;
}

// ---                                          ---
// --- get_Dspecific_heat_capacity_Dtemperature ---
// ---                                          ---

void
NAME::PureGas::get_Dspecific_heat_capacity_Dtemperature(const std::vector<double>& temperature,
                                                        std::vector<double>&       dst) const
{
    AssertThrow( temperature.size() == dst.size() , ExcDimensionMismatch(temperature.size(), dst.size()) );

    for(unsigned int q = 0; q < temperature.size(); ++q)
        dst[q] = ( 1.0*c_1
                   +
                   2.0*c_2 * (temperature[q]/1000.0)
                   +
                   3.0*c_3 * (temperature[q]/1000.0) * (temperature[q]/1000.0) ) / 1000.0;
}

// ---                                            ---
// --- get_D2specific_heat_capacity_Dtemperature2 ---
// ---                                            ---

const double
NAME::PureGas::get_D2specific_heat_capacity_Dtemperature2(const double& temperature) const
{
    return ( 2.0*c_2 + 6.0*c_3 * (temperature/1000.0) ) / (1000.0*1000.0);
}

// ---                                            ---
// --- get_D2specific_heat_capacity_Dtemperature2 ---
// ---                                            ---

void
NAME::PureGas::get_D2specific_heat_capacity_Dtemperature2(const std::vector<double>& temperature,
                                                          std::vector<double>&       dst) const
{
    AssertThrow( temperature.size() == dst.size() , ExcDimensionMismatch(temperature.size(), dst.size()) );

    for(unsigned int q = 0; q < temperature.size(); ++q)
        dst[q] = ( 2.0*c_2 + 6.0*c_3 * (temperature[q]/1000.0) ) / (1000.0*1000.0);
}