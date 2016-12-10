// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: PureGas.h
// - Description: This class describes properties of pure gases
// - Developers: Marc Secanell and Valentin N. Zingan, University of Alberta
//
// ----------------------------------------------------------------------------

#ifndef _FCST_FUELCELLSHOP_MATERIAL_PUREGAS_H_
#define _FCST_FUELCELLSHOP_MATERIAL_PUREGAS_H_

#define _DUMMY_ 1.e300

#include <materials/base_material.h>
#include <utils/fcst_constants.h>

enum enMaterialID {matNone,
    matOxygen,
    matNitrogen,
    matHydrogen,
    matWaterVapor,
    matAir,
    matHelium,
    matArgon,
    matNeon,
    matAcetone,
    matMethanol};
    
namespace FuelCellShop
{
    namespace Material
    {
        /**
         * This class is a base class for all pure gases used in OpenFCST.
         *
         * This class contains data provided by its children
         * and implements the methods by means of which
         * different properties of pure gases can be computed.
         *
         * These methods include:
         *
         * - Pressure from the ideal gases,     \f$ \quad \left[ \text{Pa} \right]                   \f$   \f$ \quad p   = \rho \frac{R}{M} T \quad \f$ and partial derivatives \f$ \quad \frac{\partial p}{\partial \rho} \quad \f$
         *                                                                                                   and \f$ \quad \frac{\partial p}{\partial T} \f$
         * - Sutherland dynamic viscosity,        \f$ \quad \left[ \text{Pa sec} \right]               \f$   \f$ \quad \mu = \frac{A T^{3/2}}{T+B} \quad \f$ and partial derivatives \f$ \quad \frac{\partial \mu}{\partial T} \f$
         * - Chapman Enskog dynamic viscosity,    \f$ \quad \left[ \text{Pa sec} \right]               \f$   \f$ \quad \mu = 2.6693 \cdot 10^{-6} \frac{\sqrt{M T}}{\sigma^2 \Omega_{\mu}} \quad \f$ and partial derivatives
         *                                                                                                   \f$ \quad \frac{\partial \mu}{\partial T} \quad \f$
         *                                                                                                   where the collision integral \f$ \quad \Omega_{\mu} \quad \f$ is given by
         *                                                                                                   \f$ \quad \Omega_{\mu} = \frac{A_{vk}}{\left(\frac{kT}{\epsilon}\right)^{B_{vk}}} +
         *                                                                                                   \frac{C_{vk}}{e^{D_{vk}\frac{kT}{\epsilon}}} +
         *                                                                                                   \frac{E_{vk}}{e^{F_{vk}\frac{kT}{\epsilon}}} \f$
         * - Zero bulk viscosity,                 \f$ \quad \left[ \text{Pa sec} \right]               \f$   \f$ \quad \lambda = 0 \quad \f$ and partial derivatives \f$ \quad \frac{\partial \lambda}{\partial T} \f$
         * - Stokes bulk viscosity,               \f$ \quad \left[ \text{Pa sec} \right]               \f$   \f$ \quad \lambda = -\frac{2}{3} \mu \quad \f$ and partial derivatives \f$ \quad \frac{\partial \lambda}{\partial T} \f$
         * - Sutherland thermal conductivity,     \f$ \quad \left[ \frac{\text{W}}{\text{m K}} \right] \f$   \f$ \quad \kappa  = \frac{\mu c_p}{\text{Pr}} \quad \f$ and partial derivatives \f$ \quad \frac{\partial \kappa}
         *                                                                                                   {\partial T} \quad \f$
         *                                                                                                   where the specific heat capacity at constant pressure \f$ \quad c_p, \quad \left[ \frac{\text{J}}{\text{kg K}} \right] \quad \f$ is given by
         *                                                                                                   \f$ \quad c_p = c_0 + c_1 \left(\frac{T}{1000}\right)^1 + c_2 \left(\frac{T}{1000}\right)^2 +
         *                                                                                                   c_3 \left(\frac{T}{1000}\right)^3 \quad \f$ and \f$ \quad \mu \quad \f$ is given by the
         *                                                                                                   Sutherland dynamic viscosity formula
         * - Chapman Enskog thermal conductivity, \f$ \quad \left[ \frac{\text{W}}{\text{m K}} \right] \f$   \f$ \quad \kappa  = 8.3224 \cdot 10^{-3} \frac{\sqrt{T/M}}{\sigma^2 \Omega_k} \quad \f$ and partial derivatives
         *                                                                                                   \f$ \quad \frac{\partial \kappa}{\partial T} \quad \f$
         *                                                                                                   where the collision integral \f$ \quad \Omega_k \quad \f$ is given by
         *                                                                                                   \f$ \quad \Omega_k = \Omega_{\mu} \f$
         * - Molar enthalpy,                      \f$ \quad \left[ \frac{\text{J}}{\text{mol}} \right] \f$   \f$ \quad H = H^{\text{ref}} + c_p M \left( T - T^{\text{ref}} \right) \quad \f$
         *                                                                                                   and partial derivatives \f$ \quad \frac{\partial H}{\partial T} \quad \f$ and \f$ \quad \frac{\partial^2
         *                                                                                                   H}{\partial T^2} \quad \f$
         *                                                                                                   where \f$ \quad c_p \quad \f$ is given by the specific heat capacity at constant pressure formula
         * - Water vapor saturation pressure,     \f$ \quad \left[ \text{Pa} \right]                   \f$   \f$ \quad p_{\text{sat}}^{H_2O} = 10^{b_0 + b_1 \left( T-273.15 \right) + b_2 \left( T-273.15 \right)^2 +
         *                                                                                                   b_3 \left( T-273.15 \right)^3} \quad \f$ and partial derivatives
         *                                                                                                   \f$ \quad \frac{\partial p_{\text{sat}}^{H_2O}}{\partial T} \f$
         *
         * All those methods receive data in SI units.
         * All those methods return the results in SI units.
         *
         * \warning Do not create an object of this class. Use its children instead.
         *
         * \note For developers: Please update this info appropriately if you add a new method.
         *
         * \author Valentin N. Zingan, Chad Balen and Marc Secanell, 2013-16
         */
        
        class PureGas : public BaseMaterial
        {
            public:
                
                ///@name Constructors, destructor, and initialization
                //@{
                /**
                * Destructor.
                */
                virtual ~PureGas();
                /**
                * Declare parameters.
                */
                virtual void declare_parameters(ParameterHandler& param) const;
                
                /**
                * Initialize parameters.
                */
                virtual void initialize(ParameterHandler& param);
                
                //@}
                
                ///@name Accessors and info
                //@{
                
                /**
                * This function returns \p molar_mass [kg/mol].
                */
                const double& get_molar_mass() const
                {
                    return molar_mass;
                }
                
                /**
                * This function returns \p collision_diameter [Angstrom].
                */
                const double& get_collision_diameter() const
                {
                    return collision_diameter;
                }
                
                /**
                * This function returns \p eps_BY_k (K).
                */
                const double& get_eps_BY_k() const
                {
                    return eps_BY_k;
                }
                
                /**
                * This function returns \p Prandtl.
                */
                const double& get_Prandtl() const
                {
                    return Prandtl;
                }
                
                /**
                * This function returns a number from \p enMaterialID enumeration.
                *
                * Such that:
                * - Oxygen      = 1
                * - Nitrogen    = 2
                * - Hydrogen    = 3
                * - Water Vapor = 4
                * - Air         = 5
                * - Helium      = 6
                * - Argon       = 7
                * - Neon        = 8
                * - Acetone     = 9
                * - Methanol    = 10
                */
                virtual enMaterialID get_ID() const = 0;
                
                /**
                * This function returns \p chemical_formula.
                */
                const std::string& get_chemical_formula() const
                {
                    return chemical_formula;
                }
                
                /**
                * This function returns
                * \p dynamic_viscosity_mode.
                */
                const std::string& get_dynamic_viscosity_mode() const
                {
                    return dynamic_viscosity_mode;
                }
                
                /**
                * This function returns
                * \p bulk_viscosity_mode.
                */
                const std::string& get_bulk_viscosity_mode() const
                {
                    return bulk_viscosity_mode;
                }
                
                /**
                * This function returns
                * \p thermal_conductivity_mode.
                */
                const std::string& get_thermal_conductivity_mode() const
                {
                    return thermal_conductivity_mode;
                }
                
                //@}
                
                ///@name Service functions. EoS.
                //@{
                
                /**
                * This function returns pressure (Pa) of a pure ideal gas.
                *
                * @param density     - density [kg/m^3],
                * @param temperature - temperature [K].
                */
                const double get_pressure(const double& density, const double& temperature) const;
                
                /**
                * This function returns pressure (Pa) of a pure ideal gas in the quadrature points of a mesh entity
                * at a constant temperature (isothermal case).
                *
                * @param density     - density in the quadrature points of a mesh entity [kg/m^3],
                * @param temperature - temperature [K],
                * @param pressure    - pressure in the quadrature points of a mesh entity [Pa].
                */
                void get_pressure(const std::vector<double>& density, const double& temperature, std::vector<double>& pressure) const;
                
                /**
                * This function returns pressure (Pa) of a pure ideal gas in the quadrature points of a mesh entity
                * at a variable temperature (non-isothermal case).
                *
                * @param density     - density in the quadrature points of a mesh entity [kg/m^3],
                * @param temperature - temperature in the quadrature points of a mesh entity [K],
                * @param pressure    - pressure in the quadrature points of a mesh entity [Pa].
                */
                void get_pressure(const std::vector<double>& density, const std::vector<double>& temperature, std::vector<double>& pressure) const;
                
                /**
                * This function returns
                * \f$ \quad \frac{\partial p}{\partial \rho} \quad \f$ of a pure ideal gas.
                *
                * @param temperature - temperature [K].
                */
                const double get_Dpressure_Ddensity(const double& temperature) const;
                
                /**
                * This function returns
                * \f$ \quad \frac{\partial p}{\partial \rho} \quad \f$ of a pure ideal gas
                * in the quadrature points of a mesh entity.
                *
                * @param temperature - temperature in the quadrature points of a mesh entity [K],
                * @param dst         - \f$ \frac{\partial p}{\partial \rho} \quad \f$ in the quadrature points of a mesh entity.
                */
                void get_Dpressure_Ddensity(const std::vector<double>& temperature, std::vector<double>& dst) const;
                
                /**
                * This function returns
                * \f$ \quad \frac{\partial p}{\partial T} \quad \f$ of a pure ideal gas.
                *
                * @param density - density [kg/m^3].
                */
                const double get_Dpressure_Dtemperature(const double& density) const;
                
                /**
                * This function returns
                * \f$ \quad \frac{\partial p}{\partial T} \quad \f$ of a pure ideal gas
                * in the quadrature points of a mesh entity.
                *
                * @param density - density in the quadrature points of a mesh entity,
                * @param dst     - \f$ \frac{\partial p}{\partial T} \quad \f$ in the quadrature points of a mesh entity.
                */
                void get_Dpressure_Dtemperature(const std::vector<double>& density, std::vector<double>& dst) const;
                
                //@}
                ///@name Service functions. Density.
                //@{
                
                /**
                * This function returns density [kg/m^3] of a pure ideal gas.
                *
                * @param temperature - temperature [K],
                * @param pressure - pressure [Pa],
                * @param molarMass - molar mass [kg/mol].
                */
                const double get_density(const double& temperature, const double& pressure, const double& molarMass) const;
                
                //@}
                ///@name Service functions. Sutherland dynamic viscosity.
                //@{
                
                /**
                * This function returns
                * Sutherland dynamic viscosity [kg/(m s)]
                * of a pure gas.
                *
                * @param temperature - temperature [K].
                */
                const double get_Sutherland_dynamic_viscosity(const double& temperature) const;
                
                /**
                * This function returns Sutherland dynamic viscosity [kg/(m s)] of a pure gas
                * in the quadrature points of a mesh entity.
                *
                * @param temperature       - temperature in the quadrature points of a mesh entity [K],
                * @param dynamic_viscosity - dynamic viscosity in the quadrature points of a mesh entity [kg/(m s)].
                */
                void get_Sutherland_dynamic_viscosity(const std::vector<double>& temperature, std::vector<double>& dynamic_viscosity) const;
                
                /**
                * This function returns the first derivative \f$ \quad \frac{\partial \mu}{\partial T} \quad \f$ of the Sutherland dynamic viscosity
                * of a pure gas.
                *
                * @param temperature - temperature [K].
                */
                const double get_DSutherland_dynamic_viscosity_Dtemperature(const double& temperature) const;
                
                /**
                * This function returns the first derivative \f$ \quad \frac{\partial \mu}{\partial T} \quad \f$ of the Sutherland dynamic viscosity
                * of a pure gas in the quadrature points of a mesh entity.
                *
                * @param temperature - temperature in the quadrature points of a mesh entity [K],
                * @param dst         - \f$ \frac{\partial \mu}{\partial T} \quad \f$ in the quadrature points of a mesh entity.
                */
                void get_DSutherland_dynamic_viscosity_Dtemperature(const std::vector<double>& temperature, std::vector<double>& dst) const;
                
                //@}
                ///@name Service functions. Chapman Enskog dynamic viscosity.
                //@{
                /**
                * This function returns Chapman Enskog dynamic viscosity [kg/(m s)] of a pure gas. Parameters:
                * @param temperature - temperature [K].
                */
                const double get_ChapmanEnskog_dynamic_viscosity(const double& temperature) const;
                
                /**
                * This function returns Chapman Enskog dynamic viscosity [kg/(m s)] of a pure gas in the quadrature points of a mesh entity.
                * Parameters:
                * @param temperature       - temperature in the quadrature points of a mesh entity [K],
                * @param dynamic_viscosity - dynamic viscosity in the quadrature points of a mesh entity [kg/(m s)].
                */
                void get_ChapmanEnskog_dynamic_viscosity(const std::vector<double>& temperature, std::vector<double>& dynamic_viscosity) const;
                
                /**
                * This function returns the first derivative \f$ \quad \frac{\partial \mu}{\partial T} \quad \f$ of the Chapman Enskog dynamic viscosity
                * of a pure gas.
                *
                * @param temperature - temperature [K].
                */
                const double get_DChapmanEnskog_dynamic_viscosity_Dtemperature(const double& temperature) const;
                
                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial \mu}{\partial T} \quad \f$ of the Chapman Enskog dynamic viscosity
                * of a pure gas
                * in the quadrature points of a mesh entity.
                *
                * @param temperature - temperature in the quadrature points of a mesh entity [K],
                * @param dst         - \f$ \frac{\partial \mu}{\partial T} \quad \f$ in the quadrature points of a mesh entity.
                */
                void get_DChapmanEnskog_dynamic_viscosity_Dtemperature(const std::vector<double>& temperature, std::vector<double>& dst) const;
                
                //@}
                ///@name Service functions. Dynamic viscosity.
                //@{
                
                /**
                * This function returns
                * desired dynamic viscosity [kg/(m s)]
                * of a pure gas.
                *
                * @param temperature - temperature [K].
                */
                const double get_dynamic_viscosity(const double& temperature) const;
                
                /**
                * This function returns
                * desired dynamic viscosity [kg/(m s)]
                * of a pure gas
                * in the quadrature points of a mesh entity.
                *
                * @param temperature       - temperature in the quadrature points of a mesh entity [K],
                * @param dynamic_viscosity - dynamic viscosity in the quadrature points of a mesh entity [kg/(m s)].
                */
                void get_dynamic_viscosity(const std::vector<double>& temperature, std::vector<double>& dynamic_viscosity) const;
                
                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial \mu}{\partial T} \quad \f$ of the desired dynamic viscosity
                * of a pure gas.
                *
                * @param temperature - temperature [K].
                */
                const double get_Ddynamic_viscosity_Dtemperature(const double& temperature) const;
                
                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial \mu}{\partial T} \quad \f$ of the desired dynamic viscosity
                * of a pure gas
                * in the quadrature points of a mesh entity.
                *
                * @param temperature - temperature in the quadrature points of a mesh entity [K],
                * @param dst         - \f$ \frac{\partial \mu}{\partial T} \quad \f$ in the quadrature points of a mesh entity.
                */
                void get_Ddynamic_viscosity_Dtemperature(const std::vector<double>& temperature, std::vector<double>& dst) const;
                
                //@}
                
                ///@name Service functions. Bulk viscosity.
                //@{
                
                /**
                * This function returns desired bulk viscosity of a pure gas.
                *
                * @param dynamic_viscosity - dynamic viscosity [kg/(m s)].
                */
                const double get_bulk_viscosity(const double& dynamic_viscosity) const;
                
                /**
                * This function returns desired bulk viscosity of a pure gas
                * in the quadrature points of a mesh entity.
                *
                * @param dynamic_viscosity - dynamic viscosity in the quadrature points of a mesh entity [kg/(m s)],
                * @param bulk_viscosity    - bulk viscosity in the quadrature points of a mesh entity [kg/(m s)].
                */
                void get_bulk_viscosity(const std::vector<double>& dynamic_viscosity, std::vector<double>& bulk_viscosity) const;
                
                /**
                * This function returns the first derivative \f$ \quad \frac{\partial \lambda}{\partial T} \quad \f$ of the desired bulk viscosity
                * of a pure gas.
                *
                * @param src - the first derivative \f$ \quad \frac{\partial \mu}{\partial T} \quad \f$ of the desired dynamic viscosity
                *              of a pure gas.
                */
                const double get_Dbulk_viscosity_Dtemperature(const double& src) const;
                
                /**
                * This function returns the first derivative \f$ \quad \frac{\partial \lambda}{\partial T} \quad \f$ of the desired bulk viscosity
                * of a pure gas in the quadrature points of a mesh entity.
                *
                * @param src - the first derivative \f$ \quad \frac{\partial \mu}{\partial T} \quad \f$ of the desired dynamic viscosity
                *              of a pure gas in the quadrature points of a mesh entity,
                * @param dst - the first derivative \f$ \quad \frac{\partial \lambda}{\partial T} \quad \f$ of the desired bulk viscosity
                *              of a pure gas in the quadrature points of a mesh entity.
                */
                void get_Dbulk_viscosity_Dtemperature(const std::vector<double>& src, std::vector<double>& dst) const;
                
                //@}
                ///@name Service functions. Partial viscosity of the mixture.
                //@{
                
                //@}
                ///@name Service functions. Sutherland thermal conductivity.
                //@{
                
                /**
                * This function returns Sutherland thermal conductivity of a pure gas.
                *
                * @param temperature - temperature [K].
                */
                const double get_Sutherland_thermal_conductivity(const double& temperature) const;
                
                /**
                * This function returns Sutherland thermal conductivity of a pure gas
                * in the quadrature points of a mesh entity.
                *
                * @param temperature          - temperature in the quadrature points of a mesh entity [K],
                * @param thermal_conductivity - thermal conductivity in the quadrature points of a mesh entity [W/(m K)].
                */
                void get_Sutherland_thermal_conductivity(const std::vector<double>& temperature, std::vector<double>& thermal_conductivity) const;
                
                /**
                * This function returns the first derivative \f$ \quad \frac{\partial \kappa}{\partial T} \quad \f$ of the Sutherland thermal conductivity
                * of a pure gas.
                *
                * @param temperature - temperature [K].
                */
                const double get_DSutherland_thermal_conductivity_Dtemperature(const double& temperature) const;
                
                /**
                * This function returns the first derivative \f$ \quad \frac{\partial \kappa}{\partial T} \quad \f$ of the Sutherland thermal conductivity
                * of a pure gas in the quadrature points of a mesh entity.
                *
                * @param temperature - temperature in the quadrature points of a mesh entity [K],
                * @param dst         - \f$ \frac{\partial \kappa}{\partial T} \quad \f$ in the quadrature points of a mesh entity.
                */
                void get_DSutherland_thermal_conductivity_Dtemperature(const std::vector<double>& temperature, std::vector<double>& dst) const;
                
                //@}
                
                ///@name Service functions. Chapman Enskog thermal conductivity.
                //@{
                
                /**
                * This function returns
                * Chapman Enskog thermal conductivity [W/(m K)]
                * of a pure gas.
                *
                * @param temperature - temperature [K].
                */
                const double get_ChapmanEnskog_thermal_conductivity(const double& temperature) const;
                
                /**
                * This function returns
                * Chapman Enskog thermal conductivity
                * of a pure gas
                * in the quadrature points of a mesh entity.
                *
                * @param temperature          - temperature in the quadrature points of a mesh entity [K],
                * @param thermal_conductivity - thermal conductivity in the quadrature points of a mesh entity [W/(m K)].
                */
                void get_ChapmanEnskog_thermal_conductivity(const std::vector<double>& temperature, std::vector<double>& thermal_conductivity) const;
                
                /**
                * This function returns the first derivative \f$ \quad \frac{\partial \kappa}{\partial T} \quad \f$ of the Chapman Enskog thermal conductivity
                * of a pure gas.
                *
                * @param temperature - temperature [K].
                */
                const double get_DChapmanEnskog_thermal_conductivity_Dtemperature(const double& temperature) const;
                
                /**
                * This function returns the first derivative \f$ \quad \frac{\partial \kappa}{\partial T} \quad \f$ of the Chapman Enskog thermal conductivity
                * of a pure gas in the quadrature points of a mesh entity.
                *
                * @param temperature - temperature in the quadrature points of a mesh entity [K],
                * @param dst         - \f$ \frac{\partial \kappa}{\partial T} \quad \f$ in the quadrature points of a mesh entity.
                */
                void get_DChapmanEnskog_thermal_conductivity_Dtemperature(const std::vector<double>& temperature, std::vector<double>& dst) const;
                
                //@}
                
                ///@name Service functions. Thermal conductivity.
                //@{
                
                /**
                * This function returns
                * desired thermal conductivity [W/(m K)]
                * of a pure gas.
                *
                * @param temperature - temperature [K].
                */
                const double get_thermal_conductivity(const double& temperature) const;
                
                /**
                * This function returns
                * desired thermal conductivity
                * of a pure gas
                * in the quadrature points of a mesh entity.
                *
                * @param temperature          - temperature in the quadrature points of a mesh entity [K],
                * @param thermal_conductivity - thermal conductivity in the quadrature points of a mesh entity [W/(m K)].
                */
                void get_thermal_conductivity(const std::vector<double>& temperature, std::vector<double>& thermal_conductivity) const;
                
                /**
                * This function returns the first derivative \f$ \quad \frac{\partial \kappa}{\partial T} \quad \f$ of the desired thermal conductivity
                * of a pure gas.
                *
                * @param temperature - temperature [K].
                */
                const double get_Dthermal_conductivity_Dtemperature(const double& temperature) const;
                
                /**
                * This function returns the first derivative \f$ \quad \frac{\partial \kappa}{\partial T} \quad \f$ of the desired thermal conductivity
                * of a pure gas in the quadrature points of a mesh entity.
                *
                * @param temperature - temperature in the quadrature points of a mesh entity [K],
                * @param dst         - \f$ \frac{\partial \kappa}{\partial T} \quad \f$ in the quadrature points of a mesh entity.
                */
                void get_Dthermal_conductivity_Dtemperature(const std::vector<double>& temperature, std::vector<double>& dst) const;
                
                //@}
                
                ///@name Service functions. Molar enthalpy.
                //@{
                
                /**
                * This function returns
                * molar enthalpy [J/mol]
                * of a pure gas.
                *
                * @param temperature - temperature [K].
                */
                const double get_molar_enthalpy(const double& temperature) const;
                
                /**
                * This function returns
                * molar enthalpy
                * of a pure gas
                * in the quadrature points of a mesh entity.
                *
                * @param temperature    - temperature in the quadrature points of a mesh entity [K],
                * @param molar_enthalpy - molar enthalpy in the quadrature points of a mesh entity [J/mol].
                */
                void get_molar_enthalpy(const std::vector<double>& temperature, std::vector<double>& molar_enthalpy) const;
                
                /**
                * This function returns the first derivative \f$ \quad \frac{\partial H}{\partial T} \quad \f$ of the molar enthalpy
                * of a pure gas.
                *
                * @param temperature - temperature [K].
                */
                const double get_Dmolar_enthalpy_Dtemperature(const double& temperature) const;
                
                /**
                * This function returns the first derivative \f$ \quad \frac{\partial H}{\partial T} \quad \f$ of the molar enthalpy
                * of a pure gas in the quadrature points of a mesh entity.
                *
                * @param temperature - temperature in the quadrature points of a mesh entity [K],
                * @param dst         - \f$ \frac{\partial H}{\partial T} \quad \f$ in the quadrature points of a mesh entity.
                */
                void get_Dmolar_enthalpy_Dtemperature(const std::vector<double>& temperature, std::vector<double>& dst) const;
                
                /**
                * This function returns the second derivative \f$ \quad \frac{\partial^2 H}{\partial T^2} \quad \f$ of the molar enthalpy
                * of a pure gas.
                *
                * @param temperature - temperature [K].
                */
                const double get_D2molar_enthalpy_Dtemperature2(const double& temperature) const;
                
                /**
                * This function returns
                * the second derivative \f$ \quad \frac{\partial^2 H}{\partial T^2} \quad \f$ of the molar enthalpy
                * of a pure gas
                * in the quadrature points of a mesh entity.
                *
                * @param temperature - temperature in the quadrature points of a mesh entity [K],
                * @param dst         - \f$ \frac{\partial^2 H}{\partial T^2} \quad \f$ in the quadrature points of a mesh entity.
                */
                void get_D2molar_enthalpy_Dtemperature2(const std::vector<double>& temperature,  std::vector<double>& dst) const;
                
                //@}
                
                ///@name Service functions. Water vapor saturation pressure.
                //@{
                
                /**
                * This function returns
                * saturation pressure [Pa]
                * of water vapor.
                *
                * @param temperature - temperature [K].
                */
                const double get_water_vapor_saturation_pressure(const double& temperature) const;
                
                /**
                * This function returns
                * saturation pressure
                * of water vapor
                * in the quadrature points of a mesh entity.
                *
                * @param temperature                     - temperature in the quadrature points of a mesh entity [K],
                * @param water_vapor_saturation_pressure - water vapor saturation pressure in the quadrature points of a mesh entity [Pa].
                */
                void get_water_vapor_saturation_pressure(const std::vector<double>& temperature, std::vector<double>& water_vapor_saturation_pressure) const;
                
                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial p_{\text{sat}}^{H_2O}}{\partial T} \quad \f$ of the saturation pressure
                * of water vapor.
                *
                * @param temperature - temperature [K].
                */
                const double get_Dwater_vapor_saturation_pressure_Dtemperature(const double& temperature) const;
                
                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial p_{\text{sat}}^{H_2O}}{\partial T} \quad \f$ of the saturation pressure
                * of water vapor
                * in the quadrature points of a mesh entity.
                *
                * @param temperature - temperature in the quadrature points of a mesh entity [K],
                * @param dst         - \f$ \frac{\partial p_{\text{sat}}^{H_2O}}{\partial T} \quad \f$ in the quadrature points of a mesh entity.
                */
                void get_Dwater_vapor_saturation_pressure_Dtemperature(const std::vector<double>& temperature, std::vector<double>& dst) const;
                
                //@}
                
                ///@name Service functions. Collision integral.
                //@{
                
                /**
                * This function returns
                * collision integral
                * of a pure gas.
                *
                * @param temperature - temperature [K].
                * 
                * @htmlonly
                *   <ul>
                *       <li> Neufeld, Philip D., A. R. Janzen, and R. A. Aziz. "Empirical Equations to Calculate 16 of the Transport Collision Integrals Ω (l, s)* for the Lennard‐Jones (12–6) Potential." The Journal of Chemical Physics 57.3 (1972): 1100-1102. </li>
                *   </ul>
                * @endhtmlonly
                */
                const double get_collision_integral(const double& temperature) const;
                
                /**
                * This function returns
                * collision integral
                * of a pure gas
                * in the quadrature points of a mesh entity.
                *
                * @param temperature        - temperature in the quadrature points of a mesh entity [K],
                * @param collision_integral - collision integral in the quadrature points of a mesh entity.
                * 
                * @htmlonly
                *   <ul>
                *       <li> Neufeld, Philip D., A. R. Janzen, and R. A. Aziz. "Empirical Equations to Calculate 16 of the Transport Collision Integrals Ω (l, s)* for the Lennard‐Jones (12–6) Potential." The Journal of Chemical Physics 57.3 (1972): 1100-1102. </li>
                *   </ul>
                * @endhtmlonly
                */
                void get_collision_integral(const std::vector<double>& temperature, std::vector<double>& collision_integral) const;
                
                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial \Omega_{\mu, \kappa}}{\partial T} \quad \f$ of the collision integral
                * of a pure gas.
                *
                * @param temperature - temperature [K].
                */
                const double get_Dcollision_integral_Dtemperature(const double& temperature) const;
                
                /**
                * This function returns the first derivative \f$ \quad \frac{\partial \Omega_{\mu, \kappa}}{\partial T} \quad \f$ of the collision integral
                * of a pure gas in the quadrature points of a mesh entity.
                *
                * @param temperature - temperature in the quadrature points of a mesh entity [K],
                * @param dst         - \f$ \frac{\partial \Omega_{\mu, \kappa}}{\partial T} \quad \f$ in the quadrature points of a mesh entity.
                */
                void get_Dcollision_integral_Dtemperature(const std::vector<double>& temperature, std::vector<double>& dst) const;
                
                //@}
                
                ///@name Service functions. Specific heat capacity at constant pressure.
                //@{
                
                /**
                * This function returns specific heat capacity [J/(g C)] at constant pressure of a pure gas.
                *
                * @param temperature - temperature [K].
                */
                const double get_specific_heat_capacity(const double& temperature) const;
                
                /**
                * This function returns specific heat capacity at constant pressure
                * of a pure gas in the quadrature points of a mesh entity.
                *
                * @param temperature            - temperature in the quadrature points of a mesh entity [K],
                * @param specific_heat_capacity - specific heat capacity at constant pressure in the quadrature points of a mesh entity [J/(g C)].
                */
                void get_specific_heat_capacity(const std::vector<double>& temperature, std::vector<double>& specific_heat_capacity) const;
                
                /**
                * This function returns the first derivative \f$ \quad \frac{\partial c_p}{\partial T} \quad \f$ of the specific heat capacity at constant pressure
                * of a pure gas.
                *
                * @param temperature - temperature [K].
                */
                const double get_Dspecific_heat_capacity_Dtemperature(const double& temperature) const;
                
                /**
                * This function returns the first derivative \f$ \quad \frac{\partial c_p}{\partial T} \quad \f$ of the specific heat capacity at constant pressure
                * of a pure gas in the quadrature points of a mesh entity.
                *
                * @param temperature - temperature in the quadrature points of a mesh entity [K],
                * @param dst         - \f$ \frac{\partial c_p}{\partial T} \quad \f$ in the quadrature points of a mesh entity.
                */
                void get_Dspecific_heat_capacity_Dtemperature(const std::vector<double>& temperature, std::vector<double>& dst) const;
                
                /**
                * This function returns the second derivative \f$ \quad \frac{\partial^2 c_p}{\partial T^2} \quad \f$ of the specific heat capacity at constant pressure
                * of a pure gas.
                *
                * @param temperature - temperature [K].
                */
                const double get_D2specific_heat_capacity_Dtemperature2(const double& temperature) const;
                
                /**
                * This function returns the second derivative \f$ \quad \frac{\partial^2 c_p}{\partial T^2} \quad \f$ of the specific heat capacity at constant pressure
                * of a pure gas in the quadrature points of a mesh entity.
                *
                * @param temperature - temperature in the quadrature points of a mesh entity,
                * @param dst         - \f$ \frac{\partial^2 c_p}{\partial T^2} \quad \f$ in the quadrature points of a mesh entity.
                */
                void get_D2specific_heat_capacity_Dtemperature2(const std::vector<double>& temperature, std::vector<double>& dst) const;
                //@}
                
            protected:
                ///@name Constructors, destructor, and initialization
                //@{
                /**
                * Constructor.
                */
                PureGas(const std::string& name);
                
                //@}
                
                //////////
                // DATA //
                //////////
                
                ///@name Fluid properties
                //@{
                
                /**
                * Molar mass [kg/mol], \f$ M \quad \left[ \frac{\text{kg}}{\text{mol}} \right] \f$.
                */
                double molar_mass;
                
                /**
                * Collision diameter [Angstrom], \f$ \sigma \quad \left[ \text{Angstrom} \right] \f$.
                */
                double collision_diameter;
                
                /**
                * The maximum energy of attraction divided by the Boltzmann constant, \f$ \frac{\epsilon}{k} \quad \left[ \text{K} \right] \f$.
                */
                double eps_BY_k;
                
                /**
                * Prandtl number, \f$ \text{Pr} \f$.
                */
                double Prandtl;
                
                //@}
                
                ///@name Coefficients
                //@{
                
                /**
                * Coefficient of the Sutherland
                * dynamic viscosity formula, \f$ A \quad \left[ \frac{\text{Pa sec}}{\text{K}^{1/2}} \right] \f$.
                */
                double A_Sutherland;
                
                /**
                * Coefficient of the Sutherland
                * dynamic viscosity formula, \f$ B \quad \left[ \text{K} \right] \f$.
                */
                double B_Sutherland;
                
                /**
                * Coefficient of the specific heat capacity at constant pressure formula,
                * \f$ c_0 \quad \left[ \frac{\text{J}}{\text{kg K}} \right] \f$.
                */
                double c_0;
                
                /**
                * Coefficient of the specific heat capacity at constant pressure formula,
                * \f$ c_1 \quad \left[ \frac{\text{J}}{\text{kg } \text{K}^2} \right] \f$.
                */
                double c_1;
                
                /**
                * Coefficient of the specific heat capacity at constant pressure formula,
                * \f$ c_2 \quad \left[ \frac{\text{J}}{\text{kg } \text{K}^3} \right] \f$.
                */
                double c_2;
                
                /**
                * Coefficient of the specific heat capacity at constant pressure formula,
                * \f$ c_3 \quad \left[ \frac{\text{J}}{\text{kg } \text{K}^4} \right] \f$.
                */
                double c_3;
                
                /**
                * Coefficient of the molar enthalpy formula
                * \f$ H^{\text{ref}} \quad \left[ \frac{\text{J}}{\text{mol}} \right] \f$.
                */
                double H_ref;
                
                /**
                * Coefficient of the molar enthalpy formula
                * \f$ T^{\text{ref}} \quad \left[ \text{K} \right] \f$.
                */
                double T_ref;
                
                //@}
                
                ///@name Names and modes
                //@{
                
                /**
                * Chemical formula of the pure gas.
                */
                std::string chemical_formula;
                
                /**
                * @param dynamic_viscosity_mode=Sutherland    - Sutherland dynamic viscosity formula is used,
                * @param dynamic_viscosity_mode=ChapmanEnskog - Chapman Enskog dynamic viscosity formula is used.
                */
                std::string dynamic_viscosity_mode;
                
                /**
                * @param bulk_viscosity_mode=zero   - zero bulk viscosity formula is used,
                * @param bulk_viscosity_mode=Stokes - Stokes bulk viscosity formula is used.
                */
                std::string bulk_viscosity_mode;
                
                /**
                * @param thermal_conductivity_mode=Sutherland    - Sutherland thermal conductivity formula is used,
                * @param thermal_conductivity_mode=ChapmanEnskog - Chapman Enskog thermal conductivity formula is used.
                */
                std::string thermal_conductivity_mode;
                //@} 
        };
        
        /**
         * This class describes properties of a dummy gas, with all properties of 1.0
         * Use this name in the parameters file if needed.
         * 
         * @param molar_mass         - [kg/mol]; 
         * @param collision_diameter - [Angstrom]; 
         * @param eps_BY_k           - [K]; 
         * @param SpecificHeat       - [J/kg K]; 
         * 
         * \author Chad Balen, 2015
         */
        class DummyGas : public PureGas
        {
        public:
            ///@name Constructors, destructor, and initialization
            //@{
            
            /**
             * Constructor.
             */
            DummyGas()
            :
            PureGas("DummyGas")
            {
                this->molar_mass         = 1.0;
                this->collision_diameter = 1.0;
                this->eps_BY_k           = 1.0;
                this->Prandtl            = 1.0;
                
                this->A_Sutherland = 1.0;
                this->B_Sutherland = 1.0;
                
                this->c_0 =   1.0;
                this->c_1 = - 1.0;
                this->c_2 =   1.0;
                this->c_3 = - 1.0;
                
                this->H_ref = 1.0;
                this->T_ref = 1.0;
                
                this->chemical_formula = "None";
            }
            
            /**
             * Destructor.
             */
            virtual ~DummyGas() { }
            
            //@}
            
            ///@name Accessors and info
            //@{
            
            virtual enMaterialID get_ID() const
            {
                return matNone;
            }
            
            //@}
            
        };
        
        /**
         * This class describes properties of pure \p oxygen.
         * Use this name in the parameters file if needed.
         * 
         * @param molar_mass         - [kg/mol]; Poling, Bruce E., et al. The properties of gases and liquids. Vol. 5. New York: McGraw-Hill, 2001.
         * @param collision_diameter - [Angstrom]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
         * @param eps_BY_k           - [K]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
         * @param SpecificHeat       - [J/kg K]; 
         * 
         * \author M. Secanell and Valentin N. Zingan, 2013
         */
        class Oxygen : public PureGas
        {
            public:
                ///@name Constructors, destructor, and initialization
                //@{
                /** Constructor */
                Oxygen()
                :
                PureGas("oxygen")
                {
                    this->molar_mass         = 31.999e-3; // [kg/mol]; Poling, Bruce E., et al. The properties of gases and liquids. Vol. 5. New York: McGraw-Hill, 2001.
                    this->collision_diameter = 3.4330000; // [Angstrom]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
                    this->eps_BY_k           = 113.00000; // [K]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
                    this->Prandtl            = 0.7130000;
                    
                    this->A_Sutherland = 1.693411300e-6;
                    this->B_Sutherland = 127.0000000000;
                    
                    this->c_0 =   0.8800e3;
                    this->c_1 = - 0.0001e3;
                    this->c_2 =   0.5400e3;
                    this->c_3 = - 0.3300e3;
                    
                    this->H_ref = 0.000;
                    this->T_ref = 298.0;
                    
                    this->chemical_formula = "O2";
                }
                
                /** Destructor */
                virtual ~Oxygen() { }
                //@}
                ///@name Accessors and info
                //@{
                virtual enMaterialID get_ID() const
                {
                    return matOxygen;
                }
                //@}
        };
        
        /**
         * This class describes properties of pure \p nitrogen.
         * Use this name in the parameters file if needed.
         * 
         * @param molar_mass         - [kg/mol]; Mattern, Glenn W. "Industrial Gases." Riegel’s Handbook of Industrial Chemistry. Springer US, 1992. 442-457.
         * @param collision_diameter - [Angstrom]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
         * @param eps_BY_k           - [K]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
         * @param SpecificHeat       - [J/kg K]; 
         * 
         * \author M. Secanell and Valentin N. Zingan, 2013
         */
        
        class Nitrogen : public PureGas
        {
            public:
                ///@name Constructors, destructor, and initialization
                //@{
                /** Constructor */
                Nitrogen()
                :
                PureGas("nitrogen")
                {
                    this->molar_mass         = 28.0134e-3; // [kg/mol]; Mattern, Glenn W. "Industrial Gases." Riegel’s Handbook of Industrial Chemistry. Springer US, 1992. 442-457.
                    this->collision_diameter = 3.66700000; // [Angstrom]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
                    this->eps_BY_k           = 99.8000000; // [K]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
                    this->Prandtl            = 0.70700000;
                    
                    this->A_Sutherland = 1.406732195e-6;
                    this->B_Sutherland = 111.0000000000;
                    
                    this->c_0 =   1.11e3;
                    this->c_1 = - 0.48e3;
                    this->c_2 =   0.96e3;
                    this->c_3 = - 0.42e3;
                    
                    this->H_ref = 0.000;
                    this->T_ref = 298.0;
                    
                    this->chemical_formula = "N2";
                }
                
                /** Destructor */
                virtual ~Nitrogen() { }
                //@}
                ///@name Accessors and info
                //@{
                
                virtual enMaterialID get_ID() const
                {
                    return matNitrogen;
                }
                //@}
        };
        
        /**
         * This class describes properties of pure \p hydrogen.
         * Use this name in the parameters file if needed.
         * 
         * @param molar_mass         - [kg/mol]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
         * @param collision_diameter - [Angstrom]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
         * @param eps_BY_k           - [K]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
         * @param SpecificHeat       - [J/kg K]; 
         * 
         * \author M. Secanell and Valentin N. Zingan, 2013
         */
        
        class Hydrogen : public PureGas
        {
            public:
                ///@name Constructors, destructor, and initialization
                //@{
                /** Constructor */
                Hydrogen()
                :
                PureGas("hydrogen")
                {
                    this->molar_mass         = 2.016e-3; // [kg/mol]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
                    this->collision_diameter = 2.915000; // [Angstrom]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
                    this->eps_BY_k           = 38.00000; // [K]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
                    this->Prandtl            = 0.685000;
                    
                    this->A_Sutherland = 6.362365620e-7;
                    this->B_Sutherland = 72.00000000000;
                    
                    this->c_0 =   13.46e3;
                    this->c_1 =   4.600e3;
                    this->c_2 = - 6.850e3;
                    this->c_3 =   3.790e3;
                    
                    this->H_ref = 0.000;
                    this->T_ref = 298.0;
                    
                    this->chemical_formula = "H2";
                }
                
                /** Destructor */
                virtual ~Hydrogen() { }
                //@}
                ///@name Accessors and info
                //@{
                virtual enMaterialID get_ID() const
                {
                    return matHydrogen;
                }
                //@}
        };
        
        /**
         * This class describes properties of pure \p WaterVapor.
         * Use this name in the parameters file if needed.
         * 
         * @param molar_mass         - [kg/mol]; Poling, Bruce E., et al. The properties of gases and liquids. Vol. 5. New York: McGraw-Hill, 2001.
         * @param collision_diameter - [Angstrom]; Poling, Bruce E., et al. The properties of gases and liquids. Vol. 5. New York: McGraw-Hill, 2001.
         * @param eps_BY_k           - [K]; Poling, Bruce E., et al. The properties of gases and liquids. Vol. 5. New York: McGraw-Hill, 2001.
         * @param SpecificHeat       - [J/kg K]; 
         * 
         * \author M. Secanell and Valentin N. Zingan, 2013
         */
        
        class WaterVapor : public PureGas
        {
            public:
                ///@name Constructors, destructor, and initialization
                //@{
                /** Constructor */
                WaterVapor()
                :
                PureGas("water")
                {
                    this->molar_mass         = 18.015e-3; // [kg/mol]; Poling, Bruce E., et al. The properties of gases and liquids. Vol. 5. New York: McGraw-Hill, 2001.
                    this->collision_diameter = 2.6410000; // [Angstrom]; Poling, Bruce E., et al. The properties of gases and liquids. Vol. 5. New York: McGraw-Hill, 2001.
                    this->eps_BY_k           = 809.10000; // [K]; Poling, Bruce E., et al. The properties of gases and liquids. Vol. 5. New York: McGraw-Hill, 2001.
                    this->Prandtl            = 0.9500000;
                    
                    this->A_Sutherland = _DUMMY_;
                    this->B_Sutherland = _DUMMY_;
                    
                    this->c_0 =   1.790e3;
                    this->c_1 =   0.107e3;
                    this->c_2 =   0.586e3;
                    this->c_3 = - 0.200e3;
                    
                    this->H_ref = - 2.41826e5;
                    this->T_ref =   298.00000;
                    
                    this->chemical_formula = "H2O Vapor";
                }
                
                /** Destructor */
                virtual ~WaterVapor() { }
                //@}
                ///@name Accessors and info
                //@{
                virtual enMaterialID get_ID() const
                {
                    return matWaterVapor;
                }
                //@}
        };
        
        /**
         * This class describes properties of pure \p air.
         * Use this name in the parameters file if needed.
         * 
         * @param molar_mass         - [kg/mol]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
         * @param collision_diameter - [Angstrom]; The Properties of Gases and Liquids by Poling et al. New York: McGraw-Hill, 2001.
         * @param eps_BY_k           - [K]; The Properties of Gases and Liquids by Poling et al. New York: McGraw-Hill, 2001.
         * @param SpecificHeat       - [J/kg K]; 
         * 
         * \author M. Secanell and Valentin N. Zingan, 2013
         */
        
        class Air : public PureGas
        {
            public:
                ///@name Constructors, destructor, and initialization
                //@{
                /** Constructor */
                Air()
                :
                PureGas("air")
                {
                    this->molar_mass         = 28.97e-3;  // [kg/mol]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
                    this->collision_diameter = 3.7110000; // [Angstrom]; The Properties of Gases and Liquids by Poling et al. New York: McGraw-Hill, 2001.
                    this->eps_BY_k           = 78.600000; // [K]; The Properties of Gases and Liquids by Poling et al. New York: McGraw-Hill, 2001.
                    this->Prandtl            = 0.7020000;
                    
                    this->A_Sutherland = 1.512041288e-6;
                    this->B_Sutherland = 120.0000000000;
                    
                    this->c_0 =   1.050e3;
                    this->c_1 = - 0.365e3;
                    this->c_2 =   0.850e3;
                    this->c_3 = - 0.390e3;
                    
                    this->H_ref = _DUMMY_;
                    this->T_ref = _DUMMY_;
                    
                    this->chemical_formula = "78.1_N2__20.9_O2__0.9_Ar__0.03_CO2__0.002_Ne__0.0005_He__0.0002_CH4__0.0001_Kr__0.00005_H2__0.000009_He";
                }
                
                /** Destructor */
                virtual ~Air() { }
                //@}
                ///@name Accessors and info
                //@{
                virtual enMaterialID get_ID() const
                {
                    return matAir;
                }
                //@}
        };
        
        /**
         * This class describes properties of pure \p helium.
         * Use this name in the parameters file if needed.
         * 
         * @param molar_mass         - [kg/mol]; Coplen, Tyler B. "Atomic weights of the elements 1999." Journal of Physical and Chemical Reference Data 30.3 (2001): 701-712.
         * @param collision_diameter - [Angstrom]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
         * @param eps_BY_k           - [K]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
         * @param SpecificHeat       - [J/kg K]; 
         * 
         * \author M. Secanell and Valentin N. Zingan, 2013
         */
        
        class Helium : public PureGas
        {
            public:
                ///@name Constructors, destructor, and initialization
                //@{
                /** Constructor */
                Helium()
                :
                PureGas("helium")
                {
                    this->molar_mass         = 4.002602e-3; // [kg/mol]; Coplen, Tyler B. "Atomic weights of the elements 1999." Journal of Physical and Chemical Reference Data 30.3 (2001): 701-712.
                    this->collision_diameter = 2.576000000; // [Angstrom]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
                    this->eps_BY_k           = 10.20000000; // [K]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
                    this->Prandtl            = 0.710000000;
                    
                    this->A_Sutherland = 1.484381490e-6;
                    this->B_Sutherland = 79.40000000000;
                    
                    this->c_0 = 5.193e3;
                    this->c_1 = 0.00000;
                    this->c_2 = 0.00000;
                    this->c_3 = 0.00000;
                    
                    this->H_ref = _DUMMY_;
                    this->T_ref = _DUMMY_;
                    
                    this->chemical_formula = "He";
                }
                
                /** Destructor */
                virtual ~Helium() { }
                //@}
                ///@name Accessors and info
                //@{
                virtual enMaterialID get_ID() const
                {
                    return matHelium;
                }
                //@}
        };
        
        /**
         * This class describes properties of pure \p Argon.
         * Use this name in the parameters file if needed.
         * 
         * @param molar_mass         - [kg/mol]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
         * @param collision_diameter - [Angstrom]; Hirschfelder, Joseph O., et al. Molecular theory of gases and liquids. New York: Wiley, 1964.
         * @param eps_BY_k           - [K]; Hirschfelder, Joseph O., et al. Molecular theory of gases and liquids. New York: Wiley, 1964.
         * @param Prandtl            - at 300K; http://www.yutopian.com/Yuan/prop/Ar.html
         * @param A_Sutherland       - Reference at T = 300K; data from http://www.nist.gov/pml/div685/grp02/srd_134_helium.cfm, L2 norm from 100K - 1000K is 1.28639595331963e-6
         * @param B_Sutherland       - Reference at T = 300K; data from http://www.nist.gov/pml/div685/grp02/srd_134_helium.cfm, L2 norm from 100K - 1000K is 1.28639595331963e-6
         * @param SpecificHeat       - [J/kg K]; http://www.engineeringtoolbox.com/specific-heat-capacity-gases-d_159.html
         *
         * \author C. Balen, 2015
         */
        
        class Argon : public PureGas
        {
            public:
                ///@name Constructors, destructor, and initialization
                //@{
                /** Constructor */
                Argon()
                :
                PureGas("argon")
                {
                    this->molar_mass         = 39.94800e-3; // [kg/mol], Transport Phenomena by Bird
                    this->collision_diameter = 3.418000000; // [Angstrom], Molecular Theory of Gases and Liquids by Hirschfelder 
                    this->eps_BY_k           = 124.0000000; // [K], Molecular Theory of Gases and Liquids by Hirschfelder 
                    this->Prandtl            = 0.680000000; // @ 300K from http://www.yutopian.com/Yuan/prop/Ar.html
                    
                    //Reference at T = 300K, data from http://www.nist.gov/pml/div685/grp02/srd_134_helium.cfm, L2 norm from 100K - 1000K is 1.28639595331963e-6
                    this->A_Sutherland = 2.03988705213516e-6; // [Pa s K^-1/2], lambda
                    this->B_Sutherland = 164.89316; // [K], C, Sutherland's constant for the gaseous material in question
                    
                    this->c_0 = 0.520e3; // [J/(kg K)], http://www.engineeringtoolbox.com/specific-heat-capacity-gases-d_159.html
                    this->c_1 = 0.00000;
                    this->c_2 = 0.00000;
                    this->c_3 = 0.00000;
                    
                    this->H_ref = _DUMMY_;
                    this->T_ref = _DUMMY_;
                    
                    this->chemical_formula = "Ar";
                }
                
                /** Destructor */
                virtual ~Argon() { }
                //@}
                ///@name Accessors and info
                //@{
                virtual enMaterialID get_ID() const
                {
                    return matArgon;
                }
                //@}
        };
        
        /**
         * This class describes properties of pure \p Neon.
         * Use this name in the parameters file if needed.
         * 
         * @param molar_mass         - [kg/mol]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
         * @param collision_diameter - [Angstrom]; Hirschfelder, Joseph O., et al. Molecular theory of gases and liquids. New York: Wiley, 1964.
         * @param eps_BY_k           - [K]; Hirschfelder, Joseph O., et al. Molecular theory of gases and liquids. New York: Wiley, 1964.
         * @param Prandtl            - at ~100K; http://www.yutopian.com/Yuan/prop/Ar.html
         * @param A_Sutherland       - Reference at T = 373.15K; data from http://www.kayelaby.npl.co.uk/general_physics/2_2/2_2_3.html, L2 norm from 273.15K - 873.15K is 0.85079117060209e-6
         * @param B_Sutherland       - Reference at 373.15K; data from http://www.kayelaby.npl.co.uk/general_physics/2_2/2_2_3.html, L2 norm from 273.15K - 873.15K is 0.85079117060209e-6
         * @param SpecificHeat       - [J/kg K]; http://www.engineeringtoolbox.com/specific-heat-capacity-gases-d_159.html
         * 
         * \author C. Balen, 2015
         */
        
        class Neon : public PureGas
        {
            public:
                ///@name Constructors, destructor, and initialization
                //@{
                /** Constructor */
                Neon()
                :
                PureGas("neon")
                {
                    this->molar_mass         = 20.18300e-3; // [kg/mol], Transport Phenomena by Bird
                    this->collision_diameter = 2.789000000; // [Angstrom], Molecular Theory of Gases and Liquids by Hirschfelder 
                    this->eps_BY_k           = 35.70000000; // [K], Molecular Theory of Gases and Liquids by Hirschfelder 
                    this->Prandtl            = 0.776171079; // @ ~100C
                    
                    //Reference at T = 373.15K, data from http://www.kayelaby.npl.co.uk/general_physics/2_2/2_2_3.html, L2 norm from 273.15K - 873.15K is 0.85079117060209e-6
                    this->A_Sutherland = 2.3955224007011e-6;  // [Pa s K^-1/2], lambda
                    this->B_Sutherland = 93.534827; // [K], C, Sutherland's constant for the gaseous material in question
                    
                    this->c_0 = 1.030e3; // [J/(kg K)], http://www.engineeringtoolbox.com/specific-heat-capacity-gases-d_159.html
                    this->c_1 = 0.00000;
                    this->c_2 = 0.00000;
                    this->c_3 = 0.00000;
                    
                    this->H_ref = _DUMMY_;
                    this->T_ref = _DUMMY_;
                    
                    this->chemical_formula = "Ne";
                }
                
                /** Destructor */
                virtual ~Neon() { }
                //@}
                ///@name Accessors and info
                //@{
                virtual enMaterialID get_ID() const
                {
                    return matNeon;
                }
                //@}
        };
        
        /**
         * This class describes properties of pure \p Acetone.
         * Use this name in the parameters file if needed.
         * 
         * @param molar_mass         - [kg/mol]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
         * @param collision_diameter - [Angstrom]; The Properties of Gases and Liquids by Poling et al. New York: McGraw-Hill, 2001.
         * @param eps_BY_k           - [K]; The Properties of Gases and Liquids by Poling et al. New York: McGraw-Hill, 2001.
         * @param Prandtl            - 
         * @param A_Sutherland       - 
         * @param B_Sutherland       - 
         * @param SpecificHeat       - 
         * 
         * \author C. Balen, 2015
         */
        
        class Acetone : public PureGas
        {
            public:
                ///@name Constructors, destructor, and initialization
                //@{
                /** Constructor */
                Acetone()
                :
                PureGas("acetone")
                {
                    this->molar_mass         = 58.07800e-3; // [kg/mol], Transport Phenomena by Bird
                    this->collision_diameter = 4.600000000; // [Angstrom], The Properties of Gases and Liquids by Poling et al.
                    this->eps_BY_k           = 560.2000000; // [K], The Properties of Gases and Liquids by Poling et al.
                    this->Prandtl            = 0.000000000; //
                    
                    //
                    this->A_Sutherland = 0.0;  // [Pa s K^-1/2], lambda
                    this->B_Sutherland = 0.0; // [K], C, Sutherland's constant for the gaseous material in question
                    
                    this->c_0 = 0.00000; // 
                    this->c_1 = 0.00000;
                    this->c_2 = 0.00000;
                    this->c_3 = 0.00000;
                    
                    this->H_ref = _DUMMY_;
                    this->T_ref = _DUMMY_;
                    
                    this->chemical_formula = "C3H6O";
                }
                
                /** Destructor */
                virtual ~Acetone() { }
                //@}
                ///@name Accessors and info
                //@{
                virtual enMaterialID get_ID() const
                {
                    return matAcetone;
                }
                //@}
        };
        
        /**
         * This class describes properties of pure \p Methanol.
         * Use this name in the parameters file if needed.
         * 
         * @param molar_mass         - [kg/mol]; Bird, R. Byron. "Transport phenomena." Applied Mechanics Reviews 55.1 (2002): R1-R4.
         * @param collision_diameter - [Angstrom]; The Properties of Gases and Liquids by Poling et al. New York: McGraw-Hill, 2001.
         * @param eps_BY_k           - [K]; The Properties of Gases and Liquids by Poling et al. New York: McGraw-Hill, 2001.
         * @param Prandtl            - 
         * @param A_Sutherland       - 
         * @param B_Sutherland       - 
         * @param SpecificHeat       - 
         * 
         * \author C. Balen, 2015
         */
        
        class Methanol : public PureGas
        {
            public:
                ///@name Constructors, destructor, and initialization
                //@{
                /** Constructor */
                Methanol()
                :
                PureGas("methanol")
                {
                    this->molar_mass         = 32.04200e-3; // [kg/mol], Transport Phenomena by Bird
                    this->collision_diameter = 3.626000000; // [Angstrom], The Properties of Gases and Liquids by Poling et al.
                    this->eps_BY_k           = 481.8000000; // [K], The Properties of Gases and Liquids by Poling et al.
                    this->Prandtl            = 0.000000000; //
                    
                    //
                    this->A_Sutherland = 0.0;  // [Pa s K^-1/2], lambda
                    this->B_Sutherland = 0.0; // [K], C, Sutherland's constant for the gaseous material in question
                    
                    this->c_0 = 0.00000; //
                    this->c_1 = 0.00000;
                    this->c_2 = 0.00000;
                    this->c_3 = 0.00000;
                    
                    this->H_ref = _DUMMY_;
                    this->T_ref = _DUMMY_;
                    
                    this->chemical_formula = "CH3OH";
                }
                
                /** Destructor */
                virtual ~Methanol() { }
                //@}
                ///@name Accessors and info
                //@{
                virtual enMaterialID get_ID() const
                {
                    return matMethanol;
                }
                //@}
        };
    } // Material             
} // FuelCellShop         

#endif