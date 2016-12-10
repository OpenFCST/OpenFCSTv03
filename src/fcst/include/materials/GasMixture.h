// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: GasMixture.h
// - Description: This class describes properties of gas mixtures
// - Developers: Valentin N. Zingan, University of Alberta
//
// ----------------------------------------------------------------------------

#ifndef _FCST_FUELCELLSHOP_MATERIAL_GASMIXTURE_H_
#define _FCST_FUELCELLSHOP_MATERIAL_GASMIXTURE_H_

#include <materials/PureGas.h>

namespace FuelCellShop
{
    namespace Material
    {

        /**
        * This class describes properties of gas mixtures.
        * 
        * \note All those methods receive data in SI units.
        * \note All those methods return the results in SI units.
        *
        * This class contains the following data:
        *
        * - \p std::vector of pointers to the \p PureGas objects called \p gases,
        * - total pressure of the whole gas mixture called \p total_pressure,
        * - temperature of the whole gas mixture called \p temperature.
        *
        * The whole gas mixture is supposed to be isobaric if some concrete value of \p total_pressure
        * is assigned by either using the \p set_total_pressure() function or defining it in the parameters file.
        *
        * If nothing happens, then the value of \p total_pressure gets the \p _DUMMY_ number equal
        * \p 1.e300 and the whole gas mixture is treated as non-isobaric.
        * In this case, the total pressure of the whole gas mixture is
        * one of the solution variables.
        *
        * The whole gas mixture is supposed to be isothermal if some concrete value of \p temperature 
        * is assigned by either using the \p set_temperature() function or defining it in the parameters file.
        *
        * If nothing happens, then the value of \p temperature gets the \p _DUMMY_ number equal
        * \p 1.e300 and the whole gas mixture is treated as non-isothermal.
        * In this case, the temperature of the whole gas mixture is one of the solution variables.
        *
        * NOTE: after calling \p initialize() please call \p set_total_pressure() 
        *       and \p set_temperature() to insert the pressure and temperature
        *       from the Operating Conditions class if fluid flow is NOT isobaric
        *       or isothermal respectively. See these function calls below for 
        *       more detail.
        * 
        * The following methods are used to compute the gas mixture properties:
        *
        * - Maxwell-Stefan isobaric diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
        *   written in the Chapman Enskog form (binary gas mixture only), \f$ \quad \left[ \frac{\text{Pa } \text{m}^2}{\text{sec}} \right] \quad \f$
        *   \f$ D_{12} \equiv p \mathscr{D}_{12} = 1.8829 \cdot 10^{-2} \frac{\sqrt{ \left( \frac{1}{M_1} + \frac{1}{M_2} \right) T^3}}{\sigma_{12}^2 \Omega_{\mathscr{D}, 12}} \quad \f$
        *   and partial derivative \f$ \quad \frac{\partial D_{12}}{\partial T} \quad \f$
        *   where the binary collision integral \f$ \quad \Omega_{\mathscr{D}, 12} \quad \f$ is given by
        *   \f$ \quad \Omega_{\mathscr{D}, 12} = \frac{A_{\text{diff}}}{\left(\frac{kT}{\epsilon_{12}}\right)^{B_{\text{diff}}}} +
        *                                        \frac{C_{\text{diff}}}{e^{D_{\text{diff}}\frac{kT}{\epsilon_{12}}}} +
        *                                        \frac{E_{\text{diff}}}{e^{F_{\text{diff}}\frac{kT}{\epsilon_{12}}}} +
        *                                        \frac{G_{\text{diff}}}{e^{H_{\text{diff}}\frac{kT}{\epsilon_{12}}}} \f$
        *
        * - Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
        *   written in the Chapman Enskog form (binary gas mixture only), \f$ \quad \left[ \frac{\text{m}^2}{\text{sec}} \right] \quad \f$
        *   \f$ \mathscr{D}_{12} = 1.8583 \cdot 10^{-7} \frac{\sqrt{ \left( \frac{1}{M_1} + \frac{1}{M_2} \right) T^3}}{p_{\text{total}} \sigma_{12}^2 \Omega_{\mathscr{D}, 12}} \quad \f$
        *   and partial derivatives \f$ \quad \frac{\partial \mathscr{D}_{12}}{\partial T} \quad \f$ and \f$ \quad \frac{\partial \mathscr{D}_{12}}{\partial p} \quad \f$
        *   where the binary collision integral \f$ \quad \Omega_{\mathscr{D}, 12} \quad \f$ is given above
        *
        * - The table \f$ N_{\text{gases}} \times N_{\text{gases}} \f$ containing
        *   Maxwell-Stefan isobaric diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
        *   written in the Chapman Enskog form (ternary and more complicated gas mixtures), \f$ \quad \left[ \frac{\text{Pa } \text{m}^2}{\text{sec}} \right] \quad \f$
        *   \f$ D_{ij} \equiv p \mathscr{D}_{ij} = 1.8829 \cdot 10^{-2} \frac{\sqrt{ \left( \frac{1}{M_i} + \frac{1}{M_j} \right) T^3}}{\sigma_{ij}^2 \Omega_{\mathscr{D}, ij}} \quad \f$
        *   and partial derivatives \f$ \quad \frac{\partial D_{ij}}{\partial T} \quad \f$
        *   where the binary collision integral \f$ \quad \Omega_{\mathscr{D}, ij} \quad \f$ is given by
        *   \f$ \quad \Omega_{\mathscr{D}, ij} = \frac{A_{\text{diff}}}{\left(\frac{kT}{\epsilon_{ij}}\right)^{B_{\text{diff}}}} +
        *                                        \frac{C_{\text{diff}}}{e^{D_{\text{diff}}\frac{kT}{\epsilon_{ij}}}} +
        *                                        \frac{E_{\text{diff}}}{e^{F_{\text{diff}}\frac{kT}{\epsilon_{ij}}}} +
        *                                        \frac{G_{\text{diff}}}{e^{H_{\text{diff}}\frac{kT}{\epsilon_{ij}}}} \f$
        *
        * - The table \f$ N_{\text{gases}} \times N_{\text{gases}} \f$ containing
        *   Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
        *   written in the Chapman Enskog form (ternary and more complicated gas mixtures), \f$ \quad \left[ \frac{\text{m}^2}{\text{sec}} \right] \quad \f$
        *   \f$ \mathscr{D}_{ij} = 1.8583 \cdot 10^{-7} \frac{\sqrt{ \left( \frac{1}{M_i} + \frac{1}{M_j} \right) T^3}}{p_{\text{total}} \sigma_{ij}^2 \Omega_{\mathscr{D}, ij}} \quad \f$
        *   and partial derivatives \f$ \quad \frac{\partial \mathscr{D}_{ij}}{\partial T} \quad \f$ and \f$ \quad \frac{\partial \mathscr{D}_{ij}}{\partial p} \quad \f$
        *   where the binary collision integral \f$ \quad \Omega_{\mathscr{D}, ij} \quad \f$ is given above
        *
        *
        * For developers: please update this info appropriately if you add a new method.
        *
        * \author Valentin N. Zingan, 2013
        * \author M Secanell, 2013-15
        */

        class GasMixture : public BaseMaterial
        {
            public:

                ///@name Constructors, destructor, and initialization
                //@{
                /**
                * Constructor.
                */
                GasMixture();
                /**
                * Constructor.
                */
                GasMixture(const std::string& name);

                /**
                * Destructor.
                */
                virtual ~GasMixture();

                /**
                * Declare parameters.
                */
                virtual void declare_parameters(ParameterHandler& param) const;

                /**
                * Initialize parameters.
                */
                virtual void initialize(ParameterHandler& param);

                /**
                * This function sets
                * \p gases.
                */
                void set_gases(const std::vector< PureGas* >& rgases)
                {
                    gases.clear();
                    gases = rgases;
                }

                /**
                * This function takes a pressure (Pa) and if \p Isobaric \p fluid \p flow is set to \p true 
                * in the data file then the inputted pressure will be used else the \p _DUMMY_ value 
                * is used.
                */
                void set_total_pressure(const double& rtotal_pressure)
                {
                    total_pressure = rtotal_pressure;
                }

                /**
                * This function takes a temperature (K) and if \p Isothermal \p fluid \p flow is set to \p true 
                * in the data file then the inputted temperature will be used else the \p _DUMMY_ value 
                * is used.
                */
                void set_temperature(const double& rtemperature)
                {
                    temperature = rtemperature;
                }

                //@}

                ///@name Accessors and info
                //@{
                /**
                * Function returning the number of gases in the mixture
                */
                inline unsigned int n_gases() const
                {
                    return gases.size();
                }
                
                /**
                 * Return gas stored in index ind
                 */
                inline PureGas* get_gas(unsigned int& ind) const
                {
                    return gases[ind];
                }

                /**
                * This function returns
                * \p gases.
                */
                const std::vector< PureGas* >& get_gases() const
                {
                    return gases;
                }

                /**
                * This function returns
                * \p total_pressure.
                */
                const double& get_total_pressure() const
                {
                    return total_pressure;
                }

                /**
                * This function returns
                * \p temperature.
                */
                const double& get_temperature() const
                {
                    return temperature;
                }

                /**
                * This function prints out
                * the material properties.
                */
                virtual void print_material_properties() const;

                //@}

                ///@name Service functions. Chapman Enskog isobaric diffusion coefficient. Binary gas mixture only.
                //@{

                /**
                * This function returns
                * Maxwell-Stefan isobaric diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
                * written in the Chapman Enskog form (binary gas mixture only) at a constant temperature.
                */
                const double get_ChapmanEnskog_isobaric_diffusion_coefficient() const;

                /**
                * This function returns
                * Maxwell-Stefan isobaric diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
                * written in the Chapman Enskog form (binary gas mixture only) at a constant temperature
                * in the quadrature points of a mesh entity.
                *
                * @param diffusion_coefficient - Chapman Enskog isobaric diffusion coefficient
                *                                at a constant temperature in the quadrature points of a mesh entity.
                */
                void get_ChapmanEnskog_isobaric_diffusion_coefficient(std::vector<double>& diffusion_coefficient) const;

                /**
                * This function returns
                * Maxwell-Stefan isobaric diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
                * written in the Chapman Enskog form (binary gas mixture only) at a variable temperature.
                *
                * @param temperature - temperature.
                */
                const double get_ChapmanEnskog_isobaric_diffusion_coefficient(const double& temperature) const;

                /**
                * This function returns
                * Maxwell-Stefan isobaric diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
                * written in the Chapman Enskog form (binary gas mixture only) at a variable temperature
                * in the quadrature points of a mesh entity.
                *
                * @param temperature           - temperature in the quadrature points of a mesh entity,
                * @param diffusion_coefficient - Chapman Enskog isobaric diffusion coefficient
                *                                at a variable temperature in the quadrature points of a mesh entity.
                */
                void get_ChapmanEnskog_isobaric_diffusion_coefficient(const std::vector<double>& temperature,
                                                                        std::vector<double>&       diffusion_coefficient) const;

                //@}

                ///@name Service functions. Derivatives of Chapman Enskog isobaric diffusion coefficient. Binary gas mixture only.
                //@{

                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial D_{12}}{\partial T} \quad \f$ of the
                * Maxwell-Stefan isobaric diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
                * written in the Chapman Enskog form (binary gas mixture only) at a variable temperature.
                *
                * @param temperature - temperature.
                */
                const double get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature(const double& temperature) const;

                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial D_{12}}{\partial T} \quad \f$ of the
                * Maxwell-Stefan isobaric diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
                * written in the Chapman Enskog form (binary gas mixture only) at a variable temperature
                * in the quadrature points of a mesh entity.
                *
                * @param temperature - temperature in the quadrature points of a mesh entity,
                * @param dst         - \f$ \frac{\partial D_{12}}{\partial T} \quad \f$
                *                      at a variable temperature in the quadrature points of a mesh entity.
                */
                void get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature(const std::vector<double>& temperature,
                                                                                    std::vector<double>&       dst) const;

                //@}

                ///@name Service functions. Chapman Enskog diffusion coefficient. Binary gas mixture only.
                //@{

                /**
                * This function returns
                * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
                * written in the Chapman Enskog form (binary gas mixture only) at a constant total pressure and temperature.
                */
                const double get_ChapmanEnskog_diffusion_coefficient() const;

                /**
                * This function returns
                * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
                * written in the Chapman Enskog form (binary gas mixture only) at a constant total pressure and temperature
                * in the quadrature points of a mesh entity.
                *
                * @param diffusion_coefficient - Chapman Enskog diffusion coefficient
                *                                at a constant total pressure and temperature in the quadrature points of a mesh entity.
                */
                void get_ChapmanEnskog_diffusion_coefficient(std::vector<double>& diffusion_coefficient) const;

                /**
                * This function returns
                * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
                * written in the Chapman Enskog form (binary gas mixture only) at a constant total pressure and variable temperature.
                *
                * @param temperature - temperature.
                */
                const double get_ChapmanEnskog_diffusion_coefficient_at_constant_pressure(const double& temperature) const;

                /**
                * This function returns
                * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
                * written in the Chapman Enskog form (binary gas mixture only) at a constant total pressure and variable temperature
                * in the quadrature points of a mesh entity.
                *
                * @param temperature           - temperature in the quadrature points of a mesh entity,
                * @param diffusion_coefficient - Chapman Enskog diffusion coefficient
                *                                at a constant total pressure and variable temperature in the quadrature points of a mesh entity.
                */
                void get_ChapmanEnskog_diffusion_coefficient_at_constant_pressure(const std::vector<double>& temperature,
                                                                                    std::vector<double>&       diffusion_coefficient) const;

                /**
                * This function returns
                * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
                * written in the Chapman Enskog form (binary gas mixture only) at a variable total pressure and constant temperature.
                *
                * @param total_pressure - total pressure.
                */
                const double get_ChapmanEnskog_diffusion_coefficient_at_constant_temperature(const double& total_pressure) const;

                /**
                * This function returns
                * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
                * written in the Chapman Enskog form (binary gas mixture only) at a variable total pressure and constant temperature
                * in the quadrature points of a mesh entity.
                *
                * @param total_pressure        - total pressure in the quadrature points of a mesh entity,
                * @param diffusion_coefficient - Chapman Enskog diffusion coefficient
                *                                at a variable total pressure and constant temperature in the quadrature points of a mesh entity.
                */
                void get_ChapmanEnskog_diffusion_coefficient_at_constant_temperature(const std::vector<double>& total_pressure,
                                                                                    std::vector<double>&       diffusion_coefficient) const;

                /**
                * This function returns
                * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
                * written in the Chapman Enskog form (binary gas mixture only) at a variable total pressure and temperature.
                *
                * @param total_pressure - total pressure,
                * @param temperature    - temperature.
                */
                const double get_ChapmanEnskog_diffusion_coefficient(const double& total_pressure,
                                                                    const double& temperature) const;
                                                                    
                /**
                * This function returns
                * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
                * written in the Chapman Enskog form (binary gas mixture only) at a variable total pressure and temperature
                * in the quadrature points of a mesh entity.
                *
                * @param total_pressure        - total pressure in the quadrature points of a mesh entity,
                * @param temperature           - temperature in the quadrature points of a mesh entity,
                * @param diffusion_coefficient - Chapman Enskog diffusion coefficient
                *                                at a variable total pressure and temperature in the quadrature points of a mesh entity.
                */
                void get_ChapmanEnskog_diffusion_coefficient(const std::vector<double>& total_pressure,
                                                            const std::vector<double>& temperature,
                                                            std::vector<double>&       diffusion_coefficient) const;

                //@}

                ///@name Service functions. Derivatives of Chapman Enskog diffusion coefficient. Binary gas mixture only.
                //@{

                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{12}}{\partial p} \quad \f$ of the
                * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
                * written in the Chapman Enskog form (binary gas mixture only) at a variable total pressure and constant temperature.
                *
                * @param total_pressure - total pressure.
                */
                const double get_DChapmanEnskog_diffusion_coefficient_Dpressure(const double& total_pressure) const;

                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{12}}{\partial p} \quad \f$ of the
                * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
                * written in the Chapman Enskog form (binary gas mixture only) at a variable total pressure and constant temperature
                * in the quadrature points of a mesh entity.
                *
                * @param total_pressure - total pressure in the quadrature points of a mesh entity,
                * @param dst            - \f$ \frac{\partial \mathscr{D}_{12}}{\partial p} \quad \f$
                *                         at a variable total pressure and constant temperature in the quadrature points of a mesh entity.
                */
                void get_DChapmanEnskog_diffusion_coefficient_Dpressure(const std::vector<double>& total_pressure,
                                                                        std::vector<double>&       dst) const;

                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{12}}{\partial p} \quad \f$ of the
                * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
                * written in the Chapman Enskog form (binary gas mixture only) at a variable total pressure and temperature.
                *
                * @param total_pressure - total pressure,
                * @param temperature    - temperature.
                */
                const double get_DChapmanEnskog_diffusion_coefficient_Dpressure(const double& total_pressure,
                                                                                const double& temperature) const;

                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{12}}{\partial p} \quad \f$ of the
                * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
                * written in the Chapman Enskog form (binary gas mixture only) at a variable total pressure and temperature
                * in the quadrature points of a mesh entity.
                *
                * @param total_pressure - total pressure in the quadrature points of a mesh entity,
                * @param temperature    - temperature in the quadrature points of a mesh entity,
                * @param dst            - \f$ \frac{\partial \mathscr{D}_{12}}{\partial p} \quad \f$
                *                         at a variable total pressure and temperature in the quadrature points of a mesh entity.
                */
                void get_DChapmanEnskog_diffusion_coefficient_Dpressure(const std::vector<double>& total_pressure,
                                                                        const std::vector<double>& temperature,
                                                                        std::vector<double>&       dst) const;

                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{12}}{\partial T} \quad \f$ of the
                * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
                * written in the Chapman Enskog form (binary gas mixture only) at a constant total pressure and variable temperature.
                *
                * @param temperature - temperature.
                */
                const double get_DChapmanEnskog_diffusion_coefficient_Dtemperature(const double& temperature) const;

                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{12}}{\partial T} \quad \f$ of the
                * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
                * written in the Chapman Enskog form (binary gas mixture only) at a constant total pressure and variable temperature
                * in the quadrature points of a mesh entity.
                *
                * @param temperature - temperature in the quadrature points of a mesh entity,
                * @param dst         - \f$ \frac{\partial \mathscr{D}_{12}}{\partial T} \quad \f$
                *                      at a constant total pressure and variable temperature in the quadrature points of a mesh entity.
                */
                void get_DChapmanEnskog_diffusion_coefficient_Dtemperature(const std::vector<double>& temperature,
                                                                            std::vector<double>&       dst) const;

                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{12}}{\partial T} \quad \f$ of the
                * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
                * written in the Chapman Enskog form (binary gas mixture only) at a variable total pressure and temperature.
                *
                * @param total_pressure - total pressure,
                * @param temperature    - temperature.
                */
                const double get_DChapmanEnskog_diffusion_coefficient_Dtemperature(const double& total_pressure,
                                                                                    const double& temperature) const;

                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{12}}{\partial T} \quad \f$ of the
                * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
                * written in the Chapman Enskog form (binary gas mixture only) at a variable total pressure and temperature
                * in the quadrature points of a mesh entity.
                *
                * @param total_pressure - total pressure in the quadrature points of a mesh entity,
                * @param temperature    - temperature in the quadrature points of a mesh entity,
                * @param dst            - \f$ \frac{\partial \mathscr{D}_{12}}{\partial T} \quad \f$
                *                         at a variable total pressure and temperature in the quadrature points of a mesh entity.
                */
                void get_DChapmanEnskog_diffusion_coefficient_Dtemperature(const std::vector<double>& total_pressure,
                                                                            const std::vector<double>& temperature,
                                                                            std::vector<double>&       dst) const;

                //@}

                ///@name Service functions. Chapman Enskog isobaric diffusion coefficients. Ternary and more complicated gas mixtures.
                //@{

                /**
                * This function returns
                * Maxwell-Stefan isobaric diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
                * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a constant temperature.
                */
                const Table< 2, double > get_ChapmanEnskog_isobaric_diffusion_coefficients() const;

                /**
                * This function returns
                * Maxwell-Stefan isobaric diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
                * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a constant temperature
                * in the quadrature points of a mesh entity.
                *
                * @param diffusion_coefficients - Chapman Enskog isobaric diffusion coefficients
                *                                 at a constant temperature in the quadrature points of a mesh entity.
                */
                void get_ChapmanEnskog_isobaric_diffusion_coefficients(std::vector< Table< 2, double > >& diffusion_coefficients) const;

                /**
                * This function returns
                * Maxwell-Stefan isobaric diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
                * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable temperature.
                *
                * @param temperature - temperature.
                */
                const Table< 2, double > get_ChapmanEnskog_isobaric_diffusion_coefficients(const double& temperature) const;

                /**
                * This function returns
                * Maxwell-Stefan isobaric diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
                * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable temperature
                * in the quadrature points of a mesh entity.
                *
                * @param temperature            - temperature in the quadrature points of a mesh entity,
                * @param diffusion_coefficients - Chapman Enskog isobaric diffusion coefficients
                *                                 at a variable temperature in the quadrature points of a mesh entity.
                */
                void get_ChapmanEnskog_isobaric_diffusion_coefficients(const std::vector<double>&         temperature,
                                                                        std::vector< Table< 2, double > >& diffusion_coefficients) const;

                //@}

                ///@name Service functions. Derivatives of Chapman Enskog isobaric diffusion coefficients. Ternary and more complicated gas mixtures.
                //@{

                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial D_{ij}}{\partial T} \quad \f$ of the
                * Maxwell-Stefan isobaric diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
                * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable temperature.
                *
                * @param temperature - temperature.
                */
                const Table< 2, double > get_DChapmanEnskog_isobaric_diffusion_coefficients_Dtemperature(const double& temperature) const;

                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial D_{ij}}{\partial T} \quad \f$ of the
                * Maxwell-Stefan isobaric diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
                * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable temperature
                * in the quadrature points of a mesh entity.
                *
                * @param temperature - temperature in the quadrature points of a mesh entity,
                * @param dst         - \f$ \frac{\partial D_{ij}}{\partial T} \quad \f$
                *                      at a variable temperature in the quadrature points of a mesh entity.
                */
                void get_DChapmanEnskog_isobaric_diffusion_coefficients_Dtemperature(const std::vector<double>&         temperature,
                                                                                    std::vector< Table< 2, double > >& dst) const;

                //@}

                ///@name Service functions. Chapman Enskog diffusion coefficients. Ternary and more complicated gas mixtures.
                //@{

                /**
                * This function returns
                * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
                * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a constant total pressure and temperature.
                */
                const Table< 2, double > get_ChapmanEnskog_diffusion_coefficients() const;

                /**
                * This function returns
                * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
                * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a constant total pressure and temperature
                * in the quadrature points of a mesh entity.
                *
                * @param diffusion_coefficients - Chapman Enskog diffusion coefficients
                *                                 at a constant total pressure and temperature in the quadrature points of a mesh entity.
                */
                void get_ChapmanEnskog_diffusion_coefficients(std::vector< Table< 2, double > >& diffusion_coefficients) const;

                /**
                * This function returns
                * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
                * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a constant total pressure and variable temperature.
                *
                * @param temperature - temperature.
                */
                const Table< 2, double > get_ChapmanEnskog_diffusion_coefficients_at_constant_pressure(const double& temperature) const;

                /**
                * This function returns
                * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
                * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a constant total pressure and variable temperature
                * in the quadrature points of a mesh entity.
                *
                * @param temperature            - temperature in the quadrature points of a mesh entity,
                * @param diffusion_coefficients - Chapman Enskog diffusion coefficients
                *                                 at a constant total pressure and variable temperature in the quadrature points of a mesh entity.
                */
                void get_ChapmanEnskog_diffusion_coefficients_at_constant_pressure(const std::vector<double>&         temperature,
                                                                                    std::vector< Table< 2, double > >& diffusion_coefficients) const;

                /**
                * This function returns
                * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
                * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable total pressure and constant temperature.
                *
                * @param total_pressure - total pressure.
                */
                const Table< 2, double > get_ChapmanEnskog_diffusion_coefficients_at_constant_temperature(const double& total_pressure) const;

                /**
                * This function returns
                * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
                * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable total pressure and constant temperature
                * in the quadrature points of a mesh entity.
                *
                * @param total_pressure         - total pressure in the quadrature points of a mesh entity,
                * @param diffusion_coefficients - Chapman Enskog diffusion coefficients
                *                                 at a variable total pressure and constant temperature in the quadrature points of a mesh entity.
                */
                void get_ChapmanEnskog_diffusion_coefficients_at_constant_temperature(const std::vector<double>&         total_pressure,
                                                                                        std::vector< Table< 2, double > >& diffusion_coefficients) const;

                /**
                * This function returns
                * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
                * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable total pressure and temperature.
                *
                * @param total_pressure - total pressure,
                * @param temperature    - temperature.
                */
                const Table< 2, double > get_ChapmanEnskog_diffusion_coefficients(const double& total_pressure,
                                                                                    const double& temperature) const;

                /**
                * This function returns
                * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
                * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable total pressure and temperature
                * in the quadrature points of a mesh entity.
                *
                * @param total_pressure         - total pressure in the quadrature points of a mesh entity,
                * @param temperature            - temperature in the quadrature points of a mesh entity,
                * @param diffusion_coefficients - Chapman Enskog diffusion coefficients
                *                                 at a variable total pressure and temperature in the quadrature points of a mesh entity.
                */
                void get_ChapmanEnskog_diffusion_coefficients(const std::vector<double>&         total_pressure,
                                                                const std::vector<double>&         temperature,
                                                                std::vector< Table< 2, double > >& diffusion_coefficients) const;

                //@}

                ///@name Service functions. Derivatives of Chapman Enskog diffusion coefficients. Ternary and more complicated gas mixtures.
                //@{

                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{ij}}{\partial p} \quad \f$ of the
                * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
                * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable total pressure and constant temperature.
                *
                * @param total_pressure - total pressure.
                */
                const Table< 2, double > get_DChapmanEnskog_diffusion_coefficients_Dpressure(const double& total_pressure) const;

                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{ij}}{\partial p} \quad \f$ of the
                * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
                * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable total pressure and constant temperature
                * in the quadrature points of a mesh entity.
                *
                * @param total_pressure - total pressure in the quadrature points of a mesh entity,
                * @param dst            - \f$ \frac{\partial \mathscr{D}_{ij}}{\partial p} \quad \f$
                *                         at a variable total pressure and constant temperature in the quadrature points of a mesh entity.
                */
                void get_DChapmanEnskog_diffusion_coefficients_Dpressure(const std::vector<double>&         total_pressure,
                                                                        std::vector< Table< 2, double > >& dst) const;

                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{ij}}{\partial p} \quad \f$ of the
                * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
                * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable total pressure and temperature.
                *
                * @param total_pressure - total pressure,
                * @param temperature    - temperature.
                */
                const Table< 2, double > get_DChapmanEnskog_diffusion_coefficients_Dpressure(const double& total_pressure,
                                                                                            const double& temperature) const;

                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{ij}}{\partial p} \quad \f$ of the
                * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
                * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable total pressure and temperature
                * in the quadrature points of a mesh entity.
                *
                * @param total_pressure - total pressure in the quadrature points of a mesh entity,
                * @param temperature    - temperature in the quadrature points of a mesh entity,
                * @param dst            - \f$ \frac{\partial \mathscr{D}_{ij}}{\partial p} \quad \f$
                *                         at a variable total pressure and temperature in the quadrature points of a mesh entity.
                */
                void get_DChapmanEnskog_diffusion_coefficients_Dpressure(const std::vector<double>&         total_pressure,
                                                                        const std::vector<double>&         temperature,
                                                                        std::vector< Table< 2, double > >& dst) const;

                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{ij}}{\partial T} \quad \f$ of the
                * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
                * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a constant total pressure and variable temperature.
                *
                * @param temperature - temperature.
                */
                const Table< 2, double > get_DChapmanEnskog_diffusion_coefficients_Dtemperature(const double& temperature) const;

                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{ij}}{\partial T} \quad \f$ of the
                * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
                * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a constant total pressure and variable temperature
                * in the quadrature points of a mesh entity.
                *
                * @param temperature - temperature in the quadrature points of a mesh entity,
                * @param dst         - \f$ \frac{\partial \mathscr{D}_{ij}}{\partial T} \quad \f$
                *                      at a constant total pressure and variable temperature in the quadrature points of a mesh entity.
                */
                void get_DChapmanEnskog_diffusion_coefficients_Dtemperature(const std::vector<double>&         temperature,
                                                                            std::vector< Table< 2, double > >& dst) const;

                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{ij}}{\partial T} \quad \f$ of the
                * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
                * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable total pressure and temperature.
                *
                * @param total_pressure - total pressure,
                * @param temperature    - temperature.
                */
                const Table< 2, double > get_DChapmanEnskog_diffusion_coefficients_Dtemperature(const double& total_pressure,
                                                                                                const double& temperature) const;

                /**
                * This function returns
                * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{ij}}{\partial T} \quad \f$ of the
                * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
                * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable total pressure and temperature
                * in the quadrature points of a mesh entity.
                *
                * @param total_pressure - total pressure in the quadrature points of a mesh entity,
                * @param temperature    - temperature in the quadrature points of a mesh entity,
                * @param dst            - \f$ \frac{\partial \mathscr{D}_{ij}}{\partial T} \quad \f$
                *                         at a variable total pressure and temperature in the quadrature points of a mesh entity.
                */
                void get_DChapmanEnskog_diffusion_coefficients_Dtemperature(const std::vector<double>&         total_pressure,
                                                                            const std::vector<double>&         temperature,
                                                                            std::vector< Table< 2, double > >& dst) const;
                //@}
                
                ///@name Service functions. Calculate partial viscosity mixture
                //@{
                
                /**
                 * This function calculates the partial viscosity based on value stored in mixture_viscosity_mode and \b ASSUMES \b ISOTHERMAL and \b NONISOBARIC fluid flow. 
                 * Developers should call this function instead of directly calling Wilke, OmegaKG, or Dynamic directly for the partial viscosity model.
                 * This function requires the following parameters:
                 * 
                 * @param tempOfMixture     - temperature of the element; [K],
                 * @param density           - density in element (i.e. density[s][q], s=species, q=quadrature); [g/cm^3],
                 * @param molarMass         - vector of molar masses of each species; [g/mol],
                 * @param dynamicViscosity  - vector of dynamic viscosities of each species; [g/cm*s],
                 * @param collisionDiameter - vector of collision diameters of each species; [m],
                 * @param porosity          - vector of quadrature points of porosity in element,
                 * @param paramMatrix       - this parameter provides a necessary matrix and is outputed by reference to be used for calculating the variation in
                 *                            partial viscosity. This is done to decrease computations, and used strictly internally,
                 *                            in get_delta_partial_viscosity(). Value is cleared and resized internally. 
                 *                            Vector format is paramMatrix[q][s1][s2], q=quadrature, s1=species1, s2=species2,
                 * @param PInv              - this parameter provides the inverse of P for the calculation of OmegaKG viscosity and variation of.
                 * 
                 * NOTE: you have to set std::vector< NAME::PureGas* > gases in your GasMixture object to be the gases used, before calling this function
                 * 
                 */
                const std::vector< std::vector<double> > get_isothermal_nonisobaric_partial_viscosity(const double&                                      tempOfMixture,
                                                                                                      const std::vector< std::vector<double> >&          density,
                                                                                                      const std::vector<double>&                         molarMass,
                                                                                                      const std::vector<double>&                         dynamicViscosity,
                                                                                                      const std::vector<double>&                         collisionDiameter,
                                                                                                      const std::vector<double>&                         porosity,
                                                                                                      std::vector< std::vector< std::vector<double> > >& paramMatrix,
                                                                                                      std::vector< FullMatrix<double> >&                 PInv) const;
                
                /**
                 * This function calculates the partial viscosity based on value stored in mixture_viscosity_mode and \b ASSUMES \b NONISOTHERMAL and \b NONISOBARIC fluid flow. 
                 * Developers should call this function instead of directly calling Wilke, OmegaKG, or Dynamic directly for the partial viscosity model.
                 * This function requires the following parameters:
                 * 
                 * @param temperature       - temperature in element (i.e. temperature[q], q=quadrature); [K],
                 * @param density           - density in element (i.e. density[s][q], s=species, q=quadrature); [g/cm^3],
                 * @param molarMass         - vector of molar masses of each species; [g/mol],
                 * @param dynamicViscosity  - dynamic viscosities of element (i.e. dynamicViscosity[s][q], s=species, q=quadrature); [g/cm*s],
                 * @param collisionDiameter - vector of collision diameters of each species; [m],
                 * @param porosity          - vector of quadrature points of porosity in element,
                 * @param paramMatrix       - this parameter provides a necessary matrix and is outputed by reference to be used for calculating the variation in
                 *                            partial viscosity. This is done to decrease computations, and used strictly internally,
                 *                            in get_delta_partial_viscosity(). Value is cleared and resized internally. 
                 *                            Vector format is paramMatrix[q][s1][s2], q=quadrature, s1=species1, s2=species2,
                 * @param PInv              - this parameter provides the inverse of P for the calculation of OmegaKG viscosity and variation of.
                 * 
                 * \b NOTE: you have to set std::vector< NAME::PureGas* > gases in your GasMixture object to be the gases used, before calling this function
                 * 
                 */
                const std::vector< std::vector<double> > get_nonisothermal_nonisobaric_partial_viscosity(const std::vector<double>&                         temperature,
                                                                                                         const std::vector< std::vector<double> >&          density,
                                                                                                         const std::vector<double>&                         molarMass,
                                                                                                         const std::vector< std::vector<double> >&          dynamicViscosity,
                                                                                                         const std::vector<double>&                         collisionDiameter,
                                                                                                         const std::vector<double>&                         porosity,
                                                                                                         std::vector< std::vector< std::vector<double> > >& paramMatrix,
                                                                                                         std::vector< FullMatrix<double> >&                 PInv) const;
                
                /**
                 * This is an overloaded function that calculates the Wilke partial viscosity model and \b ASSUMES \b ISOTHERMAL and \b NONISOBARIC fluid flow, 
                 * but doesn't return the \f$ \xi_{ij} \f$ matrix. The Wilke partial viscosity model equations are:
                 * 
                 * - \f$ \eta_i = \frac{\rho_i \eta_i^0}{ \varepsilon M_i \sum_{j=1}^n (\rho_j / M_j) \xi_{ij}} \f$,
                 * 
                 * - \f$ \xi_{ij} = \frac{\left( 1 + (\eta_i^0 / \eta_j^0)^{1/2} (M_j / M_i)^{1/4} \right)^2}{\left( 8 ( 1 + M_i / M_j ) \right)^{1/2}} \f$
                 * 
                 * @param density (\f$ \rho_i \f$)            - density; [g/cm^3],
                 * @param dynamicViscosity (\f$ \eta^0_i \f$) - dynamic viscosity [g/cm*s],
                 * @param molarMass (\f$ M_i \f$)             - molar mass [g/mol],
                 * @param porosixty (\f$ \varepsilon \f$)     - porosity.
                 * 
                 * @htmlonly
                 *   <ul>
                 *       <li> Wilke, C. R. "A viscosity equation for gas mixtures." The journal of chemical physics 18.4 (1950): 517-519. </li>
                 *   </ul>
                 * @endhtmlonly
                 * 
                 */
                std::vector<double> get_isothermal_nonisobaric_Wilke_partial_viscosity(const std::vector<double>& density,
                                                                                       const std::vector<double>& dynamicViscosity,
                                                                                       const std::vector<double>& molarMass,
                                                                                       const double&              porosity) const;
                
                /**
                 * Same as another overloaded function (\b ASSUMES \b ISOTHERMAL and \b NONISOBARIC fluid flow). The difference is that this one will return
                 * by reference \f$ \xi_{ij} \f$. This
                 * is for convience if later calculations are needed so \f$ \xi_{ij} \f$ does not need to be re-computed
                 * (like when calculating the variation in the partial viscosity of the mixture).
                 * No need to size \f$ \xi_{ij} \f$ properly as all values are cleared at the 
                 * beginning and then values pushed in after. The volume averaged Wilke partial viscosity model equations are:
                 * 
                 * - \f$ \eta_i = \frac{\rho_i \eta_i^0}{ \varepsilon M_i \sum_{j=1}^n (\rho_j / M_j) \xi_{ij}} \f$,
                 * 
                 * - \f$ \xi_{ij} = \frac{\left( 1 + (\eta_i^0 / \eta_j^0)^{1/2} (M_j / M_i)^{1/4} \right)^2}{\left( 8 ( 1 + M_i / M_j ) \right)^{1/2}} \f$
                 * 
                 * @param density (\f$ \rho_i \f$)            - phase average of density; [g/cm^3],
                 * @param dynamicViscosity (\f$ \eta^0_i \f$) - dynamic viscosity [g/cm*s],
                 * @param molarMass (\f$ M_i \f$)             - molar mass [g/mol],
                 * @param porosity (\f$ \varepsilon_i \f$)    - porosity,
                 * @param xi (\f$ \xi_{ij} \f$)               - returned by reference for possible later calculations (i.e. xi[s1][s2], s1=species1, s2=species2).
                 * 
                 * @htmlonly
                 *   <ul>
                 *       <li> Wilke, C. R. "A viscosity equation for gas mixtures." The journal of chemical physics 18.4 (1950): 517-519. </li>
                 *   </ul>
                 * @endhtmlonly
                 * 
                 */
                std::vector<double> get_isothermal_nonisobaric_Wilke_partial_viscosity(const std::vector<double>&                density,
                                                                                       const std::vector<double>&                dynamicViscosity,
                                                                                       const std::vector<double>&                molarMass,
                                                                                       const double&                             porosity,
                                                                                             std::vector< std::vector<double> >& xi) const;
                
                /**
                 * This function calculates the partial viscosity using the Wilke model based on value stored in mixture_viscosity_mode and 
                 * \b ASSUMES \b NONISOTHERMAL and \b NONISOBARIC fluid flow. 
                 */
                std::vector<double> get_nonisothermal_nonisobaric_Wilke_partial_viscosity(const std::vector<double>&                density,
                                                                                          const std::vector<double>&                dynamicViscosity,
                                                                                          const std::vector<double>&                molarMass,
                                                                                          const double&                             porosity,
                                                                                                std::vector< std::vector<double> >& xi) const;
                
                /**
                 * This is an overloaded function that calculates the Kerkhof and Geboers partial viscosity model based on the Omega integrals
                 * and \b ASSUMES \b ISOTHERMAL and \b NONISOBARIC fluid flow,
                 * but doesn't return the \f$\boldsymbol{\hat{P}}^{-1}\f$ matrix or a table of the Omega integrals. The equations for the OmegaKG partial viscosity model are:
                 * 
                 * - \f$ \boldsymbol{\hat{P}} \boldsymbol{\eta} = \boldsymbol{1} \f$,
                 * 
                 * - \f$ P_{ii}
                 * = 
                 * \frac{2}{kT} \left(
                 * \frac{4}{5} \Omega_{ii}^{(2,2)} + 
                 * \frac{16 M_i}{15 \rho_i} \sum_{j \ne i}^{n} \frac{\rho_j}{(M_i + M_j)^2} \left( 5 M_i \Omega_{ij}^{(1,1)} + \frac{3}{2} M_j \Omega_{ij}^{(2,2)}  \right) \right) \f$ ,
                 * 
                 * - \f$ P_{ij}
                 * = 
                 * -\frac{2}{kT} \left(
                 * \frac{16}{15} \left( \frac{M_i M_j}{(M_i + M_j)^2} \right) \left( 5 \Omega_{ij}^{(1,1)} - \frac{3}{2} \Omega_{ij}^{(2,2)} \right) \right) \f$ ,
                 * 
                 * - \f$ \Omega_{ij}^{(\ell,s)}
                 * =
                 * \Omega_{ij}^{*(\ell,s)} \left[\Omega_{ij}^{(\ell,s)}\right]_{rs} \f$ ,
                 * 
                 * - \f$ \left[Q_{ij}^{\ell}\right]_{rs} 
                 * =
                 * \left( 1 - \frac{1}{2} \left( \frac{1+(-1)^\ell}{1+\ell}\right) \right)\pi \sigma_{ij}^2 \f$ ,
                 * 
                 * - \f$ \left[\Omega_{ij}^{(\ell,s)}\right]_{rs}
                 * =
                 * \sqrt{\frac{kT}{2 \pi \mu_{ij}}} \left( \frac{(s+1)!}{2} \right) \left[Q_{ij}^{\ell}\right]_{rs} \f$ ,
                 * 
                 * - \f$ \sigma_{ij} 
                 * =
                 * \frac{1}{2} \left( \sigma_i + \sigma_j \right) \f$ ,
                 * 
                 * - \f$ \mu_{ij}
                 * =
                 * \frac{M_i M_j}{N_A (M_i + M_j)} \f$
                 * 
                 * @param density (\f$ \rho_i \f$ )             - density; [g/cm^3],
                 * @param collisionDiameter (\f$ \sigma_i \f$ ) - collision diameter [m],
                 * @param molarMass (\f$ M_i \f$)               - molar mass [g/mol],
                 * @param temperature (\f$ T \f$)               - temperature [K].
                 * @param porosity (\f$ \varepsilon_i \f$)      - porosity,
                 * 
                 * \b NOTE: that the global variable gases is used in this function. So set before functional call by using GasMixture object set_gases() function.
                 * 
                 * @htmlonly
                 *   <ul>
                 *       <li> Kerkhof, Piet JAM, and Marcel AM Geboers. "Toward a unified theory of isotropic molecular transport phenomena." AIChE journal 51.1 (2005): 79-121. </li>
                 *   </ul>
                 * @endhtmlonly
                 * 
                 */
                std::vector<double> get_isothermal_nonisobaric_OmegaKG_partial_viscosity(const std::vector<double>& density,
                                                                                         const std::vector<double>& collisionDiameter,
                                                                                         const std::vector<double>& molarMass,
                                                                                         const double&              temperature,
                                                                                         const double&              porosity) const;
                
                /**
                 * Same as other overloaded function (\b ASSUMES \b ISOTHERMAL and \b NONISOBARIC fluid flow). The difference is that this one will return
                 * by reference \f$ \boldsymbol{\hat{P}}^{-1} \f$ and you can specify a porosity. This
                 * is for convience if later calculations are needed so \f$ \boldsymbol{\hat{P}}^{-1} \f$ does not need to be re-computed
                 * (like when calculating the variation in the partial viscosity of the mixture).
                 * No need to size \f$ \boldsymbol{\hat{P}}^{-1} \f$ properly as all values are cleared at the 
                 * beginning and then values pushed in after. The volume averaged OmegaKG partial viscosity model equations are:
                 * 
                 * - \f$ \boldsymbol{\hat{P}} \boldsymbol{\eta} = \boldsymbol{1} \f$,
                 * 
                 * - \f$ P_{ii}
                 * = 
                 * \frac{2}{kT} \left(
                 * \frac{4}{5} \Omega_{ii}^{(2,2)} + 
                 * \frac{16 M_i}{15 \rho_i} \sum_{j \ne i}^{n} \frac{\rho_j}{(M_i + M_j)^2} \left( 5 M_i \Omega_{ij}^{(1,1)} + \frac{3}{2} M_j \Omega_{ij}^{(2,2)}  \right) \right) \f$ ,
                 * 
                 * - \f$ P_{ij}
                 * = 
                 * -\frac{2}{kT} \left(
                 * \frac{16}{15} \left( \frac{M_i M_j}{(M_i + M_j)^2} \right) \left( 5 \Omega_{ij}^{(1,1)} - \frac{3}{2} \Omega_{ij}^{(2,2)} \right) \right) \f$ ,
                 * 
                 * - \f$ \Omega_{ij}^{(\ell,s)}
                 * =
                 * \Omega_{ij}^{*(\ell,s)} \left[\Omega_{ij}^{(\ell,s)}\right]_{rs} \f$ ,
                 * 
                 * - \f$ \left[Q_{ij}^{\ell}\right]_{rs} 
                 * =
                 * \left( 1 - \frac{1}{2} \left( \frac{1+(-1)^\ell}{1+\ell}\right) \right)\pi \sigma_{ij}^2 \f$ ,
                 * 
                 * - \f$ \left[\Omega_{ij}^{(\ell,s)}\right]_{rs}
                 * =
                 * \sqrt{\frac{kT}{2 \pi \mu_{ij}}} \left( \frac{(s+1)!}{2} \right) \left[Q_{ij}^{\ell}\right]_{rs} \f$ ,
                 * 
                 * - \f$ \sigma_{ij} 
                 * =
                 * \frac{1}{2} \left( \sigma_i + \sigma_j \right) \f$ ,
                 * 
                 * - \f$ \mu_{ij}
                 * =
                 * \frac{M_i M_j}{N_A (M_i + M_j)} \f$
                 * 
                 * @param density (\f$ \rho_i \f$ )                 - phase average of density; [g/cm^3],
                 * @param collisionDiameter (\f$ \sigma_i \f$ )     - collision diameter [m],
                 * @param molarMass (\f$ M_i \f$)                   - molar mass [g/mol],
                 * @param temperature (\f$ T \f$)                   - temperature [K],
                 * @param porosity (\f$ \varepsilon \f$ )           - porosity,
                 * @param omegaIntegralTable                        - returned by reference to reduce computations when calculation variation. Contains a 2n x n matrix (n is species number) with the 
                 *                                                    \f$ \Omega_{ij}^{(1,1)} \f$ on top followed by \f$ \Omega_{ij}^{(2,2)} \f$
                 * @param PInv (\f$ \boldsymbol{\hat{P}}^{-1} \f$ ) - returned by reference for possible later calculations (i.e. PInv[s1][s2], s1=species1, s2=species2).
                 * 
                 * \b NOTE: that the global variable gases is used in this function. So set before functional call by using GasMixture object set_gases() function.
                 * 
                 * @htmlonly
                 *   <ul>
                 *       <li> Kerkhof, Piet JAM, and Marcel AM Geboers. "Toward a unified theory of isotropic molecular transport phenomena." AIChE journal 51.1 (2005): 79-121. </li>
                 *   </ul>
                 * @endhtmlonly
                 * 
                 */
                std::vector<double> get_isothermal_nonisobaric_OmegaKG_partial_viscosity(const std::vector<double>&                density,
                                                                                         const std::vector<double>&                collisionDiameter,
                                                                                         const std::vector<double>&                molarMass,
                                                                                         const double&                             temperature,
                                                                                         const double&                             porosity,
                                                                                               std::vector< std::vector<double> >& omegaIntegralTable,
                                                                                               FullMatrix<double>&                 PInv) const;
                /**
                 * This function calculates the partial viscosity using the OmegaKG model based on value stored in mixture_viscosity_mode and 
                 * \b ASSUMES \b NONISOTHERMAL and \b NONISOBARIC fluid flow. 
                 */
                std::vector<double> get_nonisothermal_nonisobaric_OmegaKG_partial_viscosity(const std::vector<double>&                density,
                                                                                            const std::vector<double>&                collisionDiameter,
                                                                                            const std::vector<double>&                molarMass,
                                                                                            const double&                             temperature,
                                                                                            const double&                             porosity,
                                                                                                  std::vector< std::vector<double> >& omegaIntegralTable,
                                                                                                  FullMatrix<double>&                 PInv) const;
                
                /**
                 * This function calculates the variation in partial viscosity based on value stored in mixture_viscosity_mode and \b ASSUMES \b ISOTHERMAL and \b NONISOBARIC fluid flow. 
                 * Developers should call this function
                 * instead of directly calling Wilke, Omega, or Dynamic directly for the delta partial viscosity model. This function requires the following 
                 * parameters:
                 * 
                 * @param tempOfMixture         - temperature of the element [K],
                 * @param paramMatrix           - this parameter provides a necessary matrix and is outputed by reference from get_partial_viscosity. This is done to decrease computations,
                 * @param PInv (\f$ \boldsymbol{\hat{P}}^{-1} \f$) - calculated previously (numerically) and returned by get_partial_viscosity,
                 * @param porosity              - vector of quadrature points of porosity in element,
                 * @param molarMass             - vector of molar masses of each species; [g/mol],
                 * @param dynamicViscosity      - vector of dynamic viscosities of each species; [g/cm*s],
                 * @param collisionDiameter     - vector of collision diameters of each species; [m],
                 * @param deltaDensity          - variation of density in element (i.e. density[s][q][k], s=species, q=quadrature, k=dof); same unit standard as density,
                 * @param density               - density in element (i.e. density[s][q], s=species, q=quadrature); [g/cm^3],
                 * @param deltaPartialViscosity - variation is partial viscosity returned by reference from function in element (i.e. deltaPartialViscosity[s][q][k], s=species, q=quadrature, k=dof).
                 * 
                 */
                void get_isothermal_nonisobaric_delta_partial_viscosity(const double&                                            tempOfMixture,
                                                                        const std::vector< std::vector< std::vector<double> > >& paramMatrix,
                                                                        const std::vector< FullMatrix<double> >&                 PInv,
                                                                        const std::vector<double>&                               porosity,
                                                                        const std::vector<double>&                               molarMass,
                                                                        const std::vector<double>&                               dynamicViscosity,
                                                                        const std::vector<double>&                               collisionDiameter,
                                                                        const std::vector< std::vector< std::vector<double> > >& deltaDensity,
                                                                        const std::vector< std::vector<double> >&                density,
                                                                        std::vector< std::vector< std::vector<double> > >&       deltaPartialViscosity) const;
                
                /**
                 * This function calculates the variation in partial viscosity based on value stored in mixture_viscosity_mode and \b ASSUMES \b NONISOTHERMAL and \b NONISOBARIC fluid flow. 
                 * Developers should call this function instead of directly calling Wilke, Omega, or Dynamic directly for the delta partial viscosity model. This function requires the following 
                 * parameters:
                 * 
                 * @param temperature           - temperature of the element (i.e. temperature[q], q=quadrature); [K],
                 * @param paramMatrix           - this parameter provides a necessary matrix and is outputed by reference from get_partial_viscosity. This is done to decrease computations,
                 * @param PInv (\f$ \boldsymbol{\hat{P}}^{-1} \f$) - calculated previously (numerically) and returned by get_partial_viscosity,
                 * @param porosity              - vector of quadrature points of porosity in element,
                 * @param molarMass             - vector of molar masses of each species; [g/mol],
                 * @param dynamicViscosity      - vector of dynamic viscosities of each species (i.e. dynamicViscosity[s][q], s=species, q=quadrature); [g/cm*s],
                 * @param collisionDiameter     - vector of collision diameters of each species; [m],
                 * @param deltaDensity          - variation of density in element (i.e. density[s][q][k], s=species, q=quadrature, k=dof); same unit standard as density,
                 * @param density               - density in element (i.e. density[s][q], s=species, q=quadrature); [g/cm^3],
                 * @param deltaPartialViscosity - variation is partial viscosity returned by reference from function in element (i.e. deltaPartialViscosity[s][q][k], s=species, q=quadrature, k=dof).
                 * 
                 */
                void get_nonisothermal_nonisobaric_delta_partial_viscosity(const std::vector<double>&                               temperature,
                                                                           const std::vector< std::vector< std::vector<double> > >& paramMatrix,
                                                                           const std::vector< FullMatrix<double> >&                 PInv,
                                                                           const std::vector<double>&                               porosity,
                                                                           const std::vector<double>&                               molarMass,
                                                                           const std::vector< std::vector<double> >&                dynamicViscosity,
                                                                           const std::vector<double>&                               collisionDiameter,
                                                                           const std::vector< std::vector< std::vector<double> > >& deltaDensity,
                                                                           const std::vector< std::vector<double> >&                density,
                                                                           std::vector< std::vector< std::vector<double> > >&       deltaPartialViscosity) const;
                
                /**
                 * Calculates the variation in Wilke's partial viscosity model w.r.t. density (\f$ \rho_{}^{} \f$). The equation it solves is:
                 * 
                 * - \f$ \delta_{\rho} \eta_i 
                 * = 
                 * \frac{\eta_i^0}{\varepsilon M_i} \left( 
                 * \delta \rho_i \left( \sum_j \frac{\rho_j \xi_{ij}}{M_j} \right)^{-1}
                 * -\frac{\rho_i}{\varepsilon^2} \left( \sum_j \frac{\rho_j \xi_{ij}}{M_j} \right)^{-2} \left( \sum_j \frac{\delta \rho_j \xi_{ij}}{M_j} \right)
                 * \right) \f$
                 * 
                 * @param xi (\f$ \xi_{ij} \f$)                         - calculated previously and returned by get_partial_viscosity,
                 * @param porosity (\f$ \varepsilon \f$)                - porosity,
                 * @param molarMass (\f$ M_{}^{} \f$)                   - vector of molar masses of each species; [g/mol],
                 * @param dynamicViscosity (\f$ \eta^{0}_{} \f$)        - vector of dynamic viscosities of each species; [g/cm*s],
                 * @param deltaDensity (\f$ \delta \rho_i \f$)          - varation in phase average of density; same unit standard as density,
                 * @param density (\f$ \rho_i \f$ )                     - phase average of density; [g/cm^3],
                 * @param deltaPartialViscosity (\f$ \delta \eta_i \f$) - variation is partial viscosity returned by reference from function in element (i.e. deltaPartialViscosity[s][q][k], s=species, q=quadrature, k=dof).
                 * 
                 * @htmlonly
                 *   <ul>
                 *       <li> Wilke, C. R. "A viscosity equation for gas mixtures." The journal of chemical physics 18.4 (1950): 517-519. </li>
                 *   </ul>
                 * @endhtmlonly
                 * 
                 */
                void get_Wilke_delta_partial_viscosity_wrt_density(const std::vector< std::vector< std::vector<double> > >& xi,
                                                                   const std::vector<double>&                               porosity,
                                                                   const std::vector<double>&                               molarMass,
                                                                   const std::vector<double>&                               dynamicViscosity,
                                                                   const std::vector< std::vector< std::vector<double> > >& deltaDensity,
                                                                   const std::vector< std::vector<double> >&                density,
                                                                   std::vector< std::vector< std::vector<double> > >&       deltaPartialViscosity) const;
                
                /**
                 * Calculates the variation in Wilke's partial viscosity model w.r.t. temperature (\f$ T_{}^{} \f$). The equation it solves is:
                 * 
                 * - \f$ \delta_{T} \eta_i 
                 * = 
                 * \frac{\rho_i}{M_i} \left( 
                 * \delta_{T} \eta_i^0 \left( \sum_j \frac{\rho_j \xi_{ij}}{M_j} \right)^{-1}
                 * -\eta_i^0 \left( \sum_j \frac{\rho_j \xi_{ij}}{M_j} \right)^{-2} \left( \sum_j \frac{\rho_j \delta_{T} \xi_{ij}}{M_j} \right)
                 * \right) \f$
                 * 
                 * - \f$ \delta_{T} \xi_{ij}
                 * =
                 * \left( 8(1 + \frac{M_i}{M_j}) \right)^{1/2} 
                 * \left( 1 + \left( \frac{\eta_i^0}{\eta_j^0} \right)^{1/2} \left( \frac{M_j}{M_i} \right)^{1/4} \right)
                 * \left( \frac{M_j}{M_i} \right)^{1/4} \left( \delta_T \eta_i^0 \left( \frac{1}{\eta_i^0 \eta_j^0} \right)^{1/2} - \delta_T \eta_j^0 \left( \frac{\eta_i^0}{(\eta_j^0)^3} \right)^{1/2} \right)
                 * \f$
                 * 
                 * \b NOTE: these equations are not volume averaged, which will need to be done.
                 * 
                 * @param xi (\f$ \xi_{ij} \f$)                         - calculated previously and returned by get_partial_viscosity,
                 * @param porosity (\f$ \varepsilon \f$)                - porosity,
                 * @param molarMass (\f$ M_{}^{} \f$)                   - vector of molar masses of each species; [g/mol],
                 * @param dynamicViscosity (\f$ \eta^{0}_{} \f$)        - vector of dynamic viscosities of each species; [g/cm*s],
                 * @param deltaTemperature (\f$ \delta_{T} T \f$)              - varation in phase average of temperature; same unit standard as temperature,
                 * @param density (\f$ \rho_i \f$ )                     - phase average of density; [g/cm^3],
                 * @param deltaPartialViscosity (\f$ \delta \eta_i \f$) - variation is partial viscosity returned by reference from function in element (i.e. deltaPartialViscosity[s][q][k], 
                 *                                                        s=species, q=quadrature, k=dof).
                 * 
                 * @htmlonly
                 *   <ul>
                 *       <li> Wilke, C. R. "A viscosity equation for gas mixtures." The journal of chemical physics 18.4 (1950): 517-519. </li>
                 *   </ul>
                 * @endhtmlonly
                 * 
                 */
                void get_Wilke_delta_partial_viscosity_wrt_temperature(const std::vector< std::vector< std::vector<double> > >& xi,
                                                                       const std::vector<double>&                               porosity,
                                                                       const std::vector<double>&                               molarMass,
                                                                       const std::vector<double>&                               dynamicViscosity,
                                                                       const std::vector< std::vector< std::vector<double> > >& deltaTemperature,
                                                                       const std::vector< std::vector<double> >&                density,
                                                                       std::vector< std::vector< std::vector<double> > >&       deltaPartialViscosity) const;
                
                /**
                 * Calculates the variation in Kerkhof and Geboers partial viscosity model using Omega integrals w.r.t. density (\f$ \rho_{}^{} \f$). The equations are:
                 * 
                 * - \f$ \delta \boldsymbol{\eta} = -\hat{\boldsymbol{P}}^{-1} \delta\hat{\boldsymbol{P}} \hat{\boldsymbol{P}}^{-1} \boldsymbol{1} \f$,
                 * 
                 * - \f$ \delta P_{ii}
                 * =
                 * \frac{32 M_i}{15kT \varepsilon} \left( 
                 * \frac{1}{\rho_i} \sum_{j \ne i}^{n} \frac{\delta \rho_j}{(M_i + M_j)^2} \left( 5 M_i \Omega_{ij}^{(1,1)} + \frac{3}{2} M_j \Omega_{ij}^{(2,2)}
                 * \right) 
                 * -\frac{\delta \rho_i}{\varepsilon^2 \rho_i^2} \sum_{j \ne i}^{n} \frac{\rho_j}{(M_i + M_j)^2} \left( 5 M_i \Omega_{ij}^{(1,1)} + \frac{3}{2} M_j \Omega_{ij}^{(2,2)}  \right) \right) \f$,
                 * 
                 * - \f$ \delta P_{ij}^{} = 0 \f$.
                 * 
                 * @param tempOfMixture (\f$ T \f$)                     - temperature of the element [K],
                 * @param omegaIntegralTable                            - returned by reference to reduce computations when calculation variation. Contains a 2n x n matrix (n is species number) with the 
                 *                                                        \f$ \Omega_{ij}^{(1,1)} \f$ on top followed by \f$ \Omega_{ij}^{(2,2)} \f$; (i.e. omegaIntegralTable[q][s1][s2], q=quadrature, 
                 *                                                        s1=species1, s2=species2).
                 * @param PInv (\f$ \boldsymbol{\hat{P}}^{-1} \f$)      - calculated previously (numerically) and returned by get_partial_viscosity,
                 * @param porosity (\f$ \varepsilon \f$)                - porosity,
                 * @param molarMass (\f$ M_{}^{} \f$)                   - vector of molar masses of each species; [g/mol],
                 * @param collisionDiameter (\f$ \sigma^{ij}_{} \f$)    - vector of collision diameters of each species; [m],
                 * @param deltaDensity (\f$ \delta \rho_i \f$)          - varation in phase average of density; same unit standard as density,
                 * @param density (\f$ \rho_i \f$ )                     - phase average of density; [g/cm^3],
                 * @param deltaPartialViscosity (\f$ \delta \eta_i \f$) - variation is partial viscosity returned by reference from function in element (i.e. deltaPartialViscosity[s][q][k], s=species, q=quadrature, k=dof).
                 * 
                 * @htmlonly
                 *   <ul>
                 *       <li> Kerkhof, Piet JAM, and Marcel AM Geboers. "Toward a unified theory of isotropic molecular transport phenomena." AIChE journal 51.1 (2005): 79-121. </li>
                 *   </ul>
                 * @endhtmlonly
                 * 
                 */
                void get_OmegaKG_delta_partial_viscosity_wrt_density(const double&                                            tempOfMixture,
                                                                     const std::vector< std::vector< std::vector<double> > >& omegaIntegralTable,
                                                                     const std::vector< FullMatrix<double> >&                 PInv,
                                                                     const std::vector<double>&                               porosity,
                                                                     const std::vector<double>&                               molarMass,
                                                                     const std::vector<double>&                               collisionDiameter,
                                                                     const std::vector< std::vector< std::vector<double> > >& deltaDensity,
                                                                     const std::vector< std::vector<double> >&                density,
                                                                     std::vector< std::vector< std::vector<double> > >&       deltaPartialViscosity) const;
                
                /**
                 * Calculates the variation in Kerkhof and Geboers partial viscosity model using Omega integrals w.r.t. temperature (\f$ T_{}^{} \f$). The equations are:
                 * 
                 * - Need to be determined.
                 * 
                 * \b NOTE: these equations are not volume averaged, which will need to be done.
                 * 
                 * @param temperature (\f$ T \f$)                       - temperature of the element [K],
                 * @param omegaIntegralTable                            - returned by reference to reduce computations when calculation variation. Contains a 2n x n matrix (n is species number) with the 
                 *                                                        \f$ \Omega_{ij}^{(1,1)} \f$ on top followed by \f$ \Omega_{ij}^{(2,2)} \f$; (i.e. omegaIntegralTable[q][s1][s2], q=quadrature, 
                 *                                                        s1=species1, s2=species2).
                 * @param PInv (\f$ \boldsymbol{\hat{P}}^{-1} \f$)      - calculated previously (numerically) and returned by get_partial_viscosity,
                 * @param porosity (\f$ \varepsilon \f$)                - porosity,
                 * @param molarMass (\f$ M_{}^{} \f$)                   - vector of molar masses of each species; [g/mol],
                 * @param collisionDiameter (\f$ \sigma^{ij}_{} \f$)    - vector of collision diameters of each species; [m],
                 * @param deltaDensity (\f$ \delta \rho_i \f$)          - varation in phase average of density; same unit standard as density,
                 * @param density (\f$ \rho_i \f$ )                     - phase average of density; [g/cm^3],
                 * @param deltaPartialViscosity (\f$ \delta \eta_i \f$) - variation is partial viscosity returned by reference from function in element (i.e. deltaPartialViscosity[s][q][k], s=species, q=quadrature, k=dof).
                 * 
                 * @htmlonly
                 *   <ul>
                 *       <li> Kerkhof, Piet JAM, and Marcel AM Geboers. "Toward a unified theory of isotropic molecular transport phenomena." AIChE journal 51.1 (2005): 79-121. </li>
                 *   </ul>
                 * @endhtmlonly
                 * 
                 */
                void get_OmegaKG_delta_partial_viscosity_wrt_temperature(const std::vector<double>&                               temperature,
                                                                         const std::vector< std::vector< std::vector<double> > >& omegaIntegralTable,
                                                                         const std::vector< FullMatrix<double> >&                 PInv,
                                                                         const std::vector<double>&                               porosity,
                                                                         const std::vector<double>&                               molarMass,
                                                                         const std::vector<double>&                               collisionDiameter,
                                                                         const std::vector< std::vector< std::vector<double> > >& deltaDensity,
                                                                         const std::vector< std::vector<double> >&                density,
                                                                         std::vector< std::vector< std::vector<double> > >&       deltaPartialViscosity) const;
                
                /**
                 * Returns a std::vector< std::vector< std::vector<double> > > of 0s. This is used for the dynamic viscosity
                 * partial viscosity model as it is not a function of density or velocity. As such, the variation is 0.
                 * Note that deltaViscosity must be properly sized before sending to this function.
                 */
                void get_Null_delta_viscosity(std::vector< std::vector< std::vector<double> > >& deltaViscosity) const;
                
                /**
                 * Loops through all variations of partial viscosity and calculates the varation in bulk viscosity according to
                 * bulk_viscosity_mode. All vectors must be sized correctly before using this function.
                 */
                void get_delta_bulk_viscosity(const std::vector< PureGas* >&                           gases,
                                              const std::vector< std::vector< std::vector<double> > >& deltaPartialViscosity,
                                                    std::vector< std::vector< std::vector<double> > >& deltaBulkViscosity) const;
                
                //@}

            protected:
                
                ///@name Service functions. Binary collision integral.
                //@{

                /**
                * This function returns binary collision integral at a constant temperature.
                *
                * @param N1 - number of the first  gas from \p gases,
                * @param N2 - number of the second gas from \p gases.
                * 
                * @htmlonly
                *   <ul>
                *       <li> Neufeld, Philip D., A. R. Janzen, and R. A. Aziz. "Empirical Equations to Calculate 16 of the Transport Collision Integrals Ω (l, s)* for the Lennard‐Jones (12–6) Potential." The Journal of Chemical Physics 57.3 (1972): 1100-1102. </li>
                *   </ul>
                * @endhtmlonly
                */
                const double get_binary_collision_integral(const unsigned int& N1 = 0,
                                                            const unsigned int& N2 = 1) const;

                /**
                * This function returns binary collision integral
                * at a constant temperature in the quadrature points of a mesh entity.
                *
                * @param binary_collision_integral - binary collision integral at a constant temperature
                *                                    in the quadrature points of a mesh entity,
                * @param N1                        - number of the first  gas from \p gases,
                * @param N2                        - number of the second gas from \p gases.
                * 
                * @htmlonly
                *   <ul>
                *       <li> Neufeld, Philip D., A. R. Janzen, and R. A. Aziz. "Empirical Equations to Calculate 16 of the Transport Collision Integrals Ω (l, s)* for the Lennard‐Jones (12–6) Potential." The Journal of Chemical Physics 57.3 (1972): 1100-1102. </li>
                *   </ul>
                * @endhtmlonly
                */
                void get_binary_collision_integral(std::vector<double>& binary_collision_integral,
                                                    const unsigned int&  N1 = 0,
                                                    const unsigned int&  N2 = 1) const;

                /**
                * This function returns binary collision integral at a variable temperature.
                *
                * @param temperature - temperature,
                * @param N1          - number of the first  gas from \p gases,
                * @param N2          - number of the second gas from \p gases.
                * 
                * @htmlonly
                *   <ul>
                *       <li> Neufeld, Philip D., A. R. Janzen, and R. A. Aziz. "Empirical Equations to Calculate 16 of the Transport Collision Integrals Ω (l, s)* for the Lennard‐Jones (12–6) Potential." The Journal of Chemical Physics 57.3 (1972): 1100-1102. </li>
                *   </ul>
                * @endhtmlonly
                */
                const double get_binary_collision_integral(const double&       temperature,
                                                            const unsigned int& N1 = 0,
                                                            const unsigned int& N2 = 1) const;

                /**
                * This function returns binary collision integral
                * at a variable temperature in the quadrature points of a mesh entity.
                *
                * @param temperature               - temperature in the quadrature points of a mesh entity,
                * @param binary_collision_integral - binary collision integral at a variable temperature
                *                                    in the quadrature points of a mesh entity,
                * @param N1                        - number of the first  gas from \p gases,
                * @param N2                        - number of the second gas from \p gases.
                * 
                * @htmlonly
                *   <ul>
                *       <li> Neufeld, Philip D., A. R. Janzen, and R. A. Aziz. "Empirical Equations to Calculate 16 of the Transport Collision Integrals Ω (l, s)* for the Lennard‐Jones (12–6) Potential." The Journal of Chemical Physics 57.3 (1972): 1100-1102. </li>
                *   </ul>
                * @endhtmlonly
                */
                void get_binary_collision_integral(const std::vector<double>& temperature,
                                                    std::vector<double>&       binary_collision_integral,
                                                    const unsigned int&        N1 = 0,
                                                    const unsigned int&        N2 = 1) const;

                /**
                * This function returns the first derivative 
                * \f$ \quad \frac{\partial \Omega_{\mathscr{D}, ij}}{\partial T} \quad \f$ of the
                * binary collision integral at a variable temperature.
                *
                * @param temperature - temperature,
                * @param N1          - number of the first  gas from \p gases,
                * @param N2          - number of the second gas from \p gases.
                */
                const double get_Dbinary_collision_integral_Dtemperature(const double&       temperature,
                                                                        const unsigned int& N1 = 0,
                                                                        const unsigned int& N2 = 1) const;

                /**
                * This function returns the first derivative 
                * \f$ \quad \frac{\partial \Omega_{\mathscr{D}, ij}}{\partial T} \quad \f$ of the
                * binary collision integral at a variable temperature in the quadrature points of a mesh entity.
                *
                * @param temperature - temperature in the quadrature points of a mesh entity,
                * @param dst         - \f$ \frac{\partial \Omega_{\mathscr{D}, ij}}{\partial T} \quad \f$
                *                      at a variable temperature in the quadrature points of a mesh entity,
                * @param N1          - number of the first  gas from \p gases,
                * @param N2          - number of the second gas from \p gases.
                */
                void get_Dbinary_collision_integral_Dtemperature(const std::vector<double>& temperature,
                                                                std::vector<double>&       dst,
                                                                const unsigned int&        N1 = 0,
                                                                const unsigned int&        N2 = 1) const;

                //@}

                ///@name Service functions. Omega* integrals.
                //@{
                
                /**
                * This function returns \f$ \Omega_{}^{*(1,1)} \f$ integral at a variable temperature.
                * This equation is only valid in range \f$ 0.3  \leq T^* \leq 400 \f$ and empirical constants
                * were taken from "Toward a Unified Theory of Isotropic Molecular Transport Phenomena" 
                * by Kerkhof and Geboers, which used Table C2 from the data in "Molecular Theory of Gases and Liquids"
                * by Hirschfelder.
                *
                * @param temperature - temperature [K],
                * @param N1          - number of the first  gas from \p gases,
                * @param N2          - number of the second gas from \p gases.
                * 
                * @htmlonly
                *   <ul>
                *       <li> Kerkhof, Piet JAM, and Marcel AM Geboers. "Toward a unified theory of isotropic molecular transport phenomena." AIChE journal 51.1 (2005): 79-121. </li>
                *       <li> Hirschfelder, Joseph O., et al. Molecular theory of gases and liquids. New York: Wiley, 1964. </li>
                *   </ul>
                * @endhtmlonly
                * 
                */
                const double get_omega_star_11_integral(const double&       temperature,
                                                        const unsigned int& N1 = 0,
                                                        const unsigned int& N2 = 1) const;
                
                /**
                * This function returns \f$ \Omega_{}^{*(1,1)} \f$ integral
                * at a variable temperature in the quadrature points of a mesh entity.
                * this equation is only valid in range \f$ 0.3  \leq T^* \leq 400 \f$ and empirical constants
                * were taken from "Toward a Unified Theory of Isotropic Molecular Transport Phenomena" 
                * by Kerkhof and Geboers, which used Table C2 from the data in "Molecular Theory of Gases and Liquids"
                * by Hirschfelder.
                * 
                * @param temperature    - temperature in the quadrature points of a mesh entity [K],
                * @param omega_integral - \f$ \Omega_{}^{*(1,1)} \f$ integral at a variable temperature in the quadrature points of a mesh entity,
                * @param N1             - index of the first  gas from \p gases,
                * @param N2             - index of the second gas from \p gases.
                * 
                * @htmlonly
                *   <ul>
                *       <li> Kerkhof, Piet JAM, and Marcel AM Geboers. "Toward a unified theory of isotropic molecular transport phenomena." AIChE journal 51.1 (2005): 79-121. </li>
                *       <li> Hirschfelder, Joseph O., et al. Molecular theory of gases and liquids. New York: Wiley, 1964. </li>
                *   </ul>
                * @endhtmlonly
                * 
                */
                void get_omega_star_11_integral(const double &       temperature,
                                                std::vector<double>& omega_integral,
                                                const unsigned int&  N1 = 0,
                                                const unsigned int&  N2 = 1) const;
                
                /**
                * This function returns \f$ \Omega_{}^{*(2,2)} \f$ integral at a variable temperature.
                * this equation is only valid in range \f$ 0.3  \leq T^* \leq 400 \f$ and empirical constants
                * were taken from "Toward a Unified Theory of Isotropic Molecular Transport Phenomena" 
                * by Kerkhof and Geboers, which used Table C2 from the data in "Molecular Theory of Gases and Liquids"
                * by Hirschfelder.
                *
                * @param temperature - temperature [K],
                * @param N1          - number of the first  gas from \p gases,
                * @param N2          - number of the second gas from \p gases.
                * 
                * @htmlonly
                *   <ul>
                *       <li> Kerkhof, Piet JAM, and Marcel AM Geboers. "Toward a unified theory of isotropic molecular transport phenomena." AIChE journal 51.1 (2005): 79-121. </li>
                *       <li> Hirschfelder, Joseph O., et al. Molecular theory of gases and liquids. New York: Wiley, 1964. </li>
                *   </ul>
                * @endhtmlonly
                * 
                */
                const double get_omega_star_22_integral(const double&       temperature,
                                                        const unsigned int& N1 = 0,
                                                        const unsigned int& N2 = 1) const;
                
                /**
                * This function returns \f$ \Omega_{}^{*(2,2)} \f$ integral
                * at a variable temperature in the quadrature points of a mesh entity.
                * this equation is only valid in range \f$ 0.3  \leq T^* \leq 400 \f$ and empirical constants
                * were taken from "Toward a Unified Theory of Isotropic Molecular Transport Phenomena" 
                * by Kerkhof and Geboers, which used Table C2 from the data in "Molecular Theory of Gases and Liquids"
                * by Hirschfelder.
                * 
                * @param temperature    - temperature in the quadrature points of a mesh entity [K],
                * @param omega_integral - \f$ \Omega_{}^{*(2,2)} \f$ integral at a variable temperature in the quadrature points of a mesh entity,
                * @param N1             - index of the first  gas from \p gases,
                * @param N2             - index of the second gas from \p gases.
                * 
                * @htmlonly
                *   <ul>
                *       <li> Kerkhof, Piet JAM, and Marcel AM Geboers. "Toward a unified theory of isotropic molecular transport phenomena." AIChE journal 51.1 (2005): 79-121. </li>
                *       <li> Hirschfelder, Joseph O., et al. Molecular theory of gases and liquids. New York: Wiley, 1964. </li>
                *   </ul>
                * @endhtmlonly
                * 
                */
                void get_omega_star_22_integral(const double&        temperature,
                                                std::vector<double>& omega_integral,
                                                const unsigned int&  N1 = 0,
                                                const unsigned int&  N2 = 1) const;
                
                //@}
                
                //////////
                // DATA //
                //////////
                    
                ///@name Fluid properties
                //@{
                /**
                * This \p std::vector contains all pure gases which form the whole gas mixture of a problem
                * at hand.
                */
                std::vector< PureGas* > gases;
                
                /**
                * This \p bool tells GasMisture if it should use isothermal assumption.
                */
                bool tempIsoTherm  = false;
                
                /**
                 * This \p bool tells GasMisture if it should use isobaric assumption.
                */
                bool pressIsoBaric = false;

                /**
                * Total pressure of the whole mixture, 
                * \f$ p_{\text{total}} \quad \left[ \text{Pa} \right] \f$.
                */
                double total_pressure = 1.0;

                /**
                * Temperature of the whole
                * gas mixture, \f$ T \quad \left[ \text{K} \right] \f$.
                */
                double temperature = 1.0;

                //@}
                
                ///@name Names and modes
                //@{
                
                /**
                * @param mixture_viscosity_mode=Wilke   - Wilke's formula for partial viscosity of each species is used,
                * @param mixture_viscosity_mode=OmegaKG - Kerkhof and Geboers formula for partial viscosity of each species, using the Omega integrals, is used,
                * @param mixture_viscosity_mode=Dynamic - The dynamic viscosity of each species is used as calculated by dynamic_viscosity_mode.
                */
                std::string mixture_viscosity_mode;
                
                //@} 

        };

    } // Material

} // FuelCellShop

#endif