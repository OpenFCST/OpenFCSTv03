//---------------------------------------------------------------------------
//    $Id: operating_conditions.h 2605 2014-08-15 03:36:44Z secanell $
//
//    Copyright (C) 2009 by Marc Secanell
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef _FUELCELL_OPERATING_CONDITIONS__H
#define _FUELCELL_OPERATING_CONDITIONS__H

// Include deal.II classes
#include <deal.II/base/parameter_handler.h>


#include <boost/lexical_cast.hpp>

#include <utils/logging.h>
#include <utils/fcst_constants.h>
#include <application_core/initial_and_boundary_data.h>
#include <grid/geometry.h>
#include <materials/PureGas.h>
#include <materials/GasMixture.h>
#include <utils/fcst_units.h>

using namespace dealii;

namespace FuelCell
{
    
    /**
     * Class used to store, read from file and define the operating conditions for a fuel cell.
     * It stores the following:
     * - cell temperature
     * - cell voltage
     * - oxygen mole fraction at the inlet of the channel
     * - relative humidity at anode and cathode
     * - total pressure at anode and cathode
     * 
     * This information, in conjunction with the water saturation equation is used to compute the
     * mole fractions for each component in the mixture.
     * 
     * If the inlet is humidified air, please set the oxygen mole fraction to 0.21 (default).
     * 
     * In the input file, the following parameters can be specified (see declare_parameters ):
     * @code
     * subsection Fuel cell data
     * (...)
     *   subsection Operating conditions
     *     set Adjust initial solution and boundary conditions = false
     *     set Temperature cell [K] = 353
     *     set Cathode pressure [Pa] = 101325
     *     set Cathode initial oxygen mole fraction (prior to humidification) = 0.21
     *     set Cathode relative humidity = 0.7
     *     set Anode pressure [Pa] = 101325
     *     set Anode relative humidity = 0.7
     *     set Voltage cell [V] = 0.6
     *     set Voltage drop in the anode [V] = 0.015
     *     set Open circuit voltage [V] = 1.23
     *   end
     * end
     * subsection Application
     *   subsection Electrode KG
     *     set Anode simulation = false
     *   end
     * end
     * @endcode 
     * 
     * <h3>Usage Details:</h3>
     * 
     * In order to use this class, first an object of the class needs to be created. Usually, one such objects exists in every application. To create the object,
     * include the .h file in the include application file and in the application data member section add the object. For example:
     * @code
     * #include "operating_conditions.h"
     * 
     * // Then in the data member declaration (usually a private member)
     * FuelCell::OperatingConditions OC;
     * @endcode
     * 
     * Once the object is created, the section where the input data will be specified in the input file needs to be declared. To do so, in the declare_parameters section 
     * of your application call the following:
     * 
     * @code
     * //--------- IN DECLARE_PARAMETERS ------------------------------------------------------
     * template <int dim>
     * void 
     * NAME::AppCathode<dim>::declare_parameters(ParameterHandler& param)
     * {
     *   (...)
     *   OC.declare_parameters(param);
     *   (...)
     * }
     * @endcode
     *          
     * 
     * Finally, once the input file has been read by our application, your class needs to be initialized. This is achieved using the function initialize()
     * @code
     * //--------- IN INITIALIZE ------------------------------------------------------
     * template <int dim>
     * void
     * NAME::AppCathode<dim>::_initialize(ParameterHandler& param)
     * {   
     *  (...) 
     *  OC.initialize(param);
     * }
     * @endcode
     * 
     * You are now ready to use your OperatingConditions object!.
     * 
     * 
     * @author M. Secanell, 2009-2013
     * 
     */
    class OperatingConditions
    {
    public:
        /**
         * Constructor
         */
        OperatingConditions();
        
        /**
         * Destructor
         */
        ~OperatingConditions();
        
        /**
         * Declare all necessary parameters in order to compute the coefficients
         * 
         * The parameters that can be specified in the input file are as follows:
         * 
         * @code
         * subsection Fuel cell data
         * (...)
         *   subsection Operating conditions
         *     set Adjust initial solution and boundary conditions = false  
         *     set Temperature cell = 353
         *     set Cathode pressure = 101325
         *     set Cathode initial oxygen mole fraction (prior to humidification) = 0.21
         *     set Cathode relative humidity = 0.7
         *     set Anode pressure = 101325
         *     set Anode relative humidity = 0.7
         *     set Voltage cell = 0.6
         *     set Voltage drop in the anode = 0.015
         *     set Open circuit voltage = 1.23
         *   end
         * end
         * @endcode 
         */
        void declare_parameters (ParameterHandler &param) const;

        /**
         * Class used to read in data and initialize the necessary data
         * to compute the coefficients.
         */
        void initialize (ParameterHandler& param);
        
        /**
         * Initialize the temperature and pressure for anode gas mixture and store a copy of GasMixture class
         * inside operating conditions.
         * Note: The gas mixture species should already be specified.
         */
        void set_gas_mixtures(FuelCellShop::Material::GasMixture& c_mix, FuelCellShop::Material::GasMixture& a_mix)
        {
            c_mix.set_total_pressure(get_pc_Pa());
            c_mix.set_temperature(get_T());
            
            cathode_mix = &c_mix;
            
            a_mix.set_total_pressure(get_pa_Pa());
            a_mix.set_temperature(get_T());
            
            anode_mix = &a_mix;
        }
            
        /**
         * Compute the density of each species in the anode based on GasMixture specified in set_gas_mixture_anode().
         */
        std::vector<double> get_rho_anode() const;
        
        /**
         * Compute the density of each species in the anode based on GasMixture specified in set_gas_mixture_cathode().
         */
        std::vector<double> get_rho_cathode() const;
        
        
        /**
         * Get the water vapour saturation pressure in atmospheres (atm) using cell temperature. 
         */
        double saturation_pressure() const;
        
        /**
         * Get the mole fraction of water at the selected electrode channel from pressure, temperature, and relative humidity
         */
        double get_x_wv(int) const;
        
        /**
         * Get the mole fraction of water at cathode channel from pressure, temperature, and relative humidity
         */
        double get_x_wv() const;

        /**
         * Get the mole fraction of water at anode channel from pressure, temperature, and relative humidity
         */
        double get_x_wv_anode() const;        
        
        /**
         * Get the mole fraction of oxygen at the cathode channel from pressure, temperature, and relative humidity
         */
        double get_x_o2() const;
        
        /**
         * Get the mole fraction of hydrogen at the anode channel from pressure, temperature, and relative humidity
         */       
        double get_x_h2() const;
        
        /**
         * Get the mole fraction of nitrogen at the anode channel from pressure, temperature, and relative humidity
         */       
        double get_x_n2(int) const;        
        
        /**
         * Get the density of water vapour (wv) at the selected electrode channel from pressure, temperature, mole fraction, and total concentration
         */  
        double get_rho_wv(int) const;
        
        /**
         * Get the density of water vapour (wv) at the cathode channel from pressure, temperature, mole fraction, and total concentration
         */  
        double get_rho_wv() const;

        /**
         * Get the density of water vapour (wv) at the anode channel from pressure, temperature, mole fraction, and total concentration
         */  
        double get_rho_wv_anode() const;        
        
        /**
         * Get the density of oxygen (O2) at the cathode channel from pressure, temperature, mole fraction, and total concentration
         */ 
        double get_rho_o2() const;
        /**
         * Get the density of hydrogen (O2) at the anode channel from pressure, temperature, mole fraction, and total concentration
         */ 
        double get_rho_h2() const;        
        /**
         * Get the density of nitrogen (N2) depending on which electrode we are solving
         */ 
        double get_rho_n2(int) const;
        /**
         * Get the density of nitrogen (N2) at the cathode channel from pressure, temperature, mole fraction, and total concentration
         */ 
        double get_rho_n2() const;
        /**
         * Get the density of nitrogen (N2) at the anode channel from pressure, temperature, mole fraction, and total concentration
         */ 
        double get_rho_n2_anode() const;        
        
        /**
         * This function is meant to be called by an application to give Operating Conditions the map relating species number (key) to 
         * gas species (value). The map is set to the private map gasSpeciesMap. At this time this is only needed by appCathodeKG as it 
         * uses density for fluid flow
         */ 
        void set_gas_map(std::map<unsigned int, std::string> tmp);
        
        /**
         * Get the species map relating species number (key) to species material (value).
         */ 
        std::map<unsigned int, std::string> get_gas_map() const;
        
        /**
         * Submit a string of an unsigned integer value and it will determine, what the respective species material is from the 
         * gasSpeciesMap and return the density at Operating Conditions. Currently, it can only return density for oxygen, 
         * water, and nitrogen.
         */ 
        double get_density_from_map(std::string species) const;
        
        /**
         * NOTE: This function is redefined in base_kinetics class, considering variable temperature and gas pressures. The function here should only be used for defining Initial condition values.
         * Get the theoretical cell voltage using the cell temperature, and the reactant gas pressures
         */
        double voltage_cell_th();
        
        /**
         * Output operating conditions to screen
         */
        void print_operating_conditions() const;
        
        /**
         * 
         */
        void adjust_initial_solution(std::vector< component_materialID_value_map >& maps,
                                     const boost::shared_ptr< FuelCellShop::Geometry::GridBase<dim> > grid) const;
        /**
         * Routine used in order to adjust boundary conditions for O2 and cell voltage. This routine should be called after
         * the component_boundaryID_value_map has been initialized.
         * 
         * In order for the application to work well, the boundary IDs for the channel/GDL and land/GDL interfaces
         * need to be properly identified in the grid section of the input file, i.e.
         * @code
         * subsection Grid generation
         *   subsection Internal mesh generator parameters
         *     subsection Boundary ID
         *       set c_Ch/GDL = 100
         *       set c_BPP/GDL = 200
         *       set a_Ch/GDL = 100
         *       set a_BPP/GDL = 200
         *     end  
         *   end
         * end
         * @endcode
         */
        void adjust_boundary_conditions(std::vector< component_boundaryID_value_map >& maps,
                                        const boost::shared_ptr< FuelCellShop::Geometry::GridBase<dim> > grid) const;
        /**
         * Get the total gas concentration in the cathode. 
         * NOTE: This is a constant value. Use only if the total gas concentration
         * is assumed to be constant
         */
        inline double get_c_c() const
        {	return c_c; }
        /** Return cell temperture as input in Operating Conditions subsection */
        inline double get_T() const
        { return T_cell; }
        /** Return cell voltage as input in Operating Conditions subsection */
        inline double get_V() const
        { return V_cell; }
        /** Return the voltage drop in the anode 
         * 
         * \note Can be used as boundary condition for anode model or initial condition
         * 		in a full MEA model
         */
        inline double get_dV_a() const
        { return dV_a; }
        /** Return cathode pressure as input in Operating Conditions subsection */
        inline double get_pc_Pa() const
        { return p_c; }
        /** Return cathode pressure as input in Operating Conditions subsection */
        inline double get_pc_atm() const
        { return p_c/Units::convert(1.,Units::ATM_to_PA); }
        /** Get the total gas concentration in the anode. */
        inline double get_c_a() const
        {	return c_a; }
        /** Return anode pressure as input in Operating Conditions subsection */
        inline double get_pa_Pa() const
        { return p_a; }
        /** Return anode pressure as input in Operating Conditions subsection */
        inline double get_pa_atm() const
        { return p_a/Units::convert(1.,Units::ATM_to_PA); }
        /** Return anode relative humidity as input in Operating Conditions subsection */
        inline double get_RH_a() const
        {return RH_a;}
        /** Return cathode relative humidity as input in Operating Conditions subsection */
        inline double get_RH_c() const
        {return RH_c;}
        /** Get the open circuit voltage for the cell */
        inline double get_OCV() const
        {return OCV;}
        
    private:
        /**
         * Private member function used to extract the species index value from name.
         * 
         * Note that density_species_1 will have index 0.
         */
        unsigned int get_species_index(std::string name) const
        {
            AssertThrow(name.compare(0,16,"density_species_") == 0, ExcMessage("This function is only meant to be used to get species number from density_species_*, please use another function."));
            
            std::string speciesStr     = name.substr(16); //get substring from position 16 to end which should be the full species number 
            unsigned int species_index = atoi(speciesStr.c_str()) - 1; //convert species string number to int to use in map
            
            return species_index;
        }
        /**
         * Private function used in adjust_initial_solution and adjust_boundary_conditions that estimates
         * the density of each species in the anode compartment.
         * 
         * The parameter species_index is the index obtained using get_species_index. Note that for the anode,
         * the index is obtained as get_species_index - n_gases in cathode mix.
         */
        double get_anode_gas_density(unsigned int species_index) const
        {                    
            double density(0.0);
            double p_sat_Pa = Units::convert(saturation_pressure(),Units::ATM_to_PA);
            
            unsigned int anode_gas_index = species_index - cathode_mix->n_gases();
            
            const std::string gas_name = anode_mix->get_gas(anode_gas_index)->name_material ();
            double M = 1e3*anode_mix->get_gas(anode_gas_index)->get_molar_mass();  // g/cm3     
            
            //-- Get pV if water is solved for only:
            double pv_a = 0.0;
            for (unsigned int ind_g = 0; ind_g<anode_mix->n_gases(); ind_g++) {
                std::string name = anode_mix->get_gas(ind_g)->name_material ();
                if (name.compare("water") == 0)
                    pv_a = p_sat_Pa*RH_a;
            }
            double p_H2 = channel_dry_hydrogen_mole_fraction*(p_a-pv_a);
            double p_N2 = p_a - p_H2 - pv_a;
            
            if (gas_name.compare("water") == 0)
                density = 1E-6*M*pv_a/(R*T_cell);         // g/cm^3
            else if (gas_name.compare("hydrogen") == 0)
                density = 1E-6*M*p_H2/(R*T_cell);         // g/cm^3
            else if (gas_name.compare("nitrogen") == 0)
                density = 1E-6*M*p_N2/(R*T_cell);         // g/cm^3
            else
                AssertThrow(false, ExcMessage("Gas not implemented"));
            
            return density;
                        
        }
        /**
         * Private function used in adjust_initial_solution and adjust_boundary_conditions that estimates
         * the density of each species in the cathode compartment of the fuel cell.
         * 
         * The parameter species_index is the index obtained using get_species_index. Note that for the anode,
         * the index is obtained as get_species_index in cathode mix.
         */
        double get_cathode_gas_density(unsigned int species_index) const
        {
            double density(0.0);
            double p_sat_Pa = Units::convert(saturation_pressure(),Units::ATM_to_PA);
            
            const std::string gas_name = cathode_mix->get_gas(species_index)->name_material ();
            double M = 1e3*cathode_mix->get_gas(species_index)->get_molar_mass();  // g/cm3  
            
            double pv_c = 0.0;
            for (unsigned int ind_g = 0; ind_g<cathode_mix->n_gases(); ind_g++) 
            {
                std::string name = cathode_mix->get_gas(ind_g)->name_material ();
                if (name.compare("water") == 0)
                    pv_c = p_sat_Pa*RH_c;
            }   
            double p_O2 = channel_oxygen_mole_fraction*(p_c-pv_c);
            double p_N2 = p_a - p_O2 - pv_c;
            
            if (gas_name.compare("water") == 0)
                density = 1E-6*M*pv_c/(R*T_cell);         // g/cm^3
            else if (gas_name.compare("oxygen") == 0)
                density = 1E-6*M*p_O2/(R*T_cell);         // g/cm^3
            else if (gas_name.compare("nitrogen") == 0)
                density = 1E-6*M*p_N2/(R*T_cell);         // g/cm^3
            else
                AssertThrow(false, ExcMessage("Gas not implemented"));
            
            return density;
        }
        
        //------------ BOUNDARY CONDITIONS -------------------------------
        /** Bool set to true if you want to modify boundary conditions */
        bool adjust_BC;
        //------------ CONSTANTS -----------------------------------------
        double R; //Gas constant 8.3144 J / (mol K)
        
        /** Initial amount of oxygen in channel prior to humidification */
        double channel_oxygen_mole_fraction;
        
        /** Initial amount of hydrogen in channel prior to humidification */
        double channel_dry_hydrogen_mole_fraction;
        
        //------------ CELL DATA -----------------------------------------
        /** Operating temperature of the cell */
        double T_cell;
        /** Operating voltage of the cell */
        double V_cell;
        /** Voltage drop in the anode */
        double dV_a;
        /** Open circuit voltage for the cell */
        double OCV;
        /** Theoretical voltage for the cell */
        double E_th;
        
        //------------- ANODE DATA -------------------------------------
        /** Pressure of the gas mixture in the anode B.C. */
        double p_a; 
        /** Concentration of the gas mixture in the anode B.C. */
        double c_a;
        /** Relative humidity of the gas mixture in the anode B.C */
        double RH_a;
        
        //------------- CATHODE DATA -------------------------------------
        /** Pressure of the gas mixture in the cathode B.C. */
        double p_c;
        /** Concentration of the gas mixture in the cathode B.C. */
        double c_c;
        /** Relative humidity of the gas mixture in the anode B.C */
        double RH_c;
        
        //------------- GAS DATA -------------------------------------
        /** Get molar mass of water vapour*/
        FuelCellShop::Material::WaterVapor water;
        const double M_water = water.get_molar_mass() * 1000; // [g/mol]
        /** Get molar mass of oxygen*/
        FuelCellShop::Material::Oxygen oxygen;
        const double M_oxygen = oxygen.get_molar_mass() * 1000; // [g/mol]
        /** Get molar mass of nitrogen*/
        FuelCellShop::Material::Nitrogen nitrogen;
        const double M_nitrogen = nitrogen.get_molar_mass() * 1000; // [g/mol]
        /** Get molar mass of nitrogen*/
        FuelCellShop::Material::Hydrogen hydrogen;
        const double M_hydrogen = hydrogen.get_molar_mass() * 1000; // [g/mol]
        /**
         * Container storing the gas mixture in the anode compartment of the cell
         */
        FuelCellShop::Material::GasMixture* anode_mix;
        /**
         * Container storing the gas mixture in the anode compartment of the cell
         */
        FuelCellShop::Material::GasMixture* cathode_mix;

        /** map relating species number (key) to species material (value)*/
        std::map<unsigned int, std::string> gasSpeciesMap;
        
        /**
        * Run an anode instead of a cathode:
        */
        bool anode = false;
        
        /**
        * Variable to switch between electrodes in molar fraction functions:
        */
        enum Electrode {noElectrode, ANODE, CATHODE, MEMBRANE};
        Electrode compartment;
    };
    
}

#endif // _FUELCELL_OPERATING_CONDITIONS__H
