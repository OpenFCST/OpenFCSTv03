//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: design_MPL.cc
//    - Description: Class used to represent a design MPL where effective properties are computed based on the porosity etc.
//    - Developers: Madhur Bhaiya
//    - $Id: design_MPL.cc 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#include <layers/design_MPL.h>

namespace NAME = FuelCellShop::Layer;

template <int dim>
const std::string NAME::DesignMPL<dim>::concrete_name ("DesignMPL");

template <int dim>
NAME::DesignMPL<dim> const* NAME::DesignMPL<dim>::PROTOTYPE = new NAME::DesignMPL<dim>();

//---------------------------------------------------------------------------
template <int dim>
NAME::DesignMPL<dim>::DesignMPL(std::string name)
  : NAME::MicroPorousLayer<dim>(name)
{
}

//---------------------------------------------------------------------------
template <int dim>
NAME::DesignMPL<dim>::DesignMPL()
  : NAME::MicroPorousLayer<dim>()
{
  //FcstUtilities::log<<" Register DesignMPL MPL to FactoryMap"<<std::endl;
  this->get_mapFactory()->insert(std::pair<std::string, MicroPorousLayer<dim>* > (this->concrete_name, this) ); 
}
//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::declare_parameters (const std::string& mpl_section_name,
                                          ParameterHandler &param) const
{
  
  NAME::MicroPorousLayer<dim>::declare_parameters(mpl_section_name, param);
  
  param.enter_subsection("Fuel cell data"); 
  {
    param.enter_subsection(mpl_section_name); 
    { 
        param.enter_subsection(concrete_name); 
        { 
            // General data
            param.declare_entry ("Porosity",  // volume fraction of void space in the GDL
                                 "0.3", // [-]
                                 Patterns::Double(0., 1.),
                                 "Volume fraction of void space in the GDL");
            
            // Method to compute effective properties
            param.declare_entry ("Method effective transport properties in pores",
                                 "Bruggemann",
                                 Patterns::Selection("Bruggemann|Percolation"),
                                 "Method used to compute effective transport properties in the void phase.");
            param.declare_entry ("Method effective transport properties in solid phase",
                                 "Bruggemann",
                                 Patterns::Selection("Given | Bruggemann | Percolation"),
                                 "Method used to compute effective transport properties in solid phase");
            param.declare_entry ("Method effective thermal conductivity",
                                 "Given",
                                 Patterns::Selection("Given"),
                                 "Method used to compute effective thermal conductivity");
            
            param.declare_entry ("Method relative liquid permeability",
                                 "Kumbur07",
                                 Patterns::Selection("Kumbur07|Wyllie"),
                                 "Method used to compute relative liquid permeability as a function of saturation.");
            param.declare_entry ("Irreducible liquid water saturation",
                                 "0.0",
                                 Patterns::Double(0., 1.),
                                 "Irreducible liquid water saturation in the MPL.");
            
            param.declare_entry ("Method capillary pressure - saturation function",
                                 "Kumbur07-corrected",
                                 Patterns::Selection("Kumbur07-corrected"),
                                 "Method used to compute capillary pressure as a function of saturation.");
            param.declare_entry ("Compaction pressure [MPa]",
                                 "0.0",
                                 Patterns::Double(0.),
                                 "Compaction pressure on the MPL, Units: [MPa].");
            param.declare_entry ("PTFE loading [% wt]",
                                 "5.0",
                                 Patterns::Double(5., 20.),
                                 "PTFE loading (% wt) in the MPL; Accepted range: 5 - 20 %.");
            
            //-- Diffusivity for the isotropic case:
            param.declare_entry ("Porosity threshold",  // volume fraction of void space in the GDL
                                 "0.12", // [-]
                                 Patterns::Double(),
                                 "Threshold value of the volume fraction of void space in the GDL."
                                 "If the porosity is less than this value transport does not occur");
            param.declare_entry ("Porosity network constant",
                                 "2.0",
                                 Patterns::Double(),
                                 "Parameter used when using percolation theory");
            param.declare_entry ("Porosity gamma network constant",
                                 "0.0",
                                 Patterns::Double(),
                                 "Parameter used when using percolation theory to account for extra diffusion");
            param.declare_entry ("Oxygen diffusion coefficient",
                                 "0.2741", //1atm 353K
                                 Patterns::Double());
            param.declare_entry ("Water vapour diffusion coefficient",
                                 "0.29676", 
                                 Patterns::Double());
            param.declare_entry ("Solid network threshold",  // volume fraction of void space in the GDL
                                 "0.12", // [-]
                                 Patterns::Double(),
                                 "Threshold value of the volume fraction of solid (electron conductive) phase in the MPL."
                                 "If the solid phase is less than this value transport in the fibre network does not occur");
            param.declare_entry ("Solid network constant",
                                 "2.0",
                                 Patterns::Double(),
                                 "Parameter used when using percolation theory");
            param.declare_entry ("Electric conductivity", 
                                 "10",
                                 Patterns::Double(),
                                 "Either effective (if Given method used) or solid electric conductivity [S/cm]");
            param.declare_entry ("Thermal conductivity", 
                                 "10",
                                 Patterns::Double(),
                                 "Either effective (if Given method used) or empirical methods [W/cm-K]");
            param.declare_entry ("Absolute permeability [cm^2]",
                                 "1.5e-9", // cm^2
                                 Patterns::Double(0.),
                                 "Absolute permeability of the layer, Units [cm^2]");
            
            //-- Diffusivity for the anisotropic case:
            param.declare_entry ("Anisotropic transport",
                                 "false", // [S/cm]
                                 Patterns::Bool(),
                                 "Boolean variable. Set to true if we want to account for anisotropy of the MPL");
            //--- XX
            param.declare_entry ("Electric conductivity X", 
                                 "10",
                                 Patterns::Double(),
                                 "Either effective (if Given method used) or solid electric conductivity [S/cm]");
            param.declare_entry ("Thermal conductivity X", 
                                 "10",
                                 Patterns::Double(),
                                 "Either effective (if Given method used) or empirical methods [W/cm-K]");
            param.declare_entry ("Porosity threshold X",  // volume fraction of void space in the GDL
                                 "0.12", // [-]
                                 Patterns::Double(),
                                 "Threshold value of the volume fraction of void space in the GDL."
                                 "If the porosity is less than this value transport does not occur");
            param.declare_entry ("Porosity network constant X",
                                 "2.0",
                                 Patterns::Double(),
                                 "Parameter used when using percolation theory");
            param.declare_entry ("Porosity gamma network constant X",
                                 "0.0",
                                 Patterns::Double(),
                                 "Parameter used when using percolation theory to account for extra diffusion");
            param.declare_entry ("Oxygen diffusion coefficient X",
                                 "0.2741", //1atm 353K
                                 Patterns::Double());
            param.declare_entry ("Water vapour diffusion coefficient X",
                                 "0.29676", 
                                 Patterns::Double());
            param.declare_entry ("Fibre network threshold X",  // volume fraction of void space in the GDL
                                 "0.12", // [-]
                                 Patterns::Double(),
                                 "Threshold value of the volume fraction of fibres in the GDL."
                                 "If the solid phase is less than this value transport in the fibre network does not occur");
            param.declare_entry ("Fibre network constant X",
                                 "2.0",
                                 Patterns::Double(),
                                 "Parameter used when using percolation theory");
            param.declare_entry ("Absolute permeability X [cm^2]",
                                 "1.5e-9", // cm^2
                                 Patterns::Double(0.),
                                 "X component - Absolute permeability of the layer, Units [cm^2]");
            
            // YY
            param.declare_entry ("Electric conductivity Y", 
                                 "10",
                                 Patterns::Double(),
                                 "Either effective (if Given method used) or solid electric conductivity [S/cm]");
            param.declare_entry ("Thermal conductivity Y", 
                                 "10",
                                 Patterns::Double(),
                                 "Either effective (if Given method used) or empirical methods [W/cm-K]");
            param.declare_entry ("Porosity threshold Y",  // volume fraction of void space in the GDL
                                 "0.12", // [-]
                                 Patterns::Double(),
                                 "Threshold value of the volume fraction of void space in the GDL."
                                 "If the porosity is less than this value transport does not occur");
            param.declare_entry ("Porosity network constant Y",
                                 "2.0",
                                 Patterns::Double(),
                                 "Parameter used when using percolation theory");
            param.declare_entry ("Porosity gamma network constant Y",
                                 "0.0",
                                 Patterns::Double(),
                                 "Parameter used when using percolation theory to account for extra diffusion");
            param.declare_entry ("Oxygen diffusion coefficient Y",
                                 "0.2741", //1atm 353K
                                 Patterns::Double());
            param.declare_entry ("Water vapour diffusion coefficient Y",
                                 "0.29676", 
                                 Patterns::Double());
            param.declare_entry ("Fibre network threshold Y",  // volume fraction of void space in the GDL
                                 "0.12", // [-]
                                 Patterns::Double(),
                                 "Threshold value of the volume fraction of fibres in the GDL."
                                 "If the solid phase is less than this value transport in the fibre network does not occur");
            param.declare_entry ("Fibre network constant Y",
                                 "2.0",
                                 Patterns::Double(),
                                 "Parameter used when using percolation theory");
            param.declare_entry ("Absolute permeability Y [cm^2]",
                                 "1.5e-9", // cm^2
                                 Patterns::Double(0.),
                                 "Y component - Absolute permeability of the layer, Units [cm^2]");
            
            // ZZ
            param.declare_entry ("Electric conductivity Z", 
                                 "10",
                                 Patterns::Double(),
                                 "Either effective (if Given method used) or solid electric conductivity [S/cm]");
            param.declare_entry ("Thermal conductivity Z", 
                                 "10",
                                 Patterns::Double(),
                                 "Either effective (if Given method used) or empirical methods [W/cm-K]");
            param.declare_entry ("Porosity threshold Z",  // volume fraction of void space in the GDL
                                 "0.12", // [-]
                                 Patterns::Double(),
                                 "Threshold value of the volume fraction of void space in the GDL."
                                 "If the porosity is less than this value transport does not occur");
            param.declare_entry ("Porosity network constant Z",
                                 "2.0",
                                 Patterns::Double(),
                                 "Parameter used when using percolation theory");
            param.declare_entry ("Porosity gamma network constant Z",
                                 "0.0",
                                 Patterns::Double(),
                                 "Parameter used when using percolation theory to account for extra diffusion");
            param.declare_entry ("Oxygen diffusion coefficient Z",
                                 "0.2741", //1atm 353K
                                 Patterns::Double());
            param.declare_entry ("Water vapour diffusion coefficient Z",
                                 "0.29676", 
                                 Patterns::Double());
            param.declare_entry ("Fibre network threshold Z",  // volume fraction of void space in the GDL
                                 "0.12", // [-]
                                 Patterns::Double(),
                                 "Threshold value of the volume fraction of fibres in the GDL."
                                 "If the solid phase is less than this value transport in the fibre network does not occur");
            param.declare_entry ("Fibre network constant Z",
                                 "2.0",
                                 Patterns::Double(),
                                 "Parameter used when using percolation theory");
            param.declare_entry ("Absolute permeability Z [cm^2]",
                                 "1.5e-9", // cm^2
                                 Patterns::Double(0.),
                                 "Z component - Absolute permeability of the layer, Units [cm^2]");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
  }
  param.leave_subsection();
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::set_parameters (const std::vector<std::string>& name_dvar,
                                      const std::vector<double>& value_dvar,
                                      ParameterHandler &param)
{
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection(this->name);
        {
            param.enter_subsection(concrete_name);
            {
                for (unsigned int i=0; i<name_dvar.size(); ++i)
                {
                    if (name_dvar[i] == "D_O2_CMPL")
                    {
                        if(this->name == "Cathode microporous layer")
                            param.set("Oxygen diffusion coefficient", value_dvar[i]);
                    }
                    else if (name_dvar[i] == "D_H2O_CMPL")
                    {
                        if(this->name == "Cathode microporous layer")
                            param.set("Water vapour diffusion coefficient", value_dvar[i]);
                    }
                    else if (name_dvar[i] == "D_H2O_AMPL")
                    {
                        if(this->name == "Anode microporous layer")
                            param.set("Water vapour diffusion coefficient", value_dvar[i]);
                    }
                    else if (name_dvar[i] == "sigma_e_MPL_X")
                    {
                        param.set("Electric conductivity X", value_dvar[i]);
                    }
                    else if (name_dvar[i] == "sigma_e_MPL_Y")
                    {
                        param.set("Electric conductivity Y", value_dvar[i]);
                    }
                    else if (name_dvar[i] == "sigma_e_MPL")
                    {
                        param.set("Electric conductivity", value_dvar[i]);
                    }
                    else if (name_dvar[i] == "thermal_conductivity_MPL")
                    {
                        param.set("Thermal conductivity", value_dvar[i]);
                    }
                    else if (name_dvar[i] == "thermal_conductivity_MPL_X")
                    {
                        param.set("Thermal conductivity X", value_dvar[i]);
                    }
                    else if (name_dvar[i] == "thermal_conductivity_MPL_Y")
                    {
                        param.set("Thermal conductivity Y", value_dvar[i]);
                    }
                    else if (name_dvar[i] == "thermal_conductivity_MPL_Z")
                    {
                        param.set("Thermal conductivity Z", value_dvar[i]);
                    }
                }
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::initialize (ParameterHandler &param)
{
    NAME::MicroPorousLayer<dim>::initialize(param);
    
    param.enter_subsection("Fuel cell data"); 
    {
        param.enter_subsection(this->name); 
        {
            param.enter_subsection(concrete_name); 
            { 
                porosity = param.get_double("Porosity");
                electrical_conductivity = param.get_double("Electric conductivity");
                thermal_conductivity = param.get_double("Thermal conductivity");
                solid_phase = 1 - porosity;
                method_eff_property_pores = param.get("Method effective transport properties in pores");
                method_eff_property_fibres = param.get("Method effective transport properties in solid phase");
                method_eff_thermal_conductivity = param.get("Method effective thermal conductivity");
                
                method_rel_liquid_permeability = param.get("Method relative liquid permeability");
                s_irr = param.get_double("Irreducible liquid water saturation");
                
                method_capillary_function = param.get("Method capillary pressure - saturation function");
                compaction_pressure = param.get_double("Compaction pressure [MPa]");
                PTFE_loading = param.get_double("PTFE loading [% wt]");
                
                // Anisotropy
                anisotropy = param.get_bool ("Anisotropic transport");
                
                matrix_electrical_conductivity.clear();
                matrix_thermal_conductivity.clear();
                porosity_th.resize(dim);
                porosity_mu.resize(dim);
                porosity_gamma.resize(dim);
                D_O2.resize(dim);
                D_wv.resize(dim);
                fibre_th.resize(dim);
                fibre_mu.resize(dim);
                abs_permeability.resize(dim);
                
                if (anisotropy == true)
                {
                    if (dim == 1)
                    {
                        // X
                        matrix_electrical_conductivity[0][0] = param.get_double("Electric conductivity X");
                        matrix_thermal_conductivity[0][0] = param.get_double("Thermal conductivity X");
                        porosity_th[0] = param.get_double("Porosity threshold X");
                        porosity_mu[0] = param.get_double("Porosity network constant X");
                        porosity_gamma[0] = param.get_double("Porosity gamma network constant X");
                        D_O2[0] = param.get_double("Oxygen diffusion coefficient X");
                        D_wv[0] = param.get_double("Water vapour diffusion coefficient X"); 
                        fibre_th[0] = param.get_double("Fibre network threshold X");
                        fibre_mu[0] = param.get_double("Fibre network constant X");
                        abs_permeability[0] = param.get_double("Absolute permeability X [cm^2]");
                    }
                    else if (dim == 2)
                    {
                        // X
                        matrix_electrical_conductivity[0][0] = param.get_double("Electric conductivity X");
                        matrix_thermal_conductivity[0][0] = param.get_double("Thermal conductivity X");
                        porosity_th[0] = param.get_double("Porosity threshold X");
                        porosity_mu[0] = param.get_double("Porosity network constant X");
                        porosity_gamma[0] = param.get_double("Porosity gamma network constant X");
                        D_O2[0] = param.get_double("Oxygen diffusion coefficient X");
                        D_wv[0] = param.get_double("Water vapour diffusion coefficient X"); 
                        fibre_th[0] = param.get_double("Fibre network threshold X");
                        fibre_mu[0] = param.get_double("Fibre network constant X");
                        abs_permeability[0] = param.get_double("Absolute permeability X [cm^2]");
                        
                        // Y
                        matrix_electrical_conductivity[1][1] = param.get_double("Electric conductivity Y");
                        matrix_thermal_conductivity[1][1] = param.get_double("Thermal conductivity Y");
                        porosity_th[1] = param.get_double("Porosity threshold Y");
                        porosity_mu[1] = param.get_double("Porosity network constant Y");
                        porosity_gamma[1] = param.get_double("Porosity gamma network constant Y");
                        D_O2[1] = param.get_double("Oxygen diffusion coefficient Y");
                        D_wv[1] = param.get_double("Water vapour diffusion coefficient Y"); 
                        fibre_th[1] = param.get_double("Fibre network threshold Y");
                        fibre_mu[1] = param.get_double("Fibre network constant Y");
                        abs_permeability[1] = param.get_double("Absolute permeability Y [cm^2]");
                    }
                    else if (dim == 3)
                    {
                        // X
                        matrix_electrical_conductivity[0][0] = param.get_double("Electric conductivity X");
                        matrix_thermal_conductivity[0][0] = param.get_double("Thermal conductivity X");
                        porosity_th[0] = param.get_double("Porosity threshold X");
                        porosity_mu[0] = param.get_double("Porosity network constant X");
                        porosity_gamma[0] = param.get_double("Porosity gamma network constant X");
                        D_O2[0] = param.get_double("Oxygen diffusion coefficient X");
                        D_wv[0] = param.get_double("Water vapour diffusion coefficient X"); 
                        fibre_th[0] = param.get_double("Fibre network threshold X");
                        fibre_mu[0] = param.get_double("Fibre network constant X");
                        abs_permeability[0] = param.get_double("Absolute permeability X [cm^2]");
                        
                        // Y
                        matrix_electrical_conductivity[1][1] = param.get_double("Electric conductivity Y");
                        matrix_thermal_conductivity[1][1] = param.get_double("Thermal conductivity Y");
                        porosity_th[1] = param.get_double("Porosity threshold Y");
                        porosity_mu[1] = param.get_double("Porosity network constant Y");
                        porosity_gamma[1] = param.get_double("Porosity gamma network constant Y");
                        D_O2[1] = param.get_double("Oxygen diffusion coefficient Y");
                        D_wv[1] = param.get_double("Water vapour diffusion coefficient Y"); 
                        fibre_th[1] = param.get_double("Fibre network threshold Y");
                        fibre_mu[1] = param.get_double("Fibre network constant Y");
                        abs_permeability[1] = param.get_double("Absolute permeability Y [cm^2]");
                        
                        // Z
                        matrix_electrical_conductivity[2][2] = param.get_double("Electric conductivity Z");
                        matrix_thermal_conductivity[2][2] = param.get_double("Thermal conductivity Z");
                        porosity_th[2] = param.get_double("Porosity threshold Z");
                        porosity_mu[2] = param.get_double("Porosity network constant Z");
                        porosity_gamma[2] = param.get_double("Porosity gamma network constant Z");
                        D_O2[2] = param.get_double("Oxygen diffusion coefficient Z");
                        D_wv[2] = param.get_double("Water vapour diffusion coefficient Z"); 
                        fibre_th[2] = param.get_double("Fibre network threshold Z");
                        fibre_mu[2] = param.get_double("Fibre network constant Z");
                        abs_permeability[2] = param.get_double("Absolute permeability Z [cm^2]");
                    }
                }
                else
                {
                    for (unsigned int i =0; i < porosity_th.size(); i++)
                    {
                        matrix_electrical_conductivity[i][i] = param.get_double("Electric conductivity");
                        matrix_thermal_conductivity[i][i] = param.get_double("Thermal conductivity");
                        porosity_th[i] = param.get_double("Porosity threshold");
                        porosity_mu[i] = param.get_double("Porosity network constant");
                        porosity_gamma[i] = param.get_double("Porosity gamma network constant");
                        D_O2[i] = param.get_double("Oxygen diffusion coefficient");
                        D_wv[i] = param.get_double("Water vapour diffusion coefficient"); 
                        fibre_th[i] = param.get_double("Solid network threshold");
                        fibre_mu[i] = param.get_double("Solid network constant");
                        abs_permeability[i] = param.get_double("Absolute permeability [cm^2]");
                    }
                }
                
                // ----------- Computing Kumbur factor (may be required later for capillary pressure calculations)
                double s_tr = ((-0.0083)*compaction_pressure*compaction_pressure) + (0.0911*compaction_pressure);
                double compressed_porosity = porosity * ( (0.9/(1.0 + s_tr)) + 0.1 );
                kumbur_factor = (std::pow(2.0, (0.4*compaction_pressure))) * (std::pow((compressed_porosity/(abs_permeability[0]*0.0001)), 0.5));
                // -----------
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();  
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::effective_gas_diffusivity(Table< 2, Tensor< 2, dim > > &D_eff) const
{
    D_eff.reinit(this->gases.size(),this->gases.size());    
    
    for (unsigned int i=0; i<this->gases.size(); i++)
    {
        for (unsigned int j=0; j<this->gases.size(); j++)
        {
            if( i != j )        // Not computing self diffusion coefficients
            {
                this->effective_gas_diffusivity( this->D_ECtheory(i,j), 0.0, D_eff(i,j) );
            }
        }
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::effective_gas_diffusivity(const double& prop,
                                                const double& saturation,
                                                double& prop_eff ) const
{
    // Assert if anisotropy is false
    AssertThrow(anisotropy == false, ExcMessage("DesignMPL::effective_transport_property_pores with double return value can only be used for isotropic GDLs"));
    Assert( saturation >= 0.0, ExcMessage("Non-negative saturation values must be passed as an argument in DesignMPL::effective_gas_diffusivity.") );
    
    double effective_porosity = porosity * (1.0 - saturation);
    
    if (method_eff_property_pores == "Bruggemann")
    {
        // Most common method
        prop_eff = prop * std::pow(effective_porosity, 1.5);
    }
    
    else if (method_eff_property_pores == "Percolation")
    {
        // Method used by M. Eikerling 
        // Note: For porosity_th[0] I use index 0 for convenience. They all contain the same value
        double aux = (std::fabs(effective_porosity - porosity_th[0])) / (1.0 - porosity_th[0]);
        double step = 1.0e-3;
        
        if (effective_porosity >= porosity_th[0])
            step = 1.0;
        
        prop_eff = prop * ( std::pow(aux, porosity_mu[0]) * step );
    }
    
    else
        AssertThrow( false, ExcNotImplemented() );
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::effective_gas_diffusivity(const double& prop,
                                                const double& saturation,
                                                Tensor<2,dim>& prop_eff) const
{
    Assert( saturation >= 0.0, ExcMessage("Non-negative saturation values must be passed as an argument in DesignMPL::effective_gas_diffusivity.") );
    
    // Initialize variable
    prop_eff.clear();
    
    double effective_porosity = porosity * (1.0 - saturation);
    
    if (method_eff_property_pores == "Bruggemann")
    {
        // Most common method
        for (unsigned int i=0; i<dim; i++)
            prop_eff[i][i] = prop * std::pow(effective_porosity, 1.5);
    }
    
    else if (method_eff_property_pores == "Percolation")
    {
        // Method used by M. Eikerling
        for (unsigned int i=0; i<dim; i++)
        {
            double aux = (std::fabs(effective_porosity - porosity_th[i])) / (1.0 - porosity_th[i]);
            double step = 1.0e-3;
            
            if (effective_porosity >= porosity_th[i])
                step = 1.0;
            
            prop_eff[i][i] = prop * ( std::pow(aux, porosity_mu[i]) * step );
        }
    }
    
    else
        AssertThrow( false, ExcNotImplemented() );
}

//---------------------------------------------------------------------------
template<int dim>
void
NAME::DesignMPL<dim>::effective_gas_diffusivity(std::vector< Tensor<2,dim> >& prop_eff_vec) const
{
    Assert( this->D_bulk.size() != 0, ExcMessage("compute_gas_diffusion not called before DesignMPL::effective_gas_diffusivity.") );
    
    prop_eff_vec.resize( this->D_bulk.size() );
    
    for (unsigned int i = 0; i < this->D_bulk.size(); ++i)
    {
        if ( this->s_vector.is_initialized() )
        {
            if ( this->s_vector[i] >= 0.0 )
                this->effective_gas_diffusivity( this->D_bulk[i], this->s_vector[i], prop_eff_vec[i] );
            else
                this->effective_gas_diffusivity( this->D_bulk[i], 0.0, prop_eff_vec[i] );
        }
        else if (this->capillary_pressure_vector.is_initialized())
        {
            std::vector<double> saturation;
            this->psd_pointer->set_capillary_pressure( this->capillary_pressure_vector );
            this->psd_pointer->set_saturation();
            this->psd_pointer->get_saturation(saturation);
            
            for (unsigned int j=0; j < dim; ++j)
                this->effective_gas_diffusivity( this->D_bulk[i], saturation[i], prop_eff_vec[i] );
        }
            
        else                
            this->effective_gas_diffusivity( this->D_bulk[i], 0.0, prop_eff_vec[i] );
    }
}

//---------------------------------------------------------------------------
template<int dim>
void
NAME::DesignMPL<dim>::derivative_effective_gas_diffusivity(std::map< VariableNames, std::vector< Tensor<2,dim> > >& dprop_eff) const
{
    Assert( this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called in DesignMPL::derivative_effective_gas_diffusivity."));
    Assert(((this->dD_bulk_dT.size()!=0) && (this->D_bulk.size()!=0)), ExcMessage("compute_gas_diffusion not called before DesignMPL::derivative_effective_gas_diffusivity."));
    
    for (unsigned int i=0; i < this->derivative_flags.size(); ++i)
    {
        std::vector< Tensor<2,dim> > Dprop( this->dD_bulk_dT.size(), Tensor<2,dim>() );
        
        if ( this->derivative_flags[i] == temperature_of_REV )
        {
            for (unsigned int j=0; j < this->dD_bulk_dT.size(); ++j)
            {
                if ( this->s_vector.is_initialized() )
                {
                    if ( this->s_vector[j] >= 0.0 )
                        this->effective_gas_diffusivity( this->dD_bulk_dT[j], this->s_vector[j], Dprop[j] );
                    else
                        this->effective_gas_diffusivity( this->dD_bulk_dT[j], 0.0, Dprop[j] );
                }
                
                else
                    this->effective_gas_diffusivity( this->dD_bulk_dT[j], 0.0, Dprop[j] );
            }
            //-----
            dprop_eff[ this->derivative_flags[i] ] = Dprop;
        }
        
        else if ( this->derivative_flags[i] == liquid_water_saturation )
        {
            Assert( this->s_vector.is_initialized(), ExcInternalError() );
            
            if (method_eff_property_pores == "Bruggemann")
            {
                for (unsigned int j=0; j < this->D_bulk.size(); ++j)
                {
                    if ( this->s_vector[j] >= 0.0 )
                    {
                        for (unsigned int k=0; k<dim; ++k)
                            Dprop[j][k][k] = this->D_bulk[j] * (-1.5) * (std::pow(porosity,1.5)) * (std::pow((1.0 - this->s_vector[j]),0.5));
                    }
                    // ELSE; for non-physical negative saturation values, 's' is taken as zero hence Deff = constant thus dDeff/ds = 0.0
                }
            }
            
            else if (method_eff_property_pores == "Percolation")
            {
                for (unsigned int j=0; j < this->D_bulk.size(); ++j)
                {
                    if (this->s_vector[j] >= 0.0)
                    {
                        for (unsigned int k=0; k<dim; ++k)
                        {
                            if ( (porosity*(1.0 - this->s_vector[j])) >= porosity_th[k] )
                            {
                                Dprop[j][k][k] = this->D_bulk[j] * (-1.0) * porosity * porosity_mu[k] * (std::pow( (1.0-porosity_th[k]), ((-1.0)*(porosity_mu[k])) )) * 
                                                 (std::pow( ((porosity*(1.0-this->s_vector[j])) - porosity_th[k]), (porosity_mu[k]-1.0) ));
                            }
                        }
                    }
                }
            }
            
            else
                AssertThrow( false, ExcNotImplemented() );
            
            //-----
            dprop_eff[ this->derivative_flags[i] ] = Dprop;
        }
        
        else if ( this->derivative_flags[i] == capillary_pressure )
        {
            Assert( this->capillary_pressure_vector.is_initialized(), ExcInternalError() );
            
            std::vector<double> saturation;
            std::vector<double> ds_dp;
            this->psd_pointer->set_capillary_pressure( this->capillary_pressure_vector );
            this->psd_pointer->set_saturation();
            this->psd_pointer->get_saturation(saturation);
            this->psd_pointer->get_derivative_saturation(ds_dp);
            
            if (method_eff_property_pores == "Bruggemann")
            {
                for (unsigned int j=0; j < this->D_bulk.size(); ++j)
                {
                    if ( saturation[j] >= 0.0 )
                    {
                        for (unsigned int k=0; k<dim; ++k)
                            Dprop[j][k][k] = this->D_bulk[j] * ds_dp[j] * (-1.5) * (std::pow(porosity,1.5)) * (std::pow((1.0 - saturation[j]),0.5));
                    }
                    // ELSE; for non-physical negative saturation values, 's' is taken as zero hence Deff = constant thus dDeff/ds = 0.0
                }
            }
            
            else if (method_eff_property_pores == "Percolation")
            {
                for (unsigned int j=0; j < this->D_bulk.size(); ++j)
                {
                    if (saturation[j] >= 0.0)
                    {
                        for (unsigned int k=0; k<dim; ++k)
                        {
                            if ( (porosity*(1.0 - saturation[j])) >= porosity_th[k] )
                            {
                                Dprop[j][k][k] = this->D_bulk[j] * ds_dp[j] * (-1.0) * porosity * porosity_mu[k] * (std::pow( (1.0-porosity_th[k]), ((-1.0)*(porosity_mu[k])) )) * 
                                (std::pow( ((porosity*(1.0-saturation[j])) - porosity_th[k]), (porosity_mu[k]-1.0) ));
                            }
                        }
                    }
                }
            }
            
            else
                AssertThrow( false, ExcNotImplemented() );
            
            //-----
            dprop_eff[ this->derivative_flags[i] ] = Dprop;
        }
        
        else
            dprop_eff[ this->derivative_flags[i] ] = Dprop;
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::effective_electron_conductivity(double& prop_eff) const
{
  effective_transport_property_solid(electrical_conductivity, 
                                     prop_eff);
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::effective_electron_conductivity(Tensor<2,dim>& prop_eff) const
{
  effective_transport_property_solid(matrix_electrical_conductivity, 
                                     prop_eff);
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::effective_thermal_conductivity(double& prop_eff) const
{
    if (method_eff_thermal_conductivity == "Given")
        prop_eff = thermal_conductivity;
    
    else
        AssertThrow( false, ExcNotImplemented() );

}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::effective_thermal_conductivity(Tensor<2,dim>& prop_eff) const
{
    if (method_eff_thermal_conductivity == "Given")
        prop_eff = matrix_thermal_conductivity;     

    else
        AssertThrow( false, ExcNotImplemented() );
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::effective_thermal_conductivity(std::vector< Tensor<2,dim> >& prop_eff) const
{
    Assert(prop_eff.size()!=0, ExcMessage("Effective thermal conductivity vector not initialized."));
    
    if (method_eff_thermal_conductivity == "Given")
        for (unsigned int i=0; i<prop_eff.size(); ++i)
            prop_eff[i] = matrix_thermal_conductivity;
    
    else
        AssertThrow( false, ExcNotImplemented() );
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::effective_transport_property_solid(const double& prop, 
                                                         double& prop_eff) const
{
    // Asset if anisotropy is false
    Assert(anisotropy == false, ExcMessage("DesignMPL::effective_transport_property_pores with double return value"
    " can only be used for isotropic GDLs"));
    // Initialize variable
    prop_eff = 0.0;
    
    // In some papers prop_eff is give directly
    if (method_eff_property_fibres == "Given")
    {
        prop_eff = prop;
    }
    else if (method_eff_property_fibres == "Bruggemann")
        // Most common method
    {
        prop_eff = prop*pow(solid_phase, 1.5);
    }
    else if (method_eff_property_fibres == "Percolation")
        // Method used by M. Eikerling
    {
        unsigned int i = 0;
        double aux = fabs(solid_phase - fibre_th[i])/(1 - fibre_th[i]);
        double step = 0.0;
        if (solid_phase >= fibre_th[i])
            step = 1.0;
        prop_eff = prop*pow(aux,fibre_mu[i])*step;
    }
    else
    {
        FcstUtilities::log << "Unknown method to compute effective conduction in "<<__FILE__ <<" line "<<__LINE__<<std::endl;
        abort();
    }
}
                                                         
//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::effective_transport_property_solid(const Tensor<2,dim>& prop, 
                                                         Tensor<2,dim>& prop_eff) const
{
    // Initialize variable
    prop_eff = 0;
    
    // In some papers prop_eff is give directly
    if (method_eff_property_fibres == "Given")
    {
        for (unsigned int i=0; i<dim; i++)
            prop_eff[i][i] = prop[i][i];
    }
    else if (method_eff_property_fibres == "Bruggemann")
    {
        // Most common method
        for (unsigned int i=0; i<dim; i++)
            prop_eff[i][i] = prop[i][i]*pow(solid_phase, 1.5);
    }
    // Method used by M. Eikerling
    else if (method_eff_property_fibres == "Percolation")
    {
        for (unsigned int i=0; i<dim; i++)
        {
            double aux = fabs(solid_phase - fibre_th[i])/(1 - fibre_th[i]);
            double step = 0.0;
            if (solid_phase >= fibre_th[i])
                step = 1.0;
            prop_eff[i][i] = prop[i][i]*pow(aux,fibre_mu[i])*step;
        }
    }
    else
    {
        FcstUtilities::log << "Unknown method to compute effective conduction in "<<__FILE__ <<" line "<<__LINE__<<std::endl;
        abort();
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::liquid_permeablity(std::vector< Tensor<2,dim> >& k_l ) const
{
    Assert( this->s_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in DesignMPL::liquid_permeablity.") );
    k_l.clear();
    k_l.resize(this->s_vector.size());
    
    if (method_rel_liquid_permeability == "Kumbur07")
    // Ref: Kumbur E, Sharp K and Mench M. On the effectiveness of Leverett approach for describing the water transport in fuel cell diffusion media. Journal of Power Sources, 168(2):356-368, 2007.
    {
        for (unsigned int i=0; i<this->s_vector.size(); ++i)
        {
            for (unsigned int j=0; j<dim; ++j)
            {
                if (this->s_vector[i] < s_irr)
                    k_l[i][j][j] = 0.0;
                
                else
                {
                    double eff_saturation = (this->s_vector[i] - s_irr)/(1.0 - s_irr);
                    k_l[i][j][j] = abs_permeability[j] * ( std::pow(eff_saturation, 2.16) );
                }
            }
        }
    }
    
    else if (method_rel_liquid_permeability == "Wyllie")
    // Ref: Nam J and Kaviany M. Effective diffusivity and water-saturation distribution in single and two-layer PEMFC diffusion medium. International Journal of Heat and Mass Transfer, 46(24):4595-4611, 2003.
    {
        for (unsigned int i=0; i<this->s_vector.size(); ++i)
        {
            for (unsigned int j=0; j<dim; ++j)
            {
                if (this->s_vector[i] < s_irr)
                    k_l[i][j][j] = 0.0;
                
                else
                {
                    double eff_saturation = (this->s_vector[i] - s_irr)/(1.0 - s_irr);
                    k_l[i][j][j] = abs_permeability[j] * ( std::pow(eff_saturation, 3.0) );
                }
            }
        }
    }
    
    else
        AssertThrow( false, ExcNotImplemented() );

}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::derivative_liquid_permeablity(std::map< VariableNames, std::vector< Tensor<2,dim> > >& deriv_k_l) const
{
    Assert( this->s_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in DesignMPL::derivative_liquid_permeablity.") );
    Assert (this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called in DesignMPL::derivative_liquid_permeablity."));
    
    for (unsigned int i=0; i<this->derivative_flags.size(); ++i)
    {
        std::vector< Tensor<2,dim> > dk_l( this->s_vector.size(), Tensor<2,dim>() );
        
        if (this->derivative_flags[i] == liquid_water_saturation)
        {
            if (method_rel_liquid_permeability == "Kumbur07")
            {
                for (unsigned int j=0; j<this->s_vector.size(); ++j)
                {
                    for (unsigned int k=0; k<dim; ++k)
                    {
                        if (this->s_vector[j] < s_irr)
                            dk_l[j][k][k] = 0.0;
                        
                        else
                            dk_l[j][k][k] = abs_permeability[k] * 2.16 * ( std::pow((this->s_vector[j] - s_irr), 1.16) ) * ( std::pow((1.0 - s_irr), -2.16) );
                    }
                }
            }
            
            else if (method_rel_liquid_permeability == "Wyllie")
            {
                for (unsigned int j=0; j<this->s_vector.size(); ++j)
                {
                    for (unsigned int k=0; k<dim; ++k)
                    {
                        if (this->s_vector[j] < s_irr)
                            dk_l[j][k][k] = 0.0;
                        
                        else
                            dk_l[j][k][k] = abs_permeability[k] * 3.0 * ( std::pow((this->s_vector[j] - s_irr), 2.0) ) * ( std::pow((1.0 - s_irr), -3.0) );
                    }
                }
            }
            
            else
                AssertThrow( false, ExcNotImplemented() );
            
            //---------------
            deriv_k_l[ this->derivative_flags[i] ] = dk_l;
        }
        
        else
            deriv_k_l[ this->derivative_flags[i] ] = dk_l;
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::saturated_liquid_permeablity_PSD(double& r_lp) const
{
    Assert( this->capillary_pressure_vector.is_initialized(), ExcMessage("Liquid water capillary values have not been set using set_capillary_pressure method in ConventionalCL::liquid_permeablity.") );
    
    this->psd_pointer->set_capillary_pressure( this->capillary_pressure_vector );
    this->psd_pointer->set_critical_radius();
    this->psd_pointer->set_saturation();
    this->psd_pointer->get_global_saturated_permeability(r_lp);
    

}
//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::relative_liquid_permeability_PSD(std::vector< Tensor<2,dim> >& r_lp ) const
{
    Assert( this->capillary_pressure_vector.is_initialized(), ExcMessage("Liquid water capillary values have not been set using set_capillary_pressure method in ConventionalCL::liquid_permeablity.") );
    r_lp.clear();
    r_lp.resize(this->capillary_pressure_vector.size(), Tensor<2,dim>());
    
    this->psd_pointer->set_capillary_pressure( this->capillary_pressure_vector );
    
    this->psd_pointer->set_critical_radius();
    this->psd_pointer->set_saturation();
    
    std::vector<double> rlp;
    
    this->psd_pointer->get_relative_liquid_permeability(rlp);
    
    for (unsigned int i=0; i<this->capillary_pressure_vector.size(); ++i)
    {
        for (unsigned int j=0; j<dim; ++j)
        {
            for (unsigned int k=0; k<dim; ++k)
            {
                if (j == k)
                    r_lp[i][j][k] = rlp[i];
                else 
                    r_lp[i][j][k] = 0.0;
            }
        }
    }

}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::derivative_relative_liquid_permeablity_PSD(std::vector<double> & derivate_r_lp ) const
{
    Assert( this->capillary_pressure_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in ConventionalCL::derivative_liquid_permeablity.") );
    Assert (this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called in ConventionalCL::derivative_liquid_permeablity."));
    
    this->psd_pointer->set_capillary_pressure(this->capillary_pressure_vector);
    this->psd_pointer->set_critical_radius();
    this->psd_pointer->set_saturation();
    
    this->psd_pointer->get_derivative_relative_liquid_permeability(derivate_r_lp);
    
}
//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::derivative_relative_liquid_permeablity_PSD(std::map< VariableNames, std::vector< Tensor<2,dim> > >& derivate_r_lp ) const
{
    Assert( this->capillary_pressure_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in ConventionalCL::derivative_liquid_permeablity.") );
    Assert (this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called in ConventionalCL::derivative_liquid_permeablity."));
    
    for (unsigned int i=0; i<this->derivative_flags.size(); ++i)
    {
        std::vector< Tensor<2,dim> > dk_l( this->capillary_pressure_vector.size(), Tensor<2,dim>() );
        
        if (this->derivative_flags[i] == capillary_pressure)
        {
            this->psd_pointer->set_capillary_pressure(this->capillary_pressure_vector);
            this->psd_pointer->set_critical_radius();
            this->psd_pointer->set_saturation();
            
            std::vector<double> derivate_rlp;
            this->psd_pointer->get_derivative_relative_liquid_permeability(derivate_rlp);
            
            for (unsigned int i=0; i<this->capillary_pressure_vector.size(); ++i)
            {
                for (unsigned int j=0; j<dim; ++j)
                {
                    for (unsigned int k=0; k<dim; ++k)
                    {
                        if (j == k)
                            dk_l[i][j][k] = derivate_rlp[i];
                        else 
                            dk_l[i][j][k] = 0.0;
                    }
                }
            }
            
            derivate_r_lp[ this->derivative_flags[i] ] = dk_l;
        }
        
        else 
            derivate_r_lp[ this->derivative_flags[i] ] = dk_l;
    }
    
}
//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::pcapillary(std::vector<double>& pc) const
{
    Assert( this->s_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in DesignMPL::pcapillary.") );
    Assert( this->T_vector.is_initialized(), ExcMessage("Temperature values have not been set using set_temperature method in DesignMPL::pcapillary.") );
    
    pc.clear();
    pc.resize(this->s_vector.size(), 0.0);
    
    if (method_capillary_function == "Kumbur07-corrected")
    {
        for (unsigned int i=0; i<this->s_vector.size(); ++i)
        {
            if (this->s_vector[i] < 1.0e-4)
            {
                pc[i] = 10.0 * (std::pow((293.0/this->T_vector[i]), 6.0)) * (0.1247 - (1.78e-4*this->T_vector[i])) * kumbur_factor * 
                        ( PTFE_loading*(1.8158 - 0.0328*PTFE_loading + 9.1235*std::pow(1.0e-4,2.0) - 2.5089*std::pow(1.0e-4,3.0)) + 0.7126*std::log(1.0e-4) );
            }
            else
            {
                pc[i] = 10.0 * (std::pow((293.0/this->T_vector[i]), 6.0)) * (0.1247 - (1.78e-4*this->T_vector[i])) * kumbur_factor * 
                        ( PTFE_loading*(1.8158 - 0.0328*PTFE_loading + 9.1235*std::pow(this->s_vector[i],2.0) - 2.5089*std::pow(this->s_vector[i],3.0)) + 0.7126*std::log(this->s_vector[i]) );
            }
        }
    }
    
    else
        AssertThrow(false, ExcNotImplemented() );
}


//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::saturation_from_capillary_equation(std::vector<double>& saturation) const
{
    Assert( this->capillary_pressure_vector.is_initialized(), ExcMessage("Liquid water capillary values have not been set using set_capillary_pressure method in design_MPL::liquid_permeablity.") );
    
    this->psd_pointer->set_capillary_pressure( this->capillary_pressure_vector );
    this->psd_pointer->set_critical_radius();
    this->psd_pointer->set_saturation();
    this->psd_pointer->get_saturation(saturation);
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::derivative_saturation_from_capillary_equation_PSD(std::vector<double>& derivative_saturation) const
{
    Assert( this->capillary_pressure_vector.is_initialized(), ExcMessage("Liquid water capillary values have not been set using set_capillary_pressure method in design_MPL::liquid_permeablity.") );
    
    this->psd_pointer->set_capillary_pressure( this->capillary_pressure_vector );
    this->psd_pointer->set_critical_radius();
    this->psd_pointer->set_saturation();
    this->psd_pointer->get_derivative_saturation(derivative_saturation);
}


//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::dpcapillary_dsat(std::vector<double> & dpc_ds) const
{
    Assert( this->s_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in DesignMPL::dpcapillary_dsat.") );
    Assert( this->T_vector.is_initialized(), ExcMessage("Temperature values have not been set using set_temperature method in DesignMPL::dpcapillary_dsat.") );
    
    dpc_ds.clear();
    dpc_ds.resize(this->s_vector.size(), 0.0);
    
    if (method_capillary_function == "Kumbur07-corrected")
    {
        for (unsigned int i=0; i<this->s_vector.size(); ++i)
        {
            if (this->s_vector[i] < 1.0e-4)
            {
                dpc_ds[i] =  0.0; // dyne/cm^2
            }
            else
            {
                dpc_ds[i] =  10.0 * (std::pow((293.0/this->T_vector[i]), 6.0)) * (0.1247 - (1.78e-4 * this->T_vector[i])) * kumbur_factor * 
                             ( (PTFE_loading*(2.0*9.1235*this->s_vector[i] - 3.0*2.5089*this->s_vector[i]*this->s_vector[i])) + (0.7126/this->s_vector[i]) ); // dyne/cm^2
            }
        }
    }
    
    else
        AssertThrow( false, ExcNotImplemented() );
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::derivative_dpcapillary_dsat(std::map< VariableNames, std::vector<double> > & d2pc_dsdu) const
{
    Assert( this->s_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in DesignMPL::dpcapillary_dsat.") );
    Assert( this->T_vector.is_initialized(), ExcMessage("Temperature values have not been set using set_temperature method in DesignMPL::derivative_dpcapillary_dsat.") );
    Assert( this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called in DesignMPL::derivative_dpcapillary_dsat.") );
    
    for (unsigned int i=0; i<this->derivative_flags.size(); ++i)
    {
        std::vector<double> d2pc(this->s_vector.size(), 0.0);
        
        if (this->derivative_flags[i] == liquid_water_saturation)
        {
            if (method_capillary_function == "Kumbur07-corrected")
            {
                for (unsigned int j=0; j<this->s_vector.size(); ++j)
                {
                    if (this->s_vector[j] < 1.0e-4)
                    {
                        d2pc[j] =  0.0; // dyne/cm^2                        
                    }
                    else
                    {
                        d2pc[j] =  10.0 * (std::pow((293.0/this->T_vector[j]), 6.0)) * (0.1247 - (1.78e-4 * this->T_vector[j])) * kumbur_factor * 
                             ( (PTFE_loading*(2.0*9.1235 - 6.0*2.5089*this->s_vector[j])) + ((-0.7126)/(this->s_vector[j]*this->s_vector[j])) ); // dyne/cm^2
                    }
                }
            }
            
            else
                AssertThrow( false, ExcNotImplemented() );
            
            //---------------
            d2pc_dsdu[ this->derivative_flags[i] ] = d2pc;
        }
        
        else if (this->derivative_flags[i] == temperature_of_REV)
        {
            if (method_capillary_function == "Kumbur07-corrected")
            {
                for (unsigned int j=0; j<this->s_vector.size(); ++j)
                {
                    if (this->s_vector[j] < 1.0e-4)
                    {
                        d2pc[j] = 0.0; // dyne/(cm^2-K)        
                    }
                    else
                    {
                        d2pc[j] = 10.0*kumbur_factor*std::pow(293.0, 6.0)*( PTFE_loading*(2.0*9.1235*this->s_vector[j]-3.0*2.5089*std::pow(this->s_vector[j],2.0))+(0.7126/this->s_vector[j]) ) * 
                                  ( (((-6.0)*(0.1247-(1.78e-4*this->T_vector[j])))/std::pow(this->T_vector[j],7.0)) + ((-1.78e-4)/std::pow(this->T_vector[j],6.0)) ); // dyne/(cm^2-K)
                    }
                }
            }
            
            else
                AssertThrow( false, ExcNotImplemented() );
            
            //---------------
            d2pc_dsdu[ this->derivative_flags[i] ] = d2pc;
        }
        
        else
            d2pc_dsdu[ this->derivative_flags[i] ] = d2pc;
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::interfacial_surface_area(std::vector<double>& a_lv) const
{
    Assert( this->s_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in DesignMPL::interfacial_surface_area.") );
    
    a_lv.clear();
    a_lv.resize( this->s_vector.size(), 0.0 );
    
    for (unsigned int i=0; i < this->s_vector.size(); ++i)
    {
        if ( (this->s_vector[i] > 0.0) && (this->s_vector[i] < 1.0) )
            a_lv[i] = 90000.0 * (std::pow(this->s_vector[i], 5.0)) * (std::pow((1.0-this->s_vector[i]),1.25));  // cm^2/cm^3
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::derivative_interfacial_surface_area(std::map< VariableNames, std::vector<double> >& deriv_a_lv) const
{
    Assert( this->s_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in DesignMPL::derivative_interfacial_surface_area.") );
    Assert( this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called in DesignMPL::derivative_interfacial_surface_area.") );
    
    for (unsigned int i=0; i < this->derivative_flags.size(); ++i)
    {
        std::vector<double> da_lv(this->s_vector.size(), 0.0);
        
        if (this->derivative_flags[i] == liquid_water_saturation)
        {
            for (unsigned int j=0; j<this->s_vector.size(); ++j)
            {
                if ( (this->s_vector[j] > 0.0) && (this->s_vector[j] < 1.0) )
                {
                    da_lv[j] = 90000.0 * ( (5.0*(std::pow(this->s_vector[j],4.0))*(std::pow((1.0-this->s_vector[j]),1.25))) - 
                                           ((std::pow(this->s_vector[j],5.0))*1.25*(std::pow((1.0-this->s_vector[j]),0.25))) );
                }
            }
            
            //---------------
            deriv_a_lv[ this->derivative_flags[i] ] = da_lv;
        }
        
        else
            deriv_a_lv[ this->derivative_flags[i] ] = da_lv;
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::interfacial_surface_area_PSD(std::vector<double>& a_lv) const
{
    Assert( this->capillary_pressure_vector.is_initialized(), ExcMessage("Liquid water capillary pressure values have not been set using set_saturation method in MPL::interfacial_surface_area.") );
    a_lv.clear();
    a_lv.resize( this->capillary_pressure_vector.size(), 0.0 );
    
    this->psd_pointer->set_capillary_pressure( this->capillary_pressure_vector );
    
    this->psd_pointer->set_critical_radius();
    this->psd_pointer->set_saturation();
    this->psd_pointer->get_liquid_gas_interfacial_surface(a_lv);
    
    
}
//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::derivative_interfacial_surface_area_PSD(std::vector<double>& deriv_a_lv) const
{
    Assert( this->capillary_pressure_vector.is_initialized(), ExcMessage("Liquid water capillary pressure values have not been set using set_saturation method in MPL::interfacial_surface_area.") );
    deriv_a_lv.clear();
    deriv_a_lv.resize( this->capillary_pressure_vector.size(), 0.0 );
    
    this->psd_pointer->set_capillary_pressure( this->capillary_pressure_vector );
    this->psd_pointer->set_critical_radius();
    this->psd_pointer->set_saturation();
    this->psd_pointer->get_derivative_liquid_gas_interfacial_surface(deriv_a_lv);

    
}
//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignMPL<dim>::derivative_interfacial_surface_area_PSD(std::map< VariableNames, std::vector<double> >& deriv_a_lv) const
{
    Assert( this->capillary_pressure_vector.is_initialized(), ExcMessage("Liquid water capillary pressure values have not been set using set_saturation method in MPL::derivative_interfacial_surface_area.") );
    Assert( this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called in MPL::derivative_interfacial_surface_area_PSD.") );
    
    for (unsigned int i=0; i < this->derivative_flags.size(); ++i)
    {
        std::vector<double> da_lv(this->capillary_pressure_vector.size(), 0.0);
        
        if (this->derivative_flags[i] == capillary_pressure)
        {
            std::vector<double> da_lv(this->capillary_pressure_vector.size(), 0.0);
            
            this->psd_pointer->set_capillary_pressure( this->capillary_pressure_vector );
            
            this->psd_pointer->set_critical_radius();
            this->psd_pointer->set_saturation();
            this->psd_pointer->get_derivative_liquid_gas_interfacial_surface(da_lv);  
            
            deriv_a_lv[ this->derivative_flags[i] ] = da_lv;
        }
        
        else
            deriv_a_lv[ this->derivative_flags[i] ] = da_lv;
    }
    
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations. 
template class NAME::DesignMPL<deal_II_dimension>;
