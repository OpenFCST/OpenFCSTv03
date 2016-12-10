//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: design_fibrous_GDL.cc
//    - Description: Class used to represent a fibrous GDL where effective properties are computed based on the porosity etc.
//    - Developers: M. Secanell (2009) and Madhur Bhaiya (2012)
//    - Id: $Id: design_fibrous_GDL.cc 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#include <layers/design_fibrous_GDL.h>

namespace NAME = FuelCellShop::Layer;

template <int dim>
const std::string NAME::DesignFibrousGDL<dim>::concrete_name ("DesignFibrousGDL");

template <int dim>
NAME::DesignFibrousGDL<dim> const* NAME::DesignFibrousGDL<dim>::PROTOTYPE = new NAME::DesignFibrousGDL<dim>();

//---------------------------------------------------------------------------
template <int dim>
NAME::DesignFibrousGDL<dim>::DesignFibrousGDL()
  : NAME::GasDiffusionLayer<dim> ()
{
    //FcstUtilities::log<<" Register DesignFibrousGDL GDL to FactoryMap"<<std::endl;
    this->get_mapFactory()->insert(std::pair<std::string, GasDiffusionLayer<dim>* > (this->concrete_name, this) ); 
}

//---------------------------------------------------------------------------
template <int dim>
NAME::DesignFibrousGDL<dim>::DesignFibrousGDL(const std::string& name)
  : NAME::GasDiffusionLayer<dim> (name)
{
}


//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignFibrousGDL<dim>::declare_parameters (const std::string& name, 
                                                 ParameterHandler &param) const
{            
    
    FuelCellShop::Layer::GasDiffusionLayer<dim>::declare_parameters(name, param);
    
    /// Nothing extra to declare since all properties are hard-coded.
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection(name);
        {
            param.enter_subsection(concrete_name); //-- Transport for the anisotropic case:
            {
                // Add generic data:
                //FuelCellShop::Layer::PorousLayer<dim>::declare_parameters(name,param);
                
                // General data
                param.declare_entry ("Porosity",  // volume fraction of void space in the GDL
                                     "0.6", // [-]
                                     Patterns::Double(0., 1.),
                                     "Volume fraction of void space in the GDL");
                
                param.declare_entry ("Method effective transport properties in pores",
                                     "Bruggemann",
                                     Patterns::Selection("Bruggemann|Percolation|Tomadakis|Mezedur"),
                                     "Method used to compute effective transport properties in the void phase."
                                     "Note the method defines the type of network, therefore the same method"
                                     "is used even in the GDL is anisotropic");
                                     
                param.declare_entry ("Method effective transport properties in solid",
                                     "Bruggemann",
                                     Patterns::Selection("Given|Bruggemann|Percolation"),
                                     "Method used to compute effective transport properties in solid"
                                     "Note the method defines the type of network, therefore the same method"
                                     "is used even in the GDL is anisotropic");
                                     
                param.declare_entry ("Method effective thermal conductivity",
                                     "Zamel",
                                     Patterns::Selection("Zamel|Given"),
                                     "Method used to compute effective thermal conductivity. Units [W/(cm-K)].");
                
                param.declare_entry ("Method relative liquid permeability",
                                     "Kumbur07",
                                     Patterns::Selection("Kumbur07|Wyllie"),
                                     "Method used to compute relative liquid permeability as a function of saturation.");
                param.declare_entry ("Irreducible liquid water saturation",
                                     "0.0",
                                     Patterns::Double(0., 1.),
                                     "Irreducible liquid water saturation in the GDL.");
                
                param.declare_entry ("Method capillary pressure - saturation function",
                                     "Kumbur07-corrected",
                                     Patterns::Selection("Kumbur07-corrected"),
                                     "Method used to compute capillary pressure as a function of saturation.");
                param.declare_entry ("Compaction pressure [MPa]",
                                     "0.0",
                                     Patterns::Double(0.),
                                     "Compaction pressure on the GDL, Units: [MPa].");
                param.declare_entry ("PTFE loading [% wt]",
                                     "5.0",
                                     Patterns::Double(5., 20.),
                                     "PTFE loading (% wt) in the GDL; Accepted range: 5 - 20 %.");
                
                //-- Anisotropic transport boolean flag:
                param.declare_entry ("Anisotropic transport",
                                     "false", // [S/cm]
                                     Patterns::Bool(),
                                     "Boolean variable. Set to true if we want to account for anisotropy of the GDL");
                
                // ---------------------------------------------------------------------------------------------------------
                //-- Transport for the isotropic case:
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
                
                param.declare_entry ("Solid network threshold",
                                     "0.12", // [-]
                                     Patterns::Double(),
                                     "Threshold value of the volume fraction of solid in the GDL."
                                     "If the solid phase is less than this value transport in the solid network does not occur");
                param.declare_entry ("Solid network constant",
                                     "2.0",
                                     Patterns::Double(),
                                     "Parameter used when using percolation theory");
                
                param.declare_entry ("Electrical conductivity [S/cm]",
                                     "100.", // [S/cm]
                                     Patterns::Double(0.),
                                     "Effective cond. if given is used, otherwise conducitivity of the raw material. Units [S/cm]");
                
                param.declare_entry ("Thermal conductivity [W/(cm-K)]",
                                     "0.01", // [W/(cm-K)]
                                     Patterns::Double(0.),
                                     "Effective thermal conductivity of the layer. Units [W/(cm-K)]");
                
                param.declare_entry ("Absolute permeability [cm^2]",
                                     "1.8e-7", // cm^2
                                     Patterns::Double(0.),
                                     "Absolute permeability of the layer, Units [cm^2]");                                     
                
                // ------------------------------------------------------------------------------------------------------------
                //-- Transport for the anisotropic case:
                //--- XX
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
                
                param.declare_entry ("Solid network threshold X",
                                     "0.12", // [-]
                                     Patterns::Double(),
                                     "Threshold value of the volume fraction of solid in the GDL."
                                     "If the solid phase is less than this value transport in the solid network does not occur");
                param.declare_entry ("Solid network constant X",
                                     "2.0",
                                     Patterns::Double(),
                                     "Parameter used when using percolation theory");
                
                param.declare_entry ("Electrical conductivity X [S/cm]",
                                     "100.", // [S/cm]
                                     Patterns::Double(0.),
                                     "X component - Effective cond. if given is used, otherwise conducitivity of the raw material. Units [S/cm]");
                
                param.declare_entry ("Thermal conductivity X [W/(cm-K)]",
                                     "0.01", // [W/(cm-K)]
                                     Patterns::Double(0.),
                                     "Component X of the thermal conductivity tensor. Units [W/(cm-K)]");
                
                param.declare_entry ("Absolute permeability X [cm^2]",
                                     "1.8e-7", // cm^2
                                     Patterns::Double(0.),
                                     "X component - Absolute permeability of the layer, Units [cm^2]");
                
                // YY
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
                
                param.declare_entry ("Solid network threshold Y",
                                     "0.12", // [-]
                                     Patterns::Double(),
                                     "Threshold value of the volume fraction of solid in the GDL."
                                     "If the solid phase is less than this value transport in the solid network does not occur");
                param.declare_entry ("Solid network constant Y",
                                     "2.0",
                                     Patterns::Double(),
                                     "Parameter used when using percolation theory");
                
                param.declare_entry ("Electrical conductivity Y [S/cm]",
                                     "100.", // [S/cm]
                                     Patterns::Double(0.),
                                     "Y component - Effective cond. if given is used, otherwise conducitivity of the raw material. Units [S/cm]");
                
                param.declare_entry ("Thermal conductivity Y [W/(cm-K)]",
                                     "0.03", // [W/(cm-K)]
                                     Patterns::Double(0.),
                                     "Component Y of the thermal conductivity tensor. Units [W/(cm-K)]");
                
                param.declare_entry ("Absolute permeability Y [cm^2]",
                                     "1.8e-7", // cm^2
                                     Patterns::Double(0.),
                                     "Y component - Absolute permeability of the layer, Units [cm^2]");
                
                // ZZ
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
                
                param.declare_entry ("Solid network threshold Z",
                                     "0.12", // [-]
                                     Patterns::Double(),
                                     "Threshold value of the volume fraction of solid in the GDL."
                                     "If the solid phase is less than this value transport in the solid network does not occur");
                param.declare_entry ("Solid network constant Z",
                                     "2.0",
                                     Patterns::Double(),
                                     "Parameter used when using percolation theory");
                
                param.declare_entry ("Electrical conductivity Z [S/cm]",
                                     "100.", // [S/cm]
                                     Patterns::Double(0.),
                                     "Z component - Effective cond. if given is used, otherwise conducitivity of the raw material. Units [S/cm]");
                
                param.declare_entry ("Thermal conductivity Z [W/(cm-K)]",
                                     "0.03", // [W/(cm-K)]
                                     Patterns::Double(0.),
                                     "Component Z of the thermal conductivity tensor. Units [W/(cm-K)]");
                
                param.declare_entry ("Absolute permeability Z [cm^2]",
                                     "1.8e-7", // cm^2
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
NAME::DesignFibrousGDL<dim>::initialize (ParameterHandler &param)
{
    NAME::GasDiffusionLayer<dim>::initialize(param);
    
    param.enter_subsection("Fuel cell data"); 
    {
        param.enter_subsection(this->name); 
        { 
            param.enter_subsection(concrete_name);
            {
                porosity = param.get_double("Porosity");
                solid_phase = 1.0 - porosity;
                method_eff_property_pores = param.get("Method effective transport properties in pores");
                method_eff_property_solid = param.get("Method effective transport properties in solid");
                method_eff_property_thermal = param.get("Method effective thermal conductivity");
                
                method_rel_liquid_permeability = param.get("Method relative liquid permeability");
                s_irr = param.get_double("Irreducible liquid water saturation");
                
                method_capillary_function = param.get("Method capillary pressure - saturation function");
                compaction_pressure = param.get_double("Compaction pressure [MPa]");
                PTFE_loading = param.get_double("PTFE loading [% wt]");
                
                // Anisotropy
                anisotropy = param.get_bool ("Anisotropic transport");
                
                // -----------
                porosity_th.resize(3);
                porosity_mu.resize(3);
                porosity_gamma.resize(3);
                
                solid_th.resize(3);
                solid_mu.resize(3);
                
                sigma_e.resize(3);
                
                k_thermal.resize(3);
                
                abs_permeability.resize(3);
                
                if (anisotropy == true)
                {
                    // X
                    porosity_th[0] = param.get_double("Porosity threshold X");
                    porosity_mu[0] = param.get_double("Porosity network constant X");
                    porosity_gamma[0] = param.get_double("Porosity gamma network constant X");
                    
                    solid_th[0] = param.get_double("Solid network threshold X");
                    solid_mu[0] = param.get_double("Solid network constant X");
                    
                    sigma_e[0] = param.get_double("Electrical conductivity X [S/cm]");
                    
                    k_thermal[0] = param.get_double("Thermal conductivity X [W/(cm-K)]");
                    
                    abs_permeability[0] = param.get_double("Absolute permeability X [cm^2]");
                    
                    // Y
                    porosity_th[1] = param.get_double("Porosity threshold Y");
                    porosity_mu[1] = param.get_double("Porosity network constant Y");
                    porosity_gamma[1] = param.get_double("Porosity gamma network constant Y");
                    
                    solid_th[1] = param.get_double("Solid network threshold Y");
                    solid_mu[1] = param.get_double("Solid network constant Y");
                    
                    sigma_e[1] = param.get_double("Electrical conductivity Y [S/cm]");
                    
                    k_thermal[1] = param.get_double("Thermal conductivity Y [W/(cm-K)]");
                    
                    abs_permeability[1] = param.get_double("Absolute permeability Y [cm^2]");
                    
                    // Z
                    porosity_th[2] = param.get_double("Porosity threshold Z");
                    porosity_mu[2] = param.get_double("Porosity network constant Z");
                    porosity_gamma[2] = param.get_double("Porosity gamma network constant Z");
                    
                    solid_th[2] = param.get_double("Solid network threshold Z");
                    solid_mu[2] = param.get_double("Solid network constant Z");
                    
                    sigma_e[2] = param.get_double("Electrical conductivity Z [S/cm]");
                    
                    k_thermal[2] = param.get_double("Thermal conductivity Z [W/(cm-K)]");
                    
                    abs_permeability[2] = param.get_double("Absolute permeability Z [cm^2]");
                }
                else
                {
                    for (unsigned int i =0; i < porosity_th.size(); i++)
                    {
                        porosity_th[i] = param.get_double("Porosity threshold");
                        porosity_mu[i] = param.get_double("Porosity network constant");
                        porosity_gamma[i] = param.get_double("Porosity gamma network constant");
                        
                        solid_th[i] = param.get_double("Solid network threshold");
                        solid_mu[i] = param.get_double("Solid network constant");
                        
                        sigma_e[i] = param.get_double("Electrical conductivity [S/cm]");
                        
                        k_thermal[i] = param.get_double("Thermal conductivity [W/(cm-K)]");
                        
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
NAME::DesignFibrousGDL<dim>::effective_gas_diffusivity(const double& prop,
                                                       const double& saturation,
                                                       double& prop_eff) const
{
    // Asset if anisotropy is false
    AssertThrow(anisotropy == false, ExcMessage("DesignFibrousGDL::effective_transport_property_pores with double return value can only be used for isotropic GDLs"));
    Assert( saturation >= 0.0, ExcMessage("Non-negative saturation values must be passed as an argument in DesignFibrousGDL::effective_gas_diffusivity.") );
    
    double effective_porosity = porosity * (1.0 - saturation);
    
    if (method_eff_property_pores == "Bruggemann")
    {
        // Most common method
        prop_eff = prop * std::pow(effective_porosity, 1.5);
    }
    
    else if (method_eff_property_pores == "Percolation")
    {
        double aux = (std::fabs(effective_porosity - porosity_th[0])) / (1.0 - porosity_th[0]);
        double step = 1.0e-3;
        
        if (effective_porosity >= porosity_th[0])
                step = 1.0;
        
        prop_eff = prop * ( std::pow(aux, porosity_mu[0]) * step );
    }
    
    else if (method_eff_property_pores == "Tomadakis")
    {
        // Method  reported by Tomadakis and Sotrichos; Correction for saturation proposed by Nam and Kaviany, 2003
        double aux = (std::fabs(porosity - porosity_th[0])) / (1.0 - porosity_th[0]);
        double step = 1.0e-3;
        
        if (porosity >= porosity_th[0])
            step = 1.0;
        
        prop_eff = prop * porosity * (std::pow(aux, porosity_mu[0])) * step * (std::pow((1.0 - saturation), 2.0));
    }
    
    else if (method_eff_property_pores == "Mezedur")
    {
        // Method reported in Mezedur et al., "Effects of pore structure, randomness
        // and size on effective mass diffusivity", AIChE J. 48 (2002) pp15-24
        prop_eff = prop * (1.0 - (std::pow( (1.0 - effective_porosity), 0.46)));
    }
    
    else
        AssertThrow( false, ExcNotImplemented() );
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignFibrousGDL<dim>::effective_gas_diffusivity(const double& prop,
                                                       const double& saturation,
                                                       Tensor<2,dim>& prop_eff) const
{
    Assert( saturation >= 0.0, ExcMessage("Non-negative saturation values must be passed as an argument in DesignFibrousGDL::effective_gas_diffusivity.") );
    
    prop_eff.clear(); // Initialize all entries to zero in the Tensor
    
    double effective_porosity = porosity * (1.0 - saturation);
    
    // Most common method
    if (method_eff_property_pores == "Bruggemann")
    {
        for (unsigned int i=0; i<dim; i++)
            prop_eff[i][i] = prop * std::pow(effective_porosity, 1.5);
    }

    else if (method_eff_property_pores == "Percolation")
    {
        for (unsigned int i=0; i<dim; i++)
        {
            double aux = (std::fabs(effective_porosity - porosity_th[i])) / (1.0 - porosity_th[i]);
            double step = 1.0e-3;
            
            if (effective_porosity >= porosity_th[i])
                step = 1.0;
            
            prop_eff[i][i] = prop * ( std::pow(aux, porosity_mu[i]) * step );
        }
    }
    
    // Method  reported by Tomadakis and Sotrichos; Correction for saturation proposed by Nam and Kaviany, 2003
    else if (method_eff_property_pores == "Tomadakis")
    {
        for (unsigned int i=0; i<dim; i++)
        {
            double aux = (std::fabs(porosity - porosity_th[i])) / (1.0 - porosity_th[i]);
            double step = 1.0e-3;
            
            if (porosity >= porosity_th[i])
                step = 1.0;
            
            prop_eff[i][i] = prop * porosity * (std::pow(aux, porosity_mu[i])) * step * (std::pow((1.0 - saturation), 2.0));
        }
    }
    
    // Method reported in Mezedur et al., "Effects of pore structure, randomness
    // and size on effective mass diffusivity", AIChE J. 48 (2002) pp15-24
    else if (method_eff_property_pores == "Mezedur")
    {
        for (unsigned int i=0; i<dim; i++)
            prop_eff[i][i] = prop * (1.0 - (std::pow( (1.0 - effective_porosity), 0.46)));
    }
    
    else
        AssertThrow( false, ExcNotImplemented() );
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignFibrousGDL<dim>::effective_gas_diffusivity(std::vector< Tensor<2,dim> >& prop_eff_vec) const
{
    Assert ( this->D_bulk.size() != 0, ExcMessage("compute_gas_diffusion not called before DesignFibrousGDL::effective_gas_diffusivity") );
    
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
template <int dim>
void
NAME::DesignFibrousGDL<dim>::effective_gas_diffusivity(Table< 2, Tensor<2,dim> > &D_eff) const
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
template<int dim>
void
NAME::DesignFibrousGDL<dim>::derivative_effective_gas_diffusivity(std::map< VariableNames, std::vector< Tensor<2,dim> > >& dprop_eff) const
{
    Assert( this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called in DesignFibrousGDL::derivative_effective_gas_diffusivity."));
    Assert(((this->dD_bulk_dT.size()!=0) && (this->D_bulk.size()!=0)), ExcMessage("compute_gas_diffusion not called before DesignFibrousGDL::derivative_effective_gas_diffusivity."));
    
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
            
            else if (method_eff_property_pores == "Tomadakis")
            {
                for (unsigned int j=0; j < this->D_bulk.size(); ++j)
                {
                    if (this->s_vector[j] >= 0.0)
                    {
                        for (unsigned int k=0; k<dim; ++k)
                        {
                            if ( porosity >= porosity_th[k] )
                            {
                                Dprop[j][k][k] = this->D_bulk[j] * porosity * (std::pow( ((porosity-porosity_th[k])/(1.0-porosity_th[k])), porosity_mu[k] )) * 
                                                 (-2.0) * (1.0 - this->s_vector[j]);
                            }
                        }
                    }
                }
            }
    
            else if (method_eff_property_pores == "Mezedur")
            {
                for (unsigned int j=0; j < this->D_bulk.size(); ++j)
                {
                    if (this->s_vector[j] >= 0.0)
                    {
                        for (unsigned int k=0; k<dim; ++k)
                            Dprop[j][k][k] = this->D_bulk[j] * (-0.46) * porosity * (std::pow( (1.0 - (porosity*(1.0 - this->s_vector[j]))), (-0.54) ));
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
                            Dprop[j][k][k] = this->D_bulk[j] * ds_dp[j] *(-1.5) * (std::pow(porosity,1.5)) * (std::pow((1.0 -  saturation[j]),0.5));
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
                                Dprop[j][k][k] = this->D_bulk[j] * ds_dp[j] *(-1.0) * porosity * porosity_mu[k] * (std::pow( (1.0-porosity_th[k]), ((-1.0)*(porosity_mu[k])) )) * 
                                (std::pow( ((porosity*(1.0-saturation[j])) - porosity_th[k]), (porosity_mu[k]-1.0) ));
                            }
                        }
                    }
                }
            }
            
            else if (method_eff_property_pores == "Tomadakis")
            {
                for (unsigned int j=0; j < this->D_bulk.size(); ++j)
                {
                    if (saturation[j] >= 0.0)
                    {
                        for (unsigned int k=0; k<dim; ++k)
                        {
                            if ( porosity >= porosity_th[k] )
                            {
                                Dprop[j][k][k] = this->D_bulk[j] * ds_dp[j] *porosity * (std::pow( ((porosity-porosity_th[k])/(1.0-porosity_th[k])), porosity_mu[k] )) * 
                                (-2.0) * (1.0 - saturation[j]);
                            }
                        }
                    }
                }
            }
            
            else if (method_eff_property_pores == "Mezedur")
            {
                for (unsigned int j=0; j < this->D_bulk.size(); ++j)
                {
                    if (saturation[j] >= 0.0)
                    {
                        for (unsigned int k=0; k<dim; ++k)
                            Dprop[j][k][k] = this->D_bulk[j] * ds_dp[j] *(-0.46) * porosity * (std::pow( (1.0 - (porosity*(1.0 - saturation[j]))), (-0.54) ));
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
NAME::DesignFibrousGDL<dim>::effective_electron_conductivity(double& prop_eff) const
{
    // Asset if anisotropy is false
    AssertThrow(anisotropy == false, ExcMessage("FibrousGDL::effective_transport_property_pores with double return value"
    " can only be used for isotropic GDLs"));
    // Initialize variable
    prop_eff = 0.0;
    const double prop = sigma_e[0];
    
    // In some papers prop_eff is give directly
    if (method_eff_property_solid == "Given")
    {
        prop_eff = prop;
    }
    // Most common method
    else if (method_eff_property_solid == "Bruggemann")
    {
        prop_eff = prop*pow(solid_phase, 1.5);
        
    }
    // Method used by M. Eikerling
    else if (method_eff_property_solid == "Percolation")
    {
        double aux = fabs(solid_phase - solid_th[0])/(1 - solid_th[0]);
        double step = 0.0;
        if (solid_phase >= solid_th[0])
            step = 1.0;
        prop_eff = prop*pow(aux,solid_mu[0])*step;
    }
    else
        AssertThrow( false, ExcNotImplemented() );
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::DesignFibrousGDL<dim>::effective_electron_conductivity(Tensor<2,dim>& prop_eff) const
{
    // Initialize variable
    prop_eff.clear(); // Initialize all entries to zero in the Tensor
    
    const std::vector<double> prop = sigma_e;
    Assert(prop.size() == 3, ExcMessage("Prop should contain the conductivity in the X, Y and Z dimension"));
    
    // In some papers prop_eff is give directly
    if (method_eff_property_solid == "Given")
    {
        for (unsigned int i=0; i<dim; i++)
            prop_eff[i][i] = prop[i];
    }
    // Most common method
    else if (method_eff_property_solid == "Bruggemann")
    {
        for (unsigned int i=0; i<dim; i++)
            prop_eff[i][i] = prop[i]*pow(solid_phase, 1.5);
    }
    // Method used by M. Eikerling
    else if (method_eff_property_solid == "Percolation")
    {
        for (unsigned int i=0; i<dim; i++)
        {
            double aux = fabs(solid_phase - solid_th[i])/(1 - solid_th[i]);
            double step = 0.0;
            if (solid_phase >= solid_th[i])
                step = 1.0;
            prop_eff[i][i] = prop[i]*pow(aux,solid_mu[i])*step;
        }
    }
    else
        AssertThrow( false, ExcNotImplemented() );
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignFibrousGDL<dim>::effective_thermal_conductivity(std::vector< Tensor<2,dim> >& prop_eff) const
{
    Assert (this->T_vector.size() != 0, ExcMessage("Temperature values have not been set. set_temperature has not been probably called."));
    prop_eff.clear();
    prop_eff.resize(this->T_vector.size());
    
    if (method_eff_property_thermal == "Zamel")
        // Ref: N. Zamel, E. Litovsky, X. Li, and J. Kleiman. Measurement of the through-plane thermal conductivity of carbon paper diffusion media for the temperature range from -50 to 120 C. International Journal of Hydrogen Energy, 36(19):12618-12625, 2011.
    {
        for ( unsigned int i = 0; i<this->T_vector.size(); ++i)
        {
            double Temp = this->T_vector[i] - 273.15;		// Converting into degree Celsius
            double K_in = (-7.166e-6)*std::pow(Temp,3.0) + (2.24e-3)*std::pow(Temp,2.0) - (0.237)*Temp + 20.1;		// In-plane thermal conductivity in W/m-C
            // Heat Barrier Resistance Coefficient, M
            double M = (-1.495e-11)*std::pow(Temp,5.0) + (2.601e-9)*std::pow(Temp,4.0) - (6.116e-8)*std::pow(Temp,3.0) - (9.829e-6)*std::pow(Temp,2.0) + (8.754e-4)*Temp + 0.0664;
            
            for ( unsigned int j = 0; j<dim; ++j)
            {
                // Factor of 0.01 is used to convert W/m-C into W/cm-K
                if ( j == 0)
                { prop_eff[i][j][j] = M*K_in*0.01;}
                else
                { prop_eff[i][j][j] = K_in*0.01;}
            }
        }
    }
    
    else if (method_eff_property_thermal == "Given")
    {
        for (unsigned int i = 0; i < prop_eff.size(); ++i)
        {
            for (unsigned int j = 0; j < dim; ++j)
                prop_eff[i][j][j] = k_thermal[j];
        }
    }
    
    else
    {
        FcstUtilities::log << "Unknown method to compute temperature dependent effective thermal conductivity in "<<__FILE__ <<" line "<<__LINE__<<std::endl;
        abort();
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignFibrousGDL<dim>::derivative_effective_thermal_conductivity(std::vector< Tensor<2,dim> >& dK_dT) const
{
    Assert (this->T_vector.size() != 0, ExcMessage("Temperature values have not been set. set_temperature has not been probably called."));
    dK_dT.clear();
    dK_dT.resize(this->T_vector.size());
    
    if (method_eff_property_thermal == "Zamel")
    // Ref: N. Zamel, E. Litovsky, X. Li, and J. Kleiman. Measurement of the through-plane thermal conductivity of carbon paper diffusion media for the temperature range from -50 to 120 C. International Journal of Hydrogen Energy, 36(19):12618-12625, 2011.
    {
        for ( unsigned int i = 0; i<this->T_vector.size(); ++i)
        {
            double Temp = this->T_vector[i] - 273.15;   // Converting into degree Celsius
            double K_in = (-7.166e-6)*std::pow(Temp,3.0) + (2.24e-3)*std::pow(Temp,2.0) - (0.237)*Temp + 20.1;		// In-plane thermal conductivity in W/m-C
            // Heat Barrier Resistance Coefficient, M
            double M = (-1.495e-11)*std::pow(Temp,5.0) + (2.601e-9)*std::pow(Temp,4.0) - (6.116e-8)*std::pow(Temp,3.0) - (9.829e-6)*std::pow(Temp,2.0) + (8.754e-4)*Temp + 0.0664;
            
            double dKin_dT = (-7.166e-6)*(std::pow(Temp,2))*3.0 + (2.24e-3)*Temp*2.0 - 0.237;
            double dM_dT = (-1.495e-11)*(std::pow(Temp,4))*5.0 + (2.601e-9)*(std::pow(Temp,3.0))*4.0 - (6.116e-8)*(std::pow(Temp,2.0))*3.0 - (9.829e-6)*Temp*2.0 + (8.754e-4);
            
            for ( unsigned int j = 0; j<dim; ++j)
            {
                // Factor of 0.01 is used to convert W/m-c into W/cm-K
                if ( j == 0)
                { dK_dT[i][j][j] = 0.01*((M*dKin_dT) + (dM_dT*K_in));}
                else
                { dK_dT[i][j][j] = 0.01*dKin_dT;}
            }
        }
    }
    
    else if (method_eff_property_thermal == "Given")
    {
        for (unsigned int i=0; i<dK_dT.size(); ++i)
        {
            for (unsigned int j=0; j<dim; ++j)
                dK_dT[i][j][j] = 0.0;
        }
    }
    
    else
    {
        FcstUtilities::log << "Unknown method to compute temperature dependent effective thermal conductivity in "<<__FILE__ <<" line "<<__LINE__<<std::endl;
        abort();
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignFibrousGDL<dim>::liquid_permeablity(std::vector< Tensor<2,dim> >& k_l ) const
{
    Assert( this->s_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in DesignFibrousGDL::liquid_permeablity.") );
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
NAME::DesignFibrousGDL<dim>::derivative_liquid_permeablity(std::map< VariableNames, std::vector< Tensor<2,dim> > >& deriv_k_l) const
{
    Assert( this->s_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in DesignFibrousGDL::derivative_liquid_permeablity.") );
    Assert (this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called in DesignFibrousGDL::derivative_liquid_permeablity."));
    
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
NAME::DesignFibrousGDL<dim>::saturated_liquid_permeablity_PSD(double & r_lp ) const
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
NAME::DesignFibrousGDL<dim>::relative_liquid_permeability_PSD(std::vector< Tensor<2,dim> >& r_lp ) const
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
NAME::DesignFibrousGDL<dim>::derivative_relative_liquid_permeablity_PSD(std::vector<double> & derivate_r_lp ) const
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
NAME::DesignFibrousGDL<dim>::derivative_relative_liquid_permeablity_PSD(std::map< VariableNames, std::vector< Tensor<2,dim> > >& derivate_r_lp ) const
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
NAME::DesignFibrousGDL<dim>::pcapillary(std::vector<double>& pc) const
{
    Assert( this->s_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in DesignFibrousGDL::pcapillary.") );
    Assert( this->T_vector.is_initialized(), ExcMessage("Temperature values have not been set using set_temperature method in DesignFibrousGDL::pcapillary.") );
    
    pc.clear();
    pc.resize(this->s_vector.size(), 0.0);
    
    if (method_capillary_function == "Kumbur07-corrected")
    {
        for (unsigned int i=0; i<this->s_vector.size(); ++i)
        {
            if (this->s_vector[i] < 1.0e-4)
            {
                pc[i] = 10.0 * (std::pow((293.0/this->T_vector[i]), 6.0)) * (0.1247 - (1.78e-4*this->T_vector[i])) * kumbur_factor * 
                        ( PTFE_loading*(0.0444 - 0.0014*PTFE_loading - 0.0275*std::pow(1.0e-4,2.0) + 0.0769*std::pow(1.0e-4,3.0)) + 0.0564*std::log(1.0e-4) );
            }
            else
            {
                pc[i] = 10.0 * (std::pow((293.0/this->T_vector[i]), 6.0)) * (0.1247 - (1.78e-4*this->T_vector[i])) * kumbur_factor * 
                        ( PTFE_loading*(0.0444 - 0.0014*PTFE_loading - 0.0275*std::pow(this->s_vector[i],2.0) + 0.0769*std::pow(this->s_vector[i],3.0)) + 0.0564*std::log(this->s_vector[i]) );
            }
        }
    }
    
    else
        AssertThrow(false, ExcNotImplemented() );
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignFibrousGDL<dim>::saturation_from_capillary_equation(std::vector<double>& saturation) const
{
    Assert( this->capillary_pressure_vector.is_initialized(), ExcMessage("Liquid water capillary values have not been set using set_capillary_pressure method in design_fibrous_GDL::liquid_permeablity.") );
    
    this->psd_pointer->set_capillary_pressure( this->capillary_pressure_vector );
    this->psd_pointer->set_critical_radius();
    this->psd_pointer->set_saturation();
    this->psd_pointer->get_saturation(saturation);
}


//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignFibrousGDL<dim>::derivative_saturation_from_capillary_equation_PSD(std::vector<double>& derivative_saturation) const
{
    Assert( this->capillary_pressure_vector.is_initialized(), ExcMessage("Liquid water capillary values have not been set using set_capillary_pressure method in design_fibrous_GDL::liquid_permeablity.") );
    
    this->psd_pointer->set_capillary_pressure( this->capillary_pressure_vector );
    this->psd_pointer->set_critical_radius();
    this->psd_pointer->set_saturation();
    this->psd_pointer->get_derivative_saturation(derivative_saturation);
}
//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignFibrousGDL<dim>::dpcapillary_dsat(std::vector<double> & dpc_ds) const
{
    Assert( this->s_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in DesignFibrousGDL::dpcapillary_dsat.") );
    Assert( this->T_vector.is_initialized(), ExcMessage("Temperature values have not been set using set_temperature method in DesignFibrousGDL::dpcapillary_dsat.") );
    
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
                             ( (PTFE_loading*(3.0*0.0769*this->s_vector[i]*this->s_vector[i] - 2.0*0.0275*this->s_vector[i])) + (0.0564/this->s_vector[i]) ); // dyne/cm^2
            }
        }
    }
    
    else
        AssertThrow( false, ExcNotImplemented() );
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignFibrousGDL<dim>::derivative_dpcapillary_dsat(std::map< VariableNames, std::vector<double> > & d2pc_dsdu) const
{
    Assert( this->s_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in DesignFibrousGDL::dpcapillary_dsat.") );
    Assert( this->T_vector.is_initialized(), ExcMessage("Temperature values have not been set using set_temperature method in DesignFibrousGDL::derivative_dpcapillary_dsat.") );
    Assert( this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called in DesignFibrousGDL::derivative_dpcapillary_dsat.") );
    
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
                             ( (PTFE_loading*(6.0*0.0769*this->s_vector[j] - 2.0*0.0275)) + ((-0.0564)/(this->s_vector[j]*this->s_vector[j])) ); // dyne/cm^2
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
                        d2pc[j] = 10.0*kumbur_factor*std::pow(293.0, 6.0)*( PTFE_loading*(3.0*0.0769*std::pow(this->s_vector[j],2.0)-2.0*0.0275*this->s_vector[j])+(0.0564/this->s_vector[j]) ) * 
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
NAME::DesignFibrousGDL<dim>::interfacial_surface_area(std::vector<double>& a_lv) const
{
    Assert( this->s_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in DesignFibrousGDL::interfacial_surface_area.") );
    
    a_lv.clear();
    a_lv.resize( this->s_vector.size(), 0.0 );
    
    for (unsigned int i=0; i < this->s_vector.size(); ++i)
    {
        if ( (this->s_vector[i] > 0.0) && (this->s_vector[i] < 1.0) )
            a_lv[i] = 140.0 * (std::pow(this->s_vector[i], 1.15)) * (std::pow((1.0-this->s_vector[i]),1.75));  // cm^2/cm^3
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DesignFibrousGDL<dim>::derivative_interfacial_surface_area(std::map< VariableNames, std::vector<double> >& deriv_a_lv) const
{
    Assert( this->s_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in DesignFibrousGDL::derivative_interfacial_surface_area.") );
    Assert( this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called in DesignFibrousGDL::derivative_interfacial_surface_area.") );
    
    for (unsigned int i=0; i < this->derivative_flags.size(); ++i)
    {
        std::vector<double> da_lv(this->s_vector.size(), 0.0);
        
        if (this->derivative_flags[i] == liquid_water_saturation)
        {
            for (unsigned int j=0; j<this->s_vector.size(); ++j)
            {
                if ( (this->s_vector[j] > 0.0) && (this->s_vector[j] < 1.0) )
                {
                    da_lv[j] = 140.0 * ( (1.15*(std::pow(this->s_vector[j],0.15))*(std::pow((1.0-this->s_vector[j]),1.75))) - 
                                           ((std::pow(this->s_vector[j],1.15))*1.75*(std::pow((1.0-this->s_vector[j]),0.75))) );
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
NAME::DesignFibrousGDL<dim>::interfacial_surface_area_PSD(std::vector<double>&  a_lv)  const
{
    Assert( this->capillary_pressure_vector.is_initialized(), ExcMessage("Liquid water capillary pressure values have not been set using set_saturation method in GDL::interfacial_surface_area.") );
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
NAME::DesignFibrousGDL<dim>::derivative_interfacial_surface_area_PSD(std::vector<double>& deriv_a_lv) const
{
    Assert( this->capillary_pressure_vector.is_initialized(), ExcMessage("Liquid water capillary pressure values have not been set using set_saturation method in GDL::interfacial_surface_area.") );
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
NAME::DesignFibrousGDL<dim>::derivative_interfacial_surface_area_PSD(std::map< VariableNames, std::vector<double> >& deriv_a_lv) const
{
    Assert( this->capillary_pressure_vector.is_initialized(), ExcMessage("Liquid water capillary pressure values have not been set using set_saturation method in ConventionalCL::derivative_interfacial_surface_area.") );
    Assert( this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called in ConventionalCL::derivative_interfacial_surface_area_PSD.") );
    
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
            
            //---------------
            deriv_a_lv[ this->derivative_flags[i] ] = da_lv;
        }
        
        else
            deriv_a_lv[ this->derivative_flags[i] ] = da_lv;
    }
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations. 
template class NAME::DesignFibrousGDL<deal_II_dimension>;
