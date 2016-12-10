//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: conventional_cl.cc
//    - Description: Class characterizing the conventional catalyst layer and methods for computing effective properties. It also provides interface to various material classes used in catalyst layer.
//    - Developers: Marc Secanell (2006-13), Peter Dobson (2011) and Madhur Bhaiya (2013)
//    - Id: $Id: conventional_CL.cc 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#include <layers/conventional_CL.h>
#include <microscale/agglomerate_base.h>

namespace NAME = FuelCellShop::Layer;

typedef std::map< VariableNames, std::vector<double> >::iterator dmap_iter; 

template <int dim>
const std::string NAME::ConventionalCL<dim>::concrete_name ("ConventionalCL");

//---------------------------------------------------------------------------
template <int dim>
NAME::ConventionalCL<dim>::ConventionalCL(std::string name)
  : CatalystLayer<dim>(name)
{}

//---------------------------------------------------------------------------
template <int dim>
NAME::ConventionalCL<dim>::ConventionalCL()
  : CatalystLayer<dim>()
{  }

//---------------------------------------------------------------------------
template <int dim>
NAME::ConventionalCL<dim>::~ConventionalCL()
{}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::declare_parameters (const std::string& cl_section_name,
                                               ParameterHandler &param) const
{
    NAME::CatalystLayer<dim>::declare_parameters(cl_section_name, param);
    
    param.enter_subsection("Fuel cell data"); 
    {
        param.enter_subsection(cl_section_name); 
        {
            param.enter_subsection(this->concrete_name); 
            {
                //-- Composition
                param.declare_entry ("Platinum loading on support (%wt)", // mass percentage of platinum catalyst on the support carbon black
                                    "", // [-]
                                    Patterns::Map( Patterns::Integer(0,255), Patterns::Double(0) ),
                                    "Mass percentage of platinum catalyst on the support carbon black eg: 2:0.20, 3:0.40, 4:0.467");
                param.declare_entry ("Platinum loading per unit volume (mg/cm3)",// catalyst platinum mass loading per unit volume of CL
                                    "",
                                    Patterns::Map( Patterns::Integer(0,255), Patterns::Double(0) ),
                                    "Catalyst platinum mass loading per unit volume of CL");

                param.declare_entry ("Method to compute porosity",
                                    "marc",
                                    Patterns::Selection("marc|NafionLoading|ICRatio"),
                                    "Method used to calculate the porosity based on the input parameter being either Electrolyte loading (%wt) or Ionomer Carbon Ratio");
                param.declare_entry ("Electrolyte loading (%wt)",  //
                                    "", // [-]
                                    Patterns::Map( Patterns::Integer(0,255), Patterns::Double(0) ),
                                    "Electrode loading is the weight percentage of ionomer per gram of CL eg: 2:0.20, 3:0.30, 4:0.40");
                param.declare_entry ("Ionomer to Carbon Ratio",  //
                                    "", // [-]
                                    Patterns::Map( Patterns::Integer(0,255), Patterns::Double(0) ),
                                    "grams of ionomer per gram of Carbon");
            
                //-- Network characteristics
                param.declare_entry ("Method effective transport properties in pores",
                                    "Bruggemann",
                                    Patterns::Selection("Given|Bruggemann|Percolation"),
                                    "Method used to compute effective transport properties in the void phase.");
                param.declare_entry ("Porosity threshold",  // volume fraction of void space in the CL
                                    "0.12", // [-]
                                    Patterns::Double(),
                                    "Threshold value of the volume fraction of void space in the CL."
                                    "If the porosity is less than this value transport does not occur");
                param.declare_entry ("Porosity network constant",
                                    "2.0",
                                    Patterns::Double(),
                                    "Parameter used when using percolation theory");
                param.declare_entry ("Porosity gamma network constant",
                                    "0.0",
                                    Patterns::Double(),
                                    "Parameter used when using percolation theory to account for extra diffusion");
                // ---
                param.declare_entry ("Method effective transport properties in solid phase",
                                    "Bruggemann",
                                    Patterns::Selection("Given|Bruggemann|Percolation"),
                                    "Method used to compute effective transport properties in pores");
                param.declare_entry ("Solid network threshold",  // volume fraction of solid in the CL
                                    "0.12", // [-]
                                    Patterns::Double(),
                                    "Threshold value of the volume fraction of solid (electron conductive) phase in the CL."
                                    "If the solid phase is less than this value transport in the fibre network does not occur");
                param.declare_entry ("Solid network constant",
                                    "2.0",
                                    Patterns::Double(),
                                    "Parameter used when using percolation theory");
                
                param.declare_entry ("Method effective transport properties in electrolyte phase",
                                    "Bruggemann",
                                    Patterns::Selection("Given|Bruggemann|Percolation|Iden11"),
                                    "Method used to compute effective transport properties in pores");
                param.declare_entry ("Electrolyte network threshold",  // volume fraction of electrolyte in the CL
                                    "0.12", // [-]
                                    Patterns::Double(),
                                    "Threshold value of the volume fraction of electrolyte (proton conductive) phase in the CL."
                                    "If the electrolyte phase is less than this value transport in the network does not occur");
                param.declare_entry ("Electrolyte network constant",
                                    "2.0",
                                    Patterns::Double(),
                                    "Parameter used when using percolation theory");
                
                param.declare_entry ("Method to compute active area",
                                    "given",
                                    Patterns::Selection("given|Marr|ETEK06|ETEK07"));
                param.declare_entry ("Active area [cm^2/cm^3]",
                                    "",
                                    Patterns::Map( Patterns::Integer(0,255), Patterns::Double(0) ),
                                    "Parameter used to represent how active the CL is. Area [cm^2] of Pt per volume of CL eg: 2:1.0e5, 3:2.0e5, 4:2.0e5");
            
                //----
                param.declare_entry ("Method effective thermal conductivity",
                                    "Given",
                                    Patterns::Selection("Given"),
                                    "Method used to compute effective thermal conductivity");
                param.declare_entry ("Thermal conductivity, [W/(cm K)]",
                                    "",
                                    Patterns::Map( Patterns::Integer(0,255), Patterns::Double(0) ),
                                    "Thermal conductivity of each CL layer eg: 2:0.015, 3:0.0027, 4:0.015");
                
                //----
                param.declare_entry ("Method relative liquid permeability",
                                    "Kumbur07",
                                    Patterns::Selection("Kumbur07|Wyllie"),
                                    "Method used to compute relative liquid permeability as a function of saturation.");
                param.declare_entry ("Irreducible liquid water saturation",
                                    "",
                                    Patterns::Map( Patterns::Integer(0,255), Patterns::Double(0., 1.) ),
                                    "Irreducible liquid water saturation in the CL.");
                param.declare_entry ("Absolute permeability [cm^2]",
                                    "", // cm^2
                                    Patterns::Map( Patterns::Integer(0,255), Patterns::Double(0) ),
                                    "Absolute permeability of the layer, Units [cm^2]");
                
                //----
                param.declare_entry ("Method capillary pressure - saturation function",
                                    "Ye07",
                                    Patterns::Selection("Ye07"),
                                    "Method used to compute capillary pressure as a function of saturation.");
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
NAME::ConventionalCL<dim>::initialize (ParameterHandler &param)
{
    NAME::CatalystLayer<dim>::initialize(param);
    
    rho_N = this->electrolyte->get_density();
    rho_c = this->catalyst_support->get_density();
    rho_Pt = this->catalyst->get_density();
    
    param.enter_subsection("Fuel cell data"); 
    {
        param.enter_subsection(this->name);
        {
            param.enter_subsection(this->concrete_name); 
            {
                prc_Pt = FcstUtilities::string_to_map<unsigned int, double>( Utilities::split_string_list( param.get("Platinum loading on support (%wt)") ) );
                V_Pt = FcstUtilities::string_to_map<unsigned int, double>( Utilities::split_string_list( param.get("Platinum loading per unit volume (mg/cm3)") ) );

                method_porosity = param.get("Method to compute porosity");

                loading_N = FcstUtilities::string_to_map<unsigned int, double>( Utilities::split_string_list( param.get("Electrolyte loading (%wt)") ) );
                IC_ratio = FcstUtilities::string_to_map<unsigned int, double>( Utilities::split_string_list( param.get("Ionomer to Carbon Ratio") ) );

                method_eff_property_pores = param.get("Method effective transport properties in pores");
                method_eff_property_solid = param.get("Method effective transport properties in solid phase");
                method_eff_property_electrolyte = param.get("Method effective transport properties in electrolyte phase");
                
                porosity_th = param.get_double("Porosity threshold");
                porosity_mu = param.get_double("Porosity network constant");
                porosity_gamma = param.get_double("Porosity gamma network constant");
                
                solid_th = param.get_double("Solid network threshold");
                solid_mu = param.get_double("Solid network constant");
                
                electrolyte_th = param.get_double("Electrolyte network threshold");
                electrolyte_mu = param.get_double("Electrolyte network constant");
                
                method_Av = param.get("Method to compute active area");
                
                // Initializing thermal conductivity value provided in parameter file
                method_eff_thermal = param.get("Method effective thermal conductivity");
                k_T = FcstUtilities::string_to_map<unsigned int, double>( Utilities::split_string_list( param.get("Thermal conductivity, [W/(cm K)]") ) );
                
                method_rel_liquid_permeability = param.get("Method relative liquid permeability");
                s_irr = FcstUtilities::string_to_map<unsigned int, double>( Utilities::split_string_list(param.get("Irreducible liquid water saturation") ) );
                abs_permeability = FcstUtilities::string_to_map<unsigned int, double>( Utilities::split_string_list( param.get("Absolute permeability [cm^2]") ) );
                
                method_capillary_function = param.get("Method capillary pressure - saturation function");
                
                // Some mandatory checks
                AssertThrow( this->material_ids.size() == prc_Pt.size(), ExcMessage("Sizes do not match for prc_Pt in ConventionalCL::initialize function.") );
                AssertThrow( this->material_ids.size() == V_Pt.size(), ExcMessage("Sizes do not match for V_Pt in ConventionalCL::initialize function.") );

                // Need to improve this here..
                //AssertThrow( this->material_ids.size() == k_T.size(), ExcMessage("Sizes do not match for k_T in ConventionalCL::initialize function.") );
                //AssertThrow( this->material_ids.size() == abs_permeability.size(), ExcMessage("Sizes do not match for abs_permeability in ConventionalCL::initialize function.") );

                for (unsigned int i=0; i<this->material_ids.size(); ++i)
                {
                    AssertThrow( prc_Pt.find(this->material_ids.at(i))!= prc_Pt.end(), ExcMessage("Material id not found in prc_Pt corresponding to define material ids in ConventionalCL.") );
                    AssertThrow( V_Pt.find(this->material_ids.at(i))!= V_Pt.end(), ExcMessage("Material id not found in V_Pt corresponding to define material ids in ConventionalCL.") );
                    // Need to improve this here..
                    //AssertThrow( k_T.find(this->material_ids.at(i))!= k_T.end(), ExcMessage("Material id not found in k_T corresponding to define material ids in ConventionalCL.") );
                    //AssertThrow( abs_permeability.find(this->material_ids.at(i))!= abs_permeability.end(), ExcMessage("Material id not found in abs_permeability corresponding to define material ids in ConventionalCL.") );
                }
                
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                if (method_porosity.compare("marc") == 0 || method_porosity.compare("NafionLoading") == 0)
                {
                    AssertThrow( this->material_ids.size()==loading_N.size(), ExcMessage("Sizes do not match for loading_N in ConventionalCL::initialize function.") );
                    IC_ratio.clear();
                    
                    for (unsigned int i=0; i<this->material_ids.size(); ++i)
                    {
                        AssertThrow( loading_N.find(this->material_ids.at(i))!= loading_N.end(), ExcMessage("Material id not found in loading_N corresponding to define material ids in ConventionalCL.") );
                        
                        unsigned int current_id = this->material_ids.at(i);
                        IC_ratio[current_id] = loading_N.at(current_id) * ( 1.0 + prc_Pt.at(current_id) ) / ( 1.0 - loading_N.at(current_id) );
                    }
                }
                
                else if (method_porosity.compare("ICRatio") == 0)
                {
                    AssertThrow( this->material_ids.size()==IC_ratio.size(), ExcMessage("Sizes do not match for IC_ratio in ConventionalCL::initialize function.") );
                    loading_N.clear();
                    
                    for (unsigned int i=0; i<this->material_ids.size(); ++i)
                    {
                        AssertThrow( IC_ratio.find(this->material_ids.at(i))!= IC_ratio.end(), ExcMessage("Material id not found in IC_ratio corresponding to define material ids in ConventionalCL.") );
                        
                        unsigned int current_id = this->material_ids.at(i);
                        loading_N[current_id] = IC_ratio.at(current_id) / (IC_ratio.at(current_id) + prc_Pt.at(current_id) + 1.0);
                    }
                }
                
                else
                {
                    FcstUtilities::log<<"Method to compute porosity not defined properly"<<std::endl;
                    abort();//throw(std::exception e);
                }
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                if(method_Av == "given")
                    Av = FcstUtilities::string_to_map<unsigned int, double>( Utilities::split_string_list( param.get("Active area [cm^2/cm^3]") ) );
                else
                    compute_Av();
                
                AssertThrow( this->material_ids.size()==Av.size(), ExcMessage("Sizes do not match for Av in ConventionalCL::initialize function.") );
                for (unsigned int i=0; i<this->material_ids.size(); ++i)
                    AssertThrow( Av.find(this->material_ids.at(i))!= Av.end(), ExcMessage("Material id not found in IC_ratio corresponding to define material ids in ConventionalCL.") );
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                
            }
            param.leave_subsection();
        }
        param.leave_subsection();   
    }
    param.leave_subsection();
    
    compute_volume_fraction();
        
    for (unsigned int i=0; i<this->material_ids.size(); ++i)
    {
        unsigned int current_id = this->material_ids.at(i);
        
        // Checking physical properties of cell. If volume fractions are unphysical code will not run.
        AssertThrow((epsilon_S.at(current_id) < 1.) or  (epsilon_V.at(current_id) > 0.), ExcMessage("ERROR: Adjust volume fractions to physical values.") );
        AssertThrow((epsilon_S.at(current_id) > 0.) or  (epsilon_V.at(current_id) < 1.), ExcMessage("ERROR: Adjust volume fractions to physical values.") );
        
        // Checking physical properties of cell. If using Percolation and void or Solid is lower than threshold, code will not run.
        if (method_eff_property_pores == "Percolation")
        {
            AssertThrow((epsilon_V.at(current_id) > porosity_th), ExcMessage("ERROR: Adjust volume fractions or CL thickness. Porosity below threshold.") );
            AssertThrow((epsilon_S.at(current_id) > solid_th), ExcMessage("ERROR: Adjust volume fractions or CL thickness. Solid below threshold.") );
        }
        
        /*
        // Previously physical checks corrected and read out a warning. This maybe useful when running an optimization.
        //Sometimes bounds can be slightly overshoot depending on the optimization algorithm.
         if (epsilon_S.at(current_id) > 1. ||  epsilon_V.at(current_id) < 0.) {
             epsilon_V.at(current_id) = 1e-4;
             epsilon_S.at(current_id) = 1. - epsilon_N.at(current_id) - epsilon_V.at(current_id);
             FcstUtilities::log << "WARNING: Adjusted volume fractions to physical values" << std::endl;
             this->print_layer_properties();//this->print_volume_fractions();
        }
        else if (epsilon_S.at(current_id) < 0. ||  epsilon_V.at(current_id) > 1.) {
             epsilon_S.at(current_id) = 1e-4;
             epsilon_V.at(current_id) = 1. - epsilon_N - epsilon_S;
             FcstUtilities::log << "WARNING: Adjusted volume fractions to physical values" << std::endl;
             this->print_layer_properties();//this->print_volume_fractions();
         }
         */
        
        // Apply new porosity to PSD
        this->porosity = epsilon_V.at(current_id);
        if(this->PSD_is_used == true) 
            this->psd_pointer->set_porosity (this->porosity);
    }
    
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::effective_gas_diffusivity(const double& prop,
                                                     const double& saturation,
                                                     double& prop_eff) const
{
    Assert( saturation >= 0.0, ExcMessage("Non-negative saturation values must be passed as an argument in ConventionalCL::effective_gas_diffusivity.") );
    
    // Initialize variable
    prop_eff = 0.0;
    
    double effective_porosity = epsilon_V.at(this->local_material_id()) * (1.0 - saturation);
    
    if (method_eff_property_pores == "Bruggemann")
    {
        // Most common method
        prop_eff = prop * std::pow(effective_porosity, 1.5);
    }
    
    else if (method_eff_property_pores == "Percolation")
    {
        // Method used by M. Eikerling 
        // Note: For porosity_th[0] I use index 0 for convenience. They all contain the same value
        double aux = (std::fabs(effective_porosity - porosity_th)) / (1.0 - porosity_th);
        double step = 1.0e-3;
        
        if (effective_porosity >= porosity_th)
            step = 1.0;
        
        prop_eff = prop * ( std::pow(aux, porosity_mu) * step );
    }
    
    else
        AssertThrow( false, ExcNotImplemented() );
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::effective_gas_diffusivity(std::vector< Tensor<2,dim> >& prop_eff_vec) const
{
    Assert ( this->D_bulk.size() != 0, ExcMessage("compute_gas_diffusion not called before ConventionalCL::effective_gas_diffusivity.") );
    
    prop_eff_vec.resize( this->D_bulk.size() );
    

    
    for (unsigned int i=0; i < this->D_bulk.size(); ++i)
    {
        if ( this->s_vector.is_initialized() )
        {
            if ( this->s_vector[i] >= 0.0 )
            {
                for (unsigned int j=0; j < dim; ++j)
                    this->effective_gas_diffusivity( this->D_bulk[i], this->s_vector[i], prop_eff_vec[i][j][j] );
            }
            else
            {
                for (unsigned int j=0; j < dim; ++j)
                    this->effective_gas_diffusivity( this->D_bulk[i], 0.0, prop_eff_vec[i][j][j] );
            }
        }
        
        else if (this->capillary_pressure_vector.is_initialized())
        {
            std::vector<double> saturation;
            this->psd_pointer->set_capillary_pressure( this->capillary_pressure_vector );
            this->psd_pointer->set_saturation();
            this->psd_pointer->get_saturation(saturation);
            
            for (unsigned int j=0; j < dim; ++j)
                this->effective_gas_diffusivity( this->D_bulk[i], saturation[i], prop_eff_vec[i][j][j] );
        }
        
        else
        {
            for (unsigned int j=0; j < dim; ++j)
                this->effective_gas_diffusivity( this->D_bulk[i], 0.0, prop_eff_vec[i][j][j] );
        }
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::derivative_effective_gas_diffusivity(std::map< VariableNames, std::vector< Tensor<2,dim> > >& dprop_eff) const
{
    Assert( this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called in ConventionalCL::derivative_effective_gas_diffusivity."));
    Assert(((this->dD_bulk_dT.size()!=0) && (this->D_bulk.size()!=0)), ExcMessage("compute_gas_diffusion not called before ConventionalCL::derivative_effective_gas_diffusivity."));
    
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
                    {
                        for (unsigned int k=0; k<dim; ++k)
                            this->effective_gas_diffusivity( this->dD_bulk_dT[j], this->s_vector[j], Dprop[j][k][k] );
                    }
                    
                    else if (this->capillary_pressure_vector.is_initialized())
                    {
                        std::vector<double> saturation;
                        this->psd_pointer->set_capillary_pressure( this->capillary_pressure_vector );
                        this->psd_pointer->set_saturation();
                        this->psd_pointer->get_saturation(saturation);
                        
                        for (unsigned int k=0; k < dim; ++k)
                            this->effective_gas_diffusivity( this->D_bulk[i], saturation[i], Dprop[j][k][k] );
                    }
                    
                    else
                    {
                        for (unsigned int k=0; k<dim; ++k)
                            this->effective_gas_diffusivity( this->dD_bulk_dT[j], 0.0, Dprop[j][k][k] );
                    }
                }
                
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
                            Dprop[j][k][k] = this->D_bulk[j] * (-1.5) * (std::pow(epsilon_V.at(this->local_material_id()),1.5)) * (std::pow((1.0 - this->s_vector[j]),0.5));
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
                            if ( (epsilon_V.at(this->local_material_id())*(1.0 - this->s_vector[j])) >= porosity_th )
                            {
                                Dprop[j][k][k] = this->D_bulk[j] * (-1.0) * epsilon_V.at(this->local_material_id()) * porosity_mu * (std::pow( (1.0-porosity_th), ((-1.0)*porosity_mu) )) *
                                                 (std::pow( ((epsilon_V.at(this->local_material_id())*(1.0-this->s_vector[j])) - porosity_th), (porosity_mu-1.0) ));
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
                            Dprop[j][k][k] = this->D_bulk[j] * ds_dp[j] *(-1.5) * (std::pow(epsilon_V.at(this->local_material_id()),1.5)) * (std::pow((1.0 - saturation[j]),0.5));
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
                            if ( (epsilon_V.at(this->local_material_id())*(1.0 - saturation[j])) >= porosity_th )
                            {
                                Dprop[j][k][k] = this->D_bulk[j] * ds_dp[j] * (-1.0) * epsilon_V.at(this->local_material_id()) * porosity_mu * (std::pow( (1.0-porosity_th), ((-1.0)*porosity_mu) )) *
                                (std::pow( ((epsilon_V.at(this->local_material_id())*(1.0-saturation[j])) - porosity_th), (porosity_mu-1.0) ));
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
template<int dim>
void
NAME::ConventionalCL<dim>::effective_gas_diffusivity(Table< 2, Tensor< 2, dim > >& prop_eff) const
{
    prop_eff.reinit(this->gases.size(),this->gases.size());
    
    for (unsigned int i = 0; i<this->gases.size(); ++i)
    {
        for (unsigned int j = 0; j<this->gases.size(); ++j)
        {
            if ( i != j)        // Not computing self diffusion coefficients
            {
                for (unsigned int d = 0; d<dim; ++d)
                {
                    this->effective_gas_diffusivity( this->D_ECtheory(i,j), 0.0, prop_eff(i,j)[d][d] );
                }
            }
        }
    }
}
//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::effective_electron_conductivity(double& prop_eff) const
{
    // Initialize variable
    prop_eff = this->catalyst_support->get_electrical_conductivity();
    
    if (method_eff_property_solid == "Bruggemann")
        // Most common method
    {
        prop_eff = prop_eff*pow(epsilon_S.at(this->local_material_id()), 1.5);
    }
    else if (method_eff_property_solid == "Percolation")
        // Method used by M. Eikerling
    {
        double aux = fabs(epsilon_S.at(this->local_material_id()) - solid_th)/(1.0 - solid_th);
        double step = 0.0;
        if (epsilon_S.at(this->local_material_id()) >= solid_th)
            step = 1.0;
        prop_eff = prop_eff*pow(aux,solid_mu)*step;
    }
    else
    {
        FcstUtilities::log << "Unknown method to compute effective transport in the solid phase in "<<__FILE__ <<" line "<<__LINE__<<std::endl;
        abort();
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::effective_electron_conductivity(Tensor<2,dim>& prop_eff) const
{
    double iso_prop_eff;
    this->effective_electron_conductivity(iso_prop_eff);
    
    for (unsigned int i=0; i<dim; ++i)
        prop_eff[i][i] = iso_prop_eff;
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::derivative_effective_electron_conductivity(std::vector<double>& ) const
{
    const std::type_info& info = typeid(*this);
    FcstUtilities::log << "Pure function " << __FUNCTION__
    << " called in Class "
    << info.name() << std::endl;
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::effective_proton_conductivity(double& prop_eff) const
{
    // Initialize variable
    prop_eff = 0.0;
    this->electrolyte->proton_conductivity(prop_eff);
    
    if (method_eff_property_electrolyte == "Bruggemann")
        // Most common method
    {
        prop_eff = prop_eff*pow(epsilon_N.at(this->local_material_id()), 1.5);
    }    
    else if (method_eff_property_electrolyte == "Percolation")
        // Method used by M. Eikerling
    {
        double aux = fabs(epsilon_N.at(this->local_material_id()) - electrolyte_th)/(1.0 - electrolyte_th);
        double step = 0.0;
        if (epsilon_N.at(this->local_material_id()) >= electrolyte_th)
            step = 1.0;
        prop_eff = prop_eff*pow(aux,electrolyte_mu)*step;
    }
    
    else if (method_eff_property_electrolyte == "Iden11")
    {
        // Mesh structure
        prop_eff = prop_eff*pow(epsilon_N.at(this->local_material_id()), 1.6);
    }
    else
    {
        FcstUtilities::log << "Unknown method to compute effective transport in the electrolyte in "<<__FILE__ <<" line "<<__LINE__<<std::endl;
        abort();
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::effective_proton_conductivity(std::vector<double>& prop_eff) const
{	
    this->electrolyte->proton_conductivity(prop_eff);
    
    if (method_eff_property_electrolyte == "Bruggemann")
        // Most common method
    {
        for (unsigned int i=0; i<prop_eff.size(); ++i)
        {
            prop_eff[i] = prop_eff[i]*pow(epsilon_N.at(this->local_material_id()), 1.5);
        }
    }    
    else if (method_eff_property_electrolyte == "Percolation")
        // Method used by M. Eikerling
    {
        for (unsigned int i=0; i<prop_eff.size(); ++i)
        {
            double aux = fabs(epsilon_N.at(this->local_material_id()) - electrolyte_th)/(1.0 - electrolyte_th);
            double step = 0.0;
            if (epsilon_N.at(this->local_material_id()) >= electrolyte_th)
                step = 1.0;
            prop_eff[i] = prop_eff[i]*pow(aux,electrolyte_mu)*step;
        }
    }
    else if (method_eff_property_electrolyte == "Iden11")
    {
        for (unsigned int i=0; i<prop_eff.size(); ++i)
        {	prop_eff[i] = prop_eff[i] * (std::pow(epsilon_N.at(this->local_material_id()), 1.6));}
    }
    else
    {
        FcstUtilities::log << "Unknown method to compute effective transport in the electrolyte in "<<__FILE__ <<" line "<<__LINE__<<std::endl;
        abort();
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::derivative_effective_proton_conductivity(std::map< VariableNames, std::vector<double> >& Dsigma_eff) const
{
    Assert (this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called in ConventionalCL::derivative_effective_proton_conductivity."));
    
    this->electrolyte->proton_conductivity_derivative(Dsigma_eff);
    
    if (method_eff_property_electrolyte == "Bruggemann")
        for ( dmap_iter i=Dsigma_eff.begin(); i!=Dsigma_eff.end(); ++i)
            for (unsigned int j=0; j < i->second.size(); ++j)
                i->second[j] = i->second[j]*pow(epsilon_N.at(this->local_material_id()), 1.5);

    else if (method_eff_property_electrolyte == "Percolation")
        for ( dmap_iter i=Dsigma_eff.begin(); i!=Dsigma_eff.end(); ++i)
            for (unsigned int j=0; j < i->second.size(); ++j)
            {
                double aux = fabs(epsilon_N.at(this->local_material_id()) - electrolyte_th)/(1. - electrolyte_th);
                double step = 0.0;
                if (epsilon_N.at(this->local_material_id()) >= electrolyte_th)
                    step = 1.0;
                i->second[j] = i->second[j]*pow(aux,electrolyte_mu)*step;
            }

    else if (method_eff_property_electrolyte == "Iden11")
        for ( dmap_iter i=Dsigma_eff.begin(); i!=Dsigma_eff.end(); ++i)
            for (unsigned int j=0; j < i->second.size(); ++j)
                i->second[j] = i->second[j]*(std::pow(epsilon_N.at(this->local_material_id()), 1.6));
   
    else
    {
        FcstUtilities::log << "Unknown method to compute effective transport in the electrolyte in "<<__FILE__ <<" line "<<__LINE__<<std::endl;
        abort();
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::effective_water_diffusivity(double& prop_eff) const
{
    // Initialize variable
    prop_eff = 0.0;
    this->electrolyte->water_diffusivity(prop_eff);
    
    if (method_eff_property_electrolyte == "Bruggemann")
        // Most common method
    {
        prop_eff = prop_eff*pow(epsilon_N.at(this->local_material_id()), 1.5);
    }    
    else if (method_eff_property_electrolyte == "Percolation")
        // Method used by M. Eikerling
    {
        double aux = fabs(epsilon_N.at(this->local_material_id()) - electrolyte_th)/(1.0 - electrolyte_th);
        double step = 0.0;
        if (epsilon_N.at(this->local_material_id()) >= electrolyte_th)
            step = 1.0;
        prop_eff = prop_eff*pow(aux,electrolyte_mu)*step;
    }
    
    else if (method_eff_property_electrolyte == "Iden11")
    {
        // Mesh structure
        prop_eff = prop_eff*pow(epsilon_N.at(this->local_material_id()), 1.6);
    }
    else
    {
        FcstUtilities::log << "Unknown method to compute effective transport in the electrolyte in "<<__FILE__ <<" line "<<__LINE__<<std::endl;
        abort();
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::effective_water_diffusivity(std::vector<double>& prop_eff) const
{
    this->electrolyte->water_diffusivity(prop_eff);
    
    if (method_eff_property_electrolyte == "Bruggemann")
        // Most common method
    {
        for (unsigned int i=0; i<prop_eff.size(); ++i)
        {	prop_eff[i] = prop_eff[i]*pow(epsilon_N.at(this->local_material_id()), 1.5);}
    }    
    else if (method_eff_property_electrolyte == "Percolation")
        // Method used by M. Eikerling
    {
        for (unsigned int i=0; i<prop_eff.size(); ++i)
        {
            double aux = fabs(epsilon_N.at(this->local_material_id()) - electrolyte_th)/(1.0 - electrolyte_th);
            double step = 0.0;
            if (epsilon_N.at(this->local_material_id()) >= electrolyte_th)
                step = 1.0;
            prop_eff[i] = prop_eff[i]*pow(aux,electrolyte_mu)*step;
        }
    }
    else if (method_eff_property_electrolyte == "Iden11")
    {
        for (unsigned int i=0; i<prop_eff.size(); ++i)
        {	prop_eff[i] = prop_eff[i] * (std::pow(epsilon_N.at(this->local_material_id()), 1.6));}
    }
    else
    {
        FcstUtilities::log << "Unknown method to compute effective transport in the electrolyte in "<<__FILE__ <<" line "<<__LINE__<<std::endl;
        abort();
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::derivative_effective_water_diffusivity(std::map< VariableNames, std::vector<double> >& Dprop_eff) const
{
    Assert (this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called in ConventionalCL::derivative_effective_water_diffusivity."));
    
    this->electrolyte->water_diffusivity_derivative(Dprop_eff);
    
    if (method_eff_property_electrolyte == "Bruggemann")
        for ( dmap_iter i=Dprop_eff.begin(); i!=Dprop_eff.end(); ++i)
            for (unsigned int j=0; j < i->second.size(); ++j)
                i->second[j] = i->second[j]*pow(epsilon_N.at(this->local_material_id()), 1.5);
   
    else if (method_eff_property_electrolyte == "Percolation")
        for ( dmap_iter i=Dprop_eff.begin(); i!=Dprop_eff.end(); ++i)
            for (unsigned int j=0; j < i->second.size(); ++j)
            {
                double aux = fabs(epsilon_N.at(this->local_material_id()) - electrolyte_th)/(1.0 - electrolyte_th);
                double step = 0.0;
                if (epsilon_N.at(this->local_material_id()) >= electrolyte_th)
                    step = 1.0;
                i->second[j] = i->second[j]*pow(aux,electrolyte_mu)*step;
            }

    else if (method_eff_property_electrolyte == "Iden11")
        for ( dmap_iter i=Dprop_eff.begin(); i!=Dprop_eff.end(); ++i)
            for (unsigned int j=0; j < i->second.size(); ++j)
                i->second[j] = i->second[j]*(std::pow(epsilon_N.at(this->local_material_id()), 1.6));

    else
    {
        FcstUtilities::log << "Unknown method to compute effective transport in the electrolyte in "<<__FILE__ <<" line "<<__LINE__<<std::endl;
        abort();
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::effective_thermoosmotic_diffusivity(std::vector<double>& prop_eff) const
{
    this->electrolyte->thermoosmotic_coeff(prop_eff);
    
    if (method_eff_property_electrolyte == "Bruggemann")
        // Most common method
    {
        for (unsigned int i=0; i<prop_eff.size(); ++i)
        {       prop_eff[i] = prop_eff[i]*pow(epsilon_N.at(this->local_material_id()), 1.5);}
    }    
    else if (method_eff_property_electrolyte == "Percolation")
        // Method used by M. Eikerling
    {
        for (unsigned int i=0; i<prop_eff.size(); ++i)
        {
            double aux = fabs(epsilon_N.at(this->local_material_id()) - electrolyte_th)/(1.0 - electrolyte_th);
            double step = 0.0;
            if (epsilon_N.at(this->local_material_id()) >= electrolyte_th)
                step = 1.0;
            prop_eff[i] = prop_eff[i]*pow(aux,electrolyte_mu)*step;
        }
    }
    else if (method_eff_property_electrolyte == "Iden11")
    {
        for (unsigned int i=0; i<prop_eff.size(); ++i)
        {       prop_eff[i] = prop_eff[i] * (std::pow(epsilon_N.at(this->local_material_id()), 1.6));}
    }
    else
    {
        FcstUtilities::log << "Unknown method to compute effective transport in the electrolyte in "<<__FILE__ <<" line "<<__LINE__<<std::endl;
        abort();
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::derivative_effective_thermoosmotic_diffusivity(std::map< VariableNames,
                                                                          std::vector<double> >& Dprop_eff) const
{
    Assert (this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called in ConventionalCL::derivative_effective_thermoosmotic_diffusivity."));
    
    this->electrolyte->thermoosmotic_coeff_derivative(Dprop_eff);
    
    if (method_eff_property_electrolyte == "Bruggemann")
        for ( dmap_iter i=Dprop_eff.begin(); i!=Dprop_eff.end(); ++i)
            for (unsigned int j=0; j < i->second.size(); ++j)
                i->second[j] = i->second[j]*pow(epsilon_N.at(this->local_material_id()), 1.5);
  
    else if (method_eff_property_electrolyte == "Percolation")
        for ( dmap_iter i=Dprop_eff.begin(); i!=Dprop_eff.end(); ++i)
            for (unsigned int j=0; j < i->second.size(); ++j)
            {
                double aux = fabs(epsilon_N.at(this->local_material_id()) - electrolyte_th)/(1 - electrolyte_th);
                double step = 0.0;
                if (epsilon_N.at(this->local_material_id()) >= electrolyte_th)
                    step = 1.0;
                i->second[j] = i->second[j]*pow(aux,electrolyte_mu)*step;
            }

    else if (method_eff_property_electrolyte == "Iden11")
        for ( dmap_iter i=Dprop_eff.begin(); i!=Dprop_eff.end(); ++i)
            for (unsigned int j=0; j < i->second.size(); ++j)
                i->second[j] = i->second[j]*(std::pow(epsilon_N.at(this->local_material_id()), 1.6));
  
    else
    {
        FcstUtilities::log << "Unknown method to compute effective transport in the electrolyte in "<<__FILE__ <<" line "<<__LINE__<<std::endl;
        abort();
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::effective_thermal_conductivity(double& prop_eff) const
{
    if (method_eff_thermal == "Given")
        prop_eff = k_T.at(this->local_material_id());  //Isotropic case
        
    else
    {
        FcstUtilities::log << "Unknown method to compute effective transport in the electrolyte in "<<__FILE__ <<" line "<<__LINE__<<std::endl;
        abort();
    }  
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::effective_thermal_conductivity(std::vector< Tensor<2,dim> >& prop_eff) const
{
    if (method_eff_thermal == "Given")
    {
        Assert( prop_eff.size()!=0, ExcMessage("Vector should be sized before passing over to ConventionalCL::effective_thermal_conductivity as an argument.") );
        for (unsigned int i=0; i<prop_eff.size(); ++i)
            for (unsigned int j=0; j<dim; ++j)
                prop_eff[i][j][j] = k_T.at(this->local_material_id());
    }
    
    else
    {
        FcstUtilities::log << "Unknown method to compute effective transport in the electrolyte in "<<__FILE__ <<" line "<<__LINE__<<std::endl;
        abort();
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::derivative_effective_thermal_conductivity(std::vector< Tensor<2,dim> >& dK_dT) const
{
    if (method_eff_thermal == "Given")
    {
        Assert (this->n_quad != 0, ExcMessage("set_solution has not been probably called in ConventionalCL."));
        dK_dT.clear();
        dK_dT.resize(this->n_quad);
        for (unsigned int i=0; i<this->n_quad; ++i)
        {
            for (unsigned int j=0; j<dim; ++j)
            {
                dK_dT[i][j][j] = 0.;
            }
        }
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::liquid_permeablity(std::vector< Tensor<2,dim> >& k_l ) const
{
    Assert( this->s_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in ConventionalCL::liquid_permeablity.") );
    k_l.clear();
    k_l.resize(this->s_vector.size(), Tensor<2,dim>());
    
    const double s_irr_local(s_irr.at(this->local_material_id()));
            
    
    if (method_rel_liquid_permeability == "Kumbur07")
    // Ref: Kumbur E, Sharp K and Mench M. On the effectiveness of Leverett approach for describing the water transport in fuel cell diffusion media. Journal of Power Sources, 168(2):356-368, 2007.
    {
        for (unsigned int i=0; i<this->s_vector.size(); ++i)
        {
            for (unsigned int j=0; j<dim; ++j)
            {
                if (this->s_vector[i] < s_irr_local)
                    k_l[i][j][j] = 0.0;
                
                else
                {
                    double eff_saturation = (this->s_vector[i] - s_irr_local)/(1.0 - s_irr_local);
                    k_l[i][j][j] = abs_permeability.at(this->local_material_id()) * ( std::pow(eff_saturation, 2.16) );
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
                if (this->s_vector[i] < s_irr_local)
                    k_l[i][j][j] = 0.0;
                
                else
                {
                    double eff_saturation = (this->s_vector[i] - s_irr_local)/(1.0 - s_irr_local);
                    k_l[i][j][j] = abs_permeability.at(this->local_material_id()) * ( std::pow(eff_saturation, 3.0) );
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
NAME::ConventionalCL<dim>::derivative_liquid_permeablity(std::map< VariableNames, std::vector< Tensor<2,dim> > >& deriv_k_l) const
{
    Assert( this->s_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in ConventionalCL::derivative_liquid_permeablity.") );
    Assert (this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called in ConventionalCL::derivative_liquid_permeablity."));
    
    const double s_irr_local(s_irr.at(this->local_material_id()));
                
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
                        if (this->s_vector[j] < s_irr_local)
                            dk_l[j][k][k] = 0.0;
                        
                        else
                            dk_l[j][k][k] = abs_permeability.at(this->local_material_id()) * 2.16 * ( std::pow((this->s_vector[j] - s_irr_local), 1.16) ) * ( std::pow((1.0 - s_irr_local), -2.16) );
                    }
                }
            }
            
            else if (method_rel_liquid_permeability == "Wyllie")
            {
                for (unsigned int j=0; j<this->s_vector.size(); ++j)
                {
                    for (unsigned int k=0; k<dim; ++k)
                    {
                        if (this->s_vector[j] < s_irr_local)
                            dk_l[j][k][k] = 0.0;
                        
                        else
                            dk_l[j][k][k] = abs_permeability.at(this->local_material_id()) * 3.0 * ( std::pow((this->s_vector[j] - s_irr_local), 2.0) ) * ( std::pow((1.0 - s_irr_local), -3.0) );
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
NAME::ConventionalCL<dim>::saturated_liquid_permeablity_PSD(double & r_lp ) const
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
NAME::ConventionalCL<dim>::relative_liquid_permeability_PSD(std::vector< Tensor<2,dim> >& r_lp ) const
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
NAME::ConventionalCL<dim>::derivative_relative_liquid_permeablity_PSD(std::vector<double> & derivate_r_lp ) const
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
NAME::ConventionalCL<dim>::derivative_relative_liquid_permeablity_PSD(std::map< VariableNames, std::vector< Tensor<2,dim> > >& derivate_r_lp ) const
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
NAME::ConventionalCL<dim>::pcapillary(std::vector<double>& pc) const
{
    Assert( this->s_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in ConventionalCL::pcapillary.") );
    pc.clear();
    pc.resize(this->s_vector.size(), 0.0);
    
    if (method_capillary_function == "Ye07")
    // Ref: Ye Q and Van Nguyen T. Three-dimensional simulation of liquid water distribution in a PEMFC with experimentally measured capillary functions. JECS, 154(12):B1242-B1251, 2007.
    {
        for (unsigned int i=0; i<this->s_vector.size(); ++i)
        {
            if (this->s_vector[i] < 0.0)
                pc[i] = ( 2395.0 + 2431.0*( std::exp(92.36*(0.0-0.567)) - std::exp((-0.0088)*(0.0-0.567)) ) ) * 10.0; // dyne/cm^2
            else
                pc[i] = ( 2395.0 + 2431.0*( std::exp(92.36*(this->s_vector[i]-0.567)) - std::exp((-0.0088)*(this->s_vector[i]-0.567)) ) ) * 10.0; // dyne/cm^2
        }
    }
    
    else
        AssertThrow( false, ExcNotImplemented() );
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::saturation_from_capillary_equation(std::vector<double>& saturation) const
{
    Assert( this->capillary_pressure_vector.is_initialized(), ExcMessage("Liquid water capillary values have not been set using set_capillary_pressure method in ConventionalCL::liquid_permeablity.") );
    
    this->psd_pointer->set_capillary_pressure( this->capillary_pressure_vector );
    this->psd_pointer->set_critical_radius();
    this->psd_pointer->set_saturation();
    this->psd_pointer->get_saturation(saturation);
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::derivative_saturation_from_capillary_equation_PSD(std::vector<double>& derivative_saturation) const
{
    Assert( this->capillary_pressure_vector.is_initialized(), ExcMessage("Liquid water capillary values have not been set using set_capillary_pressure method in ConventionalCL::liquid_permeablity.") );
    
    this->psd_pointer->set_capillary_pressure( this->capillary_pressure_vector );
    this->psd_pointer->set_critical_radius();
    this->psd_pointer->set_saturation();
    this->psd_pointer->get_derivative_saturation(derivative_saturation);
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::dpcapillary_dsat(std::vector<double> & dpc_ds) const
{
    Assert( this->s_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in ConventionalCL::dpcapillary_dsat.") );
    dpc_ds.clear();
    dpc_ds.resize(this->s_vector.size(), 0.0);
    
    if (method_capillary_function == "Ye07")
    // Ref: Ye Q and Van Nguyen T. Three-dimensional simulation of liquid water distribution in a PEMFC with experimentally measured capillary functions. JECS, 154(12):B1242-B1251, 2007.
    {
        for (unsigned int i=0; i<this->s_vector.size(); ++i)
        {
            if (this->s_vector[i] < 0.0)
                dpc_ds[i] = 0.0;
            else
                dpc_ds[i] =  2431.0 * 10.0 * ( (std::exp(92.36*(this->s_vector[i]-0.567)))*92.36 + (std::exp((-0.0088)*(this->s_vector[i]-0.567)))*0.0088 ); // [dyne/cm^2]
        }
    }
    
    else
        AssertThrow( false, ExcNotImplemented() );
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::derivative_dpcapillary_dsat(std::map< VariableNames, std::vector<double> > & d2pc_dsdu) const
{
    Assert( this->s_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in ConventionalCL::derivative_dpcapillary_dsat.") );
    Assert( this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called in ConventionalCL::derivative_dpcapillary_dsat.") );
    
    for (unsigned int i=0; i<this->derivative_flags.size(); ++i)
    {
        std::vector<double> d2pc(this->s_vector.size(), 0.0);
        
        if (this->derivative_flags[i] == liquid_water_saturation)
        {
            if (method_capillary_function == "Ye07")
            {
                for (unsigned int j=0; j<this->s_vector.size(); ++j)
                {
                    if (this->s_vector[j] < 0.0)
                        d2pc[j] = 0.0;
                    else
                        d2pc[j] =  2431.0*10.0*( (std::exp(92.36*(this->s_vector[j]-0.567)))*92.36*92.36 - (std::exp((-0.0088)*(this->s_vector[j]-0.567)))*0.0088*0.0088 );//[dyne/cm^2]
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
NAME::ConventionalCL<dim>::interfacial_surface_area(std::vector<double>& a_lv) const
{
    Assert( this->s_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in ConventionalCL::interfacial_surface_area.") );
    
    a_lv.clear();
    a_lv.resize( this->s_vector.size(), 0.0 );
    
    for (unsigned int i=0; i < this->s_vector.size(); ++i)
    {
        if ( (this->s_vector[i] > 0.0) && (this->s_vector[i] < 1.0) )
            a_lv[i] = 44000.0 * (std::pow(this->s_vector[i], 1.25)) * (std::pow((1.0-this->s_vector[i]),2.5));  // cm^2/cm^3
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::derivative_interfacial_surface_area(std::map< VariableNames, std::vector<double> >& deriv_a_lv) const
{
    Assert( this->s_vector.is_initialized(), ExcMessage("Liquid water saturation values have not been set using set_saturation method in ConventionalCL::derivative_interfacial_surface_area.") );
    Assert( this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called in ConventionalCL::derivative_interfacial_surface_area.") );
    
    for (unsigned int i=0; i < this->derivative_flags.size(); ++i)
    {
        std::vector<double> da_lv(this->s_vector.size(), 0.0);
        
        if (this->derivative_flags[i] == liquid_water_saturation)
        {
            for (unsigned int j=0; j<this->s_vector.size(); ++j)
            {
                if ( (this->s_vector[j] > 0.0) && (this->s_vector[j] < 1.0) )
                {
                    da_lv[j] = 44000.0 * ( (1.25*(std::pow(this->s_vector[j],0.25))*(std::pow((1.0-this->s_vector[j]),2.5))) - 
                                           ((std::pow(this->s_vector[j],1.25))*2.5*(std::pow((1.0-this->s_vector[j]),1.5))) );
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
NAME::ConventionalCL<dim>::interfacial_surface_area_PSD(std::vector<double>& a_lv) const
{
    Assert( this->capillary_pressure_vector.is_initialized(), ExcMessage("Liquid water capillary pressure values have not been set using set_saturation method in ConventionalCL::interfacial_surface_area.") );
    a_lv.clear();
    a_lv.resize( this->capillary_pressure_vector.size(), 0.0 );
    
    this->psd_pointer->set_capillary_pressure( this->capillary_pressure_vector );
    this->psd_pointer->set_critical_radius();
    this->psd_pointer->set_saturation();
    this->psd_pointer->get_liquid_gas_interfacial_surface(a_lv);

    // cm^2/cm^3
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::derivative_interfacial_surface_area_PSD(std::vector<double>& deriv_a_lv) const
{
    Assert( this->capillary_pressure_vector.is_initialized(), ExcMessage("Liquid water capillary pressure values have not been set using set_saturation method in ConventionalCL::interfacial_surface_area.") );
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
NAME::ConventionalCL<dim>::derivative_interfacial_surface_area_PSD(std::map< VariableNames, std::vector<double> >& deriv_a_lv) const
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
//------------ PROTECTED MEMBER FUNCTIONS:
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::compute_volume_fraction()
{
    epsilon_S.clear();
    epsilon_N.clear();
    epsilon_V.clear();

    if (method_porosity == "marc" || method_porosity == "NafionLoading" || method_porosity == "ICRatio")
    {
        for(unsigned int i=0; i<this->material_ids.size(); ++i)
        {
            unsigned int current_id = this->material_ids.at(i);

            // Solid phase volume fraction:
            epsilon_S[current_id] = (1.0/rho_Pt + (1.0-prc_Pt.at(current_id))/(prc_Pt.at(current_id)*rho_c))*(V_Pt.at(current_id)*1.e-3);
        
            // Electrolyte volume fraction:
            epsilon_N[current_id] = loading_N.at(current_id)/(1.0-loading_N.at(current_id))*(1.0/prc_Pt.at(current_id))*(1.0/rho_N)*(V_Pt.at(current_id)*1.e-3);
        
            // Porosity
            epsilon_V[current_id] = 1.0 - epsilon_S.at(current_id) - epsilon_N.at(current_id);
        }
    }
}


//---------------------------------------------------------------------------
template <int dim>
void 
NAME::ConventionalCL<dim>::compute_Av()
{
    Av.clear();

    for(unsigned int i=0; i<this->material_ids.size(); ++i)
    {
        unsigned int current_id = this->material_ids.at(i);

        // Compute the active area of catalyst. Depends on V_Pt and %Pt
        double A_0 = 0.0;
        if (method_Av == "Marr")
        {
            A_0 = 2.2779e6*pow(prc_Pt.at(current_id),3.0) - 1.5857e6*pow(prc_Pt.at(current_id),2.0) - 2.0153e6*prc_Pt.at(current_id) + 1.5950e6;
            Av[current_id] = A_0*(V_Pt.at(current_id)*1e-3);
        }
        else if (method_Av == "ETEK06")
        {
            A_0 = -4.5646e5*pow(prc_Pt.at(current_id),3.0) + 1.0618e6*pow(prc_Pt.at(current_id),2.0) - 1.8564e6*prc_Pt.at(current_id) + 1.5955e6;
            Av[current_id] = A_0*(V_Pt.at(current_id)*1e-3);
        }
        else if (method_Av == "ETEK07")
        {
            A_0 = 0.7401e7*pow(prc_Pt.at(current_id),4.0) - 1.8105e7*pow(prc_Pt.at(current_id),3.0) + 1.5449e7*pow(prc_Pt.at(current_id),2.0) - 0.6453e7*prc_Pt.at(current_id) + 0.2054e7;
            Av[current_id] = A_0*(V_Pt.at(current_id)*1e-3);
        }
    }
}


//---------------------------------------------------------------------------
template <int dim>
void
NAME::ConventionalCL<dim>::derivative_effective_proton_conductivity_wrt_electrolyte_loading(double& Dsigma_eff) const
{
//    if (method_eff_property_electrolyte == "Bruggemann")
//    {
//        Dsigma_eff = electrolyte_proton_conductivity*1.5*pow(epsilon_N, 0.5);
//    }
//    else if (method_eff_property_electrolyte == "Percolation")
//    {
//        double aux = fabs(epsilon_N - electrolyte_th)/(1.0 - electrolyte_th);
//        double Daux_Depsilon = 1/(1 - electrolyte_th);
//        double step = 0.0;
//        if (epsilon_N >= electrolyte_th)
//            step = 1.0;
//        Dsigma_eff = electrolyte_proton_conductivity*electrolyte_mu*pow(aux,electrolyte_mu-1)*Daux_Depsilon*step;
//    }
//    else
//    {
//        FcstUtilities::log << "Unknown method to compute effective conduction in "<<__FILE__ <<" line "<<__LINE__<<std::endl;
//        abort();
//    }
}


//---------------------------------------------------------------------------
template <int dim>
inline void
NAME::ConventionalCL<dim>::derivative_volume_fractions(double &Depsilon_S,
                                                        double &Depsilon_V,
                                                        double &Depsilon_N) const
{
// for (unsigned int i=0; i<this->derivative_flags.size(); ++i)
// {
// 	if (this->name == "Anode catalyst layer")
// 	{
// 		if (this->derivative_flags[i] == "epsilon_N_cat_a")
// 		{
// 			Depsilon_V = depsilon_V_cat_depsilon_N_cat();
// 			Depsilon_S = 0.0;
// 			Depsilon_N = 1.0;
// 		}
// 		else if (this->derivative_flags[i] == "V_Pt_a")
// 		{
// 			Depsilon_V = depsilon_V_cat_depsilon_S_cat()*depsilon_S_cat_dVPt(prc_Pt);
// 			Depsilon_S = depsilon_S_cat_dVPt(prc_Pt);
// 			Depsilon_N = 0.0;
// 		}
// 		else if (this->derivative_flags[i] == "prc_Pt_a")
// 		{
// 			Depsilon_V = depsilon_V_cat_depsilon_S_cat()*depsilon_S_cat_dprc_Pt(V_Pt, prc_Pt);
// 			Depsilon_S = depsilon_S_cat_dprc_Pt(V_Pt, prc_Pt);
// 			Depsilon_N = 0.0;
// 		}
// 		else if (this->derivative_flags[i] == "epsilon_V_gdl_a")
// 		{
// 			Depsilon_V = 0.0;
// 			Depsilon_S = 0.0;
// 			Depsilon_N = 0.0;
// 		}
// 		else
// 		{
// 			Depsilon_V = 0.0;
// 			Depsilon_S = 0.0;
// 			Depsilon_N = 0.0;
// 		}
// 	}
// 	if (this->name == "Cathode catalyst layer")
// 	{
// 		if (this->derivative_flags[i] == "epsilon_N_cat_c")
// 		{
// 			Depsilon_V = depsilon_V_cat_depsilon_N_cat();
// 			Depsilon_S = 0.0;
// 			Depsilon_N = 1.0;
// 		}
// 		else if (this->derivative_flags[i] == "V_Pt_c")
// 		{
// 			Depsilon_V = depsilon_V_cat_depsilon_S_cat()*depsilon_S_cat_dVPt(prc_Pt);
// 			Depsilon_S = depsilon_S_cat_dVPt(prc_Pt);
// 			Depsilon_N = 0.0;
// 		}
// 		else if (this->derivative_flags[i] == "prc_Pt_c")
// 		{
// 			Depsilon_V = depsilon_V_cat_depsilon_S_cat()*depsilon_S_cat_dprc_Pt(V_Pt, prc_Pt);
// 			Depsilon_S = depsilon_S_cat_dprc_Pt(V_Pt, prc_Pt);
// 			Depsilon_N = 0.0;
// 		}
// 		else if (this->derivative_flags[i] == "epsilon_V_gdl_c")
// 		{
// 			Depsilon_V = 0.0;
// 			Depsilon_S = 0.0;
// 			Depsilon_N = 0.0;
// 		}
// 		else
// 		{
// 			Depsilon_V = 0.0;
// 			Depsilon_S = 0.0;
// 			Depsilon_N = 0.0;
// 		}
// 	}
// }
}

//---------------------------------------------------------------------------
template <int dim>
inline void
NAME::ConventionalCL<dim>::print_layer_properties() const
{

	for (unsigned int i=0; i<this->material_ids.size(); ++i)
	{
		unsigned int current_id = this->material_ids.at(i);

		FcstUtilities::log<<"======= CATALYST LAYER STRUCTURE ======"<<std::endl;
		FcstUtilities::log<<"=== "<<this->name<<" === #:"<< (i + 1) <<" ====="<<std::endl;
		FcstUtilities::log<<"Solid phase volume fraction:"<<epsilon_S.at(current_id)<<std::endl;
		FcstUtilities::log<<"Ionomer in the CL:"<<epsilon_N.at(current_id)<<std::endl;
		FcstUtilities::log<<"Porosity in the CL:"<<epsilon_V.at(current_id)<<std::endl;
		FcstUtilities::log<<"Active area in the CL:"<<Av.at(current_id)<<std::endl;
		FcstUtilities::log<<"---------------------------------------"<<std::endl;
		FcstUtilities::log<<"Method to compute eff_property_electrolyte: "<<method_eff_property_electrolyte<<std::endl;
		FcstUtilities::log<<"======================================="<<std::endl;
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations. 
template class NAME::ConventionalCL<deal_II_dimension>;
