//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2012, 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: dummy_GDL.cc
//    - Description: Implementation of a GDL class that setup us all properties from file
//    - Developers: M. Secanell
//
//---------------------------------------------------------------------------

#include <layers/dummy_GDL.h>

namespace NAME = FuelCellShop::Layer; 

template <int dim>
const std::string NAME::DummyGDL<dim>::concrete_name ("DummyGDL");

template <int dim>
NAME::DummyGDL<dim> const* NAME::DummyGDL<dim>::PROTOTYPE = new NAME::DummyGDL<dim>();


//---------------------------------------------------------------------------
template <int dim>
NAME::DummyGDL<dim>::DummyGDL()
: FuelCellShop::Layer::GasDiffusionLayer<dim> ()
{
    //FcstUtilities::log<<" Register DummyGDL GDL to FactoryMap"<<std::endl;
    this->get_mapFactory()->insert(std::pair<std::string, FuelCellShop::Layer::GasDiffusionLayer<dim>* >(concrete_name, this));
}
            
            
//---------------------------------------------------------------------------
template <int dim>
NAME::DummyGDL<dim>::DummyGDL(std::string name)
  : NAME::GasDiffusionLayer<dim> (name)
{
  FcstUtilities::log<<" of type DummyGDL"<<std::endl;
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyGDL<dim>::declare_parameters (const std::string& name, 
                                         ParameterHandler &param) const
{
    
    FuelCellShop::Layer::GasDiffusionLayer<dim>::declare_parameters(name, param);
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection(name);
        {
            param.enter_subsection(concrete_name); //-- Transport for the anisotropic case:
            {
                
                // ISOTROPIC PROPERTIES:
                param.declare_entry ("Method effective gas transport properties",
                                     "Diffusibility",
                                     Patterns::Selection("Given|Diffusibility"),
                                     "Method used to compute effective transport properties in the void phase.");
                param.declare_entry ("Diffusivity ratio, [-]",
                                     "0.3", //1atm, 353K
                                     Patterns::Double(),
                                     "Ratio of effective diffusivity to molecular diffusivity (also known as diffusibility)");
                param.declare_entry ("Oxygen diffusion coefficient, [cm^2/s]",
                                     "0.2741", //1atm, 353K
                                     Patterns::Double(),
                                     "Oxygen diffusion coefficient given by experiment");
                param.declare_entry ("Water vapour diffusion coefficient, [cm^2/s]",
                                     "0.29646", 
                                     Patterns::Double(),
                                     "Water vapour diffusion coefficient given by experiment");
                param.declare_entry ("Electrical conductivity, [S/cm]",
                                     "100", // [S/cm]
                                     Patterns::Double(),
                                     "Effective cond. if given is used, otherwise conducitivity of the raw material. Units [S/cm]");
                param.declare_entry ("Thermal conductivity, [W/(cm K)]",
                                     "16", // [S/cm]
                                     Patterns::Double(),
                                     "Effective thermal cond. if given is used, otherwise conducitivity of the raw material. Units [W/(cm K)]");
                // ANISOTROPIC PROPERTIES IF NEEDED:    
                param.declare_entry ("Anisotropic transport",
                                     "false", // [S/cm]
                                     Patterns::Bool(),
                                     "Boolean variable. Set to true if we want to account for anisotropy of the GDL");
                //--- XX
                param.declare_entry ("Diffusivity ratio X, [-]",
                                     "0.3", //1atm, 353K
                                     Patterns::Double(),
                                     "Ratio of effective diffusivity to molecular diffusivity (also known as diffusibility)");
                param.declare_entry ("Oxygen diffusion coefficient X, [cm^2/s]",
                                     "0.2741", //1atm, 353K
                                     Patterns::Double(),
                                     "Oxygen diffusion coefficient given by experiment");
                param.declare_entry ("Water vapour diffusion coefficient X, [cm^2/s]",
                                     "0.29646", 
                                     Patterns::Double(),
                                     "Water vapour diffusion coefficient given by experiment");
                param.declare_entry ("Electrical conductivity X, [S/cm]",
                                     "100", // [S/cm]
                                     Patterns::Double(),
                                     "Component X of the electrical conductivity tensor. Units [S/cm]");
                param.declare_entry ("Thermal conductivity X, [W/(cm K)]",
                                     "16", // [S/cm]
                                     Patterns::Double(),
                                     "Effective thermal cond. if given is used, otherwise conducitivity of the raw material. Units [W/(cm K)]");
                // YY
                param.declare_entry ("Diffusivity ratio Y, [-]",
                                     "0.3", //1atm, 353K
                                     Patterns::Double(),
                                     "Ratio of effective diffusivity to molecular diffusivity (also known as diffusibility)");                
                param.declare_entry ("Oxygen diffusion coefficient Y, [cm^2/s]",
                                     "0.2741", //1atm, 353K
                                     Patterns::Double(),
                                     "Oxygen diffusion coefficient given by experiment");
                param.declare_entry ("Water vapour diffusion coefficient Y, [cm^2/s]",
                                     "0.29646", 
                                     Patterns::Double(),
                                     "Water vapour diffusion coefficient given by experiment");
                param.declare_entry ("Electrical conductivity Y, [S/cm]",
                                     "100", // [S/cm]
                                     Patterns::Double(),
                                     "Component Y of the electrical conductivity tensor. Units [S/cm]");
                param.declare_entry ("Thermal conductivity Y, [W/(cm K)]",
                                     "16", // [S/cm]
                                     Patterns::Double(),
                                     "Effective thermal cond. if given is used, otherwise conducitivity of the raw material. Units [W/(cm K)]");
                // ZZ
                param.declare_entry ("Diffusivity ratio Z, [-]",
                                     "0.3", //1atm, 353K
                                     Patterns::Double(),
                                     "Ratio of effective diffusivity to molecular diffusivity (also known as diffusibility)");                
                param.declare_entry ("Oxygen diffusion coefficient Z, [cm^2/s]",
                                     "0.2741", //1atm, 353K
                                     Patterns::Double(),
                                     "Oxygen diffusion coefficient given by experiment");
                param.declare_entry ("Water vapour diffusion coefficient Z, [cm^2/s]",
                                     "0.29646", 
                                     Patterns::Double(),
                                     "Water vapour diffusion coefficient given by experiment");
                param.declare_entry ("Electrical conductivity Z, [S/cm]",
                                     "100", 
                                     Patterns::Double(),
                                     "Component Z of the electrical conductivity tensor. Units [S/cm]");
                param.declare_entry ("Thermal conductivity Z, [W/(cm K)]",
                                     "16", // [S/cm]
                                     Patterns::Double(),
                                     "Effective thermal cond. if given is used, otherwise conducitivity of the raw material. Units [W/(cm K)]");
                
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
NAME::DummyGDL<dim>::initialize (ParameterHandler &param)
{
    
    NAME::GasDiffusionLayer<dim>::initialize(param);
    
    param.enter_subsection("Fuel cell data"); 
    {
        param.enter_subsection(this->name); 
        { 
            param.enter_subsection(concrete_name); 
            { 
                // Anisotropy
                anisotropy = param.get_bool ("Anisotropic transport");
                
                method_eff_property_pores = param.get("Method effective gas transport properties");
                
                D_D0.resize(dim);
                D_O2.resize(dim);
                D_wv.resize(dim);
                sigma_e.resize(dim);
                k_T.resize(dim);
                
                if (anisotropy == true)
                {
                    switch (dim)
                    {
                    // X
                    case 1:
                        D_D0[0] = param.get_double("Diffusivity ratio X, [-]");
                        D_O2[0] = param.get_double("Oxygen diffusion coefficient X, [cm^2/s]");
                        D_wv[0] = param.get_double("Water vapour diffusion coefficient X, [cm^2/s]");
                        sigma_e[0] = param.get_double("Electrical conductivity X, [S/cm]");	  
                        k_T[0] = param.get_double("Thermal conductivity X, [W/(cm K)]");                        
                        break;

                    case 2:
                        D_D0[0] = param.get_double("Diffusivity ratio X, [-]");                      
                        D_O2[0] = param.get_double("Oxygen diffusion coefficient X, [cm^2/s]");
                        D_wv[0] = param.get_double("Water vapour diffusion coefficient X, [cm^2/s]");
                        sigma_e[0] = param.get_double("Electrical conductivity X, [S/cm]");   
                        k_T[0] = param.get_double("Thermal conductivity X, [W/(cm K)]");
                        // Y
                        D_D0[1] = param.get_double("Diffusivity ratio Y, [-]");                       
                        D_O2[1] = param.get_double("Oxygen diffusion coefficient Y, [cm^2/s]");
                        D_wv[1] = param.get_double("Water vapour diffusion coefficient Y, [cm^2/s]");
                        sigma_e[1] = param.get_double("Electrical conductivity Y, [S/cm]");
                        k_T[1] = param.get_double("Thermal conductivity Y, [W/(cm K)]");
                        break;
                        
                    case 3:
                        D_D0[0] = param.get_double("Diffusivity ratio X, [-]");                        
                        D_O2[0] = param.get_double("Oxygen diffusion coefficient X, [cm^2/s]");
                        D_wv[0] = param.get_double("Water vapour diffusion coefficient X, [cm^2/s]");
                        sigma_e[0] = param.get_double("Electrical conductivity X, [S/cm]");   
                        k_T[0] = param.get_double("Thermal conductivity X, [W/(cm K)]");
                        // Y
                        D_D0[1] = param.get_double("Diffusivity ratio Y, [-]");                       
                        D_O2[1] = param.get_double("Oxygen diffusion coefficient Y, [cm^2/s]");
                        D_wv[1] = param.get_double("Water vapour diffusion coefficient Y, [cm^2/s]");
                        sigma_e[1] = param.get_double("Electrical conductivity Y, [S/cm]");
                        k_T[1] = param.get_double("Thermal conductivity Y, [W/(cm K)]");
                        // Z
                        D_D0[2] = param.get_double("Diffusivity ratio Z, [-]");                       
                        D_O2[2] = param.get_double("Oxygen diffusion coefficient Z, [cm^2/s]");	  
                        D_wv[2] = param.get_double("Water vapour diffusion coefficient Z, [cm^2/s]");
                        sigma_e[2] = param.get_double("Electrical conductivity Z, [S/cm]");
                        k_T[2] = param.get_double("Thermal conductivity Z, [W/(cm K)]");
                        break;
                    }
                }
                else
                {
                    for (unsigned int i =0; i < D_O2.size(); i++)
                    {
                        D_D0[0] = param.get_double("Diffusivity ratio, [-]");
                        D_O2[i] = param.get_double("Oxygen diffusion coefficient, [cm^2/s]");
                        D_wv[i] = param.get_double("Water vapour diffusion coefficient, [cm^2/s]");
                        sigma_e[i] = param.get_double("Electrical conductivity, [S/cm]");
                        k_T[i] = param.get_double("Thermal conductivity, [W/(cm K)]");
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
NAME::DummyGDL<dim>::effective_gas_diffusivity(std::vector< Tensor<2,dim> >& prop_eff_vec) const
{
    Assert( this->T_vector.is_initialized(), ExcMessage("Temperature is not set inside the porous layer. Use set_temperature method before calling this function.") );

    prop_eff_vec.resize(this->T_vector.size());
    
    std::vector<double> saturation(this->T_vector.size(), 0.0);    
    if ( this->s_vector.is_initialized() )
    {
        for (unsigned int i; i<this->s_vector.size(); i++)
            saturation[i] = this->s_vector[i];
    }
    else if ( this->capillary_pressure_vector.is_initialized() )
    {
        this->psd_pointer->set_capillary_pressure( this->capillary_pressure_vector );
        this->psd_pointer->set_saturation();
        this->psd_pointer->get_saturation(saturation);
    }
    
    if (method_eff_property_pores.compare("Given") == 0)
    {
        
        int solute_index, solvent_index;
        
        this->get_gas_index(this->solute_gas, solute_index);
        this->get_gas_index(this->solvent_gas, solvent_index);
        
        std::string solute_name(this->gases[solute_index]->name_material());
        std::string solvent_name(this->gases[solvent_index]->name_material());
        
        if (solute_name.compare("oxygen") == 0)   {
            if (solvent_name.compare("nitrogen") == 0)
                for (unsigned int q = 0; q < this->T_vector.size(); q++)
                    for (unsigned int d = 0; d<dim; d++)
                        prop_eff_vec[q][d][d] = Units::convert(D_O2[d], Units::UNIT2, Units::C_UNIT2);
        }
        else if (solute_name.compare("water") == 0) {
            if (solvent_name.compare("nitrogen") == 0)
                for (unsigned int q = 0; q < this->T_vector.size(); q++)
                    for (unsigned int d = 0; d<dim; d++)
                        prop_eff_vec[q][d][d] = Units::convert(D_wv[d], Units::UNIT2, Units::C_UNIT2);
        }
        else {
            FcstUtilities::log << "DummyGDL<dim>::effective_gas_diffusivity not implemented"<< std::endl;
            exit(1);
        }
    }
    else if (method_eff_property_pores.compare("Diffusibility") == 0) 
    {
        // -- Loop over quad points
        for (unsigned int q=0; q < this->D_bulk.size(); ++q)
            // -- Loop over dimension
            for (unsigned int k=0; k<dim; ++k) 
            {
                prop_eff_vec[q][k][k] = (1-saturation[q])*this->D_bulk[q]*D_D0[k];
            }
    }
    else
    {
        AssertThrow(false, ExcMessage ("method_eff_property_pores parameter not implemented"));
    }
}
//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyGDL<dim>::derivative_effective_gas_diffusivity(std::map< VariableNames, std::vector< Tensor<2,dim> > >&dprop_eff)  const
{   
    Assert( this->T_vector.is_initialized(), ExcMessage("Temperature is not set inside the porous layer. Use set_temperature method before calling this function.") );
    Assert( this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called in DummyGDL::derivative_effective_gas_diffusivity."));
    Assert(((this->dD_bulk_dT.size()!=0) && (this->D_bulk.size()!=0)), ExcMessage("compute_gas_diffusion not called before DummyGDL::derivative_effective_gas_diff   ")); 
   
    
    std::vector< Tensor<2,dim> > Dprop( this->dD_bulk_dT.size(), Tensor<2,dim>() );
    
    
    std::vector<double> saturation(this->T_vector.size(), 0.0);    
    if ( this->s_vector.is_initialized() )
    {
        for (unsigned int i; i<this->s_vector.size(); i++)
            saturation[i] = this->s_vector[i];
    }
    else if ( this->capillary_pressure_vector.is_initialized() )
    {
        this->psd_pointer->set_capillary_pressure( this->capillary_pressure_vector );
        this->psd_pointer->set_saturation();
        this->psd_pointer->get_saturation(saturation);
    }
    
    if (method_eff_property_pores.compare("Given") == 0)
    {
        for (unsigned int i=0; i < this->derivative_flags.size(); ++i)
        {
            if ( this->derivative_flags[i] == temperature_of_REV )
            {
                std::vector< Tensor<2,dim> > Dprop( this->T_vector.size(), Tensor<2,dim>() );            
                dprop_eff[this->derivative_flags[i]] = Dprop;
            }
        }
    }
    else if (method_eff_property_pores.compare("Diffusibility") == 0) 
    {
        for (unsigned int i=0; i < this->derivative_flags.size(); ++i)
        {
            if ( this->derivative_flags[i] == temperature_of_REV )
            {
                for (unsigned int q=0; q < this->D_bulk.size(); ++q)
                    for (unsigned int d=0; d<dim; ++d) 
                    {
                        Dprop[q][d][d] = (1-saturation[q])*this->dD_bulk_dT[q]*D_D0[d];
                    }
            }
            dprop_eff[ this->derivative_flags[i] ] = Dprop;            
        }
    }
    else
    {
        AssertThrow(false, ExcMessage ("method_eff_property_pores parameter not implemented"));
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyGDL<dim>::effective_gas_diffusivity(Table< 2, double> & prop_eff) const
{
    prop_eff.reinit(this->gases.size(),this->gases.size());
    
    unsigned int dimension = 0;

    for (unsigned int i = 0; i<this->gases.size(); i++)
    {
        std::string solute_name(this->gases[i]->name_material());

        for (unsigned int j = i+1; j<this->gases.size(); j++)
        {  
            std::string solvent_name(this->gases[j]->name_material());
            
            if (solute_name.compare("oxygen") == 0)
            {
                
                if (solvent_name.compare("nitrogen") == 0)
                {
                    prop_eff(i,j) = Units::convert(D_O2[dimension], Units::UNIT2, Units::C_UNIT2);
                    prop_eff(j,i) = prop_eff(i,j);
                    break;
                }
                else if (solvent_name.compare("water") == 0)
                {
                    prop_eff(i,j) = 1e200;
                    prop_eff(j,i) = prop_eff(i,j);
                    break;
                }
            }
            else if (solute_name.compare("water") == 0)
            {
                
                if (solvent_name.compare("nitrogen") == 0)
                {
                    prop_eff(i,j) = Units::convert(D_wv[dimension], Units::UNIT2, Units::C_UNIT2);
                    prop_eff(j,i) = prop_eff(i,j);
                    break;
                }
            }
        }
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyGDL<dim>::effective_gas_diffusivity(Table< 2, Tensor< 2, dim > > &prop_eff) const
{  
    
    prop_eff.reinit(this->gases.size(),this->gases.size());
    
    for (unsigned int i = 0; i<this->gases.size(); i++)
    {
        std::string name_i;
        name_i = this->gases[i]->name_material();
        
        for (unsigned int j = i+1; j<this->gases.size(); j++)
        {  
            std::string name_j;
            name_j = this->gases[j]->name_material();
            
            if (name_i.compare("oxygen") == 0)
            {
                
                if (name_j.compare("nitrogen") == 0)
                {
                    for (unsigned int d = 0; d<dim; d++)
                    {
                        prop_eff(i,j)[d][d] = Units::convert(D_O2[d], Units::UNIT2, Units::C_UNIT2);
                        prop_eff(j,i)[d][d] = prop_eff(i,j)[d][d];                       
                    }
                }
                else
                {
                    FcstUtilities::log << "Species "<<name_i.c_str() <<" and "<< name_j.c_str() << " diffusivity requested in "
                    << "DummyGDL<dim>::effective_gas_diffusivity not implemented"<< std::endl;
                    exit(1);
                } 
                
            }
            else if (name_i.compare("water") == 0)
            {
                
                if (name_j.compare("nitrogen") == 0)
                {
                    for (unsigned int d = 0; d<dim; d++)
                    {
                        prop_eff(i,j)[d][d] = Units::convert(D_wv[d], Units::UNIT2, Units::C_UNIT2);
                        prop_eff(j,i)[d][d] = prop_eff(i,j)[d][d];
                        break;
                    }
                }
                else
                {
                    FcstUtilities::log << "Species "<<name_i.c_str() <<" and "<< name_j.c_str() << " diffusivity requested in "
                    << "DummyGDL<dim>::effective_gas_diffusivity not implemented"<< std::endl;
                    exit(1);
                } 
            }
            else
            {
                FcstUtilities::log << "Species "<<name_i.c_str() <<" and "<< name_j.c_str() << " diffusivity requested in "
                << "DummyGDL<dim>::effective_gas_diffusivity not implemented"<< std::endl;
                exit(1);
            } 
        }
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyGDL<dim>::effective_electron_conductivity(double& prop_eff) const
{
    if (anisotropy == false)
    {
        prop_eff = sigma_e[0];
    }
    else
    {
        FcstUtilities::log <<"The member function " << __FUNCTION__
        <<" called in Class DummyGDL can only be used for isotropic materials. "
        <<" Set anisotropy to false"<<std::endl;
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyGDL<dim>::effective_electron_conductivity(Tensor<2,dim>& prop_eff) const
{
     for (unsigned int i=0; i<dim; i++)
       prop_eff[i][i] = this->sigma_e[i]; // Include in declare paramters and initialize.
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyGDL<dim>::effective_thermal_conductivity(std::vector< Tensor<2,dim> >& prop_eff) const
{
    Assert (this->T_vector.size() != 0, ExcMessage("Temperature values have not been set. set_temperature has not been probably called."));
    prop_eff.clear();
    prop_eff.resize(this->T_vector.size());
    
    for (unsigned int i = 0; i < prop_eff.size(); ++i)
    {
        for (unsigned int j = 0; j < dim; ++j)
            prop_eff[i][j][j] = k_T[j];
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyGDL<dim>::derivative_effective_thermal_conductivity(std::vector< Tensor<2,dim> >& dK_dT) const
{
    Assert (this->T_vector.size() != 0, ExcMessage("Temperature values have not been set. set_temperature has not been probably called."));
    dK_dT.clear();
    dK_dT.resize(this->T_vector.size());    
    
    for (unsigned int i=0; i<dK_dT.size(); ++i)
    {
        for (unsigned int j=0; j<dim; ++j)
            dK_dT[i][j][j] = 0.0;
    }
    
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations. 
template class NAME::DummyGDL<deal_II_dimension>;