//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2012, 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: SGL_24_BC.cc
//    - Description: Header file for a specific type of microporous layer, i.e. SIGRACET 24 BC
//    - Developers: M. Secanell and Madhur Bhaiya
//
//---------------------------------------------------------------------------

#include <layers/dummy_MPL.h>

namespace NAME = FuelCellShop::Layer;

template <int dim>
const std::string NAME::DummyMPL<dim>::concrete_name ("DummyMPL");

template <int dim>
NAME::DummyMPL<dim> const* NAME::DummyMPL<dim>::PROTOTYPE = new NAME::DummyMPL<dim>();


//---------------------------------------------------------------------------
template <int dim>
NAME::DummyMPL<dim>::DummyMPL()
  : NAME::MicroPorousLayer<dim> ()
{
    //FcstUtilities::log<<" Register DummyMPL MPL to FactoryMap"<<std::endl;
    this->get_mapFactory()->insert(std::pair<std::string, MicroPorousLayer<dim>* > (this->concrete_name, this) );           
}

//---------------------------------------------------------------------------
template <int dim>
NAME::DummyMPL<dim>::DummyMPL(std::string name)
: NAME::MicroPorousLayer<dim> (name)
{
    FcstUtilities::log<<" Created a DummyMPL MPL"<<std::endl;
    
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::DummyMPL<dim>::declare_parameters (std::string name, ParameterHandler &param) const
{   
    FuelCellShop::Layer::MicroPorousLayer<dim>::declare_parameters(name,param);
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection(name);
        {
            param.enter_subsection(concrete_name); //-- Transport for the anisotropic case:
            {
                
                // ISOTROPIC PROPERTIES:
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
NAME::DummyMPL<dim>::initialize (ParameterHandler &param)
{
    FuelCellShop::Layer::MicroPorousLayer<dim>::initialize(param);
    
    param.enter_subsection("Fuel cell data"); 
    {
        param.enter_subsection(this->name); 
        { 
            param.enter_subsection(concrete_name); 
            { 
                // Anisotropy
                this->anisotropy = param.get_bool ("Anisotropic transport");
                
                this->oxygen_diffusivity.clear();
                this->water_diffusivity.clear();
                this->electrical_conductivity.clear();
                this->thermal_conductivity.clear();
                
                if (this->anisotropy == true)
                {
                    
                    this->oxygen_diffusivity[0][0] = param.get_double("Oxygen diffusion coefficient X, [cm^2/s]");
                    this->water_diffusivity[0][0] = param.get_double("Water vapour diffusion coefficient X, [cm^2/s]");
                    this->electrical_conductivity[0][0] = param.get_double("Electrical conductivity X, [S/cm]");       
                    this->thermal_conductivity[0][0] = param.get_double("Thermal conductivity X, [W/(cm K)]");
                    // Y
                    this->oxygen_diffusivity[1][1] = param.get_double("Oxygen diffusion coefficient Y, [cm^2/s]");
                    this->water_diffusivity[1][1] = param.get_double("Water vapour diffusion coefficient Y, [cm^2/s]");
                    this->electrical_conductivity[1][1] = param.get_double("Electrical conductivity Y, [S/cm]");       
                    this->thermal_conductivity[1][1] = param.get_double("Thermal conductivity Y, [W/(cm K)]");
                    // z
                    this->oxygen_diffusivity[2][2] = param.get_double("Oxygen diffusion coefficient Z, [cm^2/s]");
                    this->water_diffusivity[2][2] = param.get_double("Water vapour diffusion coefficient Z, [cm^2/s]");
                    this->electrical_conductivity[2][2] = param.get_double("Electrical conductivity Z, [S/cm]");       
                    this->thermal_conductivity[2][2] = param.get_double("Thermal conductivity Z, [W/(cm K)]");
                    
                }
                else
                {
                    for (unsigned int i =0; i < dim; i++)
                    {
                        this->oxygen_diffusivity[i][i] = param.get_double("Oxygen diffusion coefficient, [cm^2/s]");
                        this->water_diffusivity[i][i] = param.get_double("Water vapour diffusion coefficient, [cm^2/s]");
                        this->electrical_conductivity[i][i] = param.get_double("Electrical conductivity, [S/cm]");
                        this->thermal_conductivity[i][i] = param.get_double("Thermal conductivity, [W/(cm K)]");
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
NAME::DummyMPL<dim>::effective_gas_diffusivity(Table<2, Tensor<2,dim> >&prop_eff) const
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
                        prop_eff(i,j)[d][d] = Units::convert(this->oxygen_diffusivity[d][d], Units::UNIT2, Units::C_UNIT2);
                        prop_eff(j,i)[d][d] = prop_eff(i,j)[d][d];
                    }
                }
                else
                {
                    FcstUtilities::log << "Species "<< name_j << " diffusivity requested in"
                    << "DummyCL<dim>::effective_gas_diffusivity"
                    <<"not implemented"<< std::endl;
                    exit(1);
                }

            }
            else if (name_i.compare("water") == 0)
            {

                if (name_j.compare("nitrogen") == 0)
                {
                    for (unsigned int d = 0; d<dim; d++)
                    {
                        prop_eff(i,j)[d][d] = Units::convert(this->water_diffusivity[d][d], Units::UNIT2, Units::C_UNIT2);
                        prop_eff(j,i)[d][d] = prop_eff(i,j)[d][d];
                        break;
                    }
                }
                else
                {
                    FcstUtilities::log << "Species "<< name_j << " diffusivity requested in"
                    << "DummyCL<dim>::effective_gas_diffusivity"
                    <<"not implemented"<< std::endl;
                    exit(1);
                }
            }
            else
            {
                FcstUtilities::log << "Species "<< name_i << " diffusivity requested in"
                << "DummyCL<dim>::effective_gas_diffusivity"
                <<"not implemented"<< std::endl;
                exit(1);
            }
        }
    }
}


//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyMPL<dim>::effective_gas_diffusivity(std::vector< Tensor<2,dim> >& prop_eff_vec) const
{
    Assert(this->gases.size() != 2, ExcMessage("Number of gases should be two in PorousLayer::set_gases_and_compute method in order to use effective_gas_diffusivity(std::vector< Tensor<2,dim> >& prop_eff_vec)."));
   
    prop_eff_vec.resize(this->T_vector.size());
    
    int solute_index, solvent_index;
    
    this->get_gas_index(this->solute_gas, solute_index);
    this->get_gas_index(this->solvent_gas, solvent_index);
    
    std::string solute_name(this->gases[solute_index]->name_material());
    std::string solvent_name(this->gases[solvent_index]->name_material());
    
    if (solute_name.compare("oxygen") == 0)   {
        if (solvent_name.compare("nitrogen") == 0)
            for (unsigned int q = 0; q < this->T_vector.size(); q++)
                for (unsigned int d = 0; d<dim; d++)
                    prop_eff_vec[q][d][d] = this->oxygen_diffusivity[d][d]*Units::convert(1.0, Units::UNIT2, Units::C_UNIT2);
    }
    else if (solute_name.compare("water") == 0) {
        if (solvent_name.compare("nitrogen") == 0)
            for (unsigned int q = 0; q < this->T_vector.size(); q++)
                for (unsigned int d = 0; d<dim; d++)
                    prop_eff_vec[q][d][d] = this->water_diffusivity[d][d]*Units::convert(1.0, Units::UNIT2, Units::C_UNIT2);
    }
    else {
        FcstUtilities::log << "DummyMPL<dim>::effective_gas_diffusivity not implemented"<< std::endl;
        exit(1);
    }
}
//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyMPL<dim>::derivative_effective_gas_diffusivity(std::map< VariableNames, std::vector< Tensor<2,dim> > >&dprop_eff)  const
{   
   
    for (unsigned int i=0; i < this->derivative_flags.size(); ++i)
    {
        if ( this->derivative_flags[i] == temperature_of_REV ) {
            std::vector< Tensor<2,dim> > Dprop( this->T_vector.size(), Tensor<2,dim>() );            
            dprop_eff[this->derivative_flags[i]] = Dprop;
        }
        else {
            FcstUtilities::log << "DummyMPL<dim>::effective_gas_diffusivity not implemented"<< std::endl;
            exit(1);
        }
    }
}
                                             
//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyMPL<dim>::effective_electron_conductivity(Tensor<2,dim>& prop_eff) const
{
    //Initialize variable
    prop_eff.clear();
    prop_eff = this->electrical_conductivity;
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyMPL<dim>::effective_thermal_conductivity(Tensor<2,dim>& prop_eff) const
{
    //Initialize variable
    prop_eff.clear();
    prop_eff = this->thermal_conductivity;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations. 
template class NAME::DummyMPL<deal_II_dimension>;
