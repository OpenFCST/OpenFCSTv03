// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: dummy_CL.cc
// - Description: This class characterizes a conventional catalyst layer
//                and defines constant effective properties
// - Developers: Marc Secanell and Madhur Bhaiya
// - Id: $Id$
//
// ----------------------------------------------------------------------------

#include <layers/dummy_CL.h>

namespace NAME = FuelCellShop::Layer;

template <int dim>
const std::string NAME::DummyCL<dim>::concrete_name ("DummyCL");

template <int dim>
NAME::DummyCL<dim> const* NAME::DummyCL<dim>::PROTOTYPE = new NAME::DummyCL<dim>();

//---------------------------------------------------------------------------
template <int dim>
NAME::DummyCL<dim>::DummyCL()
  : NAME::CatalystLayer<dim> ()
{
    //FcstUtilities::log<<" Register DummyCL CL to FactoryMap"<<std::endl;
    this->get_mapFactory()->insert(std::pair<std::string, FuelCellShop::Layer::CatalystLayer<dim>* >(concrete_name, this));
}

//---------------------------------------------------------------------------
template <int dim>
NAME::DummyCL<dim>::DummyCL(const std::string& name)
  : NAME::CatalystLayer<dim> (name)
{
  FcstUtilities::log<<" of type DummyCL"<<std::endl;
}

//---------------------------------------------------------------------------
template <int dim>
NAME::DummyCL<dim>::~DummyCL()
{}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyCL<dim>::initialize (ParameterHandler &param)
{
    NAME::CatalystLayer<dim>::initialize(param);

    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection(this->name);
        {
            param.enter_subsection(concrete_name); //-- Transport for the anisotropic case:
            {
                if( !param.get("Oxygen diffusion coefficient, [cm^2/s]").empty() )
                     D_O2 = FcstUtilities::string_to_map<unsigned int, double>( Utilities::split_string_list( param.get("Oxygen diffusion coefficient, [cm^2/s]") ) );
                if( !param.get("Water vapour diffusion coefficient, [cm^2/s]").empty() )
                     D_wv = FcstUtilities::string_to_map<unsigned int, double>( Utilities::split_string_list( param.get("Water vapour diffusion coefficient, [cm^2/s]") ) );
                if( !param.get("Electrical conductivity, [S/cm]").empty() )
                     sigma_e = FcstUtilities::string_to_map<unsigned int, double>( Utilities::split_string_list( param.get("Electrical conductivity, [S/cm]") ) );
                if( !param.get("Protonic conductivity, [S/cm]").empty() )
                     sigma_m = FcstUtilities::string_to_map<unsigned int, double>( Utilities::split_string_list( param.get("Protonic conductivity, [S/cm]") ) );
                if( !param.get("Active area [cm^2/cm^3]").empty() )
                     Av = FcstUtilities::string_to_map<unsigned int, double>( Utilities::split_string_list( param.get("Active area [cm^2/cm^3]") ) );
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
NAME::DummyCL<dim>::effective_gas_diffusivity(Table< 2, Tensor< 2, dim > >& prop_eff) const
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
                        prop_eff(i,j)[d][d] = Units::convert(D_O2.at(this->local_material_id()), Units::UNIT2, Units::C_UNIT2);
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
                        prop_eff(i,j)[d][d] = Units::convert(D_wv.at(this->local_material_id()), Units::UNIT2, Units::C_UNIT2);
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
NAME::DummyCL<dim>::effective_gas_diffusivity(std::vector< Tensor<2,dim> >& prop_eff_vec) const
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
                    prop_eff_vec[q][d][d] = Units::convert(D_O2.at(this->local_material_id()), Units::UNIT2, Units::C_UNIT2);
    }
    else if (solute_name.compare("water") == 0) {
        if (solvent_name.compare("nitrogen") == 0)
            for (unsigned int q = 0; q < this->T_vector.size(); q++)
                for (unsigned int d = 0; d<dim; d++)
                    prop_eff_vec[q][d][d] = Units::convert(D_wv.at(this->local_material_id()), Units::UNIT2, Units::C_UNIT2);
    }
    else {
        FcstUtilities::log << "DummyMPL<dim>::effective_gas_diffusivity not implemented"<< std::endl;
        exit(1);
    }
}
//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyCL<dim>::derivative_effective_gas_diffusivity(std::map< VariableNames, std::vector< Tensor<2,dim> > >&dprop_eff)  const
{   
   
    for (unsigned int i=0; i < this->derivative_flags.size(); ++i)
    {
        if ( this->derivative_flags[i] == temperature_of_REV )  {
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
NAME::DummyCL<dim>::effective_electron_conductivity(double& prop_eff) const
{
    prop_eff = sigma_e.at(this->local_material_id());
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyCL<dim>::effective_electron_conductivity(Tensor<2,dim>& prop_eff) const
{
    prop_eff.clear();
    for (unsigned int i=0; i<dim; i++)
        prop_eff[i][i] = sigma_e.at(this->local_material_id());
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyCL<dim>::effective_proton_conductivity(double& prop_eff) const
{
    prop_eff = sigma_m.at(this->local_material_id());
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyCL<dim>::effective_proton_conductivity(std::vector<double>& prop_eff) const
{
    Assert( prop_eff.size() != 0, ExcMessage("Input vector should be initialized before passing to this method.") );
    for (unsigned int q = 0; q < prop_eff.size(); ++q)
        prop_eff[q] = sigma_m.at(this->local_material_id());
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyCL<dim>::derivative_effective_proton_conductivity(std::map< VariableNames, std::vector<double> >& dprop_eff) const
{
    for (unsigned int i=0; i < this->derivative_flags.size(); ++i)
    {
        if ( this->derivative_flags[i] == temperature_of_REV ) {
            std::vector< double > Dprop( this->T_vector.size() , 0.0);            
            dprop_eff[this->derivative_flags[i]] = Dprop;
        }
        else if (this->derivative_flags[i] == membrane_water_content ) {
            std::vector< double > Dprop( this->T_vector.size() , 0.0 );            
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
NAME::DummyCL<dim>::current_density(std::vector<double>& coef)
{
    this->kinetics->current_density(coef);
    for (unsigned int i = 0; i<coef.size(); ++i)
        coef[i] *= Av.at(this->local_material_id());
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyCL<dim>::derivative_current_density(std::map< VariableNames, std::vector<double> >& dcoef_du)
{
    this->kinetics->derivative_current(dcoef_du);
    for (std::map< VariableNames, std::vector<double> >::iterator iter=dcoef_du.begin(); iter!=dcoef_du.end(); ++iter)
        for (unsigned int q=0; q<iter->second.size(); ++q)
            iter->second[q] *= Av.at(this->local_material_id());
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations.
template class NAME::DummyCL<deal_II_dimension>;