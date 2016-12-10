//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: PSD_base.h
//    - Description: Base class for pore size distribution model.
//    - Developers: 2009-13 by Marc Secanell, University of Alberta
//                  2013-14 by Jie Zhou, University of Alberta
//
//---------------------------------------------------------------------------

#include <microscale/PSD_base.h>

namespace NAME = FuelCellShop::MicroScale;

//---------------------------------------------------------------------------
template <int dim>
NAME::BasePSD<dim>::BasePSD(const std::string& name)
: 
name(name)
{}

//---------------------------------------------------------------------------
 
template <int dim>
void
NAME::BasePSD<dim>::declare_parameters (ParameterHandler &param) const
{
    
    std::string list_of_PSDs;
    /**
    for (typename FuelCellShop::MicroScale::BasePSD<dim>::_mapFactory::iterator iterator = FuelCellShop::MicroScale::BasePSD<dim>::get_mapFactory()->begin(); 
         iterator != FuelCellShop::MicroScale::BasePSD<dim>::get_mapFactory()->end(); iterator++)         
    {
        list_of_PSDs += iterator->first;
        if (iterator != FuelCellShop::MicroScale::BasePSD<dim>::get_mapFactory()->end())
            list_of_PSDs += "|";
        
    }
    
    list_of_PSDs.erase(list_of_PSDs.end()-1);*/
    list_of_PSDs+="DualPSD|HIPSD|HOPSD|NonePSD";
    
    param.enter_subsection("PSD parameters");
    {
        param.enter_subsection("BasePSD");
        {     
            param.declare_entry("psd type",
                                "NonePSD",
                                Patterns::Selection(list_of_PSDs),
                                "Pore size distribution type");

            param.declare_entry("Gamma",
                                "0.0",
                                Patterns::Double(0.0),
                                "Water-air interface surface tension the unit is J/m^2");
            
            param.declare_entry("lambda",
                                "0.0",
                                Patterns::Double(0.0),
                                "interconnection of the pores");            
            
            param.declare_entry("Volume fraction Hydrophilic",
                                "0.0",
                                Patterns::Double(0.0),
                                "Volume fraction of hydrophilic pores in the sample currently set to be 1.0 at any circumstance");
            
            param.declare_entry("Volume fraction Hydrophobic",
                                "0.0",
                                Patterns::Double(0.0),
                                "Volume fraction of hydrophilic pores in the sample currently set to be 1.0 at any circumstance");
            
            param.declare_entry("Mode probability global",
                                "0.0",
                                Patterns::List(Patterns::Double(0.0)),
                                "Contribution of the distribution mode into the PSD ");
            
            param.declare_entry("Mode characteristic radius global",
                                "0.0",
                                Patterns::List(Patterns::Double(0.0)),
                                "Characteristic pore size of the distribution mode into the PSD ");
            
            param.declare_entry("Mode width global",
                                "0.0",
                                Patterns::List(Patterns::Double(0.0)),
                                "Characteristic pore size of the distribution mode into the PSD ");    
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}


//---------------------------------------------------------------------------
template <int dim>
void
NAME::BasePSD<dim>::_initialize (ParameterHandler &param) 
{
    param.enter_subsection("PSD parameters");
    {
        param.enter_subsection("BasePSD");
        {
            gamma =  param.get_double("Gamma");            
            lamda = param.get_double("lambda");            
            F_HI =  param.get_double("Volume fraction Hydrophilic");            
            F_HO =  param.get_double("Volume fraction Hydrophobic");            
            f_k = FcstUtilities::string_to_number<double>( Utilities::split_string_list( param.get("Mode probability global") ) );            
            r_k = FcstUtilities::string_to_number<double>( Utilities::split_string_list( param.get("Mode characteristic radius global") ) );            
            s_k = FcstUtilities::string_to_number<double>( Utilities::split_string_list( param.get("Mode width global") ) );
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations.
template class NAME::BasePSD<deal_II_dimension>;