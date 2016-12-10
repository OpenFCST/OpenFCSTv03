//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2012, 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: SGL_24_BA.cc
//    - Description: Header file for a specific type of gas diffusion layer, i.e. SIGRACET 24 BA
//    - Developers: M. Secanell
//    - Id: $Id: SGL_24_BA.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include <layers/SGL_24_BA.h>

namespace NAME = FuelCellShop::Layer;

template <int dim>
const std::string NAME::SGL24BA<dim>::concrete_name ("SGL24BA");

template <int dim>
NAME::SGL24BA<dim> const* NAME::SGL24BA<dim>::PROTOTYPE = new NAME::SGL24BA<dim>();
//boost::shared_ptr< NAME::SGL24BA<dim> > NAME::SGL24BA<dim>::PROTOTYPE = new NAME::SGL24BA<dim>();


//---------------------------------------------------------------------------
template <int dim>
NAME::SGL24BA<dim>::SGL24BA()
  : NAME::GasDiffusionLayer<dim> ()
{
    //FcstUtilities::log<<" Register SGL24BA GDL to FactoryMap"<<std::endl;
    this->get_mapFactory()->insert(std::pair<std::string, GasDiffusionLayer<dim>* > (this->concrete_name, this) );           
}

//---------------------------------------------------------------------------
// Random values are given at the moment, once the experimental data
// is compiled..proper values can be put in constructor definition

template <int dim>
NAME::SGL24BA<dim>::SGL24BA(const std::string& name)
  : NAME::GasDiffusionLayer<dim> (name)
{
  FcstUtilities::log<<" Created a SGL_24_BA GDL"<<std::endl;
  
  // Initialize internal data
  porosity_over_tortuosity_X = 0.5;
  porosity_over_tortuosity_Y = 0.5;
  porosity_over_tortuosity_Z = 0.5;
  
  electron_conductivity_X = 100;
  electron_conductivity_Y = 100;
  electron_conductivity_Z = 100;
  
  thermal_conductivity_X = 100;
  thermal_conductivity_Y = 100;
  thermal_conductivity_Z = 100;
  
  // Setup electron_conductivity tensor:
   switch (dim)
    {
        case 1:
            porosity_over_tortuosity[0][0] = porosity_over_tortuosity_X;
            electron_conductivity[0][0] = electron_conductivity_X;
            thermal_conductivity[0][0] = thermal_conductivity_X;
            break;
        case 2:
            porosity_over_tortuosity[0][0] = porosity_over_tortuosity_X;
            porosity_over_tortuosity[1][1] = porosity_over_tortuosity_Y;
            
            electron_conductivity[0][0] = electron_conductivity_X;
            electron_conductivity[1][1] = electron_conductivity_Y;
            
            thermal_conductivity[0][0] = thermal_conductivity_X;
            thermal_conductivity[1][1] = thermal_conductivity_Y;
            break;
        case 3:
            porosity_over_tortuosity[0][0] = porosity_over_tortuosity_X;
            porosity_over_tortuosity[1][1] = porosity_over_tortuosity_Y;
            porosity_over_tortuosity[2][2] = porosity_over_tortuosity_Z;

            electron_conductivity[0][0] = electron_conductivity_X;
            electron_conductivity[1][1] = electron_conductivity_Y;
            electron_conductivity[2][2] = electron_conductivity_Z;
            
            thermal_conductivity[0][0] = thermal_conductivity_X;
            thermal_conductivity[1][1] = thermal_conductivity_Y;
            thermal_conductivity[2][2] = thermal_conductivity_Z;
            break;
        default:
            AssertThrow(false, ExcNotImplemented());
    }

}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::SGL24BA<dim>::declare_parameters(const std::string& name, 
                                       ParameterHandler &param) const
{
    FuelCellShop::Layer::GasDiffusionLayer<dim>::declare_parameters(name, param);
    
    /// Nothing extra to declare since all properties are hard-coded.
    /*
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection(name);
        {
            param.enter_subsection(concrete_name); //-- Transport for the anisotropic case:
            {
                
            }
        }
        param.leave_subsection();        
    }
    param.leave_subsection();
    */
    
}





template <int dim>
void
NAME::SGL24BA<dim>::initialize (ParameterHandler &param)
{
    
    NAME::GasDiffusionLayer<dim>::initialize(param);
    
    /// Nothing extra to declare since all properties are hard-coded.
    /*
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection(name);
        {
            param.enter_subsection(concrete_name); //-- Transport for the anisotropic case:
            {
                
            }
        }
        param.leave_subsection();        
    }
    param.leave_subsection();
    */
}















            
//---------------------------------------------------------------------------
template <int dim>
void
NAME::SGL24BA<dim>::effective_gas_diffusivity(Table<2, Tensor<2,dim> >&prop_eff) const
{
    prop_eff.reinit(this->gases.size(),this->gases.size());
    
    for (unsigned int i=0; i<this->gases.size(); i++)
        for (unsigned int j=0; j<this->gases.size(); j++)
        {
            if (i!=j)
            {
                prop_eff(i,j) = porosity_over_tortuosity*this->D_ECtheory(i,j);            
            }
        }
    
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::SGL24BA<dim>::effective_electron_conductivity(Tensor<2,dim>& prop_eff) const
{
    //Initialize variable
    prop_eff.clear();
    prop_eff = electron_conductivity;
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::SGL24BA<dim>::effective_thermal_conductivity(Tensor<2,dim>& prop_eff) const
{
    //Initialize variable
    prop_eff.clear();
    prop_eff = thermal_conductivity;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations. 
template class NAME::SGL24BA<deal_II_dimension>;
