//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: gas_diffusion_layer.cc
//    - Description: Base Gas Diffusion Layer Class. It implements the interface for other gas diffusion layer classes
//        and some common methods.
//    - Developers: M. Secanell
//    - $Id: gas_diffusion_layer.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include <layers/gas_diffusion_layer.h>


namespace NAME = FuelCellShop::Layer;

template <int dim>
const std::string NAME::GasDiffusionLayer<dim>::concrete_name ("GasDiffusionLayer");

//---------------------------------------------------------------------------
template <int dim>
NAME::GasDiffusionLayer<dim>::GasDiffusionLayer()
  : NAME::PorousLayer<dim> ()
{}

//---------------------------------------------------------------------------
template <int dim>
NAME::GasDiffusionLayer<dim>::GasDiffusionLayer(const std::string& name)
  : NAME::PorousLayer<dim> (name)
{}

//---------------------------------------------------------------------------
template <int dim>
NAME::GasDiffusionLayer<dim>::~GasDiffusionLayer()
{}
  
//---------------------------------------------------------------------------
template <int dim>
void
NAME::GasDiffusionLayer<dim>::declare_parameters (const std::string& name, 
                                                  ParameterHandler &param) const
{
    FuelCellShop::Layer::PorousLayer<dim>::declare_parameters(name,param);

    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection(name);
        {   
            param.declare_entry("Gas diffusion layer type",
                                "DummyGDL",
                                Patterns::Selection("DummyGDL|SGL24BA|DesignFibrousGDL"),
                                "Select a type of gas diffusion layer. Once this value is selected, fill in "
                                "the appropriate subsection. Options: DummyGDL|SGL24BA|DesignFibrousGDL");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
}
          
//---------------------------------------------------------------------------
template <int dim>
void
NAME::GasDiffusionLayer<dim>::initialize (ParameterHandler &param)
{

  NAME::PorousLayer<dim>::initialize(param);

  param.enter_subsection("Fuel cell data");
  {
    param.enter_subsection(this->name);
    {
    }
    param.leave_subsection();
  }
  param.leave_subsection();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::GasDiffusionLayer<dim>::test_class() const
{
  double pressure(1.0); //atm 
  double temperature(273); //K
  FuelCellShop::Layer::GasDiffusionLayer<dim> test("test");
  FuelCellShop::Material::Oxygen oxygen;
  FuelCellShop::Material::Nitrogen nitrogen;
  FcstUtilities::log<<"Computing diffusion coefficient of: "<<oxygen.name_material()<<" in "<<nitrogen.name_material()<<" at "<<pressure<<" Pa and "<<temperature<<" Kelvin"<<std::endl;
  std::vector<FuelCellShop::Material::PureGas*> gases;
  gases.push_back(&oxygen);
  gases.push_back(&nitrogen);
  test.set_gases_and_compute(gases, pressure, temperature);
  Table<2, double> D_eff;
  test.effective_gas_diffusivity(D_eff);
  
  FcstUtilities::log<<"The diffusion coefficients are: D_11="<<D_eff[0][0]<<" D_12 = "<<D_eff[0][1]<<"m2/s"<<std::endl;  
  FcstUtilities::log<<"The diffusion coefficients are: D_21="<<D_eff[1][0]<<" D_22 = "<<D_eff[1][1]<<"m2/s"<<std::endl;  
  FcstUtilities::log<<"The experimental value at 273K and 1atm is: 0.181 cm2/s"<<std::endl;
  
  int index_g(0);
  int index_sol(0);
  test.get_gas_index(&oxygen, index_g);
  test.get_gas_index(&nitrogen,index_sol);
  int marc(0); //K
  FuelCellShop::Material::WaterVapor p;
  test.get_gas_index(&p,marc);
  
  FcstUtilities::log<<"The gas "<<oxygen.name_material()<<" is "<<test.get_gas_pointer(index_g)->name_material()<<std::endl;  
  FcstUtilities::log<<"The gas "<<nitrogen.name_material()<<" is "<<test.get_gas_pointer(index_sol)->name_material()<<std::endl;  
  FcstUtilities::log<<"The unknown gas returns index "<<marc<<std::endl; 
  
  
  Table<2, Tensor<2,dim> > D_eff_tensor;
  test.effective_gas_diffusivity(D_eff_tensor);
  //FcstUtilities::log<<"The diffusion coefficients are: D_11xx="<<D_eff_tensor(0,0)[0][0]<<" D_11xy= "<<D_eff_tensor[0][0][0][1]<<"D_11yx="<<D_eff_tensor[0][0][1][0]<<" D_11yy= "<<D_eff_tensor[0][0][1][1]<<std::endl;//<<" D_11xy= "<<D_eff_tensor[0][0][0][2]<<std::endl;
  //FcstUtilities::log<<"The diffusion coefficients are: D_12xx="<<D_eff_tensor[0][1][0][0]<<" D_12xy= "<<D_eff_tensor[0][1][0][1]<<"D_12yx="<<D_eff_tensor[0][1][1][0]<<" D_12yy= "<<D_eff_tensor[0][1][1][1]<<std::endl;//<<" D_12xy= "<<D_eff_tensor[0][1][0][2]<<std::endl;
  //FcstUtilities::log<<"The diffusion coefficients are: D_21xx="<<D_eff_tensor[1][0][0][0]<<" D_21xy= "<<D_eff_tensor[1][0][0][1]<<"D_21yx="<<D_eff_tensor[1][0][1][0]<<" D_21yy= "<<D_eff_tensor[1][0][1][1]<<std::endl;//<<" D_21xy= "<<D_eff_tensor[1][0][0][2]<<std::endl;
  //FcstUtilities::log<<"The diffusion coefficients are: D_22xx="<<D_eff_tensor[1][1][0][0]<<" D_22xy= "<<D_eff_tensor[1][1][0][1]<<"D_22yx="<<D_eff_tensor[1][1][1][0]<<" D_22yy= "<<D_eff_tensor[1][1][1][1]<<std::endl;//<<" D_22xy= "<<D_eff_tensor[1][1][0][2]<<std::endl;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations.
template class NAME::GasDiffusionLayer<deal_II_dimension>;
