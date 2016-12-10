//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: 
//    - Description: 
//    - Developers: 
//    - $Id: porous_layer_test.cc 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#include <porous_layer_test.h>

//---------------------------------------------------------------------------
void PorousLayerTest::setup(){
    
    ParameterHandler param;
    // Note since I cannot create a parent, I will test the functions of the parent in one of the children:
    FuelCellShop::Layer::GasDiffusionLayer<dim>::declare_GasDiffusionLayer_parameters("Sample Gas Diffusion Layer", param);
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Sample Gas Diffusion Layer");
        {
            param.set("Gas diffusion layer type", "DesignFibrousGDL");
            {
                param.enter_subsection("Generic data");
                { 
                    param.set("Use Bosanquet approx.", true);
                    param.set("Knudsen pore radius [um]",5e-2);
                }
                param.leave_subsection();
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    
    layer = FuelCellShop::Layer::GasDiffusionLayer<dim>::create_GasDiffusionLayer("Sample Gas Diffusion Layer", param);

}

//---------------------------------------------------------------------------
void PorousLayerTest::molecular_diffusion_test()
{
    
    std::vector< FuelCellShop::Material::PureGas* > cathode_gases;
    cathode_gases.push_back(&oxygen);
    cathode_gases.push_back(&nitrogen);
    
    double pressure (1);  
    unsigned int n_quad = 3;
    FuelCellShop::SolutionVariable temperature = FuelCellShop::SolutionVariable(298.15, n_quad, temperature_of_REV);
    
    layer->set_gases (cathode_gases, pressure);

    layer->set_temperature( temperature );
    
    layer->compute_gas_diffusion(&oxygen,&nitrogen);
    
    std::vector<double> D_m;
    double D_m_check(2.0304e-5); //m^2/s
    layer->molecular_gas_diffusion_coefficient(D_m);
    
    std::ostringstream streamOut;
    for (unsigned int i=0; i<n_quad; ++i){
        streamOut <<"The value of the D_m is: "<<D_m[i]<<". The expected value is: "<<D_m_check<<std::endl;
        TEST_ASSERT_MSG(std::fabs((D_m_check - D_m[i])/D_m_check)<1e-3, streamOut.str().c_str());
        streamOut.flush();
    }
}

//---------------------------------------------------------------------------
void PorousLayerTest::Knudsen_diffusion(){
    
    std::vector<double> T_dummy(8, 298.0);
    FuelCellShop::SolutionVariable sols(T_dummy,temperature_of_REV);
    
    double D_k_check(1.480150e-5);
    
    std::vector<double> D_k;
    std::vector<double> dD_k_dT;
    FuelCellShop::Material::Oxygen oxygen;
    
    layer->compute_Knudsen_diffusion(&oxygen, sols, D_k, dD_k_dT);
    
    std::cout<<"Knudsen diffusion value for "<<oxygen.name_material()<<": "<<D_k[0]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2)<<" cm^2/s"<<std::endl;
    
    //Final check
    std::ostringstream streamOut;
    streamOut <<"PorousLayerTest::Knudsen_diffusion"<<std::endl;    
    
    for (unsigned int i=0; i<T_dummy.size(); ++i){
        streamOut <<"The value of the D_k is: "<<D_k[i]<<". The expected value is: "<<D_k_check<<std::endl;
        TEST_ASSERT_MSG(std::fabs((D_k_check - D_k[i])/D_k_check)<1e-6, streamOut.str().c_str());
        streamOut.flush();
    }
}

//---------------------------------------------------------------------------
void PorousLayerTest::Knudsen_diffusion_derivatives(){
    
    double T(298.0);
    std::vector<double> T_dummy(8, T);
    FuelCellShop::SolutionVariable sols(T_dummy,temperature_of_REV);
    
    double D_k_check(1.480150e-5);
    double dD_k_dT_check(0.5*D_k_check*(1/T));
    
    std::vector<double> D_k;
    std::vector<double> dD_k_dT;
    FuelCellShop::Material::Oxygen oxygen;
    
    layer->compute_Knudsen_diffusion(&oxygen, sols, D_k, dD_k_dT);
    
    //Final check
    std::ostringstream streamOut;
    streamOut <<"PorousLayerTest::Knudsen_diffusion_derivatives"<<std::endl;
    streamOut <<"The value of the D_k is: "<<D_k[0]<<". The expected value is: "<<D_k_check<<std::endl;
    TEST_ASSERT_MSG(std::fabs((D_k_check - D_k[0])/D_k_check)<1e-6, streamOut.str().c_str());    
    streamOut.flush();
    
    streamOut <<"The value of the dD_k_dT is: "<<dD_k_dT[0]<<". The expected value is: "<<dD_k_dT_check<<std::endl;
    TEST_ASSERT_MSG(std::fabs((dD_k_dT_check - dD_k_dT[0])/dD_k_dT_check)<1e-6, streamOut.str().c_str());
    streamOut.flush();
    
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void PorousLayerTest::diffusion_test()
{
    std::vector< FuelCellShop::Material::PureGas* > cathode_gases;
    cathode_gases.push_back(&oxygen);
    cathode_gases.push_back(&nitrogen);
    
    double pressure (1);  
    unsigned int n_quad = 3;
    FuelCellShop::SolutionVariable temperature = FuelCellShop::SolutionVariable(298.15, n_quad, temperature_of_REV);
    
    layer->set_gases (cathode_gases, pressure);

    layer->set_temperature( temperature );
    
    layer->compute_gas_diffusion(&oxygen,&nitrogen);
    
    std::vector<double> D_b;
    std::vector<double> dD_b_dT;
    layer->gas_diffusion_coefficient(D_b, dD_b_dT);
    
    double D_check(1/((1/1.480150e-5) + (1/2.0304e-5)));
    
        std::ostringstream streamOut;
    for (unsigned int i=0; i<n_quad; ++i){
        streamOut <<"The value of the D_b is: "<<D_b[i]<<". The expected value is: "<<D_check<<std::endl;
        TEST_ASSERT_MSG(std::fabs((D_check - D_b[i])/D_check)<1e-3, streamOut.str().c_str());
        streamOut.str(""); streamOut.clear();
    }
    
    // Derivative:
    double step = 1e-5;
    
    temperature = FuelCellShop::SolutionVariable(298.15+step, n_quad, temperature_of_REV);
    layer->set_temperature( temperature );
    
    layer->compute_gas_diffusion(&oxygen,&nitrogen);
    std::vector<double> Db_d;
    layer->gas_diffusion_coefficient(Db_d);
    
    for (unsigned int i=0; i<n_quad; ++i){
        double dDk_dT_expected = (Db_d[i]-D_b[i])/step;
        streamOut <<"The value of the dD_dT is: "<<dD_b_dT[i]<<". The expected value is: "<<dDk_dT_expected<<std::endl;
        TEST_ASSERT_MSG(std::fabs((dDk_dT_expected - dD_b_dT[i])/dDk_dT_expected)<1e-6, streamOut.str().c_str());
        streamOut.str(""); streamOut.clear();
    }
    
}