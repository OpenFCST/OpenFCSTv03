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
//    - $Id: design_mpl_test.cc 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------
#include <design_mpl_test.h>

void
DesignMPLTest::setup()
{
    FuelCellShop::Layer::MicroPorousLayer<dim>::declare_MicroPorousLayer_parameters("Cathode microporous layer", param);
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Cathode microporous layer");
        {
            param.set("Micro porous layer type", "DesignMPL");

            param.enter_subsection("DesignMPL");
            {
                param.set("Porosity", "0.4");
                param.set("Anisotropic transport", "true");
                param.set("Method relative liquid permeability", "Kumbur07");
                param.set("Irreducible liquid water saturation","0.0");
                param.set("Absolute permeability X [cm^2]", "1.5e-9");
                param.set("Absolute permeability Y [cm^2]", "5.0e-9");
                param.set("Method capillary pressure - saturation function", "Kumbur07-corrected");
                param.set("Compaction pressure [MPa]", "0.0");
                param.set("PTFE loading [% wt]", "5.0");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    layer = FuelCellShop::Layer::MicroPorousLayer<dim>::create_MicroPorousLayer("Cathode microporous layer", param);
}

void
DesignMPLTest::testKumburLiquidPermeability()
{
    std::vector<double> s(3);
    std::vector<VariableNames> deriv_flags;
    s[0] = -0.1; s[1] = 0.; s[2] = 0.5;
    deriv_flags.push_back(liquid_water_saturation);
    
    double delta(std::pow(10.0, -7));
    std::vector<double> delta_s(3);
    delta_s[0] = s[0] + delta; delta_s[1] = s[1] + delta; delta_s[2] = s[2] + delta;
    // ------------------------
    
    layer->set_saturation( FuelCellShop::SolutionVariable(&s, liquid_water_saturation) );
    
    std::vector< Tensor<2,dim> > k_l;
    layer->liquid_permeablity(k_l);
    
    TEST_ASSERT_DELTA_MSG(k_l[0][0][0], 0., 1.e-6, "k_l[0][0][0] failure!");
    TEST_ASSERT_DELTA_MSG(k_l[0][1][1], 0., 1.e-6, "k_l[0][1][1] failure!");
    TEST_ASSERT_DELTA_MSG(k_l[1][0][0], 0., 1.e-6, "k_l[1][0][0] failure!");
    TEST_ASSERT_DELTA_MSG(k_l[1][1][1], 0., 1.e-6, "k_l[1][1][1] failure!");
    TEST_ASSERT_DELTA_MSG(k_l[2][0][0], 3.356344016e-10, 1.e-12, "k_l[2][0][0] failure!");
    TEST_ASSERT_DELTA_MSG(k_l[2][1][1], 1.118781339e-9, 1.e-12, "k_l[2][1][1] failure!");
    // ------------------------
    
    layer->set_derivative_flags(deriv_flags);
    
    std::map< VariableNames, std::vector< Tensor<2,dim> > > deriv_k_l;
    layer->derivative_liquid_permeablity(deriv_k_l);
    
    layer->set_saturation( FuelCellShop::SolutionVariable(&delta_s, liquid_water_saturation) );
    std::vector< Tensor<2,dim> > k_l_delta;
    layer->liquid_permeablity(k_l_delta);
    
    TEST_ASSERT_DELTA_MSG(deriv_k_l[liquid_water_saturation][0][0][0], (k_l_delta[0][0][0] - k_l[0][0][0])/delta, 1.e-8, "deriv_k_l[liquid_water_saturation][0][0][0] failure!");
    TEST_ASSERT_DELTA_MSG(deriv_k_l[liquid_water_saturation][0][1][1], (k_l_delta[0][1][1] - k_l[0][1][1])/delta, 1.e-8, "deriv_k_l[liquid_water_saturation][0][1][1] failure!");
    TEST_ASSERT_DELTA_MSG(deriv_k_l[liquid_water_saturation][1][0][0], (k_l_delta[1][0][0] - k_l[1][0][0])/delta, 1.e-8, "deriv_k_l[liquid_water_saturation][1][0][0] failure!");
    TEST_ASSERT_DELTA_MSG(deriv_k_l[liquid_water_saturation][1][1][1], (k_l_delta[1][1][1] - k_l[1][1][1])/delta, 1.e-8, "deriv_k_l[liquid_water_saturation][1][1][1] failure!");
    TEST_ASSERT_DELTA_MSG(deriv_k_l[liquid_water_saturation][2][0][0], (k_l_delta[2][0][0] - k_l[2][0][0])/delta, 1.e-8, "deriv_k_l[liquid_water_saturation][2][0][0] failure!");
    TEST_ASSERT_DELTA_MSG(deriv_k_l[liquid_water_saturation][2][1][1], (k_l_delta[2][1][1] - k_l[2][1][1])/delta, 1.e-8, "deriv_k_l[liquid_water_saturation][2][1][1] failure!");
}

void
DesignMPLTest::testInterfacialArea()
{
    std::vector<double> s(3);
    std::vector<VariableNames> deriv_flags;
    s[0] = -0.1; s[1] = 0.; s[2] = 0.5;
    deriv_flags.push_back(liquid_water_saturation);
    
    double delta(std::pow(10.0, -10));
    std::vector<double> delta_s(5);
    delta_s[0] = s[0] + delta; delta_s[1] = s[1] + delta; delta_s[2] = s[2] + delta;
    // ------------------------
    
    layer->set_saturation( FuelCellShop::SolutionVariable(&s, liquid_water_saturation) );
    
    std::vector<double> a_lv;
    layer->interfacial_surface_area(a_lv);
    
    TEST_ASSERT_DELTA_MSG(a_lv[0], 0.0, 1.e-6, "a_lv[0] failure!");
    TEST_ASSERT_DELTA_MSG(a_lv[1], 0.0, 1.e-6, "a_lv[1] failure!");
    TEST_ASSERT_DELTA_MSG(a_lv[2], 1182.51058395054, 1.e-6, "a_lv[2] failure!");
    // ------------------------

    layer->set_derivative_flags(deriv_flags);
    
    std::map< VariableNames, std::vector<double> > deriv_a_lv;
    layer->derivative_interfacial_surface_area(deriv_a_lv);
    
    layer->set_saturation( FuelCellShop::SolutionVariable(&s, liquid_water_saturation) );
    layer->interfacial_surface_area(a_lv);
    
    layer->set_saturation( FuelCellShop::SolutionVariable(&delta_s, liquid_water_saturation) );
    
    std::vector<double> a_lv_delta;
    layer->interfacial_surface_area(a_lv_delta);
    
    TEST_ASSERT_DELTA_MSG(deriv_a_lv[liquid_water_saturation][0], (a_lv_delta[0]-a_lv[0])/delta, 1.e-6, "dalv_ds[0] failure!");
    TEST_ASSERT_DELTA_MSG(deriv_a_lv[liquid_water_saturation][1], (a_lv_delta[1]-a_lv[1])/delta, 1.e-1, "dalv_ds[1] failure!");
    TEST_ASSERT_DELTA_MSG(deriv_a_lv[liquid_water_saturation][2], (a_lv_delta[2]-a_lv[2])/delta, 1.e-2, "dalv_ds[2] failure!");  
}

void
DesignMPLTest::testPCapillary()
{
    std::vector<double> s(3);
    std::vector<double> T(3);
    std::vector<VariableNames> deriv_flags;
    s[0] = 0.0; s[1] = 1.0e-4; s[2] = 0.2;
    T[0] = 350.0; T[1] = 360.0; T[2] = 370.0;
    deriv_flags.push_back(liquid_water_saturation); deriv_flags.push_back(temperature_of_REV);
    
    double delta(std::pow(10.0, -10));
    std::vector<double> delta_s(3);
    std::vector<double> delta_T(3);
    delta_s[0] = s[0] + delta; delta_s[1] = s[1] + delta; delta_s[2] = s[2] + delta;
    // ------------------------
    
    layer->set_saturation( FuelCellShop::SolutionVariable(&s, liquid_water_saturation) );
    layer->set_temperature( FuelCellShop::SolutionVariable(&T, temperature_of_REV) );
    
    std::vector<double> pc;
    layer->pcapillary(pc);
    
    TEST_ASSERT_DELTA_MSG(pc[0], 594728.131286015, 1.e-6, "pcapillary[0] failure!");
    TEST_ASSERT_DELTA_MSG(pc[1], 487913.781685798, 1.e-6, "pcapillary[1] failure!");
    TEST_ASSERT_DELTA_MSG(pc[2], 2093778.09637407, 1.e-6, "pcapillary[2] failure!");
    // ------------------------
    
    std::vector<double>  dpc_ds;
    layer->dpcapillary_dsat(dpc_ds);
    
    layer->set_saturation( FuelCellShop::SolutionVariable(&delta_s, liquid_water_saturation) );
    
    std::vector<double>  pc_delta;
    layer->pcapillary(pc_delta);
    
    TEST_ASSERT_DELTA_MSG(dpc_ds[0], (pc_delta[0]-pc[0])/delta, 1.e-6, "dpc_ds[0] failure!");
    //TEST_ASSERT_DELTA_MSG(dpc_ds[1], (pc_delta[1]-pc[1])/delta, 10000., streamOut.str().c_str());
    TEST_ASSERT_DELTA_MSG(dpc_ds[2], (pc_delta[2]-pc[2])/delta, 10., "dpc_ds[2] failure!");
    // ------------------------
    
    layer->set_saturation( FuelCellShop::SolutionVariable(&s, liquid_water_saturation) );
    layer->set_temperature( FuelCellShop::SolutionVariable(&T, temperature_of_REV) );
    layer->set_derivative_flags(deriv_flags);
    
    std::map< VariableNames, std::vector<double> > deriv_dpc_ds;
    layer->derivative_dpcapillary_dsat(deriv_dpc_ds);
    
    layer->set_saturation( FuelCellShop::SolutionVariable(&delta_s, liquid_water_saturation) );
    
    std::vector<double>  dpc_ds_delta;
    layer->dpcapillary_dsat(dpc_ds_delta);
     
    TEST_ASSERT_DELTA_MSG(deriv_dpc_ds[liquid_water_saturation][0], (dpc_ds_delta[0]-dpc_ds[0])/delta, 1.e-6, "d2pc_ds2[0] failure!");
    //TEST_ASSERT_DELTA_MSG(deriv_dpc_ds[liquid_water_saturation][1], (dpc_ds_delta[1]-dpc_ds[1])/delta, 1., "d2pc_ds2[1] failure!");
    TEST_ASSERT_DELTA_MSG(deriv_dpc_ds[liquid_water_saturation][2], (dpc_ds_delta[2]-dpc_ds[2])/delta, 100., "d2pc_ds2[2] failure!");
    
    delta = std::pow(10.0, -6);
    delta_T[0] = T[0] + delta; delta_T[1] = T[1] + delta; delta_T[2] = T[2] + delta;
    
    layer->set_saturation( FuelCellShop::SolutionVariable(&s, liquid_water_saturation) );
    layer->set_temperature( FuelCellShop::SolutionVariable(&delta_T, temperature_of_REV) );
    
    layer->dpcapillary_dsat(dpc_ds_delta);
    
    TEST_ASSERT_DELTA_MSG(deriv_dpc_ds[temperature_of_REV][0], (dpc_ds_delta[0]-dpc_ds[0])/delta, 1.e-6, "d2pc_dsdT[0] failure!");
    TEST_ASSERT_DELTA_MSG(deriv_dpc_ds[temperature_of_REV][1], (dpc_ds_delta[1]-dpc_ds[1])/delta, 1., "d2pc_dsdT[1] failure!" );
    TEST_ASSERT_DELTA_MSG(deriv_dpc_ds[temperature_of_REV][2], (dpc_ds_delta[2]-dpc_ds[2])/delta, 1.e-2, "d2pc_dsdT[2] failure!");
}
