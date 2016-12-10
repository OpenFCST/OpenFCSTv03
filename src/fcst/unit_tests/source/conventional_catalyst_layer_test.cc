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
//    - $Id: conventional_catalyst_layer_test.cc 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#include <conventional_catalyst_layer_test.h>

void
ConventionalCLTest::setup()
{
    FuelCellShop::Layer::CatalystLayer<dim>::declare_CatalystLayer_parameters("Cathode catalyst layer", param);
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Cathode catalyst layer");
        {
    		param.set("Material id", "4,5");
            param.set("Catalyst layer type", "HomogeneousCL");

            param.enter_subsection("ConventionalCL");
            {
                param.set("Platinum loading on support (%wt)","4:.46,5:.46");
                param.set("Platinum loading per unit volume (mg/cm3)","4:400,5:400");
                param.set("Method to compute porosity","NafionLoading");
                param.set("Electrolyte loading (%wt)","4:0.3,5:0.3");
                param.set("Ionomer to Carbon Ratio","4:1.0,5:1.0");
                param.set("Method relative liquid permeability", "Kumbur07");
                param.set("Irreducible liquid water saturation","4:0.0");
                param.set("Absolute permeability [cm^2]", "4:1.5e-9,5:1.5e-9");
                param.set("Method capillary pressure - saturation function", "Ye07");
        		param.set("Active area [cm^2/cm^3]", "4:2.0e5,5:1.0e5");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    layer = FuelCellShop::Layer::CatalystLayer<dim>::create_CatalystLayer("Cathode catalyst layer", param);
}

//---------------------------------------------
void 
ConventionalCLTest::testVolumeFraction()
{
    std::map< std::string, double > volume_fraction;

    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Cathode catalyst layer");
        {
            param.enter_subsection("ConventionalCL");
            {
                param.set("Platinum loading per unit volume (mg/cm3)","4:400,5:300");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();

    layer.reset();

    layer = FuelCellShop::Layer::CatalystLayer<dim>::create_CatalystLayer("Cathode catalyst layer", param);

    // Setting material ID to first CL (Material ID:4).
    layer->set_local_material_id(4);
    layer->get_volume_fractions (volume_fraction);

    double expectedAnswerSolidCCL1(0.253387);
    double expectedAnswerIonomerCCL1(0.186335);
    double expectedAnswerVoidCCL1(0.560277);

	TEST_ASSERT_DELTA_MSG(expectedAnswerSolidCCL1, volume_fraction.find("Solid")->second, 1.e-5, "Solid volume fraction failed!");
	TEST_ASSERT_DELTA_MSG(expectedAnswerVoidCCL1, volume_fraction.find("Void")->second, 1.e-5, "Void volume fraction failed!");
	TEST_ASSERT_DELTA_MSG(expectedAnswerIonomerCCL1, volume_fraction.find("Ionomer")->second, 1.e-5, "Ionomer volume fraction failed!");

	// Setting material ID to second CL (Material ID:5).
    layer->set_local_material_id(5);
    layer->get_volume_fractions (volume_fraction);

    double expectedAnswerSolidCCL2(0.19004044489383215);
    double expectedAnswerIonomerCCL2(0.13975155279503104);
    double expectedAnswerVoidCCL2(0.67020800231113675);

	TEST_ASSERT_DELTA_MSG(expectedAnswerSolidCCL2, volume_fraction.find("Solid")->second, 1.e-5, "Solid volume fraction failed!");
	TEST_ASSERT_DELTA_MSG(expectedAnswerVoidCCL2, volume_fraction.find("Void")->second, 1.e-5, "Void volume fraction failed!");
	TEST_ASSERT_DELTA_MSG(expectedAnswerIonomerCCL2, volume_fraction.find("Ionomer")->second, 1.e-5, "Ionomer volume fraction failed!");
}

//---------------------------------------------
void 
ConventionalCLTest::testICRatio()
{
   double expectedAnswerCCL1 = 0.6257142857;
   double expectedAnswerCCL2 = 0.55;

    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Cathode catalyst layer");
        {
    		param.set("Material id", "4,5");
            param.set("Catalyst layer type", "HomogeneousCL");

            param.enter_subsection("ConventionalCL");
            {
                param.set("Method to compute porosity","ICRatio");
                param.set("Ionomer to Carbon Ratio","4: 0.6257142857,5: 0.55");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    layer.reset();
    
    layer = FuelCellShop::Layer::CatalystLayer<dim>::create_CatalystLayer("Cathode catalyst layer", param);

    std::map<std::string, double> info;
    std::map< std::string, double > volume_fraction;

    // Setting material ID to first CL (Material ID:4).
    layer->set_local_material_id(4);
    layer->get_loadings(info);
    layer->get_volume_fractions (volume_fraction);
    
    
    double answerCCL1 = info.find("IC_ratio")->second;
    TEST_ASSERT_MSG(expectedAnswerCCL1 == answerCCL1, "ICRatio failed!");
    
    double expectedAnswerSolidCCL1(0.253387);
    double expectedAnswerIonomerCCL1(0.186335);
    double expectedAnswerVoidCCL1(0.560277);

    TEST_ASSERT_DELTA_MSG(expectedAnswerSolidCCL1, volume_fraction.find("Solid")->second, 1.e-5, "Solid volume fraction failed!");
    TEST_ASSERT_DELTA_MSG(expectedAnswerVoidCCL1, volume_fraction.find("Void")->second, 1.e-5, "Void volume fraction failed!");
    TEST_ASSERT_DELTA_MSG(expectedAnswerIonomerCCL1, volume_fraction.find("Ionomer")->second, 1.e-5, "Ionomer volume fraction failed!");


    // Setting material ID to second CL (Material ID:5).
    layer->set_local_material_id(5);
    layer->get_loadings(info);
    layer->get_volume_fractions (volume_fraction);


    double answerCCL2 = info.find("IC_ratio")->second;
    TEST_ASSERT_MSG(expectedAnswerCCL2 == answerCCL2, "ICRatio failed!");

    double expectedAnswerSolidCCL2(0.2533872598584429);
    double expectedAnswerIonomerCCL2(0.16378796902918402);
    double expectedAnswerVoidCCL2(0.58282477111237307);

    TEST_ASSERT_DELTA_MSG(expectedAnswerSolidCCL2, volume_fraction.find("Solid")->second, 1.e-5, "Solid volume fraction failed!");
    TEST_ASSERT_DELTA_MSG(expectedAnswerVoidCCL2, volume_fraction.find("Void")->second, 1.e-5, "Void volume fraction failed!");
    TEST_ASSERT_DELTA_MSG(expectedAnswerIonomerCCL2, volume_fraction.find("Ionomer")->second, 1.e-5, "Ionomer volume fraction failed!");

    
}

//---------------------------------------------
void 
ConventionalCLTest::testNafionLoading()
{
      double expectedAnswerCCL1 = 0.3;
      double expectedAnswerCCL2 = 0.2;
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Cathode catalyst layer");
        {
    		param.set("Material id", "4,5");
            param.set("Catalyst layer type", "HomogeneousCL");

            param.enter_subsection("ConventionalCL");
            {
                param.set("Method to compute porosity","NafionLoading");
                param.set("Electrolyte loading (%wt)","4:0.3,5:0.2");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    layer.reset();
    
    layer = FuelCellShop::Layer::CatalystLayer<dim>::create_CatalystLayer("Cathode catalyst layer", param);

    std::map<std::string, double> info;
    std::map< std::string, double > volume_fraction;

    // Setting material ID to first CL (Material ID:4).
    layer->set_local_material_id(4);
    layer->get_loadings(info);
    layer->get_volume_fractions (volume_fraction);
    
    double answerCCL1 = info.find("loading_N")->second;
    TEST_ASSERT_MSG(expectedAnswerCCL1 == answerCCL1, "ICRatio failed!");
    
    double expectedAnswerSolidCCL1(0.253387);
    double expectedAnswerIonomerCCL1(0.186335);
    double expectedAnswerVoidCCL1(0.560277);
    
    TEST_ASSERT_DELTA_MSG(expectedAnswerSolidCCL1, volume_fraction.find("Solid")->second, 1.e-5, "Solid volume fraction failed!");
    TEST_ASSERT_DELTA_MSG(expectedAnswerVoidCCL1, volume_fraction.find("Void")->second, 1.e-5, "Void volume fraction failed!");
    TEST_ASSERT_DELTA_MSG(expectedAnswerIonomerCCL1, volume_fraction.find("Ionomer")->second, 1.e-5, "Ionomer volume fraction failed!");

    // Setting material ID to second CL (Material ID:5).
    layer->set_local_material_id(5);
    layer->get_loadings(info);
    layer->get_volume_fractions (volume_fraction);

    double answerCCL2 = info.find("loading_N")->second;
    TEST_ASSERT_MSG(expectedAnswerCCL2 == answerCCL2, "ICRatio failed!");

    double expectedAnswerSolidCCL2(0.2533872598584429);
    double expectedAnswerIonomerCCL2(0.10869565217391304);
    double expectedAnswerVoidCCL2(0.63791708796764401);

    TEST_ASSERT_DELTA_MSG(expectedAnswerSolidCCL2, volume_fraction.find("Solid")->second, 1.e-5, "Solid volume fraction failed!");
    TEST_ASSERT_DELTA_MSG(expectedAnswerVoidCCL2, volume_fraction.find("Void")->second, 1.e-5, "Void volume fraction failed!");
    TEST_ASSERT_DELTA_MSG(expectedAnswerIonomerCCL2, volume_fraction.find("Ionomer")->second, 1.e-5, "Ionomer volume fraction failed!");
}

//---------------------------------------------
void
ConventionalCLTest::testActiveArea()
{
    // Setting material ID to first CL (Material ID:4).
	layer->set_local_material_id(4);
	double expectedAnswerActiveArea(2.0e5);
	TEST_ASSERT_DELTA_MSG(expectedAnswerActiveArea, layer->get_active_area_Pt (), 1.e-5, "Active area comparison failed!");

    // Setting material ID to second CL (Material ID:5).
	layer->set_local_material_id(5);

	double expectedAnswerActiveAreaCCL2(1.0e5);
	TEST_ASSERT_DELTA_MSG(expectedAnswerActiveAreaCCL2, layer->get_active_area_Pt (), 1.e-5, "Active area comparison failed!");
}

//---------------------------------------------
void
ConventionalCLTest::testPtLoadingSupport()
{
    std::map<std::string, double> info;
	double expectedAnswerPtCarbonSupportCCL1(0.46);
	double expectedAnswerPtCarbonSupportCCL2(0.2);

    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Cathode catalyst layer");
        {
            param.enter_subsection("ConventionalCL");
            {
                param.set("Platinum loading on support (%wt)","4:.46,5:.2");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();

    layer.reset();

    layer = FuelCellShop::Layer::CatalystLayer<dim>::create_CatalystLayer("Cathode catalyst layer", param);


    // Setting material ID to first CL (Material ID:4).
	layer->set_local_material_id(4);
    layer->get_loadings(info);

    double answerCCL1 = info.find("prc_Pt")->second;
    TEST_ASSERT_MSG(expectedAnswerPtCarbonSupportCCL1 == answerCCL1, "prc_Pt failed!");

    // Setting material ID to second CL (Material ID:5).
	layer->set_local_material_id(5);
    layer->get_loadings(info);

    double answerCCL2 = info.find("prc_Pt")->second;
    TEST_ASSERT_MSG(expectedAnswerPtCarbonSupportCCL2 == answerCCL2, "prc_Pt failed!");
}

//---------------------------------------------
void
ConventionalCLTest::testKumburLiquidPermeability()
{
    std::vector<double> s(3);
    std::vector<VariableNames> deriv_flags;
    s[0] = -0.1; s[1] = 0.; s[2] = 0.5;
    deriv_flags.push_back(liquid_water_saturation);
    layer->set_local_material_id(4);
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
    TEST_ASSERT_DELTA_MSG(k_l[2][1][1], 3.356344016e-10, 1.e-12, "k_l[2][1][1] failure!");
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
ConventionalCLTest::testInterfacialArea()
{
    std::vector<double> s(5);
    std::vector<VariableNames> deriv_flags;
    s[0] = -0.1; s[1] = 0.; s[2] = 0.5; s[3] = 1.0; s[4] = 1.1;
    deriv_flags.push_back(liquid_water_saturation);
    
    double delta(std::pow(10.0, -10));
    std::vector<double> delta_s(5);
    delta_s[0] = s[0] + delta; delta_s[1] = s[1] + delta; delta_s[2] = s[2] + delta; delta_s[3] = s[3] - delta; delta_s[4] = s[4] + delta;
    // ------------------------

    layer->set_saturation( FuelCellShop::SolutionVariable(&s, liquid_water_saturation) );
    
    std::vector<double> a_lv;
    layer->interfacial_surface_area(a_lv);
    
    TEST_ASSERT_DELTA_MSG(a_lv[0], 0.0, 1.e-6, "a_lv[0] failure!");
    TEST_ASSERT_DELTA_MSG(a_lv[1], 0.0, 1.e-6, "a_lv[1] failure!");
    TEST_ASSERT_DELTA_MSG(a_lv[2], 3270.31956625748, 1.e-6, "a_lv[2] failure!");
    TEST_ASSERT_DELTA_MSG(a_lv[3], 0.0, 1.e-6, "a_lv[3] failure!");
    TEST_ASSERT_DELTA_MSG(a_lv[4], 0.0, 1.e-6, "a_lv[4] failure!");
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
    TEST_ASSERT_DELTA_MSG(deriv_a_lv[liquid_water_saturation][2], (a_lv_delta[2]-a_lv[2])/delta, 1.e-2, "dalv_ds[2] failure!");
    TEST_ASSERT_DELTA_MSG(deriv_a_lv[liquid_water_saturation][3], (a_lv[3] - a_lv_delta[3])/delta, 1.e-6, "dalv_ds[3] failure!");
    TEST_ASSERT_DELTA_MSG(deriv_a_lv[liquid_water_saturation][4], (a_lv_delta[4]-a_lv[4])/delta, 1.e-6, "dalv_ds[4] failure!");    
}

void
ConventionalCLTest::testPCapillary()
{
    std::vector<double> s(3);
    std::vector<double> T(3);
    std::vector<VariableNames> deriv_flags;
    s[0] = -0.1; s[1] = 0.0; s[2] = 0.2;
    T[0] = 350.0; T[1] = 360.0; T[2] = 370.0;
    deriv_flags.push_back(liquid_water_saturation); deriv_flags.push_back(temperature_of_REV);

    double delta(std::pow(10.0, -7));
    std::vector<double> delta_s(3);
    std::vector<double> delta_T(3);
    delta_s[0] = s[0] + delta; delta_s[1] = s[1] + delta; delta_s[2] = s[2] + delta;
    delta_T[0] = T[0] + delta; delta_T[1] = T[1] + delta; delta_T[2] = T[2] + delta;
    // ------------------------
    
    layer->set_saturation( FuelCellShop::SolutionVariable(&s, liquid_water_saturation) );
    layer->set_temperature( FuelCellShop::SolutionVariable(&T, temperature_of_REV) );
    
    std::vector<double> pc;
    layer->pcapillary(pc);
    
    TEST_ASSERT_DELTA_MSG(pc[0], -481.600292127737, 1.e-6, "pcapillary[0] failure!");
    TEST_ASSERT_DELTA_MSG(pc[1], -481.600292127737, 1.e-6, "pcapillary[1] failure!");
    TEST_ASSERT_DELTA_MSG(pc[2], -438.6384930866, 1.e-6, "pcapillary[2] failure!");
    // ------------------------
    
    std::vector<double>  dpc_ds;
    layer->dpcapillary_dsat(dpc_ds);
    
    layer->set_saturation( FuelCellShop::SolutionVariable(&delta_s, liquid_water_saturation) );
    
    std::vector<double>  pc_delta;
    layer->pcapillary(pc_delta);
    
    TEST_ASSERT_DELTA_MSG(dpc_ds[0], (pc_delta[0]-pc[0])/delta, 1.e-6, "dpc_ds[0] failure!");
    TEST_ASSERT_DELTA_MSG(dpc_ds[1], (pc_delta[1]-pc[1])/delta, 1.e-5, "dpc_ds[1] failure!");
    TEST_ASSERT_DELTA_MSG(dpc_ds[2], (pc_delta[2]-pc[2])/delta, 1.e-5, "dpc_ds[2] failure!");
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
    TEST_ASSERT_DELTA_MSG(deriv_dpc_ds[liquid_water_saturation][1], (dpc_ds_delta[1]-dpc_ds[1])/delta, 1.e-6, "d2pc_ds2[1] failure!");
    TEST_ASSERT_DELTA_MSG(deriv_dpc_ds[liquid_water_saturation][2], (dpc_ds_delta[2]-dpc_ds[2])/delta, 1.e-6, "d2pc_ds2[2] failure!");
    
    layer->set_saturation( FuelCellShop::SolutionVariable(&s, liquid_water_saturation) );
    layer->set_temperature( FuelCellShop::SolutionVariable(&delta_T, temperature_of_REV) );
    
    layer->dpcapillary_dsat(dpc_ds_delta);
    
    TEST_ASSERT_DELTA_MSG(deriv_dpc_ds[temperature_of_REV][0], (dpc_ds_delta[0]-dpc_ds[0])/delta, 1.e-6, "d2pc_dsdT[0] failure!");
    TEST_ASSERT_DELTA_MSG(deriv_dpc_ds[temperature_of_REV][1], (dpc_ds_delta[1]-dpc_ds[1])/delta, 1.e-6, "d2pc_dsdT[1] failure!");
    TEST_ASSERT_DELTA_MSG(deriv_dpc_ds[temperature_of_REV][2], (dpc_ds_delta[2]-dpc_ds[2])/delta, 1.e-6, "d2pc_dsdT[2] failure!");
}
