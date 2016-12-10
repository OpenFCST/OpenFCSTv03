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
//    - $Id: agglomerate_catalyst_layer_test.cc 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#include <agglomerate_catalyst_layer_test.h>

void
MultiScaleCLTest::setup()
{
    
    
    FuelCellShop::Layer::CatalystLayer<dim>::declare_CatalystLayer_parameters("Cathode catalyst layer", param);
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Cathode catalyst layer");
        {
    		param.set("Material id", "4");
            param.set("Catalyst layer type", "MultiScaleCL");

    		param.enter_subsection("ConventionalCL");
    		param.set("Platinum loading on support (%wt)", "4:0.46");
    		param.set("Platinum loading per unit volume (mg/cm3)", "4:400");
    		param.set("Electrolyte loading (%wt)", "4:0.30");
    		param.set("Active area [cm^2/cm^3]", "4:2.0e5");
    		param.leave_subsection();

            param.enter_subsection("MultiScaleCL");
                {
                param.enter_subsection("MicroScale");

                    param.set("Microscale type", "IonomerAgglomerateAnalytical");

                param.leave_subsection();
                }
                param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    layer = FuelCellShop::Layer::CatalystLayer<dim>::create_CatalystLayer("Cathode catalyst layer", param);

    layer->set_reaction_kinetics(ORR);
    layer->set_constant_solution(101325., total_pressure);
    layer->set_constant_solution(353., temperature_of_REV);
    
    
    std::vector<FuelCellShop::SolutionVariable> sols;
    
    FuelCellShop::SolutionVariable phi_s(1.5, 8, electronic_electrical_potential);
    
    std::vector<double> phiM_dummy(8, -0.2);
    FuelCellShop::SolutionVariable phi_m(phiM_dummy, protonic_electrical_potential);
    
    xO2_dummy = std::vector<double>(8, 0.5);
    FuelCellShop::SolutionVariable x_O2(&xO2_dummy, oxygen_molar_fraction);
    sols.push_back( phi_s );
    sols.push_back( phi_m);
    sols.push_back( x_O2 );
    
    layer->set_solution(sols);
}

void
MultiScaleCLTest::testSetSolution()
{
    TEST_ASSERT_DELTA_MSG( layer->solutions[electronic_electrical_potential][0], 1.5, 1.e-6,  "phi_s first element failure !!");
    TEST_ASSERT_MSG( layer->solutions[electronic_electrical_potential].size() == 8, "phi_s size failure !!");
    TEST_ASSERT_DELTA_MSG( layer->solutions[electronic_electrical_potential][7], 1.5, 1.e-6,  "phi_s last element failure !!");
    
    TEST_ASSERT_DELTA_MSG( layer->solutions[protonic_electrical_potential][0], -0.2, 1.e-6,  "phi_m first element failure !!");
    TEST_ASSERT_MSG( layer->solutions[protonic_electrical_potential].size() == 8, "phi_m size failure !!");
    TEST_ASSERT_DELTA_MSG( layer->solutions[protonic_electrical_potential][7], -0.2, 1.e-6,  "phi_m last element failure !!");
    
    TEST_ASSERT_DELTA_MSG( layer->solutions[oxygen_molar_fraction][0], 0.5, 1.e-6,  "x_O2 first element failure !!");
    TEST_ASSERT_MSG( layer->solutions[oxygen_molar_fraction].size() == 8, "x_O2 size failure !!");
    TEST_ASSERT_DELTA_MSG( layer->solutions[oxygen_molar_fraction][7], 0.5, 1.e-6,  "x_O2 last element failure !!");
    
    TEST_ASSERT_DELTA_MSG( layer->solutions[membrane_water_content][0], 12., 1.e-6,  "lambda first element failure !!");
    TEST_ASSERT_MSG( layer->solutions[membrane_water_content].size() == 8, "lambda size failure !!");
    TEST_ASSERT_DELTA_MSG( layer->solutions[membrane_water_content][7], 12., 1.e-6,  "lambda last element failure !!");
}
