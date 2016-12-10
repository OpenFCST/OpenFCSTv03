//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: PSD_test.cc
//    - Description: Unit testing class for PSD
//    - Developers: Prafful Mangal
//    - Id: $Id: PSD_HO_test.cc 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#include "PSD_HO_test.h"

//-------------------------------------------------------------
void PSD_HO_Test::setup()
{
    ParameterHandler param;
    
    boost::shared_ptr<FuelCellShop::Layer::GasDiffusionLayer<dim> > CGDL;
    
    FuelCellShop::Layer::GasDiffusionLayer<dim>::declare_GasDiffusionLayer_parameters("Cathode gas diffusion layer", param);
    
//     param.enter_subsection("Fuel cell data");
//     {
//         param.enter_subsection("Cathode gas diffusion layer");
//         {
//             param.enter_subsection("Generic data");
//             {
//                 param.set("Porosity", "0.24");
//             }
//             param.leave_subsection();
//         }
//         param.leave_subsection();
//     }
//     param.leave_subsection();
//     
//     CGDL = FuelCellShop::Layer::GasDiffusionLayer<dim>::create_GasDiffusionLayer("Cathode gas diffusion layer",param);
//     
//     psd_object.declare_parameters(param);
//     
//     param.enter_subsection("PSD parameters");
//     {
//         param.enter_subsection("BasePSD");
//         {
// 
//             //param.set("porosity", "0.24");
//             param.set("Gamma", "0.0728");
//             param.set("Contact angle", "1.396");
//             param.set("lambda", "1.0");
//             param.set("Volume fraction Hydrophobic", "0.3");
//             param.set("psd type", "HOPSD"); 
//             
//             param.set("Mode probability global", "0.31, 0.18, 0.10, 0.05, 0.36");
//             param.set("Mode characteristic radius global", "38.2, 125.0, 500.0, 1200.0, 2000.0");
//             param.set("Mode width global", "0.43, 0.55, 0.65, 0.45, 0.75");
//             
//             param.enter_subsection("HOPSD");
//             {
//                 param.set("capillay pressure", "1010000");
//                 param.set("Hydrophobic Mode probability global", "0.31, 0.18, 0.10, 0.05, 0.36");
//                 param.set("Hydrophobic Mode characteristic radius global", "38.2, 125.0, 500.0, 1200.0, 2000.0");
//                 param.set("Hydrophobic Mode width global", "0.43, 0.55, 0.65, 0.45, 0.75");
//                 
//             } 
//             
//             
//             param.leave_subsection(); 
//             
//         }
//         param.leave_subsection();
//     }
//     param.leave_subsection();
//     
//     psd_object.initialize(param);
//     
//     psd_object.set_porosity (CGDL->get_porosity());
}

// ------------------------------------------------------------
// void PSD_HO_Test::testcompute_rc_HO() 
// {
//     std::vector<double> answer(0.0);
//     
//     psd_object.get_critical_radius(answer);   
//     double expectedAnswer = 25.07024e-3;
//     
//     std::ostringstream streamOut;
//     streamOut <<"The value of the rc_HO (nm) is: "<<answer[0]<<". The expected value is: "<<expectedAnswer<<std::endl;
//     TEST_ASSERT_DELTA_MSG(expectedAnswer, answer[0], 1e-7, streamOut.str().c_str());
// }

// //-------------------------------------------------------------
// void PSD_HO_Test::testcompute_k_sat_HO() 
// {
//     double answer(0.0);
//     
//     psd_object.get_global_saturated_permeability(answer);   
//     double expectedAnswer = 3802.38608;
//     
//     std::ostringstream streamOut;
//     streamOut <<"The value of the k_sat_HO (nm^2) is: "<<answer<<". The expected value is: "<<expectedAnswer<<std::endl;
//     TEST_ASSERT_DELTA_MSG(expectedAnswer, answer, 1e-8, streamOut.str().c_str());
// }
// 
// //-------------------------------------------------------------
// void PSD_HO_Test::testcompute_sat_HO() 
// {
//     std::vector<double> answer(0.0);
//     
//     psd_object.get_saturation(answer);  
//     double expectedAnswer = 0.2856832254;
//     
//     std::ostringstream streamOut;
//     streamOut <<"The value of the Saturation_HO is: "<<answer[0]<<". The expected value is: "<<expectedAnswer<<std::endl;
//     TEST_ASSERT_DELTA_MSG(expectedAnswer, answer[0], 1.1e-3, streamOut.str().c_str());
//     
// }
// 
// //-------------------------------------------------------------
// void PSD_HO_Test::testcompute_k_L_HO() 
// {
//     std::vector<double> answer(0.0);
//     
//     psd_object.get_pore_HO_liquid_saturated_permeability(answer);   
//     double expectedAnswer = 92.44698;
//     
//     std::ostringstream streamOut;
//     streamOut <<"The value of the K_L_HO(nm^2) is: "<<answer[0]<<". The expected value is: "<<expectedAnswer<<std::endl;
//     TEST_ASSERT_DELTA_MSG(expectedAnswer, answer[0], 1e-4, streamOut.str().c_str());
//     
// }
// 
// //-------------------------------------------------------------
// void PSD_HO_Test::testcompute_kr_L_HO() 
// {
//     std::vector<double> answer(0.0);
//     
//     psd_object.get_relative_liquid_permeability(answer);  
//     double expectedAnswer = 0.0243128858;
//     
//     std::ostringstream streamOut;
//     streamOut <<"The value of the Kr_L_HO is: "<<answer[0]<<". The expected value is: "<<expectedAnswer<<std::endl;
//     TEST_ASSERT_DELTA_MSG(expectedAnswer, answer[0], 1e-7, streamOut.str().c_str());
//     
// }
// 
// //-------------------------------------------------------------
// void PSD_HO_Test::testcompute_k_G_HO() 
// {
//     std::vector<double> answer(0.0);
//     
//     psd_object.get_pore_HO_gas_saturated_permeability(answer);  
//     double expectedAnswer = 0.0114221467;
//     
//     std::ostringstream streamOut;
//     streamOut <<"The value of the K_G_HO(nm^2) is: "<<answer[0]<<". The expected value is: "<<expectedAnswer<<std::endl;
//     TEST_ASSERT_DELTA_MSG(expectedAnswer, answer[0], 1e-7, streamOut.str().c_str());
// }
// 
// //-------------------------------------------------------------
// void PSD_HO_Test::testcompute_kr_G_HO() 
// {
//     std::vector<double> answer(0.0);
//     
//     psd_object.get_relative_gas_permeability(answer);   
//     double expectedAnswer = 3.00394185e-6;
//     
//     std::ostringstream streamOut;
//     streamOut <<"The value of the Kr_G_HO(nm^2) is: "<<answer[0]<<". The expected value is: "<<expectedAnswer<<std::endl;
//     TEST_ASSERT_DELTA_MSG(expectedAnswer, answer[0], 1e-10, streamOut.str().c_str());
// }
// 
// //-------------------------------------------------------------
// void PSD_HO_Test::testcompute_interfacial_area_per_volume_HO() 
// {
//     std::vector<double> answer(0.0);
//     
//     psd_object.get_liquid_gas_interfacial_surface(answer);  
//     double expectedAnswer = 0.0001166379;
//     
//     std::ostringstream streamOut;
//     streamOut <<"The value of the interfacial area per unit volume for hydrophobic pores is: "<<answer[0]<<". The expected value is: "<<expectedAnswer<<std::endl;
//     TEST_ASSERT_DELTA_MSG(expectedAnswer, answer[0], 1e-5, streamOut.str().c_str());
// }
// 
// //-------------------------------------------------------------
// void PSD_HO_Test::testcompute_avg_rk_HO() 
// {
//     std::vector<double> answer(0.0);
//     
//     psd_object.get_knudsen_radius(answer);  
//     double expectedAnswer = 9.78731659444e-9;
//     
//     std::ostringstream streamOut;
//     streamOut <<"The value of average knudsen radius (meter) for hydrophobic pores is: "<<answer[0]<<". The expected value is: "<<expectedAnswer<<std::endl;
//     TEST_ASSERT_DELTA_MSG(expectedAnswer, answer[0], 1e-13, streamOut.str().c_str());
//     
// }


//-------------------------------------------------------------
//-------------------------------------------------------------
