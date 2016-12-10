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
//    - Id: $Id: PSD_HI_test.cc 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#include "PSD_HI_test.h"
#include <boost/concept_check.hpp>

//-------------------------------------------------------------
void PSD_HI_Test::setup()
{
    ParameterHandler param;
    
    boost::shared_ptr<FuelCellShop::Layer::GasDiffusionLayer<dim> > CGDL;
    
    FuelCellShop::Layer::GasDiffusionLayer<dim>::declare_GasDiffusionLayer_parameters("Cathode gas diffusion layer", param);
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Cathode gas diffusion layer");
        {
            
            psd_pointer = &psd_object;
            
            psd_pointer->declare_parameters(param);
            
            param.enter_subsection("PSD parameters");
            {
                param.enter_subsection("BasePSD");
                {
                    
                    //param.set("porosity", "0.84");
                    param.set("Gamma", "0.0728");
                   // param.set("Contact angle", "100");
                    param.set("lambda", "1.0");
                    param.set("Volume fraction Hydrophilic", "0.5");
                    param.set("Volume fraction Hydrophobic", "0.5");
                    param.set("Mode probability global", "0.18, 0.1, 0.72, 0.62, 0.36, 0.02");
                    param.set("Mode characteristic radius global", "500e-9, 1200e-9, 2000e-9, 38.2e-9, 125e-9, 500e-9");
                    param.set("Mode width global", "0.65, 0.45, 0.75, 0.43, 0.55, 0.65");
                    param.set("psd type", "DualPSD"); 
                    
                    param.enter_subsection("HIPSD");
                    {
                        param.set("Hydrophilic Mode probability global", "0.62, 0.36, 0.02");
                        param.set("Hydrophilic Mode characteristic radius global", "38.2e-9, 125e-9, 500e-9");
                        param.set("Hydrophilic Mode width global", "0.43, 0.55, 0.65");
                        param.set("capillay pressure", "1000.0"); 
                        param.set("Static Contact Angle HI", "80.0"); 
                        
                    }
                    param.leave_subsection(); 
                    
                    param.enter_subsection("HOPSD");
                    {
                        param.set("capillay pressure", "1000.0");
                        param.set("Hydrophobic Mode probability global", "0.18, 0.1, 0.72");
                        param.set("Hydrophobic Mode characteristic radius global", "500e-9, 1200e-9, 2000e-9");
                        param.set("Hydrophobic Mode width global", "0.65, 0.45, 0.75");
                        param.set("Static Contact Angle HO", "100.0"); 
                    } 
                    param.leave_subsection(); 
                    
                    param.enter_subsection("DualPSD");
                    {
                        param.set("capillay pressure", "1000.0");
                    } 
                    param.leave_subsection(); 
                    
                }
                param.leave_subsection();
            }
            param.leave_subsection();
            
            param.enter_subsection("Generic data");
            {
                param.set("Porosity", "0.4");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Cathode gas diffusion layer");
        {
            psd_pointer->initialize(param);
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    CGDL = FuelCellShop::Layer::GasDiffusionLayer<dim>::create_GasDiffusionLayer("Cathode gas diffusion layer",param);
    
    psd_pointer->set_porosity (0.4);
    
}

//-------------------------------------------------------------
// void PSD_HI_Test::testcompute_rc_HI()
// {
//     
//     std::vector<double> answer(0.0);
//     
//     psd_pointer->set_critical_radius();
//     
//     psd_pointer->get_critical_radius(answer);
//     
//     double expectedAnswer = 2.507024;
//   
//     std::ostringstream streamOut;
//     streamOut <<"The value of the rc_HI (microns) is: "<<answer[0]<<". The expected value is: "<<expectedAnswer<<std::endl;
//     TEST_ASSERT_DELTA_MSG(expectedAnswer, answer[0], 1e-6, streamOut.str().c_str()); 
// }

//-------------------------------------------------------------
void PSD_HI_Test::testcompute_k_sat_HI()
{
//     double answer(0.0);
//     
//     psd_pointer->set_saturation();
//   
//     psd_pointer->get_global_saturated_permeability(answer); 
//     
//     double expectedAnswer = 58.1327578;
//   
//     std::ostringstream streamOut;
//     streamOut <<"The value of the k_sat_HI (microns^2) is: "<<answer<<". The expected value is: "<<expectedAnswer<<std::endl;
//     TEST_ASSERT_DELTA_MSG(expectedAnswer, answer, 1e-7, streamOut.str().c_str());
}

//-------------------------------------------------------------
void PSD_HI_Test::testcompute_sat_HI()
{
//     std::vector<double> answer(0.0);
//   
//     psd_pointer->get_saturation(answer); 
//     
//     std::cout<< "The Saturation is "<<answer[0]<<std::endl;
//     
//     double expectedAnswer = 0.0081234135;
//   
//     std::ostringstream streamOut;
//     streamOut <<"The value of the Saturation_HI is: "<<answer[0]<<std::cos(100*3.14/180)<<". The expected value is: "<<expectedAnswer<<std::endl;
//     TEST_ASSERT_DELTA_MSG(expectedAnswer, answer[0], 1e-7, streamOut.str().c_str());
    
  /*  
    for (int p_c = 100; p_c<=300000;p_c+=100)
    {
        double critical_radius_HO;
        if (p_c < 0)
        {
            critical_radius_HO = std::numeric_limits<double>::infinity();
        }
        else
        {
            critical_radius_HO =  -1.0 * 2.0 *  0.0728
            *( std::cos(100 * 3.14/180) /p_c); 
        }
        
        double SHO = 0.0;
        
        std::vector<double> rHO_k = {500e-9, 1200e-9, 2000e-9};
        std::vector<double> fHO_k = {0.18, 0.1, 0.72};
        std::vector<double> sHO_k = {0.65, 0.45, 0.75};

        for (unsigned int i= 0; i<fHO_k.size();i++)
        {
                SHO +=   0.5
                        *0.5 
                        * fHO_k[i] 
                        * ( 1.0 - std::erf ((std::log(critical_radius_HO) - std::log(rHO_k[i])) / (sHO_k[i] * std::pow(2.0,0.5))));    
        }
        
        
        
        
        
        
            
        double critical_radius_HI;
        
        if (p_c >= 0)
        {
            critical_radius_HI = std::numeric_limits<double>::infinity();
        }
        else 
        {
            critical_radius_HI =  -1.0 * 2.0 * 0.0728 
            * (   std::cos(80 * 3.14/180) /p_c );  
        }
        
        
        
        
        std::vector<double> rHI_k = {38.2e-9, 125e-9, 500e-9};
        std::vector<double> fHI_k = {0.62, 0.36, 0.02};
        std::vector<double> sHI_k = {0.43, 0.55, 0.65};
        
        double SHI=0.0;
          
            for(unsigned int i = 0; i < rHI_k.size(); ++i) 
            {
                    SHI +=   0.5 
                            * 0.5 
                            * fHI_k[i] * ( 1.0 + std::erf (   (std::log(critical_radius_HI) - std::log(rHI_k[i]) )/(  sHI_k[i] * std::pow(2.0,0.5)  )));      
            }       
        
        
        double S = 0.0;
        S = SHI + SHO;
        
//            std::cout<<p_c<<"           " <<S<<std::endl;
        
    
           
           
        double saturated_HI_permeability = 0.0;

        
        for(unsigned int i = 0; i < sHI_k.size(); ++i)
            {
                saturated_HI_permeability     +=  0.5
                
                                                  * (1.0/16.0)
                
                                              * std::pow( 0.4 * S / 1.0 ,2.0) 
                                              
                                              * std::exp ( (-2.0) * sHI_k[i] * sHI_k[i] ) 
                                              
                                              * rHI_k[i] * rHI_k[i] * fHI_k[i]
                                              
                                              * ( std::erf ( (std::log(critical_radius_HI) - std::log(rHI_k[i])) / (sHI_k[i] * std::pow(2.0,0.5)) - sHI_k[i] * std::pow(2.0,0.5) ) + 1.0);                                    
            }
           
           
           
           double saturated_HO_permeability = 0.0;
           
           
          for (unsigned int k = 0; k < rHO_k.size(); ++k)
            {
            saturated_HO_permeability     += std::pow(0.4 * S /1.0,2.0) 
            
                                             * (std::exp ((-2.0) * sHO_k[k] * sHO_k[k]) * rHO_k[k] * rHO_k[k] * fHO_k[k])
                                             
                                             * (1.0/16.0)
                                             
                                             * 0.5
                                             
                                             * ( (-1.0) * std::erf ( (std::log(critical_radius_HO) - std::log(rHO_k[k]))/(sHO_k[k] * std::sqrt(2.0)) - sHO_k[k] * std::sqrt(2.0) ) + 1.0 );
            }
           
           
           
           
           double get_global_saturated_permeability = 0.0;
           
           psd_pointer->get_global_saturated_permeability(get_global_saturated_permeability);
           
           
           double liquid_permeability=0.0;
           
           
           liquid_permeability= (saturated_HI_permeability + saturated_HO_permeability)/get_global_saturated_permeability;
           
           
           
           std::cout<<p_c<<"                    "<<liquid_permeability<<std::endl;
           
           
           
           
           
           
           
           
           */
           
           
           
           
/*           
           
        
        double HO_liquid_gas_interfacial_surface_a = 0.0;
        
        
            for (unsigned int k = 0; k < rHO_k.size(); ++k)
            {
            
                HO_liquid_gas_interfacial_surface_a   += 0.5
                                                         * fHO_k[k] 
                                                         * std::exp(sHO_k[k] * sHO_k[k] /2.0 )  
                                                         /  rHO_k[k]  
                                                         * ( 1.0 - std::erf ( (std::log(critical_radius_HO) - std::log(rHO_k[k])) / (sHO_k[k] * std::sqrt(2.0)) + sHO_k[k] * std::sqrt(2.0) / 2.0)  )
                                                         / 8.0;
            }                                               
        double HI_liquid_gas_interfacial_surface_a = 0.0;
                                                            
            
            
        for(unsigned int i = 0; i < sHI_k.size(); ++i)
        {
            HI_liquid_gas_interfacial_surface_a +=        0.5
                                                        * fHI_k[i] 
                                                        * std::exp(sHI_k[i] * sHI_k[i] /2.0 )  
                                                        * (1.0/ rHI_k[i])  
                                                        * ( 1.0 + std::erf ( (std::log(critical_radius_HI) - std::log(rHI_k[i])) / (sHI_k[i] * std::sqrt(2.0)) + sHI_k[i] * std::sqrt(2.0) / 2.0) )
                                                        / 8.0;
        }
                                                        
        double a_max = 0.0;
        
         psd_pointer->get_maximum_cross_sectional_areas(a_max);
        
        double a_max_HI = 0;
        double a_max_HO = 0;
        
        for(unsigned int i = 0; i < sHI_k.size(); ++i)
        
        {a_max_HI +=  0.5 * fHI_k[i] * std::exp (sHI_k[i] * sHI_k[i] /2.0 )/(4.0*rHI_k[i]);}
        
        for(unsigned int i = 0; i < sHO_k.size(); ++i)
        
        {a_max_HO +=  0.5 * fHO_k[i] * std::exp (sHO_k[i] * sHO_k[i] /2.0 )/(4.0*rHO_k[i])   ;}*/
//         
//         double a_max = 0.0;
//         a_max = a_max_HI + a_max_HO;
        
/*        
        double aaaaaaaaaaaa = 0;
        aaaaaaaaaaaa = HI_liquid_gas_interfacial_surface_a + HO_liquid_gas_interfacial_surface_a;
        
         double       liquid_gas_interfacial_surface = 0.0;
         double b = 0.0;
         b = aaaaaaaaaaaa/a_max;
         liquid_gas_interfacial_surface     = (b)
                                            *(1- b)
                                            *  aaaaaaaaaaaa;*/
                                            
                                            
                                            
                                            
                                            
 //         std::cout<<S<<" " <<liquid_gas_interfacial_surface<<std::endl;
        
        
        
//         double saturated_HI_permeability = 0.0;
//         
//             
//             saturated_HI_permeability  +=      0.28
//                                               *std::pow(0.6 * S  / 2.26 ,2.0) 
//                                               *std::exp ( (-2) * 1.0 * 1.0 ) 
//                                               * 14.2e-6 * 14.2e-6 * 1.0
//                                               *( std::erf ( (std::log(critical_radius_HI) - std::log(14.2e-6 )) / ( 1.0 * std::pow(2,0.5)) - 1.0 * std::pow(2,0.5) ) + 1)  
//                                               *1/16 ;
//         
//         
//                                               
//         double saturated_HO_permeability = 0.0;
// 
//             
//             saturated_HO_permeability  += std::pow(0.6 * S  /  2.26 ,2.0) 
//                                              * (std::exp ((-2) * 0.35 * 0.35) * 34.0e-6 * 34.0e-6 * 1.0)
//                                              * 1/16
//                                              * 0.72
//                                              * ( -std::erf ( (std::log(critical_radius_HO) - std::log(34.0e-6))/(0.35 * std::pow(2,0.5)) - 0.35 * std::pow(2,0.5) ) + 1 );
//         
//         
//         
//             double answer(0.0);
//             psd_pointer->set_saturation();
//             psd_pointer->get_global_saturated_permeability(answer); 
//             
//             
//         double k = (saturated_HI_permeability + saturated_HO_permeability) / answer;
        
        
        
//         std::cout << "haha   "<<S<<" "<<k<<std::endl;
        
        
/*        
        
        
        
    }
    */
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
}

//-------------------------------------------------------------
// void PSD_HI_Test::testcompute_k_L_HI()
// {
//     std::vector<double> answer(0.0);
//     psd_object->get_pore_HI_liquid_saturated_permeability(answer);  
//     
//     double expectedAnswer = 2.9317742175e-9;
//   
//     std::ostringstream streamOut;
//     streamOut <<"The value of the K_L_HI(microns^2) is: "<<answer[0]<<". The expected value is: "<<expectedAnswer<<std::endl;
//     TEST_ASSERT_DELTA_MSG(expectedAnswer, answer[0], 1e-13, streamOut.str().c_str());
// }

//-------------------------------------------------------------
void PSD_HI_Test::testcompute_kr_L_HI()
{
//     std::vector<double> answer(0.0);
//     
//     psd_pointer->get_relative_liquid_permeability(answer);  
//     
//     double expectedAnswer = 5.04323951277e-11;
//   
//     std::ostringstream streamOut;
//     streamOut <<"The value of the Kr_L_HI is: "<<answer[0]<<". The expected value is: "<<expectedAnswer<<std::endl;
//     TEST_ASSERT_DELTA_MSG(expectedAnswer, answer[0], 1e-15, streamOut.str().c_str());
}

//-------------------------------------------------------------
// void PSD_HI_Test::testcompute_k_G_HI()
// {
//     std::vector<double> answer(0.0);
//   
//     psd_object.get_pore_HI_gas_saturated_permeability(answer);  
//     
//     double expectedAnswer = 40.0344411166;
//   
//     std::ostringstream streamOut;
//     streamOut <<"The value of the K_G_HI(microns^2) is: "<<answer[0]<<". The expected value is: "<<expectedAnswer<<std::endl;
//     TEST_ASSERT_DELTA_MSG(expectedAnswer, answer[0], 1e-5, streamOut.str().c_str());
// }

//-------------------------------------------------------------
void PSD_HI_Test::testcompute_kr_G_HI()
{
//     std::vector<double> answer(0.0);
//   
//     psd_pointer->get_relative_gas_permeability(answer);  
//     
//     double expectedAnswer = 0.6886726621;
//   
//     std::ostringstream streamOut;
//     streamOut <<"The value of the Kr_G_HI is: "<<answer[0]<<". The expected value is: "<<expectedAnswer<<std::endl;
//     TEST_ASSERT_DELTA_MSG(expectedAnswer, answer[0], 1e-7, streamOut.str().c_str());
}

//-------------------------------------------------------------
void PSD_HI_Test::testcompute_interfacial_area_per_volume_HI()
{
//     std::vector<double> answer(0.0);
//   
//     psd_pointer->get_liquid_gas_interfacial_surface(answer);  
//     double expectedAnswer = 0.000113962;
//   
//     std::ostringstream streamOut;
//     streamOut <<"The value of the interfacial area per unit volume for hydrophilic pores is: "<<answer[0]<<". The expected value is: "<<expectedAnswer<<std::endl;
//     TEST_ASSERT_DELTA_MSG(expectedAnswer, answer[0], 1e-9, streamOut.str().c_str());
}

//-------------------------------------------------------------
//-------------------------------------------------------------