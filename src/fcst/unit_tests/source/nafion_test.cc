//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: nafion_test.cc
//    - Description: Unit testing class for Nafion
//    - Developers: Madhur Bhaiya
//    - Id: $Id: nafion_test.cc 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#include "nafion_test.h"

// ------------------------------------------------------------
void
NafionTest::testSetT()
{
    electrolyte = new FuelCellShop::Material::Nafion;
    
    TEST_ASSERT_DELTA_MSG(electrolyte->T, 0., 1.e-6, "correct initialization !!");
    TEST_ASSERT_DELTA_MSG(electrolyte->p_total, 0., 1.e-6, "correct initialization !!");
    TEST_ASSERT_DELTA_MSG(electrolyte->lambda, 12., 1.e-6, "correct initialization !!");
    
    electrolyte->set_T(353.0);
    TEST_ASSERT_DELTA_MSG(electrolyte->T, 353., 1.e-6, "set_T method success !");
    TEST_ASSERT_DELTA_MSG(electrolyte->p_total, 0., 1.e-6, "set_T method success !");
    TEST_ASSERT_DELTA_MSG(electrolyte->lambda, 12., 1.e-6, "set_T method success !");

    delete electrolyte;
}

// ------------------------------------------------------------
void
NafionTest::testSetPT()
{
    electrolyte = new FuelCellShop::Material::Nafion;
    
    TEST_ASSERT_DELTA_MSG(electrolyte->T, 0., 1.e-6, "correct initialization !!");
    TEST_ASSERT_DELTA_MSG(electrolyte->p_total, 0., 1.e-6, "correct initialization !!");
    TEST_ASSERT_DELTA_MSG(electrolyte->lambda, 12., 1.e-6, "correct initialization !!");

    electrolyte->set_p_t(101325.);
    TEST_ASSERT_DELTA_MSG(electrolyte->T, 0., 1.e-6, "set_p_t method success !");
    TEST_ASSERT_DELTA_MSG(electrolyte->p_total, 101325., 1.e-6, "set_p_t method success !");
    TEST_ASSERT_DELTA_MSG(electrolyte->lambda, 12., 1.e-6, "set_p_t method success !");
    
    delete electrolyte;

}

// ------------------------------------------------------------
void
NafionTest::testSetLambda()
{
    electrolyte = new FuelCellShop::Material::Nafion;
    
    TEST_ASSERT_DELTA_MSG(electrolyte->T, 0., 1.e-6, "correct initialization !!");
    TEST_ASSERT_DELTA_MSG(electrolyte->p_total, 0., 1.e-6, "correct initialization !!");
    TEST_ASSERT_DELTA_MSG(electrolyte->lambda, 12., 1.e-6, "correct initialization !!");
    
    electrolyte->set_lambda(14.);
    TEST_ASSERT_DELTA_MSG(electrolyte->T, 0., 1.e-6, "set_lambda method success !");
    TEST_ASSERT_DELTA_MSG(electrolyte->p_total, 0., 1.e-6, "set_lambda method success !");
    TEST_ASSERT_DELTA_MSG(electrolyte->lambda, 14., 1.e-6, "set_lambda method success !");
    
    delete electrolyte;
}

// ------------------------------------------------------------
void
NafionTest::testSetTemperature()
{
    electrolyte = new FuelCellShop::Material::Nafion;

    electrolyte->set_temperature( FuelCellShop::SolutionVariable(353., 1, temperature_of_REV) );
    
    TEST_ASSERT_DELTA_MSG(electrolyte->T_var[0], 353., 1.e-6, "set_temperature method success !");
    TEST_ASSERT_DELTA_MSG(electrolyte->lambda_var[0], 12., 1.e-6, "set_lambda before set_temperature method success !");
    
    delete electrolyte;
}

// ------------------------------------------------------------
void
NafionTest::testSetMembraneWaterContent()
{
    electrolyte = new FuelCellShop::Material::Nafion;
    
    electrolyte->set_T(353.);
    electrolyte->set_membrane_water_content( FuelCellShop::SolutionVariable(14., 1, membrane_water_content) );
    
    TEST_ASSERT_DELTA_MSG(electrolyte->T_var[0], 353., 1.e-6, "set_T before set_membrane_water_content method success !");
    TEST_ASSERT_DELTA_MSG(electrolyte->lambda_var[0], 14., 1.e-6, "set_membrane_water_content method success !");
    
    delete electrolyte;
}

// ------------------------------------------------------------
void
NafionTest::testSetWaterMolarFraction()
{
    electrolyte = new FuelCellShop::Material::Nafion;
    
    electrolyte->set_T(353.);
    electrolyte->set_water_molar_fraction( FuelCellShop::SolutionVariable(0.5, 1, water_molar_fraction) );
    
    TEST_ASSERT_DELTA_MSG(electrolyte->T_var[0], 353., 1.e-6, "set_T before set_water_molar_fraction method success !");
    TEST_ASSERT_DELTA_MSG(electrolyte->xwater_var[0], 0.5, 1.e-6, "set_water_molar_fraction method success !");
    
    delete electrolyte;
}

// ------------------------------------------------------------
void
NafionTest::testProtonConductivityConstant()
{
    electrolyte = new FuelCellShop::Material::Nafion;
    electrolyte->modify_parameters(1, "Constant");
        
    // Double return method
    double sigma_p;
    electrolyte->set_T(353.);
    electrolyte->proton_conductivity(sigma_p);
    TEST_ASSERT_DELTA_MSG(sigma_p, 0.1, 1.e-6, "proton_conductivity Constant double return method success !");

    // Vector return method
    std::vector<double> sigma_vec;
    electrolyte->set_membrane_water_content( FuelCellShop::SolutionVariable(12., 1, membrane_water_content) );
    electrolyte->proton_conductivity(sigma_vec);
    TEST_ASSERT_DELTA_MSG(sigma_vec[0], 0.1, 1.e-6, "proton_conductivity Constant vector return method success !");

    // Derivatives
    std::vector<VariableNames> deriv_flags;
    deriv_flags.push_back(temperature_of_REV);
    deriv_flags.push_back(membrane_water_content);
    electrolyte->set_derivative_flags(deriv_flags);
    std::map< VariableNames, std::vector<double> > dsigma;
    electrolyte->proton_conductivity_derivative(dsigma);
    
    double delta(std::pow(10.,-8));
    std::vector<double> sigma_delta_vec;
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(353.+delta, 1, temperature_of_REV) );
    electrolyte->proton_conductivity(sigma_delta_vec);
    TEST_ASSERT_DELTA_MSG((sigma_delta_vec[0]-sigma_vec[0])/delta, dsigma[temperature_of_REV][0], 1.e-6, "proton_conductivity Constant temperature derivative method success !");
    
    electrolyte->set_membrane_water_content( FuelCellShop::SolutionVariable(12.+delta, 1, membrane_water_content) );
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(353., 1, temperature_of_REV) );
    electrolyte->proton_conductivity(sigma_delta_vec);
    TEST_ASSERT_DELTA_MSG((sigma_delta_vec[0]-sigma_vec[0])/delta, dsigma[membrane_water_content][0], 1.e-6, "proton_conductivity Constant lambda derivative method success !");
    
    delete electrolyte;
}

// ------------------------------------------------------------
void
NafionTest::testProtonConductivitySpringer()
{
    electrolyte = new FuelCellShop::Material::Nafion;
    electrolyte->modify_parameters(1, "Springer");
        
    // Double return method
    double sigma_p;
    electrolyte->set_T(353.);
    electrolyte->proton_conductivity(sigma_p);
    TEST_ASSERT_DELTA_MSG(sigma_p, 0.1056575665, 1.e-6, "proton_conductivity Springer double return method success !");

    // Vector return method
    std::vector<double> sigma_vec;
    electrolyte->set_membrane_water_content( FuelCellShop::SolutionVariable(12., 1, membrane_water_content) );
    electrolyte->proton_conductivity(sigma_vec);
    TEST_ASSERT_DELTA_MSG(sigma_vec[0], 0.1056575665, 1.e-6, "proton_conductivity Springer vector return method success !");

    // Derivatives
    std::vector<VariableNames> deriv_flags;
    deriv_flags.push_back(temperature_of_REV);
    deriv_flags.push_back(membrane_water_content);
    electrolyte->set_derivative_flags(deriv_flags);
    std::map< VariableNames, std::vector<double> > dsigma;
    electrolyte->proton_conductivity_derivative(dsigma);
    
    double delta(std::pow(10.,-8));
    std::vector<double> sigma_delta_vec;
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(353.+delta, 1, temperature_of_REV) );
    electrolyte->proton_conductivity(sigma_delta_vec);
    TEST_ASSERT_DELTA_MSG((sigma_delta_vec[0]-sigma_vec[0])/delta, dsigma[temperature_of_REV][0], 1.e-6, "proton_conductivity Springer temperature derivative method success !");
    
    electrolyte->set_membrane_water_content( FuelCellShop::SolutionVariable(12.+delta, 1, membrane_water_content) );
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(353., 1, temperature_of_REV) );
    electrolyte->proton_conductivity(sigma_delta_vec);
    TEST_ASSERT_DELTA_MSG((sigma_delta_vec[0]-sigma_vec[0])/delta, dsigma[membrane_water_content][0], 1.e-6, "proton_conductivity Springer lambda derivative method success !");
    
    delete electrolyte;
}

// ------------------------------------------------------------
void
NafionTest::testProtonConductivityNRE211()
{
    electrolyte = new FuelCellShop::Material::Nafion;
    electrolyte->modify_parameters(1, "NRE211");
        
    // Double return method
    double sigma_p;
    electrolyte->set_T(353.);
    electrolyte->proton_conductivity(sigma_p);
    TEST_ASSERT_DELTA_MSG(sigma_p, 0.1293429334, 1.e-6, "proton_conductivity NRE211 double return method success !");

    // Vector return method
    std::vector<double> sigma_vec;
    electrolyte->set_membrane_water_content( FuelCellShop::SolutionVariable(12., 1, membrane_water_content) );
    electrolyte->proton_conductivity(sigma_vec);
    TEST_ASSERT_DELTA_MSG(sigma_vec[0], 0.1293429334, 1.e-6, "proton_conductivity NRE211 vector return method success !");

    // Derivatives
    std::vector<VariableNames> deriv_flags;
    deriv_flags.push_back(temperature_of_REV);
    deriv_flags.push_back(membrane_water_content);
    electrolyte->set_derivative_flags(deriv_flags);
    std::map< VariableNames, std::vector<double> > dsigma;
    electrolyte->proton_conductivity_derivative(dsigma);
    
    double delta(std::pow(10.,-8));
    std::vector<double> sigma_delta_vec;
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(353.+delta, 1, temperature_of_REV) );
    electrolyte->proton_conductivity(sigma_delta_vec);
    TEST_ASSERT_DELTA_MSG((sigma_delta_vec[0]-sigma_vec[0])/delta, dsigma[temperature_of_REV][0], 1.e-6, "proton_conductivity NRE211 temperature derivative method success !");
    
    electrolyte->set_membrane_water_content( FuelCellShop::SolutionVariable(12.+delta, 1, membrane_water_content) );
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(353., 1, temperature_of_REV) );
    electrolyte->proton_conductivity(sigma_delta_vec);
    TEST_ASSERT_DELTA_MSG((sigma_delta_vec[0]-sigma_vec[0])/delta, dsigma[membrane_water_content][0], 1.e-6, "proton_conductivity NRE211 lambda derivative method success !");
    
    delete electrolyte;
}

// ------------------------------------------------------------
void
NafionTest::testProtonConductivityIden11()
{
    electrolyte = new FuelCellShop::Material::Nafion;
    electrolyte->modify_parameters(1, "Iden11");
    
    // Double return method
    double sigma_p;
    electrolyte->set_T(353.);
    electrolyte->proton_conductivity(sigma_p);
    TEST_ASSERT_DELTA_MSG(sigma_p, 0.1869276873, 1.e-6, "proton_conductivity Iden11 method1 success !");
    
    electrolyte->set_lambda(14.);
    electrolyte->set_T(368.);
    electrolyte->proton_conductivity(sigma_p);
    TEST_ASSERT_DELTA_MSG(sigma_p, 0.2097970581, 1.e-6, "proton_conductivity Iden11 method2 success !");
    
    // Vector return method
    std::vector<double> sigma_vec;   
    std::vector<double> lambda_vec;
    std::vector<double> T_vec;
    lambda_vec.push_back(12.); lambda_vec.push_back(14.);
    T_vec.push_back(353.); T_vec.push_back(368.);
    electrolyte->set_membrane_water_content( FuelCellShop::SolutionVariable(&lambda_vec, membrane_water_content) );
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(&T_vec, temperature_of_REV) );
    electrolyte->proton_conductivity(sigma_vec);
    TEST_ASSERT_DELTA_MSG(sigma_vec[0], 0.1869276873, 1.e-6, "proton_conductivity Iden11 method1 success !");
    TEST_ASSERT_DELTA_MSG(sigma_vec[1], 0.2097970581, 1.e-6, "proton_conductivity Iden11 method2 success !");
    
    // Derivatives
    std::vector<VariableNames> deriv_flags;
    deriv_flags.push_back(temperature_of_REV);
    deriv_flags.push_back(membrane_water_content);
    electrolyte->set_derivative_flags(deriv_flags);
    std::map< VariableNames, std::vector<double> > dsigma;
    electrolyte->proton_conductivity_derivative(dsigma);
    
    double delta(std::pow(10.,-8));
    std::vector<double> sigma_delta_vec;
    std::vector<double> T_delta_vec;
    T_delta_vec.push_back(353.+delta); T_delta_vec.push_back(368.+delta);
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(&T_delta_vec, temperature_of_REV) );
    electrolyte->proton_conductivity(sigma_delta_vec);
    TEST_ASSERT_DELTA_MSG((sigma_delta_vec[0]-sigma_vec[0])/delta, dsigma[temperature_of_REV][0], 1.e-6, "proton_conductivity Iden11 temperature derivative method1 success !");
    TEST_ASSERT_DELTA_MSG((sigma_delta_vec[1]-sigma_vec[1])/delta, dsigma[temperature_of_REV][1], 1.e-6, "proton_conductivity Iden11 temperature derivative method2 success !");
    
    std::vector<double> lambda_delta_vec;
    lambda_delta_vec.push_back(12.+delta); lambda_delta_vec.push_back(14.+delta);
    electrolyte->set_membrane_water_content( FuelCellShop::SolutionVariable(&lambda_delta_vec, membrane_water_content) );
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(&T_vec, temperature_of_REV) );
    electrolyte->proton_conductivity(sigma_delta_vec);
    TEST_ASSERT_DELTA_MSG((sigma_delta_vec[0]-sigma_vec[0])/delta, dsigma[membrane_water_content][0], 1.e-6, "proton_conductivity Iden11 lambda derivative method1 success !");
    TEST_ASSERT_DELTA_MSG((sigma_delta_vec[1]-sigma_vec[1])/delta, dsigma[membrane_water_content][1], 1.e-6, "proton_conductivity Iden11 lambda derivative method2 success !");
    
    delete electrolyte;
}

// ------------------------------------------------------------
void
NafionTest::testWaterDiffusivityConstant()
{
    electrolyte = new FuelCellShop::Material::Nafion;
    electrolyte->modify_parameters(2, "Constant");
        
    // Double return method
    double Dw;
    electrolyte->set_T(360.);
    electrolyte->water_diffusivity(Dw);
    TEST_ASSERT_DELTA_MSG(Dw, 2.e-3, 1.e-6, "water_diffusivity Constant double return method success !");

    // Vector return method
    std::vector<double> Dw_vec;
    electrolyte->set_membrane_water_content( FuelCellShop::SolutionVariable(12., 1, membrane_water_content) );
    electrolyte->water_diffusivity(Dw_vec);
    TEST_ASSERT_DELTA_MSG(Dw_vec[0], 2.e-3, 1.e-6, "water_diffusivity Springer vector return method success !");

    // Derivatives
    std::vector<VariableNames> deriv_flags;
    deriv_flags.push_back(temperature_of_REV);
    deriv_flags.push_back(membrane_water_content);
    electrolyte->set_derivative_flags(deriv_flags);
    std::map< VariableNames, std::vector<double> > dDw;
    electrolyte->water_diffusivity_derivative(dDw);
    
    double delta(std::pow(10.,-8));
    std::vector<double> Dw_delta_vec;
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(360.+delta, 1, temperature_of_REV) );
    electrolyte->water_diffusivity(Dw_delta_vec);
    TEST_ASSERT_DELTA_MSG((Dw_delta_vec[0]-Dw_vec[0])/delta, dDw[temperature_of_REV][0], 1.e-6, "water_diffusivity Constant temperature derivative method success !");
    
    electrolyte->set_membrane_water_content( FuelCellShop::SolutionVariable(12.+delta, 1, membrane_water_content) );
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(360., 1, temperature_of_REV) );
    electrolyte->water_diffusivity(Dw_delta_vec);
    TEST_ASSERT_DELTA_MSG((Dw_delta_vec[0]-Dw_vec[0])/delta, dDw[membrane_water_content][0], 1.e-6, "water_diffusivity Constant lambda derivative method success !");
    
    delete electrolyte;
}

// ------------------------------------------------------------
void
NafionTest::testWaterDiffusivitySpringer()
{
    electrolyte = new FuelCellShop::Material::Nafion;
    electrolyte->modify_parameters(2, "Springer");
        
    // Double return method
    double Dw;
    electrolyte->set_T(360.);
    electrolyte->water_diffusivity(Dw);
    TEST_ASSERT_DELTA_MSG(Dw, 0.0000044004717, 1.e-6, "water_diffusivity Springer double return method success !");

    // Vector return method
    std::vector<double> Dw_vec;
    electrolyte->set_membrane_water_content( FuelCellShop::SolutionVariable(12., 1, membrane_water_content) );
    electrolyte->water_diffusivity(Dw_vec);
    TEST_ASSERT_DELTA_MSG(Dw_vec[0], 0.0000044004717, 1.e-6, "proton_conductivity Springer vector return method success !");

    // Derivatives
    std::vector<VariableNames> deriv_flags;
    deriv_flags.push_back(temperature_of_REV);
    deriv_flags.push_back(membrane_water_content);
    electrolyte->set_derivative_flags(deriv_flags);
    std::map< VariableNames, std::vector<double> > dDw;
    electrolyte->water_diffusivity_derivative(dDw);
    
    double delta(std::pow(10.,-8));
    std::vector<double> Dw_delta_vec;
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(360.+delta, 1, temperature_of_REV) );
    electrolyte->water_diffusivity(Dw_delta_vec);
    TEST_ASSERT_DELTA_MSG((Dw_delta_vec[0]-Dw_vec[0])/delta, dDw[temperature_of_REV][0], 1.e-6, "water_diffusivity Springer temperature derivative method success !");
    
    electrolyte->set_membrane_water_content( FuelCellShop::SolutionVariable(12.+delta, 1, membrane_water_content) );
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(360., 1, temperature_of_REV) );
    electrolyte->water_diffusivity(Dw_delta_vec);
    TEST_ASSERT_DELTA_MSG((Dw_delta_vec[0]-Dw_vec[0])/delta, dDw[membrane_water_content][0], 1.e-6, "water_diffusivity Springer lambda derivative method success !");
    
    delete electrolyte;
}

// ------------------------------------------------------------
void
NafionTest::testWaterDiffusivityMotupally()
{
    electrolyte = new FuelCellShop::Material::Nafion;
    electrolyte->modify_parameters(2, "Motupally");
        
    // Double return method
    double Dw;
    electrolyte->set_T(360.);
    electrolyte->water_diffusivity(Dw);
    TEST_ASSERT_DELTA_MSG(Dw, 0.0000048722066, 1.e-6, "water_diffusivity Motupally double return method1 success !");
    electrolyte->set_lambda(2.);
    electrolyte->set_T(355.);
    TEST_ASSERT_DELTA_MSG(Dw, 0.0000048722066, 1.e-6, "water_diffusivity Motupally double return method2 success !");
    
    // Vector return method
    std::vector<double> Dw_vec;
    std::vector<double> lambda_vec, T_vec;
    lambda_vec.push_back(12.); lambda_vec.push_back(2.);
    T_vec.push_back(360.), T_vec.push_back(355.);
    electrolyte->set_membrane_water_content( FuelCellShop::SolutionVariable(&lambda_vec, membrane_water_content) );
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(&T_vec, temperature_of_REV) );
    electrolyte->water_diffusivity(Dw_vec);
    TEST_ASSERT_DELTA_MSG(Dw_vec[0], 0.0000048722066, 1.e-6, "water_diffusivity Motupally vector return method1 success !");
    TEST_ASSERT_DELTA_MSG(Dw_vec[1], 0.0000048722066, 1.e-6, "water_diffusivity Motupally vector return method2 success !");

    // Derivatives
    std::vector<VariableNames> deriv_flags;
    deriv_flags.push_back(temperature_of_REV);
    deriv_flags.push_back(membrane_water_content);
    electrolyte->set_derivative_flags(deriv_flags);
    std::map< VariableNames, std::vector<double> > dDw;
    electrolyte->water_diffusivity_derivative(dDw);

    
    double delta(std::pow(10.,-8));
    std::vector<double> Dw_delta_vec;
    std::vector<double> lambda_delta_vec, T_delta_vec;
    lambda_delta_vec.push_back(12.+delta); lambda_delta_vec.push_back(2.+delta);
    T_delta_vec.push_back(360.+delta), T_delta_vec.push_back(355.+delta);
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(&T_delta_vec, temperature_of_REV) );
    electrolyte->water_diffusivity(Dw_delta_vec);    
    TEST_ASSERT_DELTA_MSG((Dw_delta_vec[0]-Dw_vec[0])/delta, dDw[temperature_of_REV][0], 1.e-6, "water_diffusivity Motupally temperature derivative method1 success !");
    TEST_ASSERT_DELTA_MSG((Dw_delta_vec[1]-Dw_vec[1])/delta, dDw[temperature_of_REV][1], 1.e-6, "water_diffusivity Motupally temperature derivative method2 success !");
    
    electrolyte->set_membrane_water_content( FuelCellShop::SolutionVariable(&lambda_delta_vec, membrane_water_content) );
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(&T_vec, temperature_of_REV) );
    electrolyte->water_diffusivity(Dw_delta_vec);    
    TEST_ASSERT_DELTA_MSG((Dw_delta_vec[0]-Dw_vec[0])/delta, dDw[membrane_water_content][0], 1.e-6, "water_diffusivity Motupally lambda derivative method1 success !");
    TEST_ASSERT_DELTA_MSG((Dw_delta_vec[1]-Dw_vec[1])/delta, dDw[membrane_water_content][1], 1.e-6, "water_diffusivity Motupally lambda derivative method2 success !");
    
    delete electrolyte;
}

// ------------------------------------------------------------
void
NafionTest::testElectroOsmoticDragConstant()
{
    electrolyte = new FuelCellShop::Material::Nafion;
    electrolyte->modify_parameters(3, "Constant");

    // Vector return method
    std::vector<double> nd_vec;
    electrolyte->set_membrane_water_content( FuelCellShop::SolutionVariable(12., 1, membrane_water_content) );
    electrolyte->electroosmotic_drag(nd_vec);
    TEST_ASSERT_DELTA_MSG(nd_vec[0], 1., 1.e-6, "electroosmotic_drag Constant vector return method success !");

    // Derivatives
    std::vector<VariableNames> deriv_flags;
    deriv_flags.push_back(temperature_of_REV);
    deriv_flags.push_back(membrane_water_content);
    electrolyte->set_derivative_flags(deriv_flags);
    std::map< VariableNames, std::vector<double> > dnd;
    electrolyte->electroosmotic_drag_derivative(dnd);
    TEST_ASSERT_DELTA_MSG(0., dnd[temperature_of_REV][0], 1.e-6, "electroosmotic_drag Constant temperature derivative method success !");
    TEST_ASSERT_DELTA_MSG(0., dnd[membrane_water_content][0], 1.e-6, "electroosmotic_drag Constant lambda derivative method success !");
    
    delete electrolyte;
}

// ------------------------------------------------------------
void
NafionTest::testElectroOsmoticDragSpringer()
{
    electrolyte = new FuelCellShop::Material::Nafion;
    electrolyte->modify_parameters(3, "Springer");

    // Vector return method
    std::vector<double> nd_vec;
    electrolyte->set_membrane_water_content( FuelCellShop::SolutionVariable(12., 1, membrane_water_content) );
    electrolyte->electroosmotic_drag(nd_vec);
    TEST_ASSERT_DELTA_MSG(nd_vec[0], 1.363636364, 1.e-6, "electroosmotic_drag Springer vector return method success !");

    // Derivatives
    std::vector<VariableNames> deriv_flags;
    deriv_flags.push_back(temperature_of_REV);
    deriv_flags.push_back(membrane_water_content);
    electrolyte->set_derivative_flags(deriv_flags);
    std::map< VariableNames, std::vector<double> > dnd;
    electrolyte->electroosmotic_drag_derivative(dnd);
    TEST_ASSERT_DELTA_MSG(0., dnd[temperature_of_REV][0], 1.e-6, "electroosmotic_drag Constant temperature derivative method success !");
    TEST_ASSERT_DELTA_MSG(0.1136363636, dnd[membrane_water_content][0], 1.e-6, "electroosmotic_drag Constant lambda derivative method success !");
    
    delete electrolyte;
}

// ------------------------------------------------------------
void
NafionTest::testThermoOsmoticCoeffConstant()
{
    electrolyte = new FuelCellShop::Material::Nafion;
    electrolyte->modify_parameters(4, "Constant");

    // Vector return method
    std::vector<double> DT_vec;
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(360., 1, temperature_of_REV) );
    electrolyte->thermoosmotic_coeff(DT_vec);
    TEST_ASSERT_DELTA_MSG(DT_vec[0], -1.3e-7, 1.e-6, "thermoosmotic_coeff Constant vector return method success !");

    // Derivatives
    std::vector<VariableNames> deriv_flags;
    deriv_flags.push_back(temperature_of_REV);
    deriv_flags.push_back(membrane_water_content);
    electrolyte->set_derivative_flags(deriv_flags);
    std::map< VariableNames, std::vector<double> > dDT;
    electrolyte->thermoosmotic_coeff_derivative(dDT);
    TEST_ASSERT_DELTA_MSG(0., dDT[temperature_of_REV][0], 1.e-6, "thermoosmotic_coeff Constant temperature derivative method success !");
    TEST_ASSERT_DELTA_MSG(0., dDT[membrane_water_content][0], 1.e-6, "thermoosmotic_coeff Constant lambda derivative method success !");
    
    delete electrolyte;
}

// ------------------------------------------------------------
void
NafionTest::testThermoOsmoticCoeffKim()
{
    electrolyte = new FuelCellShop::Material::Nafion;
    electrolyte->modify_parameters(4, "Kim09");

    // Vector return method
    std::vector<double> DT_vec;
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(360., 1, temperature_of_REV) );
    electrolyte->thermoosmotic_coeff(DT_vec);
    TEST_ASSERT_DELTA_MSG(DT_vec[0], -1.470885913e-7, 1.e-6, "thermoosmotic_coeff Kim vector return method success !");

    // Derivatives
    std::vector<VariableNames> deriv_flags;
    deriv_flags.push_back(temperature_of_REV);
    deriv_flags.push_back(membrane_water_content);
    electrolyte->set_derivative_flags(deriv_flags);
    std::map< VariableNames, std::vector<double> > dDT;
    electrolyte->thermoosmotic_coeff_derivative(dDT);
    
    double delta(std::pow(10.,-8));
    std::vector<double> DT_delta_vec;
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(360.+delta, 1, temperature_of_REV) );
    electrolyte->thermoosmotic_coeff(DT_delta_vec);
    TEST_ASSERT_DELTA_MSG((DT_delta_vec[0]-DT_vec[0])/delta, dDT[temperature_of_REV][0], 1.e-6, "thermoosmotic_coeff Kim temperature derivative method success !");
    TEST_ASSERT_DELTA_MSG(0., dDT[membrane_water_content][0], 1.e-6, "thermoosmotic_coeff Kim lambda derivative method success !");
    
    delete electrolyte;
}

// ------------------------------------------------------------
void
NafionTest::testOxygenDiffusivity()
{
    electrolyte = new FuelCellShop::Material::Nafion;
    
    // Double return method
    double Doxy;
    electrolyte->set_T(355.);
    electrolyte->oxygen_diffusivity(Doxy);
    TEST_ASSERT_DELTA_MSG(Doxy, 0.0000104973561, 1.e-6, "oxygen_diffusivity double return method success !");

    // Vector return method
    std::vector<double> Doxy_vec;
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(355., 1, temperature_of_REV) );
    electrolyte->oxygen_diffusivity(Doxy_vec);
    TEST_ASSERT_DELTA_MSG(Doxy_vec[0], 0.0000104973561, 1.e-6, "oxygen_diffusivity vector return method success !");

    // Derivatives
    std::vector<VariableNames> deriv_flags;
    deriv_flags.push_back(temperature_of_REV);
    electrolyte->set_derivative_flags(deriv_flags);
    std::map< VariableNames, std::vector<double> > dDoxy;
    electrolyte->oxygen_diffusivity_derivative(dDoxy);
    
    double delta(std::pow(10.,-8));
    std::vector<double> Doxy_delta_vec;
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(355.+delta, 1, temperature_of_REV) );
    electrolyte->oxygen_diffusivity(Doxy_delta_vec);
    TEST_ASSERT_DELTA_MSG((Doxy_delta_vec[0]-Doxy_vec[0])/delta, dDoxy[temperature_of_REV][0], 1.e-6, "oxygen_diffusivity temperature derivative method success !");
    
    delete electrolyte;
}

// ------------------------------------------------------------
void
NafionTest::testSorptionEnthalpyConstant()
{
    electrolyte = new FuelCellShop::Material::Nafion;
    electrolyte->modify_parameters(6, "Constant");

    // Vector return method
    std::vector<double> hsorp_vec;
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(360., 1, temperature_of_REV) );
    electrolyte->sorption_enthalpy(hsorp_vec);
    TEST_ASSERT_DELTA_MSG(hsorp_vec[0], 45000., 1.e-6, "sorption enthalpy Constant vector return method success !");

    // Derivatives
    std::vector<VariableNames> deriv_flags;
    deriv_flags.push_back(temperature_of_REV);
    electrolyte->set_derivative_flags(deriv_flags);
    std::map< VariableNames, std::vector<double> > dhsorp;
    electrolyte->sorption_enthalpy_derivative(dhsorp);
    TEST_ASSERT_DELTA_MSG(0., dhsorp[temperature_of_REV][0], 1.e-6, "sorption enthalpy Constant temperature derivative method success !");
    
    delete electrolyte;
}

// ------------------------------------------------------------
void
NafionTest::testSorptionIsothermHinatsu()
{
    electrolyte = new FuelCellShop::Material::Nafion;
    electrolyte->set_p_t(101325.);
    electrolyte->modify_parameters(5, "Hinatsu");
    
    // Vector return method
    std::vector<double> lambda_eq_vec;
    std::vector<double> xwater_vec;   
    std::vector<double> T_vec;
    xwater_vec.push_back(0.4); xwater_vec.push_back(0.7);
    T_vec.push_back(353.); T_vec.push_back(360.);
    electrolyte->set_water_molar_fraction( FuelCellShop::SolutionVariable(&xwater_vec, water_molar_fraction) );
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(&T_vec, temperature_of_REV) );
    electrolyte->sorption_isotherm(lambda_eq_vec);
    TEST_ASSERT_DELTA_MSG(lambda_eq_vec[0], 6.750490538, 1.e-6, "sorption_isotherm Hinatsu method1 success !");
    TEST_ASSERT_DELTA_MSG(lambda_eq_vec[1], 9.2, 1.e-6, "sorption_isotherm Hinatsu method2 success !");
    
    // Derivatives
    std::vector<VariableNames> deriv_flags;
    deriv_flags.push_back(temperature_of_REV);
    deriv_flags.push_back(water_molar_fraction);
    electrolyte->set_derivative_flags(deriv_flags);
    std::map< VariableNames, std::vector<double> > dlambda_eq;
    electrolyte->sorption_isotherm_derivative(dlambda_eq);
    
    double delta(std::pow(10.,-8));
    std::vector<double> lambdaeq_delta_vec;
    std::vector<double> T_delta_vec;
    T_delta_vec.push_back(353.+delta); T_delta_vec.push_back(360.+delta);
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(&T_delta_vec, temperature_of_REV) );
    electrolyte->sorption_isotherm(lambdaeq_delta_vec);
    TEST_ASSERT_DELTA_MSG((lambdaeq_delta_vec[0]-lambda_eq_vec[0])/delta, dlambda_eq[temperature_of_REV][0], 1.e-6, "sorption_isotherm Hinatsu temperature derivative method1 success !");
    TEST_ASSERT_DELTA_MSG((lambdaeq_delta_vec[1]-lambda_eq_vec[1])/delta, dlambda_eq[temperature_of_REV][1], 1.e-6, "sorption_isotherm Hinatsu temperature derivative method2 success !");
    
    std::vector<double> xwater_delta_vec;
    xwater_delta_vec.push_back(0.4+delta); xwater_delta_vec.push_back(0.7+delta);
    electrolyte->set_water_molar_fraction( FuelCellShop::SolutionVariable(&xwater_delta_vec, water_molar_fraction) );
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(&T_vec, temperature_of_REV) );
    electrolyte->sorption_isotherm(lambdaeq_delta_vec);
    TEST_ASSERT_DELTA_MSG((lambdaeq_delta_vec[0]-lambda_eq_vec[0])/delta, dlambda_eq[water_molar_fraction][0], 1.e-5, "sorption_isotherm Hinatsu water molar fraction derivative method1 success !");
    TEST_ASSERT_DELTA_MSG((lambdaeq_delta_vec[1]-lambda_eq_vec[1])/delta, dlambda_eq[water_molar_fraction][1], 1.e-6, "sorption_isotherm Hinatsu water molar fraction derivative method2 success !");
    
    delete electrolyte;
}

// ------------------------------------------------------------
void
NafionTest::testSorptionIsothermLiu()
{
    electrolyte = new FuelCellShop::Material::Nafion;
    electrolyte->set_p_t(101325.);
    electrolyte->modify_parameters(5, "Liu09");
    
    // Vector return method
    std::vector<double> lambda_eq_vec;
    std::vector<double> xwater_vec;   
    std::vector<double> T_vec;
    xwater_vec.push_back(0.4); xwater_vec.push_back(0.7);
    T_vec.push_back(353.); T_vec.push_back(360.);
    electrolyte->set_water_molar_fraction( FuelCellShop::SolutionVariable(&xwater_vec, water_molar_fraction) );
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(&T_vec, temperature_of_REV) );
    electrolyte->sorption_isotherm(lambda_eq_vec);
    TEST_ASSERT_DELTA_MSG(lambda_eq_vec[0], 8.52611945404268, 1.e-6, "sorption_isotherm Liu method1 success !");
    TEST_ASSERT_DELTA_MSG(lambda_eq_vec[1], 12.59208184, 1.e-6, "sorption_isotherm Liu method2 success !");
    
    // Derivatives
    std::vector<VariableNames> deriv_flags;
    deriv_flags.push_back(temperature_of_REV);
    deriv_flags.push_back(water_molar_fraction);
    electrolyte->set_derivative_flags(deriv_flags);
    std::map< VariableNames, std::vector<double> > dlambda_eq;
    electrolyte->sorption_isotherm_derivative(dlambda_eq);
    
    double delta(std::pow(10.,-8));
    std::vector<double> lambdaeq_delta_vec;
    std::vector<double> T_delta_vec;
    T_delta_vec.push_back(353.+delta); T_delta_vec.push_back(360.+delta);
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(&T_delta_vec, temperature_of_REV) );
    electrolyte->sorption_isotherm(lambdaeq_delta_vec);
    TEST_ASSERT_DELTA_MSG((lambdaeq_delta_vec[0]-lambda_eq_vec[0])/delta, dlambda_eq[temperature_of_REV][0], 1.e-5, "sorption_isotherm Liu temperature derivative method1 success !");
    TEST_ASSERT_DELTA_MSG((lambdaeq_delta_vec[1]-lambda_eq_vec[1])/delta, dlambda_eq[temperature_of_REV][1], 1.e-6, "sorption_isotherm Liu temperature derivative method2 success !");
    
    std::vector<double> xwater_delta_vec;
    xwater_delta_vec.push_back(0.4+delta); xwater_delta_vec.push_back(0.7+delta);
    electrolyte->set_water_molar_fraction( FuelCellShop::SolutionVariable(&xwater_delta_vec, water_molar_fraction) );
    electrolyte->set_temperature( FuelCellShop::SolutionVariable(&T_vec, temperature_of_REV) );
    electrolyte->sorption_isotherm(lambdaeq_delta_vec);
    TEST_ASSERT_DELTA_MSG((lambdaeq_delta_vec[0]-lambda_eq_vec[0])/delta, dlambda_eq[water_molar_fraction][0], 1.e-5, "sorption_isotherm Liu water molar fraction derivative method1 success !");
    TEST_ASSERT_DELTA_MSG((lambdaeq_delta_vec[1]-lambda_eq_vec[1])/delta, dlambda_eq[water_molar_fraction][1], 1.e-6, "sorption_isotherm Liu water molar fraction derivative method2 success !");
    
    
    delete electrolyte;
}