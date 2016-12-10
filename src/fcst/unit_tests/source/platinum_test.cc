#include "platinum_test.h"

/**
 * Setup is called before any of the tests are called.
 */
void
PlatinumTest::setup()
{
    
    
    FuelCellShop::Material::CatalystBase::declare_Catalyst_parameters(param);    
   
    catalyst = FuelCellShop::Material::CatalystBase::create_Catalyst(param, "Platinum");    
    
    
}

void
PlatinumTest::tear_down() 
{
    //delete catalyst;
}



void
PlatinumTest::testSetSolution()
{
        
    //TEST_ASSERT_MSG(false, "Test not yet implemented!")
}


void
PlatinumTest::testAlphaAnodic()
{
    //We are testing the behaviour of Calculator::add(double a, double b)
    double answer(0);
    double expectedAnswer = 0.5;
    catalyst->set_reaction_kinetics(HOR);
    catalyst->alpha_anodic (answer);  
    
    std::ostringstream streamOut;
    streamOut <<"The value of the alpha_anodic is: "<<answer<<". The expected value is: "<<expectedAnswer<<std::endl;
    streamOut <<"testAlphaAnodic failed!"<<std::endl;

    TEST_ASSERT_MSG(expectedAnswer == answer, streamOut.str().c_str());
    
}

void
PlatinumTest::testReferenceConcentration()
{
    catalyst->set_reaction_kinetics(ORR);
    
    std::vector<VariableNames > reactants;
    reactants.push_back(oxygen_concentration);
    
    std::map< VariableNames, double > answer;
    catalyst->reference_concentration (reactants, answer);

    std::vector<double> expectedAnswer;
    expectedAnswer.push_back(1.6e-05);
   
    std::ostringstream streamOut;
    streamOut <<"The value of the reference_concentration is: "<< answer[oxygen_concentration] <<". The expected value is: "<< expectedAnswer[0]<<std::endl;
    streamOut <<"testReferenceConcentration failed!"<<std::endl;

    TEST_ASSERT_DELTA_MSG(expectedAnswer[0], answer[oxygen_concentration], 1e-06, streamOut.str().c_str());
    
}

void
PlatinumTest::testNeyerlin()
{
    ParameterHandler param2;
    
    FuelCellShop::Material::CatalystBase::declare_Catalyst_parameters(param2);    
    
    param2.enter_subsection("Materials");
    param2.enter_subsection("Platinum");
    param2.set("Method for kinetics parameters (ORR)","Neyerlin");
    param2.leave_subsection();
    param2.leave_subsection();
   
    boost::shared_ptr<FuelCellShop::Material::CatalystBase > catalyst2 = FuelCellShop::Material::CatalystBase::create_Catalyst(param2, "Platinum");    
    
    double answer(0);
    double expectedAnswer = 1;
    catalyst2->set_reaction_kinetics(ORR);
    
    // Test alpha
    {
        catalyst2->alpha_cathodic (answer); 
        std::ostringstream streamOut;
        streamOut <<"The value of the cathodic alpha is: "<<answer<<". The expected value is: "<<expectedAnswer<<std::endl;
        streamOut <<"testNeyerlin failed!"<<std::endl;
        
        TEST_ASSERT_MSG(expectedAnswer == answer, streamOut.str().c_str());
    }
    
    // Test gamma
    {
        expectedAnswer = 0.54;
        std::vector<VariableNames > reactants;
        reactants.push_back(oxygen_concentration);
        std::map< VariableNames, double > reaction_order;
        catalyst2->reaction_order (reactants, reaction_order);
        answer = reaction_order[oxygen_concentration];
        std::ostringstream streamOut;
        streamOut <<"The value of the reaction order is: "<<answer<<". The expected value is: "<<expectedAnswer<<std::endl;
        streamOut <<"testNeyerlin failed!"<<std::endl;
        
        TEST_ASSERT_MSG(expectedAnswer == answer, streamOut.str().c_str());
    }
    
    // Test exchange current density
    
}


void
PlatinumTest::testParthasarathy()
{
    ParameterHandler param2;
    
    FuelCellShop::Material::CatalystBase::declare_Catalyst_parameters(param2);    
    
    param2.enter_subsection("Materials");
    param2.enter_subsection("Platinum");
    param2.set("Method for kinetics parameters (ORR)","Parthasarathy");
    param2.leave_subsection();
    param2.leave_subsection();
   
    boost::shared_ptr<FuelCellShop::Material::CatalystBase > catalyst2 = FuelCellShop::Material::CatalystBase::create_Catalyst(param2, "Platinum");    
    
    double answer(0);
    double expectedAnswer = 1.0;
    catalyst2->set_reaction_kinetics(ORR);
    
    // Test alpha
    {
        catalyst2->alpha_cathodic (answer); 
        std::ostringstream streamOut;
        streamOut <<"The value of the cathodic_anodic is: "<<answer<<". The expected value is: "<<expectedAnswer<<std::endl;
        streamOut <<"testParthasarathy failed!"<<std::endl;
        
        TEST_ASSERT_MSG(expectedAnswer == answer, streamOut.str().c_str());
    }
    
    // Test gamma
    {
        expectedAnswer = 1.0;
        std::vector<VariableNames > reactants;
        reactants.push_back(oxygen_concentration);
        std::map< VariableNames, double > reaction_order;
        catalyst2->reaction_order (reactants, reaction_order);
        answer = reaction_order[oxygen_concentration];
        std::ostringstream streamOut;
        streamOut <<"The value of the reaction order is: "<<answer<<". The expected value is: "<<expectedAnswer<<std::endl;
        streamOut <<"testParthasarathy failed!"<<std::endl;
        
        TEST_ASSERT_MSG(expectedAnswer == answer, streamOut.str().c_str());
    }
}
    
void
PlatinumTest::testParthasarathyHCD()
{
    ParameterHandler param2;
    
    FuelCellShop::Material::CatalystBase::declare_Catalyst_parameters(param2);    
    
    param2.enter_subsection("Materials");
    param2.enter_subsection("Platinum");
    param2.set("Method for kinetics parameters (ORR)","Parthasarathy_hcd");
    param2.leave_subsection();
    param2.leave_subsection();
   
    boost::shared_ptr<FuelCellShop::Material::CatalystBase > catalyst2 = FuelCellShop::Material::CatalystBase::create_Catalyst(param2, "Platinum");    
    
    double answer(0);
    double expectedAnswer = 0.5;
    catalyst2->set_reaction_kinetics(ORR);
    
    // Test alpha
    {
        catalyst2->alpha_cathodic (answer); 
        std::ostringstream streamOut;
        streamOut <<"The value of the cathodic transfer coefficient is: "<<answer<<". The expected value is: "<<expectedAnswer<<std::endl;
        streamOut <<"testParthasarathyHCD failed!"<<std::endl;
        
        TEST_ASSERT_MSG(expectedAnswer == answer, streamOut.str().c_str());
    }
    
    // Test gamma
    {
        expectedAnswer = 1.0;
        std::vector<VariableNames > reactants;
        reactants.push_back(oxygen_concentration);
        std::map< VariableNames, double > reaction_order;
        catalyst2->reaction_order (reactants, reaction_order);
        answer = reaction_order[oxygen_concentration];
        std::ostringstream streamOut;
        streamOut <<"The value of the reaction order is: "<<answer<<". The expected value is: "<<expectedAnswer<<std::endl;
        streamOut <<"testParthasarathyHCD failed!"<<std::endl;
        
        TEST_ASSERT_MSG(expectedAnswer == answer, streamOut.str().c_str());
    }    
}