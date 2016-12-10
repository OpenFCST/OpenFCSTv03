//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: utils_test.cc
//    - Description: Unit testing class for FCST Utility classes
//    - Developers: Philip Wardlaw
//    - Id: $Id: utils_test.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include "utils_test.h"

void
UtilsTest::testIsNumber(){

    //Correct string representation of numbers
    TEST_ASSERT(FcstUtilities::is_number("-1.23"));
    TEST_ASSERT(FcstUtilities::is_number("22"));
    TEST_ASSERT(FcstUtilities::is_number("1231123"));

    //Incorrect
    TEST_ASSERT(!FcstUtilities::is_number("-+2.1"));
    TEST_ASSERT(!FcstUtilities::is_number("2+.1"));
    TEST_ASSERT(!FcstUtilities::is_number(".+1"));
    TEST_ASSERT(!FcstUtilities::is_number("-2.65.23"));
    TEST_ASSERT(!FcstUtilities::is_number("fred"));

}

//================================================
//================================================
void 
UtilsTest::testModify_parameter_file_double()
{
    ParameterHandler param;
    
    param.enter_subsection("Fuel cell data"); 
    {
        param.enter_subsection("Cathode gas diffusion layer"); 
        {

            param.enter_subsection("DesignFibrousGDL"); 
            { 
                //-- Composition
                param.declare_entry ("Porosity", // mass percentage of platinum catalyst on the support carbon black
                                     "0.46", // [-]
                                     Patterns::Double());
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    std::string name_design_var("Fuel cell data>>Cathode gas diffusion layer>>DesignFibrousGDL>>Porosity"); 
    double expectedAnswer (0.2);
    
    FcstUtilities::modify_parameter_file(name_design_var, expectedAnswer, param);
    
    double answer(0);
    
    param.enter_subsection("Fuel cell data"); 
    {
        param.enter_subsection("Cathode gas diffusion layer"); 
        {
            param.enter_subsection("DesignFibrousGDL"); 
            { 
                //-- Composition
                answer = param.get_double("Porosity");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();    
    
    TEST_ASSERT_MSG(expectedAnswer == answer, "testModify_parameter_file_double failed! You loose :(");
}

//================================================
//================================================
void 
UtilsTest::testModify_parameter_file()
{
    ParameterHandler param;
    
    param.enter_subsection("Fuel cell data"); 
    {
        param.enter_subsection("Cathode catalyst layer"); 
        {
            param.enter_subsection("ConventionalCL"); 
            { 
                //-- Composition
                param.declare_entry ("Platinum loading on support (%wt)", // mass percentage of platinum catalyst on the support carbon black
                                     "4:0.46", // [-]
                                     Patterns::Map( Patterns::Integer(0,255), Patterns::Double(0) ) ,
                                     "Mass percentage of platinum catalyst on the support carbon black");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    std::string name_design_var("Fuel cell data>>Cathode catalyst layer>>ConventionalCL>>Platinum loading on support (%wt):4"); 
    double expectedAnswer (0.2);
    
    FcstUtilities::modify_parameter_file(name_design_var, expectedAnswer, param);
    
    std::map<unsigned int, double> answer;
    
    param.enter_subsection("Fuel cell data"); 
    {
        param.enter_subsection("Cathode catalyst layer"); 
        {
            param.enter_subsection("ConventionalCL"); 
            { 
                //-- Composition
                answer = FcstUtilities::string_to_map<unsigned int, double>( Utilities::split_string_list( param.get ("Platinum loading on support (%wt)")));
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();    
    
    TEST_ASSERT_MSG(expectedAnswer == answer[4], "testModify_parameter_file failed! You loose :(");
}

//================================================
//================================================
void 
UtilsTest::testModify_parameter_file_list()
{
    ParameterHandler param;
    
    param.enter_subsection("Fuel cell data"); 
    {
        param.enter_subsection("Cathode catalyst layer"); 
        {
            param.enter_subsection("ConventionalCL"); 
            { 
                //-- Composition
                param.declare_entry ("Platinum loading on support (%wt)", // mass percentage of platinum catalyst on the support carbon black
                                     "4:0.46,5:0.6", // [-]
                                     Patterns::Map( Patterns::Integer(0,255), Patterns::Double(0) ) ,
                                     "Mass percentage of platinum catalyst on the support carbon black");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    std::string name_design_var("Fuel cell data>>Cathode catalyst layer>>ConventionalCL>>Platinum loading on support (%wt):4"); 
    double expectedAnswer (0.2);
    
    FcstUtilities::modify_parameter_file(name_design_var, expectedAnswer, param);
    
    std::string name_design_var2("Fuel cell data>>Cathode catalyst layer>>ConventionalCL>>Platinum loading on support (%wt):5"); 
    double expectedAnswer2 (0.3);
    
    FcstUtilities::modify_parameter_file(name_design_var2, expectedAnswer2, param);
    
    std::map<unsigned int, double> answer;
    
    param.enter_subsection("Fuel cell data"); 
    {
        param.enter_subsection("Cathode catalyst layer"); 
        {
            param.enter_subsection("ConventionalCL"); 
            { 
                //-- Composition
                answer = FcstUtilities::string_to_map<unsigned int, double>( Utilities::split_string_list( param.get ("Platinum loading on support (%wt)")));
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();    
    
    TEST_ASSERT_MSG(expectedAnswer == answer[4], "testModify_parameter_file_list failed! You loose :(");
    TEST_ASSERT_MSG(expectedAnswer2 == answer[5], "testModify_parameter_file_list failed! You loose :(");
}