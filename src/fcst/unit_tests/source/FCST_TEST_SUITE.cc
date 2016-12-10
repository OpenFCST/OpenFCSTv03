#include <FCST_TEST_SUITE.h>



bool FcstTestSuite::run_tests()
{
    Test::Suite ts;

    //Add sub tests suites here
    ts.add(std::auto_ptr<Test::Suite>(new FcstsUnitsTestSuite));
    ts.add(std::auto_ptr<Test::Suite>(new AnalyticalAgglomerateTest));
    ts.add(std::auto_ptr<Test::Suite>(new FCSTdatabaseTest));
    ts.add(std::auto_ptr<Test::Suite>(new MultiScaleCLTest));
    ts.add(std::auto_ptr<Test::Suite>(new PlatinumTest));
    ts.add(std::auto_ptr<Test::Suite>(new SolutionVariableTest));
    ts.add(std::auto_ptr<Test::Suite>(new NafionTest));
    ts.add(std::auto_ptr<Test::Suite>(new WaterAgglomerateTest));
    ts.add(std::auto_ptr<Test::Suite>(new IonomerAgglomerateTest));
    ts.add(std::auto_ptr<Test::Suite>(new UtilsTest));
    ts.add(std::auto_ptr<Test::Suite>(new LiquidWaterTest));
    ts.add(std::auto_ptr<Test::Suite>(new ConventionalCLTest));
    ts.add(std::auto_ptr<Test::Suite>(new DesignFibrousGDLTest));
    ts.add(std::auto_ptr<Test::Suite>(new DesignMPLTest));
    ts.add(std::auto_ptr<Test::Suite>(new EnumerationTest));
    ts.add(std::auto_ptr<Test::Suite>(new PuregasoxygenTest));
    ts.add(std::auto_ptr<Test::Suite>(new DoubleTrapTest));
    ts.add(std::auto_ptr<Test::Suite>(new GasMixtureTest));
    ts.add(std::auto_ptr<Test::Suite>(new PSD_HI_Test));
    ts.add(std::auto_ptr<Test::Suite>(new PSD_HO_Test));
    //ts.add(std::auto_ptr<Test::Suite>(new WaterPoreAgglomerateTest)); //under development
    ts.add(std::auto_ptr<Test::Suite>(new NumericalAgglomerateBaseTest));
    ts.add(std::auto_ptr<Test::Suite>(new PorousLayerTest)); ///under development
    //ts.add(std::auto_ptr<Test::Suite>(new FuelCell::UnitTest::FEVectorsTest)); ///under development
    ts.add(std::auto_ptr<Test::Suite>(new FuelCell::UnitTest::ApplicationStep3Test));
    ts.add(std::auto_ptr<Test::Suite>(new FuelCell::UnitTest::ApplicationStep8Test));    
    
    Test::TextOutput output(Test::TextOutput::Verbose);
    return ts.run(output);

}

