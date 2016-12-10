/**
 * A unit test class that tests the WaterConicalPoreAgglomerate class.
 *
 * Tests the individual protected functions that ensure public function
 * compute_current(E_r) works correctly.
 *
 *
 *
 */



#ifndef _water_pore_test
#define _water_pore_test

#include<agglomerate_water_sadeghi.h>

#include <cpptest.h>

#include <boost/shared_ptr.hpp>
#include <materials/catalyst_base.h>
#include <materials/catalyst_support_base.h>
#include <materials/polymer_electrolyte_material_base.h>
#include <materials/PureLiquid.h>
#include <application_core/fcst_variables.h>
#include <materials/PureLiquid.h>
#include <utils/fcst_utilities.h>


#include<deal.II/base/parameter_handler.h>



class WaterPoreAgglomerateTest: public Test::Suite
{
    public:
    WaterPoreAgglomerateTest()
    {
       //Add a number of tests that will be called during Test::Suite.run()
       //Generic cases

        TEST_ADD(WaterPoreAgglomerateTest::testSolvePotentials);
        TEST_ADD(WaterPoreAgglomerateTest::testSolveO2);
        TEST_ADD(WaterPoreAgglomerateTest::testCurrent);


    }
    protected:
        virtual void setup(); // setup resources...
        virtual void tear_down() {} // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
    private:
        FuelCellShop::MicroScale::WaterConicalPoreAgglomerate agg;
        boost::shared_ptr<FuelCellShop::Material::CatalystBase> catalyst;
        boost::shared_ptr<FuelCellShop::Material::CatalystSupportBase> support;
        boost::shared_ptr<FuelCellShop::Material::PolymerElectrolyteBase> electrolyte;
        FuelCellShop::Material::LiquidWater water;


        void testSolvePotentials();
        void testSolveO2();
        void testCurrent();

};

#endif
