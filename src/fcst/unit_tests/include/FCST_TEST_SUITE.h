
/**
 * If you wish to add a unit testing suite please see the example fcst_units_test .cc and .h,
 * when you've finished your test suite add it to FCST_TEST_SUITE.cc. To run the unit tests
 * set Run test = true in the Simulator subsection of your param file.
 *
 * Assertions that can are provided by cpp test for you to build your tests with:
 *
 *@code TEST_FAIL(msg) @endcode
 *Unconditional failure. For instance if assesing a flag using an if statement fails one
 *can use a TEST_FAIL to inform the user.
 *
 *@code TEST_ASSERT(expr) @endcode
 *Verify an expression and issues an assertment if it fails.
 *
 *@code TEST_ASSERT_MSG(expr, msg) @endcode
 *Verify an expression and issues user message and assertment if it fails.
 *
 *@code  TEST_ASSERT_DELTA(a, b, delta) @endcode
 *Verify that two expressions are equal up to a constant, issues an assertment if it fails.
 *
 *@code  TEST_ASSERT_DELTA_MSG(a, b, delta, msg)@endcode
 *Verify that two expressions are equal up to a constant, issues an assertment and a user message if it fails.
 *
 *@code TEST_THROWS(expr, x)@endcode
 *Verify an expression and expects an exception in return. An assertment is issued if the exception is not thrown.
 *
 *@code TEST_THROWS_MSG(expr, x, msg)@endcode
 *Verify an expression and expects an exception in return. An assertment and user message are issued if the exception is not thrown.
 *
 *@code  TEST_THROWS_ANYTHING(expr)@endcode
 *Verify an expression and expects any exception in return. An assertment is issued if no exception is thrown.
 *
 *@code TEST_THROWS_ANYTHING_MSG(expr, msg)@endcode
 *Verify an expression and expects any exception in return. An assertment and user message are issued if no exception is thrown.
 *
 *@code TEST_THROWS_NOTHING(expr)@endcode
 *Verify an expression and expects no exception in return. An assertment is issued if any exception is thrown.
 *
 *@code TEST_THROWS_NOTHING_MSG(expr, msg) @endcode
 * Verify an expression and expects no exception in return. An assertment and user message are issued if any exception is thrown.
 *
 */




#ifndef _FCST_TEST_SUITE
#define _FCST_TEST_SUITE

//List of the sub suites
#include <fcst_units_test.h>
#include <analytical_agglomerate_test.h>
#include <fcst_db_test.h>
#include <agglomerate_catalyst_layer_test.h>
#include <platinum_test.h>
#include <solution_variable_test.h>
#include <nafion_test.h>
#include <water_agglomerate_test.h>
#include <ionomer_agglomerate_test.h>
#include <utils_test.h>
#include <cpptest.h>
#include <LiquidWater_test.h>
#include <conventional_catalyst_layer_test.h>
#include <design_fibrous_gdl_test.h>
#include <design_mpl_test.h>
#include <enumeration_test.h>
#include <puregas_oxygen_test.h>
#include <DT_test.h>
#include <GasMixture_test.h>
#include <water_pore_agglomerate_test.h>
#include <numerical_agglomerate_base_test.h>
#include <porous_layer_test.h>
#include <PSD_HI_test.h>
#include <PSD_HO_test.h>
//#include <fevectors_test.h>
#include <application_step3_test.h>
#include <application_step8_test.h>

namespace FcstTestSuite
{
    bool run_tests();

} //namespace FcstTestSuite
#endif
