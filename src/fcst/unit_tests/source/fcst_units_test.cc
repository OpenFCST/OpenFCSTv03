
#include <fcst_units_test.h>

void FcstsUnitsTestSuite::perBigToSmallTest()
{
	// 1/m^2 to 1/cm^2
	static double answer = Units::convert(1, Units::PER_C_UNIT2, Units::PER_UNIT2);
	static double expected_answer =  1E-4;
	TEST_ASSERT(answer == expected_answer);
	//TEST_ASSERT_DELTA(0.5, 0.7, 0.3);
}

void FcstsUnitsTestSuite::bigToSmallTest()
{
	//m^2/a to cm^2/a: 
	static double answer = Units::convert(1, Units::C_UNIT2, Units::UNIT2);
	static double expected_answer = 10000;
	
	TEST_ASSERT(answer == expected_answer);
}

void FcstsUnitsTestSuite::perSmallToBig()
{
	//1/cm^3 to 1/m^3
	static double answer = Units::convert(1, Units::PER_UNIT3, Units::PER_C_UNIT3);
	static double expected_answer = 1000000;
	TEST_ASSERT(answer == expected_answer);
}

void FcstsUnitsTestSuite::smallToBig()
{
	//cm^3/a to m^3/a
	static double answer = Units::convert(1, Units::UNIT3, Units::C_UNIT3);
	static double expected_answer =  1E-6;
	TEST_ASSERT(answer == expected_answer);
	

}

void FcstsUnitsTestSuite::btuToKwh()
{
	
	TEST_ASSERT(Units::convert(1,Units::BTU_to_KJ) == 1.054);

}

void FcstsUnitsTestSuite::kwhToBtu()
{
	TEST_ASSERT(Units::convert(1,Units::KJ_to_BTU) == 1/1.054);

}

void FcstsUnitsTestSuite::testExceptions()
{
	TEST_THROWS(Units::convert(1,Units::PER_C_UNIT, Units::C_UNIT), std::logic_error);
	TEST_THROWS(Units::convert(1,-1), std::invalid_argument);

}

