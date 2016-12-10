/*
 * DT_test.h
 *
 *  Created on: Feb 11, 2014
 *      Author: wardlawp
 */

#ifndef DT_TEST_H_
#define DT_TEST_H_



#include <cpptest.h>
#include <string>
#include <application_core/system_management.h>
#include <reactions/base_kinetics.h>
#include <materials/catalyst_base.h>
#include <materials/polymer_electrolyte_material_base.h>
#include <map>

class DoubleTrapTest: public Test::Suite
{
public:
    DoubleTrapTest()
    {
       //Add a number of tests that will be called during Test::Suite.run()
        TEST_ADD(DoubleTrapTest::testDerivatives);
        TEST_ADD(DoubleTrapTest::testCurrent);


    }
protected:
    virtual void setup(); // setup resources... called before Test::Suite.run()
    virtual void tear_down() {} // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
private:

    //Necessary objects for kinetics
    boost::shared_ptr<FuelCellShop::Kinetics::BaseKinetics > kinetics;
    boost::shared_ptr<FuelCellShop::Material::CatalystBase > catalyst;
    boost::shared_ptr<FuelCellShop::Material::PolymerElectrolyteBase > electrolyte;

    //Param object for initialization
    ParameterHandler param;

    //Typedefs for convenience
    typedef FuelCellShop::SolutionVariable SOL;
    typedef std::vector<FuelCellShop::SolutionVariable> SOLS;

    //Tests
    void testDerivatives();
    void testCurrent();


};



#endif /* DT_TEST_H_ */
