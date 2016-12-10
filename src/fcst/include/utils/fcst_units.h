//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class:  fcst_units.h
//    - Description: A cpp object that can be used by fcst programmers
//          when manipulating units in order to preserve
//          standards and prevent confusion.
//    - Developers: Philip Wardlaw and M. Secanell
//
//---------------------------------------------------------------------------

#include <string>
#include <iostream>
#include <stdexcept>

#ifndef _FCST__UNITS
#define _FCST__UNITS

/** 
 * Class used to convert units using a standard convention. 
 * All routines are static, therefore, there is no need to create an object of this class, simply
 * use the Units::convert() to call the class and add the .h file to your header file.
 * 
 * The converstion routine convert(unitToConvert, X, Y) should be read as:
 * 
 * convert quantity unitToConvert to X from Y.
 * 
 * <h3>Usage details</h3>
 * 
 * @code
 *	FcstUtilities::log << "Tests" << std::endl;
 *	FcstUtilities::log << "1/m^2 to 1/cm^2: " <<  Units::convert(1, Units::PER_C_UNIT2, Units::PER_UNIT2)<<std::endl;
 *	FcstUtilities::log << "1/cm^3 to 1/m^3: " <<Units::convert(1, Units::PER_UNIT3, Units::PER_C_UNIT3)<<  std::endl;
 *	FcstUtilities::log << "m^2/a to cm^2/a: " <<  Units::convert(1, Units::C_UNIT2, Units::UNIT2)<<std::endl;
 *	FcstUtilities::log << "cm^3/a to m^3/a: " <<Units::convert(1, Units::UNIT3, Units::C_UNIT3)<<  std::endl;
 * @endcode 
 * 
 * @author P. Wardlaw
 * 
 * @date 2013
 */

class Units {
    
public:
    /**
     *  The conversion routine convert(unitToConvert, X, Y) should be read as:
     * 
     *  convert quantity unitToConvert from Y to X
     */
    inline static double convert(double unitToConvert, double to,  double from){
        
        if (to > 0 && from > 0){
            
            return unitToConvert*(to/from);
            
        }
        else if (to < 0 && from < 0){
            
            return unitToConvert*(from/to);
            
        }
        else {
            //invalid use of the convert function. To and from must be of the same sign
            throw std::logic_error("Incorrect unit conversion.");
        }        
    }
    
    /**
     *  Convert using on a case by case basis. In this case, the function should be read as:
     * Convert unitToConvert from KJ _to_ BTU.
     * 
     */
    inline static double convert(double unitToConvert, int specificCase){
        
        switch (specificCase)
        {
            case KJ_to_BTU:
                return unitToConvert/1.054;
            case BTU_to_KJ:
                return unitToConvert*1.054;
            case ATM_to_PA:
                return unitToConvert*1.01325e5;
        }
        
        //Invalid case, throw error
        throw std::invalid_argument("Specific case not implemented");
    }
    
    //Static Integers denoting specific conversion cases (Typically Non Metric)
    static const unsigned int KJ_to_BTU =1;
    static const unsigned int BTU_to_KJ =2;
    static const unsigned int ATM_to_PA =3;
    
    //Static doubles describing generic conversions
    static double PER_K_UNIT; //1E3
    static double PER_UNIT; // 1;
    static double PER_C_UNIT; // 1E-2;
    static double PER_MILLI_UNIT; // 1E-3;
    static double PER_MICRO_UNIT; // 1E-6;
    static double PER_N_UNIT; // 1E-9;
    static double PER_P_UNIT; // 1E-12;
    
    static double PER_UNIT2; // 1;
    static double PER_C_UNIT2; // 1E-4;
    static double PER_MILLI_UNIT2; // 1E-6;
    static double PER_MICRO_UNIT2; // 1E-12;
    static double PER_N_UNIT2; // 1E-18;
    static double PER_P_UNIT2; // 1E-24;
    
    static double PER_UNIT3; // 1;
    static double PER_C_UNIT3; // 1E-6;
    static double PER_MILLI_UNIT3; // 1E-9;
    static double PER_MICRO_UNIT3; // 1E-18;
    static double PER_N_UNIT3; // 1E-27;
    static double PER_P_UNIT3; // 1E-36;
    
    //Sign is to differentiate between "Per unit" and "by unit"
    static double K_UNIT; // -1E3
    static double UNIT; // -1;
    static double C_UNIT; //- 1E-2;
    static double MILLI_UNIT; //- 1E-3;
    static double MICRO_UNIT; //- 1E-6;
    static double N_UNIT; //- 1E-9;
    static double P_UNIT; // -1E-12;
    
    static double UNIT2; //- 1;
    static double C_UNIT2; // -1E-4;
    static double MILLI_UNIT2; // -1E-6;
    static double MICRO_UNIT2; //- 1E-12;
    static double N_UNIT2; //- 1E-18;
    static double P_UNIT2; //- 1E-24;
    
    static double UNIT3; // -1;
    static double C_UNIT3; // -1E-6;
    static double MILLI_UNIT3; // -1E-9;
    static double MICRO_UNIT3; // -1E-18;
    static double N_UNIT3; //- 1E-27;
    static double P_UNIT3; // -1E-36;
};




#endif
