//---------------------------------------------------------------------------
//
//    Copyright (C) 2011 by Marc Secanell, University of Alberta
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__CARBON_FIBER_H
#define _FUELCELLSHOP__CARBON_FIBER_H

// Include FCST classes
#include <materials/fiber_base.h>


namespace FuelCellShop
{
	namespace Material
	{
		class CarbonFiber :
		public FiberBase
		{
		public:
			
			/** Constructor 
			 * The constructor initialize parameters using the default values. This is
			 * so that if I do not want to call declare_parameters and initialize, I can
			 * still use the routine with the hard coded values.
			 */
			CarbonFiber(std::string name = "Carbon Fiber");
			
			/**  Destructor  */
			~CarbonFiber();
			
			/** Declare parameters*/
			void declare_parameters(ParameterHandler &param) const;
			
			/** Initialize */
			void initialize (ParameterHandler &param);
			
			/** Obtain the electrical conductivity */
			inline double get_electrical_conductivity() const
			{ return this->electrical_conductivity; };
			
			/** Obtain the electrical conductivity */
			double get_electrical_conductivity(double) const;
			
			double get_derivative_electrical_conductivity() const
			{
				const std::type_info& info = typeid(*this);
				FcstUtilities::log << "Pure function " << __FUNCTION__
				<< " called in Class "
				<< info.name()  << std::endl;
				return 0;
			}
			
			/** Obtain the thermal conductivity */
			inline double get_thermal_conductivity() const
			{ return this->thermal_conductivity; };
			
			double get_derivative_thermal_conductivity() const
			{
				const std::type_info& info = typeid(*this);
				FcstUtilities::log << "Pure function " << __FUNCTION__
				<< " called in Class "
				<< info.name()  << std::endl;
				return 0;
			}		
			
			/** Obtain the density */
			inline double get_density() const
			{ return this->density; };
			
			
		};
	}
}

#endif
