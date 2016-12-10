//---------------------------------------------------------------------------
//
//    Copyright (C) 2009 by Marc Secanell, University of Alberta
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__FIBER_BASE_H
#define _FUELCELLSHOP__FIBER_BASE_H

// Include FCST classes
#include <materials/base_material.h>


namespace FuelCellShop
{
	namespace Material
	{
		class FiberBase
		:
		public BaseMaterial
		{
		public:
			
			/** Constructor 
			 * The constructor initialize parameters using the default values. This is
			 * so that if I do not want to call declare_parameters and initialize, I can
			 * still use the routine with the hard coded values.
			 * 
			 * \todo1 What are we going to do with the Carbon, Nafion, Platinum, etc. classes? i.e. puresolids
			 */
			FiberBase(std::string name)
			:
			BaseMaterial(name)
			{};
			
			/**  Destructor  */
			~FiberBase()
			{};
			
			/** Declare parameters*/
			virtual void declare_parameters(ParameterHandler &param) const
			{
				FuelCellShop::Material::BaseMaterial::declare_parameters(param);
			};
			
			/** Initialize */
			virtual void initialize (ParameterHandler &param)
			{
				FuelCellShop::Material::BaseMaterial::initialize(param);
			};
			
			/** Obtain the electrical conductivity */
			virtual inline double get_electrical_conductivity() const
			{
				const std::type_info& info = typeid(*this);
				FcstUtilities::log << "Pure function " << __FUNCTION__
				<< " called in Class "
				<< info.name()  << std::endl;
				return 0;
			};
			
			/** Obtain the temperature dependent electrical conductivity which is passed as the argument */
			virtual inline double get_electrical_conductivity(double) const
			{
				const std::type_info& info = typeid(*this);
				FcstUtilities::log << "Pure function " << __FUNCTION__
				<< " called in Class "
				<< info.name()  << std::endl;
				return 0;
			};
			
			virtual double get_derivative_electrical_conductivity() const
			{
				const std::type_info& info = typeid(*this);
				FcstUtilities::log << "Pure function " << __FUNCTION__
				<< " called in Class "
				<< info.name()  << std::endl;
				return 0;
			};
			
			/** Obtain the thermal conductivity */
			virtual inline double get_thermal_conductivity() const
			{
				const std::type_info& info = typeid(*this);
				FcstUtilities::log << "Pure function " << __FUNCTION__
				<< " called in Class "
				<< info.name()  << std::endl;
				return 0;
			};
			
			virtual double get_derivative_thermal_conductivity() const
			{
				const std::type_info& info = typeid(*this);
				FcstUtilities::log << "Pure function " << __FUNCTION__
				<< " called in Class "
				<< info.name()  << std::endl;
				return 0;
			};
			
			/** Obtain the density */
			virtual inline double get_density() const
			{
				const std::type_info& info = typeid(*this);
				FcstUtilities::log << "Pure function " << __FUNCTION__
				<< " called in Class "
				<< info.name()  << std::endl;
				return 0;
			};   
			
		protected:
			/** Electrical conductivity of carbon fibers extrapolated to 100% solid phase */
			double electrical_conductivity;
			/** Thermal conductivity of carbon fibers extrapolated to 100% solid phase */
			double thermal_conductivity;
			/** Density of carbon fibers*/
			double density;
		};
	}
}

#endif
