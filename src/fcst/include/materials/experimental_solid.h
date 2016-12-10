// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: experimental_solid.h
// - Description: This class describes a solid
// - Developers: Valentin N. Zingan, University of Alberta
// - Id: $Id: experimental_solid.h 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#ifndef _FCST_FUELCELLSHOP_MATERIAL_EXPERIMENTAL_SOLID_H_
#define _FCST_FUELCELLSHOP_MATERIAL_EXPERIMENTAL_SOLID_H_

#include <materials/base_material.h>

namespace FuelCellShop
{
    namespace Material
    {

        /**
        * This class describes
        * a solid.
        *
        * The functionality of
        * this class can be extended
        * if needed.
        *
        * \author Valentin N. Zingan, 2012
        */

        class ExperimentalSolid : public BaseMaterial
        {
            public:

                ///@name Constructors, destructor, and initialization
                //@{

                /**
                * Constructor.
                */
                ExperimentalSolid(const std::string& name);

                /**
                * Destructor.
                */
                virtual ~ExperimentalSolid();

                /**
                * Declare parameters.
                */
                virtual void declare_parameters(ParameterHandler& param) const;

                /**
                * Initialize parameters.
                */
                virtual void initialize(ParameterHandler& param);

                //@}

        };

    } // Material

} // FuelCellShop

#endif