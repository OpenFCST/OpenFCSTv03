//---------------------------------------------------------------------------
//    
//    Fuel Cell Simulation Toolbox (FCST) information file	
//
//    Copyright (C) 2006-2012 by Marc Secanell
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

This folder contains information regaring the pre-processing capabilities in OpenFCST. 

FCST has a built-in grid generator set of classes under namespace FuelCellShop::Geometry. The
parent class of the built-in grid genertor is in FuelCellShop::Geometry::GridBase< dim >.
If you would like to use the built-in grid generator, please take a look at the FCST user guide for an example.

FCST is also capable of reading meshes from GMesh and from Salome. Salome is the external mesh generator
of choice within the Energy Systems Design Laboratory. In this folder, several examples of Salome
meshes are shown. A tutorial on how to use the Salome mesh generator is provided in the FCST Users Guide.
Salome is an open-source mesh generator and can be downloaded free of charge at http://www.salome-platform.org/

