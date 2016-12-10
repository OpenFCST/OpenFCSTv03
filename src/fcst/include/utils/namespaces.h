//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2011-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: namespaces.h
//    - Description: This file is used to provide DOxygen with information on namespaces.
//    - Developers: P. Dobson and M. Secanell
//
//---------------------------------------------------------------------------


#ifndef _FUEL_CELL__H
#define _FUEL_CELL__H

/**
 * \file 
 * \brief General namespace for Fuel Cell Simulation Toolbox applications and utilities
 * 
 * @author M. Secanell, P. Dobson - &copy;2011
 * 
 */
namespace FuelCell
{
    /**
     * @namespace FuelCell
     * 
     * \brief General namespace for all applications to built within FCST 
     *
     * This namespace is used to declare the application class and all necessary functions
     * and data required to solve the problem.  Application classes will also typically house all
     * data handling and post-processing routines.
     * 
     * In a typical case using FEM, the virtual functions from OptimizationBlockMatrixApplication
     * and its inheritance structure will be re-implemented to specifically suit the problem.
     * See \ref AppCathode for an example of the implementation.
     */
    namespace Application
    {
        
    }
    
    /**
     * 
     * \brief General namespace for all initial solution classes (deal.II functions) 
     *
     * This namespace is used to declare the initial solutions for all applications which
     * are linear approximations of nonlinear problems.  
     * These are mainly classes that inherit from Function<dim> and that are necessary in order to use
     * initialization subroutines in deal.II. For example InitialSolution::PemfcIC is created in order to use the
     * deal.II class VectorInterpolate which in turn is used to set up the initial solution to the problem. 
     */
    namespace InitialSolution
    {
        
    }
    
}

/**
 * \brief Container Namespace for Fuel Cell data classes
 * 
 * This namespace should not hold any member classes and will characterize sub-namespaces for fuel cell data
 *   such as Layers, Materials, Mixtures, etc.
 */
namespace FuelCellShop
{
    /**
     * \brief Namespace containing all universal constants
     *
     * This namespace contains only external inline functions that return
     * the value of a universal constant.
     */
    namespace Constants
    {}

     /**
     * 
     * \brief General namespace to generate geometries and grids for FuelCell applications 
     *
     * This namespace is used to generate a grid on the triangulation object.  General geometry information
     * is stored in the \ref GridBase class while the specific information to each geometry and application is
     * stored in the inherited classes.
     */
    namespace Geometry
    {
        
    }
    
    /** 
     * \brief Namespace to hold classes that characterize materials used in fuel cells
     *
     * \todo3 Create a liquid water class (for use in the agglomerate, UTCL?, CL w/ PSD?)
     */ 
    namespace Material
    {
        
    }
    
    /**
     * \brief Namespace to hold classes that characterize fuel cell layers
     * 
     * \todo3 Create a macrohomogeneous with pore size distribution CL class (inherits from macrohomogeneous catalyst layer).
     */ 
    namespace Layer
    {
        
    }
    
    /** 
     * \brief Namespace to hold classes that define mixture properties for gases (and liquids) in fuel cells
     * 
     * This namespace should be used to take objects or data from the material classes and define mixture properties
     *	that characterize a substance, eg. humidified air.
     */ 
    namespace Mixture
    {
        
    }
    
    /** 
     * \brief Namespace to hold classes that describe physical processes in fuel cells.
     * 
     * This namespace should be used to take objects from the layer classes, compute the coefficients for the physical equations, 
     * 	and assemble the matrices and residuals.
     */ 
    namespace Equation
    {
        
    }
    
    /**
     * \brief Namespace to hold classes that implement kinetic models. Currently implemented is 
     * Tafel kinetics, with Wang's kinetics under construction
     * 
     */ 
    namespace Kinetics
    {
        
    }	
}

#endif