//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//    - Class: homogeneous_CL.h
//    - Description: Class characterizing the macro-homogeneous catalyst layer.
//    - Developers: Marc Secanell, Peter Dobson and Madhur Bhaiya
//    - Id: $Id: homogeneous_CL.h 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__HOMOGENEOUS_CL__H
#define _FUELCELLSHOP__HOMOGENEOUS_CL__H

// Include openFCST routines:
#include <utils/fcst_constants.h>
#include <layers/conventional_CL.h>

namespace FuelCellShop
{
    namespace Layer
    {
        /**
         * This class characterizes a catalyst layer and uses this information 
         * to compute effective transport properties and interfacial areas for phase
         * change or electrochemical reactions.
         * 
         * This class implements a macrohomogeneous catalyst layer.
         *
         */
        template <int dim>
        class HomogeneousCL :
        public ConventionalCL<dim>
        {
        public:
            
            /**
             * Concrete name used for objects of this class. This name is used when
             * setting up the subsection where the data is stored in the input file.
             * 
             * The data will be store under
             * \code
             * subsection name_specified_in_constructor
             *    set Material id = 2
             *    set Catalyst layer type = DummyCL # <-here I select the type of object of type CatalystLayer
             *    subsection DummyCL # <- this is the concrete_name for this class
             *       set all info relevant to this object
             *    end
             * end
             * \endcode
             */
            static const std::string concrete_name;
            
            ///@name Constructors, destructor, and initalization
            //@{    
            /**
             * Prototye Constructor
             * 
             * \warning For internal use only
             */
            HomogeneousCL();      
            
            /**
             * Destructor
             */
            ~HomogeneousCL();
                      
            //@}
            
            ///@name Accessors and info
            //@{
            /**
             * This member function will use a FuelCellShop::BaseKinetics class in order to compute the current density
             * production in the CL
             */
            virtual void current_density(std::vector<double>&);
            
            /**
             * This member function computes the current density production in the CL. First argument is <b>current density</b>, and second is <b>effectiveness</b>, at 
             * all quadrature points in the cell. Since this is a macro-homogeneous layer, effectiveness is filled as \b 1.0
             */
            virtual void current_density(std::vector<double>& current, std::vector<double>& effectiveness)
            {
                current_density(current);
                effectiveness.assign(current.size(), 1.0);
            }
            
            /**
             * This member function will use a FuelCellShop::BaseKinetics class in order to compute the derivative of the current density
             * with respect to the variables setup using #set_derivative_flags method. 
             */
            virtual void derivative_current_density(std::map< VariableNames, std::vector<double> >& );
            //@}
            
        protected:
            ///@name Constructors and declarations
            //@{
            /**
             * Constructor
             * 
             */
            HomogeneousCL(const std::string name);                                     

            /**
             * Declare parameters for a parameter file.
             * 
             * \note This member function must be virtual since it will be accessed
             * via pointers for all children.
             */
            virtual void declare_parameters (const std::string& cl_section_name, 
                                             ParameterHandler &param) const;

            /**
             * Member function used to read in data and initialize the necessary data
             * to compute the coefficients.
             */
            void initialize (ParameterHandler &param);
            //@}
            
            ///@name Instance Delivery (Replica creator)
            //@{
            /**
             * This member function is used to create an object of type gas diffusion layer
             * 
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::Layer::CatalystLayer<dim> > create_replica (const std::string &cl_section_name)
            {
                return boost::shared_ptr<FuelCellShop::Layer::CatalystLayer<dim> > (new FuelCellShop::Layer::HomogeneousCL<dim> (cl_section_name));
            }   
            //@}           
            
            ///@name Instance Delivery (Prototype)
            //@{
            /**
             *
             */            
            static HomogeneousCL<dim> const* PROTOTYPE;
            //@}
            
            /** This routine is not used for this layer */   
            void set_cell_id(const unsigned int& ) {}
        };        
    }
}

#endif