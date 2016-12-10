// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: dummy_CL.h
// - Description: This class characterizes a conventional catalyst layer
//                and defines constant effective properties
// - Developers: Marc Secanell and Madhur Bhaiya
//
// ----------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__DUMMY_CL__H
#define _FUELCELLSHOP__DUMMY_CL__H

// Include FCST classes
#include <utils/fcst_units.h>
#include <layers/catalyst_layer.h>



namespace FuelCellShop
{
    namespace Layer
    {
        /**
         * This class characterizes a macro-homogeneous catalyst layer and should
         * be used in the case of constant effective properties, \em viz., \p effective_proton_conductivity, \p effective_gas_diffusivity,
         * \p effective_electron_conductivity and \p effective_thermal_conductivity.
         *
         * \note These constant effective properties should be set using the parameter file.
         *
         * @author M. Secanell and M. Bhaiya, 2011-13
         *
         */
        template <int dim>
        class DummyCL :
        public CatalystLayer<dim>
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
             *    set Gas diffusion layer type = DummyGDL # <-here I select the type of object of type GasDiffusionLayer
             *    subsection DummyGDL # <- this is the concrete_name for this class
             *       set all info relevant to this object
             *    end
             * end
             * \endcode
             */
            static const std::string concrete_name;


            ///@name Constructors, destructor, and initalization
            //@{
            /**
             * Prototype Constructor
             *
             * \warning For internal use only
             */
            DummyCL();

            /**
            * Destructor
            */
            ~DummyCL();

            /**
            * Declare all necessary parameters in order to compute the coefficients
            *
            * \deprecated Use declare_all_CatalystLayer_parameters
            */
            virtual void declare_parameters (ParameterHandler &param) const
            {
               declare_parameters(this->name, param);
            };

            /**
            * Member function used to read in data and initialize the necessary data
            * to compute the coefficients.
            */
            virtual void initialize (ParameterHandler &param);
            //@}

            ///@name Accessors and info
            //@{
            /**
            * Compute the effective diffusivty in the CL. This routine takes the
            * gas diffusivity from FuelCellShop::BinaryDiffusion and transforms
            * it into an effective property taking into account the porosity and
            * structure of the CL
            */
            virtual void effective_gas_diffusivity(Table< 2, Tensor< 2, dim > >&) const;

            void effective_gas_diffusivity(std::vector< Tensor<2,dim> >& prop_eff_vec) const;
            
            void derivative_effective_gas_diffusivity(std::map< VariableNames, std::vector< Tensor<2,dim> > >&dprop_eff)  const;
            
            /**
            * Compute the effective electron conductivity in the CL
            */
            virtual void effective_electron_conductivity(double& ) const;

            /**
            * Compute the effective electron conductivity in the CL
            */
            virtual void effective_electron_conductivity(Tensor<2,dim>& ) const;

            /**
            * Compute the effective proton conductivity in the CL.
            * \note This class is awesome for testing the capabilites of the model since it could
            * be implemented as dependant on the solution or not...
            * In this first iteration, I've implemented it as constant
            */
            virtual void effective_proton_conductivity(double& ) const;

            /**
             * Compute the effective proton conductivity in the CL, at quadrature points in the cell.
             * \note Input vector should be initialized (size set to number of quadrature points), before passing over to this function.
             */
            virtual void effective_proton_conductivity(std::vector<double>& ) const;

            /**
             * Compute the derivative of the effective proton conductivity in the CL
             * with respect to either the solution or design parameters. The parameters with respect to
             * which the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_effective_proton_conductivity(std::map< VariableNames, std::vector<double> >& ) const;


            /** Get the active area of platinum per unit volume of CL */
            inline double get_active_area_Pt() const
            {
                return Av.at(this->local_material_id());
            }
            /** This routine is not used for this layer */
            void set_cell_id(const unsigned int& ){}

            /**
            * This member function will use a FuelCellShop::Kinetics class in order to compute the current density
            * production in the CL
            */
            virtual void current_density(std::vector<double>&);
            /**
            * This member function will use a FuelCellShop::Kinetics class in order to compute the derivative of the current density
            * with respect to the variables setup using #set_derivative_flags method.
            */
            virtual void derivative_current_density(std::map< VariableNames, std::vector<double> >& );
            /**
             * This member function computes the current density production in the CL. First argument is <b>current density</b>, and second is <b>effectiveness</b>, at
             * all quadrature points in the cell. Since this is a dummy layer, effectiveness is filled as \b 1.0
             */
            virtual void current_density(std::vector<double>& current, std::vector<double>& effectiveness)
            {
                current_density(current);
                effectiveness.assign(current.size(), 1.0);
            }

            //@}
        private:
            ///@name Constructors
            //@{
            /**
            * Constructor
            */
            DummyCL(const std::string& name);

            /**
             * Declare parameters for a parameter file.
             *
             */
            void declare_parameters (const std::string& name,
                                     ParameterHandler &param) const
            {

                FuelCellShop::Layer::CatalystLayer<dim>::declare_parameters(name, param);

                param.enter_subsection("Fuel cell data");
                {
                    param.enter_subsection(name);
                    {
                        param.enter_subsection(concrete_name); //-- Transport for the anisotropic case:
                        {
                            param.declare_entry ("Oxygen diffusion coefficient, [cm^2/s]",
                                                 "4:0.02514", //1atm, 353K
                                                 Patterns::Map( Patterns::Integer(0,255), Patterns::Double(0) ),
                                                 "Oxygen diffusion coefficient given by experiment");
                            param.declare_entry ("Water vapour diffusion coefficient, [cm^2/s]",
                                                 "4:0.29646",
                                                 Patterns::Map( Patterns::Integer(0,255), Patterns::Double(0) ),
                                                 "Water vapour diffusion coefficient given by experiment");
                            param.declare_entry ("Electrical conductivity, [S/cm]",
                                                 "4:40", // [S/cm]
                                                 Patterns::Map( Patterns::Integer(0,255), Patterns::Double(0) ),
                                                 "Effective cond. if given is used, otherwise conductivity of the raw material. Units [S/cm]");
                            param.declare_entry ("Protonic conductivity, [S/cm]",
                                                 "4:40", // [S/cm]
                                                 Patterns::Map( Patterns::Integer(0,255), Patterns::Double(0) ),
                                                 "Effective cond. if given is used, otherwise conductivity of the raw material. Units [S/cm]");
                            param.declare_entry ("Active area [cm^2/cm^3]",
                                                 "4:2.0e5",
                                                 Patterns::Map( Patterns::Integer(0,255), Patterns::Double(0)));
                        }
                        param.leave_subsection();

                    }
                    param.leave_subsection();
                }
                param.leave_subsection();

            }


            //@}

            ///@name Instance Delivery
            //@{
            /**
             * This member function is used to create an object of type gas diffusion layer
             *
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::Layer::CatalystLayer<dim> > create_replica (const std::string &name)
            {
                return boost::shared_ptr<FuelCellShop::Layer::CatalystLayer<dim> > (new FuelCellShop::Layer::DummyCL<dim> (name));
            }
            /**
             *
             */
            static DummyCL<dim> const* PROTOTYPE;
            //@}

            ///@name Internal variables
            //@{
            /** Oxygen diffusion coefficient */
            std::map< unsigned int, double > D_O2;
            /** Water vapour diffusion coefficient */
            std::map< unsigned int, double > D_wv;
            /** Solid network conductivity */
            std::map< unsigned int, double > sigma_e;
            /** Membrane phase conductivity */
            std::map< unsigned int, double > sigma_m;
            /** Active area of catalyst per unit volume of catalyst layer */
            std::map< unsigned int, double > Av;
            //@}
        };
    }
}

#endif