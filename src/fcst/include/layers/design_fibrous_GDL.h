//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: design_fibrous_GDL.h
//    - Description: Class used to represent a fibrous GDL where effective properties are computed based on the porosity etc.
//    - Developers: M. Secanell and Madhur Bhaiya
//    - $Id: design_fibrous_GDL.h 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__DESIGN_FIBROUS_GDL_H
#define _FUELCELLSHOP__DESIGN_FIBROUS_GDL_H

// FCST classes
#include <layers/gas_diffusion_layer.h>

using namespace dealii;

namespace FuelCellShop
{
    namespace Layer
    {
        /**
         * This class defines a GDL made of fibres. Due to this, the conductivity in the X and Y direction is 
         * extremely anisotropic since the contact resistance is very high. Given this situation, I decided to
         * use a different value for the bulk resistance in the X and Y direction.
         *
         * @author M. Secanell and M. Bhaiya
         * 
         */
        template <int dim>
        class DesignFibrousGDL : 
        public  GasDiffusionLayer<dim>
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
            
            ///@name Destructor, and initalization
            //@{    
            /** 
             * Replica Constructors
             * 
             * \warning For internal use only.
             * 
             * Constructor used only to create a prototype. Do not use
             * in general since this will not include the name of the section
             * in the parameter file you need. 
             * 
             */
            DesignFibrousGDL();
            
                                 
            /** Member function used to read in data and initialize the necessary data
             * to compute the coefficients. */
            virtual void initialize (ParameterHandler &param);
            //@}
            
             ///@name Effective properties
            //@{               
            
            /**
             * Compute the effective property in the pores of the GDL. This is used for example to
             * compute effective diffusivity of gases. The method takes in bulk diffusion coefficient [\p m^2/s] and liquid water saturation as the first and second 
             * argument respectively. This routine is used in the isotropic case.
             */
            virtual void effective_gas_diffusivity(const double&, const double&, double&) const;
            /**
             * Compute the effective property in the pores of the GDL. This is used for example to
             * compute effective diffusivity of gases. The method takes in bulk diffusion coefficient [\p m^2/s] and liquid water saturation as the first and second 
             * argument respectively. This routine can be used either in the isotropic or anisotripic case
             */
            virtual void effective_gas_diffusivity(const double&, const double&, Tensor<2,dim>&) const;
            
            /**
             * Return the effective diffusivity [\p m^2/s] for nonisothermal with/without two-phase case in the GDL. It takes bulk diffusivity, computed using
             * #compute_gas_diffusion method and transforms it into an effective property,
             * taking into account the porosity, saturation and GDL structure (Anisotropic case),
             * at all quadrature points of the cell.
             * 
             * \note: For two-phase case, set_saturation should be called before using this method, otherwise this method assumes saturation value to be zero.
             */
            virtual void effective_gas_diffusivity(std::vector< Tensor<2,dim> >& ) const;
            /**
             * Return the derivative of effective diffusivity w.r.t solution variables/design parameters
             * for nonisothermal with/without two-phase case in the GDL. It transforms bulk diffusion properties computed using
             * #compute_gas_diffusion method and transforms it into an effective property,
             * taking into account the porosity, saturation and GDL structure (Anisotropic case),
             * at all quadrature points of the cell.
             * 
             * \note: For two-phase case, set_saturation should be called before using this method, otherwise this method assumes saturation value to be zero.
             */
            virtual void derivative_effective_gas_diffusivity(std::map< VariableNames, std::vector< Tensor<2,dim> > >&) const;
            
            /**
             * Compute effective diffusivity of the media.
             * 
             */
            virtual void effective_gas_diffusivity(Table< 2, Tensor<2,dim> > &D_eff) const;
            
            /**
             * Compute the electrical conductivity of the media.
             */
            virtual void effective_electron_conductivity(double&) const;
            
            /**
             * Compute the electrical conductivity of the media.
             */
            virtual void effective_electron_conductivity(Tensor<2,dim>&) const;
            
            /**
             * Compute effective thermal conductivity of the media as a function of Temperature. (Anisotropic values) or return anisotropic Tensor
             * of GIVEN values set using parameter file.
             * Ref: N. Zamel, E. Litovsky, X. Li, and J. Kleiman. Measurement of the through-plane thermal conductivity of carbon paper diffusion media
             * for the temperature range from -50 to 120 c. International Journal of Hydrogen Energy, 36(19):12618-12625, 2011.
             * It return vector of Tensor for thermal conductivity (anisotropic).
             * NOTE: Before using Zamel method, temperature values need to be set in the layer.
             */
            virtual void effective_thermal_conductivity(std::vector< Tensor<2,dim> >&) const;
            /**
             * Compute derivative of thermal conductivity with respect to variables set by set_derivative_flags.
             */
            virtual void derivative_effective_thermal_conductivity(std::vector< Tensor<2,dim> >&) const;
            
            /**
             * Compute the anisotropic GDL liquid permeability \f$ \left[ cm^2 \right] \f$, at all quadrature points in the cell.
             */
            virtual void liquid_permeablity(std::vector< Tensor<2,dim> >& ) const;
            /**
             * Compute the derivative of the anisotropic liquid permeability in the GDL
             * with respect to either the solution or design parameters, at all quadrature points in the cell. The parameters with respect to
             * which the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_liquid_permeablity(std::map< VariableNames, std::vector< Tensor<2,dim> > >& ) const;
            
            /**
             * Compute the anisotropic CL liquid permeability \f$ \left[ cm^2 \right] \f$, at all quadrature points in the cell.
             */
            
            virtual void saturated_liquid_permeablity_PSD(double&) const;
            virtual void relative_liquid_permeability_PSD(std::vector< Tensor<2,dim> >&) const;
            virtual void derivative_relative_liquid_permeablity_PSD(std::vector< double >&) const;
            
            virtual void derivative_relative_liquid_permeablity_PSD(std::map< VariableNames, std::vector< Tensor<2,dim> > >&) const;
            
            /**
             * Compute \f$ p_c \quad \left[ dyne \cdot cm^{-2}\right] \f$, at all quadrature points in the cell.
             */
            virtual void pcapillary(std::vector<double>&) const;
            virtual void saturation_from_capillary_equation(std::vector<double>&) const;
            
            virtual void derivative_saturation_from_capillary_equation_PSD(std::vector<double>&) const;
            /**
             * Compute \f$ \frac{\partial p_c}{\partial s} \quad \left[ dyne \cdot cm^{-2}\right] \f$, at all quadrature points in the GDL.
             */
            virtual void dpcapillary_dsat(std::vector<double> &) const;
            /**
             * Compute the derivative of \f$ \frac{\partial p_c}{\partial s} \quad \left[ dyne \cdot cm^{-2}\right] \f$ in the GDL, with 
             * respect to either the solution or design parameters, at all quadrature points in the cell. The parameters with respect to which 
             * the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_dpcapillary_dsat(std::map< VariableNames, std::vector<double> > &) const;
            
            /**
             * Compute the liquid-gas interfacial surface area per unit volume, \f$ a_{lv} ~\left[ \frac{cm^2}{cm^3} \right] \f$, at all 
             * quadrature points in the GDL.
             */
            virtual void interfacial_surface_area(std::vector<double>&) const;
            /**
             * Compute the derivative of the liquid-gas interfacial surface area per unit volume, with respect to either the 
             * solution variables or design parameters, at all quadrature points in the GDL. The parameters with respect to which 
             * the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags().
             */
            virtual void derivative_interfacial_surface_area(std::map< VariableNames, std::vector<double> >&) const;
            
            /**
             * Compute the liquid-gas interfacial surface area per unit volume, \f$ a_{lv} ~\left[ \frac{cm^2}{cm^3} \right] \f$, at all 
             * quadrature points in the CL.
             */
            virtual void interfacial_surface_area_PSD(std::vector<double>&) const;
            /**
             * Compute the derivative of the liquid-gas interfacial surface area per unit volume, with respect to either the 
             * solution variables or design parameters, at all quadrature points in the CL. The parameters with respect to which 
             * the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags().
             */
            virtual void derivative_interfacial_surface_area_PSD(std::vector<double>&) const;
            virtual void derivative_interfacial_surface_area_PSD(std::map< VariableNames, std::vector<double> >&) const;
            
            //@}

        private:
            ///@name Constructors, destructor, and initalization
            /** 
             * DesignFibrousGDL constructor
             * 
             * This member function is used in create_replica to create a class.
             */
            DesignFibrousGDL(const std::string& name);
                        
            //@{    
            /** Declare parameters for a parameter file. The parameters that need to be declared are 
             * The parameters that need to be declared are
             * - Porosity (default : 0.6) Represents *the porosity in the GDL
             * - Method effective transport properties in pores (default: Bruggemann) Other options "Given|Bruggemann|Percolation|Mezedur"
             * - Method effective transport properties in solid (default: Bruggemann) Other options "Given|Bruggemann|Percolation|Mezedur"
             *
             * - Electrical conductivity X
             * - Electrical conductivity Y
             * - Electrical conductivity Z
             * 
             * An example input file would look as follows:
             * 
             * @code
             * subsection Cathode gas diffusion layer
             *    set Material id = 2
             *    set Gas diffusion layer type = DesignFibrousGDL
             *    subsection DesignFibrousGDL
             *      ######### Composition: #########
             *      set Porosity = 0.6
             *      ######### Gas transport #########
             *      ## Anisotropy
             *      set Anisotropic transport = true                  # (default) false
             *      set Method effective transport properties in pores = Tomadakis    # (default) Bruggemann | Given | Percolation | Tomadakis | Mezedur
             *      set Method effective transport properties in solid = Percolation  # (default) Bruggemann | Given | Percolation
             *      ## XX
             *      set Porosity threshold X = 0.11               # (default) 0.12 | 0.118 (Peter's Thesis) | 0.11
             *      set Porosity network constant X = 0.785           # (default) 2.0  |0.785 (Peter's Thesis)  [Page 69]
             *      set Porosity gamma network constant X = 0.0           # (default) 0.0 | 
             *      #
             *      set Electrical conductivity X [S/cm] = 16.03
             *      set Solid network threshold X = 0.0               # (default) 0.12 | 
             *      set Solid network constant X = 1.5                # (default) 2.0 |
             *      ## YY
             *      set Porosity threshold Y = 0.11               # (default) 0.12 | 0.118 (Peter's Thesis) | 0.11
             *      set Porosity network constant Y = 0.521           # (default) 2.0 |
             *      set Porosity gamma network constant Y = 0.0           # (default) 0.0 |
             *      #
             *      set Electrical conductivity Y [S/cm] = 272.78
             *      set Solid network threshold Y = 0.0               # (default) 0.12 |
             *      set Solid network constant Y = 1.0                # (default) 2.0 |
             *    end
             *  end
             * @endcode
             * 
             */
            virtual void declare_parameters (const std::string& name, 
                                             ParameterHandler &param) const;

            //@}                     
            ///@name Instance delivery function
            //@{
            /**
             * This member function is used to create an object of type gas diffusion layer
             * 
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::Layer::GasDiffusionLayer<dim> > create_replica (const std::string &name)
            {
                return boost::shared_ptr<FuelCellShop::Layer::GasDiffusionLayer<dim> > (new FuelCellShop::Layer::DesignFibrousGDL<dim> (name));
            }
            //@}           
            
            ///@name Internal variables
            //@{
            //--------------------------            
            /** Porosity of the GDL */
            double porosity;
            /** Volume fraction of solid phase, i.e. solid */
            double solid_phase;
            /** Method used to compute the effective properties in the pores */
            std::string method_eff_property_pores;
            /** Method used to compute the effective properties in the solid */
            std::string method_eff_property_solid;
            /** Method used to compute the effective thermal properties*/
            std::string method_eff_property_thermal;
            
            /** Method used to compute the relative liquid permeability. */
            std::string method_rel_liquid_permeability;
            /** Irreducible liquid water saturation value in the GDL. */
            double s_irr;
            
            /** Boolean flag corresponding to anisotropy. */
            bool anisotropy;
            
            /** Porosity of the GDL threshold */
            std::vector<double> porosity_th;
            /** Network constant */
            std::vector<double> porosity_mu;
            /** Network constant gamma */
            std::vector<double> porosity_gamma;
            
            /** Solid network of the GDL threshold */
            std::vector<double> solid_th;
            /** Solid network constant */
            std::vector<double> solid_mu;
            
            /** "Bulk" electrical conductivity [S/cm] of the solid material. */
            std::vector<double> sigma_e;
            
            /** "Given" thermal conductivity [W/(cm-K)] of the layer. */
            std::vector<double> k_thermal;
            
            /** Absolute permeability [cm^2] of the layer. */
            std::vector<double> abs_permeability;
            
            /** Method used to compute capillary pressure as a function of saturation. */
            std::string method_capillary_function;
            /** GDL Compaction pressure, \f$ C ~\left[ MPa \right] \f$. */
            double compaction_pressure;
            /** PTFE loading (% wt) in the GDL.  */
            double PTFE_loading;
            /** 
             * Factor calculated based on Kumbur et al (2007), to be used in capillary pressure computation, given as:
             * \f$ 2^{0.4 C} ~ \sqrt{ \frac{\epsilon_c}{\kappa} } \f$
             */
            double kumbur_factor;
            
            //@}
            
            ///@name Instance delivery object
            //@{
            /**
             * Prototype declaration
             */             
            static DesignFibrousGDL<dim> const* PROTOTYPE;
            //@}
        };
    }
}

#endif
