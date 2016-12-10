//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: design_MPL.h
//    - Description: Class used to represent a design MPL where effective properties are computed based on the porosity etc.
//    - Developers: Madhur Bhaiya
//    - $Id: design_MPL.h 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__DESIGN_MPL_H
#define _FUELCELLSHOP__DESIGN_MPL_H

// FCST classes
#include <layers/micro_porous_layer.h>

using namespace dealii;

namespace FuelCellShop
{
    namespace Layer
    {
        /**
         * This class defines an MPL where effective transport properties are computed
         * using macro-homogeneous correlations to estimate the effective properties
         * of the media. 
         *
         * Currently the following macro-homogeneous correlations are implemented for the gas phase:
         * - Bruggemann: Most commonly used in the literature. The effective gas diffusivity only depends on 
         * the porosity of the media. The equation implemented is:
         * \f[ 
         * D_{eff} = D*\epsilon^1.5
         * \f] 
         * where \f$ \epsilon \f$ is the porosity and is specified in the input file in the following section:
         * @code
         * subsection Cathode microporous layer        <- name of the subsection identified in the FuelCellShop::Layer::MicroPorousLayer< dim >::declare_all_MicroPorousLayer_parameters and FuelCellShop::Layer::MicroPorousLayer< dim >::create_MicroPorousLayer
         *    set Material id = 3
         *    set Micro porous layer type = DesignMPL <- concrete_name for this layer
         *    subsection DesignMPL
         *      set Porosity = 0.4                            # Value used for epsilon
         *    end
         * end
         * @endcode
         * - Percolation: First proposed by Eikerling, M., Kornyshev, A.A. Modelling the performance of the cathode catalyst layer of 
         * polymer electrolyte fuel cells (1998) Journal of Electroanalytical Chemistry 453 (1-2) , pp. 89-106, 181 and later modified and used
         * in M. Secanell, K. Karan, A. Suleman and N. Djilali, "Multi-Variable Optimization of PEMFC Cathodes using an 
         * Agglomerate Model ", Electrochimica Acta, 52(22):6318-6337, 2007. and P. Dobson, C. Lei, T. Navessin, M. Secanell, "Characterization of the PEM Fuel Cell 
         * Catalyst Layer Microstructure by Nonlinear Least-Squares Parameter Estimation", Journal of the Electrochemical Society, 159(5), B1-B10, 2012.
         * Effective propeties are computed based on the porosity, a porosity threshold and a coefficient that is related to the dimensionality of the lattice 
         * structure. The equation used is:
        * \f[ 
         * D_{eff} = D*(\epsilon - \epsilon_{th})/(1 - \epsilon_{th})^\mu;
         * \f] 
         * where \f$ \epsilon \f$ is the porosity, \f$ \epsilon_{th} \f$ is the threshold porosity below which
         * the gas transport is zero due to a lack of connectiveness of the pores and \f$ \mu \f$ is a parameter
         * between 2 and 3. These parameters might be specified for each direction in the media in the input file.
         * The sections where they are specified are:
         * 
         * @code
         * subsection Fuel cell data
         * (...)
         * subsection Cathode microporous layer
         *    set Material id = 3
         *    set Micro porous layer type = DesignMPL
         *    subsection DesignMPL
         *      set Porosity = 0.4                            # From experimental data (manufacturer's data) on Sigracel 24BC
         *      set Anisotropic transport = false
         *      ######### Pore network #########
         *      set Method effective transport properties in pores = Percolation          <- Specify if you would like to use Percolation or Bruggemann
         *      set Porosity threshold = 0.118                                            <- Specify threshold porosity value
         *      set Porosity network constant = 2.0                                       <- Specify mu value/
         *      ######### Solid network #########
         *      set Method effective transport properties in solid phase = Percolation
         *      set Electric conductivity = 88.84             # From experimental data (manufacturer's data) on Sigracel 24BC
         *      set Solid network threshold = 0.118
         *      set Solid network constant = 2.0
         *    end
         *  end
         *  (...)
         * @endcode
         * If you would like to define an anisotripic media, please set Anisotropic transport = true and then
         * all the properties will have an X, Y and Z at the end of their name, e.g. set Porosity threshold X = 0.118  
         * 
         * For the solid phase properties the methods above can also be selected. The properties are set using the same rules.
         * 
         * \todo This class should be divided into two classes DummyMPL and DesignMPL. DummyMPL is a
         * class where all effective properties are directly specified in the input file. DesignMPL is the
         * class used to compute the effective properties based on effective porous media approximations.
         * 
         * 
         * \author M. Secanell, 2013
         * 
         */
        template <int dim>
        class DesignMPL : 
        public  MicroPorousLayer<dim>
        {
        public:
            ///@name Declaration
            //@{    
            /**
             * Concrete name used for objects of this class. This name is used when
             * setting up the subsection where the data is stored in the input file.
             * 
             * The data will be store under
             * \code
             * subsection name_specified_in_constructor
             *    set Material id = 2
             *    set Microporous layer type = DummyMPL # <-here I select the type of object of type GasDiffusionLayer
             *    subsection DummyMPL # <- this is the concrete_name for this class
             *       set all info relevant to this object
             *    end
             * end
             * \endcode
             */
            static const std::string concrete_name;
            
            //@}         
            ///@name Constructors, destructor, and initalization
            //@{    
            /**
             * Constructor
             */
            DesignMPL(std::string name);
            
            /**
             * Replica Constructor
             */
            DesignMPL();
            
            /**
             * Destructor
             */
            ~DesignMPL()
            {}
                       
            /**
            * Declare parameters for a parameter file. 
            * 
            * \deprecated Use FuelCellShop::Layer::MicroPorousLayer<dim>::declare_all_MicroPorousLayer_parameters
            */
            void declare_parameters (ParameterHandler &param) const
            {
                FuelCellShop::Layer::DesignMPL<dim>::declare_parameters(this->name,param);
            }
            /**
             * Member function used to set new parameters values in the optimization loop
             */
            void set_parameters (const std::vector<std::string>& name_dvar,const std::vector<double>& value_dvar,ParameterHandler &param);
            
            
            /**
             * Member function used to read in data and initialize the necessary data
             * to compute the coefficients.
             */
            void initialize (ParameterHandler &param);
            //@}
            
            ///@name Accessors and info
            //@{
                
            
            //@}
            
            ///@name Effective property calculators
            //@{

            /**
             * Compute the effective property in the pores of the MPL. This is used for example to
             * compute effective diffusivity of gases. The method takes in bulk diffusion coefficient [\p m^2/s] and \b non-negative liquid water saturation as the 
             * first and second argument respectively. This routine is used in the isotropic case.
             */
            virtual void effective_gas_diffusivity(const double& property,
                                                   const double& saturation,
                                                   double& effective_property) const;                        
            /**
             * Compute the effective property in the pores of the MPL. This is used for example to
             * compute effective diffusivity of gases. The method takes in bulk diffusion coefficient [\p m^2/s] and \b non-negative liquid water saturation as the 
             * first and second argument respectively. This routine can be used either in the isotropic or anisotripic case
             */
            virtual void effective_gas_diffusivity(const double& property,
                                                   const double& saturation,
                                                   Tensor<2,dim>& effective_property) const;
                                           
            /**
             * Return the effective diffusivity [\p m^2/s] for nonisothermal with/without two-phase case in the MPL. It takes bulk diffusivity, computed using
             * #compute_gas_diffusion method and transforms it into an effective property,
             * taking into account the porosity, saturation and MPL structure (Anisotropic case),
             * at all quadrature points of the cell.
             * 
             * \note: For two-phase case, set_saturation should be called before using this method, otherwise this method assumes saturation value to be zero.
             */
            virtual void effective_gas_diffusivity(std::vector< Tensor<2,dim> >& ) const;
            /**
             * Return the derivative of effective diffusivity w.r.t solution variables/design parameters
             * for nonisothermal with/without two-phase case in the MPL. It transforms bulk diffusion properties computed using
             * #compute_gas_diffusion method and transforms it into an effective property,
             * taking into account the porosity, saturation and MPL structure (Anisotropic case),
             * at all quadrature points of the cell.
             * 
             * \note: For two-phase case, set_saturation should be called before using this method, otherwise this method assumes saturation value to be zero.
             */
            virtual void derivative_effective_gas_diffusivity(std::map< VariableNames, std::vector< Tensor<2,dim> > >& ) const;
            
            
            /**
             * Compute the effective property in the pores. This is used to compute
             * effective diffusivity of gases.
             * This routine can be used either in the isotropic or anisotropic cases.
             * Bulk diffusion coefficients or their derivatives are obtained from Mixure::BinaryDiffusion classes
             * inside this method.
             * 
             * \note The routine FuelCellShop::Layer::PorousLayer< dim >::set_gases_and_compute (std::vector< FuelCellShop::Material::PureGas * > &gases, double pressure, double temperature)
             * (in the parent class) should have been called prior to using this class. This method is to be used only for a single-phase, isothermal application.
             */                            
            virtual void effective_gas_diffusivity(Table< 2, Tensor< 2, dim > > &D_eff) const;
            
            
            /**
             * Compute the effective conductivity. To compute the effective properties it will
             * use the method specified in set_method_effective_transport_property_solid . 
             * This function is used for the isotropic MPL
             */
            virtual void effective_electron_conductivity(double& ) const;
            
            /**
             * Compute the effective conductivity. To compute the effective properties it will
             * use the method specified in set_method_effective_transport_property_solid . 
             * This function is used for the anisotropic MPL.
             */
            virtual void effective_electron_conductivity(Tensor<2,dim>& ) const;
            /**
             * Compute the effective thermal conductivity. It's computed based on the 
             * "Method effective thermal conductivity" specified in the parameter file.
             * This function is used for the isotropic MPL
             */
            virtual void effective_thermal_conductivity(double& ) const;
            
            /**
             * Compute the effective thermal conductivity. It's computed based on the 
             * "Method effective thermal conductivity" specified in the parameter file.
             * This function is used for the anisotropic MPL.
             */
            virtual void effective_thermal_conductivity(Tensor<2,dim>& ) const;
            
            /**
             * Compute the effective thermal conductivity at all quadrature points in the cell.
             * It's computed based on the "Method effective thermal conductivity" specified in
             * the parameter file. This function is used for the anisotropic MPL.
             */
            virtual void effective_thermal_conductivity(std::vector< Tensor<2,dim> >& ) const;
            
            
            /**
             * Compute the anisotropic MPL liquid permeability \f$ \left[ cm^2 \right] \f$, at all quadrature points in the cell.
             */
            virtual void liquid_permeablity(std::vector< Tensor<2,dim> >& ) const;
            /**
             * Compute the derivative of the anisotropic liquid permeability in the MPL
             * with respect to either the solution or design parameters, at all quadrature points in the cell. The parameters with respect to
             * which the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_liquid_permeablity(std::map< VariableNames, std::vector< Tensor<2,dim> > >& ) const;
            /**
             * Compute the anisotropic CL liquid permeability \f$ \left[ cm^2 \right] \f$, at all quadrature points in the cell.
             */
            
            virtual void saturated_liquid_permeablity_PSD(double&) const;
            virtual void relative_liquid_permeability_PSD(std::vector< Tensor<2,dim> >&) const;
            virtual void derivative_relative_liquid_permeablity_PSD(std::vector<double>&) const;
            virtual void derivative_relative_liquid_permeablity_PSD(std::map< VariableNames, std::vector< Tensor<2,dim> > >&) const;
            
            
            /**
             * Compute \f$ p_c \quad \left[ dyne \cdot cm^{-2}\right] \f$, at all quadrature points in the cell.
             */
            virtual void pcapillary(std::vector<double>&) const;
            virtual void saturation_from_capillary_equation(std::vector<double>&) const;
            
            virtual void derivative_saturation_from_capillary_equation_PSD(std::vector<double>&) const;
            /**
             * Compute \f$ \frac{\partial p_c}{\partial s} \quad \left[ dyne \cdot cm^{-2}\right] \f$, at all quadrature points in the MPL.
             */
            virtual void dpcapillary_dsat(std::vector<double> &) const;
            /**
             * Compute the derivative of \f$ \frac{\partial p_c}{\partial s} \quad \left[ dyne \cdot cm^{-2}\right] \f$ in the MPL, with 
             * respect to either the solution or design parameters, at all quadrature points in the cell. The parameters with respect to which 
             * the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_dpcapillary_dsat(std::map< VariableNames, std::vector<double> > &) const;
            
            /**
             * Compute the liquid-gas interfacial surface area per unit volume, \f$ a_{lv} ~\left[ \frac{cm^2}{cm^3} \right] \f$, at all 
             * quadrature points in the MPL.
             */
            virtual void interfacial_surface_area(std::vector<double>&) const;
            /**
             * Compute the derivative of the liquid-gas interfacial surface area per unit volume, with respect to either the 
             * solution variables or design parameters, at all quadrature points in the MPL. The parameters with respect to which 
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
            
            
             /**                                                                     
             * Compute the effective property of a property that is defined by the
             * network of fibres. For example this could be used to compute the effective
             * electron conductivity or heat conduction.
             * NOTE: Isotropic case
             */
            virtual void effective_transport_property_solid(const double& property, 
                                                    double& effective_property) const;
                                                    
            /**
             * Compute the effective property of a property that is defined by the
             * network of fibres. For example this could be used to compute the effective
             * electron conductivity or heat conduction.
             * Note: Anisotropic case.
             */
            virtual void effective_transport_property_solid(const Tensor<2,dim>& property, 
                                                    Tensor<2,dim>& effective_property) const;
            
            //@}

        protected:
            ///@name Declaration
          //@{    
            /**
             * Declare parameters for a parameter file. The parameters that need to be declared are
             * - Porosity (default : 0.3) Represents the porosity in t*he GDL
             * - Method effective transport properties in pores (default: Bruggemann) Other options "Given|Bruggemann|Percolation|Mezedur"
             * - Method effective transport properties in fibres (default: Bruggemann) Other options "Given|Bruggemann|Percolation|Mezedur"
             */
            void declare_parameters (const std::string& mpl_section_name, 
                                     ParameterHandler &param) const;
                                     
                                     
            /**
             * Member function used to set new parameters values in the optimization loop
             * 
             */
            void set_parameters (const std::vector<std::string>& name_dvar,
                                 const std::vector<double>& value_dvar, 
                                 const std::string& name, 
                                 ParameterHandler &param) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
                    //@}
                                     
            ///@name Instance Delivery
            //@{
            virtual boost::shared_ptr<FuelCellShop::Layer::MicroPorousLayer<dim> > create_replica (const std::string &name)
            {
                return boost::shared_ptr<FuelCellShop::Layer::MicroPorousLayer<dim> > (new FuelCellShop::Layer::DesignMPL<dim> (name));
            }    
            
            /**
             * Prototype declaration
             */             
            static DesignMPL<dim> const* PROTOTYPE;
            
            
            //@}
            
            ///@name Internal variables
            //@{
            /** General properties */
            //--------------------------
            /** Anisotropy ? */
            bool anisotropy;
            /** Porosity of the GDL */
            double porosity;
            /** Volume fraction of solid phase, i.e. fibres */
            double solid_phase;
            /** Method used to compute the effective properties in the pores */
            std::string method_eff_property_pores;
            /** Method used to compute the effective properties in the solid phase */
            std::string method_eff_property_fibres;
            /** Method used to compute effective thermal conductivity. */
            std::string method_eff_thermal_conductivity;
            /** Electrical conductivity from the input file */
            double electrical_conductivity;
            /** Thermal conductivity from the input file */
            double thermal_conductivity;
            /** Electrical conductivity from the input file in the anisotripic case */
            Tensor<2,dim> matrix_electrical_conductivity;
            /** Thermal conductivity from the input file in the anisotripic case  */
            Tensor<2,dim> matrix_thermal_conductivity;
            
            /** Method used to compute the relative liquid permeability. */
            std::string method_rel_liquid_permeability;
            /** Irreducible liquid water saturation value in the MPL. */
            double s_irr;
            
            //----------------------
            /** Anisotropic properties */
            /** Porosity of the GDL threshold */
            std::vector<double> porosity_th;
            /** Network constant */
            std::vector<double> porosity_mu;
            /** Network constant gamma */
            std::vector<double> porosity_gamma;
            /** Oxygen Diffusion coefficient */
            std::vector<double> D_O2;
            /** Water vapour diffusion coefficient */
            std::vector<double> D_wv;
            /** Solid (electron conductive) network of the MPL threshold */
            std::vector<double> fibre_th;
            /** Solid (electron conductive) network constant */
            std::vector<double> fibre_mu;
            
            /** Absolute permeability [cm^2] of the layer. */
            std::vector<double> abs_permeability;
            
            /** Method used to compute capillary pressure as a function of saturation. */
            std::string method_capillary_function;
            /** MPL Compaction pressure, \f$ C ~\left[ MPa \right] \f$. */
            double compaction_pressure;
            /** PTFE loading (% wt) in the MPL.  */
            double PTFE_loading;
            /** 
             * Factor calculated based on Kumbur et al (2007), to be used in capillary pressure computation, given as:
             * \f$ 2^{0.4 C} ~ \sqrt{ \frac{\epsilon_c}{\kappa} } \f$
             */
            double kumbur_factor;
            
            //@}
        };
    }
}  // FuelCellShop

#endif // _FUELCELLSHOP__DESIGN_MPL_H
