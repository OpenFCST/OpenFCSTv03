//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: conventional_cl.h
//    - Description: Class characterizing the conventional catalyst layer and methods for computing effective properties. It also provides interface to various material classes used in catalyst layer.
//    - Developers: Peter Dobson (2011) and Madhur Bhaiya (2013)
//    - Id: $Id: conventional_CL.h 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__CONVENTIONAL_CL__H
#define _FUELCELLSHOP__CONVENTIONAL_CL__H

// Include deal.II classes
#include <deal.II/base/types.h>

// Include FCST classes
#include <layers/catalyst_layer.h>
#include <grid/geometries.h>
#include <boost/shared_ptr.hpp>

//Include STL
#include <algorithm>
#include <math.h> 
#include <memory>

using namespace dealii;

namespace FuelCellShop
{
    namespace Layer
    {
        /**
         * This class characterizes a catalyst layer and uses this information 
         * to compute effective transport properties and interfacial areas for phase
         * change or electrochemical reactions.
         * 
         * This class implements a macrohomogeneous homogeneous or graded catalyst layer.
         *
         */
        template <int dim>
        class ConventionalCL :
        public CatalystLayer<dim>
        {
        public:
            
            /**FcstUtilities
             * Concrete name used for objects of this class. This name is used when
             * setting up the subsection where the data is stored in the input file.
             * 
             * The data will be store under
             * \code
             * subsection name_specified_in_constructor
             *    set Material id = 4
             *    set Catalyst layer type = HomogeneousCL # <-here I select the type of object of type CatalystLayer
             *    subsection HomogeneousCL # <- this is the concrete_name for this class
             *       set all info relevant to this object
             *    end
             * end
             * \endcode
             */
            static const std::string concrete_name;
            
            ///@name Constructors, destructor, and initalization
            //@{                
            /** \warning For internal use only.
             * 
             * Constructor used only to create a prototype. Do not use
             * in general since this will not include the name of the section
             * in the parameter file you need. 
             */
            ConventionalCL();
            
            /**
             * Destructor
             */
            virtual ~ConventionalCL();
                       
            //@}
            
            ///@name Accessors and info
            //@{
            /**
             * Re-implementation of the parent set_local_material_id class to initialize the porosity
             * variable to the right value once the material ID is changed.
             */
            virtual inline void set_local_material_id(const unsigned int& id)
            {
                FuelCellShop::Layer::BaseLayer<dim>::set_local_material_id(id);
                if (!epsilon_V.empty())
                    this->porosity = epsilon_V.at(this->local_material_id());
            }
            /**
             * Print out the volume fraction in the catalyst layer
             */
            virtual void print_layer_properties() const;
            
             /** Get the volume fractions in the catalyst layer */
            virtual void get_volume_fractions(std::map<std::string, double>& volume_fractions)
            {
                //Assert( mat_id != numbers::invalid_material_id, ExcMessage("Graded Catalyst Layers needs its material id.") );
                
                compute_volume_fraction();
                volume_fractions["Solid"] = epsilon_S.at(this->local_material_id());
                volume_fractions["Void"] = epsilon_V.at(this->local_material_id());
                volume_fractions["Ionomer"] = epsilon_N.at(this->local_material_id());
                
            }
            
             /**
              * Return loadings
              * - \param V_Pt = Pt loading in ug/cm3
              * - \param loading_N = ionomer loading %wt
              * - \param IC_ratio = I/C ratio
              * @note either loading_N or IC_ratio should be passed through input file, the other computed
              * - \param prc_Pt = Pt/C ratio
              */
             virtual inline void get_loadings(std::map<std::string, double> & info)
             {
             	//Assert( mat_id != numbers::invalid_material_id, ExcMessage("Graded/Homogeneous Catalyst Layers needs its material id.") );
               info["V_Pt"] = V_Pt.at(this->local_material_id());
                 info["loading_N"] = loading_N.at(this->local_material_id());
                 info["IC_ratio"] = IC_ratio.at(this->local_material_id());
                 info["prc_Pt"] = prc_Pt.at(this->local_material_id());
             };
            
             /** Return the platinum loading per cm3 catalyst layer */
             inline double get_V_Pt(const unsigned int mat_id = numbers::invalid_material_id) const
             {
             	Assert( mat_id != numbers::invalid_material_id, ExcMessage("Graded/Homogeneous Catalyst Layers needs its material id.") );
             	// this is blah-blah return, change later
                 return V_Pt.at(this->local_material_id());
             }
            
             /**
              * Get the active area of platinum per unit volume of CL
              *
              * The active area is computed using the data in the input file according to the information in #compute_Av
              */
             virtual double get_active_area_Pt() const
             {
             	//Assert( mat_id != numbers::invalid_material_id, ExcMessage("Graded Catalyst Layers needs its material id.") );
             	// this is blah-blah return, change later
                 return Av.at(this->local_material_id());
             }
            
            //@}
            ///@name Effective property calculators
            //@{           
            
            /**
             * Compute the effective property in the pores of the CL. This is used for example to
             * compute effective diffusivity of gases. The method takes in bulk diffusion coefficient [\p m^2/s] and liquid water saturation as the first and second 
             * argument respectively. This routine is used in the isotropic case.
             */
            virtual void effective_gas_diffusivity(const double&,
            		                               const double&,
            		                               double&) const;

            /**
             * Return the effective diffusivity [\p m^2/s] for nonisothermal with/without two-phase case in the CL. It takes bulk diffusivity, computed using
             * #compute_gas_diffusion method and transforms it into an effective property,
             * taking into account the porosity, saturation and CL structure (Anisotropic case),
             * at all quadrature points of the cell.
             * 
             * \note: For two-phase case, set_saturation should be called before using this method, otherwise this method assumes saturation value to be zero.
             */
            virtual void effective_gas_diffusivity(std::vector< Tensor<2,dim> >&) const;
            /**
             * Return the derivative of effective diffusivity w.r.t solution variables/design parameters
             * for nonisothermal with/without two-phase case in the CL. It transforms bulk diffusion properties computed using
             * #compute_gas_diffusion method and transforms it into an effective property,
             * taking into account the porosity, saturation and CL structure (Anisotropic case),
             * at all quadrature points of the cell.
             * 
             * \note: For two-phase case, set_saturation should be called before using this method, otherwise this method assumes saturation value to be zero.
             */
            virtual void derivative_effective_gas_diffusivity(std::map< VariableNames, std::vector< Tensor<2,dim> > >&) const;
            
            /**
             * Compute the effective property in the pores of the CL. This is used to compute effective diffusivity of gases.
             * This routine can be used either in the isotropic or anisotropic cases.
             * Bulk diffusion coefficients or their derivatives are obtained from Mixure::BinaryDiffusion classes
             * inside this method.
             * 
             * \note The routine FuelCellShop::Layer::PorousLayer< dim >::set_gases_and_compute (std::vector< FuelCellShop::Material::PureGas * > &gases, double pressure, double temperature)
             * (in the parent class) should have been called prior to using this class. This method is to be used only for a single-phase, isothermal application.
             */
            virtual void effective_gas_diffusivity(Table< 2, Tensor< 2, dim > >&) const;
            
            /**
             * Compute the effective electron conductivity in the CL
             */
            virtual void effective_electron_conductivity(double&) const;
            
            /**
             * Compute the effective electron conductivity in the CL as an
             * anisotropic tensor.
             */
            virtual void effective_electron_conductivity(Tensor<2,dim>&) const;
            
            /**
             * Compute the derivative of the effective electron conductivity in the GDL
             * with respect to either the solution or design parameters. The parameters with respect to
             * which the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_effective_electron_conductivity(std::vector<double>&) const;
            
            /**
             * Compute the effective proton conductivity in the CL.
             */
            virtual void effective_proton_conductivity(double&) const;
            virtual void effective_proton_conductivity(std::vector<double>&) const;
            /**
             * Compute the derivative of the effective proton conductivity in the CL
             * with respect to either the solution or design parameters. The parameters with respect to
             * which the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_effective_proton_conductivity(std::map< VariableNames, std::vector<double> >&) const;
            
            /**
             * Compute the effective water diffusivity (lambda diffusivity) in the CL.
             */
            virtual void effective_water_diffusivity(double&) const;
            virtual void effective_water_diffusivity(std::vector<double>&) const;
            /**
             * Compute the derivative of the effective water diffusivity (lambda diffusivity) in the CL
             * with respect to either the solution or design parameters. The parameters with respect to
             * which the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_effective_water_diffusivity(std::map< VariableNames, std::vector<double> >&) const;
            
            /**
             * Compute the effective thermal conductivity of catalyst layer
             */
            virtual void effective_thermal_conductivity(double&) const;
            /**
             * Compute the effective thermal conductivity as a Tensor at all quadrature points
             */
            virtual void effective_thermal_conductivity(std::vector< Tensor<2,dim> >&) const;
            /**
             * Compute the derivative of the effective thermal conductivity in the CL. Currently, this function returns
             * only derivatives with respect to Temperature.
             */
            virtual void derivative_effective_thermal_conductivity(std::vector< Tensor<2,dim> >&) const;
            
            /**
             * Compute the effective thermo-osmotic diffusivity of lambda (sorbed water),
             * at all quadrature points in the CL.
             */
            virtual void effective_thermoosmotic_diffusivity(std::vector<double>&) const;
            /**
             * Compute the derivative of the effective thermo-osmotic diffusivity of lambda (sorbed water) in the CL
             * with respect to either the solution or design parameters. The parameters with respect to
             * which the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_effective_thermoosmotic_diffusivity(std::map< VariableNames, std::vector<double> >&) const;
            
            /**
             * Compute the anisotropic CL liquid permeability \f$ \left[ cm^2 \right] \f$, at all quadrature points in the cell.
             */
            virtual void liquid_permeablity(std::vector< Tensor<2,dim> >&) const;
            
            /**
             * Compute the derivative of the anisotropic liquid permeability in the CL
             * with respect to either the solution or design parameters, at all quadrature points in the cell. The parameters with respect to
             * which the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_liquid_permeablity(std::map< VariableNames, std::vector< Tensor<2,dim> > >&) const;
            
            
            
            
            
            
            
            
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
             * Compute \f$ \frac{\partial p_c}{\partial s} \quad \left[ dyne \cdot cm^{-2}\right] \f$, at all quadrature points in the CL.
             */
            virtual void dpcapillary_dsat(std::vector<double> &) const;
            /**
             * Compute the derivative of \f$ \frac{\partial p_c}{\partial s} \quad \left[ dyne \cdot cm^{-2}\right] \f$ in the CL, with 
             * respect to either the solution or design parameters, at all quadrature points in the cell. The parameters with respect to which 
             * the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_dpcapillary_dsat(std::map< VariableNames, std::vector<double> > &) const;
            
            /**
             * Compute the liquid-gas interfacial surface area per unit volume, \f$ a_{lv} ~\left[ \frac{cm^2}{cm^3} \right] \f$, at all 
             * quadrature points in the CL.
             */
            virtual void interfacial_surface_area(std::vector<double>&) const;
            /**
             * Compute the derivative of the liquid-gas interfacial surface area per unit volume, with respect to either the 
             * solution variables or design parameters, at all quadrature points in the CL. The parameters with respect to which 
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
            
        protected:
            ///@name Constructors
            //@{
                /**
             * Constructor
             */
            ConventionalCL(std::string name);
            
            /**
             * Declare parameters for a parameter file.
             * 
             * @note In order to enable graded catalyst layers, some parameters such as Platinum loading on support (%wt) utilize a 
             * map with the following structure material_id1:value1, material_id2:value2 where material_id1, material_id2
             * correspond to sublayers within the catalyst layer.
             * 
             * @code 
             * subsection Fuel cell data
             *   (...)
             *   subsection Cathode Catalyst Layer #<- This is the name fo the subsection that you specify in cl_section_name
             *     (...)
             *     subsection ConventionalCL #<- This is the name in (this->concrete_name)
             *       set Platinum loading on support (%wt) = 4:0.46                    # Mass percentage of platinum catalyst on the support carbon black
             *       set Platinum loading per unit volume (mg/cm3) = 4:400             # Catalyst platinum mass loading per unit volume of CL
             *       set Electrolyte loading (%wt)  = 4:0.3                            # Electrode loading is the weight percentage of ionomer per gram of CL
             *       //-- Network characteristics
             *       set Method effective transport properties in pores = Bruggemann # OPTIONS: Given|Bruggemann|Percolation -- Method used to compute effective transport properties in the void phase.
             *       set Porosity threshold  = 0.12                                  # Threshold value of the volume fraction of void space in the CL. If the porosity is less than this value transport does not occur
             *       set Porosity network constant = 2.0                             # Parameter used when using percolation theory      
             *       set Porosity gamma network constant = 0.0                       # Parameter used when using percolation theory to account for extra diffusion
             *       //--
             *       set Method effective transport properties in solid phase = Bruggemann # OPTIONS: Given|Bruggemann|Percolation --- Method used to compute effective transport properties in pores
             *       set Solid network threshold = 0.12                              # Threshold value of the volume fraction of solid (electron conductive) phase in the CL. If the solid phase is less than this value transport in the fibre network does not occur
             *       set Solid network constant = 2.0                                # Parameter used when using percolation theory
             *       set Method effective transport properties in electrolyte phase = Bruggemann # OPTIONS: Given|Bruggemann|Percolation|Iden11 -- Method used to compute effective transport properties in pores
             *       set Electrolyte network threshold = 0.12                        # Threshold value of the volume fraction of electrolyte (proton conductive) phase in the CL. If the electrolyte phase is less than this value transport in the network does not occur
             *       set Electrolyte network constant = 2.0                          # Parameter used when using percolation theory
             *       set Method to compute active area = given                       # OPTIONS: given|Marr|ETEK06|ETEK07 -- 
             *       set Active area [cm^2/cm^3] = 4:2.0e5
             *       set Method to compute porosity = marc                           #OPTIONS: marc
             *       //----
             *       set Method effective thermal conductivity = Given               # OPTIONS: Given -- Method used to compute effective thermal conductivity
             *       set Thermal conductivity, [W/(cm K)] = 4:0.015
             *     end
             *   end 
             * end
             * @endcode
             * 
             * For an explanation for the different options, please see the appropriate member function. For example, for the options relating
             * to active area see member function FuelCellShop::Layer::ConventionalCL::compute_Av()
             */
            void declare_parameters (const std::string& cl_section_name, 
                                     ParameterHandler &param) const;
                                  
            /**
             * Member function used to read in data and initialize the necessary data
             * to compute the coefficients.
             */
            void initialize (ParameterHandler &param);
            //@}
            /**
             * Compute porosity and volume fraction of solid and ionomer in the catalyst layer
             */
            void compute_volume_fraction();
            
            /**
             * Compute the active area of catalyst in the layer by the specified method
             * 
             * There are three methods to compute the active area, namely
             * - given
             * - Marr
             * - ETEK06
             * - ETEK07
             * 
             * The given option will simply use the value provided by the user.
             * 
             * The Marr, ETEK06 and  ETEK07 options use a polyomial approximation to obtain the surface area of Pt per gram of catalyst, i.e. \f$  A_0 \f$ , 
             * based on the Pt to carbon weight ratio, i.e. Pt|C or Platinum loading on support (%wt) in the input file. 
             * Then, the cm^2 Pt per cm^3 of CL is obtained using
             * \f[
             * A_v = A_0*(V_{Pt}*1e-3)
             * \f]
             * where \f$ V_{Pt} \f$ is the input parameter "Platinum loading per unit volume (mg/cm3)".
             * 
             * The functions to obtain \f$ A_0 \f$ are
             * 
             * 1) Marr 
             * \f[
             * A_0 = 2.2779e6*(Pt|C)^3 - 1.5857e6*(Pt|C)^2 - 2.0153e6*(Pt|C) + 1.5950e6;
             * \f]
             * The equation is from the following reference: 
             * Marr, C., Li, X. Composition and performance modelling of catalyst layer in a proton exchange membrane fuel cell. 
             * Journal of Power Sources 77 (1) , pp. 17-27, 1999
             * 
             * 2) ETEK06 
             *  \f[
             * A_0 = -4.5646e5*(Pt|C)^3 + 1.0618e6*(Pt|C)^2 - 1.8564e6*(Pt|C) + 1.5955e6;
             * \f]
             * and is a curve-fit to ETEK catalyst data reported previous to 2006.
             * 
             * 3) ETEK07
             * \f[
             * A_0 = 0.7401e7*(Pt|C)^4 - 1.8105e7*(Pt|C)^3 + 1.5449e7*(Pt|C)^2 - 0.6453e7*(Pt|C) + 0.2054e7
             *  \f]
             * The equation is derived in M. Secanell, Computational Modeling and Optimization of Proton 
             * Exchange Membrane Fuel Cells, Ph.D. thesis, University of Victoria, 2008 by curve-fitting data from ETEK from their
             * 2007 catalyst. It is used in several of M. Secanell's publications.
             * 
             */ 
            void compute_Av();
             /**
              * Compute the derivative of the effective proton conductivity w.r.t. the electrolyte
              * loading.
              */
             void derivative_effective_proton_conductivity_wrt_electrolyte_loading(double&) const;
            
             /**
              * Function to compute the partial derivative of the volume fraction
              * the different phases in the catalyst layer with respect to
              * the design variables of the optimization problem
              */
             void derivative_volume_fractions(double &Depsilon_S,
                                              double &Depsilon_V,
                                              double &Depsilon_N) const;
                                             
            /**
            * Get the effective transport method in the pores
            */
            void get_method_transport_property_pores(std::string& method) const
            {
                method = method_eff_property_pores;
            };
            
            /**
            * Get the effective transport method in the electrolyte
            */
            void get_method_transport_property_electrolyte(std::string& method) const
            {
                method = method_eff_property_electrolyte;
            };
            
            /**
            * Get the effective transport method in the solid phase
            */
            void get_method_transport_property_solid(std::string& method) const
            {
                method = method_eff_property_solid;
            };

            /**
            * Inline function to compute
            * \f[
            * \frac{\partial \epsilon_S^{cat}}{\partial \%Pt}
            * \f]
            */
            inline double depsilon_S_cat_dprc_Pt(const double& V_Pt,
                                                 const double& prc_Pt) const
            {
                return -(V_Pt*1e-3)/(rho_c*pow(prc_Pt,2.0));
            }

            /**
            * Inline function to compute
            * \f[
            * \frac{\partial \epsilon_S^{cat}}{\partial m_{Pt}}
            * \f]
            */
            inline double depsilon_S_cat_dVPt(const double& prc_Pt) const
            {
                return  (1/rho_Pt + (1-prc_Pt)/(prc_Pt*rho_c))*(1e-3);
            }

            /**
            * Inline function to compute
            * \f[
            * \frac{\partial \epsilon_V^{cat}}{\partial \epsilon_S^{cat}}
            * \f]
            */
            inline double depsilon_V_cat_depsilon_S_cat() const
            {
                return -1;
            }

            /**
            * Inline function to compute
            * \f[
            * \frac{\partial \epsilon_V^{cat}}{\partial \epsilon_S^{cat}}
            * \f]
            */
            inline double depsilon_V_cat_depsilon_N_cat() const
            {
                return -1;
            }

            //-- Composition
            /** Volume fraction of Nafion in the cathode catalyst layer. Using unsigned int specifying the material id for the specific CL*/
            std::map< unsigned int, double> epsilon_N;
            /** Void volume fraction (Porosity) of the catalyst layer. Using unsigned int specifying the material id for the specific CL*/
            std::map< unsigned int, double> epsilon_V;
            /** Solid volume fraction in the catalyst layer. Using unsigned int specifying the material id for the specific CL*/
            std::map< unsigned int, double> epsilon_S;
            /** Volume fraction of water in the cathode catalyst layer. Using unsigned int specifying the material id for the specific CL*/
            std::map< unsigned int, double> epsilon_W;
            
            //-- Catalyst properties
            /** Density of platinum */
            double rho_Pt;
            /** Density of support material */
            double rho_c;
            /** Percentage of platinum per carbon on the catalyst layer. Using unsigned int specifying the material id for the specific CL*/
            std::map< unsigned int, double> prc_Pt;
            /** Platinum loading at the catalyst layer per unit volume. Using unsigned int specifying the material id for the specific CL*/
            std::map< unsigned int, double> V_Pt;
            /** Platinum loading at the catalyst layer per unit area. Using unsigned int specifying the material id for the specific CL*/
            std::map< unsigned int, double> M_Pt;
            /** Active area of catalyst per unit volume of catalyst layer. Using unsigned int specifying the material id for the specific CL*/
            std::map< unsigned int, double> Av;
            /** Method to compute active area */
            std::string method_Av;
            /** Method to compute porosity */
            std::string method_porosity;
            /** Layer thickness or thicknesses. Using unsigned int specifying the material id*/
            std::map< unsigned int, double> L_CL;
            
            //-- Electrolyte properties
            /** Density of electrolyte */
            double rho_N;
            /** Electrolyte loading. Electrode loading is the weight percentage 
            * of ionomer per gram of CL. Which may change if multiple CLs exist. Each CL is specified using unsigned int for the material id.
            * loading_N = weight electrolyte / (weight Pt + weight C + weight electrolyte)
            */
            std::map< unsigned int, double> loading_N;
            /**
             * Ionomer to carbon ratio. Which may change if multiple CLs exist. Each CL is specified using unsigned int for the material id.
             * I/C ratio = weight electrolyte / weight C
             */
            std::map< unsigned int, double> IC_ratio;
            
            /** Percentage (mass fraction) of electrolyte in the catalyst layer. Each CL is specified using unsigned int for the material id. */
            std::map< unsigned int, double> prc_N;
            
            // Network characteristics
            /** Method used to compute effective properties -- Type of network */
            std::string method_eff_property_pores;
            /** Porous network threshold */
            double porosity_th;
            /** Porous network constant */
            double porosity_mu;
            /** */
            double porosity_gamma;
            
            /** Method used to compute effective properties -- Type of network */
            std::string method_eff_property_solid;
            /** Solid phase network threshold */
            double solid_th;
            /** Solid phase network constant */
            double solid_mu;
            /** Method used to compute effective properties -- Type of network */
            std::string method_eff_property_electrolyte;
            /** Electrolyte network threshold */
            double electrolyte_th;
            /** Electrolyte network constant */
            double electrolyte_mu;

            
            /** Method used to compute effective thermal conductivity in the catalyst layer */
            std::string method_eff_thermal;
            /** Thermal Conductivity of the layer. Using unsigned int specifying the material id for the specific CL*/
            std::map< unsigned int, double> k_T;
            
            /** Method used to compute the relative liquid permeability. */
            std::string method_rel_liquid_permeability;
            /** Irreducible liquid water saturation value in the MPL. */
            std::map< unsigned int, double> s_irr;
            /** Absolute permeability [cm^2] of the layer. Using unsigned int specifying the material id for the specific CL*/
            std::map< unsigned int, double> abs_permeability;
            
            /** Method used to compute capillary pressure as a function of saturation. */
            std::string method_capillary_function;
            
            
            
            
            
            
            
            
            
            
        };
    }
}

#endif
