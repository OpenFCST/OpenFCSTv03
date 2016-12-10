//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: gas_diffusion_layer.h
//    - Description: Base Gas Diffusion Layer Class. It implements the interface for other gas diffusion layer classes
//        and some common methods.
//    - Developers: M. Secanell
//    - $Id: gas_diffusion_layer.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__GAS_DIFFUSION_LAYER_H
#define _FUELCELLSHOP__GAS_DIFFUSION_LAYER_H

// FCST classes
#include <utils/fcst_constants.h>
#include <materials/base_material.h>
#include <layers/base_layer.h>
#include <layers/porous_layer.h>

using namespace dealii;

namespace FuelCellShop
{
    namespace Layer
    {
        /**
         * Virtual class used to provide the interface for all GasDiffusionLayer children. 
         * 
         * No object of type GasDiffusionLayer should ever be created, instead this layer
         * is used to initialize pointers of type GasDiffusionLayer. The class has a database of
         * children such that it will declare all necessary parameters for all children in the
         * input file, read the input file, create the appripriate children and return a pointer
         * to GasDiffusionLayer with the children selected.
         * 
         * All public functions are virtual but the static functions used to declare parameters and
         * to initialize a pointer of GasDiffusionLayer, i.e. declare_all_GasDiffusionLayer_parameters,
         * set_all_GasDiffusionLayer_parameters and create_GasDiffusionLayer.
         * 
         * <h3>Usage Details:</h3>
         * 
         * In order to create a gas diffusion layer within an application, the following steps need to be taken.
         * 
         * First, in the application .h file, create a pointer to a GasDiffusionLayer object, i.e.
         * @code
         * boost::shared_ptr<FuelCellShop::Layer::GasDiffusionLayer<dim> > CGDL;
         * @endcode
         * 
         * This pointer object will be available anywhere inside the application. Because we do not want to
         * worry about deleting the pointer afterwards, we use a Boost pointer which has its own memory management
         * algorithms. See the <a href="http://www.boost.org/doc/libs/1_54_0/libs/smart_ptr/smart_ptr.htm">Boost website</a> for more information
         * 
         * Once the pointer is available, we need to do three things in the application
         * - Call declare_all_GasDiffusionLayer_parameters in the application declare_parameters() member function. This 
         * member function is used to define all parameters that can be read from the input file
         * 
         * - Call create_GasDiffusionLayer and initialize. The former member function will fill the pointer created above with the appropriate
         * gas diffusion layer you have selected in the input file. Then, initialize reads from the input file all the data from the file in order
         * to setup the object.
         * 
         * The object is ready for use now.
         * 
         * Here is a code example from app_cathode.cc:
         * @code
         * //--------- IN DECLARE_PARAMETERS ------------------------------------------------------
         * template <int dim>
         * void 
         * NAME::AppCathode<dim>::declare_parameters(ParameterHandler& param)
         * {
         *   (...)
         *   // Declare section on the input file where all info will be stored. In this case Fuel Cell Data > Cathode gas diffusion layer
         *   FuelCellShop::Layer::GasDiffusionLayer<dim>::declare_GasDiffusionLayer_parameters("Cathode gas diffusion layer", param);
         *   (...)
         * }
         * 
         * //--------- IN INITIALIZE ------------------------------------------------------
         * template <int dim>
         * void
         * NAME::AppCathode<dim>::_initialize(ParameterHandler& param)
         * {   
         *  (...) 
         *   // Initialize layer classes:
         *   std::vector< FuelCellShop::Material::PureGas * > gases;
         *   gases.push_back(&oxygen);
         *   gases.push_back(&nitrogen);
         *    
         *   CGDL = FuelCellShop::Layer::GasDiffusionLayer<dim>::create_GasDiffusionLayer("Cathode gas diffusion layer",param);
         *   // Here, I specify the gases that exist in the GDL and their temperature and pressure (based on operating conditions):
         *   CGDL->set_gases_and_compute(gases, OC.get_pc_atm (), OC.get_T());
         *   (...)
         * }
         * @endcode
         * 
         * 
         * \author M. Secanell, 2013
         * 
         */
        template <int dim>
        class GasDiffusionLayer : 
        public  PorousLayer<dim>
        {
        public:
            ///@name Instance Delivery (Public variables)
            //@{         
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
            //@}             
            
            ///@name Instance Delivery (Public functions)
            /**
             * Function used to declare all the data necessary in the parameter files for
             * all GasDiffusionLayer children.
             * 
             * This member function should be used instead of declare_parameters() when we want
             * to use a GasDiffusionLayer pointer that selects the type of GDL to run at runtime.
             * 
             * \param gld_section_name Name of the section that will encapuslate all the information about the GDL
             * \param param ParameterHandler object used to store all information about the simulation. Used
             * to read the parameter file.
             * 
             * The parameter file would look as follows:
             * 
             * @code
             * subsection Fuel cell data
             * (...)
             * subsection Cathode gas diffusion layer           #This is the name provided in gld_section_name
             * 
             *   set Gas diffusion layer type = DummyGDL #[ DesignFibrousGDL | DummyGDL | SGL24BA ]
             *   set Material id = 2
             *
             *     subsection DummyGDL         # This is the subsection for the first children
             *       set Oxygen diffusion coefficient, [cm^2/s] = 0.22
             *       set Electrical conductivity, [S/cm] = 40
             *     end
             *
             *     subsection DesignFibrousGDL
             *       set Porosity = 0.6
             *       ## Anisotropy
             *       set Anisotropic transport = true                    # (default) false
             *       set Method effective transport properties in pores = Tomadakis  # (default) Bruggemann | Given | Percolation | Tomadakis | Mezedur
             *       set Method effective transport properties in solid = Percolation    # (default) Bruggemann | Given | Percolation
             *       (...)
             *     end
             * end
             * (...)
             * end
             * @endcode
             * 
             */
            static void declare_GasDiffusionLayer_parameters (const std::string& gdl_section_name, ParameterHandler &param)
            {
                
                for (typename FuelCellShop::Layer::GasDiffusionLayer<dim>::_mapFactory::iterator iterator = FuelCellShop::Layer::GasDiffusionLayer<dim>::get_mapFactory()->begin(); 
                     iterator != FuelCellShop::Layer::GasDiffusionLayer<dim>::get_mapFactory()->end(); 
                iterator++)
                     {
                         iterator->second->declare_parameters(gdl_section_name, param);
                     }        
            }

             /**
              * 
              * Function used to select the appropriate GasDiffusionLayer type as specified in the ParameterHandler under
              * line 
              * @code 
              * set Gas diffusion layer type = DummyGDL 
              * @endcode
              * current options are [ DesignFibrousGDL | DummyGDL | SGL24BA ]
              * 
              * The class will read the appropriate section in the parameter file, i.e. the one with name \param gld_section_name ,
              * create an object of the desired type and return it.
              * 
              */
             static boost::shared_ptr<FuelCellShop::Layer::GasDiffusionLayer<dim> > create_GasDiffusionLayer (const std::string& gld_section_name, 
                                                                                                              ParameterHandler &param)
             {       
                 
                 boost::shared_ptr<FuelCellShop::Layer::GasDiffusionLayer<dim> > pointer; 
                 
                 std::string concrete_name;
                 param.enter_subsection("Fuel cell data");
                 {
                     param.enter_subsection(gld_section_name);
                     {
                         concrete_name = param.get("Gas diffusion layer type");
                         FcstUtilities::log << "name: "<<concrete_name.c_str()<<std::endl;
                     }
                     param.leave_subsection();
                 }
                 param.leave_subsection();
                 
                 typename FuelCellShop::Layer::GasDiffusionLayer<dim>::_mapFactory::iterator iterator = FuelCellShop::Layer::GasDiffusionLayer<dim>::get_mapFactory()->find(concrete_name);
                 
                 if (iterator != FuelCellShop::Layer::GasDiffusionLayer<dim>::get_mapFactory()->end())
                 {
                     if (iterator->second)
                     {
                         pointer = iterator->second->create_replica(gld_section_name);
                     }
                     else 
                     {
                         FcstUtilities::log<<"Pointer not initialized"<<std::endl;
                         abort();
                     }
                 }
                 else
                 {
                     FcstUtilities::log<<"Concrete name in FuelCellShop::Layer::GasDiffusionLayer<dim>::create_GasDiffusionLayer does not exist"<<std::endl;
                     abort();
                 }
                 
                 pointer->initialize(param);
                 
                 return pointer;
             }
            //@}             

            ///@name Initalization
            //@{             
            
            /**
             * Member function used by some applications such as dummyGDL in order to know which value
             * to return. For other classes this class is not used.
             * 
             * \deprecated This member function should not be used.
             */
            virtual void set_diffusion_species_name(std::string &name)
            {
                Assert((name != "oxygen" || name != "nitrogen" || name != "water" || name != "electron"),
                       ExcNotImplemented());
                diffusion_species_name = name;
            };
            //@}
            
            ///@name Accessors and info
            //@{
            /**
             * This member function returns a type_info object with the name of the 
             * base layer type the inherited class belongs to, i.e.
             * - GasDiffusionLayer
             * - MicroPorousLayer
             * - CatalystLayer
             * - MembraneLayer
             * - SolidLayer
             * - Channel
             * 
             * Note that this is necessary if we want to find out not the name of the actual class which can be obtain using
             * @code const std::type_info& name = typeid(*this) @endcode
             * but the name of the parent class.
             * 
             * @note Do not re-implement this class in children classes
             */
            const std::type_info& get_base_type() const
            {
                return typeid(GasDiffusionLayer<dim>);
            }
            
            /**
             * Class test
             */
            virtual void test_class() const;
            //@}
            
            ///@name Effective property calculators
            //@{

            /**
             * Compute the effective property in the pores of the GDL. This is used for example to
             * compute effective diffusivity of gases. The method takes in bulk diffusion coefficient [\p m^2/s] and liquid water saturation as the first and second 
             * argument respectively. This routine is used in the isotropic case.
             */
            virtual void effective_gas_diffusivity(const double&, const double&, double&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
                
            };
            /**
             * Compute the effective property in the pores of the GDL. This is used for example to
             * compute effective diffusivity of gases. The method takes in bulk diffusion coefficient [\p m^2/s] and liquid water saturation as the first and second 
             * argument respectively. This routine can be used either in the isotropic or anisotripic case
             */
            virtual void effective_gas_diffusivity(const double&, const double&, Tensor<2,dim>&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            /**
             * Return the effective diffusivity [\p m^2/s] for nonisothermal with/without two-phase case in the GDL. It takes bulk diffusivity, computed using
             * #compute_gas_diffusion method and transforms it into an effective property,
             * taking into account the porosity, saturation and GDL structure (Anisotropic case),
             * at all quadrature points of the cell.
             */
            virtual void effective_gas_diffusivity(std::vector< Tensor<2,dim> >&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            /**
             * Return the derivative of effective diffusivity w.r.t solution variables/design parameters
             * for nonisothermal with/without two-phase case in the GDL. It transforms bulk diffusion properties computed using
             * #compute_gas_diffusion method and transforms it into an effective property,
             * taking into account the porosity, saturation and GDL structure (Anisotropic case),
             * at all quadrature points of the cell.
             */
            virtual void derivative_effective_gas_diffusivity(std::map< VariableNames, std::vector< Tensor<2,dim> > >&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            /**
             * Return the effective diffusivty in the GDL for all the gases assigned to the
             * layer using set_gases_and_compute. This routine uses FuelCellShop::Mixture::ChapmanEnskog to compute
             * the binary diffusivity for each gas with respect to each other. If the layer contains three gases the 
             * vector returns D12, D13, D23. For 4 gases, it returns D12, D13, D14, D23, D24, D34.
             * 
             * \deprecated Use #compute_gas_diffusion with appropriate gases and then #effective_gas_diffusivity
             */
            virtual void effective_gas_diffusivity(Table<2, double>& D_eff ) const 
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "!!!!!!!!!!!!!!! This class will die soon. Do not use. It is a terrible interface. !!!!!!!!!!!!!!!!!!!!!!!!!! "
                << info.name() << std::endl;
            };
            /**
             * Return a tensor with the effective diffusivty in the GDL for all the gases assigned to the
             * layer using set_gases_and_compute. This routine uses FuelCellShop::Mixture::ChapmanEnskog to compute
             * the binary diffusivity for each gas with respect to each other. If the layer contains three gases the 
             * vector returns D12, D13, D23. For 4 gases, it returns D12, D13, D14, D23, D24, D34 
             * (Anisotropic case).
             * 
             * \deprecated Use #compute_gas_diffusion with appropriate gases and then #effective_gas_diffusivity
             */
            virtual void effective_gas_diffusivity(Table< 2, Tensor<2,dim> > &D_eff ) const 
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "!!!!!!!!!!!!!!!!!!!! This class will die soon. Do not use. It is a terrible interface. !!!!!!!!!!!!!!!!!!!!!!!! "
                << info.name() << std::endl;
            };
            
            /**
             * Compute the effective electron conductivity in the GDL
             */
            virtual void effective_electron_conductivity(double& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            /**
             * Compute the effective electron conductivity in the GDL
             */
            virtual void effective_electron_conductivity(const double& /*prop*/, double& /*prop_eff*/) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            /**
             * Compute the effective electron conductivity in the GDL
             */
            virtual void effective_electron_conductivity(Tensor<2,dim>& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            /**
             * Compute the effective electron conductivity in the GDL
             */
            virtual void effective_electron_conductivity(const double& /*prop*/, Tensor<2,dim>& /*prop_eff*/) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            
            /**
             * Compute the derivative of the effective electron conductivity in the GDL
             * with respect to either the solution or design parameters. The parameters with respect to
             * which the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_effective_electron_conductivity(std::vector<Tensor<2,dim> >& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the effective thermal conductivity (isotropic) in the GDL
             */
            virtual void effective_thermal_conductivity(double& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the effective thermal conductivity (anisotropic) in the GDL
             */
            virtual void effective_thermal_conductivity(Tensor<2,dim>& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            /**
             * Compute the effective thermal conductivity (anisotropic) in the GDL, dependent on various solution variables, eg: Temperature
             */
            virtual void effective_thermal_conductivity(std::vector< Tensor<2,dim> >&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the derivative of the effective thermal conductivity in the GDL
             * with respect to either the solution or design parameters. The parameters with respect to
             * which the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_effective_thermal_conductivity(std::vector<Tensor<2,dim> >& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the GDL gas permeability
             */
            virtual void gas_permeablity(double& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the GDL gas permeability
             */
            virtual void gas_permeablity(Tensor<2,dim>& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the derivative of the effective gas permeability in the GDL
             * with respect to either the solution or design parameters. The parameters with respect to
             * which the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_gas_permeablity(std::vector<Tensor<2,dim> >& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            
            /**
             * Compute the anisotropic GDL liquid permeability \f$ \left[ cm^2 \right] \f$, at all quadrature points in the cell.
             */
            virtual void liquid_permeablity(std::vector< Tensor<2,dim> >& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the derivative of the anisotropic liquid permeability in the GDL
             * with respect to either the solution or design parameters, at all quadrature points in the cell. The parameters with respect to
             * which the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_liquid_permeablity(std::map< VariableNames, std::vector< Tensor<2,dim> > >& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the derivative of the anisotropic liquid permeability in the GDL
             * with respect to either the solution or design parameters, at all quadrature points in the cell. The parameters with respect to
             * which the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void relative_liquid_permeability_PSD(std::vector< Tensor<2,dim> >&  ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            virtual void derivative_relative_liquid_permeablity_PSD(std::vector<double> & ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            virtual void derivative_relative_liquid_permeablity_PSD(std::map< VariableNames, std::vector< Tensor<2,dim> > > & ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            virtual void saturated_liquid_permeablity_PSD(double&  ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            /**
             * Compute \f$ p_c \quad \left[ dyne \cdot cm^{-2}\right] \f$, at all quadrature points in the cell.
             */
            virtual void pcapillary(std::vector<double>&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl; 
            };
            
            virtual void saturation_from_capillary_equation(std::vector<double>&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl; 
            };
            
            virtual void derivative_saturation_from_capillary_equation_PSD(std::vector<double>&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl; 
            };
            /**
             * Compute \f$ \frac{\partial p_c}{\partial s} \quad \left[ dyne \cdot cm^{-2}\right] \f$, at all quadrature points in the GDL.
             */
            virtual void dpcapillary_dsat(std::vector<double> &) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the derivative of \f$ \frac{\partial p_c}{\partial s} \quad \left[ dyne \cdot cm^{-2}\right] \f$ in the GDL, with 
             * respect to either the solution or design parameters, at all quadrature points in the cell. The parameters with respect to which 
             * the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_dpcapillary_dsat(std::map< VariableNames, std::vector<double> > &) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            /**
             * Compute the liquid-gas interfacial surface area per unit volume, \f$ a_{lv} ~\left[ \frac{cm^2}{cm^3} \right] \f$, at all 
             * quadrature points in the GDL.
             */
            virtual void interfacial_surface_area(std::vector<double>&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the derivative of the liquid-gas interfacial surface area per unit volume, with respect to either the 
             * solution variables or design parameters, at all quadrature points in the GDL. The parameters with respect to which 
             * the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags().
             */
            virtual void derivative_interfacial_surface_area(std::map< VariableNames, std::vector<double> >&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            /**
             * Compute the liquid-gas interfacial surface area per unit volume, \f$ a_{lv} ~\left[ \frac{cm^2}{cm^3} \right] \f$, at all 
             * quadrature points in the CL.
             */
            virtual void interfacial_surface_area_PSD(std::vector<double>&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the derivative of the liquid-gas interfacial surface area per unit volume, with respect to either the 
             * solution variables or design parameters, at all quadrature points in the CL. The parameters with respect to which 
             * the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags().
             */
            virtual void derivative_interfacial_surface_area_PSD(std::map< VariableNames, std::vector<double> >&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            virtual void derivative_interfacial_surface_area_PSD(std::vector<double>&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };            
            //@}

        protected:
            ///@name Instance Delivery (Types)
            //@{         
            /** 
             * This object is used to store all objects of type GasDiffusionLayer. 
             */
            typedef std::map< std::string, GasDiffusionLayer<dim>* > _mapFactory;      
            //typedef std::map< std::string, boost::shared_ptr<GasDiffusionLayer<dim> > > _mapFactory;      
            //@}       
            
            ///@name Instance Delivery (Private and static)
            //@{         
            /**
             * 
             */
            static _mapFactory * get_mapFactory()
            {
                static _mapFactory mapFactory;
                return &mapFactory;
            }  
            //@}            
            ///@name Constructors, destructor, and initalization
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
            GasDiffusionLayer();
                       
            /**
             * Destructor
             */
            ~GasDiffusionLayer();     
            /**
             * Constructor
             */
            GasDiffusionLayer(const std::string& name);
            
            /**
             * Declare parameters for a parameter file
             * 
             * \deprecated Use declare_all_GasDiffusionLayer_parameters
             * 
             */
            void declare_parameters (ParameterHandler &param) const
            {
                this->declare_parameters(this->name, param);
            }   
            
            /**
             * Declare parameters for a parameter file. 
             * 
             */            
            virtual void declare_parameters (const std::string& name, 
                                             ParameterHandler &param) const;           
            
            /**
             * Member function used to read in data and initialize the necessary data
             * to compute the coefficients.
             * 
             * \note This routine is called in create_GasDiffusionLayer.
             */
            void initialize (ParameterHandler &param);
            //@}
                       
            ///@name Instance Delivery (Private functions)
            //@{         
            /**
             * This member function is used to create an object of type gas diffusion layer
             * 
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::Layer::GasDiffusionLayer<dim> > create_replica (const std::string &name)
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            //@}
            
            ///@name Internal variables
            //@{
            /** If GDL properties are stored inside the class (e.g. DummyGDL) then, return the property
             * stored under coefficient_name name
             */
            std::string diffusion_species_name;	
            
            /** PSD class type from input file */
            std::string PSD_type;
            
            /**
             * Porosity of the GDL
             */
            double porosity;
            
            /**
             * Tortuosity tensor of the GDL
             */
            Tensor<2,dim> tortuosity_tensor;
            
            /** 
             * Tensor storing the effective thermal conductivity 
             * of the layer.
             */
            Tensor<2, dim> thermal_conductivity_tensor;
            
            /**
             * Double storing the electric conductivity of the GDL is the
             * layer is isotropic. Otherwise the value is set to zero.
             */
            double electron_conductivity;
            
            /** 
             * Tensor storing the effective electronic conductivity 
             * of the layer. This function if used if the layer is anisotropic.
             */
            Tensor<2, dim> electron_conductivity_tensor;
            
            //@}
            
        };
        
    }
    
}  // FuelCellShop

#endif // _FUELCELLSHOP__GAS_DIFFUSION_LAYER_H
