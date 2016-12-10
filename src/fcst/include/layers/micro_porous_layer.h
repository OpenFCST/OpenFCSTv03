//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: micro_porous_layer.h
//    - Description: Base Micro-Porous Layer Class. It implements the interface for other micro-porous layer classes
//        and some common methods.
//    - Developers: M. Secanell
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__MICRO_POROUS_LAYER_H
#define _FUELCELLSHOP__MICRO_POROUS_LAYER_H

// FCST classes
#include <materials/base_material.h>
#include <layers/porous_layer.h>

using namespace dealii;

namespace FuelCellShop
{
    namespace Layer
    {
        /**
         * Virtual class used to provide the interface for all MicroPorousLayer children. 
         * 
         * No object of type MicroPorousLayer should ever be created, instead this layer
         * is used to initialize pointers of type MicroPorousLayer. The class has a database of
         * children such that it will declare all necessary parameters for all children in the
         * input file, read the input file, create the appripriate children and return a pointer
         * to MicroPorousLayer with the children selected.
         * 
         * All public functions are virtual but the static functions used to declare parameters and
         * to initialize a pointer of MicroPorousLayer, i.e. declare_all_MicroPorousLayer_parameters,
         * set_all_MicroPorousLayer_parameters and create_MicroPorousLayer.
         * 
         * 
         * <h3>Usage Details:</h3>
         * 
         * In order to create a micro porous layer within an application, the following steps need to be taken.
         * 
         * First, in the application .h file, create a pointer to a MicroPorousLayer object, i.e.
         * @code
         * boost::shared_ptr<FuelCellShop::Layer::MicroPorousLayer<dim> > CMPL;
         * @endcode
         * 
         * This pointer object will be available anywhere inside the application. Because we do not want to
         * worry about deleting the pointer afterwards, we use a Boost pointer which has its own memory management
         * algorithms. See the <a href="http://www.boost.org/doc/libs/1_54_0/libs/smart_ptr/smart_ptr.htm">Boost website</a> for more information
         * 
         * Once the pointer is available, we need to do three things in the application
         * - Call declare_MicroPorousLayer_parameters in the application declare_parameters() member function. This 
         * member function is used to define all parameters that can be read from the input file
         * 
         * - Call create_MicroPorousLayer and initialize. The former member function will fill the pointer created above with the appropriate
         * microporous layer you have selected in the input file. Then, initialize reads from the input file all the data from the file in order
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
         *   // Declare section on the input file where all info will be stored. In this case Fuel Cell Data > Cathode microporous layer
         *   FuelCellShop::Layer::MicroPorousLayer<dim>::declare_MicroPorousLayer_parameters("Cathode microporous layer", param);
         *   (...)
         * }
         * 
         * //--------- IN INITIALIZE ------------------------------------------------------
         * template <int dim>
         * void
         * NAME::AppCathode<dim>::_initialize(ParameterHandler& param)
         * {    
         *   (...)
         *   // Initialize layer classes:
         *    std::vector< FuelCellShop::Material::PureGas * > gases;
         *    gases.push_back(&oxygen);
         *    gases.push_back(&nitrogen);
         *    
         *    CMPL = FuelCellShop::Layer::MicroPorousLayer<dim>::create_MicroPorousLayer("Cathode microporous layer",param);
         * 
         *    // Here, I specify the gases that exist in the MPL and their temperature and pressure (based on operating conditions):
         *    CMPL->set_gases_and_compute(gases, OC.get_pc_atm (), OC.get_T());
         *   (...)
         * }
         * @endcode
         *
         * \author M. Secanell, 2013
         */
        template <int dim>
        class MicroPorousLayer : 
        public  PorousLayer<dim>
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
             *    set Microporous layer type = DummyMPL # <-here I select the type of object of type MicroPorousLayer
             *    subsection DummyMPL # <- this is the concrete_name for this class
             *       set all info relevant to this object
             *    end
             * end
             * \endcode
             */
            static const std::string concrete_name;
            
            //@}         

            ///@name Instance Delivery
            //@{         
            /**
             * Function used to declare all the data necessary in the parameter files former
             * all MicroPorousLayer children.
             * 
             * This member function should be used instead of declare_parameters() when we want
             * to use a MicroPorousLayer pointer that selects the type of MPL to run at runtime.
             * 
             * \param mpl_section_name Name of the section that will encapuslate all the information about the CL
             * \param param ParameterHandler object used to store all information about the simulation. Used
             * to read the parameter file.
             * 
             * The parameter file would look as follows:
             * 
             * @code
             * 
             * subsection Fuel cell data
             * (...)
             *   subsection mpl_section_name  <- This is the name used in the parameter file
             *     set Micro porous layer type = SGL24BC  # Types available: SGL24BC | DesignMPL
             *   end
             * end
             * 
             * @endcode    
             * 
             */
            static void declare_MicroPorousLayer_parameters (const std::string& mpl_section_name, 
                                                                 ParameterHandler &param)
            {
                for (typename FuelCellShop::Layer::MicroPorousLayer<dim>::_mapFactory::iterator iterator = FuelCellShop::Layer::MicroPorousLayer<dim>::get_mapFactory()->begin(); 
                     iterator != FuelCellShop::Layer::MicroPorousLayer<dim>::get_mapFactory()->end(); 
                iterator++)
                     {
                         iterator->second->declare_parameters(mpl_section_name, param);
                     }   
            }
                       
            /**
             * Function used to select the appropriate MicroPorousLayer
             */
            static boost::shared_ptr<FuelCellShop::Layer::MicroPorousLayer<dim> > create_MicroPorousLayer (const std::string& mpl_section_name, 
                                                                                                           ParameterHandler &param)
            {
                boost::shared_ptr<FuelCellShop::Layer::MicroPorousLayer<dim> > pointer;
                
                std::string concrete_name;
                param.enter_subsection("Fuel cell data");
                {
                    param.enter_subsection(mpl_section_name);
                    {
                        concrete_name = param.get("Micro porous layer type");
                        FcstUtilities::log << "name: "<<concrete_name.c_str()<<std::endl;
                    }
                    param.leave_subsection();
                }
                param.leave_subsection();
                
                typename FuelCellShop::Layer::MicroPorousLayer<dim>::_mapFactory::iterator iterator = FuelCellShop::Layer::MicroPorousLayer<dim>::get_mapFactory()->find(concrete_name);
                
                if (iterator != FuelCellShop::Layer::MicroPorousLayer<dim>::get_mapFactory()->end())
                {
                    if (iterator->second)
                    {
                        pointer = iterator->second->create_replica(mpl_section_name);
                    }
                    else 
                    {
                        FcstUtilities::log<<"Pointer not initialized"<<std::endl;
                        abort();
                    }
                }
                else
                 {
                     FcstUtilities::log<<"Concrete name in FuelCellShop::Layer::MicroPorousLayer<dim>::create_MicroPorousLayer does not exist"<<std::endl;
                     abort();
                 }
                
                pointer->initialize(param);
                
                return pointer;
            }
            //@}
            
            ///@name Initalization
            //@{    
            
            /**
             * Declare parameters for a parameter file. 
             * 
             * \deprecated Use declare_MicroPorousLayer_parameters
             */
            //void declare_parameters (ParameterHandler &param) const
            //{
            //    FuelCellShop::Layer::MicroPorousLayer<dim>::declare_parameters(this->name,param);
            //}
            
            /**
             * Member function used to read in data and initialize the necessary data
             * to compute the coefficients.
             * 
             * \deprecated
             */
            void initialize (ParameterHandler &param);
            
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
                return typeid(MicroPorousLayer<dim>);
            }
            //@}
            
            ///@name Effective property calculators
            //@{

            /**
             * Compute the effective property in the pores of the MPL. This is used for example to
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
             * Compute the effective property in the pores of the MPL. This is used for example to
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
             * Return the effective diffusivity [\p m^2/s] for nonisothermal with/without two-phase case in the MPL. It takes bulk diffusivity, computed using
             * #compute_gas_diffusion method and transforms it into an effective property,
             * taking into account the porosity, saturation and MPL structure (Anisotropic case),
             * at all quadrature points of the cell.
             */
            virtual void effective_gas_diffusivity(std::vector< Tensor<2,dim> >& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Return the derivative of effective diffusivity w.r.t solution variables/design parameters
             * for nonisothermal with/without two-phase case in the MPL. It transforms bulk diffusion properties computed using
             * #compute_gas_diffusion method and transforms it into an effective property,
             * taking into account the porosity, saturation and MPL structure (Anisotropic case),
             * at all quadrature points of the cell.
             */
            virtual void derivative_effective_gas_diffusivity(std::map< VariableNames, std::vector< Tensor<2,dim> > >& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
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
            virtual void effective_gas_diffusivity(Table< 2, Tensor< 2, dim > > & ) const
            {
                 const std::type_info& info = typeid(*this);
                 FcstUtilities::log << "Pure function " << __FUNCTION__
                 << " called in Class "
                 << info.name() << std::endl;
            };
        
             /**
              * Compute the effective conductivity. To compute the effective properties it will
              * use the method specified in set_method_effective_transport_property_solid . 
              * This function is used for the isotropic MPL
              */
            virtual void effective_electron_conductivity(double& ) const
            {
                 const std::type_info& info = typeid(*this);
                 FcstUtilities::log << "Pure function " << __FUNCTION__
                 << " called in Class "
                 << info.name() << std::endl; 
            };
            
            /**
             * Compute the effective conductivity. To compute the effective properties it will
             * use the method specified in set_method_effective_transport_property_solid . 
             * This function is used for the anisotropic MPL.
             */
            virtual void effective_electron_conductivity(Tensor<2,dim>& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            /**
             * Compute the effective thermal conductivity. To compute the effective properties it will
             * use the method specified in set_method_effective_transport_property_solid . 
             * This function is used for the isotropic MPL
             */
            virtual void effective_thermal_conductivity(double& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the effective thermal conductivity. To compute the effective properties it will
             * use the method specified in set_method_effective_transport_property_solid . 
             * This function is used for the anisotropic MPL.
             */
            virtual void effective_thermal_conductivity(Tensor<2,dim>& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the effective thermal conductivity. It will compute the thermal 
             * conductivity at all quadrature points. 
             * This function is used for the anisotropic MPL.
             */
            virtual void effective_thermal_conductivity(std::vector< Tensor<2,dim> >& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the derivative of effective thermal conductivity with respect to temperature.
             * It will compute the derivative at all quadrature points.
             * This function is used for the anisotropic MPL.
             */
            virtual void derivative_effective_thermal_conductivity(std::vector< Tensor<2,dim> >& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            /**
             * Compute the anisotropic MPL liquid permeability \f$ \left[ cm^2 \right] \f$, at all quadrature points in the cell.
             */
            virtual void liquid_permeablity(std::vector< Tensor<2,dim> >& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the derivative of the anisotropic liquid permeability in the MPL
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
            virtual void relative_liquid_permeability_PSD(std::vector< Tensor<2,dim> >& ) const
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
            
            virtual void derivative_relative_liquid_permeablity_PSD(std::vector<double> & ) const
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
             * Compute \f$ \frac{\partial p_c}{\partial s} \quad \left[ dyne \cdot cm^{-2}\right] \f$, at all quadrature points in the MPL.
             */
            virtual void dpcapillary_dsat(std::vector<double> &) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the derivative of \f$ \frac{\partial p_c}{\partial s} \quad \left[ dyne \cdot cm^{-2}\right] \f$ in the MPL, with 
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
             * quadrature points in the MPL.
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
             * solution variables or design parameters, at all quadrature points in the MPL. The parameters with respect to which 
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
             * Specify the methodology to be used to compute the effective properties 
             * for the porous phase.
             * Return true if the change has been succssful. The different methods are
             * explained in declare_parameters
             */
            virtual bool set_method_effective_transport_property_solid(std::string)
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
                
                return false;
            };
            /**
             * Compute the effective property of a property that is defined by the
             * network of fibres. For example this could be used to compute the effective
             * electron conductivity or heat conduction.
             * NOTE: Isotropic case
             */
            virtual void effective_transport_property_solid(const double& property, 
                                                            double& effective_property) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
                
            };
            
            /**
             * Compute the effective property of a property that is defined by the
             * network of fibres. For example this could be used to compute the effective
             * electron conductivity or heat conduction.
             * Note: Anisotropic case.
             */
            virtual void effective_transport_property_solid(const Tensor<2,dim>& property, 
                                                            Tensor<2,dim>& effective_property) const
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
            ///@name Declaration
            //@{    
            /**
             * Replica Constructor
             * 
             * \warning For internal use only.
             * 
             * Constructor used only to create a prototype. Do not use
             * in general since this will not include the name of the section
             * in the parameter file you need. 
             */
            MicroPorousLayer();
            
            /**
             * Constructor
             * 
             * \note Eventually, I would like to make this private.
             * 
             * \deprecated Use create_MicroPorousLayer
             */
            MicroPorousLayer(const std::string& name);
            
            /**
             * Destructor
             */
            ~MicroPorousLayer();
                       
            /**
             * Declare parameters for a parameter file. The parameters that need to be declared are
             * - Porosity (default : 0.3) Represents the porosity in t*he GDL
             * - Method effective transport properties in pores (default: Bruggemann) Other options "Given|Bruggemann|Percolation|Mezedur"
             * - Method effective transport properties in fibres (default: Bruggemann) Other options "Given|Bruggemann|Percolation|Mezedur"
             */
            void declare_parameters (const std::string& name, 
                                     ParameterHandler &param) const;
                                     

            //@}
            
            ///@name Instance Delivery
            //@{         
            /** 
             * This object is used to store all objects of type MicroPorousLayer. 
             * This information in then used in layer_generator.h in order to create the
             * correct object depending on the specified concrete type of layer selected
             * such as DummyMPL or SGL24BC.
             */
            typedef std::map< std::string, MicroPorousLayer<dim>* > _mapFactory;      
            
            /**
             * Return the map library that stores all childrens of this class. The declare_parameters of
             * each one of the children that are in the map are called in declare_all_MicroPorousLayer_parameters.
             * 
             * @warning In order for children of this class to appear in the map the following four things are 
             * necessary
             * - a static PROTOTYPE object has to be created. For example, 
             * in the .h file:
             * @code
             * static DummyMPL<dim> const* PROTOTYPE;
             * @endcode
             * in the .cc file:
             * @code
             * template <int dim>
             * NAME::DummyMPL<dim> const* NAME::DummyMPL<dim>::PROTOTYPE = new NAME::DummyMPL<dim>();
             * @endcode
             * - a default constructor which creates the PROTOTYPE is needed, e.g.
             * @code
             * template <int dim>
             * NAME::DummyMPL<dim>::DummyMPL()
             * : FuelCellShop::Layer::MicroPorousLayer<dim> ()
             * {
             *    FcstUtilities::log<<" Register DummyMPL MPL to FactoryMap"<<std::endl;
             *    this->get_mapFactory()->insert(std::pair<std::string, FuelCellShop::Layer::MicroPorousLayer<dim>* >(concrete_name, this));
             * }
             * @endcode
             * - a static std::string concrete_name must be declared and initialized to the desired name of
             * the children class. This name is used as the name in the map and as the name of the subsection 
             * where its parameters are declared. For example, 
             * in the .h file:
             * @code
             * static const std::string concrete_name;
             * @endcode
             * in the .cc file:
             * @code
             * template <int dim>
             * const std::string NAME::DummyMPL<dim>::concrete_name ("DummyMPL");
             * @endcode
             * - virtual boost::shared_ptr<FuelCellShop::Layer::MicroPorousLayer<dim> > create_replica () needs to
             * be re-implemented in the child. For example, in the .h file
             * @code
             * virtual boost::shared_ptr<FuelCellShop::Layer::MicroPorousLayer<dim> > create_replica (const std::string &name)
             * {
             *   return boost::shared_ptr<FuelCellShop::Layer::MicroPorousLayer<dim> > (new FuelCellShop::Layer::DummyMPL<dim> (name));
             * }   
             * @endcode
             */
            static _mapFactory * get_mapFactory()
            {
                static _mapFactory mapFactory;
                return &mapFactory;
            }  
            /**
             * This member function is used to create an object of type micro porous layer
             * 
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::Layer::MicroPorousLayer<dim> > create_replica (const std::string &name)
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
             //@}
             ///@name Internal variables
             //@{
             /** Anisotropy variable */
             bool anisotropy;
             /** Oxygen diffusion coefficient */
            Tensor<2,dim> oxygen_diffusivity;
             /** Water diffusion coefficient */
            Tensor<2,dim> water_diffusivity;            
            /** Electrical conductivity from the input file in the anisotripic case */
            Tensor<2,dim> electrical_conductivity;
            /** Thermal conductivity from the input file in the anisotripic case  */
            Tensor<2,dim> thermal_conductivity;
             //@}
        };
        
    }
    
}  // FuelCellShop

#endif // _FUELCELLSHOP__MICRO_POROUS_LAYER_H
