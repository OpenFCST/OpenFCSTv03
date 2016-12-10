//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//    - Class: catalyst_layer.h
//    - Description: Class characterizing the base catalyst layer and defining interface methods.
//    - Developers: Marc Secanell (2011-2013) and Madhur Bhaiya (2013)
//    - Id: $Id: catalyst_layer.h 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__LAYER_CATALYST_LAYER_H
#define _FUELCELLSHOP__LAYER_CATALYST_LAYER_H

// Include deal.II classes

// Include FCST classes
#include <utils/fcst_constants.h>
#include <layers/base_layer.h>
#include <layers/porous_layer.h>
#include <reactions/base_kinetics.h>
#include <materials/polymer_electrolyte_material_base.h>
#include <materials/catalyst_support_base.h>
#include <materials/catalyst_base.h>
#include <materials/platinum.h>
#include <materials/carbon.h>
#include <materials/nafion.h>
#include <materials/PureGas.h>
#include <materials/GasMixture.h>
#include <microscale/PSD_base.h>
#include <application_core/fe_vectors.h>

//Include STL
#include <map>

class MultiScaleCLTest;

namespace FuelCellShop
{
    namespace Layer
    {
        /**
         * Virtual class used to provide the interface for all CatalystLayer children. 
         * 
         * No object of type CatalystLayer should ever be created, instead this layer
         * is used to initialize pointers of type CatalystLayer. The class has a database of
         * children such that it will declare all necessary parameters for all children in the
         * input file, read the input file, create the appripriate children and return a pointer
         * to CatalystLayer with the children selected.
         * 
         * All public functions are virtual but the static functions used to declare parameters and
         * to initialize a pointer of CatalystLayer, i.e. declare_all_CatalystLayer_parameters,
         * set_all_CatalystLayer_parameters and create_CatalystLayer.
         * 
         * <h3>Usage Details:</h3>
         * 
         * In order to create a catalyst layer within an application, the following steps need to be taken.
         * 
         * First, in the application .h file, create a pointer to a CatalystLayer object, i.e.
         * @code
         * boost::shared_ptr<FuelCellShop::Layer::CatalystLayer<dim> > CCL;
         * @endcode
         * 
         * This pointer object will be available anywhere inside the application. Because we do not want to
         * worry about deleting the pointer afterwards, we use a Boost pointer which has its own memory management
         * algorithms. See the <a href="http://www.boost.org/doc/libs/1_54_0/libs/smart_ptr/smart_ptr.htm">Boost website</a> for more information
         * 
         * Once the pointer is available, we need to do three things in the application
         * - Call declare_CatalystLayer_parameters in the application declare_parameters() member function. This 
         * member function is used to define all parameters that can be read from the input file
         * 
         * 
         * - Call create_CatalystLayer and initialize. The former member function will fill the pointer created above with the appropriate
         * catalyst layer you have selected in the input file. Then, initialize reads from the input file all the data from the file in order
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
         *   // Declare section on the input file where all info will be stored. In this case Fuel Cell Data > Cathode catalyst layer
         *   FuelCellShop::Layer::CatalystLayer<dim>::declare_CatalystLayer_parameters("Cathode catalyst layer", param);
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
         *    catalyst.initialize(param);
         *    electrolyte.initialize(param);
         *    catalyst_support.initialize(param);  
         *     
         *    // Create a catalyst layer. When you create the layer, you also specify the type of electrolyte, catalyst support and catalyst in the layer. For example,
         *    // a conventional catalyst layer will take a FuelCellShop::Material::Nafion type electrolyte, a FuelCellShop::Material::Carbon support
         *    // and a FuelCellShop::Material::Platinum catalyst.
         *    CCL = FuelCellShop::Layer::CatalystLayer<dim>::create_CatalystLayer("Cathode catalyst layer", param);
         * 
         *    // Here, I specify the gases that exist in the CCL and their temperature and pressure (based on operating conditions):
         *    CCL->set_gases_and_compute(gases, OC.get_pc_atm (), OC.get_T());
         *    // Initialise the necessary kinetics parameters in CCL.
         *   (...)
         * }
         * @endcode
         */
        template <int dim>
        class CatalystLayer :
        public PorousLayer<dim>
        {
            
        public:
            ///@name Friend class for Unit Testing
            //@{                     
            
            /**
             * Friend class for testing purposes.
             */
            friend class ::MultiScaleCLTest;
            //@}
            
            ///@name Instance Delivery (Functions)
            //@{    

            /**
             * Function used to declare all the data necessary in the parameter files former
             * all CatalystLayer children.
             * 
             * This member function should be used instead of declare_parameters() when we want
             * to use a CatalystLayer pointer that selects the type of CL to run at runtime.
             * 
             * \param cl_section_name Name of the section that will encapuslate all the information about the CL
             * \param param ParameterHandler object used to store all information about the simulation. Used
             * to read the parameter file.
             * 
             * The parameter file would look as follows:
             * 
             * subsection cl_section_name
             *   set Catalyst layer type = DummyCL        # Options: DummyCL | HomogeneousCL | MultiScaleCL
             *   set Catalyst type = Platinum             # Options: Platinum
             *   set Catalyst support type = CarbonBlack  # Options: CarbonBlack
             *   set Electrolyte type  = Nafion           # Options: Nafion
             * 
             */
            static void declare_CatalystLayer_parameters (const std::string &cl_section_name, 
                                                          ParameterHandler &param)
            {
                for (typename FuelCellShop::Layer::CatalystLayer<dim>::_mapFactory::iterator iterator = FuelCellShop::Layer::CatalystLayer<dim>::get_mapFactory()->begin(); 
                     iterator != FuelCellShop::Layer::CatalystLayer<dim>::get_mapFactory()->end(); 
                iterator++)
                     {
                         iterator->second->declare_parameters(cl_section_name, param);
                     }   
            }
            
           /**
            * 
            * Function used to select the appropriate CatalystLayer type as specified in the ParameterHandler under
            * line 
            * @code 
            * set Catalyst layer type = DummyCL 
            * @endcode
            * current options are [ DummyCL | MultiScaleCL | HomogeneousCL ]
            * 
            * The class will read the appropriate section in the parameter file, i.e. the one with name \param cl_section_name ,
            * create an object of the desired type and return it.
            * 
            */
            static boost::shared_ptr<FuelCellShop::Layer::CatalystLayer<dim> > create_CatalystLayer (const std::string& cl_section_name, 
                                                                                                     ParameterHandler &param)
            {
                boost::shared_ptr<FuelCellShop::Layer::CatalystLayer<dim> > pointer;
                
                std::string concrete_name;
                param.enter_subsection("Fuel cell data");
                {
                    param.enter_subsection(cl_section_name);
                    {
                        concrete_name = param.get("Catalyst layer type");
                        //FcstUtilities::log << "name: "<<concrete_name.c_str()<<std::endl;
                    }
                    param.leave_subsection();
                }
                param.leave_subsection();
                
                typename FuelCellShop::Layer::CatalystLayer<dim>::_mapFactory::iterator iterator = FuelCellShop::Layer::CatalystLayer<dim>::get_mapFactory()->find(concrete_name);
                
                if (iterator != FuelCellShop::Layer::CatalystLayer<dim>::get_mapFactory()->end())
                {
                    if (iterator->second)
                    {
                        pointer = iterator->second->create_replica(cl_section_name);
                    }
                    else 
                    {
                        FcstUtilities::log<<"Pointer not initialized"<<std::endl;
                        abort();    
                    }
                }
                else
                {
                    FcstUtilities::log<<"Concrete name in FuelCellShop::Layer::CatalystLayer<dim>::create_CatalystLayer does not exist"<<std::endl;
                    abort();
                }
                
                pointer->initialize(param);
                
                return pointer;
            }
            //@}
            
            ///@name Initalization
            //@{    
            /**
             * Set those solution variables which are constant in the particular application. If the effective properties in the layer depend on other variables that are usually
             * part of the solution vector but are assumed to be constant in this simulation, the const solution value should be passed to the class using this member function. This 
             * method should be called in the initialization section of the application. This function takes value to be set as the first argument and the #VariableNames as 
             * second argument. For instance, it's required to store constant temperature value for an isothermal application, in that case this method can be used. \em e.g., in 
             * order to set temperature as \p 353.0 [\p Kelvin] in the layer, you can use the following code:
             * @code
             * // In the initialization section of the application.
             * layer.set_constant_solution(353.0, VariableNames::temperature_of_REV);
             * @endcode
             * 
             * If #temperature_of_REV is passed using this method, it also sets the temperature [\p Kelvin] values in the #electrolyte object. If #total_pressure is passed 
             * using this method, it also sets the total pressure [\p Pascals] values in the #kinetics and #electrolyte object.
             */
            virtual void set_constant_solution(const double& value, const VariableNames& name)
            {
                FuelCellShop::Layer::BaseLayer<dim>::set_constant_solution(value, name);
                
                if (name == temperature_of_REV)
                {
                    this->electrolyte->set_T(value);
                }
                
                else if (name == total_pressure)
                {
                    this->kinetics->set_p_t(value);
                    this->electrolyte->set_p_t(value);
                }
                else if (name == membrane_water_content)
                {
                    this->electrolyte->set_lambda(value);
                }
            }
            
            /**
             * This method is used to set the solution variable values in the #kinetics object, at all quadrature points in the cell. It takes
             * vector of SolutionVariable structures as input argument. Each one of them corresponds to a particular solution variable, required
             * in order to compute various terms such as non-linear current source terms etc.
             * 
             * The variables that must be set are:
             * - electronic_electrical_potential
             * - protonic_electrical_potential
             * 
             * For a cathode electrode the concentration of oxygen in the gas phase is also necessary. Note that this member function will convert the
             * gas concentration to electrolyte concentration already. For convenience, oxygen molar fractions can be passed to this class.
             * 
             * For an anode electrode the concentration of hydrogen in the gas phase is needed. Again, this member function will convert the
             * gas concentration to electrolyte concentration. For convenience, hydrogen mole fractions can be pass to this class.
             * 
             * \note Use only for solution variables.
             */
            virtual void set_solution(const std::vector< SolutionVariable >&);
            
            /**
             * Method used to set the variables for which you would like to compute the derivatives in the catalyst layer. It takes vector of
             * #VariableNames as an input argument. It also sets the derivative flags in the kinetics and electrolyte object of the catalyst layer.
             */
            virtual void set_derivative_flags(const std::vector<VariableNames>& flags)
            {
                this->derivative_flags = flags;
                this->kinetics->set_derivative_flags(flags);
                this->electrolyte->set_derivative_flags(flags);
            }
            
            /**
             * Member function used to specify the reaction for which the kinetic parameters are needed, for example
             * for a Platinum catalyst, we can specify that we need the kinetic parameters for either the oxygen reduction reaction (ORR)
             * or the hydrogen oxidation reaction (HOR)
             */
            void set_reaction_kinetics(const ReactionNames rxn_name)
            {
                this->kinetics->set_reaction_kinetics(rxn_name);
                this->catalyst->set_reaction_kinetics(rxn_name);
            }
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
                return typeid(CatalystLayer<dim>);
            }

            /**
             * Compute the volume fractions of each phase
             * 
             * The map might contains a string indicating the phase and its value. There are
             * several possible phases:
             * - Solid
             * - Void
             * - Ionomer
             * - Water
             */
            virtual void get_volume_fractions(std::map<std::string, double>& )
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            /**
             * Return loadings
             * - \param V_Pt = Pt loading in ug/cm3
             * - \param loading_N = ionomer loading %wt
             * - \param IC_ratio = I/C ratio
             * @note either loading_N or IC_ratio should be passed through input file, the other computed
             * - \param prc_Pt = Pt/C ratio           
             */
            virtual inline void get_loadings(std::map<std::string, double>& )
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl; 
            };
            
            //@}
            
            /**
             * Function for setting current cell_id from applications.
             *
             * - \param unsigned int id is the id of the current cell from the application's perspective
             *
             * @note Reimplemented by classes that need to know cell_id: MultiScaleCL
             */
            virtual void set_cell_id(const unsigned int& id){
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                        << " called in Class "
                        << info.name() << std::endl;
            }

            ///@name Effective property calculators
            //@{

            /**
             * Compute the effective property in the pores of the CL. This is used for example to
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
             * Return the effective diffusivity [\p m^2/s] for nonisothermal with/without two-phase case in the CL. It takes bulk diffusivity, computed using
             * #compute_gas_diffusion method and transforms it into an effective property,
             * taking into account the porosity, saturation and CL structure (Anisotropic case),
             * at all quadrature points of the cell.
             * 
             * \note: For two-phase case, set_saturation should be called before using this method, otherwise this method assumes saturation value to be zero.
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
             * for nonisothermal with/without two-phase case in the CL. It transforms bulk diffusion properties computed using
             * #compute_gas_diffusion method and transforms it into an effective property,
             * taking into account the porosity, saturation and CL structure (Anisotropic case),
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
             * Compute the effective property in the pores of the CL. This is used to compute effective diffusivity of gases.
             * This routine can be used either in the isotropic or anisotropic cases.
             * Bulk diffusion coefficients or their derivatives are obtained from Mixure::BinaryDiffusion classes
             * inside this method.
             * 
             * \note The routine FuelCellShop::Layer::PorousLayer< dim >::set_gases_and_compute (std::vector< FuelCellShop::Material::PureGas * > &gases, double pressure, double temperature)
             * (in the parent class) should have been called prior to using this class. This method is to be used only for a single-phase, isothermal application.
             */
            virtual void effective_gas_diffusivity(Table< 2, Tensor<2,dim> >&) const
            {
                 const std::type_info& info = typeid(*this);
                 FcstUtilities::log << "Pure function " << __FUNCTION__
                 << " called in Class "
                 << info.name() << std::endl;
            };
            
            /**
             * Compute the effective electron conductivity in the CL
             */
            virtual void effective_electron_conductivity(double& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the effective electron conductivity in the CL
             */
            virtual void effective_electron_conductivity(Tensor<2,dim>& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the derivative of the effective electron conductivity in the GDL
             * with respect to either the solution or design parameters. The parameters with respect to
             * which the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_effective_electron_conductivity(std::vector<double>& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
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
             * Compute the effective proton conductivity in the CL
             */
            virtual void effective_proton_conductivity(double& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the effective proton conductivity, at all quadrature points in the cell, mainly as a function of Temperature.
             */
            virtual void effective_proton_conductivity(std::vector<double>& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the derivative of the effective proton conductivity in the CL
             * with respect to either the solution or design parameters. The parameters with respect to
             * which the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_effective_proton_conductivity(std::map< VariableNames, std::vector<double> >& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            /**
             * Compute the effective thermal conductivity in the CL
             */
            virtual void effective_thermal_conductivity(double& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the effective thermal conductivity in the CL
             */
            virtual void effective_thermal_conductivity(Tensor<2,dim>& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the effective thermal conductivity as a Tensor at all quadrature points
             */
            virtual void effective_thermal_conductivity(std::vector< Tensor<2,dim> >& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;	
            }
            /**
             * Compute the derivative of the effective thermal conductivity in the CL
             * with respect to either the solution or design parameters. The parameters with respect to
             * which the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_effective_thermal_conductivity(std::vector<double>& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the derivative of the effective thermal conductivity in the CL
             * with respect to either the solution or design parameters. The parameters with respect to
             * which the derivatives are computed are setup in \ref FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_effective_thermal_conductivity(std::vector<Tensor<2,dim> >& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            /**
             * Compute the effective water diffusivity (lambda diffusivity) in the CL.
             */
            virtual void effective_water_diffusivity(double& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the effective water diffusivity (lambda diffusivity) at all
             * quadrature points in the CL.
             */
            virtual void effective_water_diffusivity(std::vector<double>& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            
            /**
             * Compute the derivative of the effective water diffusivity (lambda diffusivity) in the CL
             * with respect to either the solution or design parameters. The parameters with respect to
             * which the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_effective_water_diffusivity(std::map< VariableNames, std::vector<double> >& ) const
            {
               const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            
            /**
             * Compute the effective thermo-osmotic diffusivity of lambda (sorbed water),
             * at all quadrature points in the CL.
             */
            virtual void effective_thermoosmotic_diffusivity(std::vector<double>& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            /**
             * Compute the derivative of the effective thermo-osmotic diffusivity of lambda (sorbed water) in the CL
             * with respect to either the solution or design parameters. The parameters with respect to
             * which the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_effective_thermoosmotic_diffusivity(std::map< VariableNames, std::vector<double> >& ) const
            {
               const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            
            /**
             * Compute the CL gas permeability
             */
            virtual void gas_permeablity(double& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the CL gas permeability
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
             * which the derivatives are computed are setup in \ref FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_gas_permeablity(std::vector<double>& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the derivative of the effective gas permeability in the GDL
             * with respect to either the solution or design parameters. The parameters with respect to
             * which the derivatives are computed are setup in \ref FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_gas_permeablity(std::vector<Tensor<2,dim> >& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            /**
             * Compute the anisotropic CL liquid permeability \f$ \left[ cm^2 \right] \f$, at all quadrature points in the cell.
             */
            virtual void liquid_permeablity(std::vector< Tensor<2,dim> >& ) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the derivative of the anisotropic liquid permeability in the CL
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
             * Compute the anisotropic CL liquid permeability \f$ \left[ cm^2 \right] \f$, at all quadrature points in the cell.
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
             * Compute \f$ \frac{\partial p_c}{\partial s} \quad \left[ dyne \cdot cm^{-2}\right] \f$, at all quadrature points in the CL.
             */
            virtual void dpcapillary_dsat(std::vector<double> &) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            /**
             * Compute the derivative of \f$ \frac{\partial p_c}{\partial s} \quad \left[ dyne \cdot cm^{-2}\right] \f$ in the CL, with 
             * respect to either the solution or design parameters, at all quadrature points in the cell. The parameters with respect to which 
             * the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags().
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
             * quadrature points in the CL.
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
             * solution variables or design parameters, at all quadrature points in the CL. The parameters with respect to which 
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
            
          ///@name Reaction terms
          //@{            
            /**
             * This member function will use a FuelCellShop::BaseKinetics class in order to compute the current density
             * production in the CL
             */
            virtual void current_density(std::vector<double>&)
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            /**
             * This member function will use a FuelCellShop::BaseKinetics class in order to compute the current density
             * production in the CL. First argument is <b>current density</b>, and second is <b>effectiveness</b> at all quadrature points in the cell.
             */
            virtual void current_density(std::vector<double>&current, std::vector<double>&effectiveness)
            {
                for (unsigned int i = 0; i<effectiveness.size(); ++i)
                    effectiveness[i] = 0;
                
                this->current_density(current);
            };
            
            /**
             * This member function will use a FuelCellShop::BaseKinetics class in order to compute the derivative of the current density
             * with respect to the variables setup using #set_derivative_flags
             */
            virtual void derivative_current_density(std::map< VariableNames, std::vector<double> >& )
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            };
            
            /** Get the active area of platinum per unit volume of CL */
            virtual double get_active_area_Pt() const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
                return 0.0;
            };
            
            /** Method to provide access to pointer of the electrolyte object of the catalyst layer. */
            virtual inline FuelCellShop::Material::PolymerElectrolyteBase* get_electrolyte() const
            {   
                return this->electrolyte.get();
            }
            
            
            /** Method to provide access to pointer of the kinetic object of the catalyst layer. */
            virtual inline FuelCellShop::Kinetics::BaseKinetics* get_kinetics() const
            {
                return this->kinetics.get();
            }
            
            //@}
 
            /**Method for getting string describing kinetics type (corresponding to kinetics class concrete names)*/
            inline std::string get_kinetics_type(){
                return kinetics_type;
            }


            /**Method for getting coverages from kinetics objects (overloaded by MultiScaleCL)*/
            virtual SolutionMap get_coverages();




        protected:
            
            ///@name Constructors, destructor, and initalization
            //@{    
            /** \warning For internal use only.
             * 
             * Constructor used only to create a prototype. Do not use
             * in general since this will not include the name of the section
             * in the parameter file you need. 
             */
            CatalystLayer();
            
            /**
             * Destructor
             */
            ~CatalystLayer();
            
            /**
             * Constructor
             */
            CatalystLayer(const std::string& name);
                       
            
            /**
             * Default virtual declare parameters for a parameter file. 
             * 
             * \note This member function must be virtual since it will be accessed
             * via pointers for all children.
             * 
             */
            virtual void declare_parameters (const std::string& name, ParameterHandler &param) const;
                                         
            /**
             * Member function used to read in data and initialize the necessary data
             * to compute the coefficients.
             */
            void initialize (ParameterHandler &param);
            
  
            //@}
            ///@name Instance Delivery (Types)
            //@{         
            /** 
             * This object is used to store all objects of type CatalystLayer. 
             * This information in then used in layer_generator.h in order to create the
             * correct object depending on the specified concrete type of layer selected
             * such as DummyCL.
             */
            typedef std::map< std::string, CatalystLayer<dim>* > _mapFactory;      
            //@}

            ///@name Instance Delivery (Function)
            //@{         
             /**
             * Return the map library that stores all childrens of this class. The declare_parameters of
             * each one of the children that are in the map are called in declare_all_CatalystLayers.
             * 
             * @warning In order for children of this class to appear in the map the following four things are 
             * necessary
             * - a static PROTOTYPE object has to be created. For example, 
             * in the .h file:
             * @code
             * static DummyCL<dim> const* PROTOTYPE;
             * @endcode
             * in the .cc file:
             * @code
             * template <int dim>
             * NAME::DummyCL<dim> const* NAME::DummyCL<dim>::PROTOTYPE = new NAME::DummyCL<dim>();
             * @endcode
             * - a default constructor which creates the PROTOTYPE is needed, e.g.
             * @code
             * template <int dim>
             * NAME::DummyCL<dim>::DummyCL()
             * : FuelCellShop::Layer::CatalystLayer<dim> ()
             * {
             *    FcstUtilities::log<<" Register DummyGDL GDL to FactoryMap"<<std::endl;
             *    this->get_mapFactory()->insert(std::pair<std::string, FuelCellShop::Layer::GasDiffusionLayer<dim>* >(concrete_name, this));
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
             * const std::string NAME::DummyCL<dim>::concrete_name ("DummyCL");
             * @endcode
             * - virtual boost::shared_ptr<FuelCellShop::Layer::CatalystLayer<dim> > create_replica () needs to
             * be re-implemented in the child. For example, in the .h file
             * @code
             * virtual boost::shared_ptr<FuelCellShop::Layer::CatalystLayer<dim> > create_replica (std::string &name, 
             *                                                                                     FuelCellShop::Material::PolymerElectrolyteBase* electrolyte1,
             *                                                                                     FuelCellShop::Material::CatalystSupportBase* catalyst_support1,
             *                                                                                     FuelCellShop::Material::CatalystBase* catalyst1)
             * {
             *   return boost::shared_ptr<FuelCellShop::Layer::CatalystLayer<dim> > (new FuelCellShop::Layer::DummyCL<dim> (name, electrolyte1, catalyst_support1, catalyst1));
             * }   
             * @endcode
             */
            static _mapFactory * get_mapFactory()
            {
                static _mapFactory mapFactory;
                return &mapFactory;
            }                               
            ///@name Instance Delivery (Private functions)
            //@{                     
            /**
             * This member function is used to create an object of type gas diffusion layer
             * 
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::Layer::CatalystLayer<dim> > create_replica (const std::string &name)
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            //@}
            
            ///@name Internal variables
            //@{                     
            /** If CL properties are stored inside the class (e.g. DummyCL) then, return the property
             * stored under coefficient_name name
             */
            std::string diffusion_species_name;
            
            /** If the default materials are used in the layer, this will be set
             * to true. If the materials are created in the application and passed down
             * this will be false. Used for the destructor.*/
            bool default_materials;
            
            /** Catalyst type from input file */
            std::string catalyst_type;
            
            /** Catalyst Support type from input file */
            std::string catalyst_support_type;

            /** Electrolyte type from input file */
            std::string electrolyte_type;

            /** Kinetic class type from input file */
            std::string kinetics_type;
            
            /** PSD class type from input file */
            std::string PSD_type;
            
            /**
             * Pointer to the electrolyte object created in the application that is used 
             * to calculate the properties of the electrolyte in the catalyst layer.
             */
            boost::shared_ptr< FuelCellShop::Material::PolymerElectrolyteBase > electrolyte;
            
            /**
             * Pointer to the catalyst support object created in the application that 
             * is used to calculate the carbon black conductivity in the catalyst layer.
             */
            boost::shared_ptr< FuelCellShop::Material::CatalystSupportBase > catalyst_support;

            /**
             * Pointer to the catalyst object created in the application that is used 
             * to store the properties of the catalyst used in the layer.
             */
            boost::shared_ptr< FuelCellShop::Material::CatalystBase > catalyst;
            
            /** Pointer to a kinetics object */
            boost::shared_ptr< FuelCellShop::Kinetics::BaseKinetics > kinetics;
            
            /**
             * Stores the number of quadrature points in the cell.
             */
            unsigned int n_quad;
            
            /** Map storing solution variables. \p Key represents the name of the variable and \p Value represents the FuelCellShop::SolutionVariable structure. */
            std::map<VariableNames ,SolutionVariable> solutions;
            
            /** Name of the reactant which is being solved for in the catalyst layer. This is specifically used by the FuelCellShop::Layer::MultiScaleCL<dim> object. */
            VariableNames reactant;
            
            #ifdef DEBUG

            /* A list of common variable names that are set to the layer
             *
             * When running in debug this list will be compared with provided list of
             * solutions, in order to check application is consistent in providing solutions.
             *
             * Inconsistent provision of solutions can lead to invalid results.
             */
            std::vector<VariableNames> common_names;

            #endif


            //@}
        };
        
    } // Layer
    
}  // FuelCellShop
    
#endif // _FUELCELLSHOP__LAYER__CATALYST_LAYER_H
    
