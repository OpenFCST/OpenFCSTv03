//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: nafion.h
//    - Description:  Class to represent and return bulk properties for Nafion (Polymer Electrolyte Membrane)
//    - Developers: M. Secanell (2009-13) and Madhur Bhaiya (2012-13)
//    - Id: $Id: nafion.h 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------


#ifndef _FUELCELLSHOP__NAFION_H
#define _FUELCELLSHOP__NAFION_H

// Include deal.II classes
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/fe/fe_values.h>

#include <materials/polymer_electrolyte_material_base.h>
#include <materials/PureGas.h>

class NafionTest;

namespace FuelCellShop
{
    namespace Material
    {
        /**
         * Class used to store data and provide information regarding the electrolyte. In this case
         * the electrolyte that is setup is Nafion. In particular the data coded in this class is obtained
         * from the literature from Nafion 1100, usually from NRE 211 properties. 
         * 
         * In the parameter file the following parameters can be specified within the subsections specified below:
         * 
         * @code
         * subsection Fuel cell data
         *   subsection Materials
         *     subsection Nafion
         *       set Method for sorption isotherm = Hinatsu        # Options are Hinatsu | Liu09 -- Method to compute equilibrium lambda from sorption isotherm - Depends on Water molar fraction and Temperature
         *       set Method to compute proton conductivity = NRE211 # Options are: Constant|Springer|NRE211|Iden11 -- Method used to compute proton conductivity inside the membrane - Depends on lambda and T
         *       set Proton conductivity [S/cm] = 0.1                      # Proton conductivity inside the membrane, S/cm
         *       set Method to compute water diffusion = Springer   # Options are: "Constant|Springer|Motupally" -- Method used to compute diffusion of water inside the membrane - Depends on lambda and T
         *       subsection Springer coefficients --- This section is used to control the coefficients if Springer is used which are of the form: (A \lambda + B) exp [ C (1/303 - 1/T)
         *          set Proton conductivity's first coefficient = 0.005139
         *          set Proton conductivity's second coefficient = -0.00326
         *          set Proton conductivity's third coefficient = 1268.0
         *       end
         *       set Water diffusion coefficient = 2e-              # Water diffusion coefficient inside the membrane
         *       set Electro-osmotic drag method = Springer         # Options are: Constant|Springer -- Method to compute nd. Springer depends on water content
         *       set Electro-osmotic drag coefficient = 1.0         # Number of water molecules dragged by one proton
         *       set Method to compute thermo-osmotic diffusion coefficient = Kim09 # Options are: Constant | Kim09 -- Method to compute thermo-osmotic diffusion coefficient - Depends on T
         *       set Thermo-osmotic diffusion coefficient [gm/(cm-s-K)] = -1.3e-7 # Thermo-osmotic diffusion coefficient inside the membrane [gm/(cm-s-K)].
         *       set Method to compute enthalpy of sorption of water = Constant  # Options are: Constant -- Method to compute enthalpy of sorption of water - may depend on lambda and T
         *       set Enthalpy of sorption of water [J/mol] = 45000.0     # Heat released when one mole of water is sorbed inside the membrane.
         *       set Oxygen diffusion coefficient [cm^2/s] = 9.726e-6    # From J. Peron et al., “Properties of Nafion NR-211 membranes for PEMFCs”
         *       set Proton diffusion coefficient [cm^2/s] = 9.2e-5      # Peter Dobsons thesis page 42
         *       set Henry's Law Constant for Oxygen [Pa cm^3/mol] = 3.1664e10
         *       set Henry's Law Constant for Hydrogen [Pa cm^3/mol] = 6.69e10
         *     end
         *   end
         * end
         * 
         * @endcode
         * 
         * <h3>Usage Details:</h3>         
         * 
         * As with most routines you need to first declare_parameters, and then initialize. After this the class is ready for use.
         * 
         * See for example FuelCell::OperatingConditions class for more details.
         * 
         * \authors M. Secanell and M. Bhaiya
         * 
         */
        class Nafion :
        public PolymerElectrolyteBase
        {
        public:
            
            /** 
             * Name of the class. This name is used to select the layer.
             */
            static const std::string concrete_name;
            
            ///@name Friend class for Unit Testing
            //@{                     
            
            /**
             * Friend class for testing purposes.
             */
            friend class ::NafionTest;
            //@}
            
            ///@name Constructors, destructor, and parameter initalization
            //@{
            /** PROTOTYE Constructor 
             * 
             * \warning For internal use only.
             */
            Nafion(const bool);
            
            /** Constructor 
             * The constructor initialize parameters using the default values. This is
             * so that if I do not want to call declare_parameters and initialize, I can
             * still use the routine with the hard coded values.
             * 
             * \todo1 Make this private and remove name.
             */
            Nafion(std::string name = "Nafion");
            
            /**
             * Destructor
             */
            ~Nafion();
            
            /** Declare parameters.
             * 
             * \todo1 Make it protected
             * 
             * In the parameter file the following parameters can be specified within the subsections specified below:
             * 
             * @code
             * subsection Fuel cell data
             *   subsection Materials
             *     subsection Nafion
             *       set Method for sorption isotherm = Previous        # Options are Previous | Liu09 -- Method to compute equilibrium lambda from sorption isotherm - Depends on Water molar fraction and Temperature
             *       set Method to compute proton conductivity = NRE211 # Options are: Constant|Springer|NRE211|Iden11 -- Method used to compute proton conductivity inside the membrane - Depends on lambda and T
             *       set Proton conductivity [S/cm] = 0.1                      # Proton conductivity inside the membrane, S/cm
             *       set Method to compute water diffusion = Springer   # Options are: "Constant|Springer|Motupally" -- Method used to compute diffusion of water inside the membrane - Depends on lambda and T
             *       set Water diffusion coefficient = 2e-              # Water diffusion coefficient inside the membrane
             *       set Electro-osmotic drag method = Springer         # Options are: Constant|Springer -- Method to compute nd. Springer depends on water content
             *       set Electro-osmotic drag coefficient = 1.0         # Number of water molecules dragged by one proton
             *       set Method to compute thermo-osmotic diffusion coefficient = Kim09 # Options are: Constant | Kim09 -- Method to compute thermo-osmotic diffusion coefficient - Depends on T
             *       set Thermo-osmotic diffusion coefficient [gm/(cm-s-K)] = -1.3e-7 # Thermo-osmotic diffusion coefficient inside the membrane [gm/(cm-s-K)].
             *       set Method to compute enthalpy of sorption of water = Constant  # Options are: Constant -- Method to compute enthalpy of sorption of water - may depend on lambda and T
             *       set Enthalpy of sorption of water [J/mol] = 45000.0     # Heat released when one mole of water is sorbed inside the membrane.
             *       set Oxygen diffusion coefficient [cm^2/s] = 9.726e-6    # From J. Peron et al., “Properties of Nafion NR-211 membranes for PEMFCs”
             *       set Proton diffusion coefficient [cm^2/s] = 9.2e-5      # Peter Dobsons thesis page 42
             *       set Henry's Law Constant for Oxygen [Pa cm^3/mol] = 3.1664e10
             *       set Henry's Law Constant for Hydrogen [Pa cm^3/mol] = 6.69e10
             *     end
             *   end
             * end
             * 
             * @endcode
             */
            virtual void declare_parameters(ParameterHandler &param) const;
            
            /**
             * Member function used to read in data and initialize the necessary data from the parameter file to compute the coefficients.
             */
            virtual void initialize (ParameterHandler &param);
            
            //@}
            
            ///@name Transport properties and derivatives accessor methods
            //@{
            
            /**
             * Compute the equilibrium water content, \f$ \lambda_{eq} \f$, inside the Nafion for vapor-equilibriated case, at
             * every quadrature point in the cell. It takes a vector as an input argument and values are passed by reference.
             * It is required to call #set_water_molar_fraction method atleast before computing sorption isotherm values.
             */
            virtual void sorption_isotherm(std::vector<double>&) const;
            /**
             * Compute the derivatives for water sorption source terms, \f$ \frac{\partial \lambda_{eq}}{\partial u} \f$, at every quadrature 
             * point in the cell. The derivatives are computed based on the flags set by the #set_derivative_flags method. It takes map as an input argument by reference, in 
             * which \p Key corresponds to the variable about which derivative is being computed, \em i.e., #VariableNames and \p Value corresponds to
             * the vector storing the derivative values. It is required to call #set_water_molar_fraction method atleast before computing sorption isotherm values. Also, this 
             * method should only be called when atleast there is one derivative flag.
             */
            virtual void sorption_isotherm_derivative(std::map < VariableNames, std::vector<double> >&) const;
            
            /**
             * Compute the proton conductivity, \f$ \sigma_{H^+} \f$ [\p S/cm], inside the nafion, at every quadrature point in the cell. It is computed \
             * as a function of \f$ \lambda \f$ and \f$ T \f$. It can be called by either constant or variable solution setting methods to avoid duplicate equations.
             * The first parameter is the output bulk proton conductivity and the two following are input parameters such as temperature, \f$ T \f$ [\p K], and water content, \f$ \lambda \f$ [\p -].
             */
            virtual void proton_conductivity(double&, const double, const double) const;
            /**
             * Compute the proton conductivity, \f$ \sigma_{H^+} \f$ [\p S/cm], inside the nafion for a case of constant \f$ \lambda \f$ and 
             * \f$ T \f$. Before calling this method, #set_T method should be used for setting the constant temperature value. If the constant \f$ \lambda \f$ 
             * value is required to be different from the default value \b 12.0, then #set_lambda method should be used. It takes a double as an input argument and value is passed by reference. 
             * \remarks This method is recommended only for multi-scale models. In case of macrohomogeneous models, constant values should always be set before the variable ones,
             * otherwise it may give erroneous results.
             */
            virtual void proton_conductivity(double&) const;
            /**
             * Compute the proton conductivity, \f$ \sigma_{H^+} \f$ [\p S/cm], inside the nafion, at every quadrature point in the cell. It is computed 
             * as a function of \f$ \lambda \f$ and \f$ T \f$, hence solution setting methods should be called accordingly for them, whether constant or variable. It is recommended to set the
             * constant ones in the initialization of the application, unless you are solving for multi-scale model. It takes vector as an input argument and values are passed by reference.
             * \warning Precaution must be used while using two different solution setting methods, \em i.e., either constant or variable for a single solution variable on the same electrolyte object. You 
             * cannot use constant solution setting methods for this object, if variable solution setting methods are already used for a single solution variable. In case 
             * of constant \f$ \lambda \f$ and \f$ T \f$, vector size should be set equal to the number of quadrature points, before passing to this function. 
             */
            virtual void proton_conductivity(std::vector<double>&) const;
            /**
             * Compute the derivatives of proton conductivity at every quadrature point in the cell. The derivatives are computed based on the flags
             * set by the #set_derivative_flags method. It takes map as an input argument by reference, in  which \p Key corresponds to the variable about 
             * which derivative is being computed, \em i.e., #VariableNames and \p Value corresponds to the vector storing the derivative values. Note that 
             * this method should not be used when both \f$ \lambda \f$ and \f$ T \f$ are constant. Hence, atleast #set_membrane_water_content or #set_temperature should be 
             * called before using this method. Also, this method should only be called when atleast there is one derivative flag.
             */
            virtual void proton_conductivity_derivative(std::map< VariableNames, std::vector<double> >&) const;
            
            /**
             * Compute the water diffusivity, \f$ D_{\lambda} \f$ [\p cm^2/s], inside the nafion at every quadrature point in the cell. It is computed 
             * as a function of \f$ \lambda \f$ and \f$ T \f$. It can be called by either constant or variable solution setting methods to avoid duplicate equations.
             * The first parameter is the output water diffusivity and the two following are input parameters such as temperature, \f$ T \f$ [\p K], and water content, \f$ \lambda \f$ [\p -].
             */   
            virtual void water_diffusivity(double&, const double, const double) const;
            /**
             * Compute the water diffusivity, \f$ D_{\lambda} \f$ [\p cm^2/s], inside the nafion for a case of constant \f$ \lambda \f$ and 
             * \f$ T \f$. Before calling this method, #set_T method should be used for setting the constant temperature value. If the constant \f$ \lambda \f$ 
             * value is required to be different from the default value \b 12.0, then #set_lambda method should be used. It takes a double as an input argument and value is passed by reference. 
             * \remarks This method is recommended only for multi-scale models. In case of macrohomogeneous models, it's not required in particular because \f$ D_{\lambda} \f$ is 
             * required only when \f$ \lambda \f$ is one of the solution variables in the application.
             */
            virtual void water_diffusivity(double&) const;
            /**
             * Compute the water diffusivity, \f$ D_{\lambda} \f$ [\p cm^2/s], inside the nafion at every quadrature point in the cell. It is computed 
             * as a function of \f$ \lambda \f$ and \f$ T \f$, hence solution setting methods should be called accordingly for them, whether constant or variable. It is recommended to set the
             * constant ones in the initialization of the application. Note that this method should not be used when both \f$ \lambda \f$ and \f$ T \f$ are constant. Hence, 
             * atleast #set_membrane_water_content or #set_temperature should be called before using this method. It takes vector as an input argument and values are passed by reference.
             */
            virtual void water_diffusivity(std::vector<double>&) const;
            /**
             * Compute the derivatives of water diffusivity at every quadrature point in the cell. The derivatives are computed based on the flags
             * set by the #set_derivative_flags method. It takes map as an input argument by reference, in  which \p Key corresponds to the variable about 
             * which derivative is being computed, \em i.e., #VariableNames and \p Value corresponds to the vector storing the derivative values. Note that 
             * this method should not be used when both \f$ \lambda \f$ and \f$ T \f$ are constant. Hence, atleast #set_membrane_water_content or #set_temperature should be 
             * called before using this method. Also, this method should only be called when atleast there is one derivative flag.
             */
            virtual void water_diffusivity_derivative(std::map< VariableNames, std::vector<double> >&) const;
            
            /**
             * Compute the electro-osmotic drag coefficient inside the nafion at every quadrature point in
             * the cell. It takes vector as an input argument and values are passed by reference. Note that the double return method is not implemented for electro-osmotic drag 
             * coefficients because it is required only when \f$ \lambda \f$ is one of the solution variables in the application. #set_membrane_water_content 
             * should be called before using this method.
             */
            virtual void electroosmotic_drag(std::vector<double>&) const;
            /**
             * Compute the derivatives of electro-osmotic drag coefficient, at every quadrature point in the cell. The derivatives are computed based on the flags
             * set by the #set_derivative_flags method. It takes map as an input argument by reference, in  which \p Key corresponds to the variable about 
             * which derivative is being computed, \em i.e., #VariableNames and \p Value corresponds to the vector storing the derivative values. #set_membrane_water_content 
             * should be called before using this method. Also, this method should only be called when atleast there is one derivative flag.
             */
            virtual void electroosmotic_drag_derivative(std::map< VariableNames, std::vector<double> >&) const;
            
            /**
             * Compute the thermo-osmotic diffusion coefficient, [\p gm/ \p (cm-s-K \p )],  inside the nafion at every quadrature 
             * point in the cell. It takes vector as an input argument and values are passed by reference. Note that the double return method is not implemented for thermo-osmotic diffusion 
             * coefficients because it is required only when \f$ T\f$ is one of the solution variables in the application. #set_temperature should be called before using this method.
             */
            virtual void thermoosmotic_coeff(std::vector<double>&) const;
            /**
             * Compute the derivatives of thermo-osmotic diffusion coefficient, at every quadrature point in the cell. The derivatives are computed based on the flags
             * set by the #set_derivative_flags method. It takes map as an input argument by reference, in  which \p Key corresponds to the variable about 
             * which derivative is being computed, \em i.e., #VariableNames and \p Value corresponds to the vector storing the derivative values. #set_temperature should 
             * be called before using this method. Also, this method should only be called when atleast there is one derivative flag.
             */
            virtual void thermoosmotic_coeff_derivative(std::map< VariableNames, std::vector<double> >&) const;
            
            /**
             * Compute the oxygen diffusivity, [\p cm^2/s], inside the nafion for a given \f$ T \f$. It can be called by either constant or variable solution setting
             * methods to avoid duplicate equations. The first parameter is the output oxygen diffusivity and the second one is temperature, \f$ T \f$ [\p K], used as input parameter.
             */
            virtual void oxygen_diffusivity(double&, const double) const;
            /**
             * Compute the oxygen diffusivity, [\p cm^2/s], inside the nafion for a constant \f$ T \f$. Before calling this method, #set_T method should be used for 
             * setting the constant temperature value. It takes a double as an input argument and value is passed by reference.
             */
            virtual void oxygen_diffusivity(double&) const;

            /**
             * Compute the hydrogen diffusivity, [\p cm^2/s], inside the nafion for a constant \f$ T \f$. Before calling this method, #set_T method should be used for
             * setting the constant temperature value. It takes a double as an input argument and value is passed by reference.
             */
            virtual void hydrogen_diffusivity(double&) const;
            /**
             * Compute the oxygen diffusivity, [\p cm^2/s], inside the nafion at every quadrature point in the cell. It takes
             * vector as an input argument and values are passed by reference. #set_temperature should be called before using this method.
             */
            virtual void oxygen_diffusivity(std::vector<double>&) const;
            /**
             * Compute the derivatives of oxygen diffusivity, at every quadrature point in the cell. The derivatives are computed based on the flags
             * set by the #set_derivative_flags method. It takes map as an input argument by reference, in  which \p Key corresponds to the variable about 
             * which derivative is being computed, \em i.e., #VariableNames and \p Value corresponds to the vector storing the derivative values. #set_temperature should 
             * be called before using this method. Also, this method should only be called when atleast there is one derivative flag.
             */
            virtual void oxygen_diffusivity_derivative(std::map< VariableNames, std::vector<double> >&) const;
            
            /**
             * Compute the proton diffusivity, [\p cm^2/s], inside the nafion for constant case. It takes a double as an 
             * input argument and value is passed by reference. It basically returns a constant value specified using the parameter file.
             */
            virtual void proton_diffusivity(double&) const;

            /**
             * Compute the enthalpy of sorption [\p J/mol] of water in Nafion, at all quadrature points in the cell. It takes
             * vector as an input argument and values are passed by reference. Note that the double return method is not implemented for thermo-osmotic diffusion 
             * coefficients because it is required only when \f$ T\f$ is one of the solution variables in the application. #set_temperature should be called before using this method.
             */
            virtual void sorption_enthalpy(std::vector<double>& ) const;
            /**
             * Compute the derivatives of enthalpy of sorption of water, at every quadrature point in the cell. The derivatives are computed based on the flags
             * set by the #set_derivative_flags method. It takes map as an input argument by reference, in  which \p Key corresponds to the variable about 
             * which derivative is being computed, \em i.e., #VariableNames and \p Value corresponds to the vector storing the derivative values. #set_temperature 
             * should be called before using this method. Also, this method should only be called when atleast there is one derivative flag.
             */
            virtual void sorption_enthalpy_derivative(std::map< VariableNames, std::vector<double> >&) const;

            //@}
            
            ///@name Accessor methods for molar enthalpy of sorbed water and its derivatives
            //@{
            
            /**
             * Compute the molar enthalpy, \f$ \bar{H}_{\lambda} ~\f$ [\p J/mol] of sorbed water in the Nafion as a function of \f$ T \f$. It takes 
             * \b Temperature, \f$ T ~\f$ [\p Kelvin] as input by reference and returns the double value corresponding to molar enthalpy.
             */
            virtual double get_Hlambda(const double&) const;
            
            /**
             * Compute \f$ \frac{\partial  \bar{H}_{\lambda}}{\partial T} \f$ of sorbed water in the Nafion as a function of \f$ T \f$. It takes 
             * \b Temperature, \f$ T ~\f$ [\p Kelvin] as input by reference and returns the double value corresponding to the derivative.
             */
            virtual double get_dHlambda_dT(const double&) const;
            
            /**
             * Compute \f$ \frac{\partial^2  \bar{H}_{\lambda}}{\partial T^2} \f$ of sorbed water in the Nafion as a function of \f$ T \f$. It takes 
             * \b Temperature, \f$ T ~\f$ [\p Kelvin] as input by reference and returns the double value corresponding to the second derivative.
             */
            virtual double get_d2Hlambda_dT2(const double&) const;
            
            //@}
            
        private:
            
            ///@name Instance Delivery
            //@{
             /**
             * This member function is used to create an object of type Nafion material.
             * 
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::Material::PolymerElectrolyteBase > create_replica ()
            {
                return boost::shared_ptr<FuelCellShop::Material::PolymerElectrolyteBase > (new FuelCellShop::Material::Nafion ());
            }
            /**
             * Create prototype for the layer
             */            
            static Nafion const* PROTOTYPE;
            //@}
            
            /**
             * Method to modify parameters, without declaring and initializing the parameter file. It takes two arguments, first one (unsigned int) corresponds 
             * to a particular type of method and second one (std::string) corresponds to various implementations for that method type. \b index corresponds to following methods:
             *  - \b 1: protonic conductivity methods
             *  - \b 2: water diffusivity methods
             *  - \b 3: electro-osmotic drag methods
             *  - \b 4: thermo-osmotic coefficient methods
             *  - \b 5: sorption isotherm methods
             *  - \b 6: sorption enthalpy methods 
             * \note This function is used by the unit testing for this class. For \b method strings, look at declare_parameters() method.
             */
            inline void modify_parameters(const unsigned int& index, const std::string& method)
            {
                Assert(index>=1 && index<=6, ExcMessage("Wrong index input in Nafion::modify_parameters."));
                
                if (index == 1)
                    method_conductivity = method;
                else if (index == 2)
                    method_diffusivity = method;
                else if (index == 3)
                    method_electroosmotic_drag = method;
                else if (index == 4)
                    method_thermoosmosis = method;
                else if (index == 5)
                    method_sorption = method;
                else if (index == 6)
                    method_enthalpy_sorption = method;
            }
            
            /**
             * Map for proton conductivity's coefficients; used with "Springer" method.
             * Conductivity is given by 
             * 
             * \f$ \quad \sigma = \left( A \lambda + B \right) exp \left[ C \left( \frac{1}{303} - \frac{1}{T} \right) \right] \quad \f$ [\p S/cm]
             * 
             */
            std::map<std::string, double> springer_coeffs;

        };
    }
}
#endif


      
      
    
