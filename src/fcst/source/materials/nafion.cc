//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: nafion.cc
//    - Description: Class to represent and return bulk properties for Nafion (Polymer Electrolyte Membrane)
//    - Developers: M. Secanell (2009) and Madhur Bhaiya (2012-13)
//    - Id: $Id: nafion.cc 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#include <materials/nafion.h>

namespace NAME = FuelCellShop::Material;

const std::string NAME::Nafion::concrete_name ("Nafion");

NAME::Nafion const* NAME::Nafion::PROTOTYPE = new NAME::Nafion(true);

//---------------------------------------------------------------------------
NAME::Nafion::Nafion(std::string name)
: FuelCellShop::Material::PolymerElectrolyteBase(name)
{
    EW = 1100.0;
    rho_M = 2.0;                // gm/cm^3
    permittivity = 20.0;
    //
    H_O2 = 3.1664e10;           // Pa-cm^3/mol
    H_H2 = 6.69e10;             // Pa-cm^3/mol
    //
    method_sorption = "Liu09";
    //
    method_conductivity = "NRE211";
    sigma_p = 0.1;              // S/cm .. called when "Constant" method is used
    springer_coeffs["a"] = 0.005139;
    springer_coeffs["b"] = -0.00326;
    springer_coeffs["c"] = 1268.0;

    //
    method_diffusivity  = "Motupally";
    diffusion_w = 2.e-3;         // cm^2/s .. called when "Constant" method is used
    //
    method_electroosmotic_drag = "Springer";
    given_n_drag = 1.0;
    //
    method_thermoosmosis = "Kim09";
    given_thermoosmotic_coeff = -1.3e-7;        // gm/(cm-s-K) .. called when "Constant" method is used
    //
    method_enthalpy_sorption = "Constant";
    given_enthalpy_sorption = 45000.0;          // J/mol .. called when "Constant" method is used
    //
    D_O2 = 9.726e-6;            // cm^2/s
    D_Protons = 9.2e-5;         // cm^2/s
}

//---------------------------------------------------------------------------
NAME::Nafion::Nafion(const bool create_replica)
:
PolymerElectrolyteBase ()
{
    if (create_replica)
        this->get_mapFactory()->insert(std::pair<std::string, PolymerElectrolyteBase* > (this->concrete_name, this) );           
}

//---------------------------------------------------------------------------
NAME::Nafion::~Nafion()
{}

//---------------------------------------------------------------------------
void
NAME::Nafion::declare_parameters (ParameterHandler &param) const
{
    param.enter_subsection("Materials");
    {
        param.enter_subsection(NAME::Nafion::concrete_name);
        {
            param.declare_entry("Method for sorption isotherm",
                                "Liu09", //1/s
                                Patterns::Selection("Hinatsu|Liu09"),
                                                    "Method to compute equilibrium lambda from sorption isotherm - Depends on Water molar fraction and Temperature");
            //-----
            param.declare_entry ("Method to compute proton conductivity",
                                 "NRE211",
                                 Patterns::Selection("Constant|Springer|NRE211|Iden11"),
                                                     "Method used to compute proton conductivity inside the membrane - Depends on lambda and T");

            param.enter_subsection("Springer coefficients");
            {
                param.declare_entry ("Proton conductivity's first coefficient",
                                "0.005139",
                                Patterns::Double(),
                                "Proton conductivity's first coefficient - Depends on T (reference value at 30C)");            
                param.declare_entry ("Proton conductivity's second coefficient",
                                "-0.00326",
                                Patterns::Double(),
                                "Proton conductivity's second coefficient - Depends on T (reference value at 30C)");            
                param.declare_entry ("Proton conductivity's third coefficient",
                                "1268.0",
                                Patterns::Double(),
                                "Proton conductivity's third coefficient");  
           }
            param.leave_subsection();            
            
                                
            param.declare_entry ("Proton conductivity [S/cm]",
                                "0.1",
                                Patterns::Double(),
                                "Proton conductivity inside the membrane [S/cm]");                      
            //-----
            param.declare_entry ("Method to compute water diffusion",
                                 "Motupally",
                                 Patterns::Selection("Constant|Springer|Motupally"),
                                                     "Method used to compute diffusion of water inside the membrane - Depends on lambda and T");
                                 param.declare_entry ("Water diffusion coefficient [cm^2/s]",
                                 "2e-3",
                                 Patterns::Double(),
                                 "Water diffusion coefficient inside the membrane [cm^2/s]");
            //-----
            param.declare_entry ("Electro-osmotic drag method",
                                 "Springer",
                                 Patterns::Selection("Constant|Springer"),
                                                     "Method to compute nd. Springer depends on water content");
                                 param.declare_entry ("Electro-osmotic drag coefficient",
                                 "1.0",
                                 Patterns::Double(),
                                 "Number of water molecules dragged by one proton");
            //-----
            param.declare_entry ("Method to compute thermo-osmotic diffusion coefficient",
                                 "Kim09",
                                 Patterns::Selection("Constant|Kim09"),
                                                     "Method to compute thermo-osmotic diffusion coefficient - Depends on T");
                                 param.declare_entry ("Thermo-osmotic diffusion coefficient [gm/(cm-s-K)]",
                                 "-1.3e-7",
                                 Patterns::Double(),
                                 "Thermo-osmotic diffusion coefficient inside the membrane [gm/(cm-s-K)].");
            //-----
            param.declare_entry ("Method to compute enthalpy of sorption of water",
                                 "Constant",
                                 Patterns::Selection("Constant"),
                                                     "Method to compute enthalpy of sorption of water - may depend on lambda and T");
                                 param.declare_entry ("Enthalpy of sorption of water [J/mol]",
                                 "45000.0",
                                 Patterns::Double(),
                                 "Heat released when one mole of water is sorbed inside the membrane.");
            //-----
            param.declare_entry ("Oxygen diffusion coefficient [cm^2/s]",
                                 "9.726e-6",//cm^2/s @ 353K  //From J. Peron et al., “Properties of Nafion NR-211 membranes for PEMFCs”
                                 Patterns::Double());
            param.declare_entry ("Proton diffusion coefficient [cm^2/s]",
                                 "9.2e-5",//Peter Dobsons thesis page 42
                                 Patterns::Double());
            param.declare_entry ("Henry's Law Constant for Oxygen [Pa cm^3/mol]",
                                 "3.1664e10",
                                 Patterns::Double());
            param.declare_entry ("Henry's Law Constant for Hydrogen [Pa cm^3/mol]",
                                 "6.69e10",
                                 Patterns::Double());
            param.declare_entry ("Hydrogen diffusion coefficient [cm^2/s]",
                                 "12.8e-6",
                                 Patterns::Double());

        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
void
NAME::Nafion::initialize (ParameterHandler &param)
{    
    param.enter_subsection("Materials");
    {
        param.enter_subsection(NAME::Nafion::concrete_name);
        {
            method_sorption = param.get("Method for sorption isotherm");
            method_conductivity = param.get("Method to compute proton conductivity");
            sigma_p = param.get_double("Proton conductivity [S/cm]");            
            method_diffusivity = param.get("Method to compute water diffusion");
            diffusion_w = param.get_double("Water diffusion coefficient [cm^2/s]");
            method_electroosmotic_drag = param.get("Electro-osmotic drag method");
            given_n_drag = param.get_double("Electro-osmotic drag coefficient");
            method_thermoosmosis = param.get("Method to compute thermo-osmotic diffusion coefficient");
            given_thermoosmotic_coeff = param.get_double("Thermo-osmotic diffusion coefficient [gm/(cm-s-K)]");
            method_enthalpy_sorption = param.get("Method to compute enthalpy of sorption of water");
            given_enthalpy_sorption = param.get_double("Enthalpy of sorption of water [J/mol]");
            D_Protons =  param.get_double("Proton diffusion coefficient [cm^2/s]");
            D_O2 = param.get_double("Oxygen diffusion coefficient [cm^2/s]");
            D_H2 = param.get_double("Hydrogen diffusion coefficient [cm^2/s]");
            H_O2 = param.get_double("Henry's Law Constant for Oxygen [Pa cm^3/mol]");
            H_H2 = param.get_double("Henry's Law Constant for Hydrogen [Pa cm^3/mol]");
            param.enter_subsection("Springer coefficients");
            {
                springer_coeffs["a"] = param.get_double("Proton conductivity's first coefficient");
                springer_coeffs["b"] = param.get_double("Proton conductivity's second coefficient");
                springer_coeffs["c"] = param.get_double("Proton conductivity's third coefficient");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    
    param.leave_subsection();
}

//---------------------------------------------------------------------------
void
NAME::Nafion::sorption_isotherm(std::vector<double>& Leq_vec) const
{
    Assert( p_total!=0., ExcMessage("Total pressure value has not been set inside the Nafion material class before Nafion::sorption_isotherm.") );
    Assert( xwater_var.is_initialized() && T_var.is_initialized(), ExcMessage("Water molar fraction and/or temperature values have not been set before Nafion::sorption_isotherm.") );
    Assert( T_var[0]!=0., ExcMessage("Temperature values have not been set in Nafion material class using set_T method before Nafion::sorption_isotherm.") );
    
    Leq_vec.resize(xwater_var.size());
    
    if (method_sorption == "Hinatsu")
    {
        const double x_tr = 0.975;
        const double sharp = 100.0;
        for (unsigned int i=0; i<Leq_vec.size(); ++i)
        {
            double aw = (p_total*xwater_var[i])/(this->water_mat->get_water_vapor_saturation_pressure(T_var[i]));
            Leq_vec[i] = ( (0.3 + 10.8*aw - 16.0*pow(aw,2.0) + 14.1*pow(aw, 3.0))*0.5*(1.0+tanh(sharp*(x_tr - aw))) ) + ( 9.2*0.5*(1.0+tanh(sharp*(aw - x_tr))) );
        }
    }
    
    else if (method_sorption == "Liu09")
    {
        const double x_tr = 0.975;
        const double sharp = 100.0;
        
        for (unsigned int i=0; i<Leq_vec.size(); ++i)
        {
            double aw = (p_total*xwater_var[i])/(this->water_mat->get_water_vapor_saturation_pressure(T_var[i]));            
            Leq_vec[i] = ( (1.0 + aw*aw*0.2352*(T_var[i]-303.15)/30.0)*(14.22*aw*aw*aw - 18.92*aw*aw + 13.41*aw)*0.5*(1.0 + tanh(sharp*(x_tr-aw))) )
            + 
            ( (8.71 + 0.0682864*(T_var[i]-303.15)) * 0.5 * (1.0 + tanh(sharp*(aw-x_tr))) );
        }
    }
    
    else
        AssertThrow (false, ExcNotImplemented ());
    
    /* In reality (for Hinatsu)
     i f *(aw <= 1.0)
    lambda = (0.3 + 10.8*aw - 16*pow(aw,2.0) + 14.1*pow(aw, 3.0));
    else
        lambda = 9.2 + 0.01*aw;
    */
}

//---------------------------------------------------------------------------
void
NAME::Nafion::sorption_isotherm_derivative(std::map < VariableNames, std::vector<double> >& dLeq) const
{
    Assert( p_total!=0., ExcMessage("Total pressure value has not been set inside the Nafion material class before Nafion::sorption_isotherm_derivative.") );
    Assert( T_var.is_initialized() && xwater_var.is_initialized(), ExcMessage("Water molar fraction and/or temperature values have not been set before Nafion::sorption_isotherm_derivative.") );
    Assert( derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called before Nafion::sorption_isotherm_derivative.") );
    Assert( T_var[0]!=0., ExcMessage("Temperature values have not been set in Nafion material class using set_T method before Nafion::sorption_isotherm_derivative.") );
    
    // Looping over derivative flags
    for (unsigned int i=0; i<derivative_flags.size(); ++i)
    {
        std::vector<double> deriv_Leq(xwater_var.size(), 0.);      // Sets all derivatives by default to zero.
        
        if (derivative_flags[i] == water_molar_fraction)
        {
            const double x_tr = 0.975;
            const double sharp = 100.0;
        
            if (method_sorption == "Hinatsu")
                for (unsigned int j = 0; j<xwater_var.size(); ++j)
                {
                    double aw = (p_total*xwater_var[j])/(this->water_mat->get_water_vapor_saturation_pressure(T_var[j]));
                    double daw_dxwater = p_total/(this->water_mat->get_water_vapor_saturation_pressure(T_var[j]));
                    
                    deriv_Leq[j] = 0.5* ( (10.8 - 16.*2.*aw + 14.1*3.*aw*aw)*(1.+tanh(sharp*(x_tr - aw))) + (0.3 + 10.8*aw - 16.*aw*aw + 14.1*aw*aw*aw)*pow(1./cosh(sharp*(x_tr - aw)),2.)*(-sharp)
                                            +
                                        9.2*pow(1./cosh(sharp*(aw - x_tr)),2.)*sharp ) * daw_dxwater;
                }
                    
            else if (method_sorption == "Liu09")
            {
                for (unsigned int j = 0; j<xwater_var.size(); ++j)
                {
                    double aw = (p_total*xwater_var[j])/(this->water_mat->get_water_vapor_saturation_pressure(T_var[j]));
                    double daw_dxwater = p_total/(this->water_mat->get_water_vapor_saturation_pressure(T_var[j]));
                    
                    deriv_Leq[j] = 0.5* ( (2.*aw*0.2352*(T_var[j]-303.15)/30.)*(14.22*aw*aw*aw - 18.92*aw*aw + 13.41*aw)*(1. + tanh(sharp*(x_tr-aw)))    + 
                                        (1. + aw*aw*0.2352*(T_var[j]-303.15)/30.)*(14.22*3.*aw*aw - 18.92*2.*aw + 13.41)*(1. + tanh(sharp*(x_tr-aw)))  + 
                                        (1. + aw*aw*0.2352*(T_var[j]-303.15)/30.)*(14.22*aw*aw*aw - 18.92*aw*aw + 13.41*aw)*pow(1./cosh(sharp*(x_tr - aw)),2.)*(-sharp)
                                             +
                                         (8.71 + 0.0682864*(T_var[j]-303.15))*pow(1./cosh(sharp*(aw - x_tr)),2.)*sharp ) * daw_dxwater;
                }

            }
            
            else
                AssertThrow (false, ExcNotImplemented ());
            
            dLeq[ derivative_flags[i] ] = deriv_Leq;
        }
        
        
        else if (derivative_flags[i] == temperature_of_REV)
        {
            const double x_tr = 0.975;
            const double sharp = 100.0;
            
            if (method_sorption == "Hinatsu")
                for (unsigned int j = 0; j<T_var.size(); ++j)
                {
                    double aw = (p_total*xwater_var[j])/(this->water_mat->get_water_vapor_saturation_pressure(T_var[j]));
                    double daw_dT = (-1.*p_total*xwater_var[j]*this->water_mat->get_Dwater_vapor_saturation_pressure_Dtemperature(T_var[j]))/(pow(this->water_mat->get_water_vapor_saturation_pressure(T_var[j]),2.));
                    
                    deriv_Leq[j] = 0.5* ( (10.8 - 16.*2.*aw + 14.1*3.*aw*aw)*(1.+tanh(sharp*(x_tr - aw))) + (0.3 + 10.8*aw - 16.*aw*aw + 14.1*aw*aw*aw)*pow(1./cosh(sharp*(x_tr - aw)),2.)*(-sharp)
                                            +
                                        9.2*pow(1./cosh(sharp*(aw - x_tr)),2.)*sharp ) * daw_dT;
                }
            
            else if (method_sorption == "Liu09")
                for (unsigned int j = 0; j<T_var.size(); ++j)
                {
                    double aw = (p_total*xwater_var[j])/(this->water_mat->get_water_vapor_saturation_pressure(T_var[j]));
                    double daw_dT = (-1.*p_total*xwater_var[j]*this->water_mat->get_Dwater_vapor_saturation_pressure_Dtemperature(T_var[j]))/(pow(this->water_mat->get_water_vapor_saturation_pressure(T_var[j]),2.));
                    
                    deriv_Leq[j] = 0.5* ( ((0.2352/30.)*(2.*aw*daw_dT*(T_var[j]-303.15) + aw*aw))*(14.22*aw*aw*aw - 18.92*aw*aw + 13.41*aw)*(1. + tanh(sharp*(x_tr-aw)))    + 
                                        (1. + aw*aw*0.2352*(T_var[j]-303.15)/30.)*(14.22*3.*aw*aw - 18.92*2.*aw + 13.41)*daw_dT*(1. + tanh(sharp*(x_tr-aw)))  + 
                                        (1. + aw*aw*0.2352*(T_var[j]-303.15)/30.)*(14.22*aw*aw*aw - 18.92*aw*aw + 13.41*aw)*pow(1./cosh(sharp*(x_tr - aw)),2.)*(-sharp)*daw_dT
                                             +
                                         (8.71 + 0.0682864*(T_var[j]-303.15))*pow(1./cosh(sharp*(aw - x_tr)),2.)*sharp*daw_dT  +
                                          0.0682864*(1.0 + tanh(sharp*(aw-x_tr))) );                
                }
            
            else
                AssertThrow (false, ExcNotImplemented ());
            
            dLeq[ derivative_flags[i] ] = deriv_Leq;
        }
        
        else
            dLeq[ derivative_flags[i] ] = deriv_Leq;
        
    }   // End loop over derivative flags
}

//---------------------------------------------------------------------------
void
NAME::Nafion::proton_conductivity(double& S, const double T, const double lambda) const
{
    if (method_conductivity == "Constant")
        S = sigma_p;
    
    else if (method_conductivity == "Springer")
        S = (springer_coeffs.at("a")*lambda + springer_coeffs.at("b"))* exp(springer_coeffs.at("c")*(1.0/303.0 - 1.0/T) ); //T in kelvin

    else if (method_conductivity == "NRE211")
        S = ((-0.00010125)*pow(lambda,2.0) + 0.01052*lambda - 0.020634)* exp(751.5412*(1.0/303.0 - 1.0/T) ); //T in kelvin

    else if (method_conductivity == "Iden11")
    {
        // If using Idens method, the water content is passed to the effective proton conductivity function via the bulk conductivity
        //Compute water activity based on lambda
        
        double RH(0.0);         // activity of water equivalent to Relative Humidity 
        if (lambda < 13.0)
            RH = ( 0.000093561521965*pow(lambda,3.0) - 0.008649693377387*pow(lambda,2.0) + 0.183227305676982*lambda - 0.125414456218017 ) * 100.0;
        else
            RH = 100.0;
        
        S = ( 1.93134146e-7*pow(RH,3.0) - 6.73473135e-6*pow(RH,2.0) + 0.000745491046859*RH - 0.007977537491969 ) * exp (751.5412* (1.0/353.0 - 1.0/T) );
    }
    
    else
        AssertThrow (false, ExcNotImplemented ());
}

//---------------------------------------------------------------------------
void
NAME::Nafion::proton_conductivity(double& S) const
{
    Assert ((T != 0), ExcMessage("Temperature values have not been set in Nafion material class using set_T method before Nafion::proton_conductivity."));
    this->proton_conductivity(S, T, lambda);
}

//---------------------------------------------------------------------------
void
NAME::Nafion::proton_conductivity(std::vector<double>& S_vec) const
{
    if ( !T_var.is_initialized() && !lambda_var.is_initialized() )
    {
        Assert(S_vec.size() != 0, ExcMessage("Vector needs to be initialized before passing to this function, when we have constant temperature and lambda inside the membrane."));
        double sigma_m;
        proton_conductivity(sigma_m);
        for (unsigned int i=0; i<S_vec.size(); ++i)
            S_vec[i] = sigma_m;
        
        return;
    }
    
    Assert( T_var.is_initialized() && lambda_var.is_initialized(), ExcInternalError() );
    Assert( T_var[0] != 0., ExcMessage("Temperature values have not been set in Nafion material class using set_T method before Nafion::proton_conductivity.") );

    S_vec.resize(T_var.size());
        
    for (unsigned int i=0; i<S_vec.size(); ++i)
        this->proton_conductivity(S_vec[i], T_var[i], lambda_var[i]);
}

//---------------------------------------------------------------------------
void
NAME::Nafion::proton_conductivity_derivative(std::map< VariableNames, std::vector<double> >& dS) const
{
    Assert( T_var.is_initialized() && lambda_var.is_initialized(), ExcMessage("Nafion::proton_conductivity_derivative should only be used when atleast one of the lambda/temperature is variable.") );
    Assert( derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called before Nafion::proton_conductivity_derivative.") );
    Assert( T_var[0]!=0., ExcMessage("Temperature values have not been set in Nafion material class using set_T method before Nafion::proton_conductivity_derivative.") );
    
    // Looping over derivative flags
    for (unsigned int i=0; i<derivative_flags.size(); ++i)
    {
        std::vector<double> derivS_vec(lambda_var.size(), 0.);      // Sets all derivatives by default to zero.
        
        if (derivative_flags[i] == temperature_of_REV)
        {
            if (method_conductivity == "Constant")
                for (unsigned int j=0; j<T_var.size(); ++j)
                    derivS_vec[j] = 0.0;

            else if (method_conductivity == "Springer")
                for (unsigned int j=0; j<T_var.size(); ++j)
                    derivS_vec[j] = (springer_coeffs.at("a")*lambda_var[j] + springer_coeffs.at("b"))* exp(springer_coeffs.at("c")*(1.0/303.0 - 1.0/(T_var[j])) ) *springer_coeffs.at("c") * (1.0/pow((T_var[j]),2.0));
                
            else if (method_conductivity == "NRE211")
                for (unsigned int j=0; j<T_var.size(); ++j)
                    derivS_vec[j] = ((-0.00010125)*pow(lambda_var[j],2.0) + 0.01052*lambda_var[j] - 0.020634) * exp(751.5412*(1.0/303.0 - 1.0/(T_var[j])) ) *751.5412*(1.0/pow((T_var[j]),2.0));
                    
            else if (method_conductivity == "Iden11")
                for (unsigned int j=0; j<T_var.size(); ++j)
                {
                    double RH;
                    if (lambda_var[j] < 13.0)
                        RH = ( 0.000093561521965*pow(lambda_var[j],3.0) - 0.008649693377387*pow(lambda_var[j],2.0) + 0.183227305676982*lambda_var[j] - 0.125414456218017 ) * 100.0;
                    else
                        RH = 100.0;
                    
                    derivS_vec[j] = ( 1.93134146e-7*pow(RH,3.0) - 6.73473135e-6*pow(RH,2.0) + 0.000745491046859*RH - 0.007977537491969 ) *exp (751.5412* (1.0/353.0 - 1.0/T_var[j]) ) * ( 751.5412/(T_var[j]*T_var[j]) );
                }
            else
                AssertThrow (false, ExcNotImplemented ());
            
            dS[ derivative_flags[i] ] = derivS_vec;
                    
        }
        
        else if (derivative_flags[i] == membrane_water_content)
        {
            if (method_conductivity == "Constant")
                for (unsigned int j=0; j<lambda_var.size(); ++j)
                    derivS_vec[j] = 0.0;

            else if (method_conductivity == "Springer")
                for (unsigned int j=0; j<lambda_var.size(); ++j)
                    derivS_vec[j] = springer_coeffs.at("a")* exp(springer_coeffs.at("c")*(1.0/303.0 - 1.0/(T_var[j])) );
                
            else if (method_conductivity == "NRE211")
                for (unsigned int j=0; j<lambda_var.size(); ++j)
                    derivS_vec[j] = ((2.0*(-0.00010125)*lambda_var[j]) + 0.01052)* exp(751.5412*(1.0/303.0 - 1.0/(T_var[j])));
                    
            else if (method_conductivity == "Iden11")
                for (unsigned int j=0; j<lambda_var.size(); ++j)
                {
                    double RH, dRH_dlambda;
                    if (lambda_var[j] < 13.0)
                    {
                        RH = ( 0.000093561521965*pow(lambda_var[j],3.0) - 0.008649693377387*pow(lambda_var[j],2.0) + 0.183227305676982*lambda_var[j] - 0.125414456218017 ) * 100.0;
                        dRH_dlambda = ( 0.000093561521965*3.0*pow(lambda_var[j],2.0) - 0.008649693377387*2.0*lambda_var[j] + 0.183227305676982 ) * 100.0;
                    }
                    else
                    {
                        RH = 100.0;
                        dRH_dlambda = 0.0;
                    }
                    
                    derivS_vec[j] = ( 1.93134146e-7*3.0*pow(RH,2.0)*dRH_dlambda - 6.73473135e-6*2.0*RH*dRH_dlambda + 0.000745491046859*dRH_dlambda ) * exp (751.5412* (1.0/353.0 - 1.0/T_var[j]) );
                }
            else
                AssertThrow (false, ExcNotImplemented ());
            
            dS[ derivative_flags[i] ] = derivS_vec;       
        }
        
        else
            dS[ derivative_flags[i] ] = derivS_vec;
        
    }   // End loop over derivative flags
}

//---------------------------------------------------------------------------
void
NAME::Nafion::water_diffusivity(double& D_w, const double Temp, const double lamb_var) const
{
    
    if (method_diffusivity == "Constant")
        D_w = diffusion_w;

    else if (method_diffusivity == "Springer")
        D_w = 1.e-6*(2.563 - 0.33*lamb_var + 0.0264*pow(lamb_var,2.0) - 0.000671*pow(lamb_var,3.0))* exp(2416.0*(1.0/303.0 - 1.0/Temp));

    else if (method_diffusivity == "Motupally")
        if (lamb_var <= 3. && lamb_var > 0.)
            D_w = 3.10e-3*lamb_var*(-1.0 + exp(0.28*lamb_var))*exp(-2436.0/Temp);
        else
            D_w = 4.17e-4*lamb_var*(1.0 + 161.0*exp(-lamb_var))*exp(-2436.0/Temp);

    else
        AssertThrow (false, ExcNotImplemented ());

}

//---------------------------------------------------------------------------
void
NAME::Nafion::water_diffusivity(double& D_w) const
{
    Assert ((T != 0), ExcMessage("Temperature values have not been set in Nafion material class using set_T method before Nafion::water_diffusivity."));
    this->water_diffusivity( D_w, T, lambda);
}

//---------------------------------------------------------------------------
void
NAME::Nafion::water_diffusivity(std::vector<double>& Dw_vec) const
{
    Assert( T_var.is_initialized() && lambda_var.is_initialized(), ExcInternalError() );
    Assert( T_var[0] != 0., ExcMessage("Temperature values have not been set in Nafion material class using set_T method before Nafion::water_diffusivity.") );

    Dw_vec.resize(T_var.size());
    
    for (unsigned int i=0; i<Dw_vec.size(); ++i)
        this->water_diffusivity(Dw_vec[i],T_var[i],lambda_var[i]);
}

//---------------------------------------------------------------------------
void
NAME::Nafion::water_diffusivity_derivative(std::map< VariableNames, std::vector<double> >& dD_w) const
{
    Assert( T_var.is_initialized() && lambda_var.is_initialized(), ExcMessage("Nafion::water_diffusivity_derivative should only be used when atleast one of the lambda/temperature is variable.") );
    Assert( derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called before Nafion::water_diffusivity_derivative.") );
    Assert( T_var[0]!=0., ExcMessage("Temperature values have not been set in Nafion material class using set_T method before Nafion::water_diffusivity_derivative.") );
    
    // Looping over derivative flags
    for (unsigned int i=0; i<derivative_flags.size(); ++i)
    {
        std::vector<double> deriv_Dvec(lambda_var.size(), 0.);          // Sets all derivatives by default to zero.
        
        if (derivative_flags[i] == temperature_of_REV)
        {
            if (method_diffusivity == "Constant")
                for (unsigned int j=0; j<T_var.size(); ++j)
                    deriv_Dvec[j] = 0.0;
                
            else if (method_diffusivity == "Springer")
                for (unsigned int j=0; j<T_var.size(); ++j)
                    deriv_Dvec[j] = 1.e-6*(2.563 - 0.33*lambda_var[j] + 0.0264*pow(lambda_var[j],2.0) - 0.000671*pow(lambda_var[j],3.0)) * (exp(2416.0*(1.0/303.0 - 1.0/T_var[j]))*2416.0*(1.0/pow(T_var[j],2.0)));
                    
            else if (method_diffusivity == "Motupally")
                for (unsigned int j=0; j<T_var.size(); ++j)
                {
                    if (lambda_var[j] <= 3. && lambda_var[j] > 0.)
                        deriv_Dvec[j] = 3.10e-3*lambda_var[j]*(-1.0 + exp(0.28*lambda_var[j]))*exp(-2436.0/T_var[j])*(2436.0/pow(T_var[j],2.0));
                    else
                        deriv_Dvec[j] = 4.17e-4*lambda_var[j]*(1.0 + 161.0*exp(-lambda_var[j]))*exp(-2436.0/T_var[j])*(2436.0/pow(T_var[j],2.0));
                }

            else
                AssertThrow (false, ExcNotImplemented ());
            
            dD_w[ derivative_flags[i] ] = deriv_Dvec;          
        }
        
        else if (derivative_flags[i] == membrane_water_content)
        {
            if (method_diffusivity == "Constant")
                for (unsigned int j=0; j<T_var.size(); ++j)
                    deriv_Dvec[j] = 0.0;
                
            else if (method_diffusivity == "Springer")
                for (unsigned int j=0; j<lambda_var.size(); ++j)
                    deriv_Dvec[j] = 1.e-6*((-0.33) + 2*0.0264*lambda_var[j] - 3*0.000671*pow(lambda_var[j],2.0))* exp(2416.0*(1.0/303.0 - 1.0/T_var[j]));
                    
            else if (method_diffusivity == "Motupally")
                for (unsigned int j=0; j<lambda_var.size(); ++j)
                {
                    if (lambda_var[j] <= 3. && lambda_var[j] > 0.)
                        deriv_Dvec[j] = (3.10e-3*(-1.0 + exp(0.28*lambda_var[j]))*exp(-2436.0/T_var[j])) + (3.10e-3*lambda_var[j]*(0.28*exp(0.28*lambda_var[j]))*exp(-2436.0/T_var[j]));
                    else
                        deriv_Dvec[j] = (4.17e-4*(1.0 + 161.0*exp(-lambda_var[j]))*exp(-2436.0/T_var[j])) + (4.17e-4*lambda_var[j]*(-161.0*exp(-lambda_var[j]))*exp(-2436.0/T_var[j]));
                }

            else
                AssertThrow (false, ExcNotImplemented ());
                    
            dD_w[ derivative_flags[i] ] = deriv_Dvec;
        }
        
        else
            dD_w[ derivative_flags[i] ] = deriv_Dvec;
        
    }   // End loop over derivative flags
}

//---------------------------------------------------------------------------
void
NAME::Nafion::electroosmotic_drag(std::vector<double>& n_d) const
{
    Assert( lambda_var.is_initialized(), ExcMessage("Lambda values have not been set in Nafion material class using set_membrane_water_content method before Nafion::electroosmotic_drag.") );
    
    n_d.resize(lambda_var.size());
    
    if (method_electroosmotic_drag == "Constant")
        for (unsigned int i=0; i<n_d.size(); ++i)
            n_d[i] = given_n_drag;

    else if (method_electroosmotic_drag == "Springer")
        for (unsigned int i=0; i<n_d.size(); ++i)
            n_d[i] = (2.5*lambda_var[i])/22.0;

    else
        AssertThrow (false, ExcNotImplemented ());
}

//---------------------------------------------------------------------------
void
NAME::Nafion::electroosmotic_drag_derivative(std::map< VariableNames, std::vector<double> >& Dn) const
{
    Assert( derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called before Nafion::electroosmotic_drag_derivative."));
    Assert( lambda_var.is_initialized(), ExcMessage("Lambda values have not been set in Nafion material class using set_membrane_water_content method before Nafion::electroosmotic_drag_derivative.") );
    
    // Looping over derivative flags
    for (unsigned int i=0; i<derivative_flags.size(); ++i)
    {
        std::vector<double> deriv_nd(lambda_var.size(), 0.);            // Sets all derivatives by default to zero.
        
        if (derivative_flags[i] == membrane_water_content)
        {
            if (method_electroosmotic_drag == "Constant")
                for (unsigned int j=0; j<lambda_var.size(); ++j)
                    deriv_nd[j] = 0.0;

            else if (method_electroosmotic_drag == "Springer")
                for (unsigned int j=0; j<lambda_var.size(); ++j)
                    deriv_nd[j] = 2.5/22.0;

            else
                AssertThrow (false, ExcNotImplemented ());
            
            Dn[ derivative_flags[i] ] = deriv_nd;
        }
        
        else
            Dn[ derivative_flags[i] ] = deriv_nd;
        
    }   // End loop over derivative flags
}

//---------------------------------------------------------------------------
void
NAME::Nafion::thermoosmotic_coeff(std::vector<double>& D_T) const
{
    Assert( T_var.is_initialized(), ExcMessage("Temperature values have not been set in Nafion material class using set_temperature method before Nafion::thermoosmotic_coeff.") );
    
    D_T.resize(T_var.size());
    
    if (method_thermoosmosis == "Constant")
        for (unsigned int i=0; i<D_T.size(); ++i)
            D_T[i] = given_thermoosmotic_coeff;
        
    else if (method_thermoosmosis == "Kim09")
        for (unsigned int i=0; i<D_T.size(); ++i)
            D_T[i] = (-1.04e-5) * exp(-2362.0/T_var[i]) * 10.0; // Conversion factor of 10 for converting from kg/(m-s-K) to gm/(cm-s-K)
    
    else
        AssertThrow(false, ExcNotImplemented());
}

//---------------------------------------------------------------------------
void
NAME::Nafion::thermoosmotic_coeff_derivative(std::map< VariableNames, std::vector<double> >& dD_T) const
{
    Assert( derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called before Nafion::thermoosmotic_drag_derivative."));
    Assert( T_var.is_initialized(), ExcMessage("Temperature values have not been set in Nafion material class using set_temperature method before Nafion::thermoosmotic_coeff_derivative.") );
    
    // Looping over derivative flags
    for (unsigned int i=0; i<derivative_flags.size(); ++i)
    {
        std::vector<double> deriv_DT(T_var.size(), 0.);            // Sets all derivatives by default to zero.
        
        if (derivative_flags[i] == temperature_of_REV)
        {
            if (method_thermoosmosis == "Constant")
                for (unsigned int j=0; j<T_var.size(); ++j)
                    deriv_DT[j] = 0.0;
                
            else if (method_thermoosmosis == "Kim09")
                for (unsigned int j=0; j<T_var.size(); ++j)
                    deriv_DT[j] = (-1.04e-5) * exp(-2362.0/T_var[j]) * 10.0 * (-2362.0) * (-1.0/(T_var[j]*T_var[j]));
                
            else
                AssertThrow(false, ExcNotImplemented());
            
            dD_T[ derivative_flags[i] ] = deriv_DT;
        }
        
        else
            dD_T[ derivative_flags[i] ] = deriv_DT;
        
    } // END loop over derivative flags
}

//---------------------------------------------------------------------------
void
NAME::Nafion::sorption_enthalpy(std::vector<double>& h_sorp) const
{
    Assert( T_var.is_initialized(), ExcMessage("Temperature values have not been set in Nafion material class using set_temperature method before Nafion::sorption_enthalpy.") );
    
    h_sorp.resize(T_var.size());
    
    if (method_enthalpy_sorption == "Constant")
        for (unsigned int i=0; i<h_sorp.size(); ++i)
            h_sorp[i] = given_enthalpy_sorption;
    
    else
        AssertThrow(false, ExcNotImplemented());
}

//---------------------------------------------------------------------------
void
NAME::Nafion::sorption_enthalpy_derivative(std::map< VariableNames, std::vector<double> >& dh_sorp) const
{
    Assert( derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called before Nafion::sorption_enthalpy_derivative."));
    Assert( T_var.is_initialized(), ExcMessage("Temperature values have not been set in Nafion material class using set_temperature method before Nafion::sorption_enthalpy_derivative.") );
    
    // Looping over derivative flags
    for (unsigned int i=0; i<derivative_flags.size(); ++i)
    {
        std::vector<double> deriv_hsorp(T_var.size(), 0.);            // Sets all derivatives by default to zero.
        
        if (derivative_flags[i] == temperature_of_REV)
        {
            if (method_enthalpy_sorption == "Constant")
                for (unsigned int j=0; j<T_var.size(); ++j)
                    deriv_hsorp[j] = 0.0;
                
            else
                AssertThrow(false, ExcNotImplemented());
            
            dh_sorp[ derivative_flags[i] ] = deriv_hsorp;
        }
        
        else
            dh_sorp[ derivative_flags[i] ] = deriv_hsorp;
        
    } // END loop over derivative flags
}

//---------------------------------------------------------------------------
double
NAME::Nafion::get_Hlambda(const double& T_in) const
{
    if (method_enthalpy_sorption == "Constant")
        return ( water_mat->get_molar_enthalpy(T_in) - given_enthalpy_sorption );
    
    else
        AssertThrow(false, ExcNotImplemented());
}

//---------------------------------------------------------------------------
double
NAME::Nafion::get_dHlambda_dT(const double& T_in) const
{
    if (method_enthalpy_sorption == "Constant")
        return water_mat->get_Dmolar_enthalpy_Dtemperature(T_in);
    
    else
        AssertThrow(false, ExcNotImplemented());
}

//---------------------------------------------------------------------------
double
NAME::Nafion::get_d2Hlambda_dT2(const double& T_in) const
{
    if (method_enthalpy_sorption == "Constant")
        return water_mat->get_D2molar_enthalpy_Dtemperature2(T_in);
    
    else
        AssertThrow(false, ExcNotImplemented());
}

//---------------------------------------------------------------------------
void
NAME::Nafion::oxygen_diffusivity(double& D_oxy, const double Temp) const
{
    D_oxy = D_O2*exp(-(39760.0/Constants::R())*(1.0/Temp - 1.0/353.0));        // D_02 = 9.726e-6 [cm^2/s] (Default value)
}

//---------------------------------------------------------------------------
void
NAME::Nafion::oxygen_diffusivity(double& D_oxy) const
{
    Assert( (T != 0), ExcMessage("Temperature values have not been set in Nafion material class using set_T method before Nafion::oxygen_diffusivity."));
    this->oxygen_diffusivity(D_oxy,T);
}

void
NAME::Nafion::hydrogen_diffusivity(double& D_h) const
{
    Assert( (T != 0), ExcMessage("Temperature values have not been set in Nafion material class using set_T method before Nafion::hydrogen_diffusivity."));
    D_h = D_H2*exp(-(39760.0/Constants::R())*(1.0/T - 1.0/353.0));        // D_02 = 9.726e-6 [cm^2/s] (Default value)
}

//---------------------------------------------------------------------------
void
NAME::Nafion::oxygen_diffusivity(std::vector<double>& Doxy_vec) const
{
    Assert( T_var.is_initialized(), ExcMessage("Temperature values have not been set in Nafion material class using set_temperature method before Nafion::oxygen_diffusivity.") );
    
    Doxy_vec.resize(T_var.size());
    
    // Looping over quadrature points in a cell
    for (unsigned int i=0; i<Doxy_vec.size(); ++i)
        this->oxygen_diffusivity(Doxy_vec[i], T_var[i]);         // D_02 = 9.726e-6 [cm^2/s] (Default value)
}

//---------------------------------------------------------------------------
void
NAME::Nafion::oxygen_diffusivity_derivative(std::map< VariableNames, std::vector<double> >& dDoxy) const
{
    Assert( derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called before Nafion::oxygen_diffusivity_derivative."));
    Assert( T_var.is_initialized(), ExcMessage("Temperature values have not been set in Nafion material class using set_temperature method before Nafion::oxygen_diffusivity_derivative.") );

    // Looping over derivative flags
    for (unsigned int i=0; i<derivative_flags.size(); ++i)
    {
        std::vector<double> deriv_Doxy(T_var.size(), 0.);            // Sets all derivatives by default to zero.
        
        if (derivative_flags[i] == temperature_of_REV)
        {
            for (unsigned int j=0; j<T_var.size(); ++j)
                deriv_Doxy[j] = D_O2*exp(-(39760.0/Constants::R())*(1.0/T_var[j] - 1.0/353.0))*(39760.0/Constants::R())*(1.0/pow(T_var[j],2.)); // D_02 = 9.726e-6 [cm^2/s] (Default value)
                
            dDoxy[ derivative_flags[i] ] = deriv_Doxy;
        }
        
        else
            dDoxy[ derivative_flags[i] ] = deriv_Doxy;
        
    }   // End loop over derivative flags
}

//---------------------------------------------------------------------------
void
NAME::Nafion::proton_diffusivity(double& D_H) const
{
    D_H = D_Protons;
}
