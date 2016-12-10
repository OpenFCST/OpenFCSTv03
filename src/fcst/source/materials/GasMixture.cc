// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: GasMixture.cc
// - Description: This class describes properties of gas mixtures
// - Developers: Valentin N. Zingan, University of Alberta
//
// ----------------------------------------------------------------------------

#include <materials/GasMixture.h>

namespace NAME = FuelCellShop::Material;

       //////////////////////////////////////////////////
       //////////////////////////////////////////////////
       // CONSTRUCTORS, DESTRUCTOR, AND INITIALIZATION //
       //////////////////////////////////////////////////
       //////////////////////////////////////////////////

// ---             ---
// --- Constructor ---
// ---             ---
NAME::GasMixture::GasMixture()
:
NAME::BaseMaterial()
{
       pressIsoBaric = false;
       tempIsoTherm  = false;
}

NAME::GasMixture::GasMixture(const std::string& name)
:
NAME::BaseMaterial(name)
{
       pressIsoBaric = false;
       tempIsoTherm  = false;
}

// ---            ---
// --- Destructor ---
// ---            ---

NAME::GasMixture::~GasMixture()
{ }

// ---                    ---
// --- declare_parameters ---
// ---                    ---

void
NAME::GasMixture::declare_parameters(ParameterHandler& param) const
{
  param.enter_subsection("Fuel cell data");
  {
    param.enter_subsection("Materials");
    {
      param.enter_subsection(this->name);
      {
          param.declare_entry("Partial viscosity mode",
                              "Wilke",
                              Patterns::Selection("Wilke | OmegaKG | Dynamic"),
                              " ");
          
          param.declare_entry("Isothermal fluid flow",
                              "true",
                              Patterns::Bool(),
                              " ");

          param.declare_entry("Isobaric fluid flow",
                              "true",
                              Patterns::Bool(),
                              " ");
      }
      param.leave_subsection();
    }
    param.leave_subsection();
  }
  param.leave_subsection();
}

// ---            ---
// --- initialize ---
// ---            ---

void
NAME::GasMixture::initialize(ParameterHandler& param)
{
  param.enter_subsection("Fuel cell data");
  {
    param.enter_subsection("Materials");
    {
      param.enter_subsection(this->name);
      {
          mixture_viscosity_mode = param.get("Partial viscosity mode");
          tempIsoTherm           = param.get_bool("Isothermal fluid flow");
          pressIsoBaric          = param.get_bool("Isobaric fluid flow");
      }
      param.leave_subsection();
    }
    param.leave_subsection();
  }
  param.leave_subsection();
}

       ///////////////////////
       ///////////////////////
       // SERVICE FUNCTIONS //
       ///////////////////////
       ///////////////////////

  /////////////////////////////////////////////////////////////////////////////
  // Chapman Enskog isobaric diffusion coefficient - Binary gas mixture only //
  /////////////////////////////////////////////////////////////////////////////

// ---                                                  ---
// --- get_ChapmanEnskog_isobaric_diffusion_coefficient ---
// ---                                                  ---

const double
NAME::GasMixture::get_ChapmanEnskog_isobaric_diffusion_coefficient() const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( gases.size() == 2, ExcMessage("The mixture is NOT binary. Use another function.") );
    AssertThrow( tempIsoTherm != false, ExcMessage("The mixture is NOT isothermal. Use another function.") );

    const double M1 = gases[0]->get_molar_mass() * 1.0e3;
    const double M2 = gases[1]->get_molar_mass() * 1.0e3;

    const double sigma12 = 0.5*( gases[0]->get_collision_diameter() + gases[1]->get_collision_diameter() );

    const double T = temperature;

    const double omega12 = get_binary_collision_integral();

    return 1.8829e-2 * std::sqrt( ( 1.0/M1 + 1.0/M2 )*T*T*T ) / ( sigma12*sigma12*omega12 );
}

// ---                                                  ---
// --- get_ChapmanEnskog_isobaric_diffusion_coefficient ---
// ---                                                  ---

void
NAME::GasMixture::get_ChapmanEnskog_isobaric_diffusion_coefficient(std::vector<double>& diffusion_coefficient) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( gases.size() == 2, ExcMessage("The mixture is NOT binary. Use another function.") );
    AssertThrow( tempIsoTherm != false, ExcMessage("The mixture is NOT isothermal. Use another function.") );

    const double result = get_ChapmanEnskog_isobaric_diffusion_coefficient();

    for(unsigned int q = 0; q < diffusion_coefficient.size(); ++q)
        diffusion_coefficient[q] = result;
}

// ---                                                  ---
// --- get_ChapmanEnskog_isobaric_diffusion_coefficient ---
// ---                                                  ---

const double
NAME::GasMixture::get_ChapmanEnskog_isobaric_diffusion_coefficient(const double& temp) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( gases.size() == 2, ExcMessage("The mixture is NOT binary. Use another function.") );

    const double M1 = gases[0]->get_molar_mass() * 1.0e3;
    const double M2 = gases[1]->get_molar_mass() * 1.0e3;

    const double sigma12 = 0.5*( gases[0]->get_collision_diameter() + gases[1]->get_collision_diameter() );

    const double T = temp;

    const double omega12 = get_binary_collision_integral(temp);

    return 1.8829e-2 * std::sqrt( ( 1.0/M1 + 1.0/M2 )*T*T*T ) / ( sigma12*sigma12*omega12 );
}

// ---                                                  ---
// --- get_ChapmanEnskog_isobaric_diffusion_coefficient ---
// ---                                                  ---

void
NAME::GasMixture::get_ChapmanEnskog_isobaric_diffusion_coefficient(const std::vector<double>& temp,
                                                                   std::vector<double>&       diffusion_coefficient) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( gases.size() == 2, ExcMessage("The mixture is NOT binary. Use another function.") );
    AssertThrow( temp.size() == diffusion_coefficient.size() , ExcDimensionMismatch(temp.size(), diffusion_coefficient.size()) );

    const double M1 = gases[0]->get_molar_mass() * 1.0e3;
    const double M2 = gases[1]->get_molar_mass() * 1.0e3;

    const double sigma12 = 0.5*( gases[0]->get_collision_diameter() + gases[1]->get_collision_diameter() );

    std::vector<double> omega12;
    omega12.resize(temp.size());

    get_binary_collision_integral(temp, omega12);

    for(unsigned int q = 0; q < temp.size(); ++q)
        diffusion_coefficient[q] = 1.8829e-2 * std::sqrt( ( 1.0/M1 + 1.0/M2 )*temp[q]*temp[q]*temp[q] ) / ( sigma12*sigma12*omega12[q] );
}

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Derivatives of Chapman Enskog isobaric diffusion coefficient - Binary gas mixture only //
  ////////////////////////////////////////////////////////////////////////////////////////////

// ---                                                                ---
// --- get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature ---
// ---                                                                ---

const double
NAME::GasMixture::get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature(const double& temp) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( gases.size() == 2, ExcMessage("The mixture is NOT binary. Use another function.") );

    const double D12         = get_ChapmanEnskog_isobaric_diffusion_coefficient(temp);
    const double T           = temp;
    const double omega12     = get_binary_collision_integral(temp);
    const double Domega12_DT = get_Dbinary_collision_integral_Dtemperature(temp);

    return D12*(1.5/T - Domega12_DT/omega12);
}

// ---                                                                ---
// --- get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature ---
// ---                                                                ---

void
NAME::GasMixture::get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature(const std::vector<double>& temp,
                                                                                 std::vector<double>&       dst) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( gases.size() == 2, ExcMessage("The mixture is NOT binary. Use another function.") );
    AssertThrow( temp.size() == dst.size() , ExcDimensionMismatch(temp.size(), dst.size()) );

    std::vector<double> D12;
    D12.resize(temp.size());

    get_ChapmanEnskog_isobaric_diffusion_coefficient(temp, D12);

    std::vector<double> omega12;
    omega12.resize(temp.size());

    get_binary_collision_integral(temp, omega12);

    std::vector<double> Domega12_DT;
    Domega12_DT.resize(temp.size());

    get_Dbinary_collision_integral_Dtemperature(temp, Domega12_DT);

    for(unsigned int q = 0; q < temp.size(); ++q)
        dst[q] = D12[q]*(1.5/temp[q] - Domega12_DT[q]/omega12[q]);
}

  ////////////////////////////////////////////////////////////////////
  // Chapman Enskog diffusion coefficient - Binary gas mixture only //
  ////////////////////////////////////////////////////////////////////

// ---                                         ---
// --- get_ChapmanEnskog_diffusion_coefficient ---
// ---                                         ---

const double
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficient() const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( gases.size() == 2, ExcMessage("The mixture is NOT binary. Use another function.") );
    AssertThrow( tempIsoTherm != false, ExcMessage("The mixture is NOT isothermal. Use another function.") );
    AssertThrow( pressIsoBaric != false, ExcMessage("The mixture is NOT isobaric. Use another function.") );

    return get_ChapmanEnskog_isobaric_diffusion_coefficient() / total_pressure;
}

// ---                                         ---
// --- get_ChapmanEnskog_diffusion_coefficient ---
// ---                                         ---

void
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficient(std::vector<double>& diffusion_coefficient) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( gases.size() == 2, ExcMessage("The mixture is NOT binary. Use another function.") );
    AssertThrow( tempIsoTherm != false, ExcMessage("The mixture is NOT isothermal. Use another function.") );
    AssertThrow( pressIsoBaric != false, ExcMessage("The mixture is NOT isobaric. Use another function.") );

    const double result = get_ChapmanEnskog_diffusion_coefficient();

    for(unsigned int q = 0; q < diffusion_coefficient.size(); ++q)
        diffusion_coefficient[q] = result;
}

// ---                                                              ---
// --- get_ChapmanEnskog_diffusion_coefficient_at_constant_pressure ---
// ---                                                              ---

const double
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficient_at_constant_pressure(const double& temp) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( gases.size() == 2, ExcMessage("The mixture is NOT binary. Use another function.") );
    AssertThrow( pressIsoBaric != false, ExcMessage("The mixture is NOT isobaric. Use another function.") );

    return get_ChapmanEnskog_isobaric_diffusion_coefficient(temp) / total_pressure;
}

// ---                                                              ---
// --- get_ChapmanEnskog_diffusion_coefficient_at_constant_pressure ---
// ---                                                              ---

void
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficient_at_constant_pressure(const std::vector<double>& temp,
                                                                               std::vector<double>&       diffusion_coefficient) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( gases.size() == 2, ExcMessage("The mixture is NOT binary. Use another function.") );
    AssertThrow( pressIsoBaric != false, ExcMessage("The mixture is NOT isobaric. Use another function.") );
    AssertThrow( temp.size() == diffusion_coefficient.size() , ExcDimensionMismatch(temp.size(), diffusion_coefficient.size()) );

    std::vector<double> D12;
    D12.resize(temp.size());

    get_ChapmanEnskog_isobaric_diffusion_coefficient(temp, D12);

    for(unsigned int q = 0; q < temp.size(); ++q)
        diffusion_coefficient[q] = D12[q] / total_pressure;
}

// ---                                                                 ---
// --- get_ChapmanEnskog_diffusion_coefficient_at_constant_temperature ---
// ---                                                                 ---

const double
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficient_at_constant_temperature(const double& total_pres) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( gases.size() == 2, ExcMessage("The mixture is NOT binary. Use another function.") );
    AssertThrow( tempIsoTherm != false, ExcMessage("The mixture is NOT isothermal. Use another function.") );

    return get_ChapmanEnskog_isobaric_diffusion_coefficient() / total_pres;
}

// ---                                                                 ---
// --- get_ChapmanEnskog_diffusion_coefficient_at_constant_temperature ---
// ---                                                                 ---

void
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficient_at_constant_temperature(const std::vector<double>& total_pres,
                                                                                  std::vector<double>&       diffusion_coefficient) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( gases.size() == 2, ExcMessage("The mixture is NOT binary. Use another function.") );
    AssertThrow( tempIsoTherm != false, ExcMessage("The mixture is NOT isothermal. Use another function.") );
    AssertThrow( total_pres.size() == diffusion_coefficient.size() , ExcDimensionMismatch(total_pres.size(), diffusion_coefficient.size()) );

    std::vector<double> D12;
    D12.resize(total_pres.size());

    get_ChapmanEnskog_isobaric_diffusion_coefficient(D12);

    for(unsigned int q = 0; q < total_pres.size(); ++q)
        diffusion_coefficient[q] = D12[q] / total_pres[q];
}

// ---                                         ---
// --- get_ChapmanEnskog_diffusion_coefficient ---
// ---                                         ---

const double
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficient(const double& total_pres,
                                                          const double& temp) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( gases.size() == 2, ExcMessage("The mixture is NOT binary. Use another function.") );

    return get_ChapmanEnskog_isobaric_diffusion_coefficient(temp) / total_pres;
}

// ---                                         ---
// --- get_ChapmanEnskog_diffusion_coefficient ---
// ---                                         ---

void
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficient(const std::vector<double>& total_pres,
                                                          const std::vector<double>& temp,
                                                          std::vector<double>&       diffusion_coefficient) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( gases.size() == 2, ExcMessage("The mixture is NOT binary. Use another function.") );
    AssertThrow( total_pres.size() == diffusion_coefficient.size() , ExcDimensionMismatch(total_pres.size(), diffusion_coefficient.size()) );
    AssertThrow( temp.size()       == diffusion_coefficient.size() , ExcDimensionMismatch(temp.size(),       diffusion_coefficient.size()) );

    std::vector<double> D12;
    D12.resize(temp.size());

    get_ChapmanEnskog_isobaric_diffusion_coefficient(temp, D12);

    for(unsigned int q = 0; q < temp.size(); ++q)
        diffusion_coefficient[q] = D12[q] / total_pres[q];
}

  ///////////////////////////////////////////////////////////////////////////////////
  // Derivatives of Chapman Enskog diffusion coefficient - Binary gas mixture only //
  ///////////////////////////////////////////////////////////////////////////////////

// ---                                                    ---
// --- get_DChapmanEnskog_diffusion_coefficient_Dpressure ---
// ---                                                    ---

const double
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficient_Dpressure(const double& total_pres) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( gases.size() == 2, ExcMessage("The mixture is NOT binary. Use another function.") );
    AssertThrow( tempIsoTherm != false, ExcMessage("The mixture is NOT isothermal. Use another function.") );

    return - get_ChapmanEnskog_isobaric_diffusion_coefficient() / (total_pres*total_pres);
}

// ---                                                    ---
// --- get_DChapmanEnskog_diffusion_coefficient_Dpressure ---
// ---                                                    ---

void
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficient_Dpressure(const std::vector<double>& total_pres,
                                                                     std::vector<double>&       dst) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( gases.size() == 2, ExcMessage("The mixture is NOT binary. Use another function.") );
    AssertThrow( tempIsoTherm != false, ExcMessage("The mixture is NOT isothermal. Use another function.") );
    AssertThrow( total_pres.size() == dst.size() , ExcDimensionMismatch(total_pres.size(), dst.size()) );

    std::vector<double> D12;
    D12.resize(total_pres.size());

    get_ChapmanEnskog_isobaric_diffusion_coefficient(D12);

    for(unsigned int q = 0; q < total_pres.size(); ++q)
        dst[q] = - D12[q] / (total_pres[q]*total_pres[q]);
}

// ---                                                    ---
// --- get_DChapmanEnskog_diffusion_coefficient_Dpressure ---
// ---                                                    ---

const double
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficient_Dpressure(const double& total_pres,
                                                                     const double& temp) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( gases.size() == 2, ExcMessage("The mixture is NOT binary. Use another function.") );

    return - get_ChapmanEnskog_isobaric_diffusion_coefficient(temp) / (total_pres*total_pres);
}

// ---                                                    ---
// --- get_DChapmanEnskog_diffusion_coefficient_Dpressure ---
// ---                                                    ---

void
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficient_Dpressure(const std::vector<double>& total_pres,
                                                                     const std::vector<double>& temp,
                                                                     std::vector<double>&       dst) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( gases.size() == 2, ExcMessage("The mixture is NOT binary. Use another function.") );
    AssertThrow( total_pres.size() == dst.size() , ExcDimensionMismatch(total_pres.size(), dst.size()) );
    AssertThrow( temp.size()       == dst.size() , ExcDimensionMismatch(temp.size(),       dst.size()) );

    std::vector<double> D12;
    D12.resize(temp.size());

    get_ChapmanEnskog_isobaric_diffusion_coefficient(temp, D12);

    for(unsigned int q = 0; q < temp.size(); ++q)
        dst[q] = - D12[q] / (total_pres[q]*total_pres[q]);
}

// ---                                                       ---
// --- get_DChapmanEnskog_diffusion_coefficient_Dtemperature ---
// ---                                                       ---

const double
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficient_Dtemperature(const double& temp) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( gases.size() == 2, ExcMessage("The mixture is NOT binary. Use another function.") );
    AssertThrow( pressIsoBaric != false, ExcMessage("The mixture is NOT isobaric. Use another function.") );

    return get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature(temp) / total_pressure;
}

// ---                                                       ---
// --- get_DChapmanEnskog_diffusion_coefficient_Dtemperature ---
// ---                                                       ---

void
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficient_Dtemperature(const std::vector<double>& temp,
                                                                        std::vector<double>&       dst) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( gases.size() == 2, ExcMessage("The mixture is NOT binary. Use another function.") );
    AssertThrow( pressIsoBaric != false, ExcMessage("The mixture is NOT isobaric. Use another function.") );
    AssertThrow( temp.size() == dst.size() , ExcDimensionMismatch(temp.size(), dst.size()) );

    std::vector<double> DD12_DT;
    DD12_DT.resize(temp.size());

    get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature(temp, DD12_DT);

    for(unsigned int q = 0; q < temp.size(); ++q)
        dst[q] = DD12_DT[q] / total_pressure;
}

// ---                                                       ---
// --- get_DChapmanEnskog_diffusion_coefficient_Dtemperature ---
// ---                                                       ---

const double
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficient_Dtemperature(const double& total_pres,
                                                                        const double& temp) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( gases.size() == 2, ExcMessage("The mixture is NOT binary. Use another function.") );

    return get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature(temp) / total_pres;
}

// ---                                                       ---
// --- get_DChapmanEnskog_diffusion_coefficient_Dtemperature ---
// ---                                                       ---

void
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficient_Dtemperature(const std::vector<double>& total_pres,
                                                                        const std::vector<double>& temp,
                                                                        std::vector<double>&       dst) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( gases.size() == 2, ExcMessage("The mixture is NOT binary. Use another function.") );
    AssertThrow( total_pres.size() == dst.size() , ExcDimensionMismatch(total_pres.size(), dst.size()) );
    AssertThrow( temp.size()       == dst.size() , ExcDimensionMismatch(temp.size(),       dst.size()) );

    std::vector<double> DD12_DT;
    DD12_DT.resize(temp.size());

    get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature(temp, DD12_DT);

    for(unsigned int q = 0; q < temp.size(); ++q)
        dst[q] = DD12_DT[q] / total_pres[q];
}

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Chapman Enskog isobaric diffusion coefficients - Ternary and more complicated gas mixtures //
  ////////////////////////////////////////////////////////////////////////////////////////////////

// ---                                                   ---
// --- get_ChapmanEnskog_isobaric_diffusion_coefficients ---
// ---                                                   ---

const Table< 2, double >
NAME::GasMixture::get_ChapmanEnskog_isobaric_diffusion_coefficients() const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( tempIsoTherm != false, ExcMessage("The mixture is NOT isothermal. Use another function.") );

    Table< 2, double > result(gases.size(), gases.size());

    for(unsigned int i = 0; i < gases.size(); ++i)
        for(unsigned int j = 0; j < gases.size(); ++j)
            if( i != j )
            {
                const double Mi = gases[i]->get_molar_mass() * 1.0e3;
                const double Mj = gases[j]->get_molar_mass() * 1.0e3;

                const double sigmaij = 0.5*( gases[i]->get_collision_diameter() + gases[j]->get_collision_diameter() );

                const double T = temperature;

                const double omegaij = get_binary_collision_integral(i,j);

                result(i,j) = 1.8829e-2 * std::sqrt( ( 1.0/Mi + 1.0/Mj )*T*T*T ) / ( sigmaij*sigmaij*omegaij );
            }

    return result;
}

// ---                                                   ---
// --- get_ChapmanEnskog_isobaric_diffusion_coefficients ---
// ---                                                   ---

void
NAME::GasMixture::get_ChapmanEnskog_isobaric_diffusion_coefficients(std::vector< Table< 2, double > >& diffusion_coefficients) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( tempIsoTherm != false, ExcMessage("The mixture is NOT isothermal. Use another function.") );

    Table< 2, double > result(gases.size(), gases.size());
    result = get_ChapmanEnskog_isobaric_diffusion_coefficients();

    for(unsigned int q = 0; q < diffusion_coefficients.size(); ++q)
    {
        diffusion_coefficients[q].reinit(gases.size(), gases.size());
        diffusion_coefficients[q] = result;
    }
}

// ---                                                   ---
// --- get_ChapmanEnskog_isobaric_diffusion_coefficients ---
// ---                                                   ---

const Table< 2, double >
NAME::GasMixture::get_ChapmanEnskog_isobaric_diffusion_coefficients(const double& temp) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );

    Table< 2, double > result(gases.size(), gases.size());

    for(unsigned int i = 0; i < gases.size(); ++i)
        for(unsigned int j = 0; j < gases.size(); ++j)
            if( i != j )
            {
                const double Mi = gases[i]->get_molar_mass() * 1.0e3; // [g/mol]
                const double Mj = gases[j]->get_molar_mass() * 1.0e3; // [g/mol]

                const double sigmaij = 0.5*( gases[i]->get_collision_diameter() + gases[j]->get_collision_diameter() ); // [A]

                const double T = temp;

                const double omegaij = get_binary_collision_integral(temp, i,j);

                result(i,j) = 1.8829e-2 * std::sqrt( ( 1.0/Mi + 1.0/Mj )*T*T*T ) / ( sigmaij*sigmaij*omegaij );
            }

    return result;
}

// ---                                                   ---
// --- get_ChapmanEnskog_isobaric_diffusion_coefficients ---
// ---                                                   ---

void
NAME::GasMixture::get_ChapmanEnskog_isobaric_diffusion_coefficients(const std::vector<double>&         temp,
                                                                    std::vector< Table< 2, double > >& diffusion_coefficients) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( temp.size() == diffusion_coefficients.size() , ExcDimensionMismatch(temp.size(), diffusion_coefficients.size()) );

    for(unsigned int q = 0; q < temp.size(); ++q)
    {
        diffusion_coefficients[q].reinit(gases.size(), gases.size());
        diffusion_coefficients[q] = get_ChapmanEnskog_isobaric_diffusion_coefficients(temp[q]);
    }
}

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Derivatives of Chapman Enskog isobaric diffusion coefficients - Ternary and more complicated gas mixtures //
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

// ---                                                                 ---
// --- get_DChapmanEnskog_isobaric_diffusion_coefficients_Dtemperature ---
// ---                                                                 ---

const Table< 2, double >
NAME::GasMixture::get_DChapmanEnskog_isobaric_diffusion_coefficients_Dtemperature(const double& temp) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );

    Table< 2, double > D(gases.size(), gases.size());
    D = get_ChapmanEnskog_isobaric_diffusion_coefficients(temp);

    Table< 2, double > result(gases.size(), gases.size());

    for(unsigned int i = 0; i < gases.size(); ++i)
        for(unsigned int j = 0; j < gases.size(); ++j)
            if( i != j )
            {
                const double T = temp;

                const double omegaij = get_binary_collision_integral(temp, i,j);

                const double Domegaij_DT = get_Dbinary_collision_integral_Dtemperature(temp, i,j);

                result(i,j) = D(i,j)*(1.5/T - Domegaij_DT/omegaij);
            }

    return result;
}

// ---                                                                 ---
// --- get_DChapmanEnskog_isobaric_diffusion_coefficients_Dtemperature ---
// ---                                                                 ---

void
NAME::GasMixture::get_DChapmanEnskog_isobaric_diffusion_coefficients_Dtemperature(const std::vector<double>&         temp,
                                                                                  std::vector< Table< 2, double > >& dst) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( temp.size() == dst.size() , ExcDimensionMismatch(temp.size(), dst.size()) );

    for(unsigned int q = 0; q < temp.size(); ++q)
    {
        dst[q].reinit(gases.size(), gases.size());
        dst[q] = get_DChapmanEnskog_isobaric_diffusion_coefficients_Dtemperature(temp[q]);
    }
}

  ///////////////////////////////////////////////////////////////////////////////////////
  // Chapman Enskog diffusion coefficients - Ternary and more complicated gas mixtures //
  ///////////////////////////////////////////////////////////////////////////////////////

// ---                                          ---
// --- get_ChapmanEnskog_diffusion_coefficients ---
// ---                                          ---

const Table< 2, double >
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficients() const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( tempIsoTherm != false, ExcMessage("The mixture is NOT isothermal. Use another function.") );
    AssertThrow( pressIsoBaric != false, ExcMessage("The mixture is NOT isobaric. Use another function.") );

    Table< 2, double > D(gases.size(), gases.size());
    D = get_ChapmanEnskog_isobaric_diffusion_coefficients();

    Table< 2, double > result(gases.size(), gases.size());

    for(unsigned int i = 0; i < gases.size(); ++i)
        for(unsigned int j = 0; j < gases.size(); ++j)
            if( i != j )
                result(i,j) = D(i,j) / total_pressure;

    return result;
}

// ---                                          ---
// --- get_ChapmanEnskog_diffusion_coefficients ---
// ---                                          ---

void
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficients(std::vector< Table< 2, double > >& diffusion_coefficients) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( tempIsoTherm != false, ExcMessage("The mixture is NOT isothermal. Use another function.") );
    AssertThrow( pressIsoBaric != false, ExcMessage("The mixture is NOT isobaric. Use another function.") );

    Table< 2, double > result(gases.size(), gases.size());
    result = get_ChapmanEnskog_diffusion_coefficients();

    for(unsigned int q = 0; q < diffusion_coefficients.size(); ++q)
    {
        diffusion_coefficients[q].reinit(gases.size(), gases.size());
        diffusion_coefficients[q] = result;
    }
}

// ---                                                               ---
// --- get_ChapmanEnskog_diffusion_coefficients_at_constant_pressure ---
// ---                                                               ---

const Table< 2, double >
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficients_at_constant_pressure(const double& temp) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( pressIsoBaric != false, ExcMessage("The mixture is NOT isobaric. Use another function.") );

    Table< 2, double > D(gases.size(), gases.size());
    D = get_ChapmanEnskog_isobaric_diffusion_coefficients(temp);

    Table< 2, double > result(gases.size(), gases.size());

    for(unsigned int i = 0; i < gases.size(); ++i)
        for(unsigned int j = 0; j < gases.size(); ++j)
            if( i != j )
                result(i,j) = D(i,j) / total_pressure;

    return result;
}

// ---                                                               ---
// --- get_ChapmanEnskog_diffusion_coefficients_at_constant_pressure ---
// ---                                                               ---

void
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficients_at_constant_pressure(const std::vector<double>&         temp,
                                                                                std::vector< Table< 2, double > >& diffusion_coefficients) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( pressIsoBaric != false, ExcMessage("The mixture is NOT isobaric. Use another function.") );
    AssertThrow( temp.size() == diffusion_coefficients.size() , ExcDimensionMismatch(temp.size(), diffusion_coefficients.size()) );

    for(unsigned int q = 0; q < temp.size(); ++q)
    {
        diffusion_coefficients[q].reinit(gases.size(), gases.size());
        diffusion_coefficients[q] = get_ChapmanEnskog_diffusion_coefficients_at_constant_pressure(temp[q]);
    }
}

// ---                                                                  ---
// --- get_ChapmanEnskog_diffusion_coefficients_at_constant_temperature ---
// ---                                                                  ---

const Table< 2, double >
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficients_at_constant_temperature(const double& total_pres) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( tempIsoTherm != false, ExcMessage("The mixture is NOT isothermal. Use another function.") );

    Table< 2, double > D(gases.size(), gases.size());
    D = get_ChapmanEnskog_isobaric_diffusion_coefficients();

    Table< 2, double > result(gases.size(), gases.size());

    for(unsigned int i = 0; i < gases.size(); ++i)
        for(unsigned int j = 0; j < gases.size(); ++j)
            if( i != j )
                result(i,j) = D(i,j) / total_pres;

    return result;
}

// ---                                                                  ---
// --- get_ChapmanEnskog_diffusion_coefficients_at_constant_temperature ---
// ---                                                                  ---

void
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficients_at_constant_temperature(const std::vector<double>&         total_pres,
                                                                                   std::vector< Table< 2, double > >& diffusion_coefficients) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( tempIsoTherm != false, ExcMessage("The mixture is NOT isothermal. Use another function.") );
    AssertThrow( total_pres.size() == diffusion_coefficients.size() , ExcDimensionMismatch(total_pres.size(), diffusion_coefficients.size()) );

    for(unsigned int q = 0; q < total_pres.size(); ++q)
    {
        diffusion_coefficients[q].reinit(gases.size(), gases.size());
        diffusion_coefficients[q] = get_ChapmanEnskog_diffusion_coefficients_at_constant_temperature(total_pres[q]);
    }
}

// ---                                          ---
// --- get_ChapmanEnskog_diffusion_coefficients ---
// ---                                          ---

const Table< 2, double >
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficients(const double& total_pres,
                                                           const double& temp) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    
    Table< 2, double > D(gases.size(), gases.size());
    D = get_ChapmanEnskog_isobaric_diffusion_coefficients(temp);

    Table< 2, double > result(gases.size(), gases.size());

    for(unsigned int i = 0; i < gases.size(); ++i)
        for(unsigned int j = 0; j < gases.size(); ++j)
            if( i != j )
                result(i,j) = D(i,j) / total_pres;

    return result;
}

// ---                                          ---
// --- get_ChapmanEnskog_diffusion_coefficients ---
// ---                                          ---

void
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficients(const std::vector<double>&         total_pres,
                                                           const std::vector<double>&         temp,
                                                           std::vector< Table< 2, double > >& diffusion_coefficients) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( total_pres.size() == diffusion_coefficients.size() , ExcDimensionMismatch(total_pres.size(), diffusion_coefficients.size()) );
    AssertThrow( temp.size()       == diffusion_coefficients.size() , ExcDimensionMismatch(temp.size(),       diffusion_coefficients.size()) );

    for(unsigned int q = 0; q < temp.size(); ++q)
    {
        diffusion_coefficients[q].reinit(gases.size(), gases.size());
        diffusion_coefficients[q] = get_ChapmanEnskog_diffusion_coefficients(total_pres[q], temp[q]);
    }
}

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Derivatives of Chapman Enskog diffusion coefficients - Ternary and more complicated gas mixtures //
  //////////////////////////////////////////////////////////////////////////////////////////////////////

// ---                                                     ---
// --- get_DChapmanEnskog_diffusion_coefficients_Dpressure ---
// ---                                                     ---

const Table< 2, double >
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficients_Dpressure(const double& total_pres) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( tempIsoTherm != false, ExcMessage("The mixture is NOT isothermal. Use another function.") );

    Table< 2, double > D(gases.size(), gases.size());
    D = get_ChapmanEnskog_isobaric_diffusion_coefficients();

    Table< 2, double > result(gases.size(), gases.size());

    for(unsigned int i = 0; i < gases.size(); ++i)
        for(unsigned int j = 0; j < gases.size(); ++j)
            if( i != j )
                result(i,j) = - D(i,j) / (total_pres*total_pres);

    return result;
}

// ---                                                     ---
// --- get_DChapmanEnskog_diffusion_coefficients_Dpressure ---
// ---                                                     ---

void
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficients_Dpressure(const std::vector<double>&         total_pres,
                                                                      std::vector< Table< 2, double > >& dst) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( tempIsoTherm != false, ExcMessage("The mixture is NOT isothermal. Use another function.") );
    AssertThrow( total_pres.size() == dst.size() , ExcDimensionMismatch(total_pres.size(), dst.size()) );

    for(unsigned int q = 0; q < total_pres.size(); ++q)
    {
        dst[q].reinit(gases.size(), gases.size());
        dst[q] = get_DChapmanEnskog_diffusion_coefficients_Dpressure(total_pres[q]);
    }
}

// ---                                                     ---
// --- get_DChapmanEnskog_diffusion_coefficients_Dpressure ---
// ---                                                     ---

const Table< 2, double >
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficients_Dpressure(const double& total_pres,
                                                                      const double& temp) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );

    Table< 2, double > D(gases.size(), gases.size());
    D = get_ChapmanEnskog_isobaric_diffusion_coefficients(temp);

    Table< 2, double > result(gases.size(), gases.size());

    for(unsigned int i = 0; i < gases.size(); ++i)
        for(unsigned int j = 0; j < gases.size(); ++j)
            if( i != j )
                result(i,j) = - D(i,j) / (total_pres*total_pres);

    return result;
}

// ---                                                     ---
// --- get_DChapmanEnskog_diffusion_coefficients_Dpressure ---
// ---                                                     ---

void
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficients_Dpressure(const std::vector<double>&         total_pres,
                                                                      const std::vector<double>&         temp,
                                                                      std::vector< Table< 2, double > >& dst) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( total_pres.size() == dst.size() , ExcDimensionMismatch(total_pres.size(), dst.size()) );
    AssertThrow( temp.size()       == dst.size() , ExcDimensionMismatch(temp.size(),       dst.size()) );

    for(unsigned int q = 0; q < temp.size(); ++q)
    {
        dst[q].reinit(gases.size(), gases.size());
        dst[q] = get_DChapmanEnskog_diffusion_coefficients_Dpressure(total_pres[q], temp[q]);
    }
}

// ---                                                        ---
// --- get_DChapmanEnskog_diffusion_coefficients_Dtemperature ---
// ---                                                        ---

const Table< 2, double >
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficients_Dtemperature(const double& temp) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( pressIsoBaric != false, ExcMessage("The mixture is NOT isobaric. Use another function.") );

    Table< 2, double > DD_DT(gases.size(), gases.size());
    DD_DT = get_DChapmanEnskog_isobaric_diffusion_coefficients_Dtemperature(temp);

    Table< 2, double > result(gases.size(), gases.size());

    for(unsigned int i = 0; i < gases.size(); ++i)
        for(unsigned int j = 0; j < gases.size(); ++j)
            if( i != j )
                result(i,j) = DD_DT(i,j) / total_pressure;

    return result;
}

// ---                                                        ---
// --- get_DChapmanEnskog_diffusion_coefficients_Dtemperature ---
// ---                                                        ---

void
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficients_Dtemperature(const std::vector<double>&         temp,
                                                                         std::vector< Table< 2, double > >& dst) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( pressIsoBaric != false, ExcMessage("The mixture is NOT isobaric. Use another function.") );
    AssertThrow( temp.size() == dst.size() , ExcDimensionMismatch(temp.size(), dst.size()) );

    for(unsigned int q = 0; q < temp.size(); ++q)
    {
        dst[q].reinit(gases.size(), gases.size());
        dst[q] = get_DChapmanEnskog_diffusion_coefficients_Dtemperature(temp[q]);
    }
}

// ---                                                        ---
// --- get_DChapmanEnskog_diffusion_coefficients_Dtemperature ---
// ---                                                        ---

const Table< 2, double >
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficients_Dtemperature(const double& total_pres,
                                                                         const double& temp) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );

    Table< 2, double > DD_DT(gases.size(), gases.size());
    DD_DT = get_DChapmanEnskog_isobaric_diffusion_coefficients_Dtemperature(temp);

    Table< 2, double > result(gases.size(), gases.size());

    for(unsigned int i = 0; i < gases.size(); ++i)
        for(unsigned int j = 0; j < gases.size(); ++j)
            if( i != j )
                result(i,j) = DD_DT(i,j) / total_pres;

    return result;
}

// ---                                                        ---
// --- get_DChapmanEnskog_diffusion_coefficients_Dtemperature ---
// ---                                                        ---

void
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficients_Dtemperature(const std::vector<double>&         total_pres,
                                                                         const std::vector<double>&         temp,
                                                                         std::vector< Table< 2, double > >& dst) const
{
    AssertThrow( !gases.empty(), ExcMessage("std::vector gases does not contain data. Fill it up !") );
    AssertThrow( total_pres.size() == dst.size() , ExcDimensionMismatch(total_pres.size(), dst.size()) );
    AssertThrow( temp.size()       == dst.size() , ExcDimensionMismatch(temp.size(),       dst.size()) );

    for(unsigned int q = 0; q < temp.size(); ++q)
    {
        dst[q].reinit(gases.size(), gases.size());
        dst[q] = get_DChapmanEnskog_diffusion_coefficients_Dtemperature(total_pres[q], temp[q]);
    }
}

  ///////////////////////////////
  // Binary collision integral //
  ///////////////////////////////

// ---                               ---
// --- get_binary_collision_integral ---
// ---                               ---

const double
NAME::GasMixture::get_binary_collision_integral(const unsigned int& N1,
                                                const unsigned int& N2) const
{
    AssertThrow( gases.size() > 1, ExcMessage("The mixture MUST contain at least 2 gases") );
    AssertThrow( std::max(N1, N2) < gases.size(), ExcMessage("The maximum of two last arguments MUST be smaller than the total number of gases") );
    AssertThrow( tempIsoTherm != false, ExcMessage("The mixture is NOT isothermal. Use another function.") );

    const double epsN1N2_BY_k = std::sqrt( gases[N1]->get_eps_BY_k() * gases[N2]->get_eps_BY_k() );

    const double tmp = temperature/epsN1N2_BY_k;

    return Constants::A_diff() / std::pow(tmp, Constants::B_diff())
            +
            Constants::C_diff() / std::exp(Constants::D_diff()*tmp)
            +
            Constants::E_diff() / std::exp(Constants::F_diff()*tmp)
            +
            Constants::G_diff() / std::exp(Constants::H_diff()*tmp);
}

// ---                               ---
// --- get_binary_collision_integral ---
// ---                               ---

void
NAME::GasMixture::get_binary_collision_integral(std::vector<double>& binary_collision_integral,
                                                const unsigned int&  N1,
                                                const unsigned int&  N2) const
{
    AssertThrow( gases.size() > 1, ExcMessage("The mixture MUST contain at least 2 gases") );
    AssertThrow( std::max(N1, N2) < gases.size(), ExcMessage("The maximum of two last arguments MUST be smaller than the total number of gases") );
    AssertThrow( tempIsoTherm != false, ExcMessage("The mixture is NOT isothermal. Use another function.") );

    const double result = get_binary_collision_integral(N1,N2);

    for(unsigned int q = 0; q < binary_collision_integral.size(); ++q)
        binary_collision_integral[q] = result;
}

// ---                               ---
// --- get_binary_collision_integral ---
// ---                               ---

const double
NAME::GasMixture::get_binary_collision_integral(const double&       temp,
                                                const unsigned int& N1,
                                                const unsigned int& N2) const
{
    AssertThrow( gases.size() > 1, ExcMessage("The mixture MUST contain at least 2 gases") );
    AssertThrow( std::max(N1, N2) < gases.size(), ExcMessage("The maximum of two last arguments MUST be smaller than the total number of gases") );

    const double epsN1N2_BY_k = std::sqrt( gases[N1]->get_eps_BY_k() * gases[N2]->get_eps_BY_k() );

    const double tmp = temp/epsN1N2_BY_k;

    return Constants::A_diff() / std::pow(tmp, Constants::B_diff())
            +
            Constants::C_diff() / std::exp(Constants::D_diff()*tmp)
            +
            Constants::E_diff() / std::exp(Constants::F_diff()*tmp)
            +
            Constants::G_diff() / std::exp(Constants::H_diff()*tmp);
}

// ---                               ---
// --- get_binary_collision_integral ---
// ---                               ---

void
NAME::GasMixture::get_binary_collision_integral(const std::vector<double>& temp,
                                                std::vector<double>&       binary_collision_integral,
                                                const unsigned int&        N1,
                                                const unsigned int&        N2) const
{
    AssertThrow( gases.size() > 1, ExcMessage("The mixture MUST contain at least 2 gases") );
    AssertThrow( std::max(N1, N2) < gases.size(), ExcMessage("The maximum of two last arguments MUST be smaller than the total number of gases") );
    AssertThrow( temp.size() == binary_collision_integral.size() , ExcDimensionMismatch(temp.size(), binary_collision_integral.size()) );

    for(unsigned int q = 0; q < temp.size(); ++q)
        binary_collision_integral[q] = get_binary_collision_integral(temp[q], N1, N2);
}

// ---                                             ---
// --- get_Dbinary_collision_integral_Dtemperature ---
// ---                                             ---

const double
NAME::GasMixture::get_Dbinary_collision_integral_Dtemperature(const double&       temp,
                                                              const unsigned int& N1,
                                                              const unsigned int& N2) const
{
    AssertThrow( gases.size() > 1, ExcMessage("The mixture MUST contain at least 2 gases") );
    AssertThrow( std::max(N1, N2) < gases.size(), ExcMessage("The maximum of two last arguments MUST be smaller than the total number of gases") );

    const double epsN1N2_BY_k = std::sqrt( gases[N1]->get_eps_BY_k() * gases[N2]->get_eps_BY_k() );

    const double tmp = temp/epsN1N2_BY_k;

    return - (1.0/epsN1N2_BY_k)*(   Constants::A_diff()*Constants::B_diff()*std::pow( tmp, (-Constants::B_diff()-1.0) )
                                    +
                                    Constants::C_diff()*Constants::D_diff()*std::exp( -Constants::D_diff()*tmp )
                                    +
                                    Constants::E_diff()*Constants::F_diff()*std::exp( -Constants::F_diff()*tmp )
                                    +
                                    Constants::G_diff()*Constants::H_diff()*std::exp( -Constants::H_diff()*tmp )   );
}

// ---                                             ---
// --- get_Dbinary_collision_integral_Dtemperature ---
// ---                                             ---

void
NAME::GasMixture::get_Dbinary_collision_integral_Dtemperature(const std::vector<double>& temp,
                                                              std::vector<double>&       dst,
                                                              const unsigned int&        N1,
                                                              const unsigned int&        N2) const
{
    AssertThrow( gases.size() > 1, ExcMessage("The mixture MUST contain at least 2 gases") );
    AssertThrow( std::max(N1, N2) < gases.size(), ExcMessage("The maximum of two last arguments MUST be smaller than the total number of gases") );
    AssertThrow( temp.size() == dst.size() , ExcDimensionMismatch(temp.size(), dst.size()) );

    for(unsigned int q = 0; q < temp.size(); ++q)
        dst[q] = get_Dbinary_collision_integral_Dtemperature(temp[q],
                                                            N1,N2);
}

////////////////////
// omega integral //
////////////////////

// ---                            ---
// --- get_omega_star_11_integral ---
// ---                            ---
const double
NAME::GasMixture::get_omega_star_11_integral(const double&       temperature,
                                             const unsigned int& N1,
                                             const unsigned int& N2) const
{
    AssertThrow( gases.size() > 0, ExcMessage("The mixture MUST contain at least 1 gas.") );
    AssertThrow( std::max(N1, N2) < gases.size(), ExcMessage("Trying to access element outside of gases vector.") );
    AssertThrow( tempIsoTherm != false, ExcMessage("The mixture is NOT isothermal. Use another function.") );

    const double epsN1N2_BY_k = std::sqrt( gases[N1]->get_eps_BY_k() * gases[N2]->get_eps_BY_k() );
    const double tempStar = temperature/epsN1N2_BY_k;
    
    if( (tempStar >= 0.3) && (tempStar < 2.5) )
    {
        return Constants::A11_visc_R1() * std::pow(tempStar, -Constants::B11_visc_R1())
               +
               Constants::C11_visc_R1() * std::exp(-Constants::D11_visc_R1() * tempStar);
    }
    else if( (tempStar >= 2.5) && (tempStar <= 400) )
    {
        return Constants::A11_visc_R2() * std::pow(tempStar, -Constants::B11_visc_R2())
               +
               Constants::C11_visc_R2() * std::exp(-Constants::D11_visc_R2() * tempStar);
    }
    else
    {
        AssertThrow( false, ExcMessage("T* is outside of range of values Omega*(1,1) equation is valid for.") );
    }
}

void
NAME::GasMixture::get_omega_star_11_integral(const double&        temperature,
                                             std::vector<double>& omega_integral,
                                             const unsigned int&  N1,
                                             const unsigned int&  N2) const
{
    AssertThrow( gases.size() > 0, ExcMessage("The mixture MUST contain at least 1 gas.") );
    AssertThrow( std::max(N1, N2) < gases.size(), ExcMessage("Trying to access element outside of gases vector.") );
    AssertThrow( tempIsoTherm != false, ExcMessage("The mixture is NOT isothermal. Use another function.") );

    const double result = get_omega_star_11_integral(temperature, N1, N2);

    for(unsigned int q = 0; q < omega_integral.size(); ++q)
        omega_integral[q] = result;
}

// ---                            ---
// --- get_omega_star_22_integral ---
// ---                            ---
const double
NAME::GasMixture::get_omega_star_22_integral(const double&       temperature,
                                             const unsigned int& N1,
                                             const unsigned int& N2) const
{
    AssertThrow( gases.size() > 0, ExcMessage("The mixture MUST contain at least 1 gas.") );
    AssertThrow( std::max(N1, N2) < gases.size(), ExcMessage("Trying to access element outside of gases vector.") );
    AssertThrow( tempIsoTherm != false, ExcMessage("The mixture is NOT isothermal. Use another function.") );

    const double epsN1N2_BY_k = std::sqrt( gases[N1]->get_eps_BY_k() * gases[N2]->get_eps_BY_k() );
    const double tempStar = temperature/epsN1N2_BY_k;
    
    if( (tempStar >= 0.3) && (tempStar < 2.5) )
    {
        return Constants::A22_visc_R1() * std::pow(tempStar, -Constants::B22_visc_R1())
               +
               Constants::C22_visc_R1() * std::exp(-Constants::D22_visc_R1() * tempStar);
    }
    else if( (tempStar >= 2.5) && (tempStar <= 400) )
    {
        return Constants::A22_visc_R2() * std::pow(tempStar, -Constants::B22_visc_R2())
               +
               Constants::C22_visc_R2() * std::exp(-Constants::D22_visc_R2() * tempStar);
    }
    else
    {
        AssertThrow( false, ExcMessage("T* is outside of range of values Omega*(2,2) equation is valid for.") );
    }
}

void
NAME::GasMixture::get_omega_star_22_integral(const double&        temperature,
                                             std::vector<double>& omega_integral,
                                             const unsigned int&  N1,
                                             const unsigned int&  N2) const
{
    AssertThrow( gases.size() > 0, ExcMessage("The mixture MUST contain at least 1 gas.") );
    AssertThrow( std::max(N1, N2) < gases.size(), ExcMessage("Trying to access element outside of gases vector.") );
    AssertThrow( tempIsoTherm != false, ExcMessage("The mixture is NOT isothermal. Use another function.") );

    const double result = get_omega_star_22_integral(temperature, N1, N2);

    for(unsigned int q = 0; q < omega_integral.size(); ++q)
        omega_integral[q] = result;
}

//////////////////////////////////
// partial viscosity of mixture //
//////////////////////////////////

// ---                               ---
// --- get_partial_viscosity ---
// ---                               ---
const std::vector< std::vector<double> >
NAME::GasMixture::get_isothermal_nonisobaric_partial_viscosity(const double&                                      tempOfMixture,
                                                               const std::vector< std::vector<double> >&          density,
                                                               const std::vector<double>&                         molarMass,
                                                               const std::vector<double>&                         dynamicViscosity,
                                                               const std::vector<double>&                         collisionDiameter,
                                                               const std::vector<double>&                         porosity,
                                                               std::vector< std::vector< std::vector<double> > >& paramMatrix,
                                                               std::vector< FullMatrix<double> >&                 PInv) const
{
    const unsigned int numOfSpecies    = molarMass.size();
    const unsigned int n_q_points_cell = density[0].size();
    
    std::vector< std::vector<double> > partialViscosity;
    paramMatrix.clear();
    PInv.clear();
    
    if( mixture_viscosity_mode == "Wilke" )
    {
        for(unsigned int q = 0; q < n_q_points_cell; ++q)
        {
            std::vector< std::vector<double> > paramMatrixTemp; //temporary matrix of xi_{ij}
            
            std::vector<double> densityTemp(numOfSpecies, 0.0); //temporary vector for holding density values at specified quadrature point
            for(unsigned int s = 0; s < numOfSpecies; ++s)
                densityTemp[s] = density[s][q];
            
            std::vector<double> temp = get_isothermal_nonisobaric_Wilke_partial_viscosity(densityTemp, dynamicViscosity, molarMass, porosity[q], paramMatrixTemp);
            
            partialViscosity.push_back(temp);
            paramMatrix.push_back(paramMatrixTemp);
        }
    }
    else if( mixture_viscosity_mode == "OmegaKG" )
    {
        std::vector< std::vector<double> > paramMatrixTemp; //temporary matrix for table of Omega integrals
        
        for(unsigned int q = 0; q < n_q_points_cell; ++q)
        {
            std::vector<double> densityTemp(numOfSpecies, 0.0); //temporary vector for holding density values at specified quadrature point
            FullMatrix<double> PInvTemp;                        //temporary FullMatrix for PInv
            for(unsigned int s = 0; s < numOfSpecies; ++s)
                densityTemp[s] = density[s][q];
            
            std::vector<double> temp = get_isothermal_nonisobaric_OmegaKG_partial_viscosity(densityTemp, collisionDiameter, molarMass, tempOfMixture, porosity[q], paramMatrixTemp, PInvTemp);
            
            partialViscosity.push_back(temp);
            PInv.push_back(PInvTemp);
        }
        paramMatrix.push_back(paramMatrixTemp);
    }
    else if( mixture_viscosity_mode == "Dynamic" )
    {
        partialViscosity.resize(n_q_points_cell, std::vector<double>(numOfSpecies, 0.0));
        
        //Now set it to partial viscosity (done this way so calculation only done once and not for each quadrature point)
        for(unsigned int q = 0; q < n_q_points_cell; ++q) //loop through quadrature points
            for(unsigned int s = 0; s < numOfSpecies; ++s) //loop through species
                partialViscosity[q][s] = dynamicViscosity[s]; // [kg/m*s]
    }
    else
    {
        FcstUtilities::log << "Mixture viscosity mode you specified does not exist" << std::endl;
        Assert(false, ExcInternalError());
    }
    return partialViscosity;
}

const std::vector< std::vector<double> >
NAME::GasMixture::get_nonisothermal_nonisobaric_partial_viscosity(const std::vector<double>&                         temperature,
                                                                  const std::vector< std::vector<double> >&          density,
                                                                  const std::vector<double>&                         molarMass,
                                                                  const std::vector< std::vector<double> >&          dynamicViscosity,
                                                                  const std::vector<double>&                         collisionDiameter,
                                                                  const std::vector<double>&                         porosity,
                                                                  std::vector< std::vector< std::vector<double> > >& paramMatrix,
                                                                  std::vector< FullMatrix<double> >&                 PInv) const
{
    
    std::vector< std::vector<double> > partialViscosity;
    
    FcstUtilities::log << "get_nonisothermal_nonisobaric_partial_viscosity() has not been implemented yet, as I don't know how the future developer will want to volume average temperature." << std::endl;
    Assert(false, ExcInternalError());

    return partialViscosity;
}

// ---                             ---
// --- get_Wilke_partial_viscosity ---
// ---                             ---
std::vector<double>
NAME::GasMixture::get_isothermal_nonisobaric_Wilke_partial_viscosity(const std::vector<double>& density,
                                                                     const std::vector<double>& dynamicViscosity,
                                                                     const std::vector<double>& molarMass,
                                                                     const double&              porosity) const
{
    std::vector< std::vector<double> > xi;
    return get_isothermal_nonisobaric_Wilke_partial_viscosity(density, dynamicViscosity, molarMass, porosity, xi);
}

std::vector<double>
NAME::GasMixture::get_isothermal_nonisobaric_Wilke_partial_viscosity(const std::vector<double>&          density,
                                                                     const std::vector<double>&          dynamicViscosity,
                                                                     const std::vector<double>&          molarMass,
                                                                     const double&                       porosity,
                                                                     std::vector< std::vector<double> >& xi) const
{
    //Check that all vectors are the same size, if not return error
    AssertThrow( dynamicViscosity.size() == density.size(),   ExcDimensionMismatch(dynamicViscosity.size(), density.size()) );
    AssertThrow( dynamicViscosity.size() == molarMass.size(), ExcDimensionMismatch(dynamicViscosity.size(), molarMass.size()) );
    
    const unsigned int numOfSpecies = dynamicViscosity.size();
    
    //Clear and resize xi just in case anything is wrong with it
    xi.clear();
    xi.resize(numOfSpecies, std::vector<double>(numOfSpecies, 0.0));

    //Create vector for partial viscosity mixture, initialize size to be the same and set initial values to 0.0
    std::vector<double> partialViscosityMixture(numOfSpecies, 0.0);
    
    //loop through species i; the species component for the partial viscosity of mixture we are calculating
    double porosityInv = std::pow(porosity, -1);
    for(unsigned int i = 0; i < numOfSpecies; ++i)
    {
        double sum = 0.0;
        for(unsigned int j = 0; j < numOfSpecies; ++j) //Loop over all species in mixture
        {
            xi[i][j] = std::pow( 1 + std::pow( dynamicViscosity[i] / dynamicViscosity[j] , 0.5 ) * std::pow( molarMass[j] / molarMass[i], 0.25), 2 )
                       /
                       std::pow( 8 * (1 + molarMass[i] / molarMass[j]), 0.5 );
            sum += density[j] * xi[i][j] / molarMass[j];
        }
        partialViscosityMixture[i] = porosityInv * density[i] * dynamicViscosity[i] / (molarMass[i] * sum);
    }
    return partialViscosityMixture;
}

std::vector<double>
NAME::GasMixture::get_nonisothermal_nonisobaric_Wilke_partial_viscosity(const std::vector<double>&          density,
                                                                        const std::vector<double>&          dynamicViscosity,
                                                                        const std::vector<double>&          molarMass,
                                                                        const double&                       porosity,
                                                                        std::vector< std::vector<double> >& xi) const
{
    std::vector<double> partialViscosityMixture;
    
    FcstUtilities::log << "get_nonisothermal_nonisobaric_Wilke_partial_viscosity() has not been implemented yet, as I don't know how the future developer will want to volume average temperature." << std::endl;
    Assert(false, ExcInternalError());
    
    return partialViscosityMixture;
}

// ---                                ---
// --- get_OmegaKG_partial_viscosity ---
// ---                                ---
std::vector<double>
NAME::GasMixture::get_isothermal_nonisobaric_OmegaKG_partial_viscosity(const std::vector<double>& density,
                                                                       const std::vector<double>& collisionDiameter,
                                                                       const std::vector<double>& molarMass,
                                                                       const double&              temperature,
                                                                       const double&              porosity) const
{
    std::vector< std::vector<double> > omegaIntegralTable;
    FullMatrix<double>                 PInv;
    return get_isothermal_nonisobaric_OmegaKG_partial_viscosity(density, collisionDiameter, molarMass, temperature, porosity, omegaIntegralTable, PInv);
}

std::vector<double>
NAME::GasMixture::get_isothermal_nonisobaric_OmegaKG_partial_viscosity(const std::vector<double>&                density,
                                                                       const std::vector<double>&                collisionDiameter,
                                                                       const std::vector<double>&                molarMass,
                                                                       const double&                             temperature,
                                                                       const double&                             porosity,
                                                                             std::vector< std::vector<double> >& omegaIntegralTable,
                                                                             FullMatrix<double>&                 PInv) const
{
    //Check that all vectors are the same size, if not return error
    AssertThrow( collisionDiameter.size() == density.size(),   ExcDimensionMismatch(collisionDiameter.size(), density.size()) );
    AssertThrow( collisionDiameter.size() == molarMass.size(), ExcDimensionMismatch(collisionDiameter.size(), molarMass.size()) );
    omegaIntegralTable.clear(); //Clear old data as precaution
    
    //Constants
    const double       N_A          = Constants::N_A(); //Avogadro's constant
    const double       k            = Constants::K_SI(); //Boltzmann's constant in SI units
    const unsigned int numOfSpecies = collisionDiameter.size();
    
    //Size such that it stores 2 nxn matrices on top of each other creating a 2n x n matrix; Omega11 is ontop of Omega22
    omegaIntegralTable.resize(numOfSpecies*2, std::vector<double>(numOfSpecies));
    
    //P matrix that will need to be inverted for partial viscosity solution
    PInv.reinit(numOfSpecies, numOfSpecies); //resizes and sets values to 0
    
    //Create vector for partial viscosity mixture, initialize size to be the same and set initial values to 0.0
    std::vector<double> partialViscosityMixture(numOfSpecies, 0.0);
    
    //Set some constants used in calculations
    const double porosityInv = std::pow(porosity, -1);
    const double omegaCoefficient = std::sqrt( 2 * Constants::Pi() * k * temperature );
    const double pCoefficient = 2 / (k * temperature);
    
    //loop through species i & j and solve for Omega^(1,1)_ij and Omega^(2,2)_ij
    for(unsigned int i = 0; i < numOfSpecies; ++i)
    {
        double p_ij_Sum = 0.0;
        for(unsigned int j = 0; j < numOfSpecies; ++j) //Loop over all species in mixture
        {
            double collisionDiameter_ij = 0.0;
            double reducedMass_ij       = 0.0;
            if(j == i)
            {
                collisionDiameter_ij = collisionDiameter[j]; // [m]
                reducedMass_ij       = 0.5 * molarMass[i] / (1000 * N_A); // [kg/molecule]
                
                const double OmegaStar_22 = get_omega_star_22_integral(temperature, i, j);
                omegaIntegralTable[i+numOfSpecies][j] = omegaCoefficient * std::pow(collisionDiameter_ij, 2) / std::sqrt(reducedMass_ij) * OmegaStar_22; // Omega_22 value
            }
            else
            {
                collisionDiameter_ij = 0.5 * (collisionDiameter[i] + collisionDiameter[j]); // [m]
                reducedMass_ij       = molarMass[i] * molarMass[j] / (1000 * N_A * (molarMass[i] + molarMass[j])); // [kg/molecule]
                
                const double OmegaStar_11 = get_omega_star_11_integral(temperature, i, j);
                omegaIntegralTable[i][j] = 0.5 * omegaCoefficient * std::pow(collisionDiameter_ij, 2) / std::sqrt(reducedMass_ij) * OmegaStar_11; //Only need components when i != j; Omega_11 value
                
                const double OmegaStar_22 = get_omega_star_22_integral(temperature, i, j);
                omegaIntegralTable[i+numOfSpecies][j] = omegaCoefficient * std::pow(collisionDiameter_ij, 2) / std::sqrt(reducedMass_ij) * OmegaStar_22; // Omega_22 value
                
                p_ij_Sum += density[j] * std::pow(molarMass[i] + molarMass[j], -2) * (5.0 * molarMass[i] * omegaIntegralTable[i][j] + 1.5 * molarMass[j] * omegaIntegralTable[i+numOfSpecies][j]);
                
                PInv(i, j) = - pCoefficient * 16.0 / 15.0 * molarMass[i] * molarMass[j] * std::pow(molarMass[i] + molarMass[j], -2) * (5.0 * omegaIntegralTable[i][j] - 1.5 * omegaIntegralTable[i+numOfSpecies][j]);
            }
        }
        
        //Now that I have all Omega*_ij for row i can calculate PMatrix_ij for row i
        PInv(i, i) = pCoefficient * ( 0.8 * omegaIntegralTable[i+numOfSpecies][i] + porosityInv * 16.0 * molarMass[i] * p_ij_Sum / (15.0 * density[i]));
    }
    
    PInv.invert(PInv);
    for(unsigned int i = 0; i < numOfSpecies; ++i)
    {
        double sum = 0.0;
        for(unsigned int j = 0; j < numOfSpecies; ++j) //Loop over columns of row i
            sum += PInv(i, j); //Need to take sum since viscosity = PInv * Identity
        partialViscosityMixture[i] = sum * 1.0e1; //Convert kg/m*s to 10 g/cm*s
    }

    return partialViscosityMixture;
}


std::vector<double>
NAME::GasMixture::get_nonisothermal_nonisobaric_OmegaKG_partial_viscosity(const std::vector<double>&                density,
                                                                          const std::vector<double>&                collisionDiameter,
                                                                          const std::vector<double>&                molarMass,
                                                                          const double&                             temperature,
                                                                          const double&                             porosity,
                                                                                std::vector< std::vector<double> >& omegaIntegralTable,
                                                                                FullMatrix<double>&                 PInv) const
{
    std::vector<double> partialViscosityMixture;
    
    FcstUtilities::log << "get_nonisothermal_nonisobaric_OmegaKG_partial_viscosity() has not been implemented yet, as I don't know how the future developer will want to volume average temperature." << std::endl;
    Assert(false, ExcInternalError());
    
    return partialViscosityMixture;
}
///////////////////////////////////////////////
// variation of partial viscosity of mixture //
///////////////////////////////////////////////

// ---                             ---
// --- get_delta_partial_viscosity ---
// ---                             ---
void
NAME::GasMixture::get_isothermal_nonisobaric_delta_partial_viscosity(const double&                                            tempOfMixture,
                                                                     const std::vector< std::vector< std::vector<double> > >& paramMatrix,
                                                                     const std::vector< FullMatrix<double> >&                 PInv,
                                                                     const std::vector<double>&                               porosity,
                                                                     const std::vector<double>&                               molarMass,
                                                                     const std::vector<double>&                               dynamicViscosity,
                                                                     const std::vector<double>&                               collisionDiameter,
                                                                     const std::vector< std::vector< std::vector<double> > >& deltaDensity,
                                                                     const std::vector< std::vector<double> >&                density,
                                                                     std::vector< std::vector< std::vector<double> > >&       deltaPartialViscosity) const
{
    if( mixture_viscosity_mode == "Wilke" )
        get_Wilke_delta_partial_viscosity_wrt_density(paramMatrix, porosity, molarMass, dynamicViscosity, deltaDensity, density, deltaPartialViscosity);
    else if( mixture_viscosity_mode == "OmegaKG" )
        get_OmegaKG_delta_partial_viscosity_wrt_density(tempOfMixture, paramMatrix, PInv, porosity, molarMass, collisionDiameter, deltaDensity, density, deltaPartialViscosity);
    else if( mixture_viscosity_mode == "Dynamic" ) //Value only dependent on temperature and simulation is isotherm so value is 0.0
        get_Null_delta_viscosity(deltaPartialViscosity);
    else
    {
        FcstUtilities::log << "Mixture viscosity mode you specified does not exist" << std::endl;
        Assert(false, ExcInternalError());
    }
}

void
NAME::GasMixture::get_nonisothermal_nonisobaric_delta_partial_viscosity(const std::vector<double>&                               temperature,
                                                                        const std::vector< std::vector< std::vector<double> > >& paramMatrix,
                                                                        const std::vector< FullMatrix<double> >&                 PInv,
                                                                        const std::vector<double>&                               porosity,
                                                                        const std::vector<double>&                               molarMass,
                                                                        const std::vector< std::vector<double> >&                dynamicViscosity,
                                                                        const std::vector<double>&                               collisionDiameter,
                                                                        const std::vector< std::vector< std::vector<double> > >& deltaDensity,
                                                                        const std::vector< std::vector<double> >&                density,
                                                                        std::vector< std::vector< std::vector<double> > >&       deltaPartialViscosity) const
{
    FcstUtilities::log << "get_nonisothermal_nonisobaric_delta_partial_viscosity() has not been implemented yet, as I don't know how the future developer will want to volume average temperature." << std::endl;
    Assert(false, ExcInternalError());
}

// ---                             ---
// --- get_Wilke_delta_partial_viscosity ---
// ---                             ---
void
NAME::GasMixture::get_Wilke_delta_partial_viscosity_wrt_density(const std::vector< std::vector< std::vector<double> > >& xi,
                                                                const std::vector<double>&                               porosity,
                                                                const std::vector<double>&                               molarMass,
                                                                const std::vector<double>&                               dynamicViscosity,
                                                                const std::vector< std::vector< std::vector<double> > >& deltaDensity,
                                                                const std::vector< std::vector<double> >&                density,
                                                                std::vector< std::vector< std::vector<double> > >&       deltaPartialViscosity) const
{
    const unsigned int numOfSpecies  = deltaPartialViscosity.size();
    const unsigned int quadraturePts = deltaPartialViscosity[0].size();
    const unsigned int dofsPerCell   = deltaPartialViscosity[0][0].size();
    
    //Calculate sums
    for(unsigned int s = 0; s < numOfSpecies; ++s)
    {
        for(unsigned int q = 0; q < quadraturePts; ++q)
        {
            const double porosityInv        = std::pow(porosity[q], -1);
            
            double sum = 0.0;
            for(unsigned int s2 = 0; s2 < numOfSpecies; ++s2)
                sum += density[s2][q] * xi[q][s][s2] / molarMass[s2]; //density summation term
            
            for(unsigned int k = 0; k < dofsPerCell; ++k)
            {
                double delta_sum = 0.0;
                for(unsigned int s2 = 0; s2 < numOfSpecies; ++s2)
                    delta_sum += deltaDensity[s2][q][k] * xi[q][s][s2] / molarMass[s2]; //delta density summation term
    
                deltaPartialViscosity[s][q][k] = porosityInv * dynamicViscosity[s] / molarMass[s]
                                                 *
                                                 ( deltaDensity[s][q][k] / sum - density[s][q] / std::pow(sum, 2) * delta_sum ); // [1/m*s]
            }
        }
    }
}

void
NAME::GasMixture::get_Wilke_delta_partial_viscosity_wrt_temperature(const std::vector< std::vector< std::vector<double> > >& xi,
                                                                    const std::vector<double>&                               porosity,
                                                                    const std::vector<double>&                               molarMass,
                                                                    const std::vector<double>&                               dynamicViscosity,
                                                                    const std::vector< std::vector< std::vector<double> > >& deltaTemperature,
                                                                    const std::vector< std::vector<double> >&                density,
                                                                    std::vector< std::vector< std::vector<double> > >&       deltaPartialViscosity) const
{
    FcstUtilities::log << "get_Wilke_delta_partial_viscosity_wrt_temperature() has not been implemented yet, as I don't know how the future developer will want to volume average temperature." << std::endl;
    Assert(false, ExcInternalError());
}

// ---                                     ---
// --- get_OmegaKG_delta_partial_viscosity ---
// ---                                     ---
void
NAME::GasMixture::get_OmegaKG_delta_partial_viscosity_wrt_density(const double&                                            tempOfMixture,
                                                                  const std::vector< std::vector< std::vector<double> > >& omegaIntegralTable,
                                                                  const std::vector< FullMatrix<double> >&                 PInv,
                                                                  const std::vector<double>&                               porosity,
                                                                  const std::vector<double>&                               molarMass,
                                                                  const std::vector<double>&                               collisionDiameter,
                                                                  const std::vector< std::vector< std::vector<double> > >& deltaDensity,
                                                                  const std::vector< std::vector<double> >&                density,
                                                                  std::vector< std::vector< std::vector<double> > >&       deltaPartialViscosity) const
{
    const unsigned int numOfSpecies  = deltaPartialViscosity.size();
    const unsigned int quadraturePts = deltaPartialViscosity[0].size();
    const unsigned int dofsPerCell   = deltaPartialViscosity[0][0].size();
    
    const double k                = Constants::K_SI(); //Boltzmann's constant in SI units
    const double N_A              = Constants::N_A();  //Avogadro's constant
    const double pcoefficient     = 32.0 / (15.0 * k * tempOfMixture);
    const double omegaCoefficient = std::sqrt( 2 * Constants::Pi() * k * tempOfMixture );
    
    //NOTE: Since it is assumed that the simulation is isotermal sumCoefficientMatrix can be done independent of the quadrature points. 
    //      But if we ever implement a non-isothermal model, this will have to be placed inside loop below.
    //Coefficients used in summation is the same for all quadrature points and dofs. So it save computations calculate 
    //the matrix of these coefficients first.
    std::vector< std::vector<double> > sumCoefficientMatrix(numOfSpecies, std::vector<double>(numOfSpecies, 0.0));
    for(unsigned int i = 0; i < numOfSpecies; ++i)
        for(unsigned int j = 0; j < numOfSpecies; ++j)
            if(i != j)
                sumCoefficientMatrix[i][j] = std::pow(molarMass[i] + molarMass[j], -2) * (5.0 * molarMass[i] * omegaIntegralTable[0][i][j] + 1.5 * molarMass[j] * omegaIntegralTable[0][i+numOfSpecies][j]);
    
    for(unsigned int q = 0; q < quadraturePts; ++q)
    {
        for(unsigned int k = 0; k < dofsPerCell; ++k)
        {
            const double porosityInv        = std::pow(porosity[q], -1);
            
            //Loop through species, calculate deltaP_ii to construct deltaP matrix
            std::vector<double> deltaP(numOfSpecies, 0.0); // Since only need main diagonal components of deltaP_ij and it is multiplied by (1 vector) only need a vector
            for(unsigned int s = 0; s < numOfSpecies; ++s)
            {
                const double densityInv      = std::pow(density[s][q], -1);
                const double p_MMCoefficient = porosityInv * pcoefficient * molarMass[s] * densityInv;
                
                double densitySum      = 0.0;
                double deltaDensitySum = 0.0;
                for(unsigned int s2 = 0; s2 < numOfSpecies; ++s2)
                {
                    //NOTE: Do not need an if statement to jump over when i +j because sumCoefficientMatrix stores 0.0 in this indice
                    densitySum      += density[s2][q] * sumCoefficientMatrix[s][s2];
                    deltaDensitySum += deltaDensity[s2][q][k] * sumCoefficientMatrix[s][s2];
                }
                deltaP[s] = p_MMCoefficient * (deltaDensitySum - densityInv * deltaDensity[s][q][k] * densitySum);
            }
            
            //Solve for deltaViscosity = P^-1 * deltaP * P^-1 * (1 vector)
            for(unsigned int i = 0; i < numOfSpecies; ++i) //Pass through species of deltaPartialViscosity
            {
                std::vector<double> PInvTimesPii(numOfSpecies, 0.0);
                for(unsigned int s = 0; s < numOfSpecies; ++s) //Pass through row i of PInv * deltaP and produce vector of values
                    PInvTimesPii[s] = PInv[q](i, s) * deltaP[s];
                
                double pSum = 0.0;
                for(unsigned int s1 = 0; s1 < numOfSpecies; ++s1) //Pass through elements of PInvTimesPii
                {
                    double sum = 0.0;
                    for(unsigned int s2 = 0; s2 < numOfSpecies; ++s2) //Pass through columns of s2 in row s2 and sum elements of PInv
                        sum += PInv[q](s1, s2);
                    pSum += sum * PInvTimesPii[s1];
                }
                
                deltaPartialViscosity[i][q][k] = pSum * 1.0e1; //Convert kg/m*s to 10 g/cm*s
            }
        }
    }
}

void
NAME::GasMixture::get_OmegaKG_delta_partial_viscosity_wrt_temperature(const std::vector<double>&                               temperature,
                                                                      const std::vector< std::vector< std::vector<double> > >& omegaIntegralTable,
                                                                      const std::vector< FullMatrix<double> >&                 PInv,
                                                                      const std::vector<double>&                               porosity,
                                                                      const std::vector<double>&                               molarMass,
                                                                      const std::vector<double>&                               collisionDiameter,
                                                                      const std::vector< std::vector< std::vector<double> > >& deltaDensity,
                                                                      const std::vector< std::vector<double> >&                density,
                                                                      std::vector< std::vector< std::vector<double> > >&       deltaPartialViscosity) const
{
    FcstUtilities::log << "get_OmegaKG_delta_partial_viscosity_wrt_temperature() has not been implemented yet, as I don't know how the future developer will want to volume average temperature." << std::endl;
    Assert(false, ExcInternalError());
}

// ---                             ---
// --- get_Null_delta_partial_viscosity ---
// ---                             ---
void
NAME::GasMixture::get_Null_delta_viscosity(std::vector< std::vector< std::vector<double> > >& deltaViscosity) const
{
    const unsigned int numOfSpecies  = deltaViscosity.size();
    const unsigned int quadraturePts = deltaViscosity[0].size();
    const unsigned int dofsPerCell   = deltaViscosity[0][0].size();
    for(unsigned int s = 0; s < numOfSpecies; ++s)
        for(unsigned int q = 0; q < quadraturePts; ++q)
            for(unsigned int k = 0; k < dofsPerCell; ++k)
                deltaViscosity[s][q][k] = 0.0;
}

/////////////////////////////////
// variation of bulk viscosity //
/////////////////////////////////

// ---                          ---
// --- get_delta_bulk_viscosity ---
// ---                          ---
void
NAME::GasMixture::get_delta_bulk_viscosity(const std::vector< PureGas* >&                           gases,
                                           const std::vector< std::vector< std::vector<double> > >& deltaPartialViscosity,
                                                 std::vector< std::vector< std::vector<double> > >& deltaBulkViscosity) const
{
    if( mixture_viscosity_mode == "Wilke" )
    {
        AssertThrow( deltaBulkViscosity.size()       == deltaPartialViscosity.size(),       ExcDimensionMismatch(deltaBulkViscosity.size(),       deltaPartialViscosity.size()) );
        AssertThrow( deltaBulkViscosity[0].size()    == deltaPartialViscosity[0].size(),    ExcDimensionMismatch(deltaBulkViscosity[0].size(),    deltaPartialViscosity[0].size()) );
        AssertThrow( deltaBulkViscosity[0][0].size() == deltaPartialViscosity[0][0].size(), ExcDimensionMismatch(deltaBulkViscosity[0][0].size(), deltaPartialViscosity[0][0].size()) );
        
        const unsigned int numOfSpecies  = deltaBulkViscosity.size();
        const unsigned int quadraturePts = deltaBulkViscosity[0].size();
        const unsigned int dofsPerCell   = deltaBulkViscosity[0][0].size();
        for(unsigned int s = 0; s < numOfSpecies; ++s)
            for(unsigned int q = 0; q < quadraturePts; ++q)
                for(unsigned int k = 0; k < dofsPerCell; ++k)
                    deltaBulkViscosity[s][q][k] = gases[s]->get_bulk_viscosity(deltaPartialViscosity[s][q][k]);
    }
    else if( mixture_viscosity_mode == "OmegaKG" )
    {
        AssertThrow( deltaBulkViscosity.size()       == deltaPartialViscosity.size(),       ExcDimensionMismatch(deltaBulkViscosity.size(),       deltaPartialViscosity.size()) );
        AssertThrow( deltaBulkViscosity[0].size()    == deltaPartialViscosity[0].size(),    ExcDimensionMismatch(deltaBulkViscosity[0].size(),    deltaPartialViscosity[0].size()) );
        AssertThrow( deltaBulkViscosity[0][0].size() == deltaPartialViscosity[0][0].size(), ExcDimensionMismatch(deltaBulkViscosity[0][0].size(), deltaPartialViscosity[0][0].size()) );
        
        const unsigned int numOfSpecies  = deltaBulkViscosity.size();
        const unsigned int quadraturePts = deltaBulkViscosity[0].size();
        const unsigned int dofsPerCell   = deltaBulkViscosity[0][0].size();
        for(unsigned int s = 0; s < numOfSpecies; ++s)
            for(unsigned int q = 0; q < quadraturePts; ++q)
                for(unsigned int k = 0; k < dofsPerCell; ++k)
                    deltaBulkViscosity[s][q][k] = gases[s]->get_bulk_viscosity(deltaPartialViscosity[s][q][k]);
    }
    else if( mixture_viscosity_mode == "Dynamic" ) //Value only dependent on temperature and simulation is isotherm so value is 0.0
        get_Null_delta_viscosity(deltaBulkViscosity);
    else
    {
        FcstUtilities::log << "Mixture viscosity mode you specified does not exist" << std::endl;
        Assert(false, ExcInternalError());
    }
}

////////////////////////
// ACCESSORS AND INFO //
////////////////////////

// ---                           ---
// --- print_material_properties ---
// ---                           ---

void
NAME::GasMixture::print_material_properties() const
{
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "------------------------------";
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;

    FcstUtilities::log << "Parameters for pure gases :";
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;

    for(unsigned int g = 0; g < gases.size(); ++g)
    {
        FcstUtilities::log << gases[g]->name_material() << " :";
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "ID = " << gases[g]->get_ID();
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Chemical formula = " << gases[g]->get_chemical_formula();
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Molar mass, [kg/mol] = " << gases[g]->get_molar_mass();
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Collision diameter, [A] = " << gases[g]->get_collision_diameter();
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "The maximum energy of attraction divided by the Boltzmann constant, [K] = " << gases[g]->get_eps_BY_k();
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Prandtl number = " << gases[g]->get_Prandtl();
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Dynamic viscosity mode = " << gases[g]->get_dynamic_viscosity_mode();
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Bulk viscosity mode = " << gases[g]->get_bulk_viscosity_mode();
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Thermal conductivity mode = " << gases[g]->get_thermal_conductivity_mode();
        FcstUtilities::log << std::endl;

        if( tempIsoTherm == false )
        {
            FcstUtilities::log << "Dynamic viscosity, [Pa sec] = " << gases[g]->get_dynamic_viscosity(temperature);
            FcstUtilities::log << std::endl;

            FcstUtilities::log << "Bulk viscosity, [Pa sec] = " << gases[g]->get_bulk_viscosity(   gases[g]->get_dynamic_viscosity(temperature)   );
            FcstUtilities::log << std::endl;

            FcstUtilities::log << "Thermal conductivity, [W/(m K)] = " << gases[g]->get_thermal_conductivity(temperature);
            FcstUtilities::log << std::endl;

            FcstUtilities::log << "Specific heat capacity at constant pressure, [J/(kg K)] = " << gases[g]->get_specific_heat_capacity(temperature);
            FcstUtilities::log << std::endl;

            FcstUtilities::log << "Molar enthalpy, [J/mol] = " << gases[g]->get_molar_enthalpy(temperature);
            FcstUtilities::log << std::endl;

            if( gases[g]->get_chemical_formula() == "H2O Vapor" )
            {
                FcstUtilities::log << "Water vapor saturation pressure, [Pa] = " << gases[g]->get_water_vapor_saturation_pressure(temperature);
                FcstUtilities::log << std::endl;
            }
        }

        FcstUtilities::log << std::endl;
    }

    FcstUtilities::log << "Parameters for the gas mixture :";
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;

    FcstUtilities::log << this->name << " :";
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;

    if( tempIsoTherm == true && pressIsoBaric == true )
        FcstUtilities::log << "The gas mixture is both ISOTHERMAL and ISOBARIC" << std::endl;

    if( tempIsoTherm == false && pressIsoBaric == true )
        FcstUtilities::log << "The gas mixture is both NON-ISOTHERMAL and ISOBARIC" << std::endl;

    if( tempIsoTherm == true && pressIsoBaric == false )
        FcstUtilities::log << "The gas mixture is both ISOTHERMAL and NON-ISOBARIC" << std::endl;

    if( tempIsoTherm == false && pressIsoBaric == false )
        FcstUtilities::log << "The gas mixture is both NON-ISOTHERMAL and NON-ISOBARIC" << std::endl;

    if( tempIsoTherm == false )
    {
        FcstUtilities::log << "Chapman Enskog isobaric diffusion coefficients, [Pa m^2/sec] :";
        FcstUtilities::log << std::endl;

        Table< 2, double > D(gases.size(), gases.size());
        D = get_ChapmanEnskog_isobaric_diffusion_coefficients();

        for(unsigned int i = 0; i < gases.size(); ++i)
        {
            for(unsigned int j = 0; j < gases.size(); ++j)
                FcstUtilities::log << D(i,j) << "   ";
            FcstUtilities::log << std::endl;
        }

        if( pressIsoBaric == false )
        {
            FcstUtilities::log << "Chapman Enskog diffusion coefficients, [m^2/sec] :";
            FcstUtilities::log << std::endl;

            D = get_ChapmanEnskog_diffusion_coefficients();

            for(unsigned int i = 0; i < gases.size(); ++i)
            {
                for(unsigned int j = 0; j < gases.size(); ++j)
                    FcstUtilities::log << D(i,j) << "   ";
                FcstUtilities::log << std::endl;
            }
        }
    }
}