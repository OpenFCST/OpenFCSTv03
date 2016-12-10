//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: catalyst_layer.cc
//    - Description: Base Catalyst Layer Class. It implements the interface for other catalyst layer class
//        and some common methods.
//    - Developers: Marc Secanell (2011-2013) and Madhur Bhaiya (2013)
//    - $Id: catalyst_layer.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include <layers/catalyst_layer.h>

namespace NAME = FuelCellShop::Layer;

//---------------------------------------------------------------------------
template <int dim>
NAME::CatalystLayer<dim>::CatalystLayer(const std::string& name)
: NAME::PorousLayer<dim>(name)
{
    this->reactant = nothing;
    //this->constant_solutions[temperature_of_REV] = 0.;     // For debug checking purposes
    this->constant_solutions[total_pressure] = 0.;     // For debug checking purposes
    electrolyte = boost::shared_ptr<FuelCellShop::Material::PolymerElectrolyteBase > ();
    catalyst_support = boost::shared_ptr< FuelCellShop::Material::CatalystSupportBase > ();
    catalyst = boost::shared_ptr< FuelCellShop::Material::CatalystBase > ();
    default_materials = true;
    kinetics = boost::shared_ptr< FuelCellShop::Kinetics::BaseKinetics > ();
}

//---------------------------------------------------------------------------
template <int dim>
NAME::CatalystLayer<dim>::CatalystLayer()
: NAME::PorousLayer<dim> ()
{
    this->reactant = nothing;
    //this->constant_solutions[temperature_of_REV] = 0.;     // For debug checking purposes
    this->constant_solutions[total_pressure] = 0.;     // For debug checking purposes
}

//---------------------------------------------------------------------------
template <int dim>
NAME::CatalystLayer<dim>::~CatalystLayer()
{
    // No need to delete boost pointers since boost manages memory allocation
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::CatalystLayer<dim>::declare_parameters (const std::string& name,
                                              ParameterHandler &param) const
{

    FuelCellShop::Layer::PorousLayer<dim>::declare_parameters(name,param);

    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection(name);
        {
            param.declare_entry("Catalyst layer type",
                                "DummyCL",
                                Patterns::Selection("DummyCL | HomogeneousCL | MultiScaleCL"),
                                " ");
            param.declare_entry("Catalyst type",
                                "Platinum",
                                Patterns::Selection("Platinum"),
                                " ");
            param.declare_entry("Catalyst support type",
                                "CarbonBlack",
                                Patterns::Selection("CarbonBlack"),
                                " ");
            param.declare_entry("Electrolyte type",
                                "Nafion",
                                Patterns::Selection("Nafion"),
                                " ");
            param.declare_entry("Kinetics type",
                                "TafelKinetics",
                                Patterns::Selection("TafelKinetics | ButlerVolmerKinetics | DoubleTrapKinetics | DualPathKinetics"),
                                " ");   
            
  
            // Note: This must be called within the section of the CL
            FuelCellShop::Material::CatalystBase::declare_Catalyst_parameters(param);
            FuelCellShop::Material::CatalystSupportBase::declare_CatalystSupport_parameters(param);
            FuelCellShop::Material::PolymerElectrolyteBase::declare_PolymerElectrolyte_parameters(param);
            FuelCellShop::Kinetics::BaseKinetics::declare_Kinetics_parameters(param);
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::CatalystLayer<dim>::initialize (ParameterHandler &param)
{
    NAME::PorousLayer<dim>::initialize(param);

    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection(this->name);
        {
            catalyst_type = param.get("Catalyst type");
            catalyst_support_type = param.get("Catalyst support type");
            electrolyte_type = param.get("Electrolyte type");
            kinetics_type = param.get("Kinetics type");
            
            catalyst = FuelCellShop::Material::CatalystBase::create_Catalyst(param, catalyst_type);
            catalyst_support = FuelCellShop::Material::CatalystSupportBase::create_CatalystSupport(param, catalyst_support_type);
            electrolyte = FuelCellShop::Material::PolymerElectrolyteBase::create_PolymerElectrolyte(param, electrolyte_type);
            kinetics = FuelCellShop::Kinetics::BaseKinetics::create_Kinetics(param, kinetics_type);

            kinetics->set_catalyst(catalyst.get());
            kinetics->set_electrolyte(electrolyte.get());
        }
        param.leave_subsection();
    }
    param.leave_subsection();


}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::CatalystLayer<dim>::set_solution(const std::vector< SolutionVariable >& sols)
{
    Assert( std::find_if(sols.begin(), sols.end(), FuelCellShop::is_phiS) != sols.end(), ExcMessage("VariableNames::electronic_electrical_potential should exist in input vector to CatalystLayer::set_solution.") );
    Assert( std::find_if(sols.begin(), sols.end(), FuelCellShop::is_phiM) != sols.end(), ExcMessage("VariableNames::protonic_electrical_potential should exist in input vector to CatalystLayer::set_solution.") );

    std::vector<SolutionVariable> reactant_concentrations;
    
    //-- Check if the problem is isothermal or non-isothermal and extract index for T:
    int index_T = -1;
    for (unsigned int s=0; s < sols.size(); ++s) {
        if (sols[s].get_variablename() == temperature_of_REV){
            index_T = s;
            break;    
        }    
    } 

    //Ensure that we only compute with "current" solutions
    this->solutions.clear();

    for (unsigned int i=0; i < sols.size(); ++i)
    {
        if (sols[i].get_variablename() == electronic_electrical_potential)
        {
            this->kinetics->set_solid_potential(sols.at(i));
            this->solutions[electronic_electrical_potential] = sols.at(i);
        }

        else if (sols[i].get_variablename() == protonic_electrical_potential)
        {
            this->kinetics->set_electrolyte_potential(sols.at(i));
            this->solutions[protonic_electrical_potential] = sols.at(i);
        }

        else if( sols[i].get_variablename() == hydrogen_concentration )
        {
               AssertThrow( this->electrolyte->get_H_H2() != 0, ExcMessage("Henry's constant for hydrogen not initialized in the electrolyte object of the catalyst layer.") );

               std::vector<double> hydrogen_concentration_old_GL(sols[i].size());
               
               if (index_T >= 0) // Non-isothermal
                   for(unsigned int q = 0; q < hydrogen_concentration_old_GL.size(); ++q)
                       hydrogen_concentration_old_GL[q] = Constants::R()*sols.at(index_T)[q]*sols.at(i)[q]/(this->electrolyte->get_H_H2()*1.0e-6);
               else
               {
                   AssertThrow( this->constant_solutions.at(temperature_of_REV) != 0., ExcMessage("Temperature not initialized in the catalyst layer using set_T method.") );
                   for(unsigned int q = 0; q < hydrogen_concentration_old_GL.size(); ++q)
                       hydrogen_concentration_old_GL[q] = Constants::R()*this->constant_solutions.at(temperature_of_REV)*sols.at(i)[q]/(this->electrolyte->get_H_H2()*1.0e-6);
               }

               this->solutions[VariableNames::hydrogen_concentration] = FuelCellShop::SolutionVariable( hydrogen_concentration_old_GL,
                                                                                                        VariableNames::hydrogen_concentration);
               this->reactant = VariableNames::hydrogen_concentration;
               reactant_concentrations.push_back( FuelCellShop::SolutionVariable( hydrogen_concentration_old_GL,
                                                                                  VariableNames::hydrogen_concentration) );
        }

        else if( sols[i].get_variablename() == oxygen_concentration )
        {

            AssertThrow( this->electrolyte->get_H_O2() != 0, ExcMessage("Henry's constant for oxygen not initialized in the electrolyte object of the catalyst layer.") );
            
            std::vector<double> oxygen_concentration_old_GL(sols[i].size());
            
            if (index_T >= 0) // Non-isothermal
            {
                for(unsigned int q = 0; q < oxygen_concentration_old_GL.size(); ++q)
                    oxygen_concentration_old_GL[q] = Constants::R()*sols.at(index_T)[q]*sols.at(i)[q]/(this->electrolyte->get_H_O2()*1.0e-6);
            }
            else
            {
                AssertThrow( this->constant_solutions.at(temperature_of_REV) != 0., ExcMessage("Temperature not initialized in the catalyst layer using set_T method.") );
                for(unsigned int q = 0; q < oxygen_concentration_old_GL.size(); ++q)
                    oxygen_concentration_old_GL[q] = Constants::R()*this->constant_solutions.at(temperature_of_REV)*sols.at(i)[q]/(this->electrolyte->get_H_O2()*1.0e-6);
            }
            
            this->solutions[VariableNames::oxygen_concentration] = FuelCellShop::SolutionVariable( oxygen_concentration_old_GL,
                                                                                                   VariableNames::oxygen_concentration);
            this->reactant  = VariableNames::oxygen_concentration;            
            reactant_concentrations.push_back( FuelCellShop::SolutionVariable( oxygen_concentration_old_GL, VariableNames::oxygen_concentration) );
        }

        else if (sols[i].get_variablename() == temperature_of_REV)
        {
            this->kinetics->set_temperature(sols.at(i));
            this->solutions[temperature_of_REV] = sols.at(i);
        }

        else if (sols[i].get_variablename() == membrane_water_content)
        {
            this->solutions[membrane_water_content] = sols.at(i);
        }

        else if (sols[i].get_variablename() == oxygen_molar_fraction)
        {
            Assert( this->constant_solutions.at(total_pressure) != 0., ExcMessage("Total pressure not initialized in the catalyst layer using set_P method.") );
            AssertThrow( this->electrolyte->get_H_O2() != 0, ExcMessage("Henry's constant for oxygen not initialized in the electrolyte object of the catalyst layer.") );

            this->solutions[oxygen_molar_fraction] = sols.at(i);
            this->reactant = oxygen_molar_fraction;

            std::vector<double> oxy(sols[i].size());
            for (unsigned int q=0; q<oxy.size(); ++q)
                oxy[q] = sols.at(i)[q] * (this->constant_solutions.at(total_pressure)/this->electrolyte->get_H_O2());

            reactant_concentrations.push_back( SolutionVariable(oxy, oxygen_concentration) );
        }

        else if (sols[i].get_variablename() == hydrogen_molar_fraction)
        {
            Assert( this->constant_solutions.at(total_pressure) != 0., ExcMessage("Total pressure not initialized in the catalyst layer using set_P method.") );
            AssertThrow( this->electrolyte->get_H_H2() != 0, ExcMessage("Henry's constant for hydrogen not initialized in the electrolyte object of the catalyst layer.") );

            this->solutions[hydrogen_molar_fraction] = sols.at(i);
            this->reactant = hydrogen_molar_fraction;

            std::vector<double> hyd(sols[i].size());
            for (unsigned int q=0; q<hyd.size(); ++q)
                hyd[q] = sols.at(i)[q] * (this->constant_solutions.at(total_pressure)/this->electrolyte->get_H_H2());

            reactant_concentrations.push_back( SolutionVariable(hyd, hydrogen_concentration) );
        }
    }

    Assert( reactant_concentrations.size() > 0, ExcMessage("At least one reactant concentration/molar fraction is required in input vector to CatalystLayer::set_solution.") );

    this->kinetics->set_reactant_concentrations(reactant_concentrations);

    this->n_quad = this->solutions[protonic_electrical_potential].size();

    //Check for temperature, initialize if needed be
    if (this->solutions.find(temperature_of_REV) == this->solutions.end())
    {
        Assert( this->constant_solutions.at(temperature_of_REV) != 0., ExcMessage("Temperature not initialized in the catalyst layer using set_T method for isothermal case.") );
        this->solutions[temperature_of_REV] = SolutionVariable(this->constant_solutions.at(temperature_of_REV), this->n_quad, temperature_of_REV);
        this->kinetics->set_temperature( this->solutions[temperature_of_REV] );
    }

    //Check for lambda, initialize if needed be
    if (this->solutions.find(membrane_water_content) == this->solutions.end())
    {
        this->solutions[membrane_water_content] = SolutionVariable(12.0, this->n_quad, membrane_water_content);
    }


    //If anything that was being supplied before is no longer being supplied,
    #ifdef DEBUG

    std::vector<VariableNames> current_names;
    for(unsigned int j =0; j < sols.size(); j++)
        current_names.push_back(sols[j].get_variablename());
/*
    Assert( std::is_permutation(common_names.begin(), common_names.end(), current_names.begin()), ExcMessage("Inconsistent provision of solutions to catalyst layer from application, "
            "be consistent and always provide the same set of solutions"));
    common_names = current_names;
*/
    #endif
}

//---------------------------------------------------------------------------
template <int dim>
FuelCellShop::SolutionMap
NAME::CatalystLayer<dim>::get_coverages(){

    SolutionMap sols;

    if(this->kinetics->has_coverage(OH_coverage)){
        std::vector<double> OH_c;
        this->kinetics->OH_coverage(OH_c);
        sols.push_back(SolutionVariable(OH_c, OH_coverage));
    }
    else
        sols.push_back(SolutionVariable(0.0, this->solutions.size(), OH_coverage));



    if(this->kinetics->has_coverage(O_coverage)){
        std::vector<double> O_c;
        this->kinetics->O_coverage(O_c);
        sols.push_back(SolutionVariable(O_c, O_coverage));
    }
    else
        sols.push_back(SolutionVariable(0.0, this->solutions.size(), O_coverage));


    return sols;

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations.
template class NAME::CatalystLayer<deal_II_dimension>;
