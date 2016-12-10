//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2011-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: agglomerate_base.cc
//    - Description: Base class for agglomerates
//    - Developers: Philip Wardlaw, Peter Dobson, Michael Moore, M. Secanell
//    - $Id: agglomerate_base.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------


#include <microscale/agglomerate_base.h>

namespace NAME = FuelCellShop::MicroScale;


void
NAME::AgglomerateBase::set_solution(const std::map<VariableNames,SolutionVariable>& sols,const VariableNames& react, const int& solution_index){


    typedef std::map<VariableNames,SolutionVariable>::const_iterator map_itr;

    //Check for key variables (x_O2/x_H2, phi_M, and phi_S)
    AssertThrow((react == oxygen_molar_fraction) or (react == hydrogen_molar_fraction) ,
            ExcMessage("Agglomerate cannot solve for this type of reactant."));

    AssertThrow(sols.find(protonic_electrical_potential) != sols.end(),
                ExcMessage("Solution is missing protonic potential!"));

    AssertThrow((sols.find(electronic_electrical_potential) != sols.end()),
                    ExcMessage("Solution is missing solid potential!"));

    this->sol_index = solution_index;
    this->solutions = sols;

    this->sol_names.clear();
    this->sol_names.push_back(react);
    this->sol_names.push_back(protonic_electrical_potential);
    this->sol_names.push_back(electronic_electrical_potential);

    this->reactant = react;

}





//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
NAME::AgglomerateBase::declare_parameters (ParameterHandler &param) const
{
    param.enter_subsection("AgglomerateBase");{
        param.declare_entry("Radius of the agglomerate [nm]", "200",
                Patterns::Double());
        param.declare_entry("Constant agglomerate parameter [Thickness | Porosity]","Porosity",
                Patterns::Selection("Thickness|Porosity"),
                "Variable used to select if thickness or porosity should be constant");
        param.declare_entry("Thickness of the agglomerate film [nm]", "15",
                Patterns::Double());
        param.declare_entry("Agglomerate porosity", "0.25",
                Patterns::Double());

    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
void
NAME::AgglomerateBase::initialize (ParameterHandler &param)
{
    param.enter_subsection("AgglomerateBase");
    {
        r_agg = param.get_double("Radius of the agglomerate [nm]");
        delta_agg = param.get_double("Thickness of the agglomerate film [nm]");
        epsilon_agg = param.get_double("Agglomerate porosity");
        fixed_agg_variable = param.get("Constant agglomerate parameter [Thickness | Porosity]");


    }
    param.leave_subsection();

    _initialize_film_porosity();

}

void
NAME::AgglomerateBase::_initialize_film_porosity(){
    if (fixed_agg_variable.compare("Thickness") == 0)
    {

        epsilon_agg = compute_epsilon_agg(this->layer);
        // n_agg depends on epsilon_agg, therefore the former needs to be computed first.
        n_agg = compute_n(this->layer);
    }
    else if (fixed_agg_variable.compare("Porosity") == 0)
    {
        n_agg = compute_n(this->layer);
        // thickenss_agg depends on n_agg, therefore it is computed last:
        delta_agg = compute_thickness_agg(this->layer);
    }
    else
    {
        FcstUtilities::log<<"Undefined execution path in AgglomerateBase::initialize "<< std::endl;
        abort();
    }

}


//---------------------------------------------------------------------------
void
NAME::AgglomerateBase::print_properties(){
    FcstUtilities::log << "=========== CL MICROSTRUCTURE =========" << std::endl;
    FcstUtilities::log <<  "Agglomerate Type: " << get_name() << std::endl;
    FcstUtilities::log <<  "Agglomerate Radius: " << std::setw(6) <<   get_radius() << " [nm]" << std::endl;
    FcstUtilities::log <<  "Agglomerate Porosity: " << std::setw(6) << epsilon_agg << std::endl;
    FcstUtilities::log <<  "Agglomerate Thin Film: " << std::setw(6) << get_film_thickness() << " [nm]" << std::endl;
    FcstUtilities::log << "=======================================" << std::endl;
}



void
NAME::AgglomerateBase::make_thread_safe(ParameterHandler &param, unsigned int thread_index){

    //Create Local copies of important objects

    catalyst = FuelCellShop::Material::CatalystBase::create_Catalyst(param, param.get("Catalyst type"));
    electrolyte = FuelCellShop::Material::PolymerElectrolyteBase::create_PolymerElectrolyte(param, param.get("Electrolyte type"));
    kinetics = FuelCellShop::Kinetics::BaseKinetics::create_Kinetics(param, param.get("Kinetics type"));
    kinetics->set_catalyst(catalyst.get());
    kinetics->set_electrolyte(electrolyte.get());

    ReactionNames rxn_name;

    /*
     *  This following code determines reaction name based on type of kinetics model
     *  which may not necessarily be correct and only will work for currently
     *  implemented kinetics models.
     *  A more robust mothod of determining the reaction name at this scope is required.
     *  Use of new kineitics models, or existing kinetics models for different reactions
     *  will result in a bug of the following comparisions.
     */
    if((typeid(*kinetics.get()) == typeid(FuelCellShop::Kinetics::TafelKinetics)) or (typeid(*kinetics.get()) == typeid(FuelCellShop::Kinetics::DoubleTrapKinetics))){
        rxn_name = ORR;
    }
    else if(typeid(*kinetics.get()) == typeid(FuelCellShop::Kinetics::DualPathKinetics)){
        rxn_name = HOR;

    }


    kinetics->set_reaction_kinetics(rxn_name);
    catalyst->set_reaction_kinetics(rxn_name);

}

//---------------------------------------------------------------------------
double
NAME::SphericalAgglomerateGeometry::compute_n(FuelCellShop::Layer::MultiScaleCL<deal_II_dimension>* layer) const
{
    double n_agg = 0.0;

    CL_Properties props = layer->get_properties();

    n_agg =props[CLPropNames::solid_fraction]/ ((4.0 / 3.0) * pi * pow(r_agg * 1e-7, 3.0)
                                * (1 - epsilon_agg)); //# of agglomerates per unit volume

    return n_agg;
}

//---------------------------------------------------------------------------
double
NAME::SphericalAgglomerateGeometry::compute_dn_depsilon_agg(FuelCellShop::Layer::MultiScaleCL<deal_II_dimension>* layer) const
{
    double dn_de = 0.0;



    dn_de = layer->get_properties()[CLPropNames::solid_fraction]
            / ((4.0 / 3.0) * pi * pow(r_agg * 1e-7, 3.0)
            * pow(1 - epsilon_agg, 2.0)); //# of agglomerates per unit volume

    return dn_de;

}

//---------------------------------------------------------------------------
double
NAME::SphericalAgglomerateGeometry::compute_thickness_agg(FuelCellShop::Layer::MultiScaleCL<deal_II_dimension>* layer)
{
    //set default and starting values
    double thickness_0 = r_agg / 10.;
    const double tol = 1e-6;
    const double eps_tol = 1e-6;
    const double iterations = 1000;

    CL_Properties props = layer->get_properties();

    double porosity =  epsilon_agg;

    double thickness_i = thickness_0;

    for (int i = 0; i < iterations; ++i)
    {
        thickness_i = thickness_0
                - ((compute_epsilon_N(thickness_0, porosity) - props[CLPropNames::ionomer_fraction])
                        / compute_depsilonN_dthickness(thickness_0));

        if (fabs((thickness_i - thickness_0) / thickness_0) < tol
                || fabs(compute_epsilon_N(thickness_i, porosity) - compute_epsilon_N(thickness_0, porosity))
                < eps_tol)
        {
            if (thickness_i < 0)
            {
                throw std::runtime_error("Agglomerate film thickness has converged to a negative value.");
            }
            else
            {
                return thickness_i;
            }
        }
        else
            thickness_0 = thickness_i;

    }

    throw std::runtime_error("Agglomerate film thickness has not converged (r_agg = "+ std::to_string(r_agg) +")" );
}

//---------------------------------------------------------------------------
double
NAME::SphericalAgglomerateGeometry::compute_epsilon_agg(FuelCellShop::Layer::MultiScaleCL<deal_II_dimension>* layer)
{
    //
    // n_agg depends on epsilon_agg and n_agg is used in epsilon_N calcs ...
    //
    CL_Properties props = layer->get_properties();

    //set default and starting values
    double epsilon_0 = 0.4;
    const double tol = 1e-6;
    const double iterations = 1000;

    double epsilon_i = epsilon_0;

    for (int i=0; i < iterations; ++i)
    {
        n_agg = compute_n(layer);
        epsilon_i = epsilon_0 - (compute_epsilon_N(delta_agg, epsilon_0) - props[CLPropNames::ionomer_fraction]) / compute_depsilonN_depsilon_agg(layer);
        epsilon_agg = epsilon_i;

        if (fabs(compute_epsilon_N(delta_agg, epsilon_i) - props[CLPropNames::ionomer_fraction]) < tol) // option: can also check for change in epsilon_N_cat function
        {
            if (epsilon_i < 0)
            {
                return 1e-6;
                FcstUtilities::log<<"Warning:: Thickness is negative -- Results to not make sense"<<std::endl;

            }  // note: simulation will still run - porosity will be very low in the CL
            else
            {
                return epsilon_i;

            }
        }
        else
            epsilon_0 = epsilon_i;
    }

    //throw std::runtime_error("Agglomerate film thickness has converged to a negative value.");
}

//---------------------------------------------------------------------------
double
NAME::SphericalAgglomerateGeometry::compute_epsilon_N(const double delta_agg,
                                                   const double porosity_agg) const
{
    double value = 0.0;

        value = n_agg * (4.0 / 3.0) * pi
        * (pow(r_agg * 1e-7, 3.0) * (porosity_agg - 1)
        + pow(r_agg * 1e-7 + delta_agg * 1e-7, 3.0));

    return value;
}

//---------------------------------------------------------------------------
double
NAME::SphericalAgglomerateGeometry::compute_depsilonN_dthickness(const double thickness_agg) const
{
    double value = 0.0;

        value = n_agg * (4.0) * pi
                * (pow(r_agg * 1e-7 + thickness_agg * 1e-7, 2.0)) * 1e-7;

    return value;

}

//---------------------------------------------------------------------------
double
NAME::SphericalAgglomerateGeometry::compute_depsilonN_depsilon_agg(FuelCellShop::Layer::MultiScaleCL<deal_II_dimension>* layer) const
{
    double value = 0.0;

        value = n_agg * (4.0 / 3.0) * pi
        * (pow(r_agg * 1e-7, 3.0))
        +
        compute_dn_depsilon_agg(layer) * (4.0 / 3.0) * pi
        * (pow(r_agg * 1e-7, 3.0) * (epsilon_agg - 1)
        + pow(r_agg * 1e-7 + delta_agg * 1e-7, 3.0));

    return value;
}

