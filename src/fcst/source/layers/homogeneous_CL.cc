//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//    - Class: homogeneous_CL.h
//    - Description: Class characterizing the macro-homogeneous catalyst layer.
//    - Developers: Marc Secanell, Peter Dobson and Madhur Bhaiya
//    - Id: $Id: homogeneous_CL.cc 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#include <layers/homogeneous_CL.h>

namespace NAME = FuelCellShop::Layer;

template <int dim>
const std::string NAME::HomogeneousCL<dim>::concrete_name ("HomogeneousCL");

template <int dim>
NAME::HomogeneousCL<dim> const* NAME::HomogeneousCL<dim>::PROTOTYPE = new NAME::HomogeneousCL<dim>();

//---------------------------------------------------------------------------
template <int dim>
NAME::HomogeneousCL<dim>::HomogeneousCL()
  : ConventionalCL<dim>()
{
    this->get_mapFactory()->insert(std::pair<std::string, FuelCellShop::Layer::CatalystLayer<dim>* >(concrete_name, this));
}

//---------------------------------------------------------------------------
template <int dim>
NAME::HomogeneousCL<dim>::HomogeneousCL(std::string name)
  : ConventionalCL<dim>(name)
{
}

//---------------------------------------------------------------------------
template <int dim>
NAME::HomogeneousCL<dim>::~HomogeneousCL()
{}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::HomogeneousCL<dim>::declare_parameters (const std::string& cl_section_name, 
                                              ParameterHandler &param) const
{
    NAME::ConventionalCL<dim>::declare_parameters(cl_section_name, param);
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::HomogeneousCL<dim>::initialize (ParameterHandler &param)
{
    NAME::ConventionalCL<dim>::initialize(param); 
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::HomogeneousCL<dim>::current_density(std::vector<double>& coef)
{
    this->kinetics->current_density(coef);
    for (unsigned int i = 0; i<coef.size(); ++i)
        coef[i] *= this->Av.at(this->local_material_id());
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::HomogeneousCL<dim>::derivative_current_density(std::map< VariableNames, std::vector<double> >& dcoef_du)
{
    this->kinetics->derivative_current(dcoef_du);
    for (std::map< VariableNames, std::vector<double> >::iterator iter=dcoef_du.begin(); iter!=dcoef_du.end(); ++iter)
        for (unsigned int q=0; q<iter->second.size(); ++q)
            iter->second[q] *= this->Av.at(this->local_material_id());
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations. 
template class NAME::HomogeneousCL<deal_II_dimension>;