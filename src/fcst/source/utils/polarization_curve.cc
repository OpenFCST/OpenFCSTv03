//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: polarization_curve.cc
//    - Description: Child of ApplicationWrapper used to implement IV curve analysis
//    - Developers: M. Secanell
//
//---------------------------------------------------------------------------

#include "utils/polarization_curve.h"

//---------------------------------------------------------------------------
template <int dim>
FuelCell::PolarizationCurve<dim>::PolarizationCurve()
{}

//---------------------------------------------------------------------------
template <int dim>
FuelCell::PolarizationCurve<dim>::~PolarizationCurve()
{
  this->myfile.close();
}

//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::PolarizationCurve<dim>::declare_parameters(ParameterHandler& param) const
{
    // Initialize to parameters needed for evaluate (note: this is provisional)
    param.enter_subsection("Simulator");
    {
        param.enter_subsection("Polarization Curve");
        {
            param.declare_entry("Polarization curve file output",
                                "polarization_curve.dat",
                                Patterns::Anything(),
                                "File where the polarization curve results should be stored");
            param.declare_entry ("Initial voltage [V]",
                                 "1",
                                 Patterns::Double(),
                                 "Voltage at which the first point in the polarization curve will be evaluated");
            param.declare_entry ("Final voltage [V]",
                                 "0.1",
                                 Patterns::Double(),
                                 "Voltage at which the polarization curve will be terminated");
            param.declare_entry ("Increment [V]",
                                 "0.1",
                                 Patterns::Double(),
                                 "Spacing between points in the polarization curve");
            param.declare_entry ("Adaptive Increment",
                                 "true",
                                 Patterns::Bool(),
                                "Set to true if you would like to reduce the voltage increment adaptively"
                                "if convergence could not be achieved with the larger value");
            param.declare_entry ("Min. Increment [V]",
                                 "0.025",
                                 Patterns::Double(),
                                 "If Adaptive Increment? is set to true, this value controls the "
                                 "minimum change in cell voltage before the polarization curve gives up"
                                 "and the voltage is updated again. Note that this value has to be more "
                                 "than 0.01 V as a value of zero would lead to an infinite loop.");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::PolarizationCurve<dim>::initialize(ParameterHandler& param)
{
    param.enter_subsection("Simulator");
    {
        param.enter_subsection("Polarization Curve");
        {
            this->parameter_filename = param.get("Polarization curve file output");
            this->p_init = param.get_double("Initial voltage [V]");
            this->p_end = param.get_double("Final voltage [V]");
            this->dp =  param.get_double("Increment [V]");
            this->adaptive = param.get_bool("Adaptive Increment");
            this->min_dp = param.get_double("Min. Increment [V]");
            if (this->min_dp < 0.01)
                this->min_dp = 0.01;                
        }
        param.leave_subsection();
    }
    param.leave_subsection();

    this->n_dpPts = floor( (this->p_init - this->p_end)/this->dp ) + 1;
    
    // Initialize convergence flag
    this->convergence = false;
    
    // Initialize coarse solution w/ empty vector:
    this->coarse_solution = FuelCell::ApplicationCore::FEVector();
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// PRIVATE:
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::PolarizationCurve<dim>::set_parameters(ParameterHandler& param,
                                                 const shared_ptr<FuelCell::ApplicationCore::AdaptiveRefinement<dim> >& solver,
                                                 const int iteration,
                                                 const std::vector<std::string> parameter_name,
                                                 const std::vector<double> param_value)
{
    //-- Set design variables to the values given by DAKOTA:
    // Make sure that boundary conditions are modified for cell voltage:
    FcstUtilities::modify_parameter_file("Fuel cell data>>Operating conditions>>Adjust initial solution and boundary conditions", true, param);
    // Modify cell voltage:
    FcstUtilities::modify_parameter_file("Fuel cell data>>Operating conditions>>Voltage cell [V]", param_value[0], param); 
    // Make sure a solution is not read from file after the first iteration:
    if (iteration != 0)
        FcstUtilities::modify_parameter_file("Initial Solution>>Read in initial solution from file", false, param);

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::PolarizationCurve<dim>::print_parameters() const
{
    FcstUtilities::log<<"== Polarization curve parameters: =="<<std::endl;
    FcstUtilities::log<<"Initial voltage [V] : "<<this->p_init<<std::endl;
    FcstUtilities::log<<"Final voltage [V] : "<<this->p_end<<std::endl;
    FcstUtilities::log<<"Increment [V] : "<<this->dp<<std::endl;
    FcstUtilities::log<<"Adaptive Increment : "<<this->adaptive<<std::endl;
    FcstUtilities::log<<"Min. Increment [V] : "<<this->min_dp<<std::endl;
    FcstUtilities::log<<"==  =="<<std::endl;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

template <int dim>
void
FuelCell::PolarizationCurve<dim>::print_parameteric_study_header() 
{
    
    this->myfile.open (this->parameter_filename);
    this->myfile<< "# ===================================================="<<std::endl;
    this->myfile<< "# OpenFCST: Fuel cell simulation toolbox "<<std::endl;
    this->myfile<< "# ===================================================="<<std::endl;
    this->myfile<< "# Polarization curve data :"<<this->parameter_filename<<std::endl;
    this->myfile<< "# "<<std::endl;
    
    std::string header;
    
    
    header.append("Cell voltage [V]");
    for (unsigned int i=0; i<this->n_resp; i++)
    {
        if (this->name_responses[i] == "current")
        {
            header.append("\t");
            header.append("Cathode current [A/cm2]");
        }
        else
        {
            header.append("\t");
            header.append(this->name_responses[i]);
        }
    }
    this->myfile<< header<<std::endl;
}
//---------------------------------------------------------------------------
// Explicit instantations
template class FuelCell::PolarizationCurve<deal_II_dimension>;