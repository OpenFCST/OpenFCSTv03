//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-2014 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: PolyAgglomerate
//    - Description: Poly agglomerate class (work in progress)
//    - Developers: Philip Wardlaw, University of Alberta
//    - $Id: poly_agglomerate.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include <microscale/numerical_agglomerate_base.h>

namespace NAME = FuelCellShop::MicroScale;

NAME::NumericalAgglomerateBase::NumericalAgglomerateBase(){
    maxRadialDimension = 1.0;
    push_next = false;
    thread_id = 0;
}

//---------------------------------------------------------
void
NAME::NumericalAgglomerateBase::update_initial_solution(){
    //Check the database OC hasn't changed outside tolerance




    FcstUtilities::DatabaseOC snap_shot = create_OC_snapshot();

    if (guess_no_longer_valid(snap_shot)){
        if(db.connect(db_address,true)){

            unsigned int i = 0;

            for(; i<3; i++){
                if (db.has_data(this->get_name(), snap_shot, tolerance*double(1+i))){
                    initial_solution = db.get_data(this->get_name(),
                            snap_shot, tolerance*double(1+i),  column_names[0]); //TODO: change to variable solution


                    database_OC = snap_shot;
                    break;
                }
            }

            if(i>0) //If it took more than one attempt to get the data
                push_next = true; //Next time when asked, we will push new data

            db.disconnect();
        }

    }
    //The child class will now solve for these new operating conditions
    //And will store it's solution to be used again, hence it will be our
    //new guess
}

//---------------------------------------------------------
void
NAME::NumericalAgglomerateBase::save_initial_solution(){

    if (thread_id != 0) //Only thread 0 writes
           return;

    if(push_next){
        if(db.connect(db_address,true)){
            database_OC = create_OC_snapshot();

            if(not db.has_data(this->get_name(),database_OC, tolerance)){
                db.commit_data(this->get_name(), database_OC, column_names, final_results);
            }

            db.disconnect();
        }

        push_next = false;
    }
}

//---------------------------------------------------------
bool
NAME::NumericalAgglomerateBase::guess_no_longer_valid(const FcstUtilities::DatabaseOC& OC){
    return not database_OC.compare(OC, tolerance);
}

//---------------------------------------------------------
FcstUtilities::DatabaseOC
NAME::NumericalAgglomerateBase::create_OC_snapshot(){

    FcstUtilities::DatabaseOC OC;
    OC.add_param("x_O2", this->solutions[this->reactant][this->sol_index]);
    OC.add_param("phi_m",this->solutions[protonic_electrical_potential][this->sol_index]);
    OC.add_param("phi_s", this->solutions[electronic_electrical_potential][this->sol_index]);
    OC.add_param("radius", get_radius());
    OC.add_param("film_thickness", get_film_thickness());
    OC.add_param("porosity", epsilon_agg);
    OC.add_param("kinetics", typeid(*this->kinetics).name());
    OC.add_param("kin_param",this->layer->get_resource<FuelCellShop::Material::CatalystBase>()->get_kinetic_parameter_method());

    return OC;
}

//---------------------------------------------------------
void
NAME::NumericalAgglomerateBase::interpolate_initial_data(double z[], const double &x, const double& left_pos, const double& right_pos, const std::vector<double>& left_data, const std::vector<double>& right_data)
{
    //start at i=1 to ignore the first element which is x
    for(unsigned int i = 1; i<left_data.size(); i++)
        z[i-1] = left_data.at(i) + ((right_data.at(i) - left_data.at(i))*(x - left_pos))/(right_pos - left_pos);

}

//---------------------------------------------------------
bool
NAME::NumericalAgglomerateBase::use_initial_data(double z[], const double &x){
    bool answer = false;

    if((initial_solution.size() != 0)){


        //Start at centre, centre values are scaled from the first value
        double left = 0;
        std::vector<double> left_data;
        for(int j =0; j<initial_solution.at(0).size(); j++)
            left_data.push_back(0.8*initial_solution.at(0)[j]);

        double right;
        for(unsigned int i = 0; i <initial_solution.size(); i++){

            right = initial_solution.at(i)[0];
            if( (left<= x) && (x <=  right)){
                interpolate_initial_data(z,x,left,right, left_data, initial_solution.at(i));
                answer = true;
                break;
            }

            left = right;
            left_data = initial_solution.at(i);
        }

        //If we went past the right most point and haven't yet interpolated
        if (!answer and x > right){
            for(int k = 0; k < initial_solution.at(0).size() -1; k++)
                z[k] = initial_solution.at(initial_solution.size()-1)[k+1];

            answer = true;
        }

    }

    return answer;

}

//---------------------------------------------------------
void
NAME::NumericalAgglomerateBase::setUpLoadings(){

    if (loadingWeigths.size() != 0){

        for(unsigned int i = 0; i <= loadingWeigths.size(); i++)
            loadingRadii.push_back( (interface/loadingWeigths.size()) * i );

        double weightingSum =0;

        for(unsigned int i = 0; i < loadingWeigths.size(); i++)
            weightingSum += loadingWeigths[i]*( std::pow(loadingRadii[i+1], 3.0) - std::pow(loadingRadii[i],3.0));

        for(unsigned int i = 0; i < loadingWeigths.size(); i++)
            actualLoadings.push_back(double(AV/weightingSum)*loadingWeigths[i]*std::pow(interface, 3.0));

    }
}

//---------------------------------------------------------
double
NAME::NumericalAgglomerateBase::getAV(double location){


    if(actualLoadings.size() != loadingWeigths.size())
        setUpLoadings();

    if (loadingWeigths.size() == 0)
        return AV;


    for(unsigned int i = 0; i < actualLoadings.size(); i++)
        if(location > loadingRadii[i] && location <= loadingRadii[i+1]){
            return actualLoadings[i];
        }

    if(location < interface)
        return actualLoadings.back();

    else if(location <= 1.0)
        return 0.0;

    //If the code has not returned by this stage then execution is exceptional
    throw std::runtime_error("NumericalAgglomerateBase::getAV invalid location (" + std::to_string(location) + ")");


}

//---------------------------------------------------------
void
NAME::NumericalAgglomerateBase::declare_parameters (ParameterHandler &param) const
{

    AgglomerateBase::declare_parameters(param);

    param.enter_subsection("NumericalAgglomerateBase");{
        param.declare_entry("Initial condition tolerance factor",
                "0.05", Patterns::Double(0.0005, 0.5));
        param.declare_entry("Database name", "main_db", Patterns::Anything(),
                "The name of the database that reside in FCST root");
        param.declare_entry("Agglomerate Loading Profile", "1",
                Patterns::List(Patterns::Double()));
    }
    param.leave_subsection();

}

//---------------------------------------------------------
void
NAME::NumericalAgglomerateBase::initialize(ParameterHandler &param){


    AgglomerateBase::initialize(param);

    param.enter_subsection("NumericalAgglomerateBase");{
        tolerance = param.get_double("Initial condition tolerance factor");
        db_address = FcstUtilities::find_fcst_root() + "databases/" + param.get("Database name");
        loadingWeigths = FcstUtilities::string_to_number<double>( Utilities::split_string_list( param.get("Agglomerate Loading Profile")));
    }
    param.leave_subsection();
}


void
NAME::NumericalAgglomerateBase::make_thread_safe(ParameterHandler &param, unsigned int thread_index){
    thread_id = thread_index;
    AgglomerateBase::make_thread_safe(param, thread_index);

    if (thread_id == 0){
        //Avoid multiple creation, create now if necessary
        if(db.connect(db_address,true));
            db.disconnect();
    }

}
