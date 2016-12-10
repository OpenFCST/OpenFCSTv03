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

#include <microscale/poly_agglomerate.h>

namespace NAME = FuelCellShop::MicroScale;


//Declare static variables
bool NAME::PolyAgglomerate::inf_loop_preventor = true;
int NAME::PolyAgglomerate::num_sub_micro = 10;

const std::string NAME::PolyAgglomerate::concrete_name("PolyAgglomerate");
NAME::PolyAgglomerate const* NAME::PolyAgglomerate::PROTOTYPE = new NAME::PolyAgglomerate(concrete_name);



//---------------------------------------------------------------------------------//
NAME::PolyAgglomerate::PolyAgglomerate(std::string concrete_name) {
    this->get_mapFactory()->insert(std::pair<std::string, NAME::PolyAgglomerate*>(
            concrete_name, this));
}

//---------------------------------------------------------------------------------//
NAME::PolyAgglomerate::~PolyAgglomerate() {
    // Auto-generated destructor stub
}

//---------------------------------------------------------------------------------//
void
NAME::PolyAgglomerate::set_solution(const std::map<VariableNames,SolutionVariable>& sols,const VariableNames& name, const int& index){

    for(unsigned int i = 0; i < micro.size(); i++)
        micro.at(i)->set_solution(sols, name, index);
}

//---------------------------------------------------------------------------------//
FuelCellShop::SolutionMap
NAME::PolyAgglomerate::compute_current (){
    //Add contributions
    double curr = 0.0;
    double eff = 0.0;
    double C_OH = 0.0;
    double C_O = 0.0;



    for(unsigned int i = 0; i < micro.size(); i++){
        double vol = micro.volAt(i);
        SolutionMap sol = micro.at(i)->compute_current();

        curr += vol*sol.at(current_density)[0];
        eff += vol*sol.at(CL_effectiveness)[0];

        if(sol.has(OH_coverage))
            C_OH += vol*sol.at(OH_coverage)[0];
        if(sol.has(O_coverage))
            C_O += vol*sol.at(O_coverage)[0];



    }


    SolutionMap sols;
    sols.push_back(SolutionVariable(curr, 1, current_density));
    sols.push_back(SolutionVariable(eff,1, CL_effectiveness));
    sols.push_back(SolutionVariable(C_OH,1, OH_coverage));
    sols.push_back(SolutionVariable(C_O,1, O_coverage));

    return sols;

}

//---------------------------------------------------------------------------------//
double
NAME::PolyAgglomerate::aux_volume_fraction(){
    //Add contributions
    double aux_v = 0.0;


        for(unsigned int i = 0; i < micro.size(); i++){
            aux_v += micro.volAt(i)*micro.at(i)->aux_volume_fraction();
        }



    return aux_v;
}




//---------------------------------------------------------------------------------//
void
NAME::PolyAgglomerate::declare_parameters (ParameterHandler &param) const{

    //An infinite loops could occur in the parameter tree since the PolyAgglomerate
    //declares the MicroScale class parameters... which intern declare the
    //PolyAgglomerate parameters. Example:
    // subsection MicroScale
    //    subsection PolyAgglomerate
    //       subsection MicroScale
    //          subsection PolyAgglomerate
    //             ... and so on
    //
    //Therefore we employ this infinite loop prevention strategy
    //that will prevent nesting of PolyAgglomerate inside PolyAgglomerate


    if(inf_loop_preventor == true)
    {
        param.enter_subsection(concrete_name);

        //Support a combination of 5 micro structures

        //Prevent other instances of micro scale from infinitely nesting
        inf_loop_preventor = false;
        for(unsigned int i= 0; i < num_sub_micro; i++){
            param.enter_subsection("MicroStructure" + std::to_string(i));
            param.declare_entry("Volume fraction", "0.0", Patterns::Double(), "Volume fraction of specific micro structure type, "
                    "will be normalized once all micro structures types have been specified");

            FuelCellShop::MicroScale::MicroScaleBase::declare_MicroScale_parameters(param);

            param.leave_subsection();

        }

        //Unlock further declaration for additional instances
        inf_loop_preventor = true;

        param.leave_subsection();


    }

}


void

NAME::PolyAgglomerate::initialize (ParameterHandler &param){

    param.enter_subsection(concrete_name);
    bool volWarning = true;
    //Create micro scale objects
    for(unsigned int i= 0; i < num_sub_micro; i++){

        param.enter_subsection("MicroStructure" + std::to_string(i));

        double vol = param.get_double("Volume fraction");

        if(vol > 0.0){
            //The micro object will normalize the volume fractions
            micro.push_back(vol,
                    FuelCellShop::MicroScale::MicroScaleBase::create_MicroStructure(param,this->layer));
        }
        else if (volWarning){
            FcstUtilities::log << "Some micro structure in PolyAgglomerate " \
                    "have 0.0 or negative volume fractions, therefore are ignored." << std::endl;
            volWarning = false; // Don't warn again
        }

        param.leave_subsection();

    }

    if(micro.size() == 0)
        throw std::runtime_error("No micro structures of PolyAgglomerate were specified (check individual volume fractions)");
    else
        micro.normalize_vols();

    param.leave_subsection();


}

//---------------------------------------------------------------------------------//
void
NAME::PolyAgglomerate::set_structure (){
    //Do nothing
}









void
NAME::PolyAgglomerate::make_thread_safe(ParameterHandler &param, unsigned int thread_index){

    for(unsigned int i = 0; i < micro.size(); i++)
        micro.at(i)->make_thread_safe(param, thread_index);


}
