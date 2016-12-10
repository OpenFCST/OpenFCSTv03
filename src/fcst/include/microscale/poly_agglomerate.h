// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: PolyAgglomerate
// - Description: Poly agglomerate class (work in progress)
// - Developers: Philip Wardlaw, University of Alberta
// - $Id: poly_agglomerate.h 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#ifndef POLYAGGLOMERATE_H_
#define POLYAGGLOMERATE_H_

#include <microscale/micro_scale_base.h>

namespace FuelCellShop
{
namespace MicroScale
{


class MicroSet{
    /*
     * A simple container for storing a set of micro structures with corresponding volume fractions.
     *
     * Push back all the items you wish to and then normalize volume fractions using normalize_vols()
     *
     * Iterate the micro set by indexs  0 through (size() -1) using the at() and volAt() accessors
     */

public:


    MicroSet(){
        size_ =0;
    }

    ~MicroSet(){};

    inline void  push_back(const double & volumeFraction, boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase> micro){


        //Add the pointer
        Assert(micro, ExcMessage("boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase> with null pointer out of range for MicroMap."));
        micros.push_back(micro);
        vols.push_back(volumeFraction);
        size_++;


    }

    inline void normalize_vols(){
        double sum = 0.0;

        //Perform normalization to 1.0
        for(double v: vols)
            sum +=v;

        for(double &v: vols)
            v/=sum;
    }

    inline unsigned int size() const{
        return size_;
    }

    inline  boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase>& at(unsigned int idx){
        Assert(idx < size_, ExcMessage("Index out of range for MicroMap."));
        return micros.at(idx);
    }

    inline  double& volAt(unsigned int idx){
        Assert(idx < size_, ExcMessage("Index out of range for MicroMap."));
        return vols.at(idx);
    }


private:
    unsigned int size_;
    std::vector<boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase>> micros;
    std::vector<double> vols;


};


/**
 *
 * @brief An agglomerate model which considers a combination of agglomerates.
 *
 *
 * <h3> Theory </h3>
 * A combination of different micro structural formations may  exist within
 * the micro structure of PEFC catalyst layer. This model allows us to
 * assess a combination of different micro structural models
 *
 * <h3> Input parameters </h3>
 * LIST OF INPUT PARAMETERS FOR THE CLASS.
 * @code
 * ...
 *   subsection MicroScale
 *     set Microscale type = PolyAgglomerate
 *     #Further subsections for  sub micro scale models
 *     subsection MicroStructure0 #from 0 to 10
 *
 *       set Volume fraction = 0.1 #The volume fraction of the micro structure we are about to describe
 *
 *         subsection MicroScale
 *            #Describe the micro structure
 *         end
 *     end
 *
 *      #Describe more  sub micro scale models
 *
 *   end
 * ...
 * @endcode
 *
 * @warning The maximum number of sub micro structures is hard coded to 10 - see num_sub_micro
 *
 * @note Default volume fraction is 0.0, models with volume fractions of 0.0 will be ignored
 * @note Combination of volume fractions are normalized to 1.0
 *
 * @author Philip Wardlaw
 * @date 2014
 */
class PolyAgglomerate: public MicroScaleBase {
public:




    static const std::string concrete_name;
    PolyAgglomerate(std::string);
    PolyAgglomerate(){};
    virtual ~PolyAgglomerate();



    ///@name Initalization
    //@{
    /**
     * Function for setting the solution map(reactant concentration, phi_s, phi_m, etc.).
     * First argument provide the map of SolutionVariable. The second argument provide
     * the name of the primary reactant. The final argument is the index of the solution
     * map that the micro scale object should solve for.
     *
     * This function should be called immediatly before <b>compute_current</b> or
     * <b>compute_derivative_current</b>.
     *
     */
    virtual void set_solution(const std::map<VariableNames,SolutionVariable>& sols,const VariableNames& name, const int& index);
    //@}

    /**
     * Function used to compute the current density produced by the micro structure.
     * Returns current density in [A/cm^3 of solid/liquid phase micro structure]
     *
     *
     * Solves for solution variables set by the last call to <b>set_solution</b>.
     */
    virtual SolutionMap compute_current ( );

    /**
     * Function to compute the derivative of the current density at the local operating conditions.
     * Returns current density derivatives
     *

     * Solves for solution variables set by the last call to <b>set_solution</b>.
     *
     * <h3> Usage details</h3>
     * Call <b>has_derivatives</b> to check if it is OK to call this function.
     *
     */
    virtual std::vector<double> compute_derivative_current ()
    {
        Assert(false, ExcPureFunctionCalled());
        return std::vector<double>(3, 0.0);
    }

    /**
     * Returns true if the class instance can calculate
     * current density derivatives.
     *
     * <h3> Usage details</h3>
     * Call this function if unsure whether or not to call <b>compute_derivative_current</b>.
     *
     */
    virtual bool has_derivatives(){
        //Compute numerically
        return false;
    }

    /**
     * Return name of class instance, i.e. concrete name.
     */
    virtual std::string get_name(){
        return concrete_name;
    }

    /**
     * Print out key micro-structural dimensions, defined by child.
     */
    virtual void print_properties(){

        FcstUtilities::log << "=--------- PolyAgglomerate -----------="  << std::endl;
        FcstUtilities::log << "= Number of sub agglomerates: " << micro.size() << std::endl;

        for(unsigned int i = 0; i < micro.size(); i++){
            FcstUtilities::log << "Sub agglomerate #" << i << std::endl;
            FcstUtilities::log << "Associated volume fraction:"  << micro.volAt(i) << std::endl;
            micro.at(i)->print_properties();
        }

    }


    /**
     * MicroScale object may have extra contribution to volume of layer, e.g. water.
     *
     */
    virtual double aux_volume_fraction();


    virtual void make_thread_safe(ParameterHandler &param, unsigned int thread_index);

protected:



    //name Instance Delivery (Prototype)
    static PolyAgglomerate const* PROTOTYPE;

    /*
     * Protected pure virtual member function for declaring parameters, must be implemented
     * by all children. Parameter structure is hierarchical, therefore children
     * should call their parent's declare parameter function from their own declare
     * parameter function.
     */
    virtual void declare_parameters (ParameterHandler &param) const;

    /*
     * Protected pure virtual member function for initializing parameters, must be implemented
     * by all children. Parameter structure is hierarchical, therefore children
     * should call their parent's initialize function from their own initialize function.
     */
    virtual void initialize (ParameterHandler &param);
    /*
     * Protected member function for setting up MicroScale structure after object
     * initialization.
     *
     * Responsible for pulling structural properties from catalyst layer.
     *
     * Responsible for setting up points to materials and kinetics.
     *
     */
    virtual void set_structure ();

    ///@name Instance Delivery (Private functions)
    //@{
    /**
     * This member function is used to create an object of MicroScaleBase
     *
     * \warning This class MUST be redeclared in every child.
     */
    virtual boost::shared_ptr<MicroScaleBase> create_replica ()
    {
        return boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase> (new FuelCellShop::MicroScale::PolyAgglomerate ());
    }

    //@}

private:
    //The micro scale objects are stored within the following convenient contianer
    MicroSet micro;

    //loop preventor to control how parameter subsections are declared
    static bool inf_loop_preventor;

    //Number of sub micro scale objects allowed, hard coded in header
    static int num_sub_micro;

};



}
} /* namespace FuelCellShop */

#endif /* POLYAGGLOMERATE_H_ */
