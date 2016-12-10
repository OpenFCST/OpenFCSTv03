//---------------------------------------------------------------------------
// C++ Interface: agglomerate_ionomer_sun.h
//
// Description: Used to solve a system of equations representing a
//              spherical ionomer-filled agglomerate.
//
// Author: Marc Secanell, Phil Wardlaw
//         University of Alberta
//
// Copyright: See COPYING file that comes with this distribution
//
//---------------------------------------------------------------------------

#ifndef FUEL_CELL__IONOMER_AGGLOMERATE_SUN__H
#define FUEL_CELL__IONOMER_AGGLOMERATE_SUN__H


//------------------------------
// FUEL CELL DECLARATIONS
//-----------------------------

#include <microscale/agglomerate_base.h>
#include <reactions/tafel_kinetics.h>
#include <application_core/fcst_variables.h>
#include <reactions/dual_path_kinetics.h>
#include <materials/catalyst_base.h>

namespace FuelCellShop
{
    namespace MicroScale
    {
    /**
     * \brief Analytical solution to an ionomer-filled agglomerate problem in 1D
     *
     * This class implements the analytical ionomer-filled agglomerate first proposed by:
     * 
     * - W. Sun, B. A. Peppley, K. Karan, An improved two-dimensional agglomerate
     * cathode model to study the influence of catalyst layer structural parameters,
     * Electrochimica Acta 50 (16-17):3347-3358, (2005).
     * 
     * and later implemented in openFCST in order to provide the results obtained in:
     * 
     * - M. Secanell, K. Karan, A. Suleman and N. Djilali, "Multi-Variable Optimization of PEMFC Cathodes 
     * using an Agglomerate Model ", Electrochimica Acta, 52(22):6318-6337, 2007.
     * 
     * The data files for this article can be found in fcts/data/articles/Secanell_EA07_MultiVariable_Optimization_PEMFC_Cathodes_Agglomerate_Model/.
     * You can find three folders there one for analyzing a single case, one for an IV curve and one for running an optimization
     * case as discussed in the article.
     * 
     * <h3> Theory </h3>
     * The volumetric current density is
     * \f[
     *    \nabla \cdot \vec{i} = 4F \frac{p_{tot}x_{O_2}}{H_{O_2,N}}\left(\frac{1}{E_r k_c (1-\epsilon^{cl}_V)} + \frac{(r_{agg} + \delta_{agg})\delta_{agg}}{a_{agg}r_{agg}D_{O_2, N}} \right)^{-1}
     * \f]
     * where
     * \f[
     *    k_c = \frac{A_v i^{ref}_0}{4F(1-\epsilon^{cl}_V)c^{ref}_{O_2}} \exp(-\frac{\alpha_c F}{RT}(\phi_s - \phi_m))
     * \f]
     * and \f$ k_c \f$ is the reaction rate constant, \f$ s^{-1} \f$, \f$ i^{ref}_0 \f$ is the exchange current density, [ \f$ A \cdot cm^{-2} \f$ ], 
     * \f$ R \f$ is the gas constant, 8.315 \f$ J \cdot K^{-1} \cdot mol^{-1} \f$ ], \f$ T \f$ is the temperature, [K], \f$ \gamma \f$ is the 
     * coefficient in Tafel equation, [-], \f$ F \f$ is the Faraday constant, 96493 [\f$ C \cdot mol \f$ ], \f$ \alpha_c \f$ is the transfer coefficient, [-],
     * \f$ c^{ref}_{O_2} \f$ is the reference oxygen concentration, [\f$ mol \cdot cm^{-3}\f$ ], the term \f$ 1-\epsilon^{cat}_V \f$ is used to 
     * transform the active area in the catalyst layer, \f$ A_v \f$ to an active area inside the agglomerate and the effectiveness factor is given by
     * \f[
     *    E_r = \frac{1}{\phi_L}\left( \frac{1}{\tanh(3\phi_L)} - \frac{1}{3\phi_L} \right)
     * \f]
     * where \f$ E_r \f$ is the effectiveness factor, [-] and the Thiele's modulus for a spherical agglomerate is given by
     * \f[
     *    \phi_L = \frac{r_{agg}}{3} \sqrt{\frac{k_c}{D^{eff}}}
     * \f] where \f$D^{eff} \f$ is the effective oxygen diffusion coefficient inside the agglomerate, [ \f$ cm^2\cdot s^{-1} \f$],
     * \f$ r_{agg} \f$ is the radius of the spherical agglomerate and \f$ D^{eff} \f$ is the effective oxygen diffusion coefficient 
     * inside the agglomerate.
     * 
     * In these equations, there are several parameters that need to be obtained: \f$ H_{O_2,N} \f$, \f$ a_{agg} \f$, \f$ D_{O_2,N} \f$, \f$ A_v \f$, 
     * \f$ i^{ref}_0 \f$, \f$ c^{ref}_{O_2} \f$ and \f$ D^{eff} \f$. 
     * Parameters \f$ H_{O_2,N} \f$, \f$ D_{O_2, N} \f$, \f$ i^{ref}_0 \f$ and \f$ c^{ref}_{O_2} \f$ are input parameters to the model 
     * and are obtained from transport and electrochemical data. In the code, they are specified by the type of kinetics and catalyst used. Regarding
     * the type of kinetic model, ONLY TAFEL can be used in this case.
     * 
     * Parameters \f$ a_{agg} \f$, \f$ A_v \f$ and \f$ D^{eff} \f$ depend on the composition of the catalyst layer.
     * 
     * The parameter \f$ a_{agg} \f$ is defined as the ratio between the effective surface area usable to dissolve oxygen into the agglomerate
     * and the catalyst layer volume. This value can be related to the catalyst layer structure by
     * \f[
     *    a_{agg} = n 4 \pi (r_{agg} + \delta_{agg})^2 \epsilon^{cl}_V
     * \f]
     * where \f$ \delta_{agg} \f$ is the thin electrolyte film surrounding the agglomerate, [\f$ \mu m \f$], 
     * \f$ n \f$ is the number of agglomerates per unit volume, the term \f$ 4 \pi (r_{agg} + \delta_{agg})^2 \f$ is the 
     * surface of a single agglomerate and finally, \f$ \epsilon^{cl}_V \f$ is the catalyst layer porosity. 
     * 
     * For more information see:
     * 
     * - M. Secanell, K. Karan, A. Suleman and N. Djilali, "Multi-Variable Optimization of PEMFC Cathodes 
     * using an Agglomerate Model ", Electrochimica Acta, 52(22):6318-6337, 2007.
     * 
     * 
     * \author M. Secanell 2004-13
     */
        class IonomerAgglomerateSun : public AgglomerateBase
        {
            public:


                static const std::string concrete_name;


                /**
                 * Main function of the class used to compute the current over the whole agglomerate
                 * at the local operating conditions
                 */
                virtual SolutionMap compute_current ();

                /**
                 * Function to compute the derivative of the current density at the local operating conditions;
                 */
                virtual std::vector<double> compute_derivative_current ();

                /**
                 * Return name of class instance, i.e. concrete name.
                 */
                virtual std::string get_name(){
                        return concrete_name;
                    }

                /**
                 * Returns extra contribution to volume of layer
                 *
                 */
                virtual double aux_volume_fraction(){
                                return 0;
                            }


            protected:
                /*
                 * Virtual function for returning film thickness in nano meters.
                 *
                 */
                virtual double get_film_thickness(){
                        return delta_agg*1e7;
                    }

                /*
                 * Virtual function for returning agglomerate radius in nano meters.
                 *
                 */
                virtual double get_radius(){
                    return r_agg*1e7;
                }

                /*
                 *  Set the composition and structure of the agglomerate
                 */
                virtual void set_structure ();

                static IonomerAgglomerateSun const* PROTOTYPE;

                virtual boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase> create_replica ()
                {
                    return boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase> (new FuelCellShop::MicroScale::IonomerAgglomerateSun ());
                }

                /** Constructors */
                IonomerAgglomerateSun ();
                IonomerAgglomerateSun(std::string concrete_name);

                /*
                 * Protected virtual member function for declaring parameters, pure
                 * in MicroScaleBase, implemented here in IonomerAgglomerateSun.
                 * Calls parent AgglomerateBase
                 */
                virtual void declare_parameters (ParameterHandler &param) const
                {
                    AgglomerateBase::declare_parameters(param);
                    param.enter_subsection(concrete_name);{
                        //No parameters at the moment
                    }
                    param.leave_subsection();
                }

                /*
                 * Protected virtual member function for initializing parameters, pure
                 * in MicroScaleBase, implemented here in IonomerAgglomerateSun.
                 * Calls parent AgglomerateBase
                 */
                virtual void initialize (ParameterHandler &param)
                {
                    AgglomerateBase::initialize(param);
                    param.enter_subsection(concrete_name);{
                        //No parameters at the moment
                    }
                    param.leave_subsection();

                }


            private:

                /* Private member function to check kinetics are appropriate for analytical formulation.
                 * Throws exception in event of inappropriate kinetic conditions.
                 */
                void check_kinetics();

                /*bool to monitor if we have checked kinetic conditions */
                bool checked_kinetics;


                /** Function to compute the effectiveness of the agglomerate core */
                double compute_Er (const double k_c, const double D);

                /** Function to compute the derivative of the effectiveness of the agglomerate core */
                double compute_dEr (const double k_c,const double dk_c, const double D);
                
                /** Catalyst layer porosity */
                double epsilon_V;


        };// class IonomerAgglomerateSun

    } // namespace Layer

} // namespace FuelCellShop

#endif

