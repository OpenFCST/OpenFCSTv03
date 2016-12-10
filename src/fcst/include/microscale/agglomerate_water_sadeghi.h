// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: todo
// - Description: todo
// - Developers: Philip Wardlaw <wardlaw@ualberta.ca>
//                University of Alberta
//
// ----------------------------------------------------------------------------

#ifndef FUEL_CELL__WATER_PORE_AGGLOMERATE
#define FUEL_CELL__WATER_PORE_AGGLOMERATE

//------------------------------
// STD LIBRARY DECLARATIONS
//------------------------------


//------------------------------
// DEAL.II DECLARATIONS
//------------------------------

//------------------------------
// FUEL CELL DECLARATIONS
//-----------------------------

#include <microscale/agglomerate_base.h>

//Test Class
class WaterPoreAgglomerateTest;

namespace FuelCellShop
{
    namespace MicroScale
    {
        /**
         * \brief Class implementing water pore agglomerate model developed by Ehsan Sadeghi
         * 
         * WORK IN PROGRESS
         */
         
        class WaterConicalPoreAgglomerate : public AgglomerateBase
        {
        public:
            

            static const std::string concrete_name;


            /** Set the composition and structure of the agglomerate */
            

            /**
             * Main function of the class used to compute the current over the whole agglomerate
             * at the local operating conditions
             */
            virtual SolutionMap compute_current ( );
            
            /**
             * Function to compute the derivative of the current density at the local operating conditions;
             */
            //virtual std::vector<double> compute_derivative_current ();
            
            
            /**
             * Friend class for testing purposes.
             */

            friend class ::WaterPoreAgglomerateTest;

            virtual std::string get_name(){
                    return concrete_name;
                }

            virtual double aux_volume_fraction(){
                                            return 0;
                                        }


        protected:


            virtual double get_film_thickness(){
                return delta_agg*1e9;
            }
            virtual double get_radius(){
                return r_agg*1e9;
            }






            virtual void set_structure ();
            static WaterConicalPoreAgglomerate const* PROTOTYPE;

            /** Constructors */
            WaterConicalPoreAgglomerate ();
            WaterConicalPoreAgglomerate(std::string concrete_name);

            virtual boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase> create_replica ()
            {
                return boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase> (new FuelCellShop::MicroScale::WaterConicalPoreAgglomerate());
            }
            

            virtual void declare_parameters (ParameterHandler &param) const
            {
                AgglomerateBase::declare_parameters(param);
                param.enter_subsection(concrete_name);{
                }
                param.leave_subsection();




            }

            virtual void initialize (ParameterHandler &param)
            {
                AgglomerateBase::initialize(param);
                param.enter_subsection(concrete_name);{

                }
                param.leave_subsection();


            }


            /*
             * Domain temperature
             */
            double T;










        private:
            
            /**
             * Function for initializing data before problem is solved
             *
             * Return false if problem cannot be initialized, passes error msg back by reference
             */
            bool initialize_problem(std::string& error_msg);
            
            /**
             * Function for interpreting solution provided by another class.
             *
             * Return false if solution is invalid, passes error msg back by reference
             */
            bool interpret_solution(std::string& error_msg);


            /**
             * Functions which solve for Proton Potential and oxygen concentration Profiles across the pore.
             * Call these function in the following order: solveProtonPotentials, solveO2.
             *
             * Return false if they fail to converge, passes error msg back by reference
             */
            bool solveProtonPotentials(std::string& error_msg);
            bool solveO2(std::string& error_msg);


            /**
             * Functions which calculates current density,
             * proton potential and o2 concentration must be solved.
             *
             * Return current density in A/cm^3, passes agglomerate effectiveness back by reference
             */
            double calculate_j(double& E_r);

            enum SolutionIteration{
                Old,
                New,
                Intermediate
            };




            std::map<SolutionIteration, std::vector<std::vector<double>>> PorePotential;
            std::map<SolutionIteration,std::vector<std::vector<double>>> cO2;
            std::vector<double> cHwall;
            std::vector<double> etaWall;
            std::vector<double> jWall;



            //u_in goes inside newton-rapson function
            //double gama, u[201, 201, 201], h, k, pii, beta, kisi, CH;
            //f_u, df_u,


            //Factor F/(R*T)
            double beta;

            //Surface OverPotential
            double eta;

            //Angular steps
            double alfa, lambda;

            //Pore angle
            double teta;

            //Inner and outer radii
            double innerRadius;
            double outerRadius;

            //Number of pores per agglomerate
            double NP;

            //Surface Area inner/outer pore
            double S_i, S_o;
            //Dimensionless source term, dimensionless initial potential
            double eva, u_0;

            //Potential difference at agglomerate surface, dimensionless omega
            double omega, del_u;

            //Potential of zero charge
            double phi_pzc;

            //Dimensionless surface potential
            double u_s;

            //Iteraion Number for r and teta
            int N, M;

            //Helmholtz Capacitance
            double Helm;

            //Radius of a water molecule
            double waterRadius;

            //Relative permittivity of water
            double rel_permittivity;


        };// class IonomerAgglomerate2
        
    } // namespace Layer
    
} // namespace FuelCellShop

#endif
