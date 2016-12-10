// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: Micro Scale Initial Solution
// - Description: INSERT SHORT DEFINITION HERE
// - Developers: Philip Wardlaw, University of Alberta
// - $Id: numerical_agglomerate_base.h 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#include <microscale/agglomerate_base.h>
#include <utils/fcst_db.h>


#ifndef NUMERICAL_AGGLOMERATE_BASE_H_
#define NUMERICAL_AGGLOMERATE_BASE_H_


//Friend class for testing
class NumericalAgglomerateBaseTest;

namespace FuelCellShop
{
    namespace MicroScale
    {
        /**
         *
         * @brief Base class for numerical agglomerates
         *
         * Base class for numerical agglomerates, implements AgglomerateBase.
         * Provide useful methods and data for numerical agglomerates, such as
         * initial solution storage and retrieval via FCSTDatabase and
         * variable active area profiles.
         *
         *
         * <h3> Input parameters </h3>
         * LIST OF INPUT PARAMETERS FOR THE CLASS.
         * @code
         * ...
         * subsection MicroScaleBase
         *   subsection NumericalAgglomerateBase
         *     set Initial condition tolerance factor = 0.05 #Factor determining how often initial solutions are store/updated
         *     set Database name = main_db #Database existing in FCST/Databases. Setting an unused name will create a new database.
         *     set Agglomerate Loading Profile = 1 #A list of doubles of the form "1,2.2,3...", used to weight the Pt. loading profile
         *   end
         * end
         * ...
         * @endcode
         *
         * @author Philip Wardlaw
         * @date 2014
         */
        class NumericalAgglomerateBase : public AgglomerateBase
        {
        public:
            friend class::NumericalAgglomerateBaseTest;

            //Nothing public
        protected:


            //Set AV function Just for testing purposes
            inline void setAV(double newAV){
                AV = newAV;
            }

            /** Constructor */
            NumericalAgglomerateBase();

            /** Destructor */
            ~NumericalAgglomerateBase(){}

            /*
             * Protected member function for updating initial solution.
             * Function will determine if initial solution should be updated or not.
             *
             * <h3> Usage </h3>
             * Should be called before solving for current density.
             */
            void update_initial_solution();

            /*
             * Protected member function for saving an initial solution.
             * Function will determine if initial solution should be saved or not.
             *
             * <h3> Usage </h3>
             * Should be called after successfully solving system of equations.
             */
            void save_initial_solution();

            /*
             * Protected virtual member function for declaring parameters,
             * calls parent class (AgglomerateBase).
             */
            virtual void declare_parameters (ParameterHandler &param) const;

            /*
             * Protected virtual member function for initializing parameters,
             * calls parent class (AgglomerateBase).
             */
            virtual void initialize(ParameterHandler &param);

            /*
             * Protected member function for obtaining active area at
             * a given location in domain. Active area may vary depending on
             * loading profile
             *
             * \param location The current location in the agglomerate domain
             * from 1.0 to 0.0 (dimensionless linear coordinates)
             *
             * Units [cm^2 of pt surface per /cm^3 of solid CL micro structure]
             */
            double getAV(double location);

            /*
             * Function to be called by COLDAE guess function to get initial solution at point x
             *

             */
            bool use_initial_data(double z[], const double &x);

            /*
             * 2d Solutions
             */
            std::vector<std::vector<double>> final_results;
            std::vector<std::vector<double>> initial_solution;

            /*
             * Names of columns for solution data
             * <b>Important:</b> Should be implemented by child.
             */
            std::vector<std::string> column_names;


            /*
             * Absolute tolerance in reactant concentration
             */
            double R_tol;


            virtual void make_thread_safe(ParameterHandler &param, unsigned int thread_index);

        private:

            /*
             * Private function to set up the active area loading profiles depending on user input from user.
             *
             * Active area is conserved.
             *
             * See section 2.4.2.3 of Philip Wardlaw's thesis
             */

            void setUpLoadings();

            /*
             * Private member function for creating a snapshot if current operating conditions
             * so that we can query database.
             */
            FcstUtilities::DatabaseOC create_OC_snapshot();

            /*
             * Private to check if current initial solution is still valid given a set of operating conditions.
             */
            bool guess_no_longer_valid(const FcstUtilities::DatabaseOC& OC);

            /*
             * Private for interpolating initial solution.
             */
            void interpolate_initial_data(double z[], const double &x, const double& left_pos, const double& right_pos,
                    const std::vector<double>& left_data, const std::vector<double>& right_data);

            /*
             * Database variables
             */
            FcstUtilities::FCSTdatabase db;
            std::string db_address;
            FcstUtilities::DatabaseOC database_OC;

            /*
             * Bool determining if we should push data on next call to save_initial_solution();
             */
            bool push_next;

            /*
             * Tolerance determining how often we update/save initial solutions.
             */
            double tolerance;

            /*
             * Parameters for variable active area profile
             */
            double maxRadialDimension;
            std::vector<double> loadingWeigths;
            std::vector<double> actualLoadings;
            std::vector<double> loadingRadii;

            unsigned int thread_id;

        }; //class
    }
} //namespace

#endif
