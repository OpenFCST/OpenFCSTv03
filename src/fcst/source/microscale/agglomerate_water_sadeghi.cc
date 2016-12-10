// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: agglomerate_water_sadeghi.cc
// - Description: in development
// - Developers: Philip Wardlaw
// - $Id: agglomerate_water_sadeghi.cc 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#include <microscale/agglomerate_water_sadeghi.h>

namespace NAME = FuelCellShop::MicroScale;


const std::string NAME::WaterConicalPoreAgglomerate::concrete_name("WaterConicalPoreAgglomerate");
NAME::WaterConicalPoreAgglomerate const* NAME::WaterConicalPoreAgglomerate::PROTOTYPE = new NAME::WaterConicalPoreAgglomerate(concrete_name);


NAME::WaterConicalPoreAgglomerate::WaterConicalPoreAgglomerate(std::string concrete_name){
    //Model not complete, therefore do not perform registration
    //FcstUtilities::log <<" Register " + concrete_name +  " to FactoryMap"<<std::endl;
    //this->get_mapFactory()->insert(
    //        std::pair<std::string, NAME::WaterConicalPoreAgglomerate*>(
    //                concrete_name, this));
}


//---------------------------------------------------------------------------
NAME::WaterConicalPoreAgglomerate::WaterConicalPoreAgglomerate()
{
    this->has_derivatives_= false;

    phi_pzc = 0.55;
    teta=6*pi/180;


    N = 45; M = 55;
    Helm = 0.2;
    waterRadius = 0.28e-9; //TODO:: move to water object



}


//---------------------------------------------------------------------------
void
NAME::WaterConicalPoreAgglomerate::set_structure()
{
    //catalyst_support = CS;
    //catalyst = Cat;




    //CL_Properties props = this->layer->get_properties();

    //electrolyte = this->layer->get_resource<FuelCellShop::Material::PolymerElectrolyteBase>();

    //AV = props[CLPropNames::active_area_scaled];
    r_agg = 50;
    r_agg *=1e-9;// convert to m
    epsilon_agg =0.3;



    //delta_agg = delta*1e-7; // convert to cm

    innerRadius = 1e-9/sin(teta);
    outerRadius = r_agg - 2e-9;

    NP = 8*((epsilon_agg* pow(r_agg,3) - ( pow(r_agg,3)  - pow(outerRadius,3) ))/((2*pow(sin(teta),2)*cos(teta) +
            (1 - cos(teta))*(3*pow(sin(teta),2) + (1 - pow(cos(teta),2))))*(pow(outerRadius,3) - pow(innerRadius,3))));


    //Make inner and outer non dimensional
    outerRadius = outerRadius - waterRadius/sin(teta);
    innerRadius = innerRadius - waterRadius/sin(teta);

    S_i=pi*(pow(outerRadius,2)-pow(innerRadius,2))*sin(teta);
    S_o= 2*pi*pow(outerRadius+waterRadius/sin(teta)+waterRadius,2)*(2.0/NP-1+cos(teta));

    //interface = r_agg/(r_agg+delta_agg);



}

//---------------------------------------------------------------------------
bool
NAME::WaterConicalPoreAgglomerate::interpret_solution(std::string& error_msg){

    if(this->solutions[this->reactant][this->sol_index]< 0.0){
        const std::type_info& info = typeid(*this);
        error_msg = "Negative reactant concentration supplied to " + std::string(info.name());
        return false; // Quit on non-real values
    }

    //Get the electrolyte effective properties for the current temperature
    T = this->solutions[temperature_of_REV][this->sol_index];

    double H_O2 = 5.146e5*exp(-500/T); //TODO: move to water, fix units
    //c_R =1000.0/H_O2; //(this->solutions[this->reactant][this->sol_index]*P)/H_R_N; //TODO: Rework
    phi_S = this->solutions[electronic_electrical_potential][this->sol_index];
    phi_M = this->solutions[protonic_electrical_potential][this->sol_index];

    if(this->solutions.find(proton_concentration) != this->solutions.end()){
        c_H = this->solutions[proton_concentration][this->sol_index];
    }
    else{
        c_H = electrolyte->get_density() / electrolyte->get_EW();
    }



    //Perform unit conversions, we use mol/cm^3, Ehsan uses mol/m^3
    c_R =1000.0/H_O2;//c_R /= pow(0.01,3); //TODO: Rework
    c_H /= pow(0.01,3);


    return true;
}

//---------------------------------------------------------------------------
bool
NAME::WaterConicalPoreAgglomerate::initialize_problem(std::string& error_msg){


    if(interpret_solution(error_msg) == false)
        return false;



    //Initialize variables that depend on the solutions
    beta = F/(R*T);
    u_s = phi_M*beta;
    omega = phi_S - phi_pzc + phi_M;
    del_u = omega*beta;

    rel_permittivity = 60;//water.get_Relative_Permittivity(); TODO: reimplement water into this class as object
    D_R_N = 5.3e-9;//TODO replace with water->get_DO2()/pow(0.01,2);

    eva = c_H*pow(F,2)*pow(outerRadius - innerRadius,2)/(rel_permittivity*100*permittivity_0*R*T); //Dimensional factor



    PorePotential.clear(); cO2.clear();

    //Create empty vectors
    for(unsigned int i =0; i<M+1; i++){
        PorePotential[Old].push_back(std::vector<double>(N+1, 0.0));
        PorePotential[New].push_back(std::vector<double>(N+1, 0.0));
        PorePotential[Intermediate].push_back(std::vector<double>(N+1, 0.0));

        cO2[Old].push_back(std::vector<double>(N+1, 0.0));
        cO2[New].push_back(std::vector<double>(N+1, 0.0));
    }

    //Create initial solutions

    for(unsigned int j = 0; j< N+1; j++){ //rows

        for(unsigned int i = 0; i< M+1; i++){ //columns

            if(i == M){
                PorePotential[Old].at(i).at(j) = u_s;
                cO2[Old].at(i).at(j) =1; //BC [39]
            }
            else{
                PorePotential[Old].at(i).at(j) = 0.01*beta;
                cO2[Old].at(i).at(j) =0.00001;
            }

        }
    }

    cO2[New] = cO2[Old];
    PorePotential[New] = PorePotential[Old];


    cHwall  = std::vector<double> (M+1, 0.0);
    etaWall = std::vector<double> (M+1, 0.0);
    jWall   = std::vector<double> (M+1, 0.0);

    return true;
}


//---------------------------------------------------------------------------
FuelCellShop::SolutionMap
NAME::WaterConicalPoreAgglomerate::compute_current ()
{
    double answer = 0.0;


    //True pressure is now available from the layer.
    //Was 0 at set structure time.
    //P= this->layer->get_properties()[CLPropNames::pressure];



    bool error_state = true;
    std::string error_msg = "";

    //Solve
    if (initialize_problem(error_msg))
        if (solveProtonPotentials(error_msg))
            if(solveO2(error_msg))
                error_state = false;


    if(error_state){
        error_msg = std::string(typeid(*this).name()) + ":" + error_msg;
        throw std::runtime_error(error_msg.c_str());
    }

    //Solution obtained, post process
    double E_r;
    answer = calculate_j(E_r);
    SolutionMap sols;
    sols.push_back(SolutionVariable(answer, 1, current_density));
    sols.push_back(SolutionVariable(E_r, 1, CL_effectiveness));


    return sols;

}

//---------------------------------------------------------------------------

double
NAME::WaterConicalPoreAgglomerate::calculate_j(double& E_r){

    double jPore;
    //double eff_new = 0.0;
    double h = 1.0/M;
    double alpha_c = 0.5;
    double j_0_agg = 3.836906e-2;//TODO: remove

    for(unsigned int i =0; i < M; i++){
        alfa = (i)*h;
        double r = (outerRadius - innerRadius)*alfa +innerRadius;
        double r2 = r + h*(outerRadius- innerRadius);


        //jWall.at(i) = (j_0_agg*(pow(cHwall.at(i), gamm_h))*
        //        +  exp(-alpha_c*F*etaWall.at(i)/(R*T)));
        jPore+=  (cO2[New].at(i).at(N)*(pow(exp(-PorePotential[New].at(i).at(N) + u_s),1.0))*exp(-alpha_c*F*etaWall.at(i)/(R*T))*r +
                (cO2[New].at(i).at(N)*pow(exp(-PorePotential[New].at(i + 1).at(N) + u_s),1.0))*exp(-alpha_c*F*etaWall.at(i + 1)/(R*T))*r2)*h*j_0_agg;


    }




        jPore = S_i*jPore/(innerRadius + outerRadius);
        double jSurface = S_o*j_0_agg*exp(-alpha_c*F*(-eta)/(R*T));
        double jUniform = (j_0_agg*S_i + j_0_agg*S_o)*exp(-alpha_c*F*(eta)/(R*T)); //eta_0???
        E_r = (jPore + jSurface)/jUniform;




    return jPore + jSurface;

}

//---------------------------------------------------------------------------

bool
NAME::WaterConicalPoreAgglomerate::solveProtonPotentials(std::string& error_msg){


    double h = 1.0/M;
    double k = 1.0/N;
    double kisi = h/k;
    double r,gama,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,f_u,df_u;
    c2 = 1.0;
    unsigned int steps = 0; unsigned int max_steps = 50000;

    bool converged = false;

    while(not converged and steps < max_steps){
        steps++;

        //Loop to solve electrolyte potential
        for(unsigned int i = 1; i< M; i++) //Columns
        {

            for(unsigned int j = 1; j< N; j++) //Rows
            {
                alfa = (i)*h;
                lambda = (j )*k;
                r = (outerRadius - innerRadius)*alfa + innerRadius;
                c1 = pow((outerRadius - innerRadius)/(teta*r),2);
                c3 = pow(outerRadius - innerRadius,2)/(teta*tan(lambda*teta)*pow(r,2));
                c4 = 2*(outerRadius - innerRadius)/r;
                c7 = 2*c1*pow(kisi,2) + 2*c2;
                c8 = c1*pow(kisi,2) + c3*h*kisi/2;
                c9 = c1*pow(kisi,2) - c3*h*kisi/2;
                c10 = c2 + c4*h/2;
                c11 = c2 - c4*h/2;

                PorePotential[New].at(i).at(j) = ((pow(h,2)*eva*exp(-PorePotential[Old].at(i).at(j) + u_s) + c8*PorePotential[Old].at(i).at(j+1) + c9*PorePotential[New].at(i).at(j-1)
                        + c10*PorePotential[Old].at(i+1).at(j) + c11*PorePotential[New].at(i-1).at(j))/c7);

                PorePotential[Intermediate].at(i).at(j) = 5.0;



                while(std::abs((PorePotential[New].at(i).at(j) - PorePotential[Intermediate].at(i).at(j)))>0.000001*std::abs( PorePotential[New].at(i).at(j))){

                    PorePotential[Intermediate].at(i).at(j) = PorePotential[New].at(i).at(j);
                    f_u = pow(h,2)*eva*exp(-PorePotential[Intermediate].at(i).at(j) + u_s) + c8*PorePotential[Old].at(i).at(j+1) + c9*PorePotential[New].at(i).at(j-1) +
                            c10*PorePotential[Old].at(i+1).at(j) + c11*PorePotential[New].at(i-1).at(j) - c7*PorePotential[Intermediate].at(i).at(j);
                    df_u = -pow(h,2)*eva*exp(-PorePotential[Intermediate].at(i).at(j) + u_s) - c7;
                    PorePotential[New].at(i).at(j) = PorePotential[Intermediate].at(i).at(j) - f_u/df_u;
                }

            }
        }

        //Set the boundary solution
        for(unsigned int j = 1; j < N +1 ; j++)//# j<N
        {
            int i =0;
            PorePotential[New].at(i).at(j) =  (18*PorePotential[New].at(i+1).at(j) - 9*PorePotential[New].at(i+2).at(j) + 2*PorePotential[New].at(i+3).at(j))/11;
        }


        //Set pore wall boundary solutions
        for(unsigned int i = 0; i < M +1; i++) //# j<M
        {
            int j = N;
            alfa=(i)*h;
            lambda=(j)*k;

            r=(outerRadius-innerRadius)*alfa+innerRadius;

            gama= (rel_permittivity*permittivity_0*100)/(Helm*r*teta);

            PorePotential[New].at(i).at(j) = (gama*(18*PorePotential[New].at(i).at(j-1)-9*PorePotential[New].at(i).at(j-2)+
                    2*PorePotential[New].at(i).at(j-3))+6*k*del_u)/(6*k+11*gama);


            j=0;

            PorePotential[New].at(i).at(j) = (18*PorePotential[New].at(i).at(j+1)-9*PorePotential[New].at(i).at(j+2) + 2*PorePotential[New].at(i).at(j+3))/11;


        }


        //Convergence check
        converged = true;

        for(unsigned int j = 0; j< N + 1; j++){ //rows
            if(not converged)
                break;
            for(unsigned int i = 0; i< M + 1; i++){//columns
                if (std::abs(1-PorePotential[Old].at(i).at(j)/PorePotential[New].at(i).at(j))> pow(k,2)){
                    converged = false;
                    PorePotential[Old] = PorePotential[New];
                    break;
                }


            }
        }
    }

    if(not  converged){
        error_msg = "Function '"+ std::string(__PRETTY_FUNCTION__) + "' failed to converge!";
    }

    return converged;

}

//---------------------------------------------------------------------------
bool
NAME::WaterConicalPoreAgglomerate::solveO2(std::string& error_msg){

    unsigned int steps = 0; unsigned int max_steps = 150000;
    double h = 1.0/M;
    double k = 1.0/N;
    double kisi = h/k;
    double r,gama,beta_0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11;
    bool converged = false;


    eta =-0.68; //TODO: remove
    c2 = 1.0;





    double j_0_agg = 3.836906e-2;//TODO: remove


    while(not converged and steps < max_steps){
        steps++;


        for(unsigned int i = 1; i< M; i++) //Columns
        {

            for(unsigned int j = 1; j< N; j++) //Rows
            {
                alfa = (i)*h;
                lambda = (j )*k;

                r = (outerRadius - innerRadius)*alfa + innerRadius;
                c1 = pow((outerRadius - innerRadius)/(teta*r),2);
                c3 = pow(outerRadius - innerRadius,2)/(teta*tan(lambda*teta)*pow(r,2));
                c4 = 2*(outerRadius - innerRadius)/r;
                c7 = 2*c1*pow(kisi,2) + 2*c2;
                c8 = c1*pow(kisi,2) + c3*h*kisi/2;
                c9 = c1*pow(kisi,2) - c3*h*kisi/2;
                c10 = c2 + c4*h/2;
                c11 = c2 - c4*h/2;

                cO2[New].at(i).at(j) = (c8*cO2[Old].at(i).at(j+1) + c9*cO2[New].at(i).at(j-1) + c10*cO2[Old].at(i+1).at(j) + c11*cO2[New].at(i-1).at(j))/c7;
            }
        }


        //Set the boundary solution
        for(unsigned int j = 1; j < N +1 ; j++)//# j<N
        {
            int i = 0;
            cO2[New].at(i).at(j) =  (18*cO2[New].at(i+1).at(j) - 9*cO2[New].at(i+2).at(j) + 2*cO2[New].at(i+3).at(j))/11;
        }




        //Set pore wall boundary solutions
        for(unsigned int i = 0; i < M +1; i++) //# j<M
        {
            int j = N;
            alfa=(i)*h;
            lambda=(j)*k;

            r=(outerRadius-innerRadius)*alfa+innerRadius;

            cHwall.at(i) = exp(-PorePotential[New].at(i).at(N) + u_s);
            etaWall.at(i) =  eta -(PorePotential[New].at(i).at(N)- u_s)*R*T/F; //Change code to not calculate eta




            beta_0 = r*teta/(4*F*D_R_N*c_R);



            jWall.at(i) = (j_0_agg*(pow(cHwall.at(i), 1.0))*
                    +  exp(-0.5*F*etaWall.at(i)/(R*T)));//TODO: Rework

            cO2[New].at(i).at(j) = (18*cO2[New].at(i).at(j-1)-9*cO2[New].at(i).at(j-2)+
                    +   2*cO2[New].at(i).at(j-3))/(11+6*k*beta_0*jWall.at(i));//TODO: Rework



            j=0;

            cO2[New].at(i).at(j)  =(18*cO2[New].at(i).at(j+1)-9*cO2[New].at(i).at(j+2)+
                    +   2*cO2[New].at(i).at(j+3))/11;

        }



        //Convergence check
        converged = true;
        for(unsigned int j = 0; j< N + 1; j++){ //rows
            if(not converged)
                break;
            for(unsigned int i = 0; i< M + 1; i++){//columns

                if (std::abs((1 -cO2[Old].at(i).at(j)/cO2[New].at(i).at(j)) >pow(k,3.5))){
                    converged = false;
                    cO2[Old] = cO2[New];
                    break;
                }


            }
        }
    }
    if(not  converged){
        error_msg = "Function '"+ std::string(__PRETTY_FUNCTION__) + "' failed to converge!";
    }

    return converged;
}

