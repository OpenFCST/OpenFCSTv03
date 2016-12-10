//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: porous_layer.h
//    - Description: Child of base layer that implements common functionality for handling gases.
//    - Developers: M. Secanell, Madhur Bhaiya and Valentin Zingan
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__POROUS__LAYER_H
#define _FUELCELLSHOP__POROUS__LAYER_H

//Include STL
#include <algorithm>
#include <boost/concept_check.hpp>

// FCST
#include <microscale/PSD_base.h>
#include <microscale/PSD_HI.h>
#include <microscale/PSD_HO.h>
#include <microscale/PSD_dual.h>
#include <microscale/PSD_none.h>
#include <layers/base_layer.h>
#include <materials/PureGas.h>
#include <materials/GasMixture.h>
#include <utils/fcst_constants.h>
#include <utils/fcst_units.h>





using namespace dealii;

namespace FuelCellShop
{
namespace Layer
{
    
    /**
     * Virtual class used to implement properties that are characteristic of a porous layer.
     * This layer is used to setup information needed to compute molecular diffusivities and Knudsen diffusivities (pending implementation)
     * for any gases. Therefore, it provides member functions to specify the gases that exist within the layer as well as their
     * temperature and pressure.
     *
     * It also inherits all properties of BaseLayer.
     *
     * Note you cannot make an object of this class since it still lacks many member functions that are necessary for a full
     * layer implementation.
     *
     * With respect to gases: you can either use \p gases or \p gas_mixture. Whatever you prefer.
     * However, \p gas_mixture offers more methods in its structure. You cannot use both.
     * 
     * <h3>Usage Details:</h3>
     * 
     * In order to use this class, create a child object of this class. Then, initialize the object.
     * At every cell use set pressure and temperature using #set_pressure and #set_temperature to specify the
     * operating conditions at that location. Then, use compute_gas_diffusion to obtain values for 
     * diffusivity.
     * 
     * A similar process would work for other functions.
     *
     * \author M. Secanell, 2013-15
     * \author M. Bhaiya,   2013
     * \author V. Zingan,   2013
     */
    
    template<int dim>
    class PorousLayer : public BaseLayer<dim>
    {
    public:
        
        ///@name Initalization
        //@{
        /**
         * Member function used to store all the gases that are in the pore space in the gas diffusion layer as well as
         * their temperature [\p Kelvin] and total pressure [\p atm].
         * Then, when computing diffusion, the class will return each one of the diffusion coefficients for the gases.
         * For example, if we pass a vector with for gases, say oxygen and nitrogen, effective_gas_diffusivity will return
         * \f$ D_{N_2,O_2} \f$
         * 
         * \deprecated I would like to deprecate this in favor of set_temperature, set_pressure, etc.
         * 
         */
        void set_gases_and_compute (std::vector<FuelCellShop::Material::PureGas*>& gases_in,
                                    const double& pressure_in,
                                    const double& temperature_in);

        /**
         * Member function used to compute bulk diffusion coefficients (and derivatives w.r.t temperature
         * for non-isothermal case and store inside the layer). 
         * 
         * This function takes solute and solvent gas as an input argument (<b>Fick's diffusion</b> model).
         * Before calling this function, pressure [\p atm] and
         * temperature [\p Kelvin] must be set using #set_pressure and #set_temperature method.
         * 
         */
        void compute_gas_diffusion (FuelCellShop::Material::PureGas* solute_gas,
                                    FuelCellShop::Material::PureGas* solvent_gas);

        /**
         * Member function used to store all the gases that are in the pore space in the porous layer.
         * Besides this, it takes total pressure [\p atm] as input.
         * This method can be particulary used for non-isothermal applications when temperature changes
         * in every cell, hence this method doesn't compute the diffusion coefficients but store the relevant constant data.
         * \note Please ensure that in the case of <b>Fick's diffusion</b> model, \p solvent \p gas should be
         * the last object inside the input vector. On not ensuring this may not give any errors, as this function may also
         * be used for multi-component gas transport but for fickian diffusion, it may give wrong results.
         */
        inline void set_gases (std::vector<FuelCellShop::Material::PureGas*>& gases_in,
                               const double& pressure_in)
        {
            Assert(gases_in.size() >= 2, ExcMessage("Number of gases should be more than or equal to two in PorousLayer::set_gases method."));
            this->gases = gases_in;
            this->pressure = pressure_in;
            
            gas_mixture = nullptr;
        }
        
        /**
         * Set \p gas_mixture.
         */
        inline void set_gas_mixture(FuelCellShop::Material::GasMixture& rgas_mixture)
        {
            gas_mixture = &rgas_mixture;
            gases.clear();
        }
        
        /**
         * Set
         * - \p porosity_is_constant,
         * - \p permeability_is_constant,
         * - \p tortuosity_is_constant.
         */
        inline void set_porosity_permeability_tortuosity_booleans(const bool& rporosity_is_constant,
                                                                  const bool& rpermeability_is_constant,
                                                                  const bool& rtortuosity_is_constant)
        {
            porosity_is_constant     = rporosity_is_constant;
            permeability_is_constant = rpermeability_is_constant;
            tortuosity_is_constant   = rtortuosity_is_constant;
        }
        
        /**
         * Member function used to set the temperature [\p Kelvin] at every quadrature point
         * inside the cell. This function should particulary be used in the case of non-isothermal application.
         */
        inline void set_pressure (const SolutionVariable& p_in)
        {
            Assert( p_in.get_variablename() == total_pressure, ExcMessage("Wrong solution variable passed in PorousLayer::set_pressure.") );
            this->p_vector = p_in;
        }
        
        /**
         * Member function used to set the temperature [\p Kelvin] at every quadrature point
         * inside the cell. This function should particulary be used in the case of non-isothermal application.
         */
        inline void set_temperature (const SolutionVariable& T_in)
        {
            Assert( T_in.get_variablename() == temperature_of_REV, ExcMessage("Wrong solution variable passed in PorousLayer::set_temperature.") );
            this->T_vector = T_in;
        }
        
        /**
         * Member function used to set the liquid water saturation at every quadrature point inside
         * the cell. This function should particulary be used in the case of multi-phase application.
         */
        inline void set_saturation (const SolutionVariable& s_in)
        {
            Assert( s_in.get_variablename() == liquid_water_saturation, ExcMessage("Wrong solution variable passed in PorousLayer::set_saturation.") );
            this->s_vector = s_in;
        }        
        
        /**
         * Member function used to set the liquid water saturation at every quadrature point inside
         * the cell. This function should particulary be used in the case of multi-phase application.
         */
        inline void set_capillary_pressure (const SolutionVariable& p_in)
        {
            Assert( p_in.get_variablename() == capillary_pressure, ExcMessage("Wrong solution variable passed in PorousLayer::set_capillary_pressure.") );
            this->capillary_pressure_vector = p_in;
        }
        //@}
        
        ///@name Accessors and info
        //@{
        /**
         * Return the FuelCellShop::Material::PureGas  pointer that is stored inside the class in the ith position.
         * The diffusion coefficient D[i][j] would be the diffusion coefficient of gas type i with j.
         */
        FuelCellShop::Material::PureGas* get_gas_pointer(int index) const
        {
            Assert(index > 0 || index < gases.size(), ExcIndexRange(index,0,gases.size()));
            return gases[index];
        };
        
        /**
         * Returns the vector of FuelCellShop::Material::PureGas pointers stored in the porous layer.
         * This method is useful to access all gases which are being transported inside a particular porous layer object.
         */
        std::vector<FuelCellShop::Material::PureGas*> get_gases() const {
            return this->gases; }
        
        /**
         * This function returns \p gas_mixture.
         */
        const FuelCellShop::Material::GasMixture* const get_gas_mixture() const {
            return this->gas_mixture; }
        
        /** 
         * Return the gas index in the GDL class. If the gas does not exist in the GDL it will return -1 as the index.  
         */
        void get_gas_index(FuelCellShop::Material::PureGas* gas_type, int& index) const;

        /** 
         * Return the constant temperature [\p Kelvin] and constant pressure [\p atm] inside the layer. 
         * 
         * \deprecated
         */
        void get_T_and_p(double &T, double &p) const {
            T = this->temperature;
            p = this->pressure;
        }

        /** 
         * Return the constant pressure [\p atm] inside the layer. 
         * 
         * \deprecated
         */
        void get_p(double& p) const {
            p = this->pressure; }

        /**
         * This function returns \p porosity_is_constant.  
         */
        const bool& get_porosity_is_constant() const 
        {
            return porosity_is_constant; 
        }
        
        /** 
         * This function returns \p permeability_is_constant.
         */
        const bool& get_permeability_is_constant() const  {
            return permeability_is_constant; }
        
        /** 
         * This function returns \p tortuosity_is_constant. 
         */
        const bool& get_tortuosity_is_constant() const  {
            return tortuosity_is_constant; }
        
        /** 
         * This function computes constant porosity in quadrature points of a mesh entity. 
         */
        inline double get_porosity() const {
            AssertThrow( porosity_is_constant , ExcInternalError() );
            return porosity;
        }
        
        /** 
         * This function computes constant porosity in quadrature points of a mesh entity. 
         */
        inline void get_porosity(std::vector<double>& dst) const {
            AssertThrow( porosity_is_constant , ExcInternalError() );
            
            for(unsigned int q = 0; q < dst.size(); ++q)
                dst[q] = this->get_porosity();
        }
        
        /** 
         * This function computes variable porosity in quadrature points of a mesh entity.  
         */
        void get_porosity(std::vector<double>&             dst,
                          const std::vector< Point<dim> >& points) const;
                          
        /** 
         * This function computes constant permeability in quadrature points of a mesh entity. 
         */
        void get_permeability(std::vector< SymmetricTensor<2,dim> >& dst) const;
        
        /**
         * This function computes variable permeability in quadrature points of a mesh entity. 
         */
        void get_permeability(std::vector< SymmetricTensor<2,dim> >& dst, 
                              const std::vector< Point<dim> >&       points) const;
        
        /** 
         * This function computes square root of constant permeability in quadrature points of a mesh entity. 
         */
        void get_SQRT_permeability(std::vector< SymmetricTensor<2,dim> >& dst) const;
        
        /** 
         * This function computes square root of variable permeability in quadrature points of a mesh entity. 
         */
        void get_SQRT_permeability(std::vector< SymmetricTensor<2,dim> >& dst, 
                                   const std::vector< Point<dim> >&       points) const;
        
        /**
         * This function computes inverse of constant permeability in quadrature points of a mesh entity.
         */
        void get_permeability_INV(std::vector< SymmetricTensor<2,dim> >& dst) const;
        
        /**
         * This function computes inverse of variable permeability in quadrature points of a mesh entity.
         */
        void get_permeability_INV(std::vector< SymmetricTensor<2,dim> >& dst, 
                                  const std::vector< Point<dim> >&       points) const;
        
        /**
         * This function computes inverse of square root of constant permeability in quadrature points of a mesh entity.
         */
        void get_SQRT_permeability_INV(std::vector< SymmetricTensor<2,dim> >& dst) const;
        
        /**
         * This function computes inverse of square root of variable permeability in quadrature points of a mesh entity.
         */
        void get_SQRT_permeability_INV(std::vector< SymmetricTensor<2,dim> >& dst, 
                                       const std::vector< Point<dim> >&       points) const;
        
        /**
         * This function computes constant Forchheimer permeability in quadrature points of a mesh entity.
         */
        void get_Forchheimer_permeability(std::vector< SymmetricTensor<2,dim> >& dst) const;
        
        /**
         * This function computes variable Forchheimer permeability in quadrature points of a mesh entity.
         */
        void get_Forchheimer_permeability(std::vector< SymmetricTensor<2,dim> >& dst, 
                                          const std::vector< Point<dim> >&       points) const;
        
        /**
         * This function computes
         * constant tortuosity in quadrature points of a mesh entity.
         */
        void get_tortuosity(std::vector< SymmetricTensor<2,dim> >& dst) const;
        
        /**
         * This function computes
         * variable tortuosity in quadrature points of a mesh entity.
         */
        void get_tortuosity(std::vector< SymmetricTensor<2,dim> >& dst, 
                            const std::vector< Point<dim> >&       points) const;
                
        /**
         * Return the effective diffusivity [\p m^2/s] in the GDL. 
         * 
         * It transforms bulk diffusivity, computed using #compute_gas_diffusion method 
         * and transforms it into an effective property,
         * taking into account the porosity, saturation and GDL structure (Anisotropic case),
         * at all quadrature points of the cell.
         * 
         * \warning Before calling this class, make sure that you have called #compute_gas_diffusion, otherwise the
         * output will be incorrect
         * 
         */
        virtual void effective_gas_diffusivity(std::vector< Tensor<2,dim> >&) const 
        {};
        
        /**
         * Return the effective diffusivty in the GDL for all the gases assigned to the
         * layer using set_gases_and_compute. This routine uses FuelCellShop::Mixture::ChapmanEnskog to compute
         * the binary diffusivity for each gas with respect to each other. If the layer contains three gases the 
         * vector returns D12, D13, D23. For 4 gases, it returns D12, D13, D14, D23, D24, D34.
         * 
         * \deprecated Use #compute_gas_diffusion with appropriate gases and then #effective_gas_diffusivity
         */
        virtual void effective_gas_diffusivity(Table<2, double>& D_eff ) const
        {
            const std::type_info& info = typeid(*this);
            FcstUtilities::log << "Pure function " << __FUNCTION__
            << " called in Class "
            << info.name() << std::endl;
        };
        
        /**
         * Return a tensor with the effective diffusivty in the GDL for all the gases assigned to the
         * layer using set_gases_and_compute. This routine uses FuelCellShop::Mixture::ChapmanEnskog to compute
         * the binary diffusivity for each gas with respect to each other. If the layer contains three gases the 
         * vector returns D12, D13, D23. For 4 gases, it returns D12, D13, D14, D23, D24, D34 
         * (Anisotropic case).
         * 
         * \deprecated Use #compute_gas_diffusion with appropriate gases and then #effective_gas_diffusivity
         */
        virtual void effective_gas_diffusivity(Table< 2, Tensor<2,dim> > &D_eff ) const
        {
            const std::type_info& info = typeid(*this);
            FcstUtilities::log << "Pure function " << __FUNCTION__
            << " called in Class "
            << info.name() << std::endl;
        };
        
        /**
         * Return the derivative of effective diffusivity w.r.t solution variables/design parameters in the GDL. 
         * 
         * It transforms bulk diffusion derivative properties computed using #compute_gas_diffusion method 
         * and transforms it into an effective property,
         * taking into account the porosity, saturation and GDL structure (Anisotropic case),
         * at all quadrature points of the cell.
         * 
         * \warning Before calling this class, make sure that you have called #compute_gas_diffusion, otherwise the
         * output will be incorrect
         */
        virtual void derivative_effective_gas_diffusivity(std::map< VariableNames, std::vector< Tensor<2,dim> > >&) const 
        {};
        
        /**
         * Member function used to compute diffusion for a solute_gas, solvent_gas combination
         * at a given temperature and pressure.
         *
         * The diffusion coefficient is returned in SI units, i.e., m^2/s.
         * 
         * The effective diffusivity returned is either the molecular diffusivity or, if the flag
         * "set Use Bosanquet approx." use set to true, using the Bousinesq approximation is used in order 
         * to compute the bulk combined diffusivity such that
         * \f[
         * \frac{1}{D} = \frac{1}{D_{i,j}} + \frac{1}{D^K_i}
         * \f]
         * where \f$ D_{i,j}\f$ is the molecular diffusion coefficient obtained using effective_molecular_gas_diffusion_coefficient
         * and \f$ D^K_i \f$ is the Knudsen diffusivity computed using Knudsen_diffusion
         * 
         * By defaul the diffusion coefficient is computed using the pore radius in:
         * @code
         * subsection Fuel Cell Data
         *   subsection layer_name
         *     set Use Bosanquet approx. = true
         *     set Knudsen pore radius [um] = 1e6
         *   end
         * end
         * @endcode
         * 
         * If you would like to neglect the diffusion coefficient, simply make this value very large (the default 
         * basically results on negligible Knudsen effects -- the value is: 161.0949 (SI) compared to 1e-5 (SI) for molecular diff for gases
         * -- Note that you add the inverse of these values).
         * 
         * @note If you would like
         * 
         * References
         * 
         * - [1] Pollard and Present, On gaseous self-diffusion in long capillary tubes, Physical Review, Vol 73, No. 7, 1948
         * - [2] Burganos, Gas diffusion in random binary media, Journal of Chemical Physics, 109(16), 1998
         * 
         * @author M. Secanell, 2014
         */
         virtual void gas_diffusion_coefficient(std::vector<double>& D_b) const
         {
             D_b = this->D_bulk;
         }
         
        /**
         * Member function used to compute diffusion for a solute_gas, solvent_gas combination
         * at a given temperature and pressure.
         *
         * The diffusion coefficient is returned in SI units, i.e., m^2/s.
         * 
         * The effective diffusivity returned is either the molecular diffusivity or, if the flag
         * "set Use Bosanquet approx." use set to true, using the Bousinesq approximation is used in order 
         * to compute the bulk combined diffusivity such that
         * \f[
         * \frac{1}{D} = \frac{1}{D_{i,j}} + \frac{1}{D^K_i}
         * \f]
         * where \f$ D_{i,j}\f$ is the molecular diffusion coefficient obtained using effective_molecular_gas_diffusion_coefficient
         * and \f$ D^K_i \f$ is the Knudsen diffusivity computed using Knudsen_diffusion
         * 
         * References
         * - [1] Pollard and Present, On gaseous self-diffusion in long capillary tubes, Physical Review, Vol 73, No. 7, 1948 
         * - [2] Burganos, Gas diffusion in random binary media, Journal of Chemical Physics, 109(16), 1998
         * 
         * @author M. Secanell, 2014
         */
         virtual void gas_diffusion_coefficient(std::vector<double>& D_b,
                                                std::vector<double>& dD_b_dT) const
         {
             D_b = this->D_bulk;
             dD_b_dT = this->dD_bulk_dT;
         }

        /**
         * Member function used to compute molecular diffusion for a solute_gas, solvent_gas combination
         * at a given temperature and pressure.
         * 
         * The molecular diffusivity is returned in units of m^2/s
         * 
         * \warning Before calling this member function the following class should have been called: #compute_gas_diffusion
         * 
         * @author M. Secanell, 2014
         */
        void molecular_gas_diffusion_coefficient(std::vector<double>& D_m) const         
        {
            D_m = this->D_molecular;
        }
        
        /**
         * Member function used to compute molecular diffusion for a solute_gas, solvent_gas combination
         * at a given temperature and pressure.
         * 
         * The molecular diffusivity is returned in units of m^2/s
         * 
         * \warning Before calling this member function the following class should have been called: #compute_gas_diffusion
         * 
         * @author M. Secanell, 2014
         */
         void molecular_gas_diffusion_coefficient(std::vector<double>& D_m,
                                                  std::vector<double>& dD_m_dT) const
         {
             D_m = this->D_molecular;
             dD_m_dT = this->dD_molecular_dT;
         }
        
        /**
         * Member function used to get the Knudsen diffusivity in the layer after calling
         * #compute_gas_diffusion
         * 
         * @author M. Secanell, 2014
         */
        void Knudsen_diffusion(std::vector<double>& D) const
        {
            D = this->D_k;
        }
                               
        /**
         * Member function used to compute the Knudsen diffusivity in the layer.after calling
         * #compute_gas_diffusion
         * 
         * @author M. Secanell, 2014
         */
        void Knudsen_diffusion(std::vector<double>& D,
                               std::vector<double>& dD_dT) const
        {
             D = this->D_k;
             dD_dT = this->dD_k_dT;
        }
        /**
         * Compute \f$ p_c \quad \left[ dyne \cdot cm^{-2}\right] \f$, at all quadrature points in the cell.
         */
        virtual void pcapillary(std::vector<double>&) const
        {
            const std::type_info& info = typeid(*this);
            FcstUtilities::log << "Calling parent porous layer"<< std::endl; 
        }
        
        virtual void saturation_from_capillary_equation(std::vector<double>&) const
        {
            const std::type_info& info = typeid(*this);
            FcstUtilities::log << "Calling parent porous layer"<< std::endl; 
        }
        
        virtual void derivative_saturation_from_capillary_equation_PSD(std::vector<double>&) const
        {
            const std::type_info& info = typeid(*this);
            FcstUtilities::log << "Calling parent porous layer"<< std::endl; 
        }
        
        virtual void saturated_liquid_permeablity_PSD(double&) const
        {
            const std::type_info& info = typeid(*this);
            FcstUtilities::log << "Calling parent porous layer"<< std::endl; 
        }
        
        virtual void interfacial_surface_area_PSD(std::vector<double>&) const
        {
            const std::type_info& info = typeid(*this);
            FcstUtilities::log << "Calling parent porous layer"<< std::endl; 
        }
        
        virtual void relative_liquid_permeability_PSD(std::vector< Tensor<2,dim> >& ) const
        {
            const std::type_info& info = typeid(*this);
            FcstUtilities::log << "Calling parent porous layer"<< std::endl; 
        }
        
        virtual void derivative_relative_liquid_permeablity_PSD(std::vector<double>& ) const
        {
            const std::type_info& info = typeid(*this);
            FcstUtilities::log << "Calling parent porous layer"<< std::endl; 
        }
        
        virtual void derivative_interfacial_surface_area_PSD(std::vector<double>& ) const
        {
            const std::type_info& info = typeid(*this);
            FcstUtilities::log << "Calling parent porous layer"<< std::endl; 
        }
        //@}     
        ///@name Computing member functions
        //@{
        /**
         * Member function used to compute the Knudsen diffusivity in the layer.
         * - solute_gas: PureGas for which the Knudsen number should be computed.
         * - T_in: temperature in Kelvin.
         * - D_k: Knudsen diffusivity in m^2/s
         * 
         * @author M. Secanell, 2014
         */
        void compute_Knudsen_diffusion(const FuelCellShop::Material::PureGas* solute_gas,
                                       const SolutionVariable& T_in,
                                       std::vector<double>& D_k) const;
                               
        /**
         * Member function used to compute the Knudsen diffusivity in the layer.
         * - solute_gas: PureGas for which the Knudsen number should be computed.
         * - T_in: temperature in Kelvin.
         * - D_k: Knudsen diffusivity in m^2/s
         * - dD_k_dT: Knudsen diffusivity derivative vs. temperature in m^2/(s K)
         * 
         * The equation used to compute Knudsen diffusivity is
         * \f[
         * D_{k,i} = \frac{2.0 r_k 10^{-6}}{3.0} \sqrt{\frac{8.0RT}{\pi M_i}}
         * \f]
         * where R is 8.3144 J/mol K in SI units, molecular weight M_{i} is expressed in units of kg/mol and temperature T has units of K.
         * 
         * @author M. Secanell, 2014
         */
        void compute_Knudsen_diffusion(const FuelCellShop::Material::PureGas* solute_gas,
                                       const SolutionVariable& T_in,
                                       std::vector<double>& D_k,
                                       std::vector<double>& dD_k_dT) const;
            
        //@}
        ///@name Debugging and info                    
        /**
         * This member function is a virtual class that can be used to output to screen information from the layer.
         *
         * This function should be re-implemented in each layer with data that is relevant in each case.
         */
        virtual void print_layer_properties() const;
        //@}
        
    protected:
        
        ///@name Constructors, destructor, and initialization
        //@{
        
        /**
         * Constructor.
         */
        PorousLayer(const std::string& name)
        :
        FuelCellShop::Layer::BaseLayer<dim>(name)
        {
            this->derivative_flags.push_back(total_pressure);
            this->derivative_flags.push_back(temperature_of_REV);
            porosity_is_constant     = true;
            permeability_is_constant = true;
            tortuosity_is_constant   = true;
            gases.clear();
            gas_mixture = nullptr;
            PSD = nullptr;
            psd_pointer = nullptr;
        }
        
        /**
         * Constructor.
         *
         * \warning For internal use only.
         */
        PorousLayer()
        :
        FuelCellShop::Layer::BaseLayer<dim>()
        {}
        
        /**
         * Constructor.
         */
        PorousLayer(const std::string&                  name,
                    FuelCellShop::Material::GasMixture& gas_mixture)
        :
        FuelCellShop::Layer::BaseLayer<dim>(name),
        gas_mixture(&gas_mixture)
        {
            this->derivative_flags.push_back(total_pressure);
            this->derivative_flags.push_back(temperature_of_REV);
            porosity_is_constant     = true;
            permeability_is_constant = true;
            tortuosity_is_constant   = true;
            gases.clear();
            PSD = nullptr;
            psd_pointer = nullptr;
        }
        
        /**
         * Destructor.
         */
        virtual ~PorousLayer()
        {}
        
        /**
         * Declare parameters for a parameter file.
         */
        virtual void declare_parameters (const std::string &name, ParameterHandler &param) const;
        
        /**
         * Declare parameters for a parameter file
         *
         * \note Since there are two declare parameters, this one is hidden by the former, so it has to be
         * redefined otherwise it cannot be called :(
         */
        virtual void declare_parameters (ParameterHandler &param) const
        {
            declare_parameters(this->name, param);
        }
        
        /**
         * Member function used to read in data and initialize the necessary data
         * to compute the coefficients.
         */
        virtual void initialize (ParameterHandler &param);
        //@}
        
        ///@name Minor functions
        //@{
        
        /**
         * This function is used
         * to print out the name of another function
         * that has been declared in the scope of this class,
         * but not yet been implemented.
         */
        void print_caller_name(const std::string& caller_name) const;
                               
        //@}
        ///@name Internal accessors and info
        //@{
        /**
         * Return the molecular diffusivty all the gases assigned to the
         * layer using set_gases_and_compute. This routine uses FuelCellShop::Mixture::ChapmanEnskog to compute
         * the binary diffusivity for each gas with respect to each other. If the layer contains three gases the
         * vector returns D12, D13, D23. For 4 gases, it returns D12, D13, D14, D23, D24, D34.
         * 
         * \deprecated DO NOT USE
         */
        virtual void gas_diffusion_coefficients(Table< 2, double > &) const;
        
        /**
         * Return the derivative of the molecular diffusion coefficient
         * with respect to the derivative flags for all the gases assigned to the
         * layer using set_gases_and_compute. This routine uses FuelCellShop::Mixture::ChapmanEnskog to compute
         * the binary diffusivity for each gas with respect to each other. If the layer contains three gases the
         * vector returns D12, D13, D23. For 4 gases, it returns D12, D13, D14, D23, D24, D34 .
         * Derivatives are provided in the vector w.r.t following parameters in the same order as below:
         * - #total_pressure
         * - #temperature_of_REV
         * 
         * \deprecated DO NOT USE
         */
        virtual void derivative_gas_diffusion_coefficients(std::vector< Table< 2, double > >&) const;
        
        //@}
        
        //////////
        // DATA //
        //////////
        
        ///@name Layer properties
        //@{
        
        /** Gas mixture. */
        FuelCellShop::Material::GasMixture* gas_mixture;
        
        /**
         * Gases inside a porous layer.
         * \deprecated. Use \p gas_mixture instead.
         */
        std::vector<FuelCellShop::Material::PureGas*> gases;
        
        /** Variable defining if the porosity is constant. */
        bool porosity_is_constant;
        
        /** Variable defining if the permeability is constant. */
        bool permeability_is_constant;
        
        /** Variable defining if the tortuosity is constant. */
        bool tortuosity_is_constant;
        
        /** User defined constant porosity. */
        double porosity;
        
        /** Boolean flag that specifies if Knudsen effects should be accounted for.
         * This value is specified in Generic>>Use Bosanquet approx. Its default value is false.
         */
        bool use_Bosanquet;
        
        /**
         * Parameter used to define Knudsen pore radius
         */
        double Knudsen_radius;
        
        /**
         * User defined constant permeability, m^2.
         */
        SymmetricTensor<2,dim> permeability;
        
        /**
         * Square root of user defined constant permeability, m.
         */
        SymmetricTensor<2,dim> SQRT_permeability;
        
        /**
         * Inverse of user defined constant permeability, 1/m^2.
         */
        SymmetricTensor<2,dim> permeability_INV;
        
        /**
         * Inverse of square root of user defined constant permeability, 1/m.
         */
        SymmetricTensor<2,dim> SQRT_permeability_INV;
        
        /**
         * User defined constant Forchheimer permeability, 1/m.
         */
        SymmetricTensor<2,dim> Forchheimer_permeability;
        
        /**
         * User defined constant tortuosity.
         */
        SymmetricTensor<2,dim> tortuosity;
        
        /**
         * If GDL properties are stored inside the class (\em e.g DummyGDL) then, return the property
         * stored under coefficient_name name.
         */
        std::string diffusion_species_name;
        
        /**
         * Temperature [\p K] used to compute gas diffusivity.
         * \Deprecated. Use \p temperature stored in \p gas_mixture instead. Also, you should be using T_vector
         */
        double temperature;
        
        /**
         * Total pressure [\p atm] used to compute gas diffusivity.
         * \Deprecated. Use \p total_pressure stored in \p gas_mixture instead. Also you should be using p_vector.
         */
        double pressure;

        /** Pressure at every quadrature point inside the cell in [Pa]. */
        SolutionVariable p_vector;
        
        /** Temperature at every quadrature point inside the cell in [K]. */
        SolutionVariable T_vector;
        
        /** Liquid water saturation at every quadrature point inside the cell [-]. */
        SolutionVariable s_vector;
        
        /** Liquid water capillary pressure at every quadrature point inside the cell in [Pa]. */
        SolutionVariable capillary_pressure_vector;
        
        /** Tensor of diffusion coefficients
         * This are computed with setting up the gas so that they do not need to be recomputed all the time.
         */
        Table< 2, double > D_ECtheory;
        
        /** Vector of tensors for the derivative of the diffusion coefficients
         *-- This are computed with setting up the gas so that they do not need to be recomputed all the time.*/
        std::vector< Table< 2, double > > dD_ECtheory_dx;
        
        /**
         * Vector of molecular diffusion coefficients at every quadrature point inside the cell in m^2/s
         */
        std::vector<double> D_molecular;
        
        /**
         * Vector of derivatives for molecular diffusion coefficients w.r.t temperature, at every quadrature  in m^2/s
         */
        std::vector<double> dD_molecular_dT;
        
        /**
         * Vector of Knudsen diffusion coefficients at every quadrature point inside the cell in m^2/s
         */
        std::vector<double> D_k;
        
        /**
         * Vector of derivatives for Knudsen diffusion coefficients w.r.t temperature, at every quadrature  in m^2/s
         */
        std::vector<double> dD_k_dT;
        
        /**
         * Vector of bulk diffusion coefficients at every quadrature point inside the cell.
         * 
         * This might store the molecular diffusivity or the diffusivity computed using Bosanquet approx.
         */
        std::vector<double> D_bulk;
        
        /**
         * Vector of derivative of bulk diffusion coefficients w.r.t temperature, at every quadrature
         * point inside the cell. These are computed using #compute_gas_diffusion method.
         * 
         * This might store the molecular diffusivity or the diffusivity computed using Bosanquet approx.
         */
        std::vector<double> dD_bulk_dT;
        
        /** Boolean flag to specify if a PSD is to be used to estimate saturation, permeability, etc.*/
        bool PSD_is_used;
        
        /** PSD class type from input file */
        std::string PSD_type;
        
        /** Pointer to the PSD object */
        boost::shared_ptr< FuelCellShop::MicroScale::BasePSD<dim> > PSD; //psd_pointer
        
        /** Pointer to the PSD object */
        FuelCellShop::MicroScale::BasePSD<dim> * psd_pointer;

        /** Pointer used to store the solute gas for computing binary diffusion coefficients. This value
         * is used in child classes and is setup using compute_gas_diffusion
         */
        FuelCellShop::Material::PureGas* solute_gas;
        
        /** Pointer used to store the solute gas for computing binary diffusion coefficients. This value
         * is used in child classes and is setup using compute_gas_diffusion
         */
        FuelCellShop::Material::PureGas* solvent_gas;
        
        //@}
        
};

} // Layer

} // FuelCellShop

#endif