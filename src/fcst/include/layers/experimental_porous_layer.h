// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: experimental_porous_layer.h
// - Description: This class describes a porous layer
// - Developers: Valentin N. Zingan, University of Alberta
// - Id: $Id: experimental_porous_layer.h 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#ifndef _FCST_FUELCELLSHOP_LAYER_EXPERIMENTAL_POROUS_LAYER_H_
#define _FCST_FUELCELLSHOP_LAYER_EXPERIMENTAL_POROUS_LAYER_H_

#include <layers/base_layer.h>
#include <materials/experimental_fluid.h>
#include <materials/GasMixture.h>
#include <materials/experimental_solid.h>

namespace FuelCellShop
{
    namespace Layer
    {

        /**
        * This class describes
        * a porous layer and stores
        * pointers to
        *
        * - \p ExperimentalFluid object,
        * - \p GasMixture object,
        * - \p ExperimentalSolid object.
        *
        * This class is created to experiment with
        * different forms and distributions of porosity,
        * permeability, and tortuosity over porous layer.
        *
        * Let us define either "porosity" or
        * "permeability" or "tortuosity" as "X".
        *
        * This class utilizes two different modes:
        *
        * - \p X_is_constant = \p true and \p X itself
        *   is explicitly specified in parameters file.
        *   In this case, use get_X(std::vector<TYPE>& dst) function
        *   to get \p X in quadrature points of a mesh entity.
        *
        * - \p X_is_constant = \p false and \p X itself
        *   is NOT explicitly specified in parameters file.
        *   In this case, use get_X(std::vector<TYPE>& dst, const std::vector< Point<dim> >& points) function
        *   to get \p X in quadrature points of a mesh entity.
        *
        * \note The function
        * get_X(std::vector<TYPE>& dst, const std::vector< Point<dim> >& points)
        * must be implemented beforehand.
        *
        * The functionality of
        * this class can be extended
        * if needed.
        *
        * \author Valentin N. Zingan, 2012
        */

        template<int dim>
        class ExperimentalPorousLayer : public BaseLayer<dim>
        {
            public:

                ///@name Constructors, destructor, and initialization
                //@{

                /**
                * Constructor.
                */
                ExperimentalPorousLayer(const std::string& name);

                /**
                * Constructor.
                */
                ExperimentalPorousLayer(const std::string&                         name,
                                        FuelCellShop::Material::ExperimentalFluid& fluid,
                                        FuelCellShop::Material::ExperimentalSolid& solid);

                /**
                * Constructor.
                */
                ExperimentalPorousLayer(const std::string&                         name,
                                        FuelCellShop::Material::GasMixture&        gas_mixture,
                                        FuelCellShop::Material::ExperimentalSolid& solid);

                /**
                * Destructor.
                */
                virtual ~ExperimentalPorousLayer();

                /**
                * Initialize
                *
                * - \p fluid,
                * - \p solid.
                */
                void initialize(FuelCellShop::Material::ExperimentalFluid& rfluid,
                                FuelCellShop::Material::ExperimentalSolid& rsolid)
                {
                    fluid       = &rfluid;
                    gas_mixture = nullptr;
                    solid       = &rsolid;
                }

                /**
                * Initialize
                *
                * - \p gas_mixture,
                * - \p solid.
                */
                void initialize(FuelCellShop::Material::GasMixture&        rgas_mixture,
                                FuelCellShop::Material::ExperimentalSolid& rsolid)
                {
                    fluid       = nullptr;
                    gas_mixture = &rgas_mixture;
                    solid       = &rsolid;
                }

                /**
                * Initialize
                *
                * - \p porosity_is_constant,
                * - \p permeability_is_constant,
                * - \p tortuosity_is_constant.
                */
                void initialize(const bool& rporosity_is_constant,
                                const bool& rpermeability_is_constant,
                                const bool& rtortuosity_is_constant)
                {
                    porosity_is_constant     = rporosity_is_constant;
                    permeability_is_constant = rpermeability_is_constant;
                    tortuosity_is_constant   = rtortuosity_is_constant;
                }

                /**
                * Declare parameters.
                */
                virtual void declare_parameters(ParameterHandler& param) const;

                /**
                * Initialize parameters.
                */
                virtual void initialize(ParameterHandler& param);

                //@}

                ///@name Accessors and info
                //@{

                /**
                * This function returns
                * \p fluid.
                */
                const FuelCellShop::Material::ExperimentalFluid* const get_fluid() const
                {
                    return fluid;
                }

                /**
                * This function returns
                * \p gas_mixture.
                */
                const FuelCellShop::Material::GasMixture* const get_gas_mixture() const
                {
                    return gas_mixture;
                }

                /**
                * This function returns
                * \p solid.
                */
                const FuelCellShop::Material::ExperimentalSolid* const get_solid() const
                {
                    return solid;
                }

                /**
                * This function returns
                * \p porosity_is_constant.
                */
                const bool& get_porosity_is_constant() const
                {
                    return porosity_is_constant;
                }

                /**
                * This function returns
                * \p permeability_is_constant.
                */
                const bool& get_permeability_is_constant() const
                {
                    return permeability_is_constant;
                }

                /**
                * This function returns
                * \p tortuosity_is_constant.
                */
                const bool& get_tortuosity_is_constant() const
                {
                    return tortuosity_is_constant;
                }

                /**
                * This function computes
                * constant porosity in quadrature points of a mesh entity.
                */
                void get_porosity(std::vector<double>& dst) const;

                /**
                * This function computes
                * variable porosity in quadrature points of a mesh entity.
                */
                void get_porosity(std::vector<double>&             dst,
                                    const std::vector< Point<dim> >& points) const;

                /**
                * This function computes
                * constant permeability in quadrature points of a mesh entity.
                */
                void get_permeability(std::vector< SymmetricTensor<2,dim> >& dst) const;

                /**
                * This function computes
                * variable permeability in quadrature points of a mesh entity.
                */
                void get_permeability(std::vector< SymmetricTensor<2,dim> >& dst,
                                        const std::vector< Point<dim> >&       points) const;

                /**
                * This function computes
                * square root of constant permeability in quadrature points of a mesh entity.
                */
                void get_SQRT_permeability(std::vector< SymmetricTensor<2,dim> >& dst) const;

                /**
                * This function computes
                * square root of variable permeability in quadrature points of a mesh entity.
                */
                void get_SQRT_permeability(std::vector< SymmetricTensor<2,dim> >& dst,
                                            const std::vector< Point<dim> >&       points) const;

                /**
                * This function computes
                * inverse of constant permeability in quadrature points of a mesh entity.
                */
                void get_permeability_INV(std::vector< SymmetricTensor<2,dim> >& dst) const;

                /**
                * This function computes
                * inverse of variable permeability in quadrature points of a mesh entity.
                */
                void get_permeability_INV(std::vector< SymmetricTensor<2,dim> >& dst,
                                            const std::vector< Point<dim> >&       points) const;

                /**
                * This function computes
                * inverse of square root of constant permeability in quadrature points of a mesh entity.
                */
                void get_SQRT_permeability_INV(std::vector< SymmetricTensor<2,dim> >& dst) const;

                /**
                * This function computes
                * inverse of square root of variable permeability in quadrature points of a mesh entity.
                */
                void get_SQRT_permeability_INV(std::vector< SymmetricTensor<2,dim> >& dst,
                                                const std::vector< Point<dim> >&       points) const;

                /**
                * This function computes
                * constant Forchheimer permeability in quadrature points of a mesh entity.
                */
                void get_Forchheimer_permeability(std::vector< SymmetricTensor<2,dim> >& dst) const;

                /**
                * This function computes
                * variable Forchheimer permeability in quadrature points of a mesh entity.
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
                * This function returns
                * \p typeid of this class.
                *
                * All classes derived
                * from this class
                * must share
                * the same
                * \p typeid
                * for dynamic
                * swapping.
                */
                const std::type_info& get_base_type() const
                {
                    return typeid(ExperimentalPorousLayer<dim>);
                }

                /**
                * This function prints out
                * the layer properties.
                */
                virtual void print_layer_properties() const;

                //@}

            protected:

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

                //////////
                // DATA //
                //////////

                ///@name Layer properties
                //@{

                /**
                * - incompressible,
                * - isothermal,
                * - single-phase,
                * - single-component
                *
                * fluid.
                */
                FuelCellShop::Material::ExperimentalFluid* fluid;

                /**
                * Gas mixture.
                */
                FuelCellShop::Material::GasMixture* gas_mixture;

                /**
                * Solid.
                */
                FuelCellShop::Material::ExperimentalSolid* solid;

                /**
                * Variable defining
                * if the porosity is constant.
                */
                bool porosity_is_constant;

                /**
                * Variable defining
                * if the permeability is constant.
                */
                bool permeability_is_constant;

                /**
                * Variable defining
                * if the tortuosity is constant.
                */
                bool tortuosity_is_constant;

                /**
                * User defined constant porosity.
                */
                double porosity;

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

                //@}

        };

    } // Layer

} // FuelCellShop

#endif