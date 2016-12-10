======================================
 Introduction to AppPemfcTPSaturation 
======================================


Introduction
============

AppPemfcTPSaturation implements a steady-state, non-isothermal, two-phase flow model based on a saturation equation. 
This class is used to solve the physical pheonoma on a complete membrane electrode assembly.
The anode hydrogen oxydation reaction is modelled using an aglomerate model with dual-pathway kinetics and the cathode oxygen reduction
reaction using an agglomerate model and the kinetics in Sun et al., EA, 2006. The membrane is modelled using a modified Springer model.

         


Governing equations
===================
 
The model is based on the following assumptions: 

  1. The fuel cell is at steady state and operates at constant pressure (pressure gradients are negligible).
  2. The gas mixtures are assumed to have ideal gas behaviour.
  3. The cathode is fed with humidified air.
  4. The anode is fed with humidified hydrogen.
  5. The gas diffusion layers are composed of a porous fibrous matrix.
  6. The catalyst layer is composed of three phases, viz, platinum supported on carbon, membrane electrolyte ionomer and open (void) space.
  7. In the case of an agglomerate catalyst layer model being considered, the electrochemical reaction occurs inside the agglomerates.
  8. The transport of reactants from the gas channels to the catalyst layer occurs only by diffusion of reactant gas to the agglomerate surface and then by dissolution and diffusion through the ionomer to the reaction site.
  9. The transport of water inside the electrolyte in the membrane and CL is modeled using Springer's model including thermal osmosis effects.
  10. The membrane and gas phase in the CL are assumed to be in equilibrium throughout the CL, therefore they are related by means of the sorption isotherm.
  11. The transport of protons takes place only through the electrolyte, i.e. the Nafion and it is governed by Ohm's law.
  12. The transport of electrons takes place only through the solid phase, i.e. the carbon fibers in the GDL and the mixture of carbon supported platinum in the catalyst layer, and is governed by Ohm's law.
  13. Very long Brinkmann number flow is assumed, lead to negligible viscous dissipation, and hence neglected in the model.
  14. Gas and solid phases are in thermal equilibrium. This assumption is fairly valid since there is very high interstitial surface area, and convective heat transfer between these two phases would be sufficiently large such that temperature of all the phases at a particular location in the fuel cell will be approximately equal.
  15. Convection effects (including enthalpy transport by convection) are negligible for the fuel cell conditions discussed here.

The governing equations are

.. math::
  R_1(\vec{u}) = - \vec{\nabla} \cdot \left( \dfrac{p_T}{RT} ~D^{eff}_{O_2,N_2} \vec{\nabla} x_{O_2} \right)  - S_{O_2}  =  0
.. math::  
  R_2(\vec{u}) = - \vec{\nabla} \cdot \left( \dfrac{p_T}{RT} ~D^{eff}_{H_2O,N_2 ~\text{or} ~H_2}\vec{\nabla} x_{H_2O} \right) - S_{H_2O} = 0
.. math::
  R_3(\vec{u}) = - \vec{\nabla} \cdot \left( \sigma^{eff}_{m} \vec{\nabla} \phi_m \right) - S_{H^+} = 0
.. math::
  R_4(\vec{u}) = \vec{\nabla} \cdot  \left( \sigma^{eff}_{s} \vec{\nabla} \phi_s \right) - S_{e^-} = 0
.. math::
  R_5(\vec{u}) = - \vec{\nabla} \cdot \left( n_d \dfrac{\sigma_m^{eff}}{F} \vec{\nabla} \phi_m +\dfrac{\rho_{dry}}{EW} D^{eff}_{\lambda} \vec{\nabla} \lambda +\dfrac{D^{eff}_T}{M_{H_2O}} \vec{\nabla} T \right) - S_{\lambda}  =  0
.. math::
  R_6(\vec{u}) = - \vec{\nabla} \cdot \left( k^{eff} \vec{\nabla} T \right) + \sum_{\text{gases}, \lambda} \left( \vec{N}_i \cdot \vec{\nabla}\bar{H}_i \right) - S_T  =  0

 R_7(\vec{u}) = - \vec{\nabla} \cdot \left( \dfrac {\rho_l k_{rl}} {\mu_l} \vec{\nabla} s \right) - S_L  =  0
  
where the unknowns are, the oxygen mole fraction, :math:`x_{O_2}`; the water mole fraction, :math:`x_{H_2O}`; the electrolyte (membrane) and electronic potentials, :math:`\phi_m` and :math:`\phi_s` respectively; the membrane water content, :math:`\lambda`; and, the temperature, :math:`T`. The effective transport parameters :math:`D^{eff}_{O_2}`, :math:`D^{eff}_{H_2O}`, :math:`\sigma^{eff}_{m}`, :math:`\sigma^{eff}_{s}`, :math:`D^{eff}_{\lambda}`,  :math:`D^{eff}_T`, and :math:`k^{eff}`,  are different in the membrane, GDL and CL and depend non-linearly on the design variables. Due to the solution methodology, all equations need to be solved in all the domains, i.e. GDL, CL and membrane. However, some equations are not necessary in some of the cell domains. This is addressed by making the unnecessary transport parameters zero. 
  
The source terms in the system of equations are given by

.. math::
  S_{O_2} = \left\{
  \begin{array}{cl}
  \frac{-\nabla \cdot \vec{i}}{4F} \quad &\text{in CCL} \\
  0 \quad &\text{otherwise}
  \end{array}
  \right.

.. math::
  S_{H_2O} = \left\{
  \begin{array}{cl}
  \frac{\nabla \cdot \vec{i}}{2F} - k_t \frac{\rho_{dry}}{EW}(\lambda_{eq}-\lambda) \quad & \text{in CCL} \\
  - k_t \frac{\rho_{dry}}{EW}(\lambda_{eq}-\lambda) \quad & \text{in ACL} \\
  0 \quad & \text{otherwise}
  \end{array}
  \right.

.. math::
  S_{H^+} = \left\{
  \begin{array}{cl}
  -\nabla \cdot \vec{i} \quad & \text{in CCL} \\
  \nabla \cdot \vec{i} \quad & \text{in ACL} \\
  0 \quad & \text{otherwise}
  \end{array}
  \right.
  
.. math::
  S_{e^-} = \left\{
  \begin{array}{cl}
  0 &\text{in GDL and membrane} \\
  -\nabla \cdot \vec{i} \quad & \text{in CCL} \\
  \nabla \cdot \vec{i} \quad & \text{in ACL} \\
  0 \quad & \text{otherwise}
  \end{array}
  \right.
  
.. math::
  S_{\lambda} = \left\{
  \begin{array}{cl}
  k_t \frac{\rho_{dry}}{EW}(\lambda_{eq}-\lambda) \quad & \text{in CLs} \\
  0 \quad & \text{otherwise}
  \end{array}
  \right.

.. math::
  S_{T} = \left\{
  \begin{array}{cl}
  \sigma^{eff}_s ( \vec{\nabla} \phi_s \cdot \vec{\nabla} \phi_s ) \quad & \text{in GDL and MPL} \\~\\ 
  \sigma^{eff}_m ( \vec{\nabla} \phi_m \cdot \vec{\nabla} \phi_m ) \quad & \text{in Membrane} \\~\\
  -(\nabla \cdot \vec{i}) \eta + \frac{\nabla \cdot \vec{i}}{2F} ( -T f_{ORR} \Delta \bar{S}_{overall} - \bar{H}_{lv} ) + \sigma^{eff}_m ( \vec{\nabla} \phi_m \cdot \vec{\nabla} \phi_m ) \\ \qquad + \sigma^{eff}_s ( \vec{\nabla} \phi_s \cdot \vec{\nabla} \phi_s )  +k_t \frac{\rho_{dry}}{EW} ( \lambda_{eq} - \lambda ) \bar{H}_{sorption} \quad & \text{in CCL} \\~\\
  (\nabla \cdot \vec{i}) \eta + \frac{\nabla \cdot \vec{i}}{2F} ( -T (1-f_{ORR}) \Delta \bar{S}_{overall} ) + \sigma^{eff}_m ( \vec{\nabla} \phi_m \cdot \vec{\nabla} \phi_m ) \\ \qquad + \sigma^{eff}_s ( \vec{\nabla} \phi_s \cdot \vec{\nabla} \phi_s )  +k_t\frac{\rho_{dry}}{EW} ( \lambda_{eq} - \lambda ) \bar{H}_{sorption} \quad & \text{in ACL}
  \end{array}
  \right.

and 

.. math::
  S_{\L} = \left\{
  \begin{array}{cl}
  k_{e/c} a_{lv} H_{lv}  \left( \dfrac {p_v - p_{sat}(s,T)} {p_{sat}(s,T)}\right) \\
  0 \quad & \text{PEM}
  \end{array}
  \right.

where :math:`\lambda_{eq}` is given by the sorption isotherm reported by Hinatsu et al. at the corresponding water vapour activity value in the specific location in the CL.
 
Directory structure
=================== 

The AppPemfcTPSaturation directory consists of the following folders:

1. template : This folder contains the default files for running all the examples in the other folders. Please **do not** modify this file as 
it will result in all tests failing. If you would like to create your own example either include this file to your simulation using the *include*
command or copy the file to a different location. 

2. Polarization_curve : This folder contains the two polarization curve studies: dry case and wet case. Both of them contain :code:`main_test.prm` and :code:`data_test.prm` files needed to run a simulation.
Note the data file includes the template find and adds the necessary modifications.


3. Jie_EA2016_article_data : This folder contains the files to reproduce the studies that will be published in JES paper, still await for the review.
Note the data file includes the template find and adds the necessary modifications.

Setting up a AppPemfcTPSaturation simulation
===============================

In order to run OpenFCST, two files are needed that provide the necessary information for OpenFCST to execute:

  a. A main file: This file is used to select the appropriate: a) type of analysis, i.e. analysis, parametric study, polarization curve and optimization study; application; b) the nonlinear solver; c) data file name; and, d) several less critical parameters.

  b. A data file: This file is used to input all the input data used for the simulation for the application selected.

Both these files can either be loaded and modified via the openFCST graphical user interface (GUI) or modified as a text file. 

Obtaining a polarization curve
==============================

Once the `main.xml` or `.prm` file has been modified, run the simulation.

References
==========

.. _Secanell07a:
  
M. Secanell, B. Carnes, A. Suleman and N. Djilali, "Numerical Optimization of Proton Exchange Membrane Fuel Cell Cathode Electrodes", Electrochimica Acta, 52(7):2668-2682, 2007.

.. _Secanell07b:

M. Secanell, K. Karan, A. Suleman and N. Djilali, "Multi-Variable Optimization of PEMFC Cathodes using an Agglomerate Model", Electrochimica Acta, 52(22):6318-6337, 2007.

