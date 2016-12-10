===================================
 Introduction to AppPemfcNIThermal
===================================

Introduction
============

A non-isothermal, single phase membrane electrode assembly model (MEA) mathematical model accounting for most applicable heat sources, 
viz., reversible, irreversile, ohmic heating, phase change, and, heat of sorption/desorption. The thermal transport equation is fully coupled with an MEA model, 
and non-isothermal effects, such as thermal osmosis through the membrane, local relative humidity variations in the catalyst layers, and water sorption into the membrane, 
are considered in detail. Different catalyst layer models can be evaluated including:
  a. macro-homogeneous catalyst layer model
  b. agglomerate catalyst layer model.

For more information about the model description, please see BhaiyaEA14_. If using this application, please cite the article BhaiyaEA14_.

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
  15. Two-phase flow and condensation are not considered in this model, therefore the gas phase is allowed to be in super-saturated form.
  16. Convection effects (including enthalpy transport by convection) are negligible for the fuel cell conditions discussed here.

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
  
where the unknowns are, the oxygen mole fraction, ``x_{O_2}``; the water mole fraction, ``x_{H_2O}``; the electrolyte (membrane) and electronic potentials, ``\phi_m`` and ``\phi_s`` respectively; the membrane water content, ``\lambda``; and, the temperature, ``T``. The effective transport parameters ``D^{eff}_{O_2}``, ``D^{eff}_{H_2O}``, ``\sigma^{eff}_{m}``, ``\sigma^{eff}_{s}``, ``D^{eff}_{\lambda}``,  ``D^{eff}_T``, and ``k^{eff}`,  are different in the membrane, GDL and CL and depend non-linearly on the design variables. Due to the solution methodology, all equations need to be solved in all the domains, i.e. GDL, CL and membrane. However, some equations are not necessary in some of the cell domains. This is addressed by making the unnecessary transport parameters zero. 
  
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

and

.. math::
  S_{T} = \left\{
  \begin{array}{cl}
  \sigma^{eff}_s ( \vec{\nabla} \phi_s \cdot \vec{\nabla} \phi_s ) \quad & \text{in GDL and MPL} \\~\vspace{-10 pt}\\ 
  \sigma^{eff}_m ( \vec{\nabla} \phi_m \cdot \vec{\nabla} \phi_m ) \quad & \text{in Membrane} \\~\vspace{-10 pt}\\
  -(\nabla \cdot \vec{i}) \eta + \frac{\nabla \cdot \vec{i}}{2F} ( -T f_{ORR} \Delta \bar{S}_{overall} - \bar{H}_{lv} ) + \sigma^{eff}_m ( \vec{\nabla} \phi_m \cdot \vec{\nabla} \phi_m ) \\ \qquad + \sigma^{eff}_s ( \vec{\nabla} \phi_s \cdot \vec{\nabla} \phi_s )  +k_t \frac{\rho_{dry}}{EW} ( \lambda_{eq} - \lambda ) \bar{H}_{sorption} \quad & \text{in CCL} \\~\vspace{-10 pt}\\
  (\nabla \cdot \vec{i}) \eta + \frac{\nabla \cdot \vec{i}}{2F} ( -T (1-f_{ORR}) \Delta \bar{S}_{overall} ) + \sigma^{eff}_m ( \vec{\nabla} \phi_m \cdot \vec{\nabla} \phi_m ) \\ \qquad + \sigma^{eff}_s ( \vec{\nabla} \phi_s \cdot \vec{\nabla} \phi_s )  +k_t\frac{\rho_{dry}}{EW} ( \lambda_{eq} - \lambda ) \bar{H}_{sorption} \quad & \text{in ACL}
  \end{array}
  \right.

where ``\lambda_{eq}`` is given by the sorption isotherm reported by Hinatsu et al. at the corresponding water vapour activity value in the specific location in the CL.

.. Example 1:

.. include:: polarization_curve/readme.rst

.. Example 2:

.. include:: Bhaiya_EA14_Nonisothermal/readme.rst



References
==========

.. _BhaiyaEA14:
  
  M. Bhaiya, A. Putz and M. Secanell, "Analysis of non-isothermal effects on polymer electrolyte fuel cell electrode assemblies", Electrochimica Acta, 147C:294-309, 2014. DOI: http://doi.org/10.1016/j.electacta.2014.09.051