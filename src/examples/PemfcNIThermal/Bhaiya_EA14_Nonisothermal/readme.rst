
Example 2: M. Bhaiya et al., Analysis of non-isothermal effects on polymer electrolyte fuel cell electrode assemblies
=====================================================================================================================

*****************
Abstract
*****************

A non-isothermal, single phase membrane electrode assembly (MEA) mathematical model accounting for most applicable heat sources, viz., reversible, irreversible, ohmic heating, 
phase change, heat of sorption/desorption, is presented. The mathematical model fully couples a thermal transport equation with an MEA model and allows the study of non-isothermal 
effects, such as thermal osmosis through the membrane, local relative humidity variations in the catalyst layers and water sorption into the membrane. A detailed breakdown of 
various heat sources in the MEA at different current densities is provided and the impact of various thermal effects previously neglected in the literature such as thermal-osmosis, 
reversible heat distribution, and heat of sorption are studied. Results show that sorption heat cannot be neglected as it contributes up to 10\% of the total heat under normal 
operating conditions. Reversible heat distribution can significantly affect the temperature distribution shifting the hottest location of the cell from anode and cathode. 
Analyzing the water transport across the membrane, results show that thermal-osmosis contributes up to 25\% of the water flux inside the membrane at moderate and high current densities. 
The microporous layer (MPL) thermal conductivity is found to have a significant effect on fuel cell performance with severe dry-out observed at low MPL thermal conductivities.


*******************
Directory structure
*******************

This folder contains the files used in order to reproduce the results reported in reference [1].

The data is subdivided as follows:

1. Directory ``template/`` contains the main main.prm, data.prm and opt.prm file used in the simulations. In each of the sub-folders these files are extended. This template directory contains separate 
files for isothermal and non-isothermal simulations. Both the simulations use different applications, namely, Pemfc, and PemfcNIThermal. However, geometry configuration, transport properties, operating 
conditions, kinetics properties etc, are all maintained same. Base case in the template file represents, cell temperature (for the isothermal case), or bipolar plates (for the non-isothermal case) at 80C, anode/cathode pressure at 1.0 atm, and anode/cathode relative humidity at 50 %.

2. Directory ``section_5/`` contains the files to run various cases discussed in the Section 5 (Results and Discussion) of the reference [1]. At base level of this directory, there are various 
named folders, which extend the template files to simulate various operating conditions. These files are used in describing results in various subsections of the section 5 in the article. The names of 
these folders are self-explanatory. For instance, ``nonisothermal_RH_50_P_2atm``, simulates a non-isothermal MEA model at 50 % anode/cathode RH, and anode/cathode pressure at 2.0 atm.

    2(i). Sub-directory ``section_5_4`` contains the additional files required for Section 5.4 of the article. This subsection primarily deals with the effects of considering heat of sorption in the 
    non-isothermal model predictions. So the additional files in here, turn heat of sorption effects OFF, besides varying operating conditions.
    
    2(ii). Sub-directory ``section_5_5`` contains the additional files required for Section 5.5 of the article. This subsection primarily deals with the effects of reversible heat distribution in the 
    non-isothermal model predictions. So the additional files in here, vary the fraction of reversible heat in the ORR, at base case (default value of the fraction is 1.0).
    
    2(iii). Sub-directory ``section_5_6`` contains the additional files required for Section 5.6 of the article. This subsection deals with the effects of considering thermal osmosis in the 
    non-isothermal model predictions and water management of the cell. So the additional files in here, turn thermal osmosis OFF, besides varying operating conditions.

*******************
Template files
*******************


