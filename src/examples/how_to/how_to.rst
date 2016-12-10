================================
 How to create your own tutorial 
================================

Introduction
============

This tutorial describes the main structure of each tutorial for OpenFCST. Each tutorial is associated with the different applications in
OpenFCST and it is intented to provide an example for users to learn now to use the application. 

The main structure of a tuturial should contain:

1. Governing equations: This section explains the equations that are solved in the application that is being showcased.

2. Directory structure: This section explains the different folders that are in the tutorial sub-folder. Each subfolder
can be used to explain a different feature available for the given application.

3. Setting up a simulation

4. Examples: This section can contain several sub-sections

Governing equations
===================

The governing equations are (...).

More information can be found in Secanell07a_.


Directory structure
===================

Each directory consists of at least following folders:

1. template : This folder contains the default files for running all the examples in the other folders. Please **do not** modify this file as 
it will result in all tests failing. If you would like to create your own example either include this file to your simulation using the *include*
command or copy the file to a different location. 

2. analysis : This folder is used to run an analysis simulation using the template file.

3. other example folders : These folders shold contain the :code:`main_test.prm` and :code:`data_test.prm` files needed to run a simulation to reproduce the example
results. Note the data file includes the template find and adds the necessary modifications.
The script to run a test to make sure the polarization curve is running correctly is in the folder *regression* together with the default data the test is compared to.
Only the examples that are used in regression tests would have a regression subfolder.


Setting up a simulation
=======================

In this section, the :code:`main.prm` and :code:`data.prm` files are described for the given application. 

a. The :code:`main.prm` / :code:`main.xml` file: This file is used to select the appropriate: a) type of analysis, i.e. analysis, parametric study, polarization curve and optimization study; application; 
b) the nonlinear solver; c) data file name; and, d) several less critical parameters.
b. The :code:`data.prm` / :code:`data.xml` file: This file is used to input all the input data used for the simulation for the application selected.

Both these files can either be loaded and modified via the openFCST graphical user interface (GUI) or modified as a text file. 

Setting up a simulation using the OpenFCST graphical user interface (GUI)
-------------------------------------------------------------------------

Explain how to setup a simulation using the GUI, in particular explain the main sections that need to be modified in each file.


Setting up a simulation using a text (.prm) file
------------------------------------------------

Explain how to setup the data file (main differences from previous applications).

Results
-------

Explain the results obtained usig the default :code:`main.prm` and :code:`data.prm` files to run an analysis case.


Examples
========
   
.. add any .rst subfolder here


References
==========

.. _Secanell07a:
  
M. Secanell, B. Carnes, A. Suleman and N. Djilali, "Numerical Optimization of Proton Exchange Membrane Fuel Cell Cathode Electrodes", Electrochimica Acta, 52(7):2668-2682, 2007.