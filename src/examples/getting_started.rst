================
Getting started
================


Installing OpenFCST
===================

To help with configuring OpenFCST with CMake we have provide a configure script, i.e., **openFCST_install**. 

For a typical installation, go to the `fcst/` folder, and enter the following:

.. code::

  $./openFCST_install --cores=<number of cores> --install-dir=path_for_installation_directory

  
where the variable **--cores** allows you to compile the program using multiple CPUs and **--install-dir** allows you to specify the
installation directory where openFCST will be installed. By default, openFCST will create an in-build installation the openfcst/ folder. 
Inside the openfcst/ folder, two new folders will appear

    - Install
    - Build  
    
The folder **Install**  contains the installation of the code. It contains a **/bin** folder where you will find the executable files
for OpenFCST, i.e. **fuel_cell-2d.bin** and **fuel_cell-3d.bin** for 2D and 3D simulations, and the GUI file, i.e. **fcst_gui**. It
also contains the folder **examples** where you will find several tutorials on how to run openFCST. The folder **doc** contains
the HTML documentation for developers. 
The **Build** folder is the folder where all object files needed during compilation are installed. Users can ignore this folder.

For more options and information about the installation script type:

.. code::

  $./openFCST_install --help
  
System requirements
-------------------

The following software needs to also be installed in your computer in order for FCST to compile:
  
  1. GNU make and C++11 support, gcc version 4.7 or later (4.7.2 Recommended)
  2. GCC
  3. BLAS and LAPACK libraries 
  4. OpenMPI compiler
  5. gfortran compiler
  6. Bison
  7. qt4-designer 
  8. For generating the documentation: DOxygen
  9. Boost; the specific packages are iostreams, serialization, system, thread, filesystem, regex, signals, programoptions  
  10. FLEX (For Dakota)
  11. Python Packages: SciPy, NumPy, ipython, Sphinx, evtk, vtk, mayavi


  
Setting up a CONDA environment for Python packages
===================================================

The PythonFCST package which can be used for generating VTK meshes for microstructure simulations and a bunch of plotting routines use additional python libraries. The python package uses Python 2.7 and therefore is not compatible with Python 3. The migration to Python 3 will be done in a later release. Therefore, it is suggested that the users create a CONDA environment to use the PythonFCST package. 

In order to set up a CONDA environment with the necessary OpenFCST dependencies, please follow these steps:

1. Download the Python 3.5 64 Bit CONDA installer from the website: https://www.continuum.io/downloads

2. To install CONDA open a Konsole/Terminal and go to the folder where the installer was downloaded and type,

.. code:: 

    bash Anaconda3-4.2.0-Linux-x86_64.sh

3. When it asks you for a path to install anaconda, preferably put some local folder. For this walkthrough, it is assumed that you install it in /home/user/anaconda3 . Also it is assumed that you prepend your .bashrc with the anaconda bin folder that is: 

.. code::

    export PATH="/home/secanell/anaconda3/bin:$PATH"

4. To create the Python environment with the necessary pacakges required for OpenFCST type,

.. code:: 

    conda create -n fcst_python python=2.7 anaconda mayavi
    
This will install a general array of scientific python packages. To install the bare-minimum packages required for PythonFCST use,

.. code::

    conda create -n fcst_python python=2.7 scipy numpy mayavi matplotlib pandas pil

5. Next make sure that QT_API is set properly. 

.. code::

    export QT_API=pyqt


6. To use the environment, open Konsole and type:

.. code::

    source activate fcst_python


This should enable the use of all Python library functions within OpenFCST. 

.. note::
    
    export QT_API=pyqt might be required everytime. Therefore, it should be added to the .bashrc. OpenFCST environment script, i.e., $FCST_DIR/Install/fcst_env.sh sets the QT_API. Therefore, sourcing the OpenFCST environment script is sufficient and no additonal changes will be required.

  
Setting up a simulation
=======================

The easiest way to start a simulation with OpenFCST is running the examples provided in the `examples` folder. Each example is a chapter in this user guide. Once
you are comfortable running the examples, then you can use the GUI to modify different parameters to perform new simulations.

We recommend that you use the OpenFCST graphical user interface. To load the examples in the GUI, first you will need to create the `.xml` files that the GUI
reads and then you will have to open them with the GUI. For an example see the Introduction to AppCathode tutorial