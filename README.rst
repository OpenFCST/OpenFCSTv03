******************************************************
OpenFCST - Open Source Fuel Cell Simulation Toolkit
******************************************************

=================
What is OpenFCST?
=================

The Open Fuel Cell Simulation Toolbox (OpenFCST) is an open-source mathematical modelling 
package for polymer electrolyte fuel cells. FCST builds on top of the open-source 
finite element libraries deal.II, therefore many of its requirements in terms 
of operating systems and such are the same as for deal.II. OpenFCST is distributed 
under the MIT License. OpenFCST has been developed as a modular toolbox from which 
you can develop your own applications. It contains a database of physical 
phenomena equations, fuel cell layers and materials, and kinetics mathematical 
models. In addition, it already contains several applications that allow you 
to simulate different fuel cell components. For example, you can simulate a cathode 
electrode (using either a macrohomogeneous or an ionomer-filled agglomerate model), 
an anode electrode or a complete membrane electrode assembly. The applications 
already provided in FCST and they have been validated with respect to experimental data 
in the literature as well as numerical results from another model implemented
in a commercial package.

OpenFCST is being developed at the Energy Systems Design Laboratory at the 
University of Alberta in collaboration with the Automotive Fuel Cell Cooperation Corp. 
that, together with the Natural Science and Engineering Research Council of Canada
has provided the majority of the funding required to developer this code. The goal
of OpenFCST is that research groups in academia and in industry use the current 
toolbox to better understand fuel cells and to develop new physics and material 
databases that can then be integrated in the current library.

================
Getting OpenFCST
================

OpenFCST can be either downloaded from the OpenFCST website, i.e., http://www.openfcst.org or copied from the developer GIT repository.
If you are a user, the easiest way to get the code is via the website. Go to Downloads and download a .tar file with the source code. You are 
then ready to install. 

If you are using Git please follow the step in the BitBucket site for checking out a copy. If you want to modify openFCST, then
please follow the steps below:

Creating a new branch
**********************

Committing changes to the development branch is not allow. A pull request of our own branch is necessary.
Every user can create their own branch of openFCST. The recommended convention for branch usage is that each user creates his/her own branch for each issue in openFCST he/she
would like to address. The naming convention is: **username/issue_name**, for example, if secanell wants to create an issue to fix a bug on postprocessing, the class would be 
named ``secanell/postprocessing``.
 
To create a branch users can either create it in their own machine and then push it to BitBucket or create the branch directly on BitBucket. If the branch is created in 
BitBucket, then in order to checkout the branch to the appropriate machine the user needs to issue the following command::

  git branch branch_name origin/branch_name  
  git checkout branch_name
  
Both steps can be performed simultaneously with:  ``git checkout -b username/issue_name origin/username/issue_name``
 
If the branch is created in the local repository first using ``git checkout -b branch_name``, then you can commit it to BitBucket, i.e. remote server using ``git push -u origin branch_name``.
The *-u* flag means that from now on your branch branch_name will track the branch in BitBucket.

Adding, changing, staging, comitting and pusing
************************************************
 
Once the branch is created users can work on that branch for as long as needed. Users can make changes and commit the changes to their local respository using::

  git add file_name
  git commit -m "message about commit"
 
Please DO NOT use ``git add *`` or ``git add -u`` as you then have little control over what you are staging to be  committed. Using ``git status`` you can see which files have changed so that you can add them
as appropriate.
 
To commit to BitBucket, you can use::

  git push origin branch_name

Request for branch to be merged into development
*************************************************

Once you have finished fixing the issue you created the branch for, you need to follow these three steps:

#. Update your origin information using: ``git remote update`` (this will update all your local information regarding the branches on BitBucket)
#. Merge your branch with the latest version of development using: ``git merge origin/development``. This is VERY important. The administration will not accept any pull requests that 
   have not been fast-forwarded to the ``origin/development branch``.
#. Issue a pull request in BitBucket
 
There are three main branches  
* Master branch: Stable version of openFCST (no pull requests will be accepted to this branch)
* Development branch: The most up-to-date version of openFCST, personal branches should be started from this branch and all pull requests should be submitted to this branch
* Release branch: Branch containing the latest release of openFCST

Workflow for new development
*****************************

If you want to develop new code, please follow this steps: 

* Clone the repository using: git clone https://your_username@bitbucket.org/ESDLab/openfcst.git
* Create a new branch related to the new component/issue you would like to work on using: ``git checkout -b name_branch``.
  Note: The command above will create a branch named name_branch and will checkout that branch so you are ready to work.
* Once you are done with the development, ask for a pull request to merge your brach to the development branch
  Note: Merges to Master will be rejected without review.

A reminder, when developing code, please work on Debug mode (the current version gives an error once the program finishes 
in debug mode, please ignore for now as we will be working on fixing this) and test on Debug and Release mode before 
issuing a pull request. I am aware that running tests in Debug more is time consuming, but the issues that we have in 
debug mode have occurred precisely because we did not test on that mode.

===================
Installing OpenFCST
===================

Requirements:
*************
 
OpenFCST is developed on a Linux operating system using the GNU GCC compiler. It uses our own CMake scripts and the contributing libraries CMake scripts,
such as the deal.II (www.dealii.org) script, to configure and compile the library. It supports at least the following platform:
  1. OpenSUSE 12.3, 13.1, LEAP 42.1
  2. Ubuntu 14.04

The following software needs to also be installed on your computer in order for FCST to compile (make sure to have the development versions as well):
  1. CMake
  2. GNU Make and C++11 support
  3. GCC version 4.7 or later (4.8.1 Recommended)
  4. BLAS and LAPACK libraries 
  5. OpenMPI compiler
  6. GNU gfortran compiler
  7. Bison
  8. qt4-designer and libqt4 (libqt4-devel if qt4-designer is not available)
  9. For generating the documentation: DOxygen and Sphinx
  10. Boost; the specific packages are iostreams, serialization, system, thread, filesystem, regex, signals, program_options
  11. FLEX (For Dakota)
  12. Python Packages: SciPy, NumPy, ipython, Sphinx, evtk, vtk, mayavi
  13. libconfig-devel and libconfig++-devel
  14. patch
    
openFCST comes with all required libraries except:
  1. The optimization library DAKOTA from Sandia National Labs (version 5.4_r2206)
  2. The parallel adaptive mesh refinement library p4est (version 3.4.2)
  3. The partition mesh for parallel computing library Metis (version 5.1)
  4. The parallel linear algebra solver library PETSc (version 3.4.4)
  
You can either download them yourself and install them yourself, place tar files in the appropriate folder (specified below) following OpenFCST 
naming convention (specified below) or allow OpenFCST to download them for you if you have an internet connection.
  
  
Configuring and installing OpenFCST
***********************************
  
To help with configuring OpenFCST with CMake we have provided a configure script, i.e., **openFCST_install**. 

For a typical installation, go to the `openfcst/` folder, and enter the following:

.. code::

  $./openFCST_install --cores=<number of cores> --install-dir=path_for_installation_directory

  
where the variable **--cores** allows you to compile the program using multiple CPUs and **--install-dir** allows you to specify the
installation directory where openFCST will be installed. By default, openFCST will create a Build and Install folder in the same 
directory as the src folder; i.e. Inside the openfcst/ folder, two new folders will appear

    - Install
    - Build  
    
The folder **Install**  contains the installation of the code. It contains a **/bin** folder where you will find the executable files
for OpenFCST, i.e. **fuel_cell-2d.bin** and **fuel_cell-3d.bin** for 2D and 3D simulations, and the GUI file, i.e. **fcst_gui**. It
also contains the folder **examples** where you will find several tutorials on how to run openFCST. The folder **doc** contains
the HTML documentation for developers. 
The **Build** folder is the folder where all object files needed during compilation are installed. Users can ignore this folder.

If you are using any of your own pre-installed packages please consult the src/README for more information on any 
necessary changes that need to be made as is the case for Metis deal.ii and Dakota. For more options and information about the 
installation script type:

.. code::

  $./openFCST_install --help

============       
Known issues
============

No mpicc found
**************

If, after running installation script, OpenFCST reports that mpicc cannot be found, execute

.. code::

  $mpi-selector-menu
  
then logout and login again. This is a known issue of openmpi package.

Boost and deal.ii
*****************

OpenFCST uses a powerful finite element library deal.ii. Current version of OpenFCST is shipped with deal.ii 8.1.1,
which is incompatible with new versions of boost library package that is used in most up-to-date linux distributions such as
openSUSE Tumbleweed. It is recommended to manually install boost 1.53 and use --boost-dir= flag in OpenFCST installation
command as shown in example below.

.. code::

  $./openFCST_install --cores=2 --boost-dir=PATH

In this code, PATH is a full path to Install/ directory of boost.

mpif90 error in OpenSUSE LEAP
**************

If you are using OpenSUSE LEAP, you might face an error during installation of OpenFCST that says "...mpif90 is not able to compile a simple test program". In case that happens, install gcc-fortran package with its dependencies.

===============       
Getting started
===============

See user guide. To launch a sample cathode simulation, go to the install folder **Install** source the environment script. Then go to **examples/cathode/analysis** and execute the 2D fcst binary::

.. code::

  cd YourInstallDir
  source ./fcst_env.sh
  cd example/cathode/analysis
  fcst2D main_app_cathode_analysis.prm
  
This will run a cathode simulation with the simulation data parameters specified in **data_app_cathode_analysis.prm**.

=======
License
=======

Please see the file src/LICENSE or doc/LICENSE for details
  
===================
Further information
===================

Visit the "OpenFCST":http://www.openfcst.org/ website

