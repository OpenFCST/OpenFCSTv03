******************************************************
OpenFCST - Open Source Fuel Cell Simulation Toolbox
******************************************************

=================
What is OpenFCST?
=================

The Open Fuel Cell Simulation Toolbox (OpenFCST) is an open-source mathematical modelling 
package for polymer electrolyte fuel cells. FCST builds on top of the open-source 
finite element libraries deal.II, therefore many of its requirements in terms 
of operating systems and such are the same as for deal.II. OpenFCST is distributed 
under the MIT License and has been developed as a modular toolbox from which 
you can develop your own applications. It contains a database of physical 
phenomena equations, fuel cell layers and materials, and kinetics mathematical 
models. In addition, it already contains several applications that allow you 
to simulate different fuel cell components. For example, you can simulate a cathode 
electrode (using either a macrohomogeneous or an ionomer-filled agglomerate model), 
an anode electrode, or a complete membrane electrode assembly. The applications are
already provided in OpenFCST and they have been validated with respect to experimental data 
in the literature as well as numerical results from other models implemented
in commercial packages.

OpenFCST is being developed at the Energy Systems Design Laboratory at the 
University of Alberta in collaboration with the Automotive Fuel Cell Cooperation Corp. 
that, together with the Natural Science and Engineering Research Council of Canada,
has provided the majority of the funding required to develop this code. The goal
of OpenFCST is that research groups in academia and in industry use the current 
toolbox to better understand fuel cells and to develop new physics and material 
databases that can then be integrated in the current library.

================
Getting OpenFCST
================

OpenFCST can be either downloaded from the `OpenFCST website <http://www.openfcst.org>`_
or copied from the developers' `GitHub repository <https://github.com/OpenFCST/>`_.
If you are a user, the easiest way to get the code is via the website. Go to Downloads
section and download a .tar file with the source code. You are then ready to install OpenFCST. 

If you are using Git, please follow instructions in GitHub for checking out a copy.
If you want to modify openFCST, then please follow the steps below.

Creating a new branch
**********************

The development branch of OpenFCST is located in a private `Bitbucket repository <https://bitbucket.org/ESDLab/openfcst/>`_.
In order to access it, please contact `Marc Secanell <mailto:secanell@ualberta.ca>`_.

Committing changes to the development branch is not allowed. A pull request of your own branch
is necessary. Every user can create their own branch of openFCST. The recommended convention
for branch usage is that each user creates their branch for each issue in openFCST they
would like to address. The naming convention is: **username/issue_name**. For example, if the
user "secanell" wants to create an issue to fix a bug on postprocessing, the branch would be 
named ``secanell/postprocessing``.
 
Users can either create a branch in their own machines and then push it to Bitbucket or create
the branch directly on Bitbucket. If the branch is created in Bitbucket, then, in order to
checkout the branch to the appropriate machine, the user needs to issue the following command::

  git fetch && git checkout username/issue_name
  
If the branch is created in the local repository first using ``git checkout -b branch_name``,
then you can commit it to the remote Bitbucket server with ``git push -u origin branch_name``.
The *-u* flag means that from now on your branch branch_name will track the branch in Bitbucket.

Adding, changing, staging, comitting and pushing
************************************************
 
Once the branch is created, users can work on that branch for as long as it is needed. Users
can make changes and commit the changes to their local respository using::

  git add file_name
  git commit -m "message about commit"
 
Please DO NOT use ``git add *`` or ``git add -u`` as you then have little control over what
you are staging to be  committed. Using ``git status``, you can see which files have been
changed so that you can add them as appropriate.
 
To push the changes to Bitbucket, you can use::

  git push origin branch_name

Request for branch to be merged into development
*************************************************

Once you have finished fixing the issue you created the branch for, you need to follow these
three steps:

#. Update your origin information using: ``git remote update`` (this will update all your local information regarding the branches on Bitbucket).
#. Merge your branch with the latest version of development using: ``git merge origin/development``. This is VERY important. The administrators will not accept any pull requests that 
   have not been fast-forwarded to the ``origin/development`` branch.
#. Issue a pull request in Bitbucket.
 
There are three main branches:

* Master branch: Stable version of OpenFCST (no pull requests will be accepted to this branch)
* Development branch: The most up-to-date version of OpenFCST, personal branches should be
started from this branch and all pull requests should be submitted to this branch
* Release branch: Branch containing the latest release of openFCST

Workflow for new development
*****************************

If you want to develop new code, please follow this steps: 

* Clone the repository using ``git clone https://your_username@bitbucket.org/ESDLab/openfcst.git``
* Create a new branch related to the new component/issue you would like to work on using ``git checkout -b name_branch``. Note: The command above will create a branch named "name_branch" and will checkout that branch so you are ready to work.
* Once you are done with the development, create a pull request to merge your branch to the
development branch. Note: Merges to Master will be rejected without review.

A reminder: when developing code, please work in Debug mode and test in both Debug and Release
modes before issuing a pull request.

===================
Installing OpenFCST
===================

Requirements:
*************
 
OpenFCST is developed in a Linux operating system using the GNU GCC compiler. It uses our own
CMake scripts and the contributing libraries' CMake scripts, such as the `deal.II <http://www.dealii.org>`_ script,
to configure and compile the library. It supports at least the following platforms:

#. OpenSUSE 12.3, 13.1, LEAP 42.1, Tumbleweed
#. Ubuntu 14.04, 16.04

The following software needs to be installed on your computer in order for OpenFCST to compile
(make sure to have the development versions of the packages as well):
  
#. CMake
#. GNU Make and C++11 support
#. GCC version 4.7 or later (4.8.1 recommended)
#. BLAS and LAPACK libraries (blas-devel and lapack-devel)
#. OpenMPI compiler
#. GNU gfortran compiler
#. Bison
#. qt4-designer and libqt4 (libqt4-devel if qt4-designer is not available)
#. For generating the documentation: DOxygen and Sphinx
#. Boost; the specific packages are iostreams, serialization, system, thread, filesystem, regex, signals, program_options
#. FLEX (for Dakota)
#. Python Packages: SciPy, NumPy, ipython, Sphinx, evtk, vtk, mayavi
#. libconfig-devel and libconfig++-devel
#. patch
    
OpenFCST comes with all required libraries except the optimization library Dakota from Sandia National
Labs (version 5.4_r2206). You can either download and install it yourself, place tar files in the
appropriate folder (specified below) following OpenFCST 
naming convention (specified below), or allow OpenFCST to download them for you if you have an
Internet connection.
  
  
Configuring and installing OpenFCST
***********************************
  
To help with configuring OpenFCST with CMake, we provide a configuration script **openFCST_install**. 

For a typical installation, go to the `openfcst/` folder, and enter the following:

.. code::

  $./openFCST_install --cores=<number of cores> --install-dir=path_for_installation_directory

  
where the variable **--cores** allows you to compile the program using multiple cores and
**--install-dir** allows you to specify the installation directory where openFCST will be installed.
By default, openFCST will create a Build and Install folder in the same directory as the src folder.
Inside the openfcst/ folder, two new folders will appear:

* Install
* Build  
    
The folder **Install**  contains the installation of the code. It contains a **/bin** folder, where
you will find the executable files for OpenFCST, **fuel_cell-2d.bin** and **fuel_cell-3d.bin** for
2D and 3D simulations, and the GUI file,**fcst_gui**. It also contains the folder **examples**, where
you will find several tutorials on how to run openFCST. The folder **doc** contains the HTML documentation
for developers. The **Build** folder is the folder where all object files needed during compilation are
installed. Users can ignore this folder.

If you are using any of your own pre-installed packages, please consult the src/README for more information
on any necessary changes that need to be made. For more options and information about the installation script, type:

.. code::

  $./openFCST_install --help

===============       
Getting started
===============

See the user guide in src/doc/RefGuide/User_Guide.pdf. To launch a sample cathode simulation, go to the installation
folder **Install** and source the environment script. Then go to **examples/cathode/analysis** and execute the 2D fcst binary:

.. code::

  cd YourInstallDir
  source ./fcst_env.sh
  cd example/cathode/analysis
  fcst2D main.prm
  
This will run a cathode simulation with the simulation data parameters specified in **data.prm**.

============       
Known issues
============

Installation
************

No mpicc found
##############

If, after running installation script, OpenFCST reports that mpicc cannot be found, execute

.. code::

  $mpi-selector-menu
  
then logout and login again. This is a known issue of openmpi package.

"C compiler cannot create executables"
######################################

If you are installing OpenFCST with PETSc and p4est fails to install with an error
"C compiler cannot create executables", perform the same operation with mpi-selector-menu as above.

PETSc error "Could not find a functional BLAS/LAPACK"
#####################################################

If OpenFCST fails to install with PETSc with an error "Could not find a functional BLAS/LAPACK",
install blas-devel and lapack-devel packages.

mpif90 error in OpenSUSE LEAP
#############################

If you are using OpenSUSE LEAP, you might face an error during installation of OpenFCST that says
"...mpif90 is not able to compile a simple test program". In case that happens, install gcc-fortran
package with its dependencies.

Running simulations
*******************

Error "While reading VTK file, unknown file type encountered"
#############################################################

This error may appear when you are trying to run an application that needs to read a 3D .vtk mesh,
but you only have a 2D binary compiled (or vice versa). Compile the code for the required dimension
of the problem using the flag

.. code::

  --openfcst-dimen=X

where X is 1 for both 2D and 3D, 2 for 2D, and 3 for 3D.
  
=======
License
=======

Please see the file src/LICENSE or doc/LICENSE for details.
  
===================
Further information
===================

Visit the `OpenFCST <http://www.openfcst.org/>`_ website.
