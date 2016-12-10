
Example 1: A simple polarization curve 
======================================

This examples demonstrates the calculation of a polarization curve based on the
non-isothermal MEA model. This calculation is based on a default selection of
material properties described below. A regression test has been performed with three
adaptive refinement steps and hundred points on the polarization curve down to an
cell voltage of 0.01 Volts. The results of the run scripts are plotted against these
regression results.


**************************************************************
Simulation Setup
**************************************************************

The Data Files
###############

The *main.prm*, *data.prm* and *opt.prm* files are identical to the global template 
files, except for small modifications through the lexicographic replacer
sed executed by the run.sh file.



*****************************************************************
Simulation Execution
*****************************************************************

Main Simulation execution
#########################


Postprocessors
###############

The simulation produces a number of result files. The most important are:

 #. **logfile.log** - The main logfile, useful for debugging purposes
 #. **data_polcurve.dat** - The result file from the optimization. This can be used to generate polcruves, water crossover curves, voltage loss breakdowns and many more.
 #. **.vtu** files - VTK files which can be read with any recent viewer supporting vtk (e.g. Paraview, ViSIT, Mayavi2, ...). OpenFCST ships with a posprocessor to generate a correctly ordered collection file (**.pvd**)


*****************************************************************
Results
*****************************************************************


The Polarization Curve
#######################

.. plot::

  
   import PythonFCST as PythonFCST
   from PythonFCST.util.parsers import read_polcurve
   import argparse
   from itertools import cycle
   import pylab as plt
   import os
   flag=0
   path=os.getcwd()
   for subdir,dirs,files in os.walk(path):
     for file in files:
         name=os.path.join(subdir,file)
         if 'polarization_curve/regression/test_results.dat' in name:
             fname=name
   
   data01 = read_polcurve(fname)
   fig = data01.plot_polcurve(label='Reference Results',
                            format_pol   = 'r-', 
                            format_power = 'r--',
                            flag_guides=False)
   fname=[]
   for subdir,dirs,files in os.walk(path):
     for file in files:
         name=os.path.join(subdir,file)
         if 'polarization_curve.dat' in name:
             fname=name
             flag=1

   if flag:
      data02 = read_polcurve(fname)
      fig = data01.plot_polcurve(label='Current Results', fig = fig,
                            format_pol   = 'g-', 
                            format_power = 'g--',
                            flag_guides=False)   
   plt.show()

