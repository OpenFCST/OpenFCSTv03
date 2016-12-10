#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: FCST Development Team
# License: TBD
# Copyright (c) 2014, MIT License
"""
Python script to create a pvd paraview file
"""

# Basic imports:

from PythonFCST.util.parsers import read_polcurve
import argparse
from itertools import cycle
import pylab as pl
import os.path

print "="*50
print "= Creates postprocessing plots from the tabular output file of the fcst polcurve sweep"
print "="*50

print "-"*50
print "= - Parse Commandline "

# Import command line option parser:
parser = argparse.ArgumentParser(
    description="""
Python executable to create postprocessing plots from
a polarization_curve.dat results file
Format:
    
      # ====================================================
      # OpenFCST: Fuel cell simulation toolbox 
      # ====================================================
      # Polarization curve data :polarization_curve.dat
      Cell voltage [V]        Cathode current [A/cm2]
      0.9     0.007674716442 
      0.8     0.09141826831 
      0.7     0.3239937606
      0.6     0.7566928344 
      0.5     1.458893632
      0.4     2.18243
      
Usage: fcst_plot_dakota --help\n'
Load a default polcurve: fcst_plot_polcurve.py polcurve.dat polcurve_reference.dat --label ['Current', 'Reference']
    """
	,formatter_class=argparse.RawDescriptionHelpFormatter)

	
parser.add_argument('files', type=str,nargs='+',
		help='Names of one ore more polcurve files files')
                      
parser.add_argument('--label', type=str, nargs='+',default=['Testcase'],
                    help='space separated list of labels (size has to match the number of processed files)')

parser.add_argument('--output', type=str,
                    default="",
                    help='Basename for output files')


args = parser.parse_args()
sweep_files = args.files

figures={
            "polcurve":     None,
        }

plot_colors = cycle(['k','b','g','r','c','m','y'])
plot_labels = cycle(args.label)

for myfile in sweep_files:
    print "\n"*3+"-"*50
    print "= - Parse file ", myfile, '\n'*1    
    
    color = plot_colors.next()
    label = plot_labels.next()
   
    if os.path.isfile(myfile):
        parser = read_polcurve(fname=myfile)
        figures["polcurve"] = parser.plot_polcurve(fig=figures["polcurve"], label=label,
                                            format_pol   = color+'-', 
                                            format_power = color+'--',
                                            flag_guides=False,flag_power=False)
                                                             
    else:
        print "  File ", myfile, 'does not exist'

pl.show()

if args.output:
    for key in figures.keys():
        if hasattr(figures[key],'savefig'):
            print 'Save figure ', key
            figures[key].savefig(   args.output+'_' + key + '.png',
                                    bbox_inches='tight',
                                    #bbox_extra_artists=(lgd,)
                                )
    
