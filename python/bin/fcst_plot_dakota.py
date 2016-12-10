#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: FCST Development Team
# License: TBD
# Copyright (c) 2014, MIT License
"""
Python script to create a pvd paraview file
"""

# Basic imports:

from PythonFCST.util.parsers import ReadOptPrm, read_dakota
import argparse
from itertools import cycle
import pylab as pl
import os.path

print "="*50
print "= Creates postprocessing plots from one ore more dakota tabular files"
print "="*50

print "-"*50
print "= - Parse Commandline "

# Import command line option parser:
parser = argparse.ArgumentParser(
    description=
	'Python executable to create postprocessing plots from ' \
    'a dakota_tabular.dat results file' \
    'Note: needs to parse a opt.prm file to find the correct column names'
    '\n\n\n'\
    'Usage: fcst_plot_dakota --help\n' \
    'Load a default polcurve: fcst_plot_dakota dakota_tabular.dat dakota_tabular_blessed.dat --opt opt_app_pemfc_parametric.prm'
	,formatter_class=argparse.RawDescriptionHelpFormatter)

	
parser.add_argument('files', type=str,nargs='+',
		help='Names of one ore more dakota tabular files')
  
parser.add_argument('--opt', type=str,
                    default="opt.prm",
                    help='Filename of the fcst optimization input file')
                    
parser.add_argument('--label', type=str, nargs='+',default=['Testcase'],
                    help='space separated list of labels (size has to match the number of processed files)')

parser.add_argument('--output', type=str,
                    default="",
                    help='Basename for outout files')


parser.add_argument('--cell_voltage', type=str,
                    default="Fuel cell data:Operating conditions:Voltage cell",
                    help='Name of the cell voltage column')

parser.add_argument('--current_cathode', type=str,
                    default="cathode_current",
                    help='Name of the total cathodic current (default="current_cathode")')

parser.add_argument('--current_anode', type=str,
                    default="anode_current",
                    help='Name of the total cathodic current (default="current_cathode")')
                    
parser.add_argument('--water_cathode', type=str,
                    default="water_cathode",
                    help='water source term in the cathode')
                    
parser.add_argument('--water_anode', type=str,
                    default="water_anode",
                    help='water source term in the anode')
                    
parser.add_argument('--max_temp', type=str,
                    default="max_temperature",
                    help='maximum temperature')

args = parser.parse_args()
dakota_files = args.files


prm=ReadOptPrm(args.opt)
design_variables = prm.find_design_variables()
objectives       = prm.find_objective_functionals()

figures={
            "polcurve":     None,
            "crossover":    None,
            "max_temp":     None,
            "thermal_losses":   None
        }

plot_colors = cycle(['k','b','g','r','c','m','y'])
plot_labels = cycle(args.label)

for myfile in dakota_files:
    print "\n"*3+"-"*50
    print "= - Parse file ", myfile, '\n'*1    
    
    color = plot_colors.next()
    label = plot_labels.next()
    
    
    if os.path.isfile(myfile):
        dakota = read_dakota(fname=myfile, parameters=design_variables, responses=objectives)
        
        dakota.set_names_polcurve(cell_voltage=args.cell_voltage,
                                  current_cathode=args.current_cathode,
                                  current_anode=args.current_anode
                                 )
        dakota.set_names_crossover( water_cathode=args.water_cathode,
                                    water_anode=args.water_anode
                                  )
        dakota.set_names_thermal( max_temp=args.max_temp
                                )                          
        
        figures["polcurve"] = dakota.plot_polcurve(fig=figures["polcurve"], label=label,
                                            format_pol   = color+'-', 
                                            format_power = color+'--',
                                            flag_guides=False)
        figures["crossover"]= dakota.plot_water_crossover(fig=figures["crossover"],
                                                          label=label,
                                                          format_cathode=color+'-',
                                                          format_anode  =color+'.')
        figures["max_temp"]= dakota.plot_max_temp(fig=figures["max_temp"],
                                                          label=label,
                                                          format_plot=color+'-')
        #figures["thermal_losses"]= dakota.plot_max_temp(fig=figures["thermal_losses"],
        #                                                  label=label,
        #                                                  format_plot=color+'-')                                                 
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
    
