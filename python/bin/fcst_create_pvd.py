#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: FCST Development Team
# License: TBD
# Copyright (c) 2014, MIT License
"""
Python script to create a pvd paraview file
"""

# Basic imports:

from PythonFCST.util.parsers import CreatePvd
import argparse

print "="*50
print "= Paraview pvd file creator"
print "="*50

print "-"*50
print "= - Parse Commandline "
# Import command line option parser:
parser = argparse.ArgumentParser(
    description=
	'Program to assemble a sequence of vtk files into a .pvd file. ' \
    'Note: this does not work with the legacy format, change to vtu!!!' \
    '\n\n\n'\
    'Usage: fcst_create_pvd.py --help'
	,formatter_class=argparse.RawDescriptionHelpFormatter)
	
parser.add_argument('--pattern', type=str,
       default='fuel\\_cell\\_solution\\_DataFile_Cycle\\_[0-9]\\_Sol\\_(?P<order>[0-9]+)\\.vtu',
	help='regex pattern for the filename (default is set)')
  
parser.add_argument('--output', type=str,
                    default="paraview.pvd",
                    help='Name for the paraview .pvd file')

args = parser.parse_args()

print "-"*50
print "= - Parse current directory"
pvdparser=CreatePvd(pattern=args.pattern,fname=args.output,group='order')

pvdparser.run()
del pvdparser
    

