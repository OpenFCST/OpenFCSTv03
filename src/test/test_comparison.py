#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
mkPolcurve.py
Displays polarisation curve created by fcst in the parametric mode

Input:
    - Datafile, e.g. "input.dat"
      Format:
      %eval_id         V_cell         obj_fn 
       1            0.3   -2.380636381 
       2           0.35   -2.223311025 
       3            0.4   -2.015362082 
       4           0.45   -1.749832606 
       5            0.5   -1.430078875 
       6           0.55   -1.078322382 
       7            0.6  -0.7347111835 
       8           0.65  -0.4420970865 
       9            0.7  -0.2278052956 
      10           0.75 -0.09479962371 
      11            0.8 -0.02921816945 
"""
#! Polarization Curve
#! ###############################################

# Basic imports:
import sys, os, shutil
import glob
import re
import commands

from StringIO import StringIO # Create stringlike file objects

# Import scientific packages:
#from pylab import *
from math import sqrt
import scipy
import scipy.io
import scipy.interpolate
import scipy.optimize

print "="*50
print "="*50
print "= Checking the accuracy of the test results against the expected results"
print "-"*50

print "-"*50
print "= - Parse Commandline "
# Import command line option parser:
from optparse import OptionParser
# Setup of the command line options:
usage = "usage: %prog [options] filename.dat"
parser = OptionParser(usage)
options, filename = parser.parse_args(sys.argv[1:])
print "-"*50

print "-"*50

print "= - Load data"
print "Checking test results from file:", filename[0]
print "by comparing with data in file:", filename[1]
myTestData = scipy.loadtxt(filename[0],skiprows=1);
myExpectedResults = scipy.loadtxt(filename[1],skiprows=1);
print "-"*50

print "-"*50
print "Error between the two files:"


"""diff=(myTestData[:,2]-myExpectedResults[:,2])**2;
diff2 = []
for d in diff:
	diff2.append(sqrt(d))

diff = diff2
diff=diff[:]/-myExpectedResults[:,2];
print diff"""


if len(myTestData) != len(myExpectedResults):
  print "Lengths of data do not match"
  exit(1) 
  
diff = []  

for i in range(len(myTestData)):
    diff.append(abs((myTestData[i,2]-myExpectedResults[i,2])/myExpectedResults[i,2]))



print diff
print "-"*50

print "-"*50
print "Largest error between between the results:"
max_error = max(diff)
print max_error

if max_error == 0:
  print "Test results match expected results"
  print "="*50
  print "="*50
  exit(0)
elif max_error < 0.01:
  print "Test results are within numerical error (1%), with the greatest being: ", max_error
  print "="*50
  print "="*50
  exit(0)
elif max_error < 0.1:
  print "Test results have an unacceptable error (between 1% and 10%), with the largest being: ", max_error
  print "="*50
  print "="*50
  exit(1)   
else:
  print "Test results differ significantly from the expected results (greater than 10%): ", max_error
  print "="*50
  print "="*50
  exit(1)  


