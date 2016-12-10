#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 16:56:58 2014

@author: Mayank Sabharwal, 2015
"""
PARAVIEW_PATH='@PARAVIEW_DIR@/Build'


paraview_path=PARAVIEW_PATH
import os,sys
import scipy as sp
import numpy as np

if os.path.exists(paraview_path):
    for x in os.walk(paraview_path):
        sys.path.append(x[0])
        os.environ["LD_LIBRARY_PATH"]=os.environ.get("LD_LIBRARY_PATH")+":"+x[0]
    try:
        import paraview.simple as s
        flag=1
    except:
        flag=0
        print "Failed to import Paraview Python Libraries"
        print "Exiting code"
        exit(3)
else: 
    print "Paraview source build not found!"
    print "Set the Paraview Path in openfcst/src/CMakeLists.txt!"
    flag=0
    exit(3)



print "="*50
print "="*50
print "= Checking the accuracy of the test results against the expected results"
print "-"*50

print "-"*50
print "= - Parse Commandline "
# Import command line option parser:
from optparse import OptionParser

# Setup of the command line options:
usage = "usage: %prog [options] fuel_cell_solution_DataFile_00001_Cycle_4.vtk reference_data.dat"
parser = OptionParser(usage)
options, filename = parser.parse_args(sys.argv[1:])
print "-"*50

print "-"*50

print "= - Load data"
print "Checking test results from file:", filename[1]
print "by comparing with data in the simulation result file:", filename[0]
tmp=sp.loadtxt(filename[1],dtype='string')
header = tmp[0]

refData = np.array(tmp[1:],dtype='float')

x=np.array(refData[:,0])
y=np.array(refData[:,1])
z=np.array(refData[:,2])

refResponses={}
for i in range(np.size(header)-3):
    refResponses[header[i+3]] = refData[:,i+3]
simResponses={}
for name in refResponses.iterkeys():
    simResponses[name]=[]

fname=os.getcwd()+'/'+filename[0]

extension = os.path.splitext(filename[0])[1]    
if extension == '.vtk':
    solution = s.LegacyVTKReader( guiName="solution", FileNames=[fname] )
elif extension == '.vtu':
    solution = s.XMLUnstructuredGridReader( guiName="solution", FileName=[fname] )
else:
    print "= - Unknown file format of type: ", extension
    
for i in range(np.size(x)):
    temp=[]
    ProbeLocation1=[]
    ProbeLocation1 = s.ProbeLocation( guiName="ProbeLocation1", ProbeType="Fixed Radius Point Source", Input = solution )
    ProbeLocation1.ProbeType.Center = [x[i],y[i],z[i]]
    temp=s.servermanager.Fetch(ProbeLocation1)
    for name in refResponses.iterkeys():
        if name == 'velocity_X':
            simResponses[name].append(temp.GetPointData().GetArray('velocity').GetValue(0))
        elif name == 'velocity_Y':
            simResponses[name].append(temp.GetPointData().GetArray('velocity').GetValue(1))
        elif name == 'velocity_Z':            
            simResponses[name].append(temp.GetPointData().GetArray('velocity').GetValue(2))
        else:
            simResponses[name].append(temp.GetPointData().GetArray(name).GetValue(0))
error={}
for name in simResponses.iterkeys():
    error[name]=np.absolute(refResponses[name]-simResponses[name])/refResponses[name]
max_error=np.max(error.values())

if max_error == 0:
  print "Test results match expected results"
  print "="*50
  print "="*50
  exit(0)
elif max_error < 0.01:
  print "Test results are within numerical error (1%), with the greatest being: ", max_error*100
  print "="*50
  print "="*50
  exit(0)
elif max_error < 0.1:
  print "Test results have an unacceptable error (between 1% and 10%), with the largest being: ", max_error*100
  print "="*50
  print "="*50
  exit(1)   
else:
  print "Test results differ significantly from the expected results (greater than 10%): ", max_error*100
  print "="*50
  print "="*50
  exit(1)  