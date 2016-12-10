#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 16:56:58 2014

@author: msabharwal
"""

from PIL import Image
import os
import numpy as np
import time
from optparse import OptionParser

#==============================================================================
#Initializing the different options for the application
#any new options should be added in this section
#==============================================================================
parser=OptionParser()
parser.add_option("-f","--file",dest="foldername",action="store",help="relative path of the images to be converted",metavar="FILE",type="string")
parser.add_option("-b", "--basename",action="store",dest="basename",help="basename of the image files. Default: 'Slice_' would mean images with name Slice_001 and so on",default="Slice_",type="string")
parser.add_option("-n", action="store",dest="num",help="Number of images to be read. Default: 200",default=200,type="int")
parser.add_option("-o", action="store",dest="output",help="Output file name. Default: test.vtk",default='test.vtk',type="string")
parser.add_option("-m","--material",dest="material",action="store",help="material id of the phase to be meshed. Default: 0",default=0,type="int")
parser.add_option("-v","--voxel",dest="voxel",action="store",help="voxel size as x y z",metavar="VOXEL",type="float",nargs=3)
parser.add_option("-e","--extension",dest="extension",action="store",help="extension for the input files. Default is .tiff",default=".tiff",type="string")
(options,args)=parser.parse_args()

if options.foldername == None:
    parser.error("Please enter the image location using the option -f or --file; For more help options check -h or --help")

foldername = options.foldername
basename = options.basename
num_images = options.num
output = options.output
start=time.time()
material = options.material
extension=options.extension
voxel=options.voxel
if not voxel:
    voxel=[1e-6,1e-6,1e-6]    
voxel=[voxel[1],voxel[0],voxel[2]]

#==============================================================================
#Reading Images from the data folder and constructing a 3D stack
#==============================================================================
name=os.getcwd()+'/'+foldername+'/'+basename+'%.3d'

tmp=[]
for i in range(num_images):
    filename=name%(i+1)+extension
    tmp.append(np.asarray(Image.open(filename),dtype=np.float))
temp=np.dstack(tmp)
image=np.ones(temp.shape)*20
image[temp==material]=0
#==============================================================================
#Generating the VTK mesh using the Mesh module in PythonFCST
#==============================================================================
import PythonFCST as fcst
g = fcst.mesh.PhaseGenerator
o=g(image,output,scale=voxel)
o.write()
print "Phase with material id:",material," is represented by material id 0 in the mesh file: ",output
print "time elapsed is ", time.time()-start, " seconds"