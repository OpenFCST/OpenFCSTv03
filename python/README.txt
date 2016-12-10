========================================================
PythonFCST - Python interface for fcst
========================================================

The fcst python module contains usefuel pre and processors
to fcst. A SWIG interface to the source code is to follow.


Usage Examples
==============

import PythonFCST

Mesh Generation
----------------
The mesh module in PythonFCST can be used to generate a VTK mesh from a 3D numpy array. In order generate a VTK file simply use the writeVTK.py file in the python/bin folder.
The following commands can be used:
>./writeVTK.py -f stack -b Slice_ -n 200 -o test.vtk
The above command will go into the "stack" folder relative to the current working directory and read images with the name "Slice_001.tiff" upto "Slice_200.tiff". The output mesh file would have a name test.vtk
Except for the -f flag the other flags have the same default value as above. For a complete list of options just use the command:
>./writeVTK.py -h
