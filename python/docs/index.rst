.. FCSTpython documentation master file, created by
   sphinx-quickstart on Thu Nov 29 11:31:33 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
****************************************
PythonFCST: High level interface to fcst 
****************************************


PythonFCST is a suite of python classes for high level pre and post processing for fcst as well as 
a wrapper for the entire fcst library (future!!!)





Documentation
=============

Contents:

.. toctree::
   :maxdepth: 1

   getting_started.rst
   developer_guide.rst
   api.rst
   docPythonFCST.rst
   
.. todolist::

Documentation is available in the docstrings and in ths sphinx documentation.

Contents
--------
The PythonFCST package imports all the functions from the top level modules.

Subpackages
-----------
::

 util	--> 	Contains the base class and little helpers and everything 
 		else feeling homeless otherwise.
 mesh	--> 	Mesh generators (using Salome)
 VISU	--> 	Visualisation methods
 
 
Import
------
>>> import PythonFCST as fcst


Inheritance diagram
-------------------
    
.. inheritance-diagram:: PythonFCST

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

