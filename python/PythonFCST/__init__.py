# Author: Andreas Putz
# Copyright (c) 2013, PyhtonFCST
# License: TBD.
r"""
:mod:`PythonFCST`: Python helpers for OpenFCST
===============================================================================

Documentation is available in the docstrings and in ths sphinx documentation.

Contents
--------
The PythonFCST package imports all the functions from the top level modules.

Subpackages
-----------

.. list-table:: `PythonFCST` module structure.
   :widths: 10 80 10
   :header-rows: 1

   * - name
     - description
     - autoload
   * - :mod:`PythonFCST.util`
     - common utilities and classes used by most of the othe rmodules
     - yes
   * - :mod:`PythonFCST.mesh`
     - Mesh generation classes
     - yes
   * - `VISU`
     - Visualisation methods
     - yes
   * - `interactive/`
     - setup of IPython-based shell `isfepy`
     -
 
Prefix convention
-----------------
::

  None as yet
 
 
Utility tools
-------------
::

 TODO                --- Todo
 
 
Import
------
>>> import PythonFCST as fcst

Inheritance diagram
-------------------
    
.. inheritance-diagram:: PythonFCST

Package Documentation
---------------------

.. automodule:: util
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: mesh
   :members:
   :undoc-members:
   :show-inheritance:

"""

__version__ = '0.0.1'

__requires__ = [
    'scipy',
    'numpy',
]

#__extras_require__ = {PyFCell
#    'app': [
#        'envisage',
#    ],
#}

# __all__ = ['misc']






import util
import mesh

import pylab as pl
def pylab_setup():
    fig_width_pt        = 246.0         # Get this from LaTeX using \showthe\columnwidth
    inches_per_pt   = 1.0/72.27                 # Convert pt to inch
    fig_width_cm        = 30.0
    inches_per_cm   = 0.393
    golden_mean = (pl.sqrt(5)-1.0)/2.0             # Aesthetic ratio
    fig_width = fig_width_cm*inches_per_cm      # width in inches
    fig_height = fig_width*golden_mean          # height in inches
    fig_size =  [fig_width,fig_height]
    params = {
                'backend': 'ps',
    #           'backend': 'svg',
                'axes.labelsize': 24,
                'axes.titlesize': 28,
#                'text.fontsize': 20,                
                'legend.fontsize': 20,
#                'ticks.font':'Helvetica',
                'xtick.labelsize': 20,
                'ytick.labelsize': 20,
                'text.usetex': True,
                'text.latex.preamble': [
                    r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
                    #r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
                    r'\usepackage{helvet}',    # set the normal font here
                    r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
                    r'\sansmath'],  # <- tricky! -- gotta actually tell tex to use!,
                'figure.figsize': fig_size,
                'font.family':'sans-serif',
                'font.serif':'Computer Modern Roman',
                'font.sans-serif':'Helvetica'
            }
    for key in params.keys():
        if pl.rcParams.has_key(key):
            pl.rcParams[key] = params[key]
        else:
            print "Your version of matplotlib does not support the parameter ", key

pylab_setup()