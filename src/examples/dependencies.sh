#!/bin/bash 

if (( $EUID != 0 ));
  then 
  echo "#################################################################"
  echo "#	This script requires root privilidges to execute	#"
  echo "#   Please run the script as root again to install dependencies #"
  echo "#################################################################"
else 

zypper in -f vtk-6.0.0-3.1.4

  echo "#################################################################"
if ! $(python -c "import traits" &> /dev/null); then
  echo "#   		 Installing traits			       #"  
  easy_install -q traits 
        
else
  echo "#   	 Python module traits already present		       #"
fi


if ! $(python -c "import matplotlib" &> /dev/null); then
  echo "#   		 Installing matplotlib			       #"  
  easy_install -q matplotlib 
        
else
  echo "#   	 Python module matplotlib already present	       #"
fi

if ! $(python -c "import sphinxcontrib.programoutput" &> /dev/null); then
  echo "#   	Installing sphinxcontrib-programoutput		       #"  
  easy_install -q sphinxcontrib-programoutput
        
else
  echo "#   Python module sphinxcontrib-programoutput already present  #"
fi

if ! $(python -c "import apptools" &> /dev/null); then
  echo "#   		 Installing apptools			       #"  
  easy_install -q apptools 
        
else
  echo "#   	 Python module apptools already present		       #"
fi

if ! $(python -c "from tvtk.api import tvtk" &> /dev/null); then
  
  echo "#   		 Installing Mayavi			       #"        
  git clone https://github.com/enthought/mayavi.git >>install.log
  cd mayavi
  python setup.py install >>install.log
  cd ..
  rm -rf mayavi
  
else
  echo "#   	 Python module tvtk already present		       #"
fi

  echo "#################################################################"

fi
