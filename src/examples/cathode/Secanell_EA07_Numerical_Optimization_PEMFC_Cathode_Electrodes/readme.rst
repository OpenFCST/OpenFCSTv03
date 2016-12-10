Example: Optimization of a macro-homogeneous cathode electrode (Ref: Secanell, Electrochimica Acta, 2007)
=========================================================================================================

This folder contains the files used in order to reproduce the results reported in reference [1].

The data is subdivided as follows:

1. Directory ``template/`` contains the main main.prm, data.prm and opt.prm file used in the simulations. In each of the sub-folders these files are extended

2. Directory ``analysis/`` contains the data to run a single point in the polarization curve. It uses the data in ``template/`` directly without modifications.

3. Directory ``parameteric/`` contains the data files used to compute a full polarization curve. To achieve this goal, only the main.prm is modified with respect to the template files by including information about the number of points needed in the polarization curve. The directory also includes an .ods file used to compare the data. The results, as of Feb. 2, 2014, OpenFCST Rev: 1817 match the article.

4. Directory ``design/`` contains two subfolders. Each folder contains the necessary files to setup a simulation in order to obtain the results in sections 4.1 and 4.2 in the artice. Each subfolder references the template data files, so only the main modifications to the original file are highlighted.
Each subdirectory also includes an .ods file used to compare the data. The results, as of Feb. 2, 2014, OpenFCST Rev: 1817 match the article. 

**Warning:** This directory is not tested with regression tests, therefore it is possible that the parameter file has to be modified slightly in order to obtain the correct solution. The folder was last validated on Feb. 3, 2015
prior to openFCST release 0.2.

References
----------
M. Secanell , B. Carnes , A. Suleman , N. Djilali, Numerical optimization of proton exchange membrane fuel cell cathodes, Electrochimica Acta 52 (2007) 2668–2682.
