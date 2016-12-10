#Test necessary python packages exits

exitCode = 0

try:
  #from pylab import *
  import scipy
  import scipy.io
  import scipy.interpolate
  import scipy.optimize
except:
  exitCode =1

exit(exitCode)