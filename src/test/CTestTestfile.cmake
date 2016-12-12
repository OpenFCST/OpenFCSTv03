################################################################################################################
#
# This file contains locations and names of the tests that will be run by CTest.
#
# Usage:
#
# Within the brackets, the first item is simply a name and can be anything that will describe the test being
# run. The second item is the location of the script that contains the test that is to be run.
#
# Authors: M. Moore, P. Wardlaw, A. Kosakian, C. Balen and M. Secanell, 2013-15
#
################################################################################################################
## -- Edit the address below
## -- file(STRINGS file.txt variable) reads lines from the file and
## -- stores them as a list in the variable. End-of-line characters are ignored.
file(STRINGS @CMAKE_INSTALL_PREFIX@/test/config.txt FLAGS)

## -- list(FIND list value variable) will return index of element
## -- in the list through variable if the value was found or -1 if it was not
list(FIND FLAGS "PETScON" PETSc_FLAG)
list(FIND FLAGS "DAKOTAON" Dakota_FLAG)
list(FIND FLAGS "--cores=" Cores_INDEX)

if(PETSc_FLAG GREATER -1)
    MATH(EXPR Cores_INDEX "${Cores_INDEX}+1") #increase index by 1 since number of cores to use is on the next line
    list(GET FLAGS ${Cores_INDEX} Cores_FLAG)
else()
    set(Cores_FLAG "1")
endif()
    
# --- Unit tests:
# NOTE: This case needs to be run first to create the log file
ADD_TEST(unit_tests "./unit_tests/unit_tests_case.sh")

# --- app_cathode ---
ADD_TEST(AppCathode>>PolarizationCurve "../examples/cathode/polarization_curve/regression/run_test.sh" COMMAND "--cores=${Cores_FLAG}")
ADD_TEST(AppAnode>>PolarizationCurve "../examples/anode/polarization_curve/regression/run_test.sh" COMMAND "--cores=${Cores_FLAG}")

# --- app_pemfc ---
ADD_TEST(AppPemfc>>Homogeneous "../examples/Pemfc/homogeneous/regression/run_test.sh" COMMAND "--cores=${Cores_FLAG}")
ADD_TEST(AppPemfc>>Agglomerate "../examples/Pemfc/agglomerate/regression/run_test.sh" COMMAND "--cores=${Cores_FLAG}")
ADD_TEST(AppPemfc>>GradedElectrode "../examples/Pemfc/graded_electrode/regression/run_test.sh" COMMAND "--cores=${Cores_FLAG}")

# --- app_nonisothermal ---
ADD_TEST(PemfcNIThermal>>PolarizationCurve "../examples/PemfcNIThermal/polarization_curve/regression/run_test.sh" COMMAND "--cores=${Cores_FLAG}")

# --- app_nonisothermal_two_phase ---
ADD_TEST(AppPemfcTPSaturation>>PolarizationCurve "../examples/PemfcTPSaturation/polarization_curve/regression/run_test.sh") #Do not run with multiple cores, else MUMPS can't solve at 0.5V

# --- app_diffusion ---
ADD_TEST(AppDiffusion3D>>Analysis "../examples/microscale/test/regression/run_test.sh")

# --- app_ohmic
ADD_TEST(AppOhmic>>Analysis "../examples/ohmic/analysis/regression/run_test.sh" COMMAND "--cores=${Cores_FLAG}")

############################
### PETSc and Dakota:
############################
## -- Parallel tests
IF(PETSc_FLAG GREATER -1)
    ADD_TEST(AppCathode>>PolarizationCurve>>Parallel "../examples/cathode/parallel/regression/run_test.sh" COMMAND "--cores=${Cores_FLAG}")
ENDIF() # PETSc_FLAG

# -- Dakota integration test:
IF(Dakota_FLAG GREATER -1)
    ADD_TEST(AppCathode>>Optimization "../examples/cathode/optimization/regression/run_test.sh")
ENDIF() # Dakota_FLAG

############################
############################
############################
