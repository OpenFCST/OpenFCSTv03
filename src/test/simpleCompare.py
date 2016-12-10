#! /usr/bin/env python2.7
"""
Brief:
  Simple comparison script which compares 1d vectors read in from text files.
  Vector members are end line character deliminated.

Usage:
  python simpleCompare.py solution_ref.dat solution.dat

Results:
  returns 0 if data matches within tolerance
  returns 1 otherwise

  by Phil Wardlaw, 2014

"""

import sys

if __name__ == "__main__":


    tolerance = 0.01 #  1%


    assert len(sys.argv) ==3, "Invalid number of arguments"

    # Arg 0 will be python file name
    file0 = open(sys.argv[1])
    file1 = open(sys.argv[2])
    data0 = []
    data1 = []

    for l in file0.readlines():
        data0.append(float(l))

    for l in file1.readlines():
        data1.append(float(l))


    if len(data0) != len(data1):
        print "Vector lengths do not match"
        exit(1)

    for i in range(len(data0)):
        if data0[i] > 10e-10:
            diff = abs(data1[i]-data0[i])/data0[i]
        else:
            diff = abs(data1[i]-data0[i])
        #print diff
        if diff > tolerance:
            print "Vector members " + str(i) + " do not match within tolerance"
            exit(1)

    #Test passed
    #print "All good"
    exit(0)