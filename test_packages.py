#!/usr/bin/env python

"""

A quick script to check for the present of the required packages, and 
print a warning for missing packages. 


"""

alldone=0

try:
    import Tkinter
except:
    print "Tkinter not installed\n"
    alldone=1


try:
    import pylab
except:
    print "pylab not installed\n"
    alldone=1


try:
    import numpy
except:
    print "numpy not installed\n"
    alldone=1


try:
    import scipy
except:
    print "scipy not installed\n"
    alldone=1


try:
    import pyfits
except:
    print "pyfits not installed\n"
    alldone=1


try:
    import pickle
except:
    print "pickle not installed\n"
    alldone=1

try:
    import matplotlib
except:
    print "matplotlib not installed\n"
    alldone=1

try:
    import os
except:
    print "os not installed\n"
    alldone=1


try:
    import Pmw
except:
    print "Pmw not installed\n"
    alldone=1
try:
    import string
except:
    print "string not installed\n"
    alldone=1
try:
    import sys
except:
    print "sys not installed\n"
    alldone=1



if (alldone==0):
    print "All packages found\n"
