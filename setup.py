#!/usr/bin/env python
"""
setup.py file for SWIG ccdoes

This file is part of Sprout.

    Sprout is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sprout is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sprout.  If not, see <http://www.gnu.org/licenses/>.
     
   Copyright 2013, Michihiro Takami, Hyosun Kim, Jennifer Karr
   All rights reserved.   
   Contact: Jennifer Karr (jkarr@asiaa.sinica.edu.tw)


   -----

   A routine for setting things up to compile the C code and wrap
   it with the python code. 

   Version 1.0 (release) Last Updated May 25 2013


"""

from distutils.core import setup, Extension

name = 'ccodes'

my_module = Extension('_'+name, sources=[name+'_wrap.c',
                                         'VectorRotate.c',
                                         'StokesMuller.c',
                                         'NewPosition_func.c',
                                         'Initialize.c',
                                         'mt19937ar.c'], )

setup (name = name,
       version     = '1',
       author      = "Hyosun Kim",
       description = """ swig module ccodes """,
       ext_modules = [my_module],
       py_modules  = [name],
       )

