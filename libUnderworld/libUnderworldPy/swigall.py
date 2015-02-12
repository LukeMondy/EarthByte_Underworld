#!/usr/bin/env python
import sys, subprocess, shutil, os, glob	

dir = os.path.join( '..', 'build')

swigcommandbase = 'swig -Wextra -python -importall -ignoremissing'
swigcommandbase +=  ' -I'+os.path.join( dir, 'include')
swigcommandbase +=  ' -I'+os.path.join( dir, 'include', 'StGermain')
swigcommandbase +=  ' -I'+os.path.join( dir, 'include', 'StgDomain')
swigcommandbase +=  ' -I'+os.path.join( dir, 'include', 'StgFEM')
swigcommandbase +=  ' -I'+os.path.join( dir, 'include', 'PICellerator')
swigcommandbase +=  ' -I'+os.path.join( dir, 'include', 'Underworld')
swigcommandbase +=  ' -I'+os.path.join( dir, 'include', 'gLucifer')
swigcommandbase +=  ' -I'+os.path.join( dir, 'include', 'ImportersToolbox')
swigcommandbase +=  ' -I'+os.path.join( dir, 'include', 'petsc')
swigcommandbase +=  ' -I'+os.path.join( '..', 'ctools')

swigfiles = [ 
              "c_arrays.i",
              "c_pointers.i",
              "petsc.i",
              "StGermain_Tools.i",
              "gLucifer.i",
              "ImportersToolbox.i",
              "Underworld.i",
              "PICellerator.i",
              "StgFEM.i",
              "StgDomain.i",
              "StGermain.i",
              "analytic.i"
             ]

for swigfile in swigfiles:
	print ""
	print "executing command:"
	print swigcommandbase + ' ' + swigfile
	print ""
	subp = subprocess.Popen(swigcommandbase + ' ' + swigfile, shell=True )
	subp.wait()

