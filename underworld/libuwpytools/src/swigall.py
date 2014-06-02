#!/usr/bin/env python
import sys, subprocess, shutil, os, glob	

dir = os.environ.get('UW_DIR',None)
if dir is None:
	print 'UW_DIR environment variable not set! this must point to your Underworld build directory.'
	sys.exit()

swigcommandbase = 'swig -Wextra -python -importall -ignoremissing'
swigcommandbase +=  ' -I'+os.path.join( os.environ.get('UW_DIR'), 'include')
swigcommandbase +=  ' -I'+os.path.join( os.environ.get('UW_DIR'), 'include', 'StGermain')
swigcommandbase +=  ' -I'+os.path.join( os.environ.get('UW_DIR'), 'include', 'StgDomain')
swigcommandbase +=  ' -I'+os.path.join( os.environ.get('UW_DIR'), 'include', 'StgFEM')
swigcommandbase +=  ' -I'+os.path.join( os.environ.get('UW_DIR'), 'include', 'PICellerator')
swigcommandbase +=  ' -I'+os.path.join( os.environ.get('UW_DIR'), 'include', 'Underworld')
swigcommandbase +=  ' -I'+os.path.join( os.environ.get('UW_DIR'), 'include', 'gLucifer')
swigcommandbase +=  ' -I'+os.path.join( os.environ.get('UW_DIR'), 'include', 'petsc')

swigfiles = [ 
              "c_arrays.i",
              "c_pointers.i",
              "petsc.i",
              "gLucifer.i",
              "Underworld.i",
              "StgFEM.i",
              "StgDomain.i",
              "StGermain_Tools.i",
              "StGermain.i",
              "PICellerator.i"
             ]

for swigfile in swigfiles:
	print ""
	print "executing command:"
	print swigcommandbase + ' ' + swigfile
	print ""
	subp = subprocess.Popen(swigcommandbase + ' ' + swigfile, shell=True )
	subp.wait()
	fname = swigfile.split(".i")[0]
	fname = fname + ".py"
	subp = subprocess.Popen('mv ' + fname + ' ../../uwpytools/.', shell=True )

