#!/usr/bin/env python
import os, sys, subprocess

guys =  [ 
"StGermain",
"StgDomain",
"StgFEM",
"PICellerator",
"Underworld",
"gLucifer",
"gLucifer_Viewer",
"Geothermal",
"ImportersToolbox",
"Solvers",
"Spherical",
"Viscoelastic",
"libUnderworldPy" 
]

guy = None
try:
	guy = sys.argv[1]
	guys = [guy.strip('/')]
except:
	pass

os.environ["BUILDNOW"]="DONOTBUILDANYTHING"
sconsBin = os.path.join('config', 'scons', 'scons.py')
subp = subprocess.Popen(
    sconsBin + ' ' + ' '.join(sys.argv[1:]), shell=True
)
subp.wait()

for guy in guys:
	os.environ["BUILDNOW"]=guy
	sconsBin = os.path.join('config', 'scons', 'scons.py')

	print "Building " + guy
	subp = subprocess.Popen(
	    sconsBin + ' ' + ' '.join(sys.argv[1:]), shell=True
	)
	
	if subp.wait() != 0:
		print "-------------------------------------------"
		print "Error while compiling " + guy
		print "-------------------------------------------"
		sys.exit(1)
