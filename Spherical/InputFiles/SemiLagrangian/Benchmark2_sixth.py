#!/usr/bin/env python
'''
This example shows how you can modify an IC
'''

import os
import subprocess
import underworld

home=os.getenv("HOME")
pwd=os.getcwd();
# init using simple models

# standard model
model_input_files = pwd+"/Benchmark2_sixth.xml "

# init model
underworld.Init(model_input_files)

# grab the dict
stgdict = underworld.dictionary.GetDictionary()

# now make changes before run time
Ra = 1e4
restartTimestep = 108
stgdict["Ra"]=Ra
stgdict["outputPath"]="./benchmark2"

# important restart bits
#stgdict["checkpointReadPath"]="/home/JulianGiordani/scratch/Ra_1e4_"+str(radial_elementRes)+"_parall/"
#stgdict["restartTimestep"]=13400

# set to initialise and solve
stgdict["maxTimeSteps"]=800
stgdict["checkpointEvery"]=1
stgdict["restartTimestep"]=restartTimestep

stgdict["pauseToAttachDebugger"]=0

# don't forget to set the dict back again to affect the above changes
#underworld.dictionary.SetDictionary(stgdict)

underworld.Construct()

from libUnderworld import c_arrays
from libUnderworld import StgDomain
from libUnderworld import StgFEM
import math
import numpy as np

#   grab the TemperatureField
tfield = underworld._stgermain.GetLiveComponent("TemperatureField")

# get the number of local nodes on the temperature mesh
nLocalNodes = underworld.libUnderworld.StgDomain.Mesh_GetLocalSize( tfield.feMesh, 0 ) # 0 represents the 0th topological element of the mesh, ie the nodes
mesh = tfield.feMesh.getAsNumpyArray()

maxR = 6.0
minR = 3.0
dx   = minR - maxR
d2x  = minR*minR - maxR*maxR

temp = np.zeros( mesh.shape[0] )

if restartTimestep == 0:
	cVal = c_arrays.DoubleArray(1)
	for ii in range( 0, nLocalNodes ):
		r        = math.sqrt(mesh[ii][0]*mesh[ii][0] + mesh[ii][1]*mesh[ii][1] + mesh[ii][2]*mesh[ii][2])
		eta      = math.atan2( mesh[ii][0], mesh[ii][2] )
		zeta     = math.atan2( mesh[ii][1], mesh[ii][2] )
		#scaledR = (maxR - r) / ( maxR - minR ) 
# calculate a perturbation on the temperature
		height   = 0.25*(math.sin(2*eta) + 1.0)*(math.sin(2*zeta) + 1.0)
	        #temp    = scaledR + 0.1*math.sin(math.pi*scaledR)*height
		a        = height
	        b        = (1.0 - a*d2x)/dx
		c        = maxR*(a*d2x - 1.0)/dx - a*maxR*maxR
	        temp[ii] = a*r*r + b*r + c
	        if( temp[ii] ) < 0.0:
			temp[ii] = 0.0

	cVal = c_arrays.DoubleArray(1)
	for ii in range( 0, nLocalNodes ):
		cVal[0] = temp[ii]
	# set the temperature
		StgFEM.FeVariable_SetValueAtNode( tfield, ii, cVal.cast() )

underworld.RunMainLoop()

underworld.Finalise()
