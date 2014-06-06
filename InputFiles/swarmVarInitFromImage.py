#!/usr/bin/env python
'''
This example shows how you can extract and/or set the values on particle swarm variables. 
Here we have use a raster image to set the density values for the simulation.

Note that Scipy / Numpy are required for this example
'''

import underworld
import underworld.c_arrays as c_arrays


# init with this vanilla model
underworld.InitWithArgs("BuoyancyDrivenVanilla.xml")

# grab the dict
stgdict = underworld.GetCurrentDictionary()

# set to initialise and solve
stgdict["maxTimeSteps"]=100
stgdict["pauseToAttachDebugger"]=0
stgdict["components"]["default_ppcManager"]["Gravity"] = 10.0
stgdict["components"]["window"]["Viewport"] = "ParticleDensityVP"

 
# don't forget to set the dict back again to affect the above changes
underworld.SetDictionary(stgdict)


# underworld.PrettyDictionaryPrint(stgdict, indent=3)


underworld.Construct()
underworld.BuildAndInitialise()


##  lets reinit swarm guys
#   grab the material swarm
swarm = underworld.GetLiveComponent("materialSwarm")
underworld.Swarm_PrintVariables(swarm)

# This is ugly, but gives us (2 way) access to the gravity array in the ppcManager
ppc_mgr = underworld.GetLiveComponent("default_ppcManager")
gravVector = c_arrays.DoubleArray_frompointer(ppc_mgr.gravity)

variables = underworld.Swarm_GetVariables(swarm)

# lets fire up smug moresi
from scipy import misc
smugMoresi = misc.imread('moresi.png')
# rotate
import numpy
smugMoresi = numpy.rot90(smugMoresi,k=3)

# setup scaling to image
xscaling = float(smugMoresi.shape[0]) / ( stgdict["maxX"] - stgdict["minX"] )
yscaling = float(smugMoresi.shape[1]) / ( stgdict["maxY"] - stgdict["minY"] )
densscaling = float(smugMoresi.ptp())  # the ptp function gets the peak-to-peak of the array

# get position variable.
posVariable  = variables[3][1]  # note that we know that the variable is the third in the list from inspecting the output to Swarm_PrintVariables
# lets grab the density variable directly from the dictionary, just to demonstrate the other option.
densVariable = underworld.GetLiveComponent("DensitySwarmVariable")

# ok, lets sweep image
for ii in range(0,swarm.particleLocalCount):
	# grab the position
 	pos = underworld.SwarmVariable_GetValueAt( posVariable, ii )
 	# scale to location in mat
 	scalposX = int( (pos[0]-stgdict["minX"]) * xscaling ) 
 	scalposY = int( (pos[1]-stgdict["minY"]) * yscaling ) 
 	densVal = ( densscaling - float(smugMoresi[scalposX][scalposY]) ) / densscaling  # note that white = largeValues, black = smallValues, hence we invert
 	# now set the density swarm variable
 	underworld.SwarmVariable_SetValueAt( densVariable, ii, [densVal] )  # note that we pass in a list of values
 	#print underworld.SwarmVariable_GetValueAt( densVariable, ii )

# and how does young moresi hold up under the weight of his own smugness? 

underworld.Step(steps=50)
for i in range(0,10):
	gravVector[1] *= -1.0
	underworld.Step(steps=100)
gravVector[1] *= -1.0

underworld.Finalise()
