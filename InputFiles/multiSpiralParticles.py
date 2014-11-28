#!/usr/bin/env python
'''
This example shows how you can add multiple set of particles to a swarm. 
Here the particles are layed out in a spiral configuration.

Note that Scipy / Numpy are not required for this example
'''

import underworld
import math
import underworld.c_arrays
import underworld.StgDomain as StgDomain
import underworld.PICellerator as PICellerator
import underworld.Underworld as Underworld

# init using simple models
underworld.Init("RayleighTaylorBenchmark.xml PICellerator/PassiveTracerSwarm.xml")

# grab the dict
stgdict = underworld.GetCurrentDictionary()

# set to initialise and solve
stgdict["maxTimeSteps"]=250
stgdict["pauseToAttachDebugger"]=0
stgdict["checkpointEvery"]=50
stgdict["dumpDataEvery"]=10

# mod domain
stgdict["minX"] = -1.
stgdict["maxX"] =  1.
stgdict["minY"] = -1.
stgdict["maxY"] =  1.

# Add new variable to the passive swarm

colourVariable = underworld.NewComponentEntryInDictionary(name="colourVariable", type="Variable", globalDict=stgdict)
colourVariable["Rank"] = "Scalar"
colourVariable["DataType"] = "Double"

colourSwarmVariable = underworld.NewComponentEntryInDictionary(name="colourSwarmVariable", type="SwarmVariable", globalDict=stgdict)
colourSwarmVariable["Swarm"] = "passiveTracerSwarm"
colourSwarmVariable["Variable"] = "colourVariable"

stgdict["components"]["particleDots"]["ColourVariable"] = "colourSwarmVariable"


# lets get ride of manual particle layout
del stgdict["components"]["passiveTracerParticleLayout"]

# mod viz
stgdict["components"]["particleDots"]["Swarm"] = "passiveTracerSwarm"   # set lucSwarmViewer to viz passive swarm
stgdict["components"]["densityTitle"]["string"] = "Teleportation"    # rename this guy
stgdict["components"]["window"]["Viewport"] = "ParticleDensityVP"       # only viz this (now misnamed) viewport
stgdict["components"]["densityColourMap"]["colours"] = "Red Yellow Green Blue Purple"
stgdict["components"]["densityColourMap"]["dynamicRange"] = "False"
stgdict["components"]["densityColourMap"]["minimum"] = "1.0"
stgdict["components"]["densityColourMap"]["maximum"] = "100.0"


# don't forget to set the dict back again to affect the above changes
underworld.SetDictionary(stgdict)

underworld.Construct()

# for y in range(-100,100):
# 	t = (0.0, y/100.0, 0.0 )
# 	print "{} - {}".format( y/100.0, StgDomain.Stg_Shape_IsCoordInside(teleporterBox, t) )
 

##  lets reinit swarm guys
#   grab the passiveTracerSwarm guy
swarm = underworld.GetLiveComponent("passiveTracerSwarm")

#   and extract the variables stored there (would be nicer to match this against what we want)
variablesOnPassiveSwarm = underworld.Swarm_GetVariables(swarm)
positionXSwarmVariable = variablesOnPassiveSwarm[1][1]
positionYSwarmVariable = variablesOnPassiveSwarm[2][1]
colourSwarmVariable = underworld.GetLiveComponent("colourSwarmVariable")

particleCount = 500
dim = 2 
revCount = 5.0
radius = 0.075

valuePtr = underworld.c_arrays.DoubleArray(particleCount*dim)  # allocate c array for particle coords

for ii in range(0,particleCount):                             # set particle coords in spiral config
	factor =  float(ii)/float(particleCount)
	valuePtr[dim*ii+0] = radius * factor * math.cos(revCount*2.*math.pi*factor)
	valuePtr[dim*ii+1] = radius * factor * math.sin(revCount*2.*math.pi*factor)

underworld.PICellerator.GeneralSwarm_AddParticlesFromCoordArray(swarm, particleCount, dim, valuePtr.cast())   # ok, add particles

for ii in range(0,swarm.particleLocalCount):
	underworld.SwarmVariable_SetValueAt( colourSwarmVariable, ii, [ 1 ] )   		


# Timestepping loop
for step in range(0,100):
	# Periodically add another spiral into the passive array 
	if (step % 10 == 0):
		oldNumberOfLocalParticles=swarm.particleLocalCount
		underworld.PICellerator.GeneralSwarm_AddParticlesFromCoordArray(swarm, particleCount, dim, valuePtr.cast()) 
		for ii in range(oldNumberOfLocalParticles,swarm.particleLocalCount):
	 		underworld.SwarmVariable_SetValueAt( colourSwarmVariable, ii, [ step ] )   		

	underworld.Step(steps=1)

		

underworld.Finalise()

quit()






