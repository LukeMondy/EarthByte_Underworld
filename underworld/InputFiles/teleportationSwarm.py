#!/usr/bin/env python
'''
This example shows how you can detect when tracers are within
a shape after the simulation has been initialised. It also
shows how to move tracers or delete them when such conditions
are met and update the swarm state appropriately

Note that Scipy / Numpy are not required for this example
'''

import uwpytools
import math
import uwpytools.c_arrays
import uwpytools.StgDomain as StgDomain
import uwpytools.PICellerator as PICellerator
import uwpytools.Underworld as Underworld

# init using simple models
uwpytools.InitWithArgs("RayleighTaylorBenchmark.xml PICellerator/PassiveTracerSwarm.xml")

# grab the dict
stgdict = uwpytools.GetCurrentDictionary()

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

colourVariable = uwpytools.NewComponentEntryInDictionary(name="passiveTracerSwarm-colourVariable", Type="Variable", globalDict=stgdict)
colourVariable["Rank"] = "Scalar"
colourVariable["DataType"] = "Double"

colourSwarmVariable = uwpytools.NewComponentEntryInDictionary(name="passiveTracerSwarm-colour", Type="SwarmVariable", globalDict=stgdict)
colourSwarmVariable["Swarm"] = "passiveTracerSwarm"
colourSwarmVariable["Variable"] = "passiveTracerSwarm-colourVariable"

stgdict["components"]["particleDots"]["ColourVariable"] = "passiveTracerSwarm-colour"


# lets get ride of manual particle layout
del stgdict["components"]["passiveTracerParticleLayout"]

# mod viz
stgdict["components"]["particleDots"]["Swarm"] = "passiveTracerSwarm"   # set lucSwarmViewer to viz passive swarm
stgdict["components"]["densityTitle"]["string"] = "Watch Me Teleport !"    # rename this guy
stgdict["components"]["window"]["Viewport"] = "ParticleDensityVP"       # only viz this (now misnamed) viewport
stgdict["components"]["densityColourMap"]["colours"] = "Red Yellow Green Blue Purple"
stgdict["components"]["densityColourMap"]["dynamicRange"] = "False"
stgdict["components"]["densityColourMap"]["minimum"] = "1.0"
stgdict["components"]["densityColourMap"]["maximum"] = "100.0"


# Particle teleporting box

teleporterBox = uwpytools.NewComponentEntryInDictionary(name="teleporterBox", Type="Box", globalDict=stgdict)
teleporterBox["startX"] = -0.2
teleporterBox["endX"]   = 0.0
teleporterBox["startY"] = 0.25
teleporterBox["endY"]   = 0.45


# don't forget to set the dict back again to affect the above changes
uwpytools.SetDictionary(stgdict)

uwpytools.Construct()
uwpytools.BuildAndInitialise()

# Components are now live !

# We need the actual Shape component so we can call its functions
teleporterBox = uwpytools.GetLiveComponent("teleporterBox")

#   grab the passiveTracerSwarm and find the variables it contains

swarm = uwpytools.GetLiveComponent("passiveTracerSwarm")
varDict = uwpytools.Swarm_GetVariablesAsDict(swarm)
positionXSwarmVariable = varDict["passiveTracerSwarm-PositionX"]["swarmVariable"]
positionYSwarmVariable = varDict["passiveTracerSwarm-PositionY"]["swarmVariable"]
colourSwarmVariable =    varDict["passiveTracerSwarm-colour"   ]["swarmVariable"]

particleCount = 1000
dim = 2 
revCount = 10.0
radius = 0.075

valuePtr = uwpytools.c_arrays.DoubleArray(particleCount*dim)  # allocate c array for particle coords

for ii in range(0,particleCount):                             # set particle coords in spiral config
	factor =  float(ii)/float(particleCount)
	valuePtr[dim*ii+0] = radius * factor * math.cos(revCount*2.*math.pi*factor)
	valuePtr[dim*ii+1] = radius * factor * math.sin(revCount*2.*math.pi*factor)

uwpytools.PICellerator.GeneralSwarm_AddParticlesFromCoordArray(swarm, particleCount, dim, valuePtr.cast())   # ok, add particles

for ii in range(0,swarm.particleLocalCount):
	uwpytools.SwarmVariable_SetValueAt( colourSwarmVariable, ii, [ 1 ] )   		


# Timestepping loop
for step in range(0,100):
	# Periodically add another spiral into the passive array 
	if (step % 10 == 0):
		oldNumberOfLocalParticles=swarm.particleLocalCount
		uwpytools.PICellerator.GeneralSwarm_AddParticlesFromCoordArray(swarm, particleCount, dim, valuePtr.cast()) 
		for ii in range(oldNumberOfLocalParticles,swarm.particleLocalCount):
	 		uwpytools.SwarmVariable_SetValueAt( colourSwarmVariable, ii, [ step ] )   		


	# teleporting box  		
	while ii < swarm.particleLocalCount:
		x = uwpytools.SwarmVariable_GetValueAt( positionXSwarmVariable , ii )
		y = uwpytools.SwarmVariable_GetValueAt( positionYSwarmVariable , ii )
		if StgDomain.Stg_Shape_IsCoordInside(teleporterBox, (x[0],y[0],0.0)):
			newY = y[0]
			while StgDomain.Stg_Shape_IsCoordInside(teleporterBox, (x[0],newY,0.0)):
				newY += 0.01
			uwpytools.SwarmVariable_SetValueAt( positionYSwarmVariable, ii, [ newY ] ) # exterminate, exterminate !!!	
			
			# OR .. teleport them into an unknown dimension (but be careful, Victor, they may clutter up the machine)
			# StgDomain.Swarm_DeleteParticle(swarm, ii)
			# Note - we use the while loop here because the Delete function will change the particle count

		ii += 1


	# in parallel we have to use this guy and not the particle-by-particle one		
	StgDomain.Swarm_UpdateAllParticleOwners(swarm) 

	uwpytools.Step(steps=1)

uwpytools.Finalise()

quit()






