#!/usr/bin/env python
'''
This example shows how 2 swarms can interact using radial basis functions
to detect their proximity. The particles in the passive tracer swarm are
teleported vertically when they come close to the particles in the RBF swarm.
(see teleportationSwarm for the more simple example of this).

This also shows how to initialise the RBF swarm from python and ensure
that the basis functions are properly constructed - this is an additional
step required for a passive Swarm.

Note that Scipy / Numpy are not required for this example
'''

import uwpytools
import math
import uwpytools.c_arrays
import uwpytools.StgDomain as StgDomain
import uwpytools.PICellerator as PICellerator
import uwpytools.Underworld as Underworld

# init using simple models
uwpytools.InitWithArgs("RayleighTaylorBenchmark.xml")

# grab the dictionary
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


# Add a passive swarm and the required advector component
# Add a new variable for the visuals

passiveTracerSwarm = uwpytools.NewComponentEntryInStgDict( globalDict=stgdict,
								 name="passiveTracerSwarm", 
								 Type="GeneralSwarm", 
								 CellLayout= "elementCellLayout",
								 FeMesh= "elementMesh",
								 ParticleCommHandlers= [ "pMovementHandler" ] )

passiveTracerAdvector = uwpytools.NewComponentEntryInStgDict( globalDict=stgdict,
								 name="passiveTracerAdvector", 
								 Type="SwarmAdvector", 
								 Swarm= "passiveTracerSwarm",
								 TimeIntegrator= "timeIntegrator",
								 VelocityField= "VelocityField",
								 PeriodicBCsManager= "periodicBCsManager" )

colourVariable = uwpytools.NewComponentEntryInStgDict( globalDict=stgdict,
								 name="passiveTracerSwarm-colourVariable", 
								 Type="Variable", 
								 Rank= "Scalar",
								 DataType= "Double" )
 
colourSwarmVariable = uwpytools.NewComponentEntryInStgDict( globalDict=stgdict,
								 name="passiveTracerSwarm-colour", 
								 Type="SwarmVariable",
								 Swarm= "passiveTracerSwarm",
								 Variable= "passiveTracerSwarm-colourVariable")

# push this into the existing viz component
stgdict["components"]["particleDots"]["ColourVariable"] = colourSwarmVariable["name"]

# Add an RBFSwarm (no advector at the moment)
# Add a variable and a shape 
# It might be useful to bundle all of this up into one call

rbfTracerSwarm = uwpytools.NewComponentEntryInStgDict( globalDict=stgdict,
								 name="rbfTracerSwarm", 
								 Type="Swarm", 
								 CellLayout = "CLLCellLayout",
								 FeMesh = "elementMesh" )
								 # ParticleCommHandlers = [ "pMovementHandler" ] )   # only needed if moving 


RBFCellLayout = uwpytools.NewComponentEntryInStgDict( globalDict=stgdict,
								 name="CLLCellLayout",
								 Type="CLLCellLayout",
								 GeometryMesh="linearMesh",
								 CellSize="0.025",
 								 MeshDim="2" )

rbfManager = uwpytools.NewComponentEntryInStgDict( globalDict=stgdict,
								 name="rbfManager", 
								 Type="RBFManager", 
								 ParticleSupportRadius="0.05",
								 RBFdim="2",
								 RBFSwarm=rbfTracerSwarm["name"] )


rbfTracerVariable = uwpytools.NewComponentEntryInStgDict( globalDict=stgdict,
								 name="rbfTracerSwarm-valueVariable", 
								 Type="Variable", 
								 Rank= "Scalar", 
								 DataType= "Double" )

rbfTracerSwarmVariable = uwpytools.NewComponentEntryInStgDict( globalDict=stgdict,
								 name="rbfTracerSwarm-value", 
								 Type="SwarmVariable", 
								 Swarm= rbfTracerSwarm["name"],
								 Variable= rbfTracerVariable["name"] )

rbfFieldVariable = uwpytools.NewComponentEntryInStgDict( globalDict=stgdict, 
								 name="rbfFieldVariable",
								 Type="RBFFieldVariable",
								 RBFManager=rbfManager["name"],
								 SwarmVariable=rbfTracerSwarmVariable["name"] 
								 # SwarmVariable="rbfTracerSwarm-PositionY"
								 )

# New heightfield shape (Does this work for a 2D model ?) 

rbfHeightFieldShape = uwpytools.NewComponentEntryInStgDict( globalDict=stgdict, 
								 name="rbfHeightFieldShape",
								 Type="BelowHeightField", 
								 HeightField=rbfFieldVariable["name"] )

# New rbf field-value shape 

rbfFieldValueShape = uwpytools.NewComponentEntryInStgDict( globalDict=stgdict, 
								 name="rbfFieldValueShape",
								 Type="FieldValueShape", 
								 ValueField=rbfFieldVariable["name"], 
								 LowerLimit="0.5", 
								 UpperLimit="1.5"	 )


rbfParticleDots = uwpytools.NewComponentEntryInStgDict( globalDict=stgdict,
			name="rbfParticleDots", Type="lucSwarmViewer", Swarm=rbfTracerSwarm["name"], 
			Colour="Black", pointSize = "2.0" )

stgdict["components"]["ParticleDensityVP"]["DrawingObject"].append(rbfParticleDots["name"])

# mod viz
stgdict["components"]["particleDots"]["Swarm"] = "passiveTracerSwarm"   # set lucSwarmViewer to viz passive swarm
stgdict["components"]["densityTitle"]["string"] = "Watch Me Teleport !"    # rename this guy
stgdict["components"]["window"]["Viewport"] = "ParticleDensityVP"       # only viz this (now misnamed) viewport
stgdict["components"]["densityColourMap"]["colours"] = "Red Yellow Green Blue Purple"
stgdict["components"]["densityColourMap"]["dynamicRange"] = "False"
stgdict["components"]["densityColourMap"]["minimum"] = "1.0"
stgdict["components"]["densityColourMap"]["maximum"] = "250.0"


# Push the dictionary back to StGermain before construct / build phase
uwpytools.SetDictionary(stgdict)

uwpytools.Construct()
uwpytools.BuildAndInitialise()



# Components are now live !

# We need the actual Shape component so we can call its functions
teleporterBox = uwpytools.GetLiveComponent(rbfFieldValueShape["name"])

# Set up rbf swarm, and put in some points

dim=2
maxjj = 8
maxii = 250 
particleCount=maxii*maxjj	

rbfValuePtr = uwpytools.c_arrays.DoubleArray(particleCount*dim)  # allocate c array for particle coords

rbfSwarm = uwpytools.GetLiveComponent(rbfTracerSwarm["name"])
varDict = uwpytools.Swarm_GetVariablesAsDict(rbfSwarm)
positionXSwarmVariable = varDict["rbfTracerSwarm-PositionX"]["swarmVariable"]
positionYSwarmVariable = varDict["rbfTracerSwarm-PositionY"]["swarmVariable"]
valueSwarmVariable     = varDict["rbfTracerSwarm-value"    ]["swarmVariable"]

# uwpytools.PrettyDictionaryPrint(varDict, indent=3)

particleNumber=0
for ii in range(0,maxii):
	for jj in range(0,maxjj):
		rbfValuePtr[dim*particleNumber+0] = 0.75 * float(ii) / maxii - 0.75 + 0.05 * float(jj) / maxjj;
		rbfValuePtr[dim*particleNumber+1] = float(jj) / maxjj - 0.25 + 0.015 * math.sin(20*math.pi*float(ii) / maxii ) 
		particleNumber += 1	

# for ii in range(0,particleCount):      # set particle coords in spiral config
# 	factor =  float(ii)/float(particleCount)
# 	rbfValuePtr[dim*ii+0] = radius * factor * math.cos(revCount*2.*math.pi*factor)
# 	rbfValuePtr[dim*ii+1] = radius * factor * math.sin(revCount*2.*math.pi*factor)  + 0.25

uwpytools.PICellerator.GeneralSwarm_AddParticlesFromCoordArray(rbfSwarm, particleCount, dim, rbfValuePtr.cast())   # ok, add particles

for ii in range(0,rbfSwarm.particleLocalCount):
	uwpytools.SwarmVariable_SetValueAt( valueSwarmVariable, ii, [ 1.0 ] )   	


# Now we have to rebuild the RBF stuff as a result of moving the particles

rbfManagerCPT = uwpytools.GetLiveComponent(rbfManager["name"])
Underworld.RBFManager_CalculateParticleDensities(rbfManagerCPT)		


#   grab the passiveTracerSwarm and find the variables it contains

swarm = uwpytools.GetLiveComponent("passiveTracerSwarm")
varDict = uwpytools.Swarm_GetVariablesAsDict(swarm)
positionXSwarmVariable = varDict["passiveTracerSwarm-PositionX"]["swarmVariable"]
positionYSwarmVariable = varDict["passiveTracerSwarm-PositionY"]["swarmVariable"]
colourSwarmVariable =    varDict["passiveTracerSwarm-colour"   ]["swarmVariable"]

particleCount = 2000
dim = 2 
revCount = 10.0
radius = 0.15

valuePtr = uwpytools.c_arrays.DoubleArray(particleCount*dim)  # allocate c array for particle coords

for ii in range(0,particleCount):                             # set particle coords in spiral config
	factor =  float(ii)/float(particleCount)
	valuePtr[dim*ii+0] = radius * factor * math.cos( revCount*2.*math.pi*factor)
	valuePtr[dim*ii+1] = radius * factor * math.sin(-revCount*2.*math.pi*factor) - 0.5 

uwpytools.PICellerator.GeneralSwarm_AddParticlesFromCoordArray(swarm, particleCount, dim, valuePtr.cast())   # ok, add particles

for ii in range(0,swarm.particleLocalCount):
	uwpytools.SwarmVariable_SetValueAt( colourSwarmVariable, ii, [ 1 ] )   		

RBFFieldVariable = uwpytools.GetLiveComponent(rbfFieldVariable["name"])

# Timestepping loop
for step in range(0,250):
	# Periodically add another spiral into the passive array 
	if (step % 25 == 0):
		oldNumberOfLocalParticles=swarm.particleLocalCount
		uwpytools.PICellerator.GeneralSwarm_AddParticlesFromCoordArray(swarm, particleCount, dim, valuePtr.cast()) 
		for ii in range(oldNumberOfLocalParticles,swarm.particleLocalCount):
	 		uwpytools.SwarmVariable_SetValueAt( colourSwarmVariable, ii, [ step ] )   		

	# teleporting box 
	ii = 0 		
	while ii < swarm.particleLocalCount:
		x = uwpytools.SwarmVariable_GetValueAt( positionXSwarmVariable , ii )
		y = uwpytools.SwarmVariable_GetValueAt( positionYSwarmVariable , ii )

		## We could use this directly ... 
		# value, result = uwpytools.FieldVariable_InterpolateValueAt( RBFFieldVariable, (x[0],y[0],0.0) )
		## etc etc

		# or use the Shape routines ... 
		if StgDomain.Stg_Shape_IsCoordInside(teleporterBox, (x[0],y[0],0.0)):
			newY = y[0]
			while StgDomain.Stg_Shape_IsCoordInside(teleporterBox, (x[0],newY,0.0)):
				newY += 0.01
			uwpytools.SwarmVariable_SetValueAt( positionYSwarmVariable, ii, [ newY ] ) 	
			
		ii += 1

	# in parallel we have to use this guy and not the particle-by-particle one		
	StgDomain.Swarm_UpdateAllParticleOwners(swarm) 

	uwpytools.Step(steps=1)

uwpytools.Finalise()

quit()






