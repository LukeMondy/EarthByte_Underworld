#!/usr/bin/env python
'''
This example shows how you can add a set of particles to a swarm. 

Here the particles are layed out in a spiral configuration.

Note that Scipy / Numpy are required for this example
'''

import underworld

# init using simple models
underworld.Init("BuoyancyDrivenVanilla.xml PICellerator/PassiveTracerSwarm.xml")

# grab the dict
stgdict = underworld.GetCurrentDictionary()

# set to initialise and solve
stgdict["maxTimeSteps"]=-1
stgdict["pauseToAttachDebugger"]=0
stgdict["checkpointEvery"]=1

# mod domain
stgdict["minX"] = -1.
stgdict["maxX"] =  1.
stgdict["minY"] = -1.
stgdict["maxY"] =  1.

# lets get ride of manual particle layout
del stgdict["components"]["passiveTracerParticleLayout"]

# mod viz
stgdict["components"]["particleDots"]["Swarm"] = "passiveTracerSwarm"   # set lucSwarmViewer to viz passive swarm
del stgdict["components"]["particleDots"]["ColourVariable"]             # get rid of these as we do not have a swarmvariable on the passive swarm
del stgdict["components"]["particleDots"]["ColourMap"]
stgdict["components"]["densityTitle"]["string"] = "Spiral Particles"    # rename this guy
stgdict["components"]["window"]["Viewport"] = "ParticleDensityVP"       # only viz this (now misnamed) viewport

# don't forget to set the dict back again to affect the above changes
underworld.SetDictionary(stgdict)

underworld.Construct()

##  lets reinit swarm guys
#   grab the passiveTracerSwarm guy
swarm = underworld.GetLiveComponent("passiveTracerSwarm")

#underworld.Swarm_PrintVariables(swarm)
#variables = underworld.Swarm_GetVariables(swarm)

particleCount = 10000
dim = 2 
revCount = 10.

import underworld.c_arrays
valuePtr = underworld.c_arrays.DoubleArray(particleCount*dim)  # allocate c array for particle coords

import math
for ii in range(0,particleCount):                             # set particle coords in spiral config
	factor =  float(ii)/float(particleCount)
	valuePtr[dim*ii+0] = factor*math.cos(revCount*2.*math.pi*factor)
	valuePtr[dim*ii+1] = factor*math.sin(revCount*2.*math.pi*factor)

underworld.PICellerator.GeneralSwarm_AddParticlesFromCoordArray(swarm, particleCount, dim, valuePtr.cast())   # ok, add particles

underworld.RunMainLoop()

underworld.Finalise()
