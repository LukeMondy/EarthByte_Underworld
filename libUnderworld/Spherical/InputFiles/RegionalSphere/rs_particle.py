#!/usr/bin/python

# coding: utf-8

import os
import subprocess
import math
import underworld

pwd=os.getcwd();
# init using simple models

# standard model
model_input_files = pwd+"/RS_ball.xml "

### command line modifications ###
# model_input_files += " -Uzawa_velSolver_ksp_type preonly -Uzawa_velSolver_pc_type lu " 

# read xml
underworld.Init(model_input_files)

# grab the dict
stgdict = underworld.dictionary.GetDictionary()


# set to global parameters
stgdict["maxTimeSteps"]=0
stgdict["checkpointEvery"]=1

stgdict["pauseToAttachDebugger"]=0
stgdict["outputPath"]="./redef_particles"

# set the dictionary again in Underworld
underworld.dictionary.SetDictionary(stgdict)

underworld.Construct()

def inSphere(xyz):
    """
     Return 1 if xyz location is inside circle of radius 1 centred at (2,1,3)
    """
    x = xyz[0]-2
    y = xyz[1]-1
    z = xyz[2]-3
    
    rad = x*x + y*y + z*z
    
    if rad < 1:
        return True
    else:
        return False
    


# now we get the material swarms and material variables
matswarm = underworld._stgermain.GetLiveComponent("materialSwarm")
matVars = underworld.swarms.tools.Swarm_GetVariables(matswarm)

# print the variables
[x[0] for x in matVars]


for x in matVars:
    if x[0] == 'materialSwarm-MaterialIndex':
        print matVars.index(x)


posVar = matVars[4][1]
matVar = matVars[5][1]



for ii in range(0,matswarm.particleLocalCount):
	# grab the position
    pos = underworld.swarms.tools.SwarmVariable_GetValueAt( posVar, ii )
    #pointVal = map_rbfs(pos, slab_funcs, 10.0, 9.75, 9.25, 9.0, 0.0)
    if(inSphere(pos)):
        underworld.swarms.tools.SwarmVariable_SetValueAt( matVar , ii, [1] )



underworld.RunMainLoop()

