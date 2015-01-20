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
model_input_files = pwd+"/4_input.xml "
# velocity BC def

model_input_files += pwd+"/kspinterface.xml -scr_ksp_type bicg "
#model_input_files += pwd+"/../../vcs/annulus.Periodic.freeSlip.xml "
model_input_files += pwd+"/../../vcs/annulus.Periodic.noSlip.xml "
# temperature BC def
model_input_files += pwd+"/../../vcs/temperatureBCs.xml "

# init model
underworld.Init(model_input_files)

# grab the dict
stgdict = underworld.dictionary.GetDictionary()

# now make changes before run time
radial_elementRes=40
Ra = 1e4
Di = 0.25
stgdict["Ra"]=Ra
stgdict["Di"]=Di
stgdict["outputPath"]="./output/"

# important restart bits
#stgdict["checkpointReadPath"]="/home/JulianGiordani/scratch/Ra_1e4_"+str(radial_elementRes)+"_parall/"
#stgdict["restartTimestep"]=13400


# set to initialise and solve
stgdict["elementResI"]=radial_elementRes
stgdict["elementResJ"]=round(10*radial_elementRes)

stgdict["maxTimeSteps"]=1
stgdict["checkpointEvery"]=1

#stgdict["pauseToAttachDebugger"]=10

# don't forget to set the dict back again to affect the above changes
underworld.dictionary.SetDictionary(stgdict)

underworld.Construct()

##  lets reinit swarm guys
#   grab the TemperatureField
tfield = underworld._stgermain.GetLiveComponent("TemperatureField")

# get the number of local nodes on the temperature mesh
nLocalNodes = underworld.libUnderworld.StgDomain.Mesh_GetLocalSize( tfield.feMesh, 0 ) # 0 represents the 0th topological element of the mesh, ie the nodes

from libUnderworld import c_arrays
from libUnderworld import StgDomain
from libUnderworld import StgFEM
import math

maxR = 2.22
minR = 1.22

cVal = c_arrays.DoubleArray(1)
for ii in range( 0, nLocalNodes ):
   result = StgDomain.Mesh_GetVertex(tfield.feMesh, ii )
   ptr = c_arrays.DoubleArray.frompointer(result)
   pos = [ptr[0], ptr[1]] # make a position vector, assume 2D
   r = math.sqrt(pos[0]*pos[0] + pos[1]*pos[1])
   scaledR = (maxR - r) / ( maxR - minR ) 
# calculate a perturbation on the temperature
   angle = math.atan2( pos[1], pos[0] )
   temp = (scaledR) + 0.1 * (math.cos(4*angle) * math.sin(math.pi*scaledR) ) 
   cVal[0] = temp
# set the temperature
   StgFEM.FeVariable_SetValueAtNode( tfield, ii, cVal.cast() )

underworld.RunMainLoop()

underworld.Finalise()
