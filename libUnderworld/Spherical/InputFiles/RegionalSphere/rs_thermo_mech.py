#!/usr/bin/python
'''
This example shows how you can modify a temperature IC
'''

import os
import subprocess
import underworld

pwd=os.getcwd();
# init using simple models

# standard model
model_input_files = pwd+"/RS_themal.xml "
### command line modifications ###
# model_input_files += " -Uzawa_velSolver_ksp_type preonly -Uzawa_velSolver_pc_type lu " 

# read xml
underworld.Init(model_input_files)

# grab the dict
stgdict = underworld.dictionary.GetDictionary()

# python variables to be used to set xml
maxR=6
minR=3
T_max = 1
T_min = 0
wallInput = { 'left' : T_max, 'right' : T_min }

# some python trickery to set the temperature BCs
# read the existing temperatureBCs xml inforomation
nTempBCs = len( stgdict["temperatureBCs"]["vcList"] )

for w_i in range(nTempBCs):
    # for each wall defined in xml 
    entry = stgdict["temperatureBCs"]["vcList"][w_i]
    # we overwrite the value as per the 'wallInput' python dictionary above
    if wallInput.has_key(entry["wall"]):
        entry['variables'][0]['value'] = wallInput[ entry["wall"] ]


# now make changes before run time
stgdict["outputPath"]="./rs_thermo_pic"
stgdict["minX"]=minR
stgdict["maxX"]=maxR

nTempWallBCs = len( stgdict["temperatureBCs"]["vcList"] )

# set to global parameters
stgdict["maxTimeSteps"]=10
stgdict["checkpointEvery"]=1

stgdict["pauseToAttachDebugger"]=0

# set the dictionary again in Underworld
underworld.dictionary.SetDictionary(stgdict)

underworld.Construct()

##  lets reinit swarm guys
#   grab the TemperatureField
tfield = underworld._stgermain.GetLiveComponent("TemperatureField")

# get the number of local nodes on the temperature mesh
nLocalNodes = underworld.libUnderworld.StgDomain.Mesh_GetLocalSize( tfield.feMesh, 0 ) # 0 represents the 0th topological element of the mesh, ie the nodes

from libUnderworld import c_arrays
import math

# parameter for building the temperature profile
grad = (T_min-T_max) / float(maxR-minR)

cVal = c_arrays.DoubleArray(1)

def inCircle( xyz ) :
    """
     Return 1 if xyz location is inside circle of radius 0.5 centred at (1,1,4.5)
    """

    x = xyz[0] - 1
    y = xyz[1] - 1
    z = xyz[2] - 4.5
    if( x*x + y*y + z*z < 0.5 ):
        return 1
    else:
        return 0

for ii in range( 0, nLocalNodes ):
   result = underworld.libUnderworld.StgDomain.Mesh_GetVertex(tfield.feMesh, ii )
   ptr = c_arrays.DoubleArray.frompointer(result)
   pos = [ptr[0], ptr[1], ptr[2]] # make a position vector, assume 3D

   """
   angle from equator-equator centre
   _long = math.degrees( math.atan2( pos[0], pos[2] ))
   _lat = math.degrees( math.atan2(pos[1], pos[2]) )
   """

   # apply linear temperature
   r = math.sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2])
   temp = grad*r-grad*maxR

   # apply perturbation if position is inCircle, see func above
   if( inCircle(pos) ) :
       temp = temp+0.3

   cVal[0] = temp
   # set the temperature
   underworld.libUnderworld.StgFEM.FeVariable_SetValueAtNode( tfield, ii, cVal.cast() )

underworld.RunMainLoop()
