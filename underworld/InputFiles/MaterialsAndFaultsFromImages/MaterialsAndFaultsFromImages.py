# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import uwpytools
uwpytools.InitWithArgs("BasicRT.xml")

# <codecell>

uwdict = uwpytools.GetCurrentDictionary()
# lets setup domain
uwdict["minX"] = 0.
uwdict["maxX"] = 4.
uwdict["minY"] = 0.
uwdict["maxY"] = 1.
# also model run parameters
uwdict["maxTimeSteps"] = -1
uwdict["checkpointEvery"] = 1

# <codecell>

# lets create materials.  we'll just use default values for now, and set everything later
uwdict["components"]["layer1"]={ "Type":"RheologyMaterial", "Shape":"backgroundShape", "density":1., "Rheology":"backgroundViscosity"}
uwdict["components"]["layer2"]={ "Type":"RheologyMaterial", "Shape":"backgroundShape", "density":1., "Rheology":"backgroundViscosity"}
uwdict["components"]["layer3"]={ "Type":"RheologyMaterial", "Shape":"backgroundShape", "density":1., "Rheology":"backgroundViscosity"}
uwdict["components"]["layer4"]={ "Type":"RheologyMaterial", "Shape":"backgroundShape", "density":1., "Rheology":"backgroundViscosity"}
uwdict["components"]["layer5"]={ "Type":"RheologyMaterial", "Shape":"backgroundShape", "density":1., "Rheology":"backgroundViscosity"}
uwdict["components"]["air"]   ={ "Type":"RheologyMaterial", "Shape":"backgroundShape", "density":1., "Rheology":"backgroundViscosity"}

# <codecell>

# ok, lets create some RBF guys.  Define a function
def NewRBFSwarmDictionary( rbfname, supportRadius ):
    # manager dict
    managerDict   = { "Type":"RBFManager", "ParticleSupportRadius":supportRadius, "RBFdim":2, "RBFSwarm":rbfname+"Swarm" }
    swarmDict     = { "Type":"Swarm"     , "CellLayout":rbfname+"CLLLayout" }
    layoutDict    = { "Type":"CLLCellLayout",  "GeometryMesh":"linearMesh", "CellSize":supportRadius, "MeshDim":2 }
    variableDict  = { "Type":"Variable",  "Rank":"Scalar", "DataType":"Double" }
    swarmVarDict  = { "Type":"SwarmVariable",  "Swarm":rbfname+"Swarm", "Variable":rbfname+"Variable" }
    return { rbfname+"Manager"  : managerDict, 
             rbfname+"Swarm"    : swarmDict, 
             rbfname+"CLLLayout": layoutDict,
             rbfname+"Variable" : variableDict,
             rbfname+"SwarmVar" : swarmVarDict }

# generate a support radius.  lets try for 10 particles per element width
supportRadius = 2.*(float(uwdict["maxX"]) - float(uwdict["minX"]))/float(uwdict["elementResI"]*10)
uwdict["components"].update( NewRBFSwarmDictionary( "faults1", supportRadius ) )
uwdict["components"].update( NewRBFSwarmDictionary( "faults2", supportRadius ) )
uwdict["components"].update( NewRBFSwarmDictionary( "faults3", supportRadius ) )

# <codecell>

uwpytools.SetDictionary(uwdict) 

# <codecell>

uwpytools.Construct()

# <codecell>

uwpytools.BuildAndInitialise()

# <codecell>

# fire up the materials image... some massaging (via Gimp) was requried to get it to a usable state
%matplotlib inline
import matplotlib.pylab as pylab
pylab.rcParams['figure.figsize'] = 16, 12
import matplotlib.pyplot as plt
from scipy import misc
matImg = misc.imread('faults2int_mats.png')
cmap = plt.get_cmap('terrain', 6)
plt.imshow(matImg, cmap=cmap)
plt.colorbar()

# <codecell>

# ok, lets build a map from values found within the image, to the material index values stgermain has set
context = uwpytools.GetLiveComponent( "context" )
layer1 = uwpytools.PICellerator.Materials_Register_GetByName( context.materials_Register, "layer1" )
layer2 = uwpytools.PICellerator.Materials_Register_GetByName( context.materials_Register, "layer2" )
layer3 = uwpytools.PICellerator.Materials_Register_GetByName( context.materials_Register, "layer3" )
layer4 = uwpytools.PICellerator.Materials_Register_GetByName( context.materials_Register, "layer4" )
layer5 = uwpytools.PICellerator.Materials_Register_GetByName( context.materials_Register, "layer5" )
air = uwpytools.PICellerator.Materials_Register_GetByName( context.materials_Register, "air" )
# by inspection of above image, create map
mapDictionary = { 4: layer1.index, 2:layer2.index, 0:layer3.index, 3:layer4.index, 1:layer5.index, 5:air.index }

# <codecell>

# lets define a function to handle the image
def ImageAsField( numpyImage, minTuple, maxTuple, queryLocationTuple ):
    # lets make the image domain a touch bigger to make life easier
    epsX = (maxTuple[0]- minTuple[0])/100000.
    epsY = (maxTuple[1]- minTuple[1])/100000.
    numPixX = float(numpyImage.shape[1])
    numPixY = float(numpyImage.shape[0])
    
    qPixX = int(numPixX*(queryLocationTuple[0] - minTuple[0]+epsX)/(maxTuple[0]-minTuple[0]+2*epsX)) 
    qPixY = int(numPixY) - int(numPixY*(queryLocationTuple[1] - minTuple[1]+epsY)/(maxTuple[1]-minTuple[1]+2*epsY)) -1
    
    return numpyImage[qPixY][qPixX]

# <codecell>

swarm = uwpytools.GetLiveComponent("materialSwarm")
uwpytools.Swarm_PrintVariables(swarm)
variables = uwpytools.Swarm_GetVariables(swarm)

# <codecell>

matIndexVariable = variables[4][1]
posVariable = variables[3][1]

# <codecell>

# ok, lets sweep image
# setup img domain as follows
imgMinDom = (0.,0.)
imgMaxDom = (4.252,1.042)

for ii in range(0,swarm.particleLocalCount):
	# grab the position
    pos = uwpytools.SwarmVariable_GetValueAt( posVariable, ii )
    imgVal = ImageAsField( matImg, imgMinDom, imgMaxDom, tuple(pos) )
    uwpytools.SwarmVariable_SetValueAt( matIndexVariable, ii, [mapDictionary[imgVal]] ) 

# <codecell>

faults1Img = misc.imread('faults2int_faults1.png')
plt.imshow(faults1Img, )

# <codecell>

faults2Img = misc.imread('faults2int_faults2.png')
plt.imshow(faults2Img, )

# <codecell>

faults3Img = misc.imread('faults2int_faults3.png')
plt.imshow(faults3Img, )

# <codecell>

# ok, lets add some particles along fault lines to RBF swarm
numpartsX = uwdict["elementResI"]*10
numpartsY = int(float(uwdict["elementResJ"]*10./4.))

minX = float(uwdict["minX"])
maxX = float(uwdict["maxX"])
minY = float(uwdict["minY"])
maxY = float(uwdict["maxY"])

f1swarm = uwpytools.GetLiveComponent("faults1Swarm")
f2swarm = uwpytools.GetLiveComponent("faults2Swarm")
f3swarm = uwpytools.GetLiveComponent("faults3Swarm")

for ii in range(0,numpartsX):
    for jj in range(0,numpartsY):
        partX = minX + (float(ii)+0.5)*(maxX-minX)/float(numpartsX)
        partY = minY + (float(jj)+0.5)*(maxY-minY)/float(numpartsY)
        f1Val = ImageAsField( faults1Img, imgMinDom, imgMaxDom, (partX,partY) )
        f2Val = ImageAsField( faults2Img, imgMinDom, imgMaxDom, (partX,partY) )
        f3Val = ImageAsField( faults3Img, imgMinDom, imgMaxDom, (partX,partY) )
        if f1Val != 0 : uwpytools.PICellerator.GeneralSwarm_AddParticle( f1swarm, 2, partX, partY, 0.)
        if f2Val != 0 : uwpytools.PICellerator.GeneralSwarm_AddParticle( f2swarm, 2, partX, partY, 0.)
        if f3Val != 0 : uwpytools.PICellerator.GeneralSwarm_AddParticle( f3swarm, 2, partX, partY, 0.)

# <codecell>

# need to call init func to calc rbf densities
uwpytools.StGermain.Stg_Component_Initialise(uwpytools.GetLiveComponent("faults1Manager"), context, True)
uwpytools.StGermain.Stg_Component_Initialise(uwpytools.GetLiveComponent("faults2Manager"), context, True)
uwpytools.StGermain.Stg_Component_Initialise(uwpytools.GetLiveComponent("faults3Manager"), context, True)

# <codecell>

uwpytools.Swarm_PrintVariables(f1swarm)
variables = uwpytools.Swarm_GetVariables(f1swarm)
variables[5][1].isCheckpointedAndReloaded = True

# <codecell>

uwpytools.RunMainLoop()
uwpytools.Finalise()

