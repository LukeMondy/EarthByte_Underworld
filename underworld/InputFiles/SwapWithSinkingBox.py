## Rayleigh Taylor benchmark with No XML at all!                                                                              

import uwpytools
import uwpytools as uw

import uwpytools.meshing._meshing as mesh
import uwpytools.fields._fields as field
import uwpytools.swarms._swarms as swarm
import uwpytools.matrix._matrix as matrix
import uwpytools.shapes._shapes as shape
import uwpytools.rheology._rheology as rheo
import uwpytools.boundary._boundary as bc
import uwpytools.equations._equations as eq
import uwpytools.geometry._geometry as geom
import uwpytools.physics._physics as phys
import uwpytools.visualisation._visualisation as vis

uw.Init()

gdict=uw.GetCurrentPythonDictionary()

pd=uw.PrettyDictionaryPrint


uw.initDefaultParameters()
# For Rayleigh-Taylor                                                                                                           
uw.setParameters(minX=0.0,  maxX=0.9142,
                 minY=0.0,  maxY=1.0,
                 gravity = 1.0,
                 resX=48, resY=48,
                 dim=2,
                 outputPath="swapBoxViscContrastNic32",
                 maxTimeSteps     = 120,
                 saveDataEvery    = 10,
                 restartTimestep  = 120,
                 checkpointEvery = 10,
                 checkpointAppendStep=0,
                 interpolateRestart = 1,
                 particlesPerCell = 35
                 )

#help(uw.importToolBox)                                                                                                         
uw.importToolBox('Solvers')
uw.importToolBox('gLucifer')

uw.addPlugin("StgFEM_FrequentOutput")
uw.addPlugin("StgFEM_CPUTime")
uw.addPlugin("StgFEM_StandardConditionFunctions")



# Set up a standard Q1P0 Mesh                                                                                                   
[meshQ1Name,        meshP0Name,
 velFeVarName,      pressFeVarName,   # This function sets the FeVariables names on the dictionary now too. Maybe easier.                         
 gaussIntSwarmName, picIntSwarmName] = geom.meshQ1P0CartesianCreate(dim="dim", pic=True)



# The vrms plugin requires the parameters 'GaussSwarm' and 'VelocityField' to be set                                            
uw.addPlugin("Underworld_Vrms", GaussSwarm=gaussIntSwarmName, VelocityField=velFeVarName)


# Need to set up Shapes and Rheology                                                                                            
# Set up default background shape                                                                                               
shape._everywhereCreate(componentName="backgroundShape")

# Set's up the shape required for Rayleigh Taylor by default. Note that the order matters with shapes!                          
shape._belowCosinePlaneCreate(componentName="belowCosinePlaneShape")

shape._boxCreate(componentName="boxShape", startList=[0.4,0.35], endList=[0.6,0.6])
# Sets up "backgroundViscosity" by default                                                                                      
rheo._isoviscousCreate( componentName="backgroundViscosity", eta0=1.0)  # or could set eta0 in the above Parameters list        
rheo._isoviscousCreate( componentName="lightLayerViscosity", eta0=1000.0)
rheo._isoviscousCreate( componentName="boxViscosity", eta0=10000.0)
# connect shapes and rheologies                                                                                                 
rheo._joinRheologyAndShape( componentName = "background",
                            rheologyName  = "backgroundViscosity",
                            shapeName     = "backgroundShape",
                            density       = 1.0
                            )
rheo._joinRheologyAndShape( componentName = "lightLayer",
                            rheologyName  = "lightLayerViscosity",
                            shapeName     = "belowCosinePlaneShape",
                            density       = 1.5
                            )
rheo._joinRheologyAndShape( componentName = "boxLayer",
                            rheologyName  = "boxViscosity",
                            shapeName     = "boxShape",
                            density       = 2.0
                            )

uw.addPlugin("Underworld_SwapRheologies", UpperRheologyMaterial="background", LowerRheologyMaterial="lightLayer", MaterialSwarm="materialSwarm", Height=0.343)

# Need a system of Equations and Solver                                                                                         
eq.stokesSystemCreate(solver="stokesblockkspinterface", buoyancy=True, pic=True)                                                         
# The uzawa by default at the moment                                                                                            
#eq.stokesSystemCreate(buoyancy=True, pic=True)
help(eq.stokesSystemCreate)

uwpytools.solvers._multigrid.multigrid(3)

uwpytools.solvers._options.options(optionsString="-A11_ksp_type fgmres -A11_ksp_converged_reason -xxscr_ksp_view -Q22_pc_type gkgdiag -log_summary", optionsFilename="options/options-scr-mg.opt")


field._operatorFeVariableCreate(feVariableName="VelocityField", operator="Magnitude")


vis.cameraCreate(centreFieldVariable="VelocityField", name="camera")
#create some arrows for visualisation
vis.fieldArrowsCreate(name="velocityArrows", fieldVariable="VelocityField")
#vis.viewPortCreate(fieldVariable="VelocityField", camera="camera")
vis.viewPortCreate(fieldVariable="PressureField", camera="camera")
vis.databaseCreate()

vis._swarmViewerCreate(name="particleDots", pointSize=2.5, swarm="materialSwarm", variable="Density")
vis.viewPortCreate(name="particleDotsVP",scalarFieldMap="particleDots", border="border", camera="camera")
vis.viewPortCreate(fieldVariable="magnitudeVelocityField", camera="camera", arrows="velocityArrows")


vis.windowCreate(viewPortList=["PressureFieldVP magnitudeVelocityFieldVP particleDotsVP"])


# Create BC's                                                                                                                   

# Free slip on side walls                                                                                                       
bc._wallFreeSlipCreate(bcEntry="velocityMeshBCs", wall="left")
bc._wallFreeSlipCreate(bcEntry="velocityMeshBCs", wall="right")
# No slip on top and bottom                                                                                                     
bc._wallNoSlipCreate(bcEntry="velocityMeshBCs", wall="top")
bc._wallNoSlipCreate(bcEntry="velocityMeshBCs", wall="bottom")

# Let particles leave box just in case                                                                                          
uw.NewComponentEntryInStgDict( gdict, name="escapedRoutine", Type="EscapedRoutine")

#pd(gdict)                                                                                                                      
#help(swarm._integrationSwarmCreate)                                                                                            
#uw.Construct()
#uw.BuildAndInitialise()
#uw.Step(steps=12)
#uw.Finalise()

# Woo!!



#pd(gdict)

uw.Construct()
uw.BuildAndInitialise()
uw.Step(steps=10)
uw.Finalise()

