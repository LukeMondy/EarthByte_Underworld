
# coding: utf-8

# In[1]:

## Thermal Convection                                                                             

from pprint import pprint

import uwpytools as uw

import uwpytools.meshing  as meshing
import uwpytools.geometry as geometry
import uwpytools.fields   as fields
import uwpytools.swarms   as swarms
import uwpytools.material as material
import uwpytools.visualisation as visual
# import uwpytools.fields._fields as fields

# import uwpytools.swarms._swarms as swarms
# import uwpytools.matrix._matrix as matrix
import uwpytools.shapes as shape
import uwpytools.rheology as rheo
import uwpytools.boundary as bc
import uwpytools.equations as eq
import uwpytools.physics as phys

uw.Init()

gdict=uw.GetCurrentPythonDictionary()

uw.PrettyDictionaryPrint(gdict)


# In[2]:

# These defaults should be split up either by problem or by module (fields etc)

uw.initDefaultParameters()

minX=0.0
maxX=1.0
minY=0.0
maxY=1.0

resX=64
resY=64
steps=6
particlesPerCell=35

dim = 2

uw.setParameters(gravity = 1.0,
                 outputPath="thermalConvection"+str(resX),
                 #restartTimestep  = 59,
                 maxTimeSteps     = steps,
                 )


# In[3]:


# Set up a standard Q1P0 Mesh (with temperature field)
geoNames=geometry.setup.meshQ1P0CartesianCreate(dim=dim, pic=True, 
                                                minX=minX, maxX=maxX, minY=minY, maxY=maxY,
                                                resX=resX, resY=resY,
                                                particlesPerCell=particlesPerCell, withTemperature=True)


# In[4]:

# Need to set up Shapes and Rheology
# Set up default background shape
shape.setup.everywhereCreate(componentName="backgroundShape")

# Sets up "backgroundViscosity" by default
rheo.setup.isoviscousCreate( componentName="backgroundViscosity", eta0=1.0)  # or could set eta0 in the above Parameters list

# connect shapes and rheologies
material.setup.materialCreate( componentName = "background",
                               rheologyName  = "backgroundViscosity",
                               shapeName     = "backgroundShape",
                               density       = 1.0
                               )


# In[5]:

# Need a system of Equations and Solver                                                                                         
#eq.stokesSystemCreate(solver="stokesblockkspinterface", buoyancy=True, pic=True)                                                         
# The uzawa by default at the moment                                                                                            
eq.setup.stokesSystemCreate(buoyancy=True, pic=True, buoyancyType="thermal")
#help(eq.setup.stokesSystemCreate)



# In[6]:

eq.setup.advectionDiffusionEquationCreate(equationName="energyEqn",                                                        
                                          diffusivity=1.0 )



# In[7]:

# Create BC's                                                                                                                   

# Free slip on side walls
bc.setup.wallFreeSlipCreate( wall="left")
bc.setup.wallFreeSlipCreate( wall="right")

# Free slip on top
bc.setup.wallFreeSlipCreate( wall="top")

# Free slip on bottom
bc.setup.wallFreeSlipCreate( wall="bottom")

bc.setup.wallTemperatureCreate(wall="bottom", value=1.0)
bc.setup.wallTemperatureCreate(wall="top", value=0.0)

bc.setup.temperatureICSinusoidalCreate(TopLayerCoord=maxY, BottomLayerCoord=minY, BottomLayerBC=1.0)

# Let particles leave box just in case                                                                                          
uw.NewComponentEntryInStgDict( gdict, name="escapedRoutine", Type="EscapedRoutine")

#pd(gdict)                                                                                                                      
#help(swarm._integrationSwarmCreate)     


# In[9]:

uw.importToolBox('gLucifer')


fields.setup.operatorFeVariableCreate(feVariableName="TemperatureField", operator="Gradient", componentName="TemperatureGradientsField")
# Create a field to visualise.
fields.setup.operatorFeVariableCreate(feVariableName="VelocityField", operator="Magnitude", componentName="VelocityMagnitudeField")
fields.setup.operatorFeVariableCreate(feVariableName="VelocityField", operator="Gradient", componentName="VelocityGradientsField")
fields.setup.operatorFeVariableCreate(feVariableName="VelocityGradientsField", operator="TensorInvariant", componentName="VelocityGradientsInvariantField")

uw.addCheckPointVariables(["TemperatureGradientsField", "VelocityGradientsField", "VelocityGradientsInvariantField"])

visual.setup.cameraCreate(centreFieldVariable="VelocityField", name="camera")
#create some arrows for visualisation
visual.setup.fieldArrowsCreate(name="velocityArrows", fieldVariable="VelocityField")
#vis.viewPortCreate(fieldVariable="VelocityField", camera="camera")
visual.setup.viewPortCreate(fieldVariable="PressureField", camera="camera")
visual.setup.viewPortCreate(fieldVariable="VelocityMagnitudeField", camera="camera", arrows="velocityArrows")
visual.setup.viewPortCreate(fieldVariable="TemperatureField", camera="camera")

visual.setup.databaseCreate()


visual.setup.windowCreate(viewPortList=["TemperatureFieldVP VelocityMagnitudeFieldVP PressureFieldVP"])



# In[10]:

uw.Construct()
uw.BuildAndInitialise()
uw.Step(steps=steps)
uw.Finalise()

# Woo!!


# In[11]:

uw.getInfo()


# In[8]:



