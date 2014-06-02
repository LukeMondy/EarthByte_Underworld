
# coding: utf-8

# In[1]:

## Lid Driven

import uwpytools
import uwpytools as uw

import uwpytools.meshing  as meshing
import uwpytools.geometry as geometry
import uwpytools.fields   as fields
import uwpytools.swarms   as swarms
import uwpytools.material as material
import uwpytools.visualisation as visual

import uwpytools.shapes as shape
import uwpytools.rheology as rheo
import uwpytools.boundary as bc
import uwpytools.equations as eq
import uwpytools.physics as phys

uw.Init()

gdict=uw.GetCurrentPythonDictionary()

uw.PrettyDictionaryPrint(gdict)


# In[2]:

minX=-2.0
maxX=2.0
minY=0.0
maxY=2.0
dim=2
resX=64
resY=32
particlesPerCell = 35
steps=3

uw.initDefaultParameters()
uw.setParameters(outputPath="liddrivenSinePIC"+str(resX), 
                 maxTimeSteps     = steps
                 )


# In[3]:


uw.importToolBox('gLucifer')


# In[4]:

# Set up a standard Q1P0 Mesh
geoNames=geometry.setup.meshQ1P0CartesianCreate(dim=dim, pic=True, 
                                                minX=minX, maxX=maxX, minY=minY, maxY=maxY,
                                                resX=resX, resY=resY,
                                                particlesPerCell=particlesPerCell)


# In[5]:

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


# In[6]:

# Need a system of Equations and Solver                                                                                         
#eq.stokesSystemCreate(solver="stokesblockkspinterface", pic=True)                                                         
# The uzawa by default at the moment                                                                                            
eq.setup.stokesSystemCreate(pic=True)
#help(eq.setup.stokesSystemCreate)



# In[7]:


# Create BC's

# Free slip on side walls
bc.setup.wallFreeSlipCreate(wall="left")
bc.setup.wallFreeSlipCreate(wall="right")

# Set function on top boundary
bc.setup.wallSetFuncCreate(wall="top", func="Velocity_SinusoidalLid")

# No slip on bottom
bc.setup.wallNoSlipCreate(wall="bottom")


# In[8]:

uw.importToolBox('gLucifer')
            
# Create a field to visualise.
fields.setup.operatorFeVariableCreate(feVariableName="VelocityField", operator="Magnitude")

visual.setup.cameraCreate(centreFieldVariable="VelocityField", name="camera")
#create some arrows for visualisation
visual.setup.fieldArrowsCreate(name="velocityArrows", fieldVariable="VelocityField")
#vis.viewPortCreate(fieldVariable="VelocityField", camera="camera")
visual.setup.viewPortCreate(fieldVariable="PressureField", camera="camera")
visual.setup.viewPortCreate(fieldVariable="magnitudeVelocityField", camera="camera", arrows="velocityArrows")

visual.setup.databaseCreate()


visual.setup.windowCreate(viewPortList=["PressureFieldVP magnitudeVelocityFieldVP"])


# In[9]:

# Let particles leave box just in case
uw.NewComponentEntryInStgDict( gdict, name="escapedRoutine", Type="EscapedRoutine")

uw.Construct()
uw.BuildAndInitialise()
uw.Step(steps=steps)
uw.Finalise()

# Woo!!


# In[ ]:




# In[ ]:



