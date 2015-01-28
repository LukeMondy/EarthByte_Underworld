
# coding: utf-8

# In[1]:

## Rayleigh Taylor benchmark with No XML at all!                                                                              

import underworld as uw

import underworld.meshing  as meshing
import underworld.geometry as geometry
import underworld.fields   as fields
import underworld.swarms   as swarms
import underworld.material as material
import underworld.visualisation as visual
import underworld.shapes as shape
import underworld.rheology as rheo
import underworld.boundary as bc
import underworld.equations as eq
import underworld.physics as phys


# In[2]:

gdict=uw.dictionary.GetDictionary()
uw.dictionary.initDefaultParameters()


# In[3]:

# For Rayleigh-Taylor
minX=0.0
maxX=0.9142
minY=0.0
maxY=1.0

resX=32
resY=32
steps=10
particlesPerCell=35

dim = 2

uw.dictionary.setParameters(gravity = 1.0,
                 outputPath="raytay",
                 #restartTimestep  = 59,
                 maxTimeSteps     = steps,
                 resX = 32,
                 resY = 32
                 )
            


# In[4]:


# Set up a standard Q1P0 Mesh
geoNames=geometry.setup.meshQ1P0CartesianCreate(dim=dim, pic=True, 
                                                minX=minX, maxX=maxX, minY=minY, maxY=maxY,
                                                resX=resX, resY=resY,
                                                particlesPerCell=particlesPerCell)


# In[6]:

# Need to set up Shapes and Rheology                                                                                            
# Set up default background shape                                                                                               
shape.setup.everywhereCreate(componentName="backgroundShape")

# Set's up the shape required for Rayleigh Taylor by default. Note that the order matters with shapes!                          
shape.setup.belowCosinePlaneCreate(componentName="belowCosinePlaneShape")

# Sets up "backgroundViscosity" by default                                                                                      
rheo.setup.isoviscousCreate( componentName="backgroundViscosity", eta0=1.0)  # or could set eta0 in the above Parameters list        
rheo.setup.isoviscousCreate( componentName="lightLayerViscosity", eta0=1.0)
# connect shapes and rheologies   

material.setup.materialCreate( componentName = "background",
                               rheologyName  = "backgroundViscosity",
                               shapeName     = "backgroundShape",
                               density       = 1.0
                               )
material.setup.materialCreate( componentName = "lightLayer",
                               rheologyName  = "lightLayerViscosity",
                               shapeName     = "belowCosinePlaneShape",
                               density       = 0.0
                               )


# In[7]:

# Need a system of Equations and Solver                                                                                         
#eq.stokesSystemCreate(solver="stokesblockkspinterface", buoyancy=True, pic=True)                                                         
# The uzawa by default at the moment                                                                                            
eq.setup.stokesSystemCreate(buoyancy=True, pic=True)
#help(eq.setup.stokesSystemCreate)



# In[8]:

# Create BC's                                                                                                                   

# Free slip on side walls                                                                                                       
bc.setup.wallFreeSlipCreate( wall="left")
bc.setup.wallFreeSlipCreate( wall="right")
# No slip on top and bottom                                                                                                     
bc.setup.wallNoSlipCreate( wall="top")
bc.setup.wallNoSlipCreate( wall="bottom")

# Let particles leave box just in case                                                                                          
uw.dictionary.UpdateDictWithComponent( gdict, name="escapedRoutine", Type="EscapedRoutine")
fields.setup.operatorFeVariableCreate(feVariableName="VelocityField", operator="Magnitude")
uw.Construct()
uw.RunMainLoop()


# In[5]:

# Create a field to visualise.

#import glucifer.pylab as plt


#visual.setup.cameraCreate(centreFieldVariable="VelocityField", name="camera")
# create some arrows for visualisation
#visual.setup.fieldArrowsCreate(name="velocityArrows", fieldVariable="VelocityField")
#visual.setup.viewPortCreate(fieldVariable="PressureField", camera="camera")
#visual.setup.databaseCreate()

# todo: get swarm etc from outputs from above
#visual.setup.swarmViewerCreate(name="particleDots", pointSize=2.5, swarm="materialSwarm", variable="Density")
#visual.setup.viewPortCreate(name="particleDotsVP",scalarFieldMap="particleDots", border="border", camera="camera")
#visual.setup.viewPortCreate(fieldVariable="magnitudeVelocityField", camera="camera", arrows="velocityArrows")


#visual.setup.windowCreate(viewPortList=["PressureFieldVP magnitudeVelocityFieldVP particleDotsVP"])



# In[6]:

# In[7]:

#fig = plt.figure(num="need3g2y")
#fig.Points(swarm="materialSwarm", colourVariable="materialSwarm-MaterialIndex", pointSize=3.)


# In[8]:

#fig.show()


# In[9]:

#fig.saveDB("somenewDB.gldb")


# In[9]:




# In[9]:



