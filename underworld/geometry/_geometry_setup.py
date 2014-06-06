# Geometry package - set up higher level Meshes, Shapes etc
import underworld as _uw
# We should import the names with underscores to limit what people are exposed to when they examine this module
import underworld.meshing as _meshing
import underworld.fields as _fields
import underworld.swarms as _swarms

#import underworld.matrix as _matrix
#import underworld.boundary as _boundary
#import underworld.equations as _equations


##############################################################################
## This code adds what is required to the python dictionary 
## We eventually pass the python dictionary back to Underworld
## and Underworld then uses this information to configure and set
## itself up.
##############################################################################

'''
This code adds what is required to the python dictionary for geometric objects
Ultimately the global Dictionary gets passed back to Underworld which then actually creates the simulation

'''


def meshQ1P0CartesianCreate(resX="resX", resY="resY", resZ="resZ",
                            dim ="dim",
                            minX="minX", minY="minY", minZ="minZ",
                            maxX="maxX", maxY="maxY", maxZ="maxZ",
                            pic=True,
                            particlesPerCell="particlesPerCell",
                            withTemperature=False
                            ):
    """
    Creates a Q1-P0 mesh system.
    
    The function creates the required FeMeshes along with
    MeshVariables and FeVariables to match.
    
    Arguments like "minX" etc can be left to their default values as string parameters.
    It is also possible to set them to numerical values here.

    In the first case, matching parameters must be set at the top level of the dictionary and given
    numerical values.
    This is easily done using the 'setParameters' function in the main module.
    
    (It might be better to have the variables created in another function?)

    """
    namesDict={}  # to be returned

    globalDict = _uw.GetCurrentPythonDictionary()

    # Set minX etc to info dictionary
    globalDict["info"]["minX"]=minX
    globalDict["info"]["maxX"]=maxX
    globalDict["info"]["minY"]=minY
    globalDict["info"]["maxY"]=maxY
    globalDict["info"]["minZ"]=minZ
    globalDict["info"]["maxZ"]=maxZ
    globalDict["info"]["resX"]=resX
    globalDict["info"]["resY"]=resY
    globalDict["info"]["resZ"]=resZ
    globalDict["info"]["dim"]=dim
    globalDict["info"]["pic"]=str(pic)
    globalDict["info"]["withTemperature"]=str(withTemperature)
    globalDict["info"]["particlesPerCell"]=particlesPerCell
    # Unfortunately we currently need these on top level of Dicitionary as well.
    globalDict["minX"]=minX
    globalDict["maxX"]=maxX
    globalDict["minY"]=minY
    globalDict["maxY"]=maxY
    globalDict["minZ"]=minZ
    globalDict["maxZ"]=maxZ
    globalDict["resX"]=resX
    globalDict["resY"]=resY
    globalDict["resZ"]=resZ
    globalDict["dim"]=dim
    globalDict["particlesPerCell"]=particlesPerCell
    
    # The velocity mesh Q1 (linear)
    meshQ1Name="linearMesh"
    meshDictQ1 = _meshing.setup.cartesianMeshCreate(meshElementType="linear",
                                                    componentName=meshQ1Name,
                                                    minX=minX, minY=minY, minZ=minZ,
                                                    maxX=maxX, maxY=maxY, maxZ=maxZ,
                                                    resX=resX, resY=resY, resZ=resZ,
                                                    dim=dim
                                                    )
    # The pressure mesh P0 (constant)
    
    # This might be the correct level for this?
    globalDict["info"]["velocityMesh"]=meshQ1Name  # This is what the BC setup functions will refer to by default
    
    meshP0Name="constantMesh"
    meshDictP0 = _meshing.setup.cartesianMeshCreate(meshElementType="constant", componentName=meshP0Name, primaryMeshName=meshQ1Name)

    globalDict["info"]["pressureMesh"]=meshP0Name

    # Let's put some useful variables on the Mesh
    # Need a Velocity FeVariable for a velocity field on the mesh
    # need a mesh variable first

    velMeshVarName="velocityMeshVariable"
    velMeshVarName="velocity"

    _fields.setup.meshVariableCreate( componentName=velMeshVarName, dim=dim, rankType="Vector", dataType="Double", meshName=meshQ1Name)

    # Name the field exactly this because stuff in Underworld assumes these names
    # e.g. the Standard Condition functions.
    # But it is annoying because the capitalisation is not consistent everywhere !! 

    velName = "VelocityField"
    velFeVar = _fields.setup.feVariableCreate( componentName=velName, meshName=meshQ1Name, meshVariableName=velMeshVarName)
    globalDict["info"]["velocityField"]=velName
        
    pressMeshVarName="pressureMeshVariable"
    pressMeshVarName="pressure"
    _fields.setup.meshVariableCreate( componentName=pressMeshVarName, dim=dim, rankType="Scalar", dataType="Double", meshName=meshP0Name)
    
    pressName  = "PressureField"
    pressFeVar = _fields.setup.feVariableCreate( componentName=pressName, meshName=meshP0Name, meshVariableName=pressMeshVarName)
    globalDict["info"]["pressureField"]=pressName

    
    if withTemperature:
        tempMeshVar = "temperatureMeshVariable"   # temp lives on velocity mesh
        tempMeshVar = "temperature"   # temp lives on velocity mesh
        _fields.setup.meshVariableCreate( componentName=tempMeshVar, dim=dim, rankType="Scalar", dataType="Double", meshName=meshQ1Name)
        tempName = "TemperatureField"
        tempFeVar = _fields.setup.feVariableCreate( componentName=tempName, meshName=meshQ1Name, meshVariableName=tempMeshVar)
        namesDict["temperatureFeVariable"]=tempName
        globalDict["info"]["temperatureField"]=tempName

    # we usually need the Gauss swarm in any case

    # this should be _swarms.integrationSwarmCreate 

    swarmType="Gauss"
    [gaussIntSwarm, integratorName, emptymapper] = _swarms.setup.integrationSwarmCreate(swarmType=swarmType, meshName=meshQ1Name)
    globalDict["info"]["gaussIntSwarm"]=gaussIntSwarm["name"]

    picIntSwarm={}
    if pic:
        swarmType="PIC"
        # need a material Swarm first if swarm is type PIC (let it create one with a default name)

        swarmDict = _swarms.setup.materialSwarmCreate(meshName=meshQ1Name, particlesPerCell=particlesPerCell)

        # This function creates what it needs for a functioning integration swarm by default

        [picIntSwarm, integratorName, mapper] = _swarms.setup.integrationSwarmCreate( materialSwarmName=swarmDict["name"], 
                                                                                      swarmType=swarmType,
                                                                                      integratorName = integratorName,
                                                                                      meshName=meshQ1Name)

        globalDict["info"]["picIntSwarm"]=picIntSwarm["name"]

        # Need an advector for the swarm
        _swarms.setup.materialSwarmAdvectorCreate( integratorName    = integratorName,
                                                   swarmName         = swarmDict["name"],
                                                   velocityFieldName = velName)

    # Set FeVariable names on the dictionary at top level.
    # Other functions can then access these names so users don't have to think about this stuff.
    # BUT NOT IN THE ROOT DICTIONARY, SURELY !!!!! 

    #globalDict["velocityFeVariable"]= velName
    #globalDict["pressureFeVariable"]= pressName
    #globalDict["gaussIntSwarm"]     = gaussIntSwarm["name"]
    #globalDict["picIntSwarm"]       = picIntSwarm["name"]
    globalDict["velocityMesh"]       = meshQ1Name
    globalDict["temperatureMesh"]       = meshQ1Name
    globalDict["pressureMesh"]       = meshP0Name

    # It is easier to just return some names here (... but we should actually be consistent)
    namesDict["linearMesh"]=meshQ1Name
    namesDict["constantMesh"]=meshP0Name
    namesDict["velocityField"]=velName
    namesDict["pressureField"]=pressName
    namesDict["gaussIntSwarm"]     =gaussIntSwarm["name"]
    namesDict["picIntSwarm"]       =picIntSwarm["name"]
    return namesDict

