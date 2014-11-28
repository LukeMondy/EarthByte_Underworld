# Solvers package - activates and configures the "new solver system"
import underworld as _uw
##############################################################################
# This code adds what is required to the python dictionary
# to set up a Mesh for Underworld.
# We eventually pass the python dictionary back to Underworld
# and Underworld then uses this information to configure and set
# itself up.
##############################################################################

'''
This code adds what is required to the python dictionary for a Mesh
Ultimately the global Dictionary gets passed back to Underworld which then actually creates the simulation

'''


def cartesianMeshCreate( componentName="",
                         resX=4, resY=4, resZ=1,
                         minX=0, minY=0, minZ=0,   # Maybe can use minCoord etc here instead of these. But this way defaults are easier?
                         maxX=1, maxY=1, maxZ=1,   # Should be able to pass in a list that has only 2 members in 2d. Will test this to make sure
                         dim=2,
                         meshElementType="linear",
                         primaryMeshName="linearMesh",
                         nameSpace="stokesEqn"):
    """
    Set up a Mesh in the dictionary.
    Possible meshElementTypes are ( constant, linear, quadratic )

    Args:
        componentName (string)         : The name of the Mesh (optional: will create a name based on meshElementType if null)
        resX (int)
        resY (int)
        resZ (int)
        minX (double)
        minY (double)
        minZ (double)
        maxX (double)
        maxY (double)
        maxZ (double)
        dim  (int)
        meshElementType (string)       : Must be one of ( constant, linear, quadratic ). Defaults to linear
        primaryMeshName (string)       : Only needed for the 'constant' Mesh.
    """
    meshAllowedTypes = ["linear", "constant", "quadratic"]

    if meshElementType not in meshAllowedTypes:
        print "Error: Need to choose one of  ( constant, linear, quadratic ) for meshElementType"
        return

    globalDict = _uw.dictionary.GetDictionary()

    if componentName == "":  # Then create a name automatically
        componentName = meshElementType + "Mesh"    # e.g. linearMesh
    baseMeshType = "FeMesh"                       # All meshes in Underworld are this basic type AFAIKATMIT

    # Want to test if we already have a Mesh of the name we are going to use here.
    componentName  = _uw.utils.uniqueComponentNameGlobalDict(componentName)

    # At this point we should have ensured a new unique name for our new Mesh component
    # e.g. Now we can have more than one linear cartesian mesh
    newComponentMeshDict = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                                   name=componentName,
                                                                   Type=baseMeshType,    # currently always "FeMesh"
                                                                   elementType=meshElementType
                                                                   )
    ################################################################
    # add this mesh name to the 'info' sub-dictionary for other components to find  #
    ################################################################
    # globalDict["info"][meshElementType+"Mesh"]="componentName"

    ############################################################################################
    # Now we need a Mesh Generator (this should already be a unique name.
    # If it's not then bad luck)
    ############################################################################################

    thisMeshName = componentName
    componentName = componentName + "Generator"

    ############################################################################################
    # It appears that the constant Mesh needs to have a so-called "elementMesh"
    # which is in fact the mesh used for any swarms, which is always the Velocity Mesh.
    # The current implementation of the constant Mesh gets its dimensions from this other Mesh
    ############################################################################################

    if( meshElementType == "constant" ):
        if primaryMeshName not in globalDict["components"]:
            print "Error: Need to create a primary Mesh first for the constant Mesh"
            print "       This is almost always the mesh used for the Velocity in the model e.g. \"linearMesh\""
            print "       The constant Mesh gets it's dimensions from this \"primary\" Mesh"
            return

        newComponentMeshGenDict = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                                          name=componentName,
                                                                          Type="C0Generator",  # This is the generator function for a constant Mesh
                                                                          mesh=thisMeshName,
                                                                          elementMesh=primaryMeshName  # The constant Mesh gets it's dimensions from primaryMesh
                                                                          )
        return [newComponentMeshDict, newComponentMeshGenDict]

    # Create the Generator component for this Mesh
    if( meshElementType == "linear"):
        genType = "CartesianGenerator"
    if( meshElementType == "quadratic"):
        genType = "C2Generator"

    minCoordList = [minX, minY, minZ]
    maxCoordList = [maxX, maxY, maxZ]
    sizeList     = [resX, resY, resZ]

    newComponentMeshGenDict = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                                      name     = componentName,
                                                                      Type     = genType,  # This is the generator function for a linear Mesh
                                                                      mesh     = thisMeshName,
                                                                      maxCoord = maxCoordList,  # passing in lists here
                                                                      minCoord = minCoordList,
                                                                      size     = sizeList,
                                                                      dims = dim,
                                                                      shadowDepth = 1  # Hard wired for now...it always seems to be set to 1 anyway
                                                                      )

    return [newComponentMeshDict, newComponentMeshGenDict]


#       newComponentMeshDict = dict()
#       newComponentMeshDict["minCoord"] = minCoordList
#       newComponentMeshDict["maxCoord"] = maxCoordList
#       newComponentMeshDict["size"]     = sizeList

# newComponentMeshDict["name"] = componentName  # useful to have this ...
#    newComponentMeshDict["mesh"] = thisMeshName
#    newComponentMeshDict["Type"] = "CartesianGenerator"
#    newComponentMeshDict["dims"] = dim
# newComponentMeshDict["shadowDepth"] = 1        # Hard wired for now...it always seems to be set to 1 anyway
# globalDict["components"][componentName] = newComponentMeshDict   # The order of operations actually doesn't matter here
