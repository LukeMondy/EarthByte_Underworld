# Fields package - configures Fields and Variables on a Mesh
import underworld as _uw
##############################################################################
# This code adds what is required to the python dictionary
# to set up Fields and Variables for Underworld.
# We eventually pass the python dictionary back to Underworld
# and Underworld then uses this information to configure and set
# itself up.
##############################################################################

'''
This code adds what is required to the python dictionary for feVariables and Variables
Ultimately the global Dictionary gets passed back to Underworld which then actually creates the simulation

'''

########################################################################################################################


def feVariableCreate( componentName="",
                      dim=2,
                      units="",
                      meshName="",
                      meshVariableName="",    # Needed by DofLayout component
                      removeBCs="True"):
    """
    Create a new FeVariable on a mesh.
    feVariable variables provide interpolation of nodal values anywhere based on the nodes on a Mesh.
    feVariable variables are required input for the Matrices and Vectors in the system.

    Args:
       meshName         :  e.g. velocityMesh
       meshVariableName :  This is needed for the DofLayout component that gets created as well.
       removeBCs        :  True or False to create matrices with BCs in or left out.
    """

    if meshName == "":
        _uw.utils.sendError("meshName=\"someMesh\" is a required argument")
        print "e.g. use cartesianMeshCreate from the meshing module to create one"
        return -1

    if meshVariableName == "":
        _uw.utils.sendError("meshVariableName=\"someMeshVariable\" is a required argument")
        print "e.g. use meshVariableCreate from the this module to create one"
        return -1

    globalDict = _uw.dictionary.GetDictionary()
    _uw.utils.warnMissingComponent(globalDict, meshName )
    _uw.utils.warnMissingComponent(globalDict, meshVariableName )

    if componentName == "":
        mName = meshName[0].upper() + meshName[1:]  # capitalize only first letter and leave others as they were.
        componentName = "feVariable" + mName     # e.g. feVariableVelocityMesh

    # Want to test if we already have an operator Fe Variable of the name we are going to use here.

    componentName = _uw.utils.checkForNewComponentName(globalDict, componentName)

    # At this point we should have ensured a new unique name for the new FeVariable component
    # Create required components for this variable

    dofLayoutName = meshName + "DofLayout"
    dofLayoutName = _uw.utils.checkForNewComponentName(globalDict, dofLayoutName)
    newDofLayout = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                           name = dofLayoutName,
                                                           Type = "DofLayout",
                                                           MeshVariable = meshVariableName
                                                           )
    # It seems we always need some BCs and ICs structs when we create a Fe Variable
    # We then need to create a struct of the same name outside of the "components" dictionary
    # at the top level dictionary to set BCs etc. This struct must have the same name as the BC structs
    # in the component dictionary <-- weird

    newBCname = _uw.utils.checkForNewComponentName(globalDict, meshName + "BCs")
    newBC = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                    name = newBCname,
                                                    Type = "CompositeVC",
                                                    Data = meshName
                                                    )
    newICname = _uw.utils.checkForNewComponentName(globalDict, meshName + "ICs")
    newIC = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                    name = newICname,
                                                    Type = "CompositeVC",
                                                    Data = meshName
                                                    )
    newComponentFeVarDict = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                                    name      = componentName,
                                                                    Type      = "FeVariable",
                                                                    FEMesh    = meshName,
                                                                    DofLayout = dofLayoutName,
                                                                    BC        = newBCname,
                                                                    IC        = newICname,
                                                                    LinkedDofInfo = meshName + "LinkedDofs",  # dummy name? appears to be not usually used.
                                                                    #removeBCs = removeBCs,
                                                                    outputUnits = units
                                                                    )

    # need to get correct BCs name for BC setup
    # Field variable name better to refer to for BCs rather then the mesh name?
    # e.g. we might have velocity and temperature but both live on the velocity mesh
    globalDict["info"][componentName + "BCs"] = newBCname
    globalDict["info"][componentName + "ICs"] = newICname

    return [newBCname, newComponentFeVarDict, newDofLayout, newBC, newIC]


########################################################################################################################


def operatorFeVariableCreate( componentName="",
                              dim=2,
                              units="",
                              feVariableName="",
                              operator=""):
    # need a list of valid units here
    """
    Args:
        feVariableName  :  input fe variable for new variable
        operator        :  operator that takes feVariableName as an argument. See list of Operators below
        units           :  optional. e.g. "cm/yr"
    Operators:
        TensorInvariant     :      Is in fact the sqrt of the 2nd Invariant; I_2 = \sqrt{ 0.5 J_{ij} J_{ij} }
        Magnitude           :
        Gradient            :
        TakeFirstComponent  :
        TakeSecondComponent :
        TensorSymmetricPart :
        TensorAntisymmetricPart   :
        SymmetricTensor_Invariant :
        Tensor_AverageTrace       :
        (These are defined in file StgDomain/Utils/src/Operator.c for more)
    """

    # Might just give a warning if not in this list. That way don't have to panic about always being up to date.

    opList = ["Magnitude", "Gradient", "TensorInvariant", "TakeFirstComponent", "TakeSecondComponent",
              "TensorSymmetricPart", "TensorAntisymmetricPart", "SymmetricTensor_Invariant", "Tensor_AverageTrace"]

    if operator == "":
        _uw.utils.sendError("operator=\"someOperator\" is a required argument")
        print "e.g. One of:"
        print opList
        return -1

    if feVariableName == "":
        _uw.utils.sendError("feVariableName=\"someFeVariable\" is a required argument")
        print "e.g. use feVariableCreate from the this module to create one"
        return -1

    if operator not in opList:
        _uw.utils.sendWarning("Operator " + operator + " not in list of known operators.")

    globalDict = _uw.dictionary.GetDictionary()
    _uw.utils.warnMissingComponent(globalDict, feVariableName )

    if componentName == "":
        fName = feVariableName[0].upper() + feVariableName[1:]  # capitalize only first letter and leave others as they were.
        opName = operator[0].lower() + operator[1:]   # lower first letter only (making a convention that all structs have names with lower-case 1st letter)
        componentName = opName + fName     # e.g. gradientVelocityFeVariable

    # Want to test if we already have an operator Fe Variable of the name we are going to use here.

    componentName = _uw.utils.checkForNewComponentName(globalDict, componentName)

    # At this point we should have ensured a new unique name for the new Fe Variable component
    newComponentFeVarDict = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                                    name     = componentName,
                                                                    Type     = "OperatorFeVariable",
                                                                    Operator = operator,
                                                                    Operand  = feVariableName,
                                                                    outputUnits = units
                                                                    )
    return newComponentFeVarDict

########################################################################################################################


def meshVariableCreate( componentName="",
                        dim=2,
                        rankType="Vector",
                        dataType="Double",
                        meshName="velocityMesh",
                        mergeType="merge"):
    """
    Set up a Mesh Variable in the dictionary.
    Possible rankTypes are ( Scalar, Vector )

    Needed by DofLayouts which are in turn needed by FeVariables (e.g. VelocityFeVariable)

    Args:
        componentName (string)         : The name of the mesh Variable.
        rankType  (string)             : type of Variable ( Scalar, Vector )
        dataType  (string)             : data type (Int, Float, Double, Possoms)
        dim  (int)                     :
        meshName (string)              : Mesh which this Variable is defined on
        mergeType (string)             : If this is "replace" then then any existing variable of the same name will be clobbered

    Returns:
        Dictionary node for the component definition

    Behaviours:
        Will define a default name of meshVector or meshScalar and will ignore mergeType if this is the case

    """

    rankAllowedTypes = ["Scalar", "Vector"]

    if (rankType not in rankAllowedTypes):
        _uw.utils.sendError("Need to choose one of  ( Scalar, Vector ) for rankType")
        return -1

    if componentName == "":
        mergeType = "merge"                  # Should not over-write auto-generated variables.
        componentName = "mesh" + rankType      # e.g. meshVector

    # Want to test if we already have a mesh Variable of the name we are going to use here.

    if mergeType != "replace":
        componentName = _uw.utils.uniqueComponentNameGlobalDict(componentName)

    # At this point we should have ensured a new unique name for the new mesh Variable component

    # Will these vx,vy,vz sub-names be unique global names in the dictionary ?

    if rankType == "Vector":
        namesList = ["vx", "vy", "vz"]
        newComponentMeshVarDict = _uw.dictionary.UpdateDictWithComponent(   name=componentName,
                                                                            Type="MeshVariable",
                                                                            Rank=rankType,
                                                                            DataType=dataType,
                                                                            VectorComponentCount=dim,
                                                                            names=namesList,
                                                                            mesh =meshName
                                                                            )

    if rankType == "Scalar":
        newComponentMeshVarDict = _uw.dictionary.UpdateDictWithComponent(   name=componentName,
                                                                            Type="MeshVariable",
                                                                            Rank=rankType,
                                                                            DataType=dataType,
                                                                            mesh =meshName
                                                                            )

    return newComponentMeshVarDict

########################################################################################################################

# Suggest we deprecate this function and use "mergeType" (or change this to be something different) in the function above


def _meshVariableClobber( componentName="myMeshVariable",
                          dim=2,
                          rankType="Vector",
                          dataType="Double",
                          meshName="velocityMesh"):
    """
    Set up a Mesh Variable in the dictionary but Clobber it if it already exists.

    This function is for people who wear steel-cap boots, smoke cigars, drink their rum straight,
    are probably armed and don't want no steenkin' warnings.

    Possible rankTypes are ( Scalar, Vector )

    Needed by DofLayouts which are in turn needed by FeVariables (e.g. VelocityFeVariable)

    Args:
        componentName (string)         : The name of the mesh Variable.
        rankType  (string)             : type of Variable ( Scalar, Vector )
        dataType  (string)             : data type (Int, Float, Double, Possums)
        dim  (int)                     :
        meshName (string)              : Mesh which this Variable is defined on
    """
    globalDict = _uw.dictionary.GetDictionary()
    newComponentMeshVarDict = -1
    if componentName in globalDict["components"]:  # Then clobber it!
        if rankType == "Vector":
            namesList = ["vx", "vy", "vz"]
            newComponentMeshVarDict = _uw.dictionary.UpdateDictWithComponent(   name=componentName,
                                                                                Type="MeshVariable",
                                                                                Rank=rankType,
                                                                                DataType=dataType,
                                                                                VectorComponentCount=dim,
                                                                                names=namesList,
                                                                                mesh =meshName
                                                                                )
        if rankType == "Scalar":
            newComponentMeshVarDict = _uw.dictionary.UpdateDictWithComponent(   name=componentName,
                                                                                Type="MeshVariable",
                                                                                Rank=rankType,
                                                                                DataType=dataType,
                                                                                mesh =meshName
                                                                                )
    else:
        newComponentMeshVarDict = meshVariableCreate( componentName=componentName,
                                                      dim=dim,
                                                      rankType=rankType,
                                                      dataType=dataType,
                                                      meshName=meshName )

    return newComponentMeshVarDict
