# Matrix package
import underworld as _uw
#import underworld as _uw
#import underworld.utils._utils as utils
##############################################################################
# This code adds what is required to the python dictionary
# to set up the Matrices for a system of equations
# We eventually pass the python dictionary back to Underworld
# and Underworld then uses this information to configure and set
# itself up.
##############################################################################

'''
This code adds what is required to the python dictionary for Matrices required by a system of equations
Ultimately the global Dictionary gets passed back to Underworld which then actually creates the simulation

'''
# Create an Element Assembly term with attached matrix and force-vector.
# The dependence is:
# Element Assembly Term
#               -> Matrix  (attached to assembly term)
#                    -> Force Vector  (attached to matrix)
# hence this function creates 3 structs in the dictionary
#
# This function attempts to make it very easy to create the usual suspects


def matrixCreate( matrixName="stokes_matrix",
                  rowFeVariable="", colFeVariable="",
                  rhsVector="", transposeRHSVector="",
                  allowZeroElementContributions="False",
                  matrixTermType="",
                  intSwarmName="",
                  createTransposeRHS=False,
                  comment="",
                  Context="context"):
    """
    Create a matrix element Term with Matrix and its (required) ForceVector

    Args:
       rowFeVariable      : An FeVariable that determines the row dimension (required)
       colFeVariable      : An FeVariable that determines the col dimension (required)
       matrixTermType     : A type for the element assembly routine         (required)
                            One of :
                               [ ConstitutiveMatrixCartesian, GradientStiffnessMatrixTerm, UzawaPreconditionerTerm,
                                 PressMassMatrixTerm .... ]

       rhsVector          : A force vector for the RHS of the matrix. (will create a default one if not supplied)

       createTransposeRHS : If True will create a default transposeRHSVector (a RHS for the transpose of the matrix)
       transposeRHSVector : A force vector for the transpose of the matrix (optional)

       intSwarmName       : An integration swarm for the matrix element assembly routine struct (required)

    """
    globalDict = _uw.dictionary.GetDictionary()

    mList = [ "ConstitutiveMatrixCartesian", "GradientStiffnessMatrixTerm",
              "UzawaPreconditionerTerm", "PressMassMatrixTerm"]
    if matrixTermType == "":
        _uw.utils.sendError("matrixTermType is a required argument")
        print "Choose one of:"
        print mList
        return -1

    if matrixTermType not in mList:
        _uw.utils.sendWarning("The matrix term type: " + matrixTermType + " is not in a list of known types")
        print "The known types are:"
        print mList

    if rowFeVariable == "":
        _uw.utils.sendError("rowFeVariable is a required argument")
        print "create one using feVariableCreate in the fields module"
        return -1
    if colFeVariable == "":
        _uw.utils.sendError("colFeVariable is a required argument")
        print "create one using feVariableCreate in the fields module"
        return -1
    if intSwarmName == "":
        _uw.utils.sendError("intSwarmName is a required argument")
        print "create one using integrationSwarmCreate in the swarms module"
        return -1

    _uw.utils.warnMissingComponent(globalDict, intSwarmName )
    _uw.utils.warnMissingComponent(globalDict, rowFeVariable )
    _uw.utils.warnMissingComponent(globalDict, colFeVariable )

    if rhsVector != "" and rhsVector != None:
        _uw.utils.warnMissingComponent(globalDict, rhsVector )
    else:  # create a default one if nothing passed in
        newRHSName = matrixName + "RHS"
        newRHSName = _uw.utils.checkForNewComponentName(globalDict, newRHSName)
        newRHS     = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                             name = newRHSName,
                                                             Type="ForceVector",
                                                             FeVariable=rowFeVariable,
                                                             ExtraInfo ="context"
                                                             )

        rhsVector = newRHSName

    # no default creation for the transpose vector.
    if transposeRHSVector != "" and transposeRHSVector != None:
        _uw.utils.warnMissingComponent(globalDict, transposeRHSVector )
    if createTransposeRHS:
        newRHSName = matrixName + "TransposeRHS"
        newRHSName = _uw.utils.checkForNewComponentName(globalDict, newRHSName)
        newRHS     = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                             name = newRHSName,
                                                             Type="ForceVector",
                                                             FeVariable=colFeVariable,
                                                             ExtraInfo ="context"
                                                             )
        transposeRHS = newRHSName

    matrixName = _uw.utils.checkForNewComponentName(globalDict, matrixName)
    newMatrixDict = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                            name     = matrixName,
                                                            Type     = "StiffnessMatrix",
                                                            RowVariable    = rowFeVariable,
                                                            ColumnVariable = colFeVariable,
                                                            RHS            = rhsVector,
                                                            transposeRHS   = transposeRHSVector,
                                                            comment        = comment,
                                                            allowZeroElementContributions = str(allowZeroElementContributions)
                                                            )

    termName = matrixName + "AssemblyTerm"
    # If the Type is ConstitutiveMatrixCartesian then we HAVE to name it
    # 'constitutiveMatrix' or the PpcManager doesn't get built unless we explicitly set it up!
    # This is too restrictive.
    # We can't build the ppc_Manager unless we have the int Swarm for the constitutive matrix
    # But we don't want to know the names of these things before they are built
    # because we want the option of creating some names automatically and also
    # the freedom to name them what we want.
    # For now we do this HACK
    commentStr = ""
    if matrixTermType == "ConstitutiveMatrixCartesian":
        termName = "constitutiveMatrix"
        commentStr = "If the Type is ConstitutiveMatrixCartesian then we HAVE to name it 'constitutiveMatrix' or the PpcManager doesn't get built unless we explicitly set it up!"
    termName = _uw.utils.checkForNewComponentName(globalDict, termName)
    newTermDict = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                          name    = termName,
                                                          Type    = matrixTermType,
                                                          Swarm   = intSwarmName,
                                                          StiffnessMatrix = matrixName,
                                                          comment         = commentStr
                                                          )
    return [newMatrixDict, newTermDict]


def vectorCreate( vectorName="solutionVelocity",
                  feVariable="",
                  vectorType="SolutionVector",
                  extraInfo="",
                  comment="",
                  Context="context"):

    globalDict = _uw.dictionary.GetDictionary()

    vectorName  = _uw.utils.checkForNewComponentName(globalDict, vectorName)
    newVectDict = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                          name     = vectorName,
                                                          Type     = vectorType,
                                                          FeVariable = feVariable,
                                                          ExtraInfo  = Context   # might be better if the name on the struct was 'Context'?
                                                          )
    return newVectDict
