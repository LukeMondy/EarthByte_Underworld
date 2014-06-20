# Solvers package - activates and configures solvers
import underworld as _uw

##############################################################################
# This code adds what is required to the python dictionary
# to set up the Solvers for Underworld.
# We eventually pass the python dictionary back to Underworld
# and Underworld then uses this information to configure and set
# itself up.
##############################################################################

'''
This code adds what is required to the python dictionary for the Underworld Solvers
Ultimately the global Dictionary gets passed back to Underworld which then actually creates the simulation

'''


def uzawaCreate(preconditionerMatrix="preconditioner", tol="1.0e-5", monitor=False, maxIts=5000, minIts=1, verbose=False):
    """
    Set up the Uzawa Solver in the dictionary

    Args:
        preconditionerMatrix  : a matrix for preconditioning.
    """
    globalDict = _uw.dictionary.GetDictionary()
    _uw.utils.warnMissingComponent(globalDict, preconditionerMatrix )

    uzawaDict = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                        name = "uzawa",
                                                        Type = "Stokes_SLE_UzawaSolver",
                                                        velocitySolver = "matrixSolver",  # not used atm?
                                                        Preconditioner = preconditionerMatrix,
                                                        tolerance      = tol,
                                                        monitor = str(monitor),
                                                        maxIterations = maxIts,
                                                        minIterations = minIts
                                                        )

    globalDict["solver"] = "UzawaSolver"
    if verbose:
        print "Now you need to set up a system of equations to attach this to:"
        print "use the stokesEquationCreate function in the equations module"

    return uzawaDict


def stokesBlockKSPInterfaceCreate(preconditionerMatrix="preconditioner", tol="1.0e-5", verbose=False):
    """
    Set up the Block Stokes SCR Solver in the dictionary

    Args:
        preconditionerMatrix  : a matrix for preconditioning.
    """
    globalDict = _uw.dictionary.GetDictionary()
    _uw.utils.warnMissingComponent(globalDict, preconditionerMatrix )

    sbkiDict = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                       name = "stokesblockkspinterface",
                                                       Type = "StokesBlockKSPInterface",
                                                       Preconditioner = preconditionerMatrix
                                                       )

    _uw.solvers.setSolverType("StokesBlockKSP")  # So the options function can check is this solver active
    # _uw.solvers._solverType = "StokesBlockKSP"  # So the options function can check is this solver active

    globalDict["solver"] = "StokesBlockInterfaceSolver"
    if verbose:
        print "Now you need to set up a system of equations to attach this to:"
        print "use the stokesEquationCreate function in the equations module"

    return sbkiDict
