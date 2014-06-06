
import underworld as _underworld 

##############################################################################
## Dictionary changes / definitions from MultigridForRegularSCR.xml or 
##  from MultigridForRegular.xml
##############################################################################

def multigrid( mgLevels ):
    """

    """

    global _solversOn

    _solverType = _underworld.solvers.solverType()

    if _solverType == "None":
        print "solvers.multigrid(): You should first call one of the solvers.setup functions"
        return


    globalDict = _underworld.GetCurrentPythonDictionary()    

    globalDict['mgLevels'] = mgLevels

    # The old-style solver requires a plugin for the MG, the new does not

    if "Uzawa" in _solverType:
        globalDict["plugins"].append( {"Type" : "StgFEM_Multigrid" , "Context":"context" } ) 


    mgsolver = _underworld.NewComponentEntryInStgDict( globalDict,
                                name="mgSolver",
                                Type="PETScMGSolver",
                                levels=mgLevels,
                                opGenerator="mgOpGenerator"
                                )

    mgOpGenerator = _underworld.NewComponentEntryInStgDict( globalDict,
                                name="mgOpGenerator",
                                Type="SROpGenerator",
                                fineVariable="VelocityField"
                                )

    print " *  Added MG Solver to Stokes Eqn (using {} levels)".format(mgLevels)


    _solverType = _solverType + "+Multigrid"
    _underworld.solvers.setSolverType( _solverType )



    return    

