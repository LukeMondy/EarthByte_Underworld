
import underworld as _underworld


def Uzawa( PIC=True ):
    """
    Creates new components required to light up the Uzawa solvers
        Previously, this would be activated by the following XML files
            Underworld/StokesFlowUzawa.xml     or
            Underworld/StokesFlowUzawaPIC.xml

    Args:
        PIC (Bool)                   : Whether to use gauss or PIC integration

    Returns:
        Nothing
    """

    global _solversOn

    if _underworld.isInitialised() == False:
        print "StGermain has not been initialised yet. You need to run one of the Init functions."
        return

    if _underworld.isConstructed() == True or _underworld.isBuilt == True():
        print "StGermain has already been Constructed/Built/Initialised."
        print "Modifications to the dictionary will possibly have no effect."
        print "You may need to restart, either by quitting or calling Finalise()"

    if PIC == True:
        integrationSwarm = "picIntegrationPoints"
        _solverType = "Uzawa+PIC"
    else:
        integrationSwarm = "gaussSwarm"
        _solverType = "Uzawa"

    globalDict = _underworld.GetCurrentDictionary()

    preconditioner = _underworld.dictionary.UpdateDictWithComponent( globalDict,
                                                                     name="preconditioner",
                                                                     Type="StiffnessMatrix",
                                                                     RowVariable="PressureField",
                                                                     ColumnVariable="PressureField",
                                                                     RHS="cont_force",
                                                                     allowZeroElementContributions="True"
                                                                     )
    preconditionerTerm = _underworld.dictionary.UpdateDictWithComponent( globalDict,
                                                                         name="preconditionerTerm",
                                                                         Type="UzawaPreconditionerTerm",
                                                                         Swarm=integrationSwarm,
                                                                         StiffnessMatrix=preconditioner["name"]
                                                                         )

    uzawaSolver = _underworld.dictionary.UpdateDictWithComponent( globalDict,
                                                                  name="uzawaSolver",
                                                                  Type="Stokes_SLE_UzawaSolver",
                                                                  velocitySolver="matrixSolver",
                                                                  Preconditioner=preconditioner["name"],
                                                                  tolerance="1.0e-5",
                                                                  monitor="false",
                                                                  maxIterations="5000",
                                                                  minIterations="1"
                                                                  )

    # Not really a solver as such !!
    stokesEqn = _underworld.dictionary.UpdateDictWithComponent( globalDict,
                                                                name="stokesEqn",
                                                                Type="Stokes_SLE",
                                                                SLE_Solver="uzawa",
                                                                Context="context",
                                                                StressTensorMatrix="k_matrix",
                                                                GradientMatrix="g_matrix",
                                                                DivergenceMatrix="",
                                                                CompressibilityMatrix="c_matrix",
                                                                VelocityVector="solutionVelocity",
                                                                PressureVector="solutionPressure",
                                                                ForceVector="mom_force",
                                                                ContinuityForceVector="cont_force",
                                                                killNonConvergent="false",
                                                                nonLinearMaxIterations="nonLinearMaxIterations",
                                                                nonLinearTolerance="nonLinearTolerance",
                                                                makeConvergenceFile="false"
                                                                )

    _solversOn = True

    return
