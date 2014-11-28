import underworld as _underworld


def StokesBlockKSP( penaltyNumber=0.0 ):
    """
    Creates new components required to light up the Stokes Block solvers

    Args:
        penaltyNumber (float)      : Augmented Lagrangian penalty number (0.0 reverts to non-penalty form)

    Returns:
        Nothing
    """

    global _solversOn

    if _underworld.isInitialised() == False:
        print "StGermain has not been initialised yet. You need to run one of the Init functions."
        return

    if _underworld.isConstructed() == True or _underworld.isBuilt() == True:
        print "StGermain has already been Constructed/Built/Initialised."
        print "Modifications to the dictionary will possibly have no effect."
        print "You may need to restart, either by quitting or calling Finalise()"

    print "Initialising Velic & Moresi robust solver toolbox (penalty = {})".format(penaltyNumber)

    # The preconditioner does not use material point information
    integrationSwarm = "gaussSwarm"

    globalDict = _underworld.dictionary.GetDictionary()


# Add the solver toolbox

    globalDict["import"].append('Solvers')

# Root level parameters for the solver (we can change where they live later !)

    globalDict['penaltyNumber'] = penaltyNumber  # OK default value
    globalDict['hFactor'] = 0.0  # What is this number ??

##############################################################################
# Dictionary changes / definitions from AugLagStokesSLE-GtMG.xml
##############################################################################

# RHS vector for mass matrix

    scratchForceVector = _underworld.dictionary.UpdateDictWithComponent( globalDict,
                                                                         name="scratchForceVector",
                                                                         Type="ForceVector",
                                                                         FeVariable= "PressureField",
                                                                         ExtraInfo= "context"
                                                                         )

    mass_matrix = _underworld.dictionary.UpdateDictWithComponent( globalDict,
                                                                  name="m_matrix",
                                                                  Type="StiffnessMatrix",
                                                                  RowVariable="PressureField",
                                                                  ColumnVariable="PressureField",
                                                                  RHS=scratchForceVector["name"],
                                                                  allowZeroElementContributions="True"
                                                                  )

    mass_matrix_assembly = _underworld.dictionary.UpdateDictWithComponent( globalDict,
                                                                           name="PressureMassMatrixTerm",
                                                                           Type="PressMassMatrixTerm",
                                                                           Swarm=integrationSwarm,
                                                                           GeometryMesh="linearMesh",
                                                                           StiffnessMatrix=mass_matrix["name"]
                                                                           )

    print " *  Initialised Pressure Mass Matrix as {}".format(mass_matrix["name"])
    print " *  Initialised Pressure Mass Matrix RHS as {}".format(scratchForceVector["name"])
    print " *  Initialised Pressure Mass Matrix Assembly Operator as {}".format(mass_matrix_assembly["name"])


# Modifications to the Stokes equation are the same (?) for GTG and GTMG
# ** Except for the m_matrix term

#    del globalDict["components"]["stokesEqn"]

# Note, since the names of some of the input parameters for this component start with
# numbers, we have to pass the arguments to the function as a dictionary explicitly.
# (Maybe this is work adopting anyway, since it is general)

    stokesEqn = globalDict["components"]["stokesEqn"] = {
        "globalDict": "stgdict",
        "name": "stokesEqn",
        "Type": "AugLagStokes_SLE",
        "SLE_Solver": "uzawa",
        "Context": "context",
        "StressTensorMatrix": "k_matrix",
        "2ndStressTensorMatrix": "",
                                 "GradientMatrix": "g_matrix",
                                 "MassMatrix": mass_matrix["name"],
                                 "JunkForceVector": "junk_force",
                                 "DivergenceMatrix": "d_matrix",
                                 "CompressibilityMatrix": "c_matrix",
                                 "VelocityVector": "solutionVelocity",
                                 "PressureVector": "solutionPressure",
                                 "ForceVector": "mom_force",
                                 "2ndForceVector": "",
                                 "ContinuityForceVector": "cont_force",
                                 "killNonConvergent": "false",
                                 "nonLinearMaxIterations": "nonLinearMaxIterations",
                                 "nonLinearTolerance": "nonLinearTolerance",
                                 "makeConvergenceFile": "false",
                                 "penaltyNumber": penaltyNumber,
                                 "hFactor": "hFactor"
    }

    print " *  Re-Initialised Stokes Equation as {}".format(stokesEqn["name"])

##############################################################################
# Dictionary changes / definitions from VelocityMassMatrixSLE.xml
##############################################################################

    scratchVelForceVector = _underworld.dictionary.UpdateDictWithComponent( globalDict,
                                                                            name="scratchVelForceVector",
                                                                            Type="ForceVector",
                                                                            FeVariable= "VelocityField",
                                                                            ExtraInfo= "context"
                                                                            )

    vel_mass_matrix = _underworld.dictionary.UpdateDictWithComponent( globalDict,
                                                                      name="vm_matrix",
                                                                      Type="StiffnessMatrix",
                                                                      RowVariable="VelocityField",
                                                                      ColumnVariable="VelocityField",
                                                                      RHS=scratchVelForceVector["name"],
                                                                      allowZeroElementContributions="True"
                                                                      )

    vel_mass_matrix_assembly = _underworld.dictionary.UpdateDictWithComponent( globalDict,
                                                                               name="VelocityMassMatrixTerm",
                                                                               Type="VelocityMassMatrixTerm",
                                                                               Swarm=integrationSwarm,
                                                                               GeometryMesh="linearMesh",
                                                                               StiffnessMatrix=vel_mass_matrix["name"]
                                                                               )

    # Should the mass matrix be PIC / Non PIC ?

    print " *  Initialised Velocity Mass Matrix as {}".format(vel_mass_matrix["name"])
    print " *  Initialised Velocity Mass Matrix RHS as {}".format(scratchVelForceVector["name"])
    print " *  Initialised Velocity Mass Matrix Assembly Operator as {}".format(vel_mass_matrix_assembly["name"])

    stokesEqn["VelocityMassMatrix"] = vel_mass_matrix["name"]
    stokesEqn["VMassForceVector"]   = scratchVelForceVector["name"]

    print " *  Added Velocity Mass Matrix to Stokes Eqn"


##############################################################################
# Dictionary changes / definitions from kspinterface.xml
##############################################################################

    if "uzawa" in globalDict["components"].keys():
        del globalDict["components"]["uzawa"]

    stokesblockkspinterface = _underworld.dictionary.UpdateDictWithComponent( globalDict,
                                                                              name="stokesblockkspinterface",
                                                                              Type="StokesBlockKSPInterface",
                                                                              Preconditioner="preconditioner"
                                                                              )

    stokesEqn["SLE_Solver"] = stokesblockkspinterface["name"]

    print " *  Added Block Solver to Stokes Eqn"

    _solversOn = True
    _solverType = "StokesBlockKSP"

    if penaltyNumber != 0.0:
        _solverType = _solverType + "+Penalty"

    _underworld.solvers.setSolverType( _solverType )
