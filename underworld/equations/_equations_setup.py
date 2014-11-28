# Equations package - activates and configures Systems of Equations
import underworld as _uw
import underworld.matrix as _matrix
import underworld.solvers as _solvers
##############################################################################
# This code adds what is required to the python dictionary
# We eventually pass the python dictionary back to Underworld
# and Underworld then uses this information to configure and set
# itself up.
##############################################################################

'''
This code adds what is required to the python dictionary for Systems of Equations
Ultimately the global Dictionary gets passed back to Underworld which then actually creates the simulation

'''

# will add compressibility etc later
# have to have a little think about the "buoyancy" stuff (see _physics.py)
# could read PIC flag of main dictionary. But like to have it explicit here at least for now

# velocityField etc have default names here.
# these names get created in _geometry_setup in meshQ1P0CartesianCreate
# we can get these names as output from meshQ1P0CartesianCreate and then pass them into here to ensure consistency


def stokesSystemCreate(equationName="stokesEqn", solver="uzawa",
                       buoyancy=False,
                       buoyancyType="compositional",
                       Rayleigh=1e6,
                       compressibility=False,
                       oneOverlambda=10.0,
                       pic=True, penaltyNumber=0.0,
                       velocityField ="",
                       pressureField ="",
                       gaussIntSwarm ="",
                       picIntSwarm   ="",
                       StoreDensityOnParticles = "False"
                       ):
    """
    Create the full Stokes system of equations.
    Creates the required matrices and solver.

    This function also sets the FeVariables and Integration Swarm names on the top level of the
    dictionary as parameters.

    Args:
       solver        : One of "stokesblockkspinterface" or "uzawa"
       buoyancy      : Adds Buoyancy force term
       buoyancyType  : One of ['compositional', 'thermal']
       Rayleigh      : Rayleigh Number for thermal buoyancy
       pic           : True for particles or False for Gauss swarms
       penaltyNumber : Penalty number for the "stokesblockinterface" Augmented Lagrangian solver
       StoreDensityOnParticles: If buoyancyType is compositional, set this to true in order to visualise temperature-dependent density on "your material swarm"-Density

    """

    globalDict = _uw.dictionary.GetDictionary()

    c_mat = ""

    if velocityField == "":
        velocityField = globalDict["info"]["velocityField"]
    if pressureField == "":
        pressureField = globalDict["info"]["pressureField"]
    if gaussIntSwarm == "":
        gaussIntSwarm = globalDict["info"]["gaussIntSwarm"]

    if not pic:
        picIntSwarm = gaussIntSwarm
    else:
        if picIntSwarm == "":
            picIntSwarm = globalDict["info"]["picIntSwarm"]

    _uw.dictionary.addCheckPointVariables([velocityField, pressureField])

    # Create Stokes Matrices and Vectors
    [kmatDict, kmatTermDict] = _matrix.setup.matrixCreate(matrixName     = "k_matrix",
                                                          rowFeVariable  = velocityField,
                                                          colFeVariable  = velocityField,
                                                          matrixTermType = "ConstitutiveMatrixCartesian",
                                                          intSwarmName   = picIntSwarm,
                                                          comment        = "momentum"
                                                          )
    mom_forceVector = kmatDict["RHS"]  # momentum force vector - needed for _stokesCreate and  buoyancy
    # Need to make a transposeRHSvector first for this matrix
    transRHS = _matrix.setup.vectorCreate(vectorName="g_matrixTransRHS",
                                          feVariable=pressureField,
                                          vectorType="ForceVector",
                                          Context   ="context" )
    [gmatDict, gmatTermDict] = _matrix.setup.matrixCreate(matrixName     = "g_matrix",
                                                          rowFeVariable  = velocityField,
                                                          colFeVariable  = pressureField,
                                                          matrixTermType = "GradientStiffnessMatrixTerm",
                                                          intSwarmName   = gaussIntSwarm,
                                                          rhsVector      = mom_forceVector,
                                                          transposeRHSVector = "g_matrixTransRHS",
                                                          comment        = "gradient"
                                                          )
    cont_forceVector = "g_matrixTransRHS"  # continuity force vector - needed for _stokesCreate
    [pmatDict, pmatTermDict] = _matrix.setup.matrixCreate(matrixName     = "preconditioner",
                                                          rowFeVariable  = pressureField,
                                                          colFeVariable  = pressureField,
                                                          matrixTermType = "UzawaPreconditionerTerm",
                                                          intSwarmName   = picIntSwarm,
                                                          comment        = "preconditioner"
                                                          )
    mass_forceVector = ""
    vmass_forceVector = ""
    vm_mat = ""
    m_mat = ""
    if solver == "stokesblockkspinterface":
        _uw.dictionary.importToolBox('Solvers')
        _solvers.setup.stokesBlockKSPInterfaceCreate(preconditionerMatrix=pmatDict["name"])
        # create some auxiliary stuff
        [mmatDict, mmatTermDict] = _matrix.setup.matrixCreate(matrixName     = "m_matrix",
                                                              rowFeVariable  = pressureField,
                                                              colFeVariable  = pressureField,
                                                              matrixTermType = "PressMassMatrixTerm",
                                                              intSwarmName   = gaussIntSwarm,
                                                              comment        = "pressure Mass matrix for Aug Lag Solver"
                                                              )
        # The PressMassMatrixTerm Assembly Term needs the 'GeometryMesh' set on it
        # Not sure if we need to make the matrixCreate func more clever or not?
        # Matrix Assembly terms don't usually have Meshes attached
        # We know that the correct Mesh is attached to the Q1(velocity) FeVariable
        # going to have to think about this a little.
        # m_matrixAssemblyTerm is the auto name the matrix Create routine makes up
        globalDict["components"]["m_matrixAssemblyTerm"]["GeometryMesh"] = globalDict["components"][velocityField]["FEMesh"]  # Attach the Q1 (Velocity) Mesh here.
        mass_forceVector = mmatDict["RHS"]
        m_mat = mmatDict["name"]
        [vmmatDict, vmmatTermDict] = _matrix.setup.matrixCreate(matrixName     = "vm_matrix",
                                                                rowFeVariable  = velocityField,
                                                                colFeVariable  = velocityField,
                                                                matrixTermType = "VelocityMassMatrixTerm",
                                                                intSwarmName   = gaussIntSwarm,
                                                                comment        = "velocity Mass matrix for Aug Lag Solver"
                                                                )
        vmass_forceVector = vmmatDict["RHS"]
        vm_mat = vmmatDict["name"]
    else:
        # create uzawa solver
        _solvers.setup.uzawaCreate(preconditionerMatrix=pmatDict["name"])

    # Now we need Solution Vectors (maybe should put these at top level dictionary too)
    velSol   = _matrix.setup.vectorCreate(vectorName="solutionVelocity", feVariable=velocityField,
                                          vectorType="SolutionVector", Context=globalDict["context"])
    pressSol = _matrix.setup.vectorCreate(vectorName="solutionPressure", feVariable=pressureField,
                                          vectorType="SolutionVector", Context=globalDict["context"])

    # Add buoyancy
    if buoyancy:
        if buoyancyType == "compositional":
            addBuoyancy(forceVector=mom_forceVector, intSwarm=picIntSwarm, StoreDensityOnParticles=StoreDensityOnParticles)  # changed from picIntSwarm..test this
        if buoyancyType == "thermal":
            addThermalBuoyancy(forceVector=mom_forceVector, intSwarm=picIntSwarm, Ra=str(Rayleigh))

    # Can now set up the stokesEqn entry which connects the Solver to the Matrix system of equations

    stokesCreate( equationName="stokesEqn",
                  solver=solver,
                  preconditionerMatrix = pmatDict["name"],
                  stressTensorMatrix   = kmatDict["name"],
                  gradientMatrix       = gmatDict["name"],
                  velocitySolution     = velSol["name"],
                  pressureSolution     = pressSol["name"],
                  momentumForceVector   = mom_forceVector,
                  continuityForceVector = cont_forceVector,
                  context               = globalDict["context"],
                  massMatrix              = m_mat,
                  junkForceVector         = mass_forceVector,
                  velocityMassMatrix      = vm_mat,
                  velocityMassForceVector = vmass_forceVector,
                  penaltyNumber           = penaltyNumber,
                  compressibilityMatrix   = c_mat
                  )
    return


def stokesCreate(equationName="stokesEqn",
                 solver="uzawa",
                 preconditionerMatrix="preconditioner",
                 stressTensorMatrix="k_matrix",    # maybe we could store these default names to ensure consistency across modules.
                 gradientMatrix="g_matrix",
                 divergenceMatrix="",
                 velocitySolution="solutionVelocity",
                 pressureSolution="solutionPressure",
                 momentumForceVector="mom_force",
                 continuityForceVector="cont_force",
                 context="context",
                 secondStressTensorMatrix="",
                 massMatrix="m_matrix",
                 junkForceVector="junk_force",
                 secondForceVector="",
                 velocityMassMatrix="vm_matrix",
                 velocityMassForceVector="vm_force",
                 compressibilityMatrix="",
                 penaltyNumber = 0.1
                 ):
    """
    Sets up a stokes equation system with provided matrices and solution vectors.
    """

    globalDict = _uw.dictionary.GetDictionary()
    missing = _uw.utils.warnMissingComponent(globalDict, solver )
    _uw.utils.warnMissingComponent(globalDict, preconditionerMatrix )
    _uw.utils.warnMissingComponent(globalDict, stressTensorMatrix )
    _uw.utils.warnMissingComponent(globalDict, gradientMatrix )
    #_uw.utils.warnMissingComponent(globalDict, divergenceMatrix )
    _uw.utils.warnMissingComponent(globalDict, velocitySolution )
    _uw.utils.warnMissingComponent(globalDict, pressureSolution )
    _uw.utils.warnMissingComponent(globalDict, momentumForceVector )
    _uw.utils.warnMissingComponent(globalDict, continuityForceVector )
    _uw.utils.warnMissingComponent(globalDict, context )

    stokesDict = {}  # empty dictionary
    if solver == "uzawa":
        stokesDict = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                             name = equationName,
                                                             Type = "Stokes_SLE",
                                                             SLE_Solver = solver,
                                                             Context    = context,
                                                             StressTensorMatrix     = stressTensorMatrix,
                                                             GradientMatrix         = gradientMatrix,
                                                             DivergenceMatrix       = divergenceMatrix,
                                                             CompressibilityMatrix  = compressibilityMatrix,
                                                             VelocityVector         = velocitySolution,
                                                             PressureVector         = pressureSolution,
                                                             ForceVector            = momentumForceVector,
                                                             ContinuityForceVector  = continuityForceVector,
                                                             killNonConvergent      = str(False),
                                                             nonLinearMaxIterations = "nonLinearMaxIterations",
                                                             nonLinearTolerance     = "nonLinearTolerance",
                                                             makeConvergenceFile    = str(False)
                                                             )
    if solver == "stokesblockkspinterface":
        stokesDict = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                             name = equationName,
                                                             Type = "AugLagStokes_SLE",
                                                             SLE_Solver = solver,
                                                             Context    = context,
                                                             StressTensorMatrix     = stressTensorMatrix,
                                                             GradientMatrix         = gradientMatrix,
                                                             DivergenceMatrix       = divergenceMatrix,
                                                             CompressibilityMatrix  = compressibilityMatrix,
                                                             VelocityVector         = velocitySolution,
                                                             PressureVector         = pressureSolution,
                                                             ForceVector            = momentumForceVector,
                                                             ContinuityForceVector  = continuityForceVector,
                                                             killNonConvergent      = str(False),
                                                             nonLinearMaxIterations = "nonLinearMaxIterations",
                                                             nonLinearTolerance     = "nonLinearTolerance",
                                                             makeConvergenceFile    = str(False),
                                                             MassMatrix              = massMatrix,
                                                             penaltyNumber           = penaltyNumber,
                                                             hFactor                 = "hFactor",
                                                             VelocityMassMatrix      = velocityMassMatrix,
                                                             VMassForceVector        = velocityMassForceVector,
                                                             JunkForceVector         = junkForceVector
                                                             )
        # have to get a bit "hacky" because of bad names we want in the dictionary. We should just change these names in UW itself.
        globalDict["components"][equationName]["2ndStressTensorMatrix"] = secondStressTensorMatrix
        globalDict["components"][equationName]["2ndForceVector"] = secondForceVector
    return stokesDict


def addBuoyancy(forceVector="mom_force", intSwarm="", temperatureField="TemperatureField", comment="", StoreDensityOnParticles="False"):

    globalDict = _uw.dictionary.GetDictionary()

    if temperatureField not in globalDict["components"].keys():
        comment = "Temperature field is dummy variable here."

    if "gravity" not in globalDict.keys():
        _uw.utils.sendWarning("The 'gravity' parameter is missing from the Dictionary: adding it with default value = 1")
        globalDict["gravity"] = 1

    buoyancy = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                       name = "buoyancyForceTerm",
                                                       Type = "BuoyancyForceTerm",
                                                       TemperatureField = temperatureField,  # optional: temp of 0.0 used if no Field
                                                       ForceVector = forceVector,
                                                       Swarm = intSwarm,
                                                       gravity = "gravity",
                                                       comment = comment,
                                                       StoreDensityOnParticles = StoreDensityOnParticles
                                                       )
    return buoyancy


def addThermalBuoyancy(forceVector="mom_force", intSwarm="", temperatureField="TemperatureField", Ra="1e6", comment=""):

    globalDict = _uw.dictionary.GetDictionary()

    if temperatureField not in globalDict["components"].keys():
        _uw.utils.sendWarning("Temperature field is missing from dictionary")
    if "gravity" not in globalDict.keys():
        _uw.utils.sendWarning("The 'gravity' parameter is missing from the Dictionary: adding it with default value = 1")
        globalDict["gravity"] = 1

    buoyancy = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                       name = "thermalBuoyancyForceTerm",
                                                       Type = "ThermalBuoyancyForceTerm",
                                                       TemperatureField = temperatureField,  # optional: temp of 0.0 used if no Field
                                                       ForceVector = forceVector,
                                                       Swarm = intSwarm,
                                                       gravity = "gravity",
                                                       Ra      = str(Ra),
                                                       comment = comment
                                                       )
    return buoyancy


def advectionDiffusionEquationCreate(equationName="energyEqn", phiField="",
                                     gaussIntSwarm="", velocityField="",
                                     diffusivity=1.0,
                                     upwindFunc="DoublyAsymptoticAssumption",
                                     diffusivityVariable="",
                                     courantFactor=0.25,
                                     context="context"):
    """
    Set up an Advection Diffusion Equation and solver.

    Based on:
    "Streamline upwind/Petrov-Galerkin formulations for convection dominated flows
            with particular emphasis on the incompressible Navier-Stokes equations"
    AN Brooks, TJR Hughes,    Computer methods in applied mechanics and engineering 32 (1), 199-259, 1982

    Args:
      phiField   : scalar field (usually the TemperatureField)
      upwindFunc : one of ["DoublyAsymptoticAssumption", "CriticalAssumption", "Exact"]
    """

    globalDict = _uw.dictionary.GetDictionary()

    if phiField == "":  # should handle case where these don't exist I suppose...
        phiField = globalDict["info"]["temperatureField"]
    if velocityField == "":
        velocityField = globalDict["info"]["velocityField"]
    if gaussIntSwarm == "":
        gaussIntSwarm = globalDict["info"]["gaussIntSwarm"]

    if phiField not in globalDict["components"].keys():
        _uw.utils.sendWarning("Required field is missing from dictionary")
    if gaussIntSwarm not in globalDict["components"].keys():
        _uw.utils.sendWarning("Gauss swarm is missing from dictionary")
    if velocityField not in globalDict["components"].keys():
        _uw.utils.sendWarning("Velocity Field is missing from dictionary")

    massMat = "massMatrix"
    residual = "residual"

    _uw.dictionary.addCheckPointVariables([phiField, phiField + "-phiDotField"])

    force = _matrix.setup.vectorCreate(vectorName=residual, feVariable=phiField, vectorType="ForceVector")
    mass  = _matrix.setup.vectorCreate(vectorName=massMat, feVariable=phiField, vectorType="ForceVector")
    predictor = _uw.dictionary.UpdateDictWithComponent( globalDict, name="predictorMulticorrector", Type="AdvDiffMulticorrector")

    lumped = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                     name = "lumpedMassMatrixForceTerm",
                                                     Type = "LumpedMassMatrixForceTerm",
                                                     ForceVector = massMat,
                                                     Swarm = gaussIntSwarm
                                                     )

    defRes = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                     name= "defaultResidualForceTerm",
                                                     Type= "AdvDiffResidualForceTerm",
                                                     Swarm=gaussIntSwarm,
                                                     ForceVector=residual,
                                                     ExtraInfo=equationName,
                                                     VelocityField=velocityField,
                                                     DiffusivityVariable=diffusivityVariable,
                                                     defaultDiffusivity=str(diffusivity),
                                                     UpwindXiFunction=upwindFunc
                                                     )

    advDiffDict = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                          name = equationName,
                                                          Type = "AdvectionDiffusionSLE",
                                                          SLE_Solver = "predictorMulticorrector",
                                                          Context    = context,
                                                          MassMatrix = massMat,
                                                          Residual   = "residual",
                                                          PhiField   = phiField,
                                                          courantFactor = str(courantFactor)
                                                          )

    return
