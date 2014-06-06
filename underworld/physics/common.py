import underworld as _underworld 


## These ones belong with meshes !

def constant_mesh( meshName="constantMesh", elementMeshName="elementMesh"):
    """
    ConstantMesh.xml !
    Args: 
        name (string)
    """

    print " *  Piecewise constant mesh name - ", meshName


    _underworld.StgComponentToGlobalDict(
         name = meshName,
         Type = "FeMesh",
         elementType = "constant" 
         )

    _underworld.StgComponentToGlobalDict(
        name = meshName + "Generator",
        Type = "C0Generator",
        mesh = meshName,
        elementMesh = elementMeshName
        )
    

    return

def linear_mesh( meshName="linearMesh", meshSize = (32,32,32), span=[(0.0,0.0,0.0),(1.0,1.0,1.0)], dimensions=None, shadowDepth=1):
    """
        LinearMesh.xml !
        Args:
            name (string)
            meshSize = (int elementResI, int elementResJ, int elementResK) 
            span = [(float minX, float minY, float minZ), (float maxX, float maxY, float maxZ) ]
            dimensions = (int) - mesh dimensions dim=2 or dim=3
            shadowDepth = (int) - whatever that really means !
    """

    ## Should do some error checking that the tuples / lists are the correct size 
    ## for the specified dimensions 

    if dimensions == None:
        print "Please specify mesh dimensions !"
        return

    print " *  Piecewise linear mesh name - ", meshName


    _underworld.StgComponentToGlobalDict(
        name = meshName,
        Type = "FeMesh",
        elementType = "linear"
        )

    _underworld.StgComponentToGlobalDict(
        name = meshName+"Generator",
        Type = "CartesianGenerator",
        mesh = meshName,
        dims = dimensions,
        shadowDepth = shadowDepth,
        size = meshSize,  # Tuple should create a list
        minCoord = span[0],
        maxCoord = span[1]
        )

## Velocity and Pressure Field

## Can we generalise this to vector fields in general 

def velocity_field( fieldName="velocity", meshName="velocityMesh", dimensions=2  ):

    print " *  Velocity field name - ", fieldName


    _underworld.StgComponentToGlobalDict(
        name = fieldName,
        Type = "MeshVariable",
        mesh = meshName,
        Rank = "Vector",
        DataType = "double",
        VectorComponentCount = dimensions,
        names = ("vx", "vy", "vz")
        )

    _underworld.StgComponentToGlobalDict(
        name = fieldName+"Field",
        Type = "FeVariable",
        FEMesh = meshName,
        DofLayout = fieldName+"DofLayout",
        BC = fieldName+"BCs",
        IC = fieldName+"ICs",
        LinkedDofInfo = fieldName+"LinkedDofs",
        outputUnits="cm/yr"
        )

    _underworld.StgComponentToGlobalDict(
        name = fieldName+"BCs",
        Type = "CompositeVC",
        Data = meshName
        )

    _underworld.StgComponentToGlobalDict(
        name = fieldName+"ICs",
        Type = "CompositeVC",
        Data = meshName
        )

    _underworld.StgComponentToGlobalDict(
        name = fieldName+"DofLayout",
        Type = "DofLayout",
        MeshVariable = fieldName
        )


    ## Pre-define a "standard collection" of operators
    

    _underworld.StgComponentToGlobalDict(
        name = fieldName+"MagnitudeField",
        Type = "OperatorFeVariable",
        Operator = "Magnitude",
        Operand = fieldName+"Field",
        outputUnits="cm/yr"

        )

    _underworld.StgComponentToGlobalDict(
        name = fieldName+"GradientsField",
        Type = "OperatorFeVariable",
        Operator = "Gradient",
        Operand = fieldName+"Field"
        )

    # Do we use this for anything ?

    _underworld.StgComponentToGlobalDict(
        name = fieldName+"GradientsInvariantField",
        Type = "OperatorFeVariable",
        Operator = "Gradient",
        Operand = fieldName+"Field"
        )

    _underworld.StgComponentToGlobalDict(
        name = fieldName+"Xcomponent",
        Type = "OperatorFeVariable",
        Operator = "TakeFirstComponent",
        Operand = fieldName+"Field",
        outputUnits="cm/yr"
        )

    _underworld.StgComponentToGlobalDict(
        name = fieldName+"Ycomponent",
        Type = "OperatorFeVariable",
        Operator = "TakeSecondComponent",
        Operand = fieldName+"Field",
        outputUnits="cm/yr"
        )

    if (dimensions == 3):
        _underworld.StgComponentToGlobalDict(
            name = fieldName+"Zcomponent",
            Type = "OperatorFeVariable",
            Operator = "TakeThirdComponent",
            Operand = fieldName+"Field",
            outputUnits="cm/yr"
        )

 
     # Usually Strain rate ... but I want to tag it with the field name so this seems more appropriate (LM)
     # Probably s^-1 would be better too

    _underworld.StgComponentToGlobalDict(
        name = fieldName+"SymmetricGradient",
        Type = "OperatorFeVariable",
        Operator = "TensorSymmetricPart",
        Operand = fieldName+"GradientsField",
        outputUnits="yr^-1"
        )

    # Spin tensor (was called vorticity, but conventionally, we should use that name for curl(field))
    _underworld.StgComponentToGlobalDict(
        name = fieldName+"AntiSymmetricGradient",
        Type = "OperatorFeVariable",
        Operator = "TensorAntiSymmetricPart",
        Operand = fieldName+"GradientsField",
        outputUnits="yr^-1"
        )

    _underworld.StgComponentToGlobalDict(
        name = fieldName+"GradientInvariantSymmetricPart",
        Type = "OperatorFeVariable",
        Operator = "SymmetricTensor_Invariant",
        Operand = fieldName+"SymmetricGradient",
        outputUnits="yr^-1"
        )

    if dimensions == 2:

        _underworld.StgComponentToGlobalDict(
            name = fieldName+"GradientXXSymmetricPart",
            Type = "OperatorFeVariable",
            Operator = "TakeFirstComponent",
            Operand = fieldName+"SymmetricGradient",
            outputUnits="yr^-1"
            )

        _underworld.StgComponentToGlobalDict(
            name = fieldName+"GradientYYSymmetricPart",
            Type = "OperatorFeVariable",
            Operator = "TakeSecondComponent",
            Operand = fieldName+"SymmetricGradient",
            outputUnits="yr^-1"
            )

        _underworld.StgComponentToGlobalDict(
            name = fieldName+"GradientXYSymmetricPart",
            Type = "OperatorFeVariable",
            Operator = "TakeThirdComponent",
            Operand = fieldName+"SymmetricGradient",
            outputUnits="yr^-1"
            )

    else:   # 3D

        _underworld.StgComponentToGlobalDict(
            name = fieldName+"GradientXXSymmetricPart",
            Type = "OperatorFeVariable",
            Operator = "TakeFirstComponent",
            Operand = fieldName+"SymmetricGradient",
            outputUnits="yr^-1"
            )

        _underworld.StgComponentToGlobalDict(
            name = fieldName+"GradientYYSymmetricPart",
            Type = "OperatorFeVariable",
            Operator = "TakeSecondComponent",
            Operand = fieldName+"SymmetricGradient",
            outputUnits="yr^-1"
            )

        _underworld.StgComponentToGlobalDict(
            name = fieldName+"GradientZZSymmetricPart",
            Type = "OperatorFeVariable",
            Operator = "TakeThirdComponent",
            Operand = fieldName+"SymmetricGradient",
            outputUnits="yr^-1"
            )

        _underworld.StgComponentToGlobalDict(
            name = fieldName+"GradientXYSymmetricPart",
            Type = "OperatorFeVariable",
            Operator = "TakeFourthComponent",
            Operand = fieldName+"SymmetricGradient",
            outputUnits="yr^-1"
            )

        _underworld.StgComponentToGlobalDict(
            name = fieldName+"GradientXZSymmetricPart",
            Type = "OperatorFeVariable",
            Operator = "TakeFifthComponent",
            Operand = fieldName+"SymmetricGradient",
            outputUnits="yr^-1"
            )

        _underworld.StgComponentToGlobalDict(
            name = fieldName+"GradientYZSymmetricPart",
            Type = "OperatorFeVariable",
            Operator = "TakeSixthComponent",
            Operand = fieldName+"SymmetricGradient",
            outputUnits="yr^-1"
            )


#   <param name="velocityMesh">linearMesh</param>
#   <param name="elementMesh">linearMesh</param>

# Note: This creates "pressureField" and not "PressureField"
   
def pressure_field( fieldName="pressure", meshName="pressureMesh", dimensions=2  ):

    # This is general (cf. velocity - we just need to change Rank, DataType, Units)

    print " *  Pressure field name - ", fieldName



    _underworld.StgComponentToGlobalDict(
        name = fieldName,
        Type = "MeshVariable",
        Rank = "Scalar",
        DataType = "Double"
        )

    _underworld.StgComponentToGlobalDict(
        name = fieldName+"BCs",
        Type = "CompositeVC",
        Data = meshName
        )

    _underworld.StgComponentToGlobalDict(
        name = fieldName+"ICs",
        Type = "CompositeVC",
        Data = meshName
        )

    _underworld.StgComponentToGlobalDict(
        name = fieldName+"DofLayout",
        Type = "DofLayout",
        MeshVariable = fieldName
        )

    _underworld.StgComponentToGlobalDict(
        name = fieldName+"Field",
        Type = "FeVariable",
        FEMesh = meshName,
        DofLayout = fieldName+"DofLayout",
        BC = fieldName+"BCs",
        IC = fieldName+"ICs",
        LinkedDofInfo = fieldName+"LinkedDofs",
        outputUnits="GPa"
        )



def gauss_point_integration_swarm( gaussSwarmName="GaussSwarm", meshName="elementMesh", timeIntegratorName="timeIntegrator"):

    print " *  Gauss point integration component: ", gaussSwarmName

    # Should these 2 be private names or public ?

    _underworld.StgComponentToGlobalDict(
        name = "cellLayout",
        Type = "SingleCellLayout"
        )

    _underworld.StgComponentToGlobalDict(
        name = "particleLayout",
        Type = "GaussParticleLayout"
        )

    _underworld.StgComponentToGlobalDict(
        name = "gaussSwarm",
        Type = "IntegrationPointsSwarm",
        CellLayout = "cellLayout",
        ParticleLayout = "particleLayout",
        FeMesh = meshName,
        TimeIntegrator = timeIntegratorName
        )


def time_integration( timeIntegratorName="timeIntegrator", timeIntegratorOrder=2, simultaneous="False", contextName="context"):

    print " *  Time integration component: ", timeIntegratorName

    _underworld.StgComponentToGlobalDict(
        name = timeIntegratorName,
        Type = "TimeIntegrator",
        order= timeIntegratorOrder,
        simultaneous = simultaneous,
        Context = contextName
        )



def material_point_swarm( picIntSwarmName="picIntegrationPoints", matSwarmName="materialSwarm", 
                          meshName="elementMesh", timeIntegratorName="timeIntegrator", 
                          pcdvcResolution=(10,10,10), particlesPerCell=20,
                          velocityFieldName="velocityField", periodicBCsManagerName="periodicBCsManager" ):

    
    print " *  Material point integration component: ", picIntSwarmName


    _underworld.StgComponentToGlobalDict(
        name = "elementCellLayout",
        Type = "ElementCellLayout",
        mesh = meshName
        )

    _underworld.StgComponentToGlobalDict(
        name = "localLayout",
        Type = "MappedParticleLayout"
        )

    # We can always drill down and change stuff, but here are some defaults
    _underworld.StgComponentToGlobalDict(
        name = "integrationWeights",
        Type = "PCDVC",
        resolutionX = pcdvcResolution[0],
        resolutionY = pcdvcResolution[1],
        resolutionZ = pcdvcResolution[2],
        upperT = 25,
        lowerT = 0.6,
        maxDeletions = 3,
        maxSplits = 3,
        MaterialPointsSwarm = matSwarmName
        )

    _underworld.StgComponentToGlobalDict(
        name = picIntSwarmName,
        Type = "IntegrationPointsSwarm",
        ParticleLayout = "localLayout",
        FeMesh = meshName,
        WeightsCalculator = "integrationWeights",
        TimeIntegrator = timeIntegratorName,
        IntegrationPointMapper = "ipMapper"
        )    

    _underworld.StgComponentToGlobalDict(
        name = "ipMapper",
        Type = "CoincidentMapper",
        IntegrationPointsSwarm = picIntSwarmName,
        MaterialPointsSwarm = matSwarmName
        )    

    _underworld.StgComponentToGlobalDict(
        name = "materialSwarmParticleLayout",
        Type = "SpaceFillerParticleLayout",
        averageInitialParticlesPerCell = particlesPerCell
        )

    _underworld.StgComponentToGlobalDict(
        name = "pMovementHandler",
        Type = "ParticleMovementHandler"
        )

    _underworld.StgComponentToGlobalDict(
        name = "pShadowSync",
        Type = "ParticleShadowSync"
        )

        
    _underworld.StgComponentToGlobalDict(
        name = matSwarmName,
        Type = "MaterialPointsSwarm",
        CellLayout = "elementCellLayout",
        ParticleLayout = "materialSwarmParticleLayout",
        FeMesh = meshName,
        ParticleCommHandlers = ["pMovementHandler", "pShadowSync"],
        SplittingRoutine = "splittingRoutine",
        RemovalRoutine = "removalRoutine",
        EscapedRoutine = "escapedRoutine"
        )

    _underworld.StgComponentToGlobalDict(
        name = "materialSwarmAdvector",
        Type = "SwarmAdvector",
        Swarm = matSwarmName,
        TimeIntegrator = timeIntegratorName,
        VelocityField = velocityFieldName,
        PeriodicBCsManager = periodicBCsManagerName,
        allowFallbackToFirstOrder = "True"
        )
    

## Stokes flow setup & Solver (which we can remove)


def stokes_equation(velocityFieldName="velocityField", 
                    pressureFieldName="pressureField", 
                    gaussSwarmName="GaussSwarm",
                    nonLinearMaxIterations = "100",
                    nonLinearTolerance = "1e-4",
                    contextName="context"):


    print " *  Stokes Equation defined as: stokesEqn"



    # Solution vectors

    _underworld.StgComponentToGlobalDict(
        name = "solutionVelocity",
        Type = "SolutionVector",
        FeVariable = velocityFieldName
        )

    _underworld.StgComponentToGlobalDict(
        name = "solutionPressure",
        Type = "SolutionVector",
        FeVariable = pressureFieldName
        )
    
    # RHS vectors

    _underworld.StgComponentToGlobalDict(  
        name = "mom_force",
        Type = "ForceVector",
        FeVariable = velocityFieldName,
        ExtraInfo = contextName
        )
        
    _underworld.StgComponentToGlobalDict(  
        name = "cont_force",
        Type = "ForceVector",
        FeVariable = pressureFieldName,
        ExtraInfo = contextName
        )
        

    # Saddle point system: matrices

    _underworld.StgComponentToGlobalDict(
        name = "k_matrix",
        Type = "StiffnessMatrix",
        RowVariable = velocityFieldName,
        ColumnVariable = velocityFieldName,
        RHS = "mom_force",
        allowZeroElementContributions = "False"
        )

    _underworld.StgComponentToGlobalDict(
        name = "constitutiveMatrix",
        Type = "ConstitutiveMatrixCartesian",
        Swarm = gaussSwarmName,
        StiffnessMatrix = "k_matrix"
        )

    _underworld.StgComponentToGlobalDict(
        name = "g_matrix",
        Type = "StiffnessMatrix",
        RowVariable = velocityFieldName,
        ColumnVariable = pressureFieldName,
        RHS = "mom_force",
        transposeRHS = "cont_force",
        allowZeroElementContributions = "False"
        )


    _underworld.StgComponentToGlobalDict(
        name = "gradientStiffnessMatrixTerm",
        Type = "GradientStiffnessMatrixTerm",
        Swarm = gaussSwarmName,
        StiffnessMatrix = "g_matrix"
        )

    # SLE 

    _underworld.StgComponentToGlobalDict(
        name = "stokesEqn",
        Type = "Stokes_SLE",
        Context = contextName,
        StressTensorMatrix = "k_matrix",
        GradientMatrix = "g_matrix",
        DivergenceMatrix = "",
        CompressibilityMatrix = "c_matrix",
        VelocityVector = "solutionVelocity",
        PressureVector = "solutionPressure",
        ForceVector = "mom_force",
        ContinuityForceVector = "cont_force",
        killNonConvergent = "False",
        nonLinearMaxIterations = nonLinearMaxIterations,
        nonLinearTolerance = nonLinearTolerance,
        makeConvergenceFile = "False"
        )


# Need to expose / exchange some of these values with other components - but what is the best way ?

def thermal_compositional_buoyancy_rhs( forceVectorName="mom_force", integrationSwarmName="picIntegrationPoints",
                                        temperatureFieldName="TemperatureField", gravity=1.0):

    _underworld.StgComponentToGlobalDict(
        name = "buoyancyForceTerm",
        Type="BuoyancyForceTerm", 
        ForceVector = forceVectorName,
        TemperatureField = temperatureFieldName,
        Swarm = integrationSwarmName,
        gravity = gravity
        )

## Boundary conditions

def velocity_boundary_conditions( top    = (None,None,0.0 ,None,None,None),  
                                  bottom = (None,None,0.0 ,None,None,None), 
                                    left   = (0.0 ,None,None,None,None,None), 
                                    right  = (0.0 ,None,None,None,None,None), 
                                    front  = (None,None,None,None,0.0 ,None), 
                                    back   = (None,None,None,None,0.0 ,None), 
                                    periodicX="False",
                                    periodicY="False",
                                    periodicZ="False",
                                    dimensions=None):

    '''
        Semi-automatic boundary conditions for Cartesian or Regional Spherical meshes
        - specify for each wall (top, bottom, left, right, front, back)
        - Tuple of (vx, Fx, vy, Fy, vz, Fz) - v values will over-ride sigma. 
          Leave Fx,Fy,Fz as None if a zero-stress condition is appropriate (sigma_s not yet implemented)
        - periodicX / periodicY / periodicZ which will over-ride the above choices (not yet implemented)

    '''

    if dimensions == None:
        print "Please specify mesh dimensions - 2D or 3D !"
        return


    print " *  Setting boundary conditions: (vx, Fx, vy, Fy, vz, Fz) - "
    print "    Top    - ", top
    print "    Bottom - ", bottom
    print "    Left   - ", left
    print "    Right  - ", right
    if dimensions == 3:
        print "    Front  - ", front
        print "    Back   - ", back



    ## Note: This is not a component so we don't have a python wrapper to stuff this in the dictionary !!
    
    globalDict = _underworld.GetCurrentPythonDictionary()        
    globalDict["velocityBCs"] = { "type" : "CompositeVC" }
          
    ## Make a list of all the boundary conditions as dictionaries;

    noSlip = [ { "name": "vx", "type": "double", "value": 0.0 },
               { "name": "vy", "type": "double", "value": 0.0 } ]
    if dimensions == 3:
        noSlip.append( {"name": "vz", "type": double, "value": 0.0 } )

    freeslipXwall = [ { "name": "vx", "type": "double", "value": 0.0 } ]
    freeslipYwall = [ { "name": "vy", "type": "double", "value": 0.0 } ]
    freeslipZwall = [ { "name": "vz", "type": "double", "value": 0.0 } ]
    

    ## Periodic - set relevant wall boundary conditions to None and then set up as required


    mapBCs = [ "vx", "Fx", "vy", "Fy", "vz", "Fz" ]
    BCdict = { "top":top, "bottom":bottom, "left":left, "right":right, "front":front, "back":back }

    vcList = []

    for wall in [ "top", "bottom", "left", "right", "front", "back" ]:
        # Velocity degrees of freedom
        for dof in range(0,dimensions*2-1 ,2):  
            if BCdict[wall][dof] != None:
            #    print wall, " -> ", mapBCs[dof], " -> ", BCdict[wall][dof]        
                vcList.append( {"type":"WallVC", "wall": wall, "variables": { "name": mapBCs[dof], "type": "double", "value": BCdict[wall][dof] } } )




    globalDict["velocityBCs"]["vcList"] = vcList

    print vcList








