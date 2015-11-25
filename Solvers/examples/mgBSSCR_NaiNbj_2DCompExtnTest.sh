#! /bin/bash

#
#   Solve using BSSCR SLE augmented lagrangian method.
#   -ksp_k2_type SLE will result in BSSCR using user-defined assembled matrix for K2.
#   In the solver K <- K+penalty*K2 where K2 is integral(Na,i*Nb,j). ( user-defined matrix for this model )
#   A vector f2 is set up in xml that will contain boundary condition terms from K2.
#   This is added to the original force vector with penalty.
#   The default is to add f2 to forceterm. This can be changed by using -forcecorrection false.
#
#   K2 is assembled in Experimental/Rheology/src/Matrix_NaiNai.c;
#        set up via Experimental/InputFiles/MVModels/solvertest/_solvers/AugLagSchurCompSolverSLE-intNaiNai.xml
#
#   There exist two penalty terms for this configuration:
#   1. --components.MatrixTerm.incompressibility_Penalty passed directly into assembly routine
#   2. --components.stokesEqn.penaltyNumber              passed into BSSCR solver and used in K <- K + penalty*K2
#   The final effective penalty number will be the product of these two.
#
#   -augmented_lagrangian option will cause BSSCR to use function that is "augmented aware" but will act like regular schur solve if
#                         penalty is zero or -ksp_k2_type NULL.
#
#   -Q22_pc_type uw       option will cause BSSCR to use Schur preconditioner matrix that is set up in xml by default
#                         K used in preconditioner matrix is vanilla K.
#   -Q22_pc_type uwscale  option will cause BSSCR to use Schur preconditioner matrix created from G'*diag(K)^(-1)*G in BSSCR solver. 
#                         (option name misleading; will change)
#                         K used in preconditioner matrix is augmented K.
#
#   The K2 matrix for the GMG setup is identically proportional (4x) to the viscous penalty matrix from here and mgSchur2DCompExtnTest.sh.
#   i.e. G'*M^(-1)*G == 4.0*integral(Na,i*Nb,j) ( at least for Q1-P0 elements ) ( relation identical for 64x32 and 128x64 meshes )
#
#   With current configuration: (output identical to mgSchur2DCompExtnTest.sh and mgBSSCR2DCompExtnTest.sh both with -Q22_pc_type uwscale option)
#     Solution time:      9.086343e+01    (on Intel(R) Core(TM)2 Quad CPU    Q9650  @ 3.00GHz : in serial) (REP deactivated)
#     Picard residual:    9.93499269e-04
#     Picard Iterations:  123
#     Total Pressure its: 755
#
#   With -forcecorrection false: (output identical to mgBSSCR_GMG_2DCompExtnTest.sh with penalty=2.5)
#     Solution time:      8.740317e+01    ( as above )
#     Picard residual:    9.49525379e-04
#     Picard Iterations:  113
#     Total Pressure its: 669
#   With -forcecorrection false and -Q22_pc_type uw: (output identical to mgBSSCR_GMG_2DCompExtnTest.sh with penalty=2.5)
#     Solution time:      8.552165e+01
#     Picard residual:    9.29550173e-04
#     Picard Iterations:  98
#     Total Pressure its: 1023
#
#   With 128x64 mesh and -forcecorrection false: (output identical to mgBSSCR_GMG_2DCompExtnTest.sh with penalty=2.5 and 128x64 mesh)
#      Solution time:      2.866074e+02
#      Picard residual:    9.85464953e-04
#      Picard Iterations:  87
#      Total Pressure its: 712
#
export UWPATH=`./getUWD.sh`
export UWEXEC=$UWPATH/build/bin/Underworld
OUT="outputMGBSSCR_NaiNbj"
mkdir $OUT >& /dev/null
$UWEXEC $UWPATH/Underworld/InputFiles/+LithosphereSandbox/Sandbox2D.xml  $UWPATH/Underworld/InputFiles/Sandbox2DOverRide.xml \
    $UWPATH/Solvers/InputFiles/AugLagStokesSLE-NaiNbj.xml \
    $UWPATH/Solvers/InputFiles/kspinterface.xml \
    $UWPATH/Solvers/InputFiles/MultigridForRegularSCR.xml \
    -options_file ./options-scr-mg-accelerating.opt \
  		--outputPath="./$OUT" \
  		--components.stokesEqn.isNonLinear=True \
  		--components.box.startY=0.0 \
  		--components.box.endY=1.0 \
  		--components.weakZone.startY=0.0 --components.weakZone.endY=0.05 --components.substrateViscosity.eta0=0.1 \
  		--saveDataEvery=1 --checkpointEvery=2 --checkpointWritePath="./$OUT/Checkpoints" --checkpointAppendStep=1 \
  		--Cohesion=50 --Cohesion2=50 --FrictionCoeff=0.25 --FrictionCoeff2=0.25 \
                --components.stokesEqn.penaltyNumber=1.0 \
                --components.MatrixTerm.incompressibility_Penalty=10.0 \
		--components.MatrixTerm.viscosity_weighting=false \
  		-ksp_type bsscr -pc_type none -ksp_k2_type SLE -augmented_lagrangian \
  		-Q22_pc_type uwscale -Xrestore_K -no_scale \
  		-Xbuild_cb_pressure_nullspace \
  		-Xbuild_const_pressure_nullspace \
  		--mgLevels=5 \
  		-Xscr_ksp_max_it 1000 \
  		-scr_ksp_rtol 1.0e-4 \
  		-A11_ksp_rtol 1.0e-4 \
                -backsolveA11_ksp_type fgmres \
                -backsolveA11_ksp_rtol 1.0e-7 \
  		--minYieldCriterion=1.0 --minViscosity=1.0 \
                --nonLinearTolerance=0.001 \
  		--elementResI=64 --elementResJ=32 \
  		--maxTimeSteps=0 -Xdump_matvec -forcecorrection false\
                < /dev/null 2> "./$OUT/errors.txt" | tee "./$OUT/output.txt"
 	
./getconv.pl < "$OUT/output.txt"

# Mutually incompatible at present: 
# -build_cb_pressure_nullspace  
# -build_const_pressure_nullspace


# Choices
# -ksp_k2_type NULL GG GMG DGMGD SLE


