#! /bin/bash

#
#   Solve using BSSCR GG/GMG augmented lagrangian method.
#   In the solver K <- K+penalty*K2 where K2 is G'*M^(-1)*G.
#   The M matrix is setup and attached to an SLE in Experimental/InputFiles/MVModels/solvertest/_solvers/AugLagSchurCompSolverSLE-DGTGD.xml.
#   K2 is calculated in the BSSCR solver.
#
#   -augmented_lagrangian option will cause BSSCR to use function that is "augmented aware" but will act like regular schur solve if
#                         penalty is zero or -ksp_k2_type NULL.
#
#   -Q22_pc_type uw       option will cause BSSCR to use Schur preconditioner matrix that is set up in xml by default
#                         K used in preconditioner matrix is vanilla K.
#   -Q22_pc_type uwscale  option will cause BSSCR to use Schur preconditioner matrix created from G'*diag(K)^(-1)*G in BSSCR solver.
#                         (option name misleading; will change)
#                         The preconditioner matrix will be based on the augmented lagrangian K matrix.
#
#   The K2 matrix for this setup is identically proportional (4x) to the viscous penalty matrix from mgBSSCR_NaiNai_2DCompExtnTest.sh and mgSchur2DCompExtnTest.sh.
#   i.e. G'*M^(-1)*G == 4.0*integral(Na,i*Nb,j) ( at least for Q1-P0 elements ) ( relation identical for 64x32 and 128x64 meshes )
#
#
#   With 64x32 mesh -Q22_pc_type uwscale: (output identical to mgBSSCR_NaiNai_2DCompExtnTest.sh with -forcecorrection false) (note that penalty=2.5 here)
#      Solution time:      7.666981e+01  (on Intel(R) Core(TM)2 Quad CPU    Q9650  @ 3.00GHz : in serial) (REP deactivated)
#      Picard residual:    9.49525379e-04
#      Picard Iterations:  113
#      Total Pressure its: 669
#   With 128x64 mesh -Q22_pc_type uwscale: (output identical to mgBSSCR_NaiNai_2DCompExtnTest.sh with -forcecorrection false and 128x64 mesh) (note that penalty=2.5 here)
#      Solution time:      2.457116e+02
#      Picard residual:    9.85464953e-04
#      Picard Iterations:  87
#      Total Pressure its: 712
#
#   With -Q22_pc_type uw: (output identical to mgBSSCR_NaiNai_2DCompExtnTest.sh with -forcecorrection false and -Q22_pc_type uw) (note that penalty=2.5 here)
#      Solution time:      7.471834e+01
#      Picard residual:    9.29550173e-04
#      Picard Iterations:  98
#      Total Pressure its: 1023
#
#   With -Q22_pc_type gtkg (with Velocity Mass Matrix) 64x32 resolution.
#      Solution time:      6.807188e+01
#      Picard residual:    9.43188306e-04
#      Picard Iterations:  106
#      Total Pressure its: 274
#   With -Q22_pc_type gtkg (with Velocity Mass Matrix) 128x64 resolution. 
#      Solution time:      2.746766e+02
#      Picard residual:    9.80036554e-04
#      Picard Iterations:  83
#      Total Pressure its: 518
#
#
export UWPATH=`./getUWD.sh`
export UWEXEC=$UWPATH/build/bin/Underworld
OUT="outputMGBSSCR_GMG"
mkdir $OUT >& /dev/null
$UWEXEC $UWPATH/Underworld/InputFiles/+LithosphereSandbox/Sandbox2D.xml $UWPATH/Underworld/InputFiles/Sandbox2DOverRide.xml \
    $UWPATH/Solvers/InputFiles/AugLagStokesSLE-GtMG.xml \
    $UWPATH/Solvers/InputFiles/VelocityMassMatrixSLE.xml \
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
		--components.constitutiveMatrix.viscosity_weighting=false \
  		-ksp_type bsscr -pc_type none -ksp_k2_type GMG -augmented_lagrangian --penaltyNumber=2.5 \
  		-XQ22_pc_type uwscale \
  		-Q22_pc_type uwscale -Xrestore_K -no_scale -Xscr_pc_gtkg_ksp_view -Xscr_pc_gtkg_ksp_rtol 1e-2 -Xscr_pc_gtkg_ksp_type cg \
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
  		--maxTimeSteps=0 -Xdump_matvec \
                < /dev/null 2> "./$OUT/errors.txt" | tee "./$OUT/output.txt"
 		
./getconv.pl < "$OUT/output.txt"

#> "./outputMGBSSCR_GMG/output.txt" 2>&1	

# Mutually incompatible at present: 
# -build_cb_pressure_nullspace  
# -build_const_pressure_nullspace


# Choices
# -ksp_k2_type NULL GG GMG DGMGD SLE
