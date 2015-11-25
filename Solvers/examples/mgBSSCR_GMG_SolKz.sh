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
#
export UWPATH=`./getUWD.sh`
export UWEXEC=$UWPATH/build/bin/Underworld
OUT="outputMGBSSCR_GMG_SolKz"
mkdir $OUT >& /dev/null
$UWEXEC $UWPATH/Underworld/SysTest/PerformanceTests/testVelicSolKz.xml \
    $UWPATH/Solvers/InputFiles/AugLagStokesSLE-GtMG.xml \
    $UWPATH/Solvers/InputFiles/VelocityMassMatrixSLE.xml \
    $UWPATH/Solvers/InputFiles/kspinterface.xml \
    $UWPATH/Solvers/InputFiles/MultigridForRegularSCR.xml \
    $UWPATH/Solvers/InputFiles/analyticVis.xml \
    -options_file ./options-scr-mg-accelerating.opt \
  		--outputPath="./$OUT" \
  		--components.stokesEqn.isNonLinear=False \
  		--saveDataEvery=1 --checkpointEvery=2 --checkpointWritePath="./$OUT/Checkpoints" --checkpointAppendStep=1 \
		--components.constitutiveMatrix.viscosity_weighting=false \
                --components.FieldTest.normaliseByAnalyticSolution=False \
  		-ksp_type bsscr -pc_type none -ksp_k2_type GMG -augmented_lagrangian --penaltyNumber=2.5 \
  		-XQ22_pc_type uwscale \
  		-Q22_pc_type gtkg -Xrestore_K -no_scale -scr_pc_gtkg_ksp_view -scr_pc_gtkg_ksp_rtol 1e-2 -scr_pc_gtkg_ksp_type cg \
  		-Xbuild_cb_pressure_nullspace \
  		-Xbuild_const_pressure_nullspace \
  		--mgLevels=5 \
  		-Xscr_ksp_max_it 1000 \
  		-scr_ksp_rtol 1.0e-4 \
  		-A11_ksp_rtol 1.0e-4 \
                -backsolveA11_ksp_type fgmres \
                -backsolveA11_ksp_rtol 1.0e-7 \
  		--elementResI=128 --elementResJ=128 \
  		--maxTimeSteps=0 -Xdump_matvec \
                < /dev/null 2> "./$OUT/errors.txt" | tee "./$OUT/output.txt"
 		
./getconv.pl < "$OUT/output.txt"

#> "./outputMGBSSCR_GMG/output.txt" 2>&1	

# Mutually incompatible at present: 
# -build_cb_pressure_nullspace  
# -build_const_pressure_nullspace


# Choices
# -ksp_k2_type NULL GG GMG DGMGD SLE
