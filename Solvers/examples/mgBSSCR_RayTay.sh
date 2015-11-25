#! /bin/bash

#
#   Solve using BSSCR GG/GMG augmented lagrangian method (optionally).
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
#
#   It seems that this model works better without the augmented lagrangian.
#   It also converges very well with just the standard Underworld Schur compliment preconditioner (-Q22_pc_type uw)
#   For a 128x128 mesh it converges in 2 iterations for the pressure ( -scr_ksp_rtol 1.0e-3 -A11_ksp_rtol 1.0e-2 ) and
#   3 final velocity iterations for the backsolve ( -backsolveA11_ksp_rtol 1.0e-7 -backsolveA11_ksp_type gmres ).
#$UWEXEC $UWPATH/Underworld/SysTest/PerformanceTests/testVelicSolCx.xml \

export UWPATH=`./getUWD.sh`
export UWEXEC="mpirun -n 4 $UWPATH/build/bin/Underworld"
OUT="rt_64x64afterVecDest_MGBSSCR_GMG_RayTay"
mkdir $OUT >& /dev/null
$UWEXEC $UWPATH/Underworld/InputFiles/RayleighTaylorBenchmark.xml \
    $UWPATH/Solvers/InputFiles/AugLagStokesSLE-GtMG.xml \
    $UWPATH/Solvers/InputFiles/kspinterface.xml \
    $UWPATH/Solvers/InputFiles/MultigridForRegularSCR.xml \
    -options_file ./options-scr-mg-accelerating.opt \
  		--outputPath="./$OUT" \
  		--components.stokesEqn.isNonLinear=False \
  		--saveDataEvery=5 --checkpointEvery=5 --checkpointWritePath="./$OUT/Checkpoints" --checkpointReadPath="./$OUT/Checkpoints" --checkpointAppendStep=1 \
  		-ksp_type bsscr -pc_type none -ksp_k2_type GMG -augmented_lagrangian 1 --penaltyNumber=0.1 \
  		-Q22_pc_type gkgdiag \
  		-XQ22_pc_type gtkg -Xrestore_K -no_scale \
                -Xscr_pc_gtkg_ksp_view -Xscr_pc_gtkg_ksp_monitor -Xscr_pc_gtkg_ksp_rtol 1e-6 -Xscr_pc_gtkg_ksp_type cg \
  		-remove_checkerboard_pressure_null_space 0 \
  		-remove_constant_pressure_null_space 1 \
  		--mgLevels=5 \
  		-Xscr_ksp_max_it 1000 \
                -scr_ksp_type minres \
  		-scr_ksp_rtol 1.0e-4 \
  		-A11_ksp_rtol 1.0e-5 \
                -backsolveA11_ksp_type gmres \
                -backsolveA11_ksp_rtol 1.0e-5 \
  		--elementResI=64 --elementResJ=64 \
 		--maxTimeSteps=100 -Xdump_matvec -log_summary \
                < /dev/null 2> "./$OUT/errors.txt" | tee "./$OUT/output.txt"

./getconv.pl < "$OUT/output.txt"

#> "./outputMGBSSCR_GMG/output.txt" 2>&1	

#                --restartTimestep=1000 \ 


# Mutually incompatible at present: 
# -build_cb_pressure_nullspace  
# -build_const_pressure_nullspace


# Choices
# -ksp_k2_type NULL GG GMG DGMGD SLE
