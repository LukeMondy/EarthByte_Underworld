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
#    -options_file ./options-scr-mg-accelerating.opt \

export UWPATH=`./getUWD.sh`
export UWEXEC=$UWPATH/build/bin/Underworld
OUT="outputMGBSSCR_GMG_SolCx-test"
mkdir $OUT >& /dev/null


$UWEXEC $UWPATH/Underworld/SysTest/PerformanceTests/testVelicSolCx.xml \
    $UWPATH/Solvers/InputFiles/AugLagStokesSLE-GtMG.xml \
    $UWPATH/Solvers/InputFiles/kspinterface.xml \
    $UWPATH/Solvers/InputFiles/MultigridForRegularSCR.xml \
    $UWPATH/Solvers/InputFiles/analyticVis.xml \
    -options_file options-scr-mg-accelerating.opt \
  		--outputPath="$OUT" \
  		--components.stokesEqn.isNonLinear=True \
  		--saveDataEvery=1 --checkpointEvery=2 --checkpointWritePath="./$OUT/Checkpoints" --checkpointAppendStep=1 \
      --components.FieldTest.normaliseByAnalyticSolution=False \
      --solCx_etaA=10000000 \
      --solCx_xc=0.6 \
      --solCx_etaB=1 \
      --solCx_n=2.0 \
      --wavenumberY=2.0 \
      -uzawastyle 0 \
      -scrPCKSP_ksp_type fgmres \
      -scrPCKSP_ksp_converged_reason \
      -scrPCKSP_ksp_view \
  		-ksp_type bsscr -pc_type none -ksp_k2_type GMG -augmented_lagrangian --penaltyNumber=0.0 \
  		-Q22_pc_type uw \
  		-XQ22_pc_type gtkg -Xrestore_K -no_scale \
      -Xscr_pc_gtkg_ksp_view -Xscr_pc_gtkg_ksp_monitor -Xscr_pc_gtkg_ksp_rtol 1e-6 -Xscr_pc_gtkg_ksp_type cg \
  		-Xbuild_cb_pressure_nullspace \
  		-build_const_pressure_nullspace \
  		--mgLevels=2 \
  		-Xscr_ksp_max_it 1000 \
      -scr_ksp_type fgmres \
  		-scr_ksp_rtol 1.0e-6 \
  		-A11_ksp_rtol 1.0e-8 \
      -backsolveA11_ksp_type gmres \
      -backsolveA11_ksp_rtol 1.0e-7 \
  		--elementResI=128 --elementResJ=128 \
  		--maxTimeSteps=10 -Xdump_matvec \
          < /dev/null 2> "./$OUT/errors.txt" | tee "./$OUT/output.txt"

./getconv.pl < "$OUT/output.txt"

#> "./outputMGBSSCR_GMG/output.txt" 2>&1	

# Mutually incompatible at present: 
# -build_cb_pressure_nullspace  
# -build_const_pressure_nullspace


# Choices
# -ksp_k2_type NULL GG GMG DGMGD SLE
