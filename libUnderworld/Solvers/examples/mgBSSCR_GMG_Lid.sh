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
#

export UWPATH=`./getUWD.sh`
export UWEXEC="cgdb --args $UWPATH/build/bin/Underworld"
#export UWEXEC="$UWPATH/build/bin/Underworld"
OUT="outputMGBSSCR_GMG_SolCx"
mkdir $OUT >& /dev/null
$UWEXEC $UWPATH/Solvers/InputFiles/liddriven.xml \
    $UWPATH/Solvers/InputFiles/AugLagStokesSLE-GtMG.xml \
    $UWPATH/Solvers/InputFiles/VelocityMassMatrixSLE.xml \
    $UWPATH/Solvers/InputFiles/kspinterface.xml \
    $UWPATH/Solvers/InputFiles/quiet.xml \
                --dim=3 \
  		--outputPath="./$OUT" \
  		--components.stokesEqn.isNonLinear=False \
  		--saveDataEvery=1 --checkpointEvery=2 --checkpointWritePath="./$OUT/Checkpoints" --checkpointAppendStep=1 \
                -uzawastyle 0 \
                -scrPCKSP_ksp_type fgmres \
                -XscrPCKSP_ksp_converged_reason \
                -XscrPCKSP_ksp_view \
  		-ksp_type bsscr -pc_type none -ksp_k2_type DGMGD -augmented_lagrangian --penaltyNumber=0.05 \
  		-Q22_pc_type uw \
  		-XQ22_pc_type gtkg -Xrestore_K -no_scale \
                -Xscr_pc_gtkg_ksp_view -Xscr_pc_gtkg_ksp_monitor -scr_pc_gtkg_ksp_rtol 1e-6 -scr_pc_gtkg_ksp_type cg \
  		-Xbuild_cb_pressure_nullspace \
  		-build_const_pressure_nullspace \
  		--mgLevels=5 \
  		-Xscr_ksp_max_it 1000 \
                -scr_ksp_type fgmres \
                -Xscr_ksp_norm_type preconditioned \
                -Xscr_ksp_left_pc \
  		-scr_ksp_rtol 1.0e-8 \
                -scr_ksp_monitor_true_residual \
              -A11_pc_type hypre -A11_pc_hypre_type boomeramg -XA11_pc_hypre_boomeramg_print_statistics \
  		-A11_ksp_rtol 1.0e-6 \
                -A11_ksp_type gmres \
              -XA11_pc_hypre_boomeramg_grid_sweeps_all 5 \
              -XA11_pc_hypre_boomeramg_tol 1e-3 \
                -A11_ksp_converged_reason \
               -XA11_ksp_norm_inf_monitor \
               -XA11_use_norm_inf_stopping_condition \
               -XA11_ksp_monitor_true_residual \
                -XA11_ksp_view \
                -backsolveA11_ksp_type fgmres \
                -backsolveA11_ksp_rtol 1.0e-6 \
  		--elementResI=16 --elementResJ=16 --elementResK=16\
  		--maxTimeSteps=0 -Xdump_matvec -Xhelp \
               < /dev/null 2> "./$OUT/errors.txt" | tee "./$OUT/output.txt"
./getconv2.pl < "$OUT/output.txt"

#                -Xscr_ksp_norm_type natural \ <-- not supported when using KSPRelativeRhsConvergenced()

#    
#
#    $UWPATH/Solvers/InputFiles/MultigridForRegularSCR.xml \
#    -options_file ./options-scr-mg-accelerating.opt \

# Mutually incompatible at present: 
# -build_cb_pressure_nullspace  
# -build_const_pressure_nullspace


# Choices
# -ksp_k2_type NULL GG GMG DGMGD SLE
