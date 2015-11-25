#! /bin/bash

#
#   Solve using Uzawa+MG.
#

export UWPATH=`./getUWD.sh`
export UWEXEC="cgdb --args $UWPATH/build/bin/Underworld"
#export UWEXEC="$UWPATH/build/bin/Underworld"
OUT="bsscrLid"
RES=16
mkdir $OUT >& /dev/null

#    -options_file ./options-scr-mg-accelerating.opt \

$UWEXEC $UWPATH/Solvers/InputFiles/liddrivenperiodicx.xml \
    $UWPATH/Solvers/InputFiles/kspinterface.xml \
    $UWPATH/Solvers/InputFiles/MultigridForRegularSCR.xml \
    -options_file ./options-scr-mg-accelerating.opt \
               --components.stokesEqn.isNonLinear=False \
                --mgLevels=3 \
                --dim=3 \
  		--outputPath="./$OUT" \
                -ksp_type bsscr -pc_type none -ksp_k2_type NULL -augmented_lagrangian --penaltyNumber=0.0 \
  		--elementResI=$RES --elementResJ=$RES --elementResK=$RES \
  		--maxTimeSteps=0 -Xdump_matvec -Xhelp \
                -A11_ksp_type richardson -XA11_ksp_converged_reason \
                -Xscr_ksp_monitor_true_residual -Xscr_ksp_norm_type unpreconditioned \
                -scr_ksp_view -help \
  		--saveDataEvery=1 --checkpointEvery=1 --checkpointWritePath="./$OUT/Checkpoints" --checkpointAppendStep=1 \

#

#                < /dev/null 2> "./$OUT/errors.txt" | tee "./$OUT/output.txt"
#
#   Sizes of matrices dumped to files from PETScMGSolver_UpdateMatrices function (line 491 $UWPATH/StgFEM/SLE/SystemSetup/src/PETScMGSolver.c)
#   for MG with 3 levels keeping BC's in. (velocity)
#   for 16^3 mesh mg level 2-> 14739 unknowns  <-- 17^3 ( 17 nodes for 16 elements)
#   for 16^3 mesh mg level 1-> 2187 unknowns   <--  9^3 (  9 nodes for 8 elements )
#   for 16^3 mesh mg level 0-> 375 unknowns    <--  5^3 (  5 nodes for 4 elements )

#   for MG with 3 levels taking BC's out. (velocity)
#   for 16^3 mesh mg level 2-> 12495 unknowns  <-- 17^3-bc's ( 17 nodes for 16 elements)
#   for 16^3 mesh mg level 1-> 2187 unknowns   <--  9^3 (  9 nodes for 8 elements )
#   for 16^3 mesh mg level 0-> 375 unknowns    <--  5^3 (  5 nodes for 4 elements )

# using  <include>Underworld/VariableConditions/velocityBC.topLidSinusoidal_NoSlipBottom.xml</include>
# with vx vy set on top and bottom walls -> 2*17^2*2 = 1156
# vz set on front and back walls         -> 2*17*17  =  578
# vx set on left and right walls         -> 2*15*17   = 510
#                                               total  2244
#                                          unknowns = 17^3 - 2244 = 12495
#
#   Direct solve via following options shows velocity system is consistent (matrix has full rank)
#    -Uzawa_velSolver_ksp_type preonly -Uzawa_velSolver_pc_type lu -Uzawa_velSolver_ksp_view -Uzawa_velSolver_ksp_converged_reason \
#   sub-matrices also have full rank via rank calculation in "octave".
#