#! /bin/bash

#
#   Solve using Uzawa+MG.
#

export UWPATH=`./getUWD.sh`
export UWEXEC="cgdb --args $UWPATH/build/bin/Underworld"
export UWEXEC="$UWPATH/build/bin/Underworld"
OUT="uzawaLidRamp"
RES=16
mkdir $OUT >& /dev/null

SOLVE=LU
SOLVE=MG
if [ "$SOLVE" = LU ]
then
    CONF="-Uzawa_velSolver_ksp_type preonly -Uzawa_velSolver_pc_type lu -Uzawa_velSolver_ksp_view -Uzawa_velSolver_ksp_converged_reason"
else
    CONF="$UWPATH/Experimental/InputFiles/Experimental_Components/MultigridForRegularUzawa.xml"
fi

$UWEXEC $UWPATH/Solvers/InputFiles/liddrivenperiodicxz.xml \
               $CONF \
               --components.stokesEqn.isNonLinear=False \
                --mgLevels=3 \
                --dim=3 \
  		--outputPath="./$OUT" \
  		--elementResI=$RES --elementResJ=$RES --elementResK=$RES \
  		--maxTimeSteps=0 -dump_matvec -Xhelp \
                -XA11_ksp_type richardson -XA11_ksp_converged_reason \
                --components.uzawa.monitor=true --components.uzawa.maxIterations=1000 \
  		--saveDataEvery=1 --checkpointEvery=1 --checkpointWritePath="./$OUT/Checkpoints" --checkpointAppendStep=1 \

#

#                < /dev/null 2> "./$OUT/errors.txt" | tee "./$OUT/output.txt"
#   b StgFEM/SLE/SystemSetup/src/PETScMGSolver.c:532
#   Sizes of matrices dumped to files from PETScMGSolver_UpdateMatrices function (line 491 $UWPATH/StgFEM/SLE/SystemSetup/src/PETScMGSolver.c)
#
#   Sinusoidal lid
#   for MG with 3 levels taking BC's out. (velocity)
#   for 8^3 mesh periodic in x and z
#   rank/size(A2)=1344/1344 (8*8*9-8*8*2)*3   svd max=0.978772 svd min=0.019030
#   rank/size(A1)=375/375  (5^3)*3   svd max=3.4873979 svd min=0.0029633
#   rank/size(A0)=81/81   (3^3)*3   svd max=6.45678 svd min=0.12870

#   Velocity_Lid_RampWithCentralMax with periodic BCs in x and z
#   13056 = 17*16*16*3
#   for 16^3 mesh periodic in x and z (sinusoidal lid and ramp lid)
#   size(A2)=11520  (=13056 - 1536) = 3*(17*16*16-16*16*2)
#   size(A1)=2187   (=3*9^3 curious what MG is doing here..as if periodic bc's are included in matrix)
#   size(A0)=375    (=3*5^3)
#

# using  Velocity_Lid_RampWithCentralMax
# periodic in x and z reduces nodes to 16 and 16 in x and z directions
# with vx vy vz set on top and bottom walls -> 3*16^2*2 = 1536
#
#
#   Direct solve via following options shows velocity system is consistent (matrix has full rank)
#    -Uzawa_velSolver_ksp_type preonly -Uzawa_velSolver_pc_type lu -Uzawa_velSolver_ksp_view -Uzawa_velSolver_ksp_converged_reason \
#   sub-matrices also have full rank via rank calculation in "octave".
#
#  Solves ok, so MG seems fine with periodic BC's here
#
# with LU
#  Uzawa its. = 0019 , Uzawa residual = 9.2990381545207e-06
#  |G^T u|/|u|               = 1.91484126e-04
#  |f - K u - G p|/|f|       = 4.80277171e-16
#  |f - K u - G p|_w/|f|_w   = 1.00000000e+00
#  |u|_{\infty} = 1.00000000e+00 , u_rms = 1.78534620e-01
#  |p|_{\infty} = 4.90057936e+00 , p_rms = 1.22868121e+00
#  min/max(u) = -1.47743588e-01 [12133] / 1.00000000e+00 [792]
#  min/max(p) = -4.90057936e+00 [1271] / 4.89703303e+00 [1791]
#  \sum_i p_i = -5.14149903e-04
#
# with MG
#  Uzawa its. = 0019 , Uzawa residual = 9.2990027179255e-06
#  |G^T u|/|u|               = 1.91484342e-04
#  |f - K u - G p|/|f|       = 2.81157852e-07
#  |f - K u - G p|_w/|f|_w   = 1.00000000e+00
#  |u|_{\infty} = 9.99999997e-01 , u_rms = 1.78534419e-01
#  |p|_{\infty} = 4.90058009e+00 , p_rms = 1.22868272e+00
#  min/max(u) = -1.47743519e-01 [12133] / 9.99999997e-01 [792]
#  min/max(p) = -4.90058009e+00 [1271] / 4.89703057e+00 [1791]
#  \sum_i p_i = -5.14060023e-04
#
