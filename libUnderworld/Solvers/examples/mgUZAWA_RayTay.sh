#! /bin/bash

#
#   Solve using Uzawa + MG
#
#   Uzawa Solver has two ksps
#   1: the preconditioner ksp with prefix -Uzawa_pcSolver_
#   2: the velocity ksp with prefix       -Uzawa_velSolver_
#   When Multigrid is used the MG functions replace the velocity ksp
#   and set it's prefix to -A11       ( line 308 in /StgFEM/SLE/SystemSetup/src/PETScMGSolver.c )
#   The preconditioner ksp is unchanged.
#
export UWPATH=`./getUWD.sh`
export UWEXEC="mpirun -n 4 $UWPATH/build/bin/Underworld"
OUT="testpng"
mkdir $OUT >& /dev/null
$UWEXEC $UWPATH/Underworld/InputFiles/RayleighTaylorBenchmark.xml \
    $UWPATH/StgFEM/Apps/StgFEM_Components/MultigridForRegular.xml \
    -options_file $UWPATH/Underworld/InputFiles/options/options-uzawa-mg.opt \
  		--outputPath="./$OUT" \
  		--components.stokesEqn.isNonLinear=False \
  		--saveDataEvery=50 --checkpointEvery=50 --checkpointWritePath="./$OUT/Checkpoints" --checkpointReadPath="./$OUT/Checkpoints" --checkpointAppendStep=1 \
  		--elementResI=128 --elementResJ=128 \
                --mgLevels=5 \
  		--maxTimeSteps=1000 -log_summary \
                < /dev/null 2> "./$OUT/errors.txt" | tee "./$OUT/output.txt"

#
#