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
export UWEXEC=$UWPATH/build/bin/Underworld
OUT="globe2"
mkdir $OUT >& /dev/null

#    $UWPATH/StgFEM/Apps/StgFEM_Components/MultigridForRegular.xml \
#    -options_file $UWPATH/Underworld/InputFiles/options/options-uzawa-mg.opt \

mpirun -n 2 $UWEXEC ./gplateVoxelGMTGlobe.xml \
    --outputPath="$OUT" \
    --particlesPerCell=30 \
    --averageInitialParticlesPerCell=30 \
    --components.weights.Inflow=True \
    --components.weights.maxSplits=25 \
    --components.uzawa.monitor=1 \
    --mgLevels=3 \
    -A11_ksp_rtol 1.0e-3 \
    --elementResLong=24 --elementResLat=16 --elementResR=12 \
    --maxTimeSteps=1 -Xdump_matvec \
    < /dev/null 2> "./$OUT/errors.txt" | tee "./$OUT/output.txt"



