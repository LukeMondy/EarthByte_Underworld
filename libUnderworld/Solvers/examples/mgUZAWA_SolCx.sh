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

#        $UWPATH/StgFEM/Apps/StgFEM_Components/MultigridForRegular.xml \
#        -options_file $UWPATH/Underworld/InputFiles/options/options-uzawa-mg.opt \
#$UWEXEC $UWPATH/Solvers/InputFiles/testVelicSolCx.xml \

export UWPATH=`./getUWD.sh`
export UWEXEC="mpirun -n 2 $UWPATH/build/bin/Underworld"
export UWEXEC=$UWPATH/build/bin/Underworld
#export UWEXEC="cgdb --args $UWPATH/build/bin/Underworld"

OUT="outputMGUZAWA_solCx"
mkdir $OUT >& /dev/null

$UWEXEC $UWPATH/Underworld/SysTest/PerformanceTests/testVelicSolCx.xml \
        $UWPATH/Solvers/InputFiles/SolversToolBox.xml \
        $UWPATH/Solvers/InputFiles/analyticVis.xml \
        -options_file $UWPATH/Underworld/InputFiles/options/uzawa-mumps-3.0.opt \
  		--outputPath="./$OUT" \
  		--components.stokesEqn.isNonLinear=False \
  		--saveDataEvery=1 --checkpointEvery=1 --checkpointWritePath="./$OUT/Checkpoints" --checkpointAppendStep=1 \
                --components.FieldTest.normaliseByAnalyticSolution=False \
                --solCx_etaA=10 \
                --solCx_xc=0.6 \
                --solCx_etaB=1 \
                --solCx_n=2.0 \
                --wavenumberY=2.0 \
                --nonLinearTolerance=0.001 \
  		--elementResI=32 --elementResJ=32 \
                --mgLevels=3 \
  		--maxTimeSteps=0 -dump_matvec \
                -Uzawa_velSolver_ksp_view -help \
#                < /dev/null 2> "./$OUT/errors.txt" | tee "./$OUT/output.txt"

#
#
