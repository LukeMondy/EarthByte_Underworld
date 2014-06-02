#!/bin/sh
cd /Users/julian/codes/underworld2/libUnderworld/Underworld/SysTest/RegressionTests
mpiexec -np 2 /Users/julian/codes/underworld2/libUnderworld/build/bin/StGermain testSLAdvectionDiffusion.xml /Users/julian/codes/underworld2/libUnderworld/Underworld/SysTest/RegressionTests/output/testSLAdvectionDiffusion-analyticTest-np2-elementResI-10-elementResJ-10-elementResK-10/testRun/credo-analysis.xml --elementResI=10 --elementResJ=10 --elementResK=10
