#!/bin/sh
cd /Users/julian/codes/underworld2/libUnderworld/Underworld/SysTest/RegressionTests
mpiexec -np 1 /Users/julian/codes/underworld2/libUnderworld/build/bin/StGermain Nightmare.xml /Users/julian/codes/underworld2/libUnderworld/Underworld/SysTest/RegressionTests/output/Nightmare-referenceTest-np1-lowres/testRun/credo-analysis.xml --elementResI=10 --elementResJ=10 --elementResK=10