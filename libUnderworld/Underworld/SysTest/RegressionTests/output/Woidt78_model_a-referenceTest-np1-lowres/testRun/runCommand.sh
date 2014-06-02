#!/bin/sh
cd /Users/julian/codes/underworld2/libUnderworld/Underworld/SysTest/RegressionTests
mpiexec -np 1 /Users/julian/codes/underworld2/libUnderworld/build/bin/StGermain Woidt78_model_a.xml /Users/julian/codes/underworld2/libUnderworld/Underworld/SysTest/RegressionTests/output/Woidt78_model_a-referenceTest-np1-lowres/testRun/credo-analysis.xml --elementResI=10 --elementResJ=25 --elementResK=10
