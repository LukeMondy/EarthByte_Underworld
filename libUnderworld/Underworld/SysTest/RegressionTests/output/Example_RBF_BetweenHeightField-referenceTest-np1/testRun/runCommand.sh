#!/bin/sh
cd /Users/julian/codes/underworld2/libUnderworld/Underworld/SysTest/RegressionTests
mpiexec -np 1 /Users/julian/codes/underworld2/libUnderworld/build/bin/StGermain Example_RBF_BetweenHeightField.xml /Users/julian/codes/underworld2/libUnderworld/Underworld/SysTest/RegressionTests/output/Example_RBF_BetweenHeightField-referenceTest-np1/testRun/credo-analysis.xml 
