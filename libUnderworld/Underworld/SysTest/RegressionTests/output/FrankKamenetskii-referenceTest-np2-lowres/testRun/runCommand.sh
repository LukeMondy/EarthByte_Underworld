#!/bin/sh
cd /Users/julian/codes/underworld2/libUnderworld/Underworld/SysTest/RegressionTests
mpiexec -np 2 /Users/julian/codes/underworld2/libUnderworld/build/bin/StGermain FrankKamenetskii.xml /Users/julian/codes/underworld2/libUnderworld/Underworld/SysTest/RegressionTests/output/FrankKamenetskii-referenceTest-np2-lowres/testRun/credo-analysis.xml --elementResI=10 --elementResJ=10 --elementResK=10
