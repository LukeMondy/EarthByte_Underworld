#!/bin/sh
cd /Users/julian/codes/underworld2/libUnderworld/Underworld/SysTest/RegressionTests
mpiexec -np 1 /Users/julian/codes/underworld2/libUnderworld/build/bin/StGermain PoiseuilleFlow.xml /Users/julian/codes/underworld2/libUnderworld/Underworld/SysTest/RegressionTests/output/PoiseuilleFlow-referenceTest-np1/testRun/credo-analysis.xml 
