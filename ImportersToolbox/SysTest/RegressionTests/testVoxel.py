#!/usr/bin/env python

from credo.systest import *

testSuite = SysTestSuite("Underworld", "RegressionTests")

# voxel tests
for np in [1,2]:
	testSuite.addStdTest(ReferenceTest, "testVoxelGeomodeller.xml",
		nproc=np, expPathPrefix="expected", runSteps=0, fieldsToTest=["PressureField"])
	testSuite.addStdTest(ReferenceTest, "testVoxelASCII.xml",
		nproc=np, expPathPrefix="expected", runSteps=0, fieldsToTest=["PressureField"])
	testSuite.addStdTest(ReferenceTest, "testVoxelGMTPPC.xml",
		nproc=np, expPathPrefix="expected", runSteps=0, fieldsToTest=["VelocityField"])

# Set to a high time, as extension test can be quite slow on older procs
testSuite.setAllTimeouts(minutes=60)

def suite():
    return testSuite

if __name__ == "__main__":
    testRunner = SysTestRunner()
    testRunner.runSuite(testSuite)
