#!/usr/bin/env python

from credo.systest import *

testSuite = SysTestSuite("Spherical", "RegressionTests")

# testing spherical-segment mesh
for np in [1,2]:
        testSuite.addStdTest(ReferenceTest, "quasi_annulus.xml", runSteps=1, nproc=np)

testSuite.addStdTest(ReferenceTest, "RTP3D_pic.xml", runSteps=1, nproc=1)
testSuite.addStdTest(ReferenceTest, "RS_pic.xml", runSteps=1, nproc=1)
testSuite.addStdTest(ReferenceTest, "RS_fem_lidDriven.xml", runSteps=1, nproc=1)
# Set to a high time, as extension test can be quite slow on older procs
testSuite.setAllTimeouts(minutes=60)

def suite():
    return testSuite

if __name__ == "__main__":
    testRunner = SysTestRunner()
    testRunner.runSuite(testSuite)
