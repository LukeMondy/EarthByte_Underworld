#!/usr/bin/env python

from credo.systest import *

testSuite = SysTestSuite("Spherical", "RegressionTests")

models = ["LidDrivenFEM.xml",
          "LidDrivenPIC.xml",
          "Periodic.annulus.xml" ]

# testing spherical-segment mesh
for np in [1,2]:
  testSuite.addStdTest(ReferenceTest, "LidDrivenFEM.xml", runSteps=0, nproc=np)

testSuite.addStdTest(ReferenceTest, "LidDrivenPIC.xml", runSteps=5, nproc=1)
testSuite.addStdTest(ReferenceTest, "Periodic.annulus.xml", runSteps=1, nproc=1)

# Set to a high time, as extension test can be quite slow on older procs
testSuite.setAllTimeouts(minutes=60)

def suite():
    return testSuite

if __name__ == "__main__":
    testRunner = SysTestRunner()
    testRunner.runSuite(testSuite)
