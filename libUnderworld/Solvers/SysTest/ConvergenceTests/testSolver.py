#!/usr/bin/env python
import os
import copy
from credo.systest import *

resSet2D = [(32,32),(48,48),(96,96)] # if only one res set here then we get a glorious divide by zero error in credo.

solvSuite = SysTestSuite("Solvers", "RegressionTests")

models = ["testVelicSolCx.xml","testVelicSolKx.xml"]
solcxOpts = {"mgLevels":3,"elementResI":64, "elementResJ":64,"components.FieldTest.normaliseByAnalyticSolution":False,"particlesPerCell":40,"penaltyNumber":0.0,"solCx_etaA":1000.0,"solCx_etaB":1.0,"solCx_n":2.0,"wavenumberY":2.0}
solkxOpts = {"mgLevels":3,"elementResI":64, "elementResJ":64,"components.FieldTest.normaliseByAnalyticSolution":False,"particlesPerCell":40,"penaltyNumber":0.25,"solKx_m":2.0,"solKx_n":2.0,"solKx_sigma":1.0,"solKx_twiceB":6.90775527898213705205,"wavenumberX":2.0,"wavenumberY":2.0}

optionsList=[]

optionsList.append(solcxOpts)
optionsList.append(solkxOpts)
num=0
for modelXML in models:
    for np in [1,2,4]:
        solvSuite.addStdTest(AnalyticMultiResTest, modelXML, nproc=np, paramOverrides=optionsList[num], resSet=resSet2D)
    num=num+1

solvSuite.setAllTimeouts(minutes=2)

mgSetupXML = "Solvers/MultigridForRegularSCR.xml"
alSetupXML = "Solvers/AugLagStokesSLE-GtMG.xml"
vmSetupXML = "Solvers/VelocityMassMatrixSLE.xml"
kiSetupXML = "Solvers/kspinterface.xml"

mgSolverOpts = "options-solv-mg.opt"

num=0 # Lets fix credo's ridiculously long path names
for sysTest in solvSuite.sysTests:
    path="output/test-"+str(num)+"-solv"
    test="test-"+str(num)+"-solv"
    sysTest.testName=test
    #sysTest.updateOutputPaths(sysTest.outputPathBase + "-solv")
    sysTest.updateOutputPaths(path)
    num=num+1
    sysTest.inputFiles.append(mgSetupXML)
    sysTest.inputFiles.append(alSetupXML)
    sysTest.inputFiles.append(vmSetupXML)
    sysTest.inputFiles.append(kiSetupXML)
    sysTest.solverOpts = mgSolverOpts

def suite():
    return mgSuite

if __name__ == "__main__":
    testRunner = SysTestRunner()
    testRunner.runSuite(solvSuite)

