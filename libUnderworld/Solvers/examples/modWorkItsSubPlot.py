#!/usr/bin/env python
import sys
import os
from subprocess import call
import math
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab # In just python needed to install python-tk for this
                                 # But worked in ipython?

def setargsmodel(out, pp, vv, uw, mglevels, scr, scrtol, scrnormtype, a11, a11tol, resx, resy, model, pen):
    args=[]
    if "sinker" in model:
        args.append("--components.isoViscosity.eta0=1.0  -Xhelp ")
        args.append("--components.circleViscosity.eta0=%(vv)s " % {"vv": vv})
    if "cx" in model:
        args.append("--components.FieldTest.normaliseByAnalyticSolution=False ")
        args.append("--solCx_etaA=%(vv)s " % {"vv": vv})
        args.append("--solCx_etaB=1.0")
        args.append("--solCx_xc=0.5 ")
        args.append("--solCx_n=2.0 ")
        args.append("--wavenumberY=2.0 ")
    if "kx" in model:
        vv = math.log(vv) # natural log
        args.append("--components.FieldTest.normaliseByAnalyticSolution=False ")
        args.append("--solKx_n=2.0   --wavenumberX=2.0 ")
        args.append("--solKx_m=2.0   --wavenumberY=2.0 ")
        args.append("--solKx_sigma=1.0 ")
        args.append("--solKx_twiceB=%(vv)s " % {"vv": vv})

    args.append("--outputPath=%(out)s " % {"out": out})
    args.append("--particlesPerCell=%(pp)s " % {"pp": pp} )
    args.append("--saveDataEvery=1 --checkpointEvery=1 --checkpointWritePath=./%(out)s/Checkpoints --checkpointAppendStep=1 " % {"out": out})
    args.append("-force_correction 1 ")
    args.append("-scr_ksp_set_min_it_converge 1 ")
    args.append("-ksp_type bsscr -pc_type none -ksp_k2_type GMG  -augmented_lagrangian 1  -rescale_equations 0 -k_scale_only 1 --penaltyNumber=%(pen)s " % {"pen": pen} )
    args.append("-Q22_pc_type %(uw)s " % {"uw": uw})
    args.append("-remove_constant_pressure_null_space 1 ")
    args.append("--mgLevels=%(mglevels)s " % {"mglevels": mglevels})
    args.append("-scr_ksp_type %(scr)s " % {"scr": scr})
    args.append("-scr_ksp_rtol %(scrtol)s " % {"scrtol": scrtol})

    args.append("%(scrnormtype)s " % {"scrnormtype": scrnormtype})
    args.append("-A11_ksp_rtol %(a11tol)s " % {"a11tol": a11tol})
    args.append("-A11_ksp_type %(a11)s " % {"a11": a11})
    args.append("-A11_ksp_converged_reason")

    args.append("-change_backsolve 1")
    args.append("-backsolveA11_ksp_type preonly -backsolveA11_pc_type lu  ")
    args.append("-change_A11rhspresolve 1 -rhsA11_ksp_type preonly -rhsA11_pc_type lu -rhsA11_ksp_rtol 1.0e-6 ")
    #args.append("-backsolveA11_ksp_rtol 1.0e-3 -backsolveA11_pc_mg_galerkin true -restore_K 1 ")
    args.append("--elementResI=%(resx)s --elementResJ=%(resy)s " % {"resx": resx, "resy": resy})
    args.append("--maxTimeSteps=0 -Xdump_matvec -Xmatsuffix _%(resx)sx%(resy)s_10e%(vc)s_%(model)s_" % {"resx": resx, "resy": resy, "vc": jump, "model": model})
    return args

def setargsmodelsinker(out, pp, vv, uw, mglevels, scr, scrnormtype, a11, a11tol, resx, resy, model, pen):
    args=[]
    args.append("--outputPath=%(out)s " % {"out": out})
    args.append("-NN %(out)s " % {"out": out})
    args.append("--particlesPerCell=%(pp)s " % {"pp": pp} )
    args.append("--saveDataEvery=1 --checkpointEvery=1 --checkpointWritePath=./%(out)s/Checkpoints --checkpointAppendStep=1 " % {"out": out})
    args.append("--components.isoViscosity.eta0=1.0  -Xhelp ")
    args.append("--components.circleViscosity.eta0=%(vv)s " % {"vv": vv})
    args.append("-force_correction 1 ")
    args.append("-ksp_type bsscr -pc_type none -ksp_k2_type GMG  -augmented_lagrangian 1  -rescale_equations 0 -k_scale_only 1 --penaltyNumber=%(pen)s " % {"pen": pen} )
    args.append("-Q22_pc_type %(uw)s " % {"uw": uw})
    args.append("-remove_constant_pressure_null_space 1 ")
    args.append("--mgLevels=%(mglevels)s " % {"mglevels": mglevels})
    args.append("-scr_ksp_type %(scr)s " % {"scr": scr})
    args.append("%(scrnormtype)s " % {"scrnormtype": scrnormtype})
    args.append("-A11_ksp_rtol %(a11tol)s " % {"a11tol": a11tol})
    args.append("-A11_ksp_type %(a11)s " % {"a11": a11})
    args.append("-A11_ksp_converged_reason")
    #args.append("-backsolveA11_ksp_type fgmres -XbacksolveA11_ksp_monitor ")
    args.append("-backsolveA11_ksp_type preonly -backsolveA11_pc_type lu -backsolveA11_ksp_monitor ")
    args.append("-change_A11_prefix 1 -RHSA11_ksp_type preonly -RHSA11_pc_type lu -RHSA11_ksp_monitor ")
    #args.append("-backsolveA11_ksp_rtol 1.0e-3 -backsolveA11_pc_mg_galerkin true -restore_K 1 ")
    args.append("--elementResI=%(resx)s --elementResJ=%(resy)s " % {"resx": resx, "resy": resy})
    args.append("--maxTimeSteps=0 -dump_matvec -matsuffix _%(resx)sx%(resy)s_10e%(vc)s_%(model)s_" % {"resx": resx, "resy": resy, "vc": jump, "model": model})
    return args

#if debugging
#uwexec="cgdb --args "+uwexec
pylab.rcParams['figure.figsize'] = 12, 14  # create a 1200x800 image
pltIts={}
pltPtime={}
pltPress={}
f, ((pltIts[0], pltIts[1]), (pltPtime[0], pltPtime[1])) = plt.subplots(2, 2)
f, ((pltIts[0], pltIts[1]), (pltPtime[0], pltPtime[1]), (pltPress[0], pltPress[1])) = plt.subplots(3, 2)

#models=["sinker","sinkerq2p1"]
models=["cx","cxq2p1"]
#models=["kx","kxq2p1"]
modelNum=0

for model in models:
    count=1

    varmg=0

    element='Q1-P0'
    matrix='K'
    fac=2

    if 'q2q1' in model:
        element='Q2-Q1'
        fac=1
    if 'q2p1' in model:
        element='Q2-P1'
        fac=1
    if 'q2p1' in model:
        element='Q2-P1'
        fac=1


    mglevels=3

    pwd=os.getcwd()
    pwd=os.path.split(pwd)[0]
    pwd=os.path.split(pwd)[0]
    #uwexec=os.path.join(pwd,"build","bin","Underworld")
    uwexec=os.path.join("/home","mvelic","uworld","underworld2","libUnderworld","build","bin","Underworld")
    uwpath=os.path.join("/home","mvelic","uworld","underworld2","libUnderworld")
    
    pltlist=[]
    pltlegend=[]
    pltlist2=[]
    pltlegend2=[]
    
    resList=[16*fac,24*fac,32*fac]
    resCount=float(len(resList))
    resAlpha=resCount

    ncols=len(resList)
    #resList=[16*fac,24*fac]
    vjumpList=[2,4,6,8]
    #pens=[0.01, 0.5, 1.0, 3, 10, 20, 100, 500, 2000, 10000]
    pens=[0.01, 1.0, 20, 100, 500, 1000, 10000]
    pensCount=float(len(pens))

    lineColorsK=['m','y','r','g']
    lineMarkers=['p','s','o','d']
    
    PEN={}
    ITS={}
    PTIME={}
    PRESSL={}
    PRESSU={}
    
    fstr={}
    
    colstr=0
    
    penAll=[]
    
    scrp="default"
    mg="gmg"
    
    #output path
    
    uw="gkgdiag"
    
    scr="fgmres"
    a11=scr
    scrnormtype=""
    a11tol=1e-4
    scrtol=1e-7
    pp=40
    
    if model == "sinker":
        modelxml=uwpath+"/Solvers/InputFiles/sinker.xml "
    if model == "sinkerq2q1":
        modelxml=uwpath+"/Solvers/InputFiles/sinkerq2q1Nearest.xml "
    if model == "cx":
        modelxml=uwpath+"/Underworld/SysTest/PerformanceTests/testVelicSolCx.xml"
    if model == "cxq2q1":
        modelxml=uwpath+"/Solvers/InputFiles/testVelicSolCxQ2Q1Nearest.xml"
    if model == "cxq2p1":
        modelxml=uwpath+"/Solvers/InputFiles/testVelicSolCxQ2P1Nearest.xml"
    if model == "kx":
        modelxml=uwpath+"/Solvers/InputFiles/testVelicSolKx.xml"
    if model == "kxq2q1":
        modelxml=uwpath+"/Solvers/InputFiles/testVelicSolKxQ2Q1Nearest.xml"
    if model == "kxq2p1":
        modelxml=uwpath+"/Solvers/InputFiles/testVelicSolKxQ2P1Nearest.xml"
    
    mgop   = uwpath+"/Solvers/InputFiles/MultigridForRegularSCR.xml -options_file "+uwpath+"/Solvers/examples/options-scr-mg-rob.opt "
    auglag = uwpath+"/Solvers/InputFiles/AugLagStokesSLE-GtMG.xml "
    vmass  = uwpath+"/Solvers/InputFiles/VelocityMassMatrixSLE.xml "
    kspint = uwpath+"/Solvers/InputFiles/kspinterface.xml "
    vis    = uwpath+"/Solvers/InputFiles/analyticVis.xml "
    quiet  = uwpath+"/Solvers/InputFiles/quiet.xml "
    
    xml=[]
    xml.append(modelxml)
    xml.append(auglag)
    xml.append(vmass)
    xml.append(kspint)
#    if "sinker" not in model:
#        xml.append(vis)
    xml.append(quiet)
    xml.append(mgop)


    for res in resList:
        if varmg==1:
            mglevels=int(math.log(res)/math.log(2))
    
        scale="no_scale"
        colstr=0
        resx=resy=res


        if(varmg==1):
            dirstr="%(model)s_%(resx)sx%(resy)s_%(varmg)s" % {"model": model, "resx": resx, "resy": resy, "varmg":varmg}
        else:
            dirstr="%(model)s_%(resx)sx%(resy)s_X" % {"model": model, "resx": resx, "resy": resy}
        print dirstr
        if not os.path.exists(dirstr):
            os.mkdir(dirstr)
    
        for jump in vjumpList:
            vv = 10**jump
            #vv = math.log(va) # natural log
            PEN[jump]=[]
            ITS[jump]=[]
            PTIME[jump]=[]
            PRESSL[jump]=[]
            PRESSU[jump]=[]

            for penalty in pens:
                out=dirstr+"/"+"%(model)s_10e%(vc)s_%(pen)s_%(scrtol)s_%(a11tol)s" % {"model": model, "vc":jump, "pen":penalty, "scrtol": str(scrtol), "a11tol": str(a11tol)}
                PEN[jump].append(penalty)
                
                #if not os.path.exists(out):
                #    os.mkdir(out)
                args=setargsmodel(out, pp, vv, uw, mglevels, scr, scrtol, scrnormtype, a11, a11tol, resx, resy, model, penalty)
                cmdstr=uwexec+" "
                for xmlstr in xml:
                    xmlstr +=" \\\n"
                    cmdstr += xmlstr
                cmdargs=""
                for arg in args:
                    arg+=" \\\n"
                    cmdargs += arg
                #print cmdargs
                cmd=cmdstr+cmdargs+" > ./%s/output.txt 2>&1" % (out)
                #print cmd
                if not os.path.exists(out): # run the cmd only if out is not already created.
                    os.mkdir(out)
                    call(cmd,shell=True)
                f=open("./%s/output.txt" % (out),"r")
                ptime=0.0
                pits=0
                vsum=0
                plow=0.0
                pup=0.0
                varr=[]
                for line in f:
                    #print line
                    p = re.compile(r"""Pressure Solve:         = (\d.*\d) secs \/ (\d+) its""")
                    m = p.search(line)
                    if m != None:
                        #print m.group(1), m.group(2)
                        ptime  +=  float(m.group(1))
    
                    p = re.compile(r"""Linear solve converged due to CONVERGED_RTOL iterations (\d+)""")
                    m = p.search(line)
                    if m != None:
                        #print "v its= "+m.group(1)
                        vsum += int(m.group(1))
                        varr.append(int(m.group(1)))
                    p = re.compile(r"""min\/max\(p\)    = (.*) \[\d+\] \/ (.*) \[\d+\]""")
                    m = p.search(line)
                    if m != None:
                        plow = float(m.group(1))
                        pup  = float(m.group(2))
                print "vel its = "+str(vsum)+" sum= "+str(sum(varr))
                PTIME[jump].append(ptime)
                ITS[jump].append(vsum)
                PRESSL[jump].append(plow)
                PRESSU[jump].append(pup)
            #end penalty loop
            print PEN[jump]
            print ITS[jump]
            print res
            print jump
            print "resAlpha = "+str(resAlpha)
            #plt.subplot(211) 
            #pltarr[0].semilogx(PEN[jump], ITS[jump], marker=lineMarkers[colstr], linestyle='--', label="10e%(jump)s %(res)sx%(res)s" % {"jump":jump, "res":res})
            pltIts[modelNum].loglog(PEN[jump], ITS[jump], alpha=resAlpha/resCount, color=lineColorsK[colstr], marker=lineMarkers[colstr], linestyle='--', label="10e%(jump)s %(res)sx%(res)s" % {"jump":jump, "res":res})
            #plt.subplot(212) 
            #pltPtime[modelNum].semilogx(PEN[jump], PTIME[jump], marker=lineMarkers[colstr], linestyle='--', label="10e%(jump)s %(res)sx%(res)s" % {"jump":jump, "res":res})
            pltPtime[modelNum].loglog(PEN[jump], PTIME[jump], alpha=resAlpha/resCount, color=lineColorsK[colstr], marker=lineMarkers[colstr], linestyle='--', label="10e%(jump)s %(res)sx%(res)s" % {"jump":jump, "res":res})
            pltPress[modelNum].semilogx(PEN[jump], PRESSL[jump], alpha=resAlpha/resCount, color=lineColorsK[colstr], marker=lineMarkers[colstr], linestyle='--', label="10e%(jump)s %(res)sx%(res)s" % {"jump":jump, "res":res})
            pltPress[modelNum].semilogx(PEN[jump], PRESSU[jump], alpha=resAlpha/resCount, color=lineColorsK[colstr], marker=lineMarkers[colstr], linestyle='--', label="10e%(jump)s %(res)sx%(res)s" % {"jump":jump, "res":res})
            colstr=colstr+1
        #end jump loop
        resAlpha=resAlpha-1.0
    #end res loop
    plt.subplots_adjust(hspace=0.35)
    
    strT="Sol%(model)s Total Inner Work (its):\n %(element)s elements (scr/a11 %(scrtol)s/%(a11tol)s)" % {"model":model,"element":element,"scrtol":scrtol,"a11tol":a11tol}
    pltIts[modelNum].set_title(strT)
    pltIts[modelNum].set_xlabel("Penalty Number")
    pltIts[modelNum].set_ylabel("Total Iterations")
    
    #plt.legend( loc='upper center', ncol=3, shadow=True, fancybox=True, fontsize='small', framealpha=0.2)
    
    #For version of matplotlib < 1.2.1
    pltIts[modelNum].legend( loc='upper center', ncol=ncols, shadow=True, fancybox=True, prop={'size':6})
    pltIts[modelNum].grid(True)
    strIm="workItssol%(model)s_totalWork_%(scale)s" % {"model":model, "scale": scale}
    #pltIts[modelNum].savefig(strIm)
    
    ###################################
    ## Plot pressure solve time.. #####
    ###################################
    
    strT="Sol%(model)s Pressure Solve Time (s):\n %(element)s elements (scr/a11 %(scrtol)s/%(a11tol)s)" % {"model":model,"element":element,"scrtol":scrtol,"a11tol":a11tol}
    pltPtime[modelNum].set_title(strT)
    pltPtime[modelNum].set_xlabel("Penalty Number")
    pltPtime[modelNum].set_ylabel("Solve Time (s)")
    
    #plt.legend( loc='upper center', ncol=3, shadow=True, fancybox=True, fontsize='small', framealpha=0.2)
    
    #For version of matplotlib < 1.2.1
    pltPtime[modelNum].legend( loc='upper center', ncol=ncols, shadow=True, fancybox=True, prop={'size':6})
    pltPtime[modelNum].grid(True)
    ###################################
    ## Plot pressure bounds.. #####
    ###################################
    
    strT="Sol%(model)s Pressure Bounds:\n %(element)s elements (scr/a11 %(scrtol)s/%(a11tol)s)" % {"model":model,"element":element,"scrtol":scrtol,"a11tol":a11tol}
    pltPress[modelNum].set_title(strT)
    pltPress[modelNum].set_xlabel("Penalty Number")
    pltPress[modelNum].set_ylabel("Pressure")
    
    #plt.legend( loc='upper center', ncol=3, shadow=True, fancybox=True, fontsize='small', framealpha=0.2)
    
    #For version of matplotlib < 1.2.1
    #leg=pltPress[modelNum].legend( loc='upper center', ncol=ncols, shadow=True, fancybox=True, prop={'size':6})
    #leg.get_frame().set_alpha(0.0)

    pltPress[modelNum].grid(True)
    modelNum=modelNum+1

strIm="workPtimesol%(model)s_totalWork_%(scale)s_[%(scrtol)s__%(a11tol)s].pdf" % {"model":model, "scale": scale, "scrtol": str(scrtol), "a11tol": str(a11tol)}

plt.savefig(strIm, format='pdf')
    
# for models
