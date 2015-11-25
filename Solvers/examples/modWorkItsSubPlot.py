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

def setargsmodel(out, pp, vv, uw, mglevels, scr, scrtol, scrnormtype, a11, a11tol, resx, resy, model, pen, gp):
    args=[]
    if "sinker" in model:
        args.append("--components.isoViscosity.eta0=1.0  -Xhelp ")
        args.append("--components.circleViscosity.eta0=%(vv)s " % {"vv": vv})
    if "cx" in model:
        args.append("--components.FieldTest.normaliseByAnalyticSolution=False ")
        args.append("--solCx_etaB=%(vv)s " % {"vv": vv})
        args.append("--solCx_etaA=1.0")
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

    args.append("-scr_ksp_max_it 2000")
    args.append("-A11_ksp_max_it 2000")

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

    #args.append("--gaussParticlesZ=2 ")
    args.append("--gaussParticlesX=%(gp)s " % {"gp": gp})
    args.append("--gaussParticlesY=%(gp)s " % {"gp": gp})


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
pylab.rcParams['figure.figsize'] = 12, 22  # create a 1200x800 image
pltIts={}
pltPtime={}
pltPress={}
pltVel={}
pltPits={}
#f, ((pltIts[0], pltIts[1]), (pltPtime[0], pltPtime[1])) = plt.subplots(2, 2)
#f, ((pltIts[0], pltIts[1]), (pltPtime[0], pltPtime[1]), (pltPress[0], pltPress[1])) = plt.subplots(3, 2)
#f, ((pltIts[0], pltIts[1]), (pltPtime[0], pltPtime[1]), (pltPress[0], pltPress[1]), (pltVel[0], pltVel[1])) = plt.subplots(4, 2)
f, ((pltIts[0], pltIts[1]), (pltPtime[0], pltPtime[1]), (pltPress[0], pltPress[1]), (pltVel[0], pltVel[1]), (pltPits[0], pltPits[1])) = plt.subplots(5, 2)

solve="LU"
solve="gmg"
basemodel="cx"
modelstr=basemodel
gp="6"
gpstr=gp+"x"+gp
#intstyle="pcdvc"
#intstyle="gauss"
intstyle="nearest" # Nearest
quadmodel="q2p1"

a11tol=1e-4
scrtol=1e-3

varmg=0
#######################################################################################################
#######################################################################################################
basemodelquad=basemodel+quadmodel+intstyle
models=[basemodel,basemodelquad]

##models=["sinker","sinkerq2p1"]
#models=["cx","cxq2p1"]
#models=["kx","kxq2p1"]
modelNum=0

for model in models:
    count=1

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
    
    resList=[16*fac,24*fac]
    #resList=[64*fac,128*fac]
    resCount=float(len(resList))
    resAlpha=resCount

    ncols=len(resList)
    #resList=[16*fac,24*fac]
    vjumpList=[2,6]
    #pens=[0.01, 0.5, 1.0, 3, 10, 20, 100, 500, 2000, 10000]
    #pens=[0.0, 0.01, 1.0, 20, 100, 500, 1000, 10000]
    pens=[0.01, 500, 10000]
    pensCount=float(len(pens))

    lineColorsK=['m','y','r','g']
    lineMarkers=['p','s','o','d']
    
    PEN={}
    ITS={}  # total velocity iterations
    PITS={} # pressure iterations
    PTIME={}
    PRESSL={}
    PRESSU={}
    VELOCITYL={}
    VELOCITYU={}

    fstr={}
    colstr=0
    penAll=[]
    scrp="default"
    mg="gmg"
    uw="gkgdiag"
    #uw="gtkg"
    
    scr="fgmres"
    a11=scr
    scrnormtype=""

    a11tolstr=str(a11tol)
    if "LU" in solve:
        a11tolstr="LU"

    pp=40
    
    if model == "sinker":
        modelxml=uwpath+"/Solvers/InputFiles/sinker.xml "
    if model == "sinkerq2q1nearest":
        modelxml=uwpath+"/Solvers/InputFiles/sinkerq2q1Nearest.xml "

    if model == "sinkerq2p1pcdvc":
        modelxml=uwpath+"/Solvers/InputFiles/sinkerq2p1pcdvc.xml "
    if model == "sinkerq2q1pcdvc":
        modelxml=uwpath+"/Solvers/InputFiles/sinkerq2q1PCDVC.xml "

    if model == "sinkerq2p1nearest":
        modelxml=uwpath+"/Solvers/InputFiles/sinkerq2p1Nearest.xml "
    if model == "sinkerq2p1gauss":
        modelxml=uwpath+"/Solvers/InputFiles/sinkerq2p1Gauss.xml "

    if model == "cx":
        modelxml=uwpath+"/Underworld/SysTest/PerformanceTests/testVelicSolCx.xml"
    if model == "cxq2q1nearest":
        modelxml=uwpath+"/Solvers/InputFiles/testVelicSolCxQ2Q1Nearest.xml"
    if model == "cxq2p1nearest":
        modelxml=uwpath+"/Solvers/InputFiles/testVelicSolCxQ2P1Nearest.xml"
    if model == "kx":
        modelxml=uwpath+"/Solvers/InputFiles/testVelicSolKx.xml"
    if model == "kxq2q1nearest":
        modelxml=uwpath+"/Solvers/InputFiles/testVelicSolKxQ2Q1Nearest.xml"
    if model == "kxq2p1nearest":
        modelxml=uwpath+"/Solvers/InputFiles/testVelicSolKxQ2P1Nearest.xml"
    

    if "LU" in solve:
        mgop   = " -options_file "+uwpath+"/Solvers/examples/options-scr-mumps-petsc3.opt "
    else:
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
    if "sinker" not in model:
        xml.append(vis)
    xml.append(quiet)
    xml.append(mgop)


    for res in resList:
        if varmg==1:
            mglevels=int(math.log(res)/math.log(2))
    
        scale="no_scale"
        colstr=0
        resx=resy=res

        if(varmg==1):
            #dirstr="%(model)s_%(resx)sx%(resy)s_%(varmg)s" % {"model": model, "resx": resx, "resy": resy, "varmg":varmg}
            dirstr="%(model)s_%(resx)sx%(resy)s_%(intstyle)s_%(gpstr)s_%(solve)s_%(varmg)s" % {"model": model, "resx": resx, "resy": resy, "intstyle": intstyle, "gpstr": gpstr, "solve": solve, "varmg":varmg}
        else:
            dirstr="%(model)s_%(resx)sx%(resy)s_%(intstyle)s_%(gpstr)s_%(solve)s_0.5_Switched" % {"model": model, "resx": resx, "resy": resy, "intstyle": intstyle, "gpstr": gpstr, "solve": solve}
        print dirstr
        if not os.path.exists(dirstr):
            os.mkdir(dirstr)
    
        for jump in vjumpList:
            vv = 10**jump
            #vv = math.log(va) # natural log
            PEN[jump]=[]
            ITS[jump]=[]
            PITS[jump]=[]
            PTIME[jump]=[]
            PRESSL[jump]=[]
            PRESSU[jump]=[]
            VELOCITYL[jump]=[]
            VELOCITYU[jump]=[]

            for penalty in pens:
                filestr="%(model)s_10e%(vc)s_%(pen)s_%(scrtol)s_%(a11tolstr)s" % {"model": model, "vc":jump, "pen":penalty, "scrtol": str(scrtol), "a11tolstr": str(a11tolstr)}
                out=dirstr+"/"+filestr
                PEN[jump].append(penalty)
                
                #if not os.path.exists(out):
                #    os.mkdir(out)
                args=setargsmodel(out, pp, vv, uw, mglevels, scr, scrtol, scrnormtype, a11, a11tol, resx, resy, model, penalty, gp)
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
                    outdir=uwpath+"/Solvers/examples/png"
                    if not os.path.exists(outdir):
                        os.mkdir(outdir)
                mvpngcmd="cp "+out+"/window.00000.png"+" "+uwpath+"/Solvers/examples/png/"+filestr+"_%(resx)sx%(resy)s.png" % {"resx": resx, "resy": resy}
                #call(mvpngcmd,shell=True)
                
                f=open("./%s/output.txt" % (out),"r")
                ptime=0.0
                pits=0
                vsum=0
                plow=0.0 #pressure bounds
                pup=0.0
                vlow=0.0 #velocity bounds
                vup=0.0
                varr=[]
                for line in f:
                    #print line
                    p = re.compile(r"""Pressure Solve:         = (\d.*\d) secs \/ (\d+) its""")
                    m = p.search(line)
                    if m != None:
                        #print m.group(1), m.group(2)
                        ptime  +=  float(m.group(1))
                        pits  +=  int(m.group(2))
                    p = re.compile(r"""Linear solve converged due to CONVERGED_RTOL iterations (\d+)""")
                    m = p.search(line)
                    if m != None:
                        #print "v its= "+m.group(1)
                        vsum += int(m.group(1))
                        varr.append(int(m.group(1)))
                    #Get pressure bounds
                    p = re.compile(r"""min\/max\(p\)    = (.*) \[\d+\] \/ (.*) \[\d+\]""")
                    m = p.search(line)
                    if m != None:
                        plow = float(m.group(1))
                        pup  = float(m.group(2))
                    #Get velocity bounds
                    p = re.compile(r"""min\/max\(u\)    = (.*) \[\d+\] \/ (.*) \[\d+\]""")
                    m = p.search(line)
                    if m != None:
                        vlow = float(m.group(1))
                        vup  = float(m.group(2))

                print "vel its = "+str(vsum)+" sum= "+str(sum(varr))
                PTIME[jump].append(ptime)
                ITS[jump].append(vsum)
                PITS[jump].append(pits)
                PRESSL[jump].append(plow)
                PRESSU[jump].append(pup)
                VELOCITYL[jump].append(vlow)
                VELOCITYU[jump].append(vup)
            #end penalty loop
            print PEN[jump]
            print ITS[jump]
            print res
            print jump
            print "resAlpha = "+str(resAlpha)
            #plt.subplot(211) 
            #pltarr[0].semilogx(PEN[jump], ITS[jump], marker=lineMarkers[colstr], linestyle='--', label="10e%(jump)s %(res)sx%(res)s" % {"jump":jump, "res":res})
            pltIts[modelNum].loglog(PEN[jump], ITS[jump], alpha=resAlpha/resCount, color=lineColorsK[colstr], marker=lineMarkers[colstr], linestyle='--', label="10e%(jump)s %(res)sx%(res)s" % {"jump":jump, "res":res})
            pltPits[modelNum].loglog(PEN[jump], PITS[jump], alpha=resAlpha/resCount, color=lineColorsK[colstr], marker=lineMarkers[colstr], linestyle='--', label="10e%(jump)s %(res)sx%(res)s" % {"jump":jump, "res":res})
            #plt.subplot(212) 
            #pltPtime[modelNum].semilogx(PEN[jump], PTIME[jump], marker=lineMarkers[colstr], linestyle='--', label="10e%(jump)s %(res)sx%(res)s" % {"jump":jump, "res":res})
            pltPtime[modelNum].loglog(PEN[jump], PTIME[jump], alpha=resAlpha/resCount, color=lineColorsK[colstr], marker=lineMarkers[colstr], linestyle='--', label="10e%(jump)s %(res)sx%(res)s" % {"jump":jump, "res":res})
            pltPress[modelNum].semilogx(PEN[jump], PRESSL[jump], alpha=resAlpha/resCount, color=lineColorsK[colstr], marker=lineMarkers[colstr], linestyle='--', label="10e%(jump)s %(res)sx%(res)s" % {"jump":jump, "res":res})
            pltPress[modelNum].semilogx(PEN[jump], PRESSU[jump], alpha=resAlpha/resCount, color=lineColorsK[colstr], marker=lineMarkers[colstr], linestyle='--', label="10e%(jump)s %(res)sx%(res)s" % {"jump":jump, "res":res})
            pltVel[modelNum].semilogx(PEN[jump], VELOCITYL[jump], alpha=resAlpha/resCount, color=lineColorsK[colstr], marker=lineMarkers[colstr], linestyle='--', label="10e%(jump)s %(res)sx%(res)s" % {"jump":jump, "res":res})
            pltVel[modelNum].semilogx(PEN[jump], VELOCITYU[jump], alpha=resAlpha/resCount, color=lineColorsK[colstr], marker=lineMarkers[colstr], linestyle='--', label="10e%(jump)s %(res)sx%(res)s" % {"jump":jump, "res":res})
            colstr=colstr+1
        #end jump loop
        resAlpha=resAlpha-1.0
    #end res loop
    plt.subplots_adjust(hspace=0.35)
    
    strT="Sol%(model)s Total Inner Work (its):\n %(element)s elements (scr/a11 %(scrtol)s/%(a11tol)s)" % {"model":model,"element":element,"scrtol":scrtol,"a11tol":a11tolstr}
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
    
    strT="Sol%(model)s Pressure Solve Time (s):\n %(element)s elements (scr/a11 %(scrtol)s/%(a11tol)s)" % {"model":model,"element":element,"scrtol":scrtol,"a11tol":a11tolstr}
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
    
    strT="Sol%(model)s Pressure Bounds:\n %(element)s elements (scr/a11 %(scrtol)s/%(a11tol)s)" % {"model":model,"element":element,"scrtol":scrtol,"a11tol":a11tolstr}
    pltPress[modelNum].set_title(strT)
    pltPress[modelNum].set_xlabel("Penalty Number")
    pltPress[modelNum].set_ylabel("Pressure")
    
    ###################################
    ## Plot velocity bounds.. #####
    ###################################
    
    strT="Sol%(model)s Velocity Bounds:\n %(element)s elements (scr/a11 %(scrtol)s/%(a11tol)s)" % {"model":model,"element":element,"scrtol":scrtol,"a11tol":a11tolstr}
    pltVel[modelNum].set_title(strT)
    pltVel[modelNum].set_xlabel("Penalty Number")
    pltVel[modelNum].set_ylabel("Velocity")

    ###################################
    ## Plot Pressure Iterations.. #####
    ###################################
    
    strT="Sol%(model)s Pressure Iterations:\n %(element)s elements (scr/a11 %(scrtol)s/%(a11tol)s)" % {"model":model,"element":element,"scrtol":scrtol,"a11tol":a11tolstr}
    pltPits[modelNum].set_title(strT)
    pltPits[modelNum].set_xlabel("Penalty Number")
    pltPits[modelNum].set_ylabel("Pressure Iterations")
    #plt.legend( loc='upper center', ncol=3, shadow=True, fancybox=True, fontsize='small', framealpha=0.2)
    
    #For version of matplotlib < 1.2.1
    #leg=pltPress[modelNum].legend( loc='upper center', ncol=ncols, shadow=True, fancybox=True, prop={'size':6})
    #leg.get_frame().set_alpha(0.0)

    pltPress[modelNum].grid(True)
    pltVel[modelNum].grid(True)
    pltPits[modelNum].grid(True)
    modelNum=modelNum+1

#strIm="work6x6-q1p0-q2p1_pcdvcsol%(model)s_totalWork_%(scale)s_[%(scrtol)s__%(a11tol)s].pdf" % {"model":model, "scale": scale, "scrtol": str(scrtol), "a11tol": str(a11tol)}
#strIm="sinker6x6-q1p0-q2p1gauss_6x6_totalWork_%(scale)s_[%(scrtol)s__%(a11tol)s].pdf" % {"model":model, "scale": scale, "scrtol": str(scrtol), "a11tol": str(a11tol)}
#strIm="hr-%(modelstr)s6x6-%(solve)s-q1p0-q2p1Nearest_6x6_totalWork_%(scale)s_[%(scrtol)s__%(a11tol)s].pdf" % {"modelstr":modelstr, "scale": scale, "scrtol": str(scrtol), "a11tol": str(a11tolstr), "solve":solve}
strIm="%(modelstr)s%(gpstr)s-%(solve)s-q1p0-%(quadmodel)s%(Intstyle)s_totalWork_%(scale)s_[%(scrtol)s__%(a11tol)s].pdf" % {"modelstr":modelstr, "gpstr": gpstr, "quadmodel":quadmodel, "Intstyle":intstyle.upper(), "scale": scale, "scrtol": str(scrtol), "a11tol": str(a11tolstr), "solve":solve}

if varmg==1:
    strIm="%(modelstr)s%(gpstr)s-%(solve)s-q1p0-%(quadmodel)s%(Intstyle)s_totalWork_%(scale)s_[%(scrtol)s__%(a11tol)s]_varMG.pdf" % {"modelstr":modelstr, "gpstr": gpstr, "quadmodel":quadmodel, "Intstyle":intstyle.upper(), "scale": scale, "scrtol": str(scrtol), "a11tol": str(a11tolstr), "solve":solve}

plt.savefig(strIm, format='pdf')
    
# for models
