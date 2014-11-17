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

#if debugging
#uwexec="cgdb --args "+uwexec
pylab.rcParams['figure.figsize'] = 12, 10  # create a 1200x800 image
pltK={}
pltS={}
f, ((pltK[0], pltK[1]), (pltS[0], pltS[1])) = plt.subplots(2, 2)

models=["sinker","sinkerq2q1"]
modelNum=0

minPen=0.01
maxPen=10000.0

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


    #pwd=os.getcwd()
    #pwd=os.path.split(pwd)[0]
    #pwd=os.path.split(pwd)[0]
    ##uwexec=os.path.join(pwd,"build","bin","Underworld")
    #uwexec=os.path.join("/home","mvelic","uworld","underworld2","libUnderworld","build","bin","Underworld")
    #uwpath=os.path.join("/home","mvelic","uworld","underworld2","libUnderworld")
    
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
    #pens=[0.01, 0.5, 1.0, 3, 10,20, 100, 500, 2000, 10000]
    #pensCount=float(len(pens))
    lineColorsK=['m','y','r','g']
    lineMarkers=['p','s','o','d']
    
    PEN={}
    K={}
    S={}
    #ITS={}
    #PTIME={}
    fstr={}

    colstr=0    
    #penAll=[]
    
    for res in resList:    
        scale="no_scale"
        colstr=0
        resx=resy=res    
        for jump in vjumpList:
            vv = 10**jump
            #vv = math.log(va) # natural log
            PEN[jump]=[]
            K[jump]=[]
            S[jump]=[]
            #ITS[jump]=[]
            #PTIME[jump]=[]
            fstr[jump]="plot_%(model)s_%(res)s_10e%(jump)s_%(scale)s.txt" % {"model": model, "res": res, "jump": jump, "scale": scale}
            ff = open(fstr[jump], "r")
            for line in ff:
                if float(line.split()[0]) >= minPen and  float(line.split()[0]) <= maxPen:
                    PEN[jump].append(line.split()[0])
                    K[jump].append(line.split()[2])
                    S[jump].append(line.split()[3])

            #penAll.extend(pen[jump])

            print PEN[jump]
            print K[jump]
            print res
            print jump
            print "resAlpha = "+str(resAlpha)
            #plt.subplot(211) 
            #pltarr[0].semilogx(PEN[jump], ITS[jump], marker=lineMarkers[colstr], linestyle='--', label="10e%(jump)s %(res)sx%(res)s" % {"jump":jump, "res":res})
            pltK[modelNum].loglog(PEN[jump], K[jump], alpha=resAlpha/resCount, color=lineColorsK[colstr], marker=lineMarkers[colstr], linestyle='--', label="10e%(jump)s %(res)sx%(res)s" % {"jump":jump, "res":res})
            #plt.subplot(212) 
            #pltS[modelNum].semilogx(PEN[jump], PTIME[jump], marker=lineMarkers[colstr], linestyle='--', label="10e%(jump)s %(res)sx%(res)s" % {"jump":jump, "res":res})
            pltS[modelNum].loglog(PEN[jump], S[jump], alpha=resAlpha/resCount, color=lineColorsK[colstr], marker=lineMarkers[colstr], linestyle='--', label="10e%(jump)s %(res)sx%(res)s" % {"jump":jump, "res":res})
            colstr=colstr+1
        #end jump loop
        resAlpha=resAlpha-1.0
    #end res loop
    plt.subplots_adjust(hspace=0.35)
    
    strT="Sol%(model)s K matrix Condition Numbers:\n %(element)s elements" % {"model":model,"element":element}
    pltK[modelNum].set_title(strT)
    pltK[modelNum].set_xlabel("Penalty Number")
    pltK[modelNum].set_ylabel("Condition Number ($\kappa_2$)")
    
    
    pltK[modelNum].grid(True)
    #For version of matplotlib < 1.2.1
    pltK[modelNum].legend( loc='upper center', ncol=ncols, shadow=True, fancybox=True, prop={'size':6})
    #else
    #pltK[modelNum].legend( loc='upper center', ncol=3, shadow=True, fancybox=True, fontsize='small', framealpha=0.2)

    strIm="condNumberssol%(model)s_%(scale)s" % {"model":model, "scale": scale}
    #pltK[modelNum].savefig(strIm)
    
    ###################################
    ## Plot pressure solve time.. #####
    ###################################
    
    strT="Sol%(model)s S matrix Condition Numbers:\n %(element)s elements" % {"model":model,"element":element}
    pltS[modelNum].set_title(strT)
    pltS[modelNum].set_xlabel("Penalty Number")
    pltS[modelNum].set_ylabel("Condition Number ($\kappa_2$)")
    

    pltS[modelNum].grid(True)
    #For version of matplotlib < 1.2.1
    pltS[modelNum].legend( loc='upper center', ncol=ncols, shadow=True, fancybox=True, prop={'size':6})
    #For version of matplotlib >= 1.2.1
    #pltS[modelNum]..legend( loc='upper center', ncol=3, shadow=True, fancybox=True, fontsize='small', framealpha=0.2)
    ##
    modelNum=modelNum+1

strIm="condNumberssol%(model)s__%(scale)s.pdf" % {"model":model, "scale": scale}

plt.savefig(strIm, format='pdf')
    
# for models
