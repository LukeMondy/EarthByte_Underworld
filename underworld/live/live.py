import underworld as uw

def getAssignFuncFromType(typ):

    fem=uw.StgFEM
    und=uw.Underworld
    pic=uw.PICellerator
    stg=uw.StGermain
    dom=uw.StgDomain
    luc=uw.gLucifer

    fl=[fem,und,pic,stg,dom,luc]

    func=0
    for lib in fl:
        #print lib
        for key in lib.__dict__.keys():
            #if typ in key and "Assign" in key:
            if typ+"_AssignFromXML" in key:
                print key, lib.__name__
                func=lib.__dict__[key]

    return func

def listComponentTypesDict():

    gd=uw.dictionary.GetDictionary()


    comps=gd["components"]

    for key in comps.keys():
        print comps[key]["Type"]

    return

def getListComponentTypesDict():

    tlist=[]

    gd=uw.dictionary.GetDictionary()


    comps=gd["components"]

    for key in comps.keys():
        #print comps[key]["Type"]
        tlist.append(comps[key]["Type"])

    return tlist

def getAssignFromXMLFuncList():
    flist=[]

    tlist=getListComponentTypesDict()

    for typ in tlist:
        func=getAssignFuncFromType(typ)
        flist.append(func)
        if func==0:
            print typ

    return flist
