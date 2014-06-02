#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscpc.h>
#include <petscsnes.h>

#ifdef HAVE_PETSCEXT
#include <petscext.h>
#include <petscext_pc.h>
#endif

#include "private/kspimpl.h"   /*I "petscksp.h" I*/

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include "Solvers/KSPSolvers/KSPSolvers.h"

#ifdef HAVE_PETSCEXT
#include "Test/TestKSP.h"
#include "BSSCR/BSSCR.h"
#endif

/************************************************************************/
/*** This function is called from _StokesBlockKSPInterface_Initialise **************/
/************************************************************************/
#undef __FUNCT__  
#define __FUNCT__ "KSPRegisterAllKSP"
PetscErrorCode KSPRegisterAllKSP(const char path[])
{

    PetscFunctionBegin;
    #ifdef HAVE_PETSCEXT
    KSPRegisterTEST(path);/* not sure if the path matters much here. everything still worked even though I had it wrong */
    KSPRegisterBSSCR (path);
    #endif
    PetscFunctionReturn(0);
}

