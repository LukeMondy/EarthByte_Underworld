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

#include <petscversion.h>
#if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >=3) )
  #include "petsc-private/kspimpl.h"   /*I "petscksp.h" I*/
#else
  #include "private/kspimpl.h"   /*I "petscksp.h" I*/
#endif

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include "Solvers/KSPSolvers/KSPSolvers.h"
#include "BSSCR/petsccompat.h"

#include "Test/TestKSP.h"
#include "BSSCR/BSSCR.h"

/************************************************************************/
/*** This function is called from _StokesBlockKSPInterface_Initialise **************/
/************************************************************************/
#undef __FUNCT__  
#define __FUNCT__ "KSPRegisterAllKSP"
PetscErrorCode KSPRegisterAllKSP(const char path[])
{

    PetscFunctionBegin;
    KSPRegisterTEST(path);/* not sure if the path matters much here. everything still worked even though I had it wrong */
    KSPRegisterBSSCR (path);
    PetscFunctionReturn(0);
}

