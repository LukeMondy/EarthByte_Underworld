#include <petsc.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>
#include <petscpc.h>

#include "common-driver-utils.h"
#include "pc_GtKG.h"
#include "pc_ScaledGtKG.h"


PetscErrorCode BSSCR_PetscExtStokesSolversInitialize( void )
{
	Stg_PCRegister( "gtkg", "Solvers/KSPSolvers/src/BSSCR", "BSSCR_PCCreate_GtKG", BSSCR_PCCreate_GtKG );
	PetscFunctionReturn(0);
}


PetscErrorCode BSSCR_PetscExtStokesSolversFinalize( void )
{
	
#if ((PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR>=4))
  PCFinalizePackage();
#else
  PCRegisterDestroy();
#endif
	
	PetscFunctionReturn(0);
}



