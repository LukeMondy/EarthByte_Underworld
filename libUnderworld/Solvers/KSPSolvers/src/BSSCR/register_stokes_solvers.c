
#include <petsc.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>
#include <petscpc.h>

#include "common-driver-utils.h"
#ifdef HAVE_PETSCEXT
#include "pc_GtKG.h"
#include "pc_ScaledGtKG.h"
#endif

PetscErrorCode BSSCR_PetscExtStokesSolversInitialize( void )
{
#ifdef HAVE_PETSCEXT	
	PCRegister( "gtkg", "Solvers/KSPSolvers/src/BSSCR", "BSSCR_PCCreate_GtKG", BSSCR_PCCreate_GtKG );
//	PCRegister( "scgtkg", "pc/impls/gtkg", "BSSCR_PCCreate_ScGtKG", BSSCR_PCCreate_ScGtKG );
#endif
	
	PetscFunctionReturn(0);
}


PetscErrorCode BSSCR_PetscExtStokesSolversFinalize( void )
{
	
	PCRegisterDestroy();
	
	PetscFunctionReturn(0);
}



