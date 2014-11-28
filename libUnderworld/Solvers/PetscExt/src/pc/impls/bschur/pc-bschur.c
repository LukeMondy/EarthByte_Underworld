/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2007 - 2008, 
 **        Dave A. May [dave.mayhem23@gmail.com]
 **        School of Mathematical Sciences, 
 **        Monash University
 **        Victoria, 3800
 **        Australia
 **
 **    Copyright (c) 2009 - 2010, 
 **        Dave A. May [dave.mayhem23@gmail.com]
 **        Geophysical Fluid Dynamics, 
 **        Department of Earth Sciences,
 **        ETH Zürich,
 **        Switzerland
 **
 **    Project:       PetscExt_v3.1
 **    Description:   A collection of extensions to the PETSc library.
 **
 **    All rights reserved.
 **    Redistribution and use in source and binary forms, with or without modification,
 **    are permitted provided that the following conditions are met:
 **
 **        * Redistributions of source code must retain the above copyright notice,
 **              this list of conditions and the following disclaimer.
 **        * Redistributions in binary form must reproduce the above copyright
 **              notice, this list of conditions and the following disclaimer in the
 **              documentation and/or other materials provided with the distribution.
 **        * Neither the name of the Monash University nor the names of its contributors
 **              may be used to endorse or promote products derived from this software
 **              without specific prior written permission.
 **
 **    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 **    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 **    THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 **    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
 **    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 **    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 **    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 **    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 **    LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
 **    OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 **
 **
 **    $Id$
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/


#include <stdlib.h>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscpc.h>
#include <petscversion.h>
#if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >=3) )
  #include <petsc-private/kspimpl.h>
  #include <petsc-private/pcimpl.h>
#else
  #include <private/kspimpl.h>
  #include <private/pcimpl.h>
#endif

#include "private/vec/petscvec-block.h"
#include "private/mat/petscmat-block.h"

#include "pc-bschur-impl.h"
#include "pc-bschur-ops.h"
#include "private/pc/pc-bschur.h"


const char *PCBSchurTypes[] = { "upper", "lower", "PCBSchurType", "PC_BSCHUR_", 0 };



PetscErrorCode Stg_PCDestroy_BSchur( PC pc )
{
	PetscErrorCode ierr;
	PC_BSchur      s = (PC_BSchur)pc->data;
	
	
	PetscFunctionBegin;
	
	if(s->schur != PETSC_NULL) {
		Stg_MatDestroy(&s->schur);
	}
	
	if(s->ksp_1 != PETSC_NULL) {
		Stg_KSPDestroy( &s->ksp_1);
	}
	if(s->ksp_2 != PETSC_NULL) {
		Stg_KSPDestroy( &s->ksp_2);
	}
	
	if( s->explicit_schur != PETSC_NULL ) {
		Stg_MatDestroy( &s->explicit_schur);
	}
	if( s->x_prime != PETSC_NULL ) {
		Stg_VecDestroy( &s->x_prime);
	}
	
	
	
	/* release implementation data pointer */
	ierr=PetscFree( s );CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


/* =====================================  Public functions  ====================================== */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "PCBSchurSetSubKSP_BSchur"
PetscErrorCode PCBSchurSetSubKSP_BSchur( PC pc, PetscInt i, KSP sub_ksp )
{
	PC_BSchur       s;
	
	s = (PC_BSchur)pc->data;
	
	if( i < 2 || i < 0) Stg_SETERRQ( PETSC_ERR_ARG_SIZ, "Only valid for 2 x 2 block matrix" );
	
	PetscObjectReference( (PetscObject)sub_ksp );
	if( i==0 ) {
		if( s->ksp_1 != PETSC_NULL ) {
			Stg_KSPDestroy( & s->ksp_1 );
		}
		s->ksp_1 = sub_ksp;
	}
	if( i==1 ) {
		if( s->ksp_2 != PETSC_NULL ) {
			Stg_KSPDestroy( & s->ksp_2 );
		}
		s->ksp_2 = sub_ksp;
	}
	
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "PCBSchurSetSubKSP"
PetscErrorCode PCBSchurSetSubKSP( PC pc, PetscInt i, KSP sub_ksp )
{
	PetscErrorCode ierr,(*f)(PC,PetscInt,KSP);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)pc,"PCBSchurSetSubKSP_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(pc,i,sub_ksp);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}

/* =========================================================================== */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "PCBSchurSetType_BSchur"
PetscErrorCode PCBSchurSetType_BSchur( PC pc, PCBSchurType bt )
{
	PC_BSchur       s;
	
	s = (PC_BSchur)pc->data;
	
	
	s->application_type = bt;
	
	/* Set the appropriate methods */
	if( s->application_type == PC_BSCHUR_UPPER ) {
		pc->ops->apply          = PCApply_BSchur_UPPER;
//		pc->ops->applytranspose = PCApplyTranspose_BSchur_UPPER;
	}	
	else if( s->application_type == PC_BSCHUR_LOWER ) {
		pc->ops->apply          = PCApply_BSchur_LOWER;
//		pc->ops->applytranspose = PCApplyTranspose_BSchur_LOWER;
	}
	else {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "PCBSchurType must be one of { PC_BSCHUR_UPPER, PC_BSCHUR_LOWER } \n" );
	}
	
	PetscFunctionReturn(0);	
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "PCBSchurSetType"
PetscErrorCode PCBSchurSetType( PC pc, PCBSchurType bt )
{
	PetscErrorCode ierr,(*f)(PC,PCBSchurType);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)pc,"PCBSchurSetType_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(pc,bt);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}


/* Constructor */
static struct _PCOps PCBSchurOps_Values = {
/* 0 */		PCSetUp_BSchur,
			PCApply_BSchur_UPPER,				/* These get specified depending on location of objects set */
			0, //PCApplyRichardson_Block,
			0, //PCApplyBA_Block,
			PCApply_BSchur_UPPER,
/* 5 */		0, //PCApplyBATranspose_Block,
			PCSetFromOptions_BSchur,
			0, //PCPreSolve_Block,
			0, //PCPostSolve_Block,
			0, //PCGetFactoredMatrix_Block,
/* 10 */	0, //PCApplySymmetricLeft_Block,
			0, //PCApplySymmetricRight_Block,
			0, //PCSetUpOnBlocks_Block,
			Stg_PCDestroy_BSchur,
			PCView_BSchur
};



/*
This will be called after PCSetOperators, i.e. during PCSetType
*/
EXTERN_C_BEGIN
PetscErrorCode PCCreate_BSchur( PC pc )
{
	PC_BSchur      s;
	PetscErrorCode ierr;
	MPI_Comm       comm;
	const char     *pc_prefix;
	
	
	
	PetscFunctionBegin;
	
	
	/* define the operations */
	ierr = PetscMemcpy( pc->ops, &PCBSchurOps_Values,sizeof(struct _PCOps) );CHKERRQ(ierr);
	
	
	/* allocate and set pointer for implememtation data */
	ierr = PetscMalloc( sizeof(struct _PC_BSchur), &s );CHKERRQ(ierr);
	ierr = PetscMemzero( s, sizeof(struct _PC_BSchur) );CHKERRQ(ierr);
	pc->data            = (void*)s;
	
	/* Init the implementation data */
	s->application_type  = PC_BSCHUR_UPPER;
	s->schur             = PETSC_NULL;
	s->explicit_operator = PETSC_FALSE;
	s->explicit_schur    = PETSC_NULL;
	s->ksp_1             = PETSC_NULL;
	s->ksp_2             = PETSC_NULL;
	s->block_factorisation = PETSC_FALSE;
	s->x_prime             = PETSC_NULL;
	
	
	PetscObjectGetComm( (PetscObject)pc, &comm );
	PCGetOptionsPrefix( pc, &pc_prefix );
	
	KSPCreate( comm, &s->ksp_1 );
	KSPSetOptionsPrefix( s->ksp_1, pc_prefix );
	KSPAppendOptionsPrefix( s->ksp_1, "pc_bschur_0_" );
	
	KSPCreate( comm, &s->ksp_2 );
	KSPSetOptionsPrefix( s->ksp_2, pc_prefix );
	KSPAppendOptionsPrefix( s->ksp_2, "pc_bschur_1_" );
	
	
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCBSchurSetType_C",
			"PCBSchurSetType_BSchur",
			PCBSchurSetType_BSchur);CHKERRQ(ierr);
	
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCBSchurSetSubKSP_C",
			"PCBSchurSetSubKSP_BSchur",
			PCBSchurSetSubKSP_BSchur);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}
EXTERN_C_END
