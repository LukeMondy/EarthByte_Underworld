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

#include "pc-block-impl.h"
#include "pc-block-ops.h"
#include "private/pc/petscpc-block.h"


/* I don't think there is a requirement to keep the cases matching with what is in PCBlockType, just the order */
/*const char *PCBlockTypes[] = { "DIAGONAL", "UPPER", "LOWER", "FULL", "PCBlocType", "PC_BLOCK_", 0 };*/
const char *PCBlockTypes[] = { "diagonal", "upper", "lower", "full", "PCBlocType", "PC_BLOCK_", 0 };



PetscErrorCode Stg_PCDestroy_Block( PC pc )
{
	PetscErrorCode ierr;
	PC_Block      s = (PC_Block)pc->data;
	PetscInt i;
	
	
	PetscFunctionBegin;
	
	/* release implemenation internals */
	for( i=0; i<s->nr; i++ ) {
		Stg_KSPDestroy( & s->diag[i] );
	}
	free( s->diag );
	
	// Free up the mem for the mat array
	for( i=0; i<s->nr; i++ ) {
		free( s->mat[i] );
	}
	free( s->mat );
	
	
	
	/* free the vector */
	/*
	Stg_VecDestroy( & s->block_vec );
	*/
	if( s->t != PETSC_NULL ) {
		for( i=0; i<s->nr; i++ ) {
			if( s->t[i] != PETSC_NULL ) {
				Stg_VecDestroy( & s->t[i] );
				s->t[i] = PETSC_NULL;
			}
		}
		free( s->t );
	}
	
	if( s->t_trans != PETSC_NULL ) {
		for( i=0; i<s->nr; i++ ) {
			if( s->t_trans[i] != PETSC_NULL ) {
				Stg_VecDestroy( & s->t_trans[i] );
				s->t_trans[i] = PETSC_NULL;
			}
		}
		free( s->t_trans );
	}
	
	
	
	
	/* release implementation data pointer */
	ierr=PetscFree( s );CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


/* =====================================  Public functions  ====================================== */

EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "PCBlockGetSubKSP_Block"
PetscErrorCode PCBlockGetSubKSP_Block( PC pc, PetscInt i, KSP *sub_ksp )
{
	PC_Block       s;
	
	
	s = (PC_Block)pc->data;
	
	if( !pc->setupcalled ) {
		PCSetUp( pc );
	}
	
	
	if( i >= s->nr || i < 0 ) Stg_SETERRQ2( PETSC_ERR_ARG_SIZ, "Global dim %D. User value %D",s->nr, i );
	
	*sub_ksp = s->diag[i];
	
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "PCBlockGetSubKSP"
PetscErrorCode PCBlockGetSubKSP( PC pc, PetscInt i, KSP *sub_ksp )
{
	PetscErrorCode ierr,(*f)(PC,PetscInt,KSP*);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)pc,"PCBlockGetSubKSP_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(pc,i,sub_ksp);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}

/* =========================================================================== */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "PCBlockSetSubKSP_Block"
PetscErrorCode PCBlockSetSubKSP_Block( PC pc, PetscInt i, KSP sub_ksp )
{
	PC_Block       s;
	
	s = (PC_Block)pc->data;
	
	if( i >= s->nr || i < 0 ) Stg_SETERRQ2( PETSC_ERR_ARG_SIZ, "Global dim %D. User value %D",s->nr, i );
	
	
	
	if( s->diag[i] != PETSC_NULL ) {
		Stg_KSPDestroy( & s->diag[i] );
	}
	PetscObjectReference( (PetscObject)sub_ksp );
	
	s->diag[i] = sub_ksp;
	
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "PCBlockSetSubKSP"
PetscErrorCode PCBlockSetSubKSP( PC pc, PetscInt i, KSP sub_ksp )
{
	PetscErrorCode ierr,(*f)(PC,PetscInt,KSP);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)pc,"PCBlockSetSubKSP_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(pc,i,sub_ksp);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}

/* =========================================================================== */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "PCBlockGetSubMatrix_Block"
PetscErrorCode PCBlockGetSubMatrix_Block( PC pc, PetscInt i, PetscInt j, Mat *sub_mat )
{
	PC_Block       s;
	
	
	s = (PC_Block)pc->data;
	
	if( i >= s->nr || i < 0 ) Stg_SETERRQ2( PETSC_ERR_ARG_SIZ, "Global dim %D (row). User value %D",s->nr, i );
	if( j >= s->nc || j < 0 ) Stg_SETERRQ2( PETSC_ERR_ARG_SIZ, "Global dim %D (col). User value %D",s->nc, j );
	
	
	
	if( i == j ) {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, 
				"No matrices formally exist on the diagonal. \n"
				"If you require the operator associated with a PC at this (i,j) \n"
				"index, use PCBlockGetSubKSP followed by KSPGetOperators. \n" );
	}
	
	*sub_mat = s->mat[i][j];
	
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "PCBlockGetSubMatrix"
PetscErrorCode PCBlockGetSubMatrix( PC pc, PetscInt i, PetscInt j, Mat *sub_mat )
{
	PetscErrorCode ierr,(*f)(PC,PetscInt,PetscInt,Mat*);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)pc,"PCBlockGetSubMatrix_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(pc,i,j,sub_mat);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}

/* =========================================================================== */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "PCBlockSetSubMatrix_Block"
PetscErrorCode PCBlockSetSubMatrix_Block( PC pc, PetscInt i, PetscInt j, Mat sub_mat )
{
	PC_Block       s;
	
	s = (PC_Block)pc->data;
	
	if( i >= s->nr || i < 0 ) Stg_SETERRQ2( PETSC_ERR_ARG_SIZ, "Global dim %D (row). User value %D",s->nr, i );
	if( j >= s->nc || j < 0 ) Stg_SETERRQ2( PETSC_ERR_ARG_SIZ, "Global dim %D (col). User value %D",s->nc, j );
	
	
	s->mat[i][j] = sub_mat;
	
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "PCBlockSetSubMatrix"
PetscErrorCode PCBlockSetSubMatrix( PC pc, PetscInt i, PetscInt j, Mat sub_mat )
{
	PetscErrorCode ierr,(*f)(PC,PetscInt,PetscInt,Mat);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)pc,"PCBlockSetSubMatrix_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(pc,i,j,sub_mat);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}


/* =========================================================================== */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "PCBlockSetBlockType_Block"
PetscErrorCode PCBlockSetBlockType_Block( PC pc, PCBlockType bt )
{
	PC_Block       s;
	
	s = (PC_Block)pc->data;
	
	
	//_PCSetUp_BlockTypeValid( pc, &s->application_type );
	s->application_type = bt;
	
	/* Set the appropriate methods */
	if( s->application_type == PC_BLOCK_DIAGONAL ) {
		/* default is already set to be DIAGONAL so do nothing */
	}
	else if( s->application_type == PC_BLOCK_UPPER ) {
		pc->ops->apply = PCApply_Block_UPPER;
		pc->ops->applytranspose = PCApplyTranspose_Block_UPPER;
	}	
	else if( s->application_type == PC_BLOCK_LOWER ) {
		pc->ops->apply = PCApply_Block_LOWER;
		pc->ops->applytranspose = PCApplyTranspose_Block_LOWER;
	}
	else if( s->application_type == PC_BLOCK_FULL ) {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "PCBLOCK_FULL not implemented \n" );
	}
	
	PetscFunctionReturn(0);	
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "PCBlockSetBlockType"
PetscErrorCode PCBlockSetBlockType( PC pc, PCBlockType bt )
{
	PetscErrorCode ierr,(*f)(PC,PCBlockType);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)pc,"PCBlockSetBlockType_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(pc,bt);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}


/* =========================================================================== */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "PCBlockGetBlockType_Block"
PetscErrorCode PCBlockGetBlockType_Block( PC pc, const char **bt )
{
	PC_Block       s;
	
	s = (PC_Block)pc->data;
	
	*bt = PCBlockTypes[ s->application_type ];
	
	PetscFunctionReturn(0);	
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "PCBlockGetBlockType"
PetscErrorCode PCBlockGetBlockType( PC pc, const char **bt )
{
	PetscErrorCode ierr,(*f)(PC,const char**);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)pc,"PCBlockGetBlockType_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(pc,bt);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}


/* Constructor */
static struct _PCOps PCBlockOps_Values = {
/* 0 */		PCSetUp_Block,
			PCApply_Block_DIAGONAL,				/* These get specified depending on location of objects set */
			0, //PCApplyRichardson_Block,
			0, //PCApplyBA_Block,
			PCApplyTranspose_Block_DIAGONAL,
/* 5 */		0, //PCApplyBATranspose_Block,
			PCSetFromOptions_Block,
			0, //PCPreSolve_Block,
			0, //PCPostSolve_Block,
			0, //PCGetFactoredMatrix_Block,
/* 10 */	0, //PCApplySymmetricLeft_Block,
			0, //PCApplySymmetricRight_Block,
			0, //PCSetUpOnBlocks_Block,
			Stg_PCDestroy_Block,
			PCView_Block
};



/*
This will be called after PCSetOperators, i.e. during PCSetType
*/
EXTERN_C_BEGIN
PetscErrorCode PCCreate_Block( PC pc )
{
	PC_Block       s;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	/* define the operations */
	ierr = PetscMemcpy( pc->ops, &PCBlockOps_Values,sizeof(struct _PCOps) );CHKERRQ(ierr);
	
	
	/* allocate and set pointer for implememtation data */
	ierr = PetscMalloc( sizeof(struct _PC_Block), &s );CHKERRQ(ierr);
	ierr = PetscMemzero( s, sizeof(struct _PC_Block) );CHKERRQ(ierr);
	pc->data            = (void*)s;
	
	/* Init the implementation data */
	s->nr = 0;
	s->nc = 0;
	s->mat               = PETSC_NULL;
	s->diag              = PETSC_NULL;
	s->block_vec         = PETSC_NULL;
	s->t                 = PETSC_NULL;
	s->apply_setupcalled         = PETSC_FALSE;
	s->applytranspose_setupcalled = PETSC_FALSE;
	s->t_trans                    = PETSC_NULL;
	
	
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCBlockGetSubKSP_C",
			"PCBlockGetSubKSP_Block",
			PCBlockGetSubKSP_Block);CHKERRQ(ierr);
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCBlockSetSubKSP_C",
			"PCBlockSetSubKSP_Block",
			PCBlockSetSubKSP_Block);CHKERRQ(ierr);
	
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCBlockGetSubMatrix_C",
			"PCBlockGetSubMatrix_Block",
			PCBlockGetSubMatrix_Block);CHKERRQ(ierr);
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCBlockSetSubMatrix_C",
			"PCBlockSetSubMatrix_Block",
			PCBlockSetSubMatrix_Block);CHKERRQ(ierr);
	
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCBlockGetBlockType_C",
			"PCBlockGetBlockType_Block",
			PCBlockGetBlockType_Block);CHKERRQ(ierr);
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCBlockSetBlockType_C",
			"PCBlockSetBlockType_Block",
			PCBlockSetBlockType_Block);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}
EXTERN_C_END
