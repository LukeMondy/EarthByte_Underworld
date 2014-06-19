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

/* contains the constructor, destructor and any interface functions for the class */


#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 3) )
  #include <petsc-private/matimpl.h>
#else
  #include <private/matimpl.h>
#endif

#include "private/vec/petscvec-block.h"
#include "private/mat/petscmat-block.h"

#include "mat-schur-impl.h"
#include "mat-schur-ops.h"
#include "private/mat/mat-schur.h"


const char *_MatSchurComplementType[] = { "MatSchur_A11", "MatSchur_A22", 0 };


/* external prototypes */
/*
PetscErrorCode MatSchurSetOperators_Schur( Mat A, Mat A11,Mat A12,Mat A21,Mat A22 );
PetscErrorCode MatSchurSetOperatorsFromBlock_Schur( Mat A, Mat BlockA );

PetscErrorCode MatSchurSetSchurComplementType_Schur( Mat A, MatSchurComplementType type );
PetscErrorCode MatSchurGetSchurComplementType_Schur( Mat A, MatSchurComplementType *type );

PetscErrorCode MatSchurApplyReductionToVec_Schur( Mat A, Vec f1, Vec f2, Vec subb );
PetscErrorCode MatSchurApplyReductionToVecFromBlock_Schur( Mat A, Vec F, Vec subb );

PetscErrorCode MatSchurSetScalar_Schur( Mat A, PetscScalar alpha );
PetscErrorCode MatSchurGetScalar_Schur( Mat A, PetscScalar *alpha );

PetscErrorCode MatSchurSetKSP_Schur( Mat A, KSP ksp );
PetscErrorCode MatSchurGetKSP_Schur( Mat A, KSP *ksp );
*/





#undef __FUNCT__
#define __FUNCT__ "MatDestroy_Schur"
PetscErrorCode MatDestroy_Schur( Mat A )
{
	PetscErrorCode ierr;
	Mat_Schur      s = (Mat_Schur)A->data;
	
	PetscFunctionBegin;
	/* release implemenation internals. Not NULL, then release */
	if (s->A11) {	Stg_MatDestroy( &s->A11 );	} /* If have A11, then free */
	if (s->A12) {	Stg_MatDestroy( &s->A12 );	}
	if (s->A21) {	Stg_MatDestroy( &s->A21 );	}
	if (s->A22) {	Stg_MatDestroy( &s->A22 );	}
	
	if (s->t1) {	Stg_VecDestroy( & s->t1 );	}
	if (s->t2) {	Stg_VecDestroy( & s->t2 );	}
	if (s->t1a) {	Stg_VecDestroy( & s->t1a );	}
	if (s->t2a) {	Stg_VecDestroy( & s->t2a );	}
	
	if (s->ksp) {	Stg_KSPDestroy( & s->ksp );	}
	
	/* release implementation data pointer */
	ierr=PetscFree(s);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


/*-----------------------------------------  Public interfaces  -------------------------------------------*/

EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "MatSchurSetOperators_Schur"
PetscErrorCode MatSchurSetOperators_Schur( Mat A, Mat A11,Mat A12,Mat A21,Mat A22 )
{
	Mat_Schur s = (Mat_Schur)A->data;
	
	PetscFunctionBegin;
	
	/* A11 */
	if (s->A11) {	Stg_MatDestroy( & s->A11 );						}
	s->A11 = A11;
	if (A11) {		PetscObjectReference( (PetscObject)A11 );	}
	
	/* A12 */
	if (s->A12) {	Stg_MatDestroy( & s->A12 );						}
	s->A12 = A12;
	if (A12) {		PetscObjectReference( (PetscObject)A12 );	}
	
	/* A21 */
	if (s->A21) {	Stg_MatDestroy( & s->A21 );						}
	s->A21 = A21;
	if (A21) {		PetscObjectReference( (PetscObject)A21 );	}
	
	/* A22 */
	if (s->A22) {	Stg_MatDestroy( & s->A22 );						}
	s->A22 = A22;
	if (A22 ) {		PetscObjectReference( (PetscObject)A22 );	}
	
	
	/* Check BOTH A12 and A21 have been specified */
	if (!s->A12) {
		Stg_SETERRQ( PETSC_ERR_ARG_NULL, "MatSchurSetOperators: You must specify A12" );
	}
	if (!s->A21) {
		Stg_SETERRQ( PETSC_ERR_ARG_NULL, "MatSchurSetOperators: You must specify A21" );
	}
	
	s->operators_set = PETSC_TRUE;
	
	
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "MatSchurSetOperators"
PetscErrorCode MatSchurSetOperators( Mat A, Mat A11,Mat A12,Mat A21,Mat A22)
{
	PetscErrorCode ierr,(*f)(Mat,Mat,Mat,Mat,Mat);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)A,"MatSchurSetOperators_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(A,A11,A12,A21,A22);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}



EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "MatSchurSetOperatorsFromBlock_Schur"
PetscErrorCode MatSchurSetOperatorsFromBlock_Schur( Mat A, Mat BlockA )
{
	Mat A11,A12,A21,A22;
	PetscTruth is_block;
	PetscInt M,N;
	
	PetscFunctionBegin;
	
	/* Check that BlockA is of type "block" */
	Stg_PetscTypeCompare( (PetscObject)BlockA, "block", &is_block );
	if (!is_block) {
		Stg_SETERRQ( PETSC_ERR_SUP, "MatSchurSetOperatorsFromBlock: Operator must be a block matrix" );
	}
	/* Check its a 2x2 block */
	MatGetSize( BlockA, &M,&N );
	if( M!=2 && N!=2 ) {
		Stg_SETERRQ2( PETSC_ERR_SUP, "MatSchurSetOperatorsFromBlock: Supplied block matrix has dimension %Dx%D. Must supply a 2x2 block matrix.", M,N );
	}
	
	
	/* Extract blocks */
	MatBlockGetSubMatrix( BlockA, 0,0, &A11 );
	MatBlockGetSubMatrix( BlockA, 0,1, &A12 );
	MatBlockGetSubMatrix( BlockA, 1,0, &A21 );
	MatBlockGetSubMatrix( BlockA, 1,1, &A22 );
	
	
	MatSchurSetOperators_Schur( A, A11,A12,A21,A22 );
	
	MatBlockRestoreSubMatrices(BlockA);
	
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "MatSchurSetOperatorsFromBlock"
PetscErrorCode MatSchurSetOperatorsFromBlock( Mat A, Mat BlockA )
{
	PetscErrorCode ierr,(*f)(Mat,Mat);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)A,"MatSchurSetOperatorsFromBlock_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(A,BlockA);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}

/* ------------- */

EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "MatSchurSetSchurComplementType_Schur"
PetscErrorCode MatSchurSetSchurComplementType_Schur( Mat A, MatSchurComplementType type )
{
	Mat_Schur s = (Mat_Schur)A->data;
	//PetscInt i;
	PetscTruth flg, valid_type;
	
	
	PetscFunctionBegin;
	
	/* Check operators set */
	if (!s->operators_set) {
		Stg_SETERRQ( PETSC_ERR_ORDER, "MatSchurSetSchurComplementType: You must call MatSchurSetOperators_Schur() first" );
	}
	
	/* Check if type is already set */
	if (s->stype) {
		Stg_SETERRQ1( PETSC_ERR_ARG_WRONGSTATE, "MatSchurSetSchurComplementType: Type already set to %s. This cannot be modified", s->stype );
	}
	
	//i = 0;
	valid_type = PETSC_FALSE;
	PetscStrcmp( type, "MatSchur_A11", &flg );
	if (flg) {
		s->stype = _MatSchurComplementType[0];
		
		/* Check we have provided the necesary operator */
		if (!s->A11) {
			Stg_SETERRQ( PETSC_ERR_ARG_NULL, "MatSchurSetSchurComplementType: You must specify A11 if build S_11" );
		}
		
		KSPSetOperators( s->ksp, s->A11, s->A11, SAME_NONZERO_PATTERN );
		MatGetVecs( s->A11, PETSC_NULL, &s->t1 );
		VecDuplicate( s->t1, &s->t1a );
		
		valid_type = PETSC_TRUE;
	}
	
	PetscStrcmp( type, "MatSchur_A22", &flg );
	if (flg) { // is TRUE
		s->stype = _MatSchurComplementType[1];
		
		/* Check we have provided the necesary operator */
		if (!s->A22) { // == PETSC_NULL 
			Stg_SETERRQ( PETSC_ERR_ARG_NULL, "MatSchurSetSchurComplementType: You must specify A22 if build S_22" );
		}
		
		KSPSetOperators( s->ksp, s->A22, s->A22, SAME_NONZERO_PATTERN );
		MatGetVecs( s->A22, PETSC_NULL, &s->t2 );
		VecDuplicate( s->t2, &s->t2a );
		
		/* set some special ops here */
		A->ops->mult          = MatMult_Schur_22;
		A->ops->multtranspose = MatMultTranspose_Schur_22;
		
		
		
		valid_type = PETSC_TRUE;
	}
	KSPSetFromOptions( s->ksp );
	
	if (!valid_type) { // is FALSE
		Stg_SETERRQ1( PETSC_ERR_ARG_UNKNOWN_TYPE, "MatSchurSetSchurComplementType: Type %s is invalid. Must be one of \"MatSchur_A11\",\"MatSchur_A22\"", type );
	}
	
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "MatSchurSetSchurComplementType"
PetscErrorCode MatSchurSetSchurComplementType( Mat A, MatSchurComplementType type )
{
	PetscErrorCode ierr,(*f)(Mat,MatSchurComplementType);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)A,"MatSchurSetSchurComplementType_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(A,type);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}

/* ------------- */

/*
"MatSchur_A11"
Solve A11 {t1} = {f1}
{b2} = A21 {t1}
{b2} = {b2} - {f2}

"MatSchur_A22"
Solve A22 {t2} = {f2}
{b1} = A12 {t2}
{b1} = {b1} - {f1}

*/
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "MatSchurApplyReductionToVec_Schur"
PetscErrorCode MatSchurApplyReductionToVec_Schur( Mat A, Vec f1,Vec f2, Vec subb )
{
	Mat_Schur s = (Mat_Schur)A->data;
	PetscTruth flg_11, flg_22;
	PetscInt subb_N, N;
	
	PetscFunctionBegin;
	
	
	/* check type is set */
	if (!s->stype) {
		Stg_SETERRQ( PETSC_ERR_ORDER, "You must call MatSchurSetSchurComplementType() first" );
	}
	
	VecGetSize( subb, &subb_N );
	
	PetscStrcmp( s->stype, "MatSchur_A11", &flg_11 );
	PetscStrcmp( s->stype, "MatSchur_A22", &flg_22 );
	if (flg_11) {
		/* check size of subb is compatable (matches f2) */
		VecGetSize( f2, &N );
		if( N != subb_N ) {
			Stg_SETERRQ( PETSC_ERR_SUP, "Output vector must have same dimensions as f2" );
		}

		/* A21 A11^{-1} f1 - f2 */		
		KSPSolve( s->ksp, f1, s->t1 );
		MatMult( s->A21, s->t1, subb );
		VecAXPY( subb, -1.0, f2 ); /* subb <- subb - f2 */
	}
	else if (flg_22) {
		/* check size of subb is compatable (matches f1) */
		VecGetSize( f1, &N );
		if( N != subb_N ) {
			Stg_SETERRQ( PETSC_ERR_SUP, "Output vector must have same dimensions as f1" );
		}
	
		/* A12 A22^{-1} f2 - f1 */	
		KSPSolve( s->ksp, f2, s->t2 );
		MatMult( s->A12, s->t2, subb );
		VecAXPY( subb, -1.0, f1 ); /* subb <- subb - f1 */
		
	}
	else {
		Stg_SETERRQ1( PETSC_ERR_ARG_WRONGSTATE, "MatSchurSetSchurComplementType %s is invalid. Must be one of \"MatSchur_A11\",\"MatSchur_A22\"", s->stype );
	}
	
	/* Scaling result if necessary */
	if (s->scale_set) {
		VecScale( subb, s->alpha );
	}
	
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "MatSchurApplyReductionToVec"
PetscErrorCode MatSchurApplyReductionToVec( Mat A, Vec f1, Vec f2, Vec subb )
{
	PetscErrorCode ierr,(*f)(Mat,Vec,Vec,Vec);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)A,"MatSchurApplyReductionToVec_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(A,f1,f2,subb);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}


EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "MatSchurApplyReductionToVecFromBlock_Schur"
PetscErrorCode MatSchurApplyReductionToVecFromBlock_Schur( Mat A, Vec F, Vec subb )
{
	PetscTruth is_block;
	Vec f1, f2;
	PetscInt M;
	
	PetscFunctionBegin;
	
	/* Check b is of type "block" */
	Stg_PetscTypeCompare( (PetscObject)F, "block", &is_block );
	if (!is_block) {
		Stg_SETERRQ( PETSC_ERR_SUP, "You can only apply reduction to a block vector. Try using MatSchurApplyReductionToVec()" );
	}
	/* Check its a 2x1 block */
	VecGetSize( F, &M );
	if( M!=2 ) {
		Stg_SETERRQ1( PETSC_ERR_SUP, "The block vector has dimension %Dx1. You must supply a 2x1 block vector", M );
	}
	
	VecBlockGetSubVector( F, 0, &f1 );
	VecBlockGetSubVector( F, 1, &f2 );
	
	MatSchurApplyReductionToVec_Schur( A, f1,f2, subb );
	
	VecBlockRestoreSubVectors( F );
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "MatSchurApplyReductionToVecFromBlock"
PetscErrorCode MatSchurApplyReductionToVecFromBlock( Mat A, Vec F, Vec subb )
{
	PetscErrorCode ierr,(*f)(Mat,Vec,Vec);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)A,"MatSchurApplyReductionToVecFromBlock_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) { /* If have f, (ie. not NULL) then call f() */
		ierr = (*f)(A,F,subb);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}

/* ------------- */

EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "MatSchurSetScalar_Schur"
PetscErrorCode MatSchurSetScalar_Schur( Mat A, PetscScalar alpha )
{
	Mat_Schur s = (Mat_Schur)A->data;
	
	PetscFunctionBegin;
	s->scale_set = PETSC_TRUE;
	s->alpha = alpha;
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "MatSchurSetScalar"
PetscErrorCode MatSchurSetScalar( Mat A, PetscScalar alpha )
{
	PetscErrorCode ierr,(*f)(Mat,PetscScalar);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)A,"MatSchurSetScalar_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(A,alpha);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}

/* ------------- */

EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "MatSchurGetScalar_Schur"
PetscErrorCode MatSchurGetScalar_Schur( Mat A, PetscScalar *alpha )
{
	Mat_Schur s = (Mat_Schur)A->data;
	
	PetscFunctionBegin;
	*alpha = s->alpha;
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "MatSchurGetScalar"
PetscErrorCode MatSchurGetScalar( Mat A, PetscScalar *alpha )
{
	PetscErrorCode ierr,(*f)(Mat,PetscScalar*);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)A,"MatSchurGetScalar_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(A,alpha);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}

/* ------------- */

EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "MatSchurGetKSP_Schur"
PetscErrorCode MatSchurGetKSP_Schur( Mat A, KSP *ksp )
{
	Mat_Schur s = (Mat_Schur)A->data;
	
	PetscFunctionBegin;
	*ksp = s->ksp;
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "MatSchurGetKSP"
PetscErrorCode MatSchurGetKSP( Mat A, KSP *ksp )
{
	PetscErrorCode ierr,(*f)(Mat,KSP*);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)A,"MatSchurGetKSP_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(A,ksp);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}

/* ---------- */

EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "MatSchurSetKSP_Schur"
PetscErrorCode MatSchurSetKSP_Schur( Mat A, KSP ksp )
{
	Mat_Schur s = (Mat_Schur)A->data;
	
	PetscFunctionBegin;
	
	if (s->ksp) {
		Stg_KSPDestroy( & s->ksp );
	}
	s->ksp = ksp;
	PetscObjectReference( (PetscObject)ksp );
	
	PetscFunctionReturn(0);
}
EXTERN_C_END


#undef __FUNCT__  
#define __FUNCT__ "MatSchurSetKSP"
PetscErrorCode MatSchurSetKSP( Mat A, KSP ksp )
{
	PetscErrorCode ierr,(*f)(Mat,KSP);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)A,"MatSchurSetKSP_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(A,ksp);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}

/* ------------- */

EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "MatSchurGetSchurComplementType_Schur"
PetscErrorCode MatSchurGetSchurComplementType_Schur( Mat A, MatSchurComplementType *type )
{
	Mat_Schur s = (Mat_Schur)A->data;
	
	PetscFunctionBegin;
	*type = s->stype;
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "MatSchurGetSchurComplementType"
PetscErrorCode MatSchurGetSchurComplementType( Mat A, MatSchurComplementType *type )
{
	PetscErrorCode ierr,(*f)(Mat,MatSchurComplementType*);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)A,"MatSchurGetSchurComplementType_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(A,type);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}


/* Constructors */
#undef __FUNCT__
#define __FUNCT__ "MatSetOps_Schur"
PetscErrorCode MatSetOps_Schur( struct _MatOps* ops )
{
	PetscFunctionBegin;

	/* 0*/
	ops->setvalues  = 0; // MatSetValues_Schur;
	ops->getrow     = 0; // MatGetRow_Schur;
	ops->restorerow = 0; // MatRestoreRow_Schur;
	ops->mult       = MatMult_Schur_11;
	ops->multadd    = MatMultAdd_Schur;
	/* 5*/
	ops->multtranspose    = MatMultTranspose_Schur_11;
	ops->multtransposeadd = MatMultTransposeAdd_Schur;
	ops->solve            = 0; // MatSolve_Schur;
	ops->solveadd         = 0; // MatSolveAdd_Schur;
	ops->solvetranspose   = 0; // MatSolveTranspose_Schur;
	/*10*/
	ops->solvetransposeadd = 0; // MatSolveTransposeAdd_Schur;
	ops->lufactor          = 0; // MatLUFactor_Schur;
	ops->choleskyfactor    = 0; // MatCholeskyFactor_Schur;
	ops->sor               = 0; // MatSOR_Schur;
	ops->transpose         = 0; // MatTranspose_Schur;
	/*15*/
	ops->getinfo       = 0; // MatGetInfo_Schur;
	ops->equal         = 0; // MatEqual_Schur;
	ops->getdiagonal   = 0; // MatGetDiagonal_Schur;
	ops->diagonalscale = 0; // MatDiagonalScale_Schur;
	ops->norm          = 0; // MatNorm_Schur;
	/*20*/
	ops->assemblybegin = MatAssemblyBegin_Schur;
	ops->assemblyend   = MatAssemblyEnd_Schur;
	ops->setoption     = 0; // MatSetOption_Schur;
	ops->zeroentries   = 0; // MatZeroEntries_Schur;
	/*24*/
	ops->zerorows               = 0; // MatZeroRows_Schur;
	ops->lufactorsymbolic       = 0; // MatLUFactorSymbolic_Schur;
	ops->lufactornumeric        = 0; // MatLUFactorNumeric_Schur;
	ops->choleskyfactorsymbolic = 0; // MatCholeskyFactorSymbolic_Schur;
	ops->choleskyfactornumeric  = 0; // MatCholeskyFactorNumeric_Schur;
	/*29*/

	ops->ilufactorsymbolic  = 0; // MatILUFactorSymbolic_Schur;
	ops->iccfactorsymbolic  = 0; // MatICCFactorSymbolic_Schur;


	/*34*/
	ops->duplicate     = 0; // MatDuplicate_Schur;
	ops->forwardsolve  = 0; // MatForwardSolve_Schur;
	ops->backwardsolve = 0; // MatBackwardSolve_Schur;
	ops->ilufactor     = 0; // MatILUFactor_Schur;
	ops->iccfactor     = 0; // MatICCFactor_Schur;
	/*39*/
	ops->axpy            = 0; // MatAXPY_Schur;
	ops->getsubmatrices  = 0; // MatGetSubMatrices_Schur;
	ops->increaseoverlap = 0; // MatIncreaseOverlap_Schur;
	ops->getvalues       = 0; // MatGetValues_Schur;
	ops->copy            = 0; // MatCopy_Schur;
	/*44*/
	ops->getrowmax   = 0; // MatGetRowMax_Schur;
	ops->scale       = MatScale_Schur;
	ops->shift       = 0; // MatShift_Schur;
	ops->diagonalset = 0; // MatDiagonalSet_Schur;
	//ops->dummy       = 0; // MatILUDTFactor_Schur;
	/*49*/

	ops->getrowij        = 0; // MatGetRowIJ_Schur;
	ops->restorerowij    = 0; // MatRestorRowIJ_Schur;
	ops->getcolumnij     = 0; // MatGetColumnIJ_Schur;
	ops->restorecolumnij = 0; // MatRestoreColumnIJ_Schur;
	/*54*/
	ops->fdcoloringcreate = 0; // MatFDColoringCreate_Schur;
	ops->coloringpatch    = 0; // MatColoringPatch_Schur;
	ops->setunfactored    = 0; // MatSetUnfactored_Schur;
	ops->permute          = 0; // MatPermute_Schur;
	ops->setvaluesblocked = 0; // MatSetValuesBlocked_Schur;
	/*59*/
	ops->getsubmatrix  = 0; // MatGetSubMatrix_Schur;
	ops->destroy       = MatDestroy_Schur;
	ops->view          = MatView_Schur;
	ops->convertfrom   = 0; // MatConvertFrom_Schur;

	/*64*/


	ops->setlocaltoglobalmapping = 0; // MatSetLocalToGlobalMapping_Schur;
	ops->setvalueslocal          = 0; // MatSetValuesLocal_Schur;
	ops->zerorowslocal           = 0; // MatZeroRowsLocal_Schur;
	/*69*/
	ops->getrowmaxabs    = 0; // MatGetRowMaxAbs_Schur;
	ops->getrowminabs    = 0; // 
	ops->convert         = 0; // MatConvert_Schur;
	ops->setcoloring     = 0; // MatSetColoring_Schur;

	/* 74 */
	ops->setvaluesadifor = 0; // MatSetValuesAdifor_Schur;
	ops->fdcoloringapply              = 0; // MatFDColoringApply_Schur;
	ops->setfromoptions               = MatSetFromOptions_Schur;
	ops->multconstrained              = 0; // MatMultConstrained_Schur;
	ops->multtransposeconstrained     = 0; // MatMultTransposeConstrained_Schur;
	/*79*/
	//ops->permutesparsify = 0; // MatPermuteSparsify_Schur;
	ops->mults           = 0; // MatMults_Schur;
	ops->solves          = 0; // MatSolves_Schur;
	ops->getinertia      = 0; // MatGetInertia_Schur;
	ops->load            = 0; // MatLoad_Schur;
	/*84*/
	ops->issymmetric             = 0; // MatIsSymmetric_Schur;
	ops->ishermitian             = 0; // MatIsHermitian_Schur;
	ops->isstructurallysymmetric = 0; // MatIsStructurallySymmetric_Schur;
	//ops->dummy4                  = 0; // MatPBRelax_Schur;
	ops->getvecs                 = MatGetVecs_Schur;
	/*89*/
	ops->matmult         = 0; // MatMatMult_Schur;
	ops->matmultsymbolic = 0; // MatMatMultSymbolic_Schur;
	ops->matmultnumeric  = 0; // MatMatMultNumeric_Schur;
	ops->ptap            = 0; // MatPtAP_Schur;
	ops->ptapsymbolic    = 0; // MatPtAPSymbolic_Schur;
	/*94*/
	ops->ptapnumeric              = 0; // MatPtAPNumeric_Schur;




	/*99*/



	ops->conjugate           = 0; // MatConjugate_Schur;

	/*104*/
	ops->setvaluesrow              = 0; // MatSetValuesRow_Schur;
	ops->realpart                  = 0; // MatRealPart_Schur;
	ops->imaginarypart             = 0; // MatImaginaryPart_Schur;
	ops->getrowuppertriangular     = 0; // MatGetRowUpperTriangular_Schur;
	ops->restorerowuppertriangular = 0; // MatRestoreRowUpperTriangular_Schur;
	/*109*/
	ops->matsolve           = 0; // MatMatSolve_Schur;
	ops->getredundantmatrix = 0; // MatGetRedundantMatrix_Schur;
	ops->getrowmin          = 0; // MatGetRowMin_Schur;
	ops->getcolumnvector    = 0;  // MatGetColumnVector_Schur;
	ops->missingdiagonal    = 0;
	/* 114 */
	ops->getseqnonzerostructure = 0;
	ops->create                 = 0;
	ops->getghosts              = 0;
	////ops->dummy2                 = 0;
	////ops->dummy3                 = 0;
	/* 119 */
	ops->multdiagonalblock = 0;
	ops->hermitiantranspose = 0;
	ops->multhermitiantranspose = 0;
	ops->multhermitiantransposeadd = 0;
	
	PetscFunctionReturn(0);
}

/*
The following operations are not to used with MatSchur
 - MatGetOwnershipRange(A,&s,&e):  Result, s=0,e=0.
 - MatGetOwnershripRanges(A,&r):   Result, error.
*/
EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "MatCreate_Schur"
PetscErrorCode MatCreate_Schur( Mat A )
{
	Mat_Schur       s;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* define the operations */
	MatSetOps_Schur( A->ops );
	
//	A->precision       = PETSC_SCALAR;
	//A->factor           = MAT_FACTOR_NONE;
	//A->mapping          = 0;
	A->spptr            = 0;
	A->same_nonzero     = PETSC_FALSE;
	A->assembled        = PETSC_FALSE;
	
	/* allocate and set pointer for implememtation data */
	ierr = PetscMalloc( sizeof(struct _Mat_Schur), &s );CHKERRQ(ierr);
	ierr = PetscMemzero( s, sizeof(struct _Mat_Schur) );CHKERRQ(ierr);
	A->data            = (void*)s;
	
	/* Init the implementation data */
	s->stype = PETSC_NULL;
	s->operators_set = PETSC_FALSE;
	s->scale_set = PETSC_FALSE;
	s->alpha = 1.0;
	s->A11 = PETSC_NULL;
	s->A12 = PETSC_NULL;
	s->A21 = PETSC_NULL;
	s->A22 = PETSC_NULL;
	s->t1 = s->t1a = PETSC_NULL;
	s->t2 = s->t2a = PETSC_NULL;
	
	KSPCreate( ((PetscObject)A)->comm, &s->ksp );
	KSPSetOptionsPrefix( s->ksp, "mat_schur_" );
	
	/* Set type */
	PetscObjectChangeTypeName((PetscObject) A, "schur" );
	
	/* define public functions */
	
	ierr = PetscObjectComposeFunctionDynamic( (PetscObject)A,"MatSchurSetOperators_C",
			"MatSchurSetOperators_Schur",
			MatSchurSetOperators_Schur );CHKERRQ(ierr);
	
	ierr = PetscObjectComposeFunctionDynamic( (PetscObject)A,"MatSchurSetOperatorsFromBlock_C",
			"MatSchurSetOperatorsFromBlock_Schur",
			MatSchurSetOperatorsFromBlock_Schur );CHKERRQ(ierr);
	
	/* ------ */
	
	ierr = PetscObjectComposeFunctionDynamic( (PetscObject)A,"MatSchurSetSchurComplementType_C",
			"MatSchurSetSchurComplementType_Schur",
			MatSchurSetSchurComplementType_Schur );CHKERRQ(ierr);
	
	ierr = PetscObjectComposeFunctionDynamic( (PetscObject)A,"MatSchurGetSchurComplementType_C",
			"MatSchurGetSchurComplementType_Schur",
			MatSchurGetSchurComplementType_Schur );CHKERRQ(ierr);
	
	/* ------ */
	
	ierr = PetscObjectComposeFunctionDynamic( (PetscObject)A,"MatSchurApplyReductionToVec_C",
			"MatSchurApplyReductionToVec_Schur",
			MatSchurApplyReductionToVec_Schur );CHKERRQ(ierr);
	
	ierr = PetscObjectComposeFunctionDynamic( (PetscObject)A,"MatSchurApplyReductionToVecFromBlock_C",
			"MatSchurApplyReductionToVecFromBlock_Schur",
			MatSchurApplyReductionToVecFromBlock_Schur );CHKERRQ(ierr);
	
	/* ------ */
	
	ierr = PetscObjectComposeFunctionDynamic( (PetscObject)A,"MatSchurSetScalar_C",
			"MatSchurSetScalar_Schur",
			MatSchurSetScalar_Schur );CHKERRQ(ierr);

	ierr = PetscObjectComposeFunctionDynamic( (PetscObject)A,"MatSchurGetScalar_C",
			"MatSchurGetScalar_Schur",
			MatSchurGetScalar_Schur );CHKERRQ(ierr);
	
	/* ------ */

	ierr = PetscObjectComposeFunctionDynamic( (PetscObject)A,"MatSchurGetKSP_C",
			"MatSchurGetKSP_Schur",
			MatSchurGetKSP_Schur );CHKERRQ(ierr);
	
	ierr = PetscObjectComposeFunctionDynamic( (PetscObject)A,"MatSchurSetKSP_C",
			"MatSchurSetKSP_Schur",
			MatSchurSetKSP_Schur );CHKERRQ(ierr);
	
	/* ------ */


	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "MatCreateSchur"
PetscErrorCode MatCreateSchur( MPI_Comm comm, Mat A11,Mat A12,Mat A21,Mat A22, PetscScalar alpha, MatSchurComplementType type, Mat *A )
{	
	PetscTruth flg;
	PetscInt m,n,M,N;
	
	
	PetscFunctionBegin;
	
	/* If don't have A12, then error */
	if (!A12) {	Stg_SETERRQ( PETSC_ERR_ARG_NULL, "MatCreateSchur: You must specify A12" );	}
	MatGetLocalSize(A12,&m,&n);
	MatGetSize(A12,&M,&N);
	
	MatCreate( comm, A );
	
	PetscStrcmp( type, "MatSchur_A11", &flg );
	if (flg) {		MatSetSizes( *A, n,n,N,N);		}
	/* If have flg, then types match */
	
	PetscStrcmp( type, "MatSchur_A22", &flg );
	if (flg) {		MatSetSizes( *A, m,m,M,M);		}
	
	MatSetType( *A, "schur" );
#if (((PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR>=3)) || (PETSC_VERSION_MAJOR>3) )
        MatSetUp( *A );
#endif
//	PetscObjectChangeTypeName((PetscObject)*A,"schur");
	
	
	MatSchurSetOperators_Schur( *A, A11,A12,A21,A22 );
	
	MatSchurSetSchurComplementType_Schur( *A, type );
	if (alpha) {
		MatSchurSetScalar_Schur( *A, alpha );
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatCreateSchurFromBlock"
PetscErrorCode MatCreateSchurFromBlock( Mat bmat, PetscScalar alpha, MatSchurComplementType type, Mat *A )
{	
	PetscTruth is_block;
	Mat A11,A12,A21,A22;
	PetscInt M,N;
	
	PetscFunctionBegin;
	
	/* Check b is of type "block" */
	Stg_PetscTypeCompare( (PetscObject)bmat, "block", &is_block );
	if (!is_block) {		Stg_SETERRQ( PETSC_ERR_SUP, "MatCreateSchurFromBlock: Must supply a block matrix" );		}
	
	/* Check its a 2x2 block */
	MatGetSize( bmat, &M,&N );
	if( M!=2 && N!=2 ) {
		Stg_SETERRQ2( PETSC_ERR_SUP, "MatCreateSchurFromBlock: Supplied block matrix has dimension %Dx%D. Must supply a 2x2 block matrix.", M,N );
	}
	
	MatBlockGetSubMatrix( bmat, 0,0, &A11 );
	MatBlockGetSubMatrix( bmat, 0,1, &A12 );
	MatBlockGetSubMatrix( bmat, 1,0, &A21 );
	MatBlockGetSubMatrix( bmat, 1,1, &A22 );
	
	MatCreateSchur( ((PetscObject)bmat)->comm, A11,A12,A21,A22, alpha,type, A );
	
	MatBlockRestoreSubMatrices(bmat);
	
	PetscFunctionReturn(0);
}
