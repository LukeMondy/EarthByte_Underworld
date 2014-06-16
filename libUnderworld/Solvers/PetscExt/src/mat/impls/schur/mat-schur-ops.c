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
#include "private/mat/mat-schur.h"
#include "mat-schur-ops.h"


extern const char *_MatSchurComplementType[];

/* for logging */
PetscLogStage __PETScExt_MatSchur_MatMult_Stage;
PetscLogStage __PETScExt_MatSchur_MatMultTranspose_Stage;

PetscLogEvent __PETScExt_MatSchur_MatMults_Event;
PetscLogEvent __PETScExt_MatSchur_A11_inverse_Event;
PetscLogEvent __PETScExt_MatSchur_A22_inverse_Event;


/* 0*/
#undef __FUNCT__
#define __FUNCT__ "MatSetValues_Schur"
PetscErrorCode MatSetValues_Schur(Mat mat,PetscInt m,const PetscInt idxm[],PetscInt n,const PetscInt idxn[],const PetscScalar v[],InsertMode addv)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatGetRow_Schur"
PetscErrorCode MatGetRow_Schur(Mat mat,PetscInt row,PetscInt *ncols,PetscInt *cols[],PetscScalar *vals[])
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatRestoreRow_Schur"
PetscErrorCode MatRestoreRow_Schur(Mat mat,PetscInt row,PetscInt *ncols,PetscInt *cols[],PetscScalar *vals[])
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );

	PetscFunctionReturn(0);
}


//MatMult_Schur
/*
"MatSchur_A11": y = S_11 x
Solve A11 {t1} = {x}
{y} = A21 {t1}
{y} <- -{y}
{y} = {y} + A22 {x}
{y} <- -{y}
*/
#undef __FUNCT__
#define __FUNCT__ "MatMult_Schur_11"
PetscErrorCode MatMult_Schur_11(Mat mat,Vec x,Vec y)
{
	Mat_Schur s = (Mat_Schur)mat->data;
	PetscScalar final_scale;
	PetscTruth apply;
	
	PetscFunctionBegin; 
	
	//PetscLogStagePush( __PETScExt_MatSchur_MatMult_Stage );
	
	//PetscLogEventBegin(__PETScExt_MatSchur_MatMults_Event,0,0,0,0);
	/* this will save us doing one extra vector scale */
	apply = PETSC_FALSE;
	if (s->A22!=PETSC_NULL && s->scale_set == PETSC_TRUE ) {
		final_scale = - s->alpha;
		apply = PETSC_TRUE;
	}
	else if (s->A22!=PETSC_NULL && s->scale_set == PETSC_FALSE ) {
		final_scale = -1.0;
		apply = PETSC_TRUE;
	}
	else if (s->A22==PETSC_NULL && s->scale_set==PETSC_TRUE ) {
		final_scale = s->alpha;
		apply = PETSC_TRUE;
	}
	
	
	MatMult( s->A12, x, s->t1a );
	//PetscLogEventEnd(__PETScExt_MatSchur_MatMults_Event,0,0,0,0);
	
	//PetscLogEventBegin(__PETScExt_MatSchur_A11_inverse_Event,0,0,0,0);
	KSPSolve( s->ksp, s->t1a, s->t1 );
	//PetscLogEventEnd(__PETScExt_MatSchur_A11_inverse_Event,0,0,0,0);
	
	//PetscLogEventBegin(__PETScExt_MatSchur_MatMults_Event,0,0,0,0);
	MatMult( s->A21, s->t1, y );
	
	if( s->A22 != PETSC_NULL ) {
		VecScale( y, -1.0 );
		
		MatMultAdd( s->A22,x, y,y );
	}
	
	if (apply) {
		VecScale( y, final_scale );
	}
	//PetscLogEventEnd(__PETScExt_MatSchur_MatMults_Event,0,0,0,0);
	
	//PetscLogStagePop();
	
	PetscFunctionReturn(0);
}

/*
"MatSchur_A22": y = S_22 x
Solve A22 {t2} = {x}
{y} = A12 {t2}

{y} <- -{y}
{y} = {y} + A11 {x}
{y} <- -{y}
*/
#undef __FUNCT__
#define __FUNCT__ "MatMult_Schur_22"
PetscErrorCode MatMult_Schur_22(Mat mat,Vec x,Vec y)
{
	Mat_Schur s = (Mat_Schur)mat->data;
	PetscScalar final_scale;
	PetscTruth apply;
	
	PetscFunctionBegin; 
	
	//PetscLogStagePush( __PETScExt_MatSchur_MatMult_Stage );
	
	//PetscLogEventBegin(__PETScExt_MatSchur_MatMults_Event,0,0,0,0);
	/* this will save us doing one extra vector scale */
	apply = PETSC_FALSE;
	if (s->A11!=PETSC_NULL && s->scale_set == PETSC_TRUE ) {
		final_scale = - s->alpha;
		apply = PETSC_TRUE;
	}
	else if (s->A11!=PETSC_NULL && s->scale_set == PETSC_FALSE ) {
		final_scale = -1.0;
		apply = PETSC_TRUE;
	}
	else if (s->A11==PETSC_NULL && s->scale_set==PETSC_TRUE ) {
		final_scale = s->alpha;
		apply = PETSC_TRUE;
	}
	
	
	MatMult( s->A21, x, s->t2a );
	//PetscLogEventEnd(__PETScExt_MatSchur_MatMults_Event,0,0,0,0);
	
	//PetscLogEventBegin(__PETScExt_MatSchur_A22_inverse_Event,0,0,0,0);
	KSPSolve( s->ksp, s->t2a, s->t2 );
	//PetscLogEventEnd(__PETScExt_MatSchur_A22_inverse_Event,0,0,0,0);
	
	//PetscLogEventBegin(__PETScExt_MatSchur_MatMults_Event,0,0,0,0);
	MatMult( s->A12, s->t2, y );
	
	if( s->A11 != PETSC_NULL ) {
		VecScale( y, -1.0 );
		
		MatMultAdd( s->A11,x, y,y );
	}
	
	
	if (apply) {
		VecScale( y, final_scale );
	}
	
	//PetscLogEventEnd(__PETScExt_MatSchur_MatMults_Event,0,0,0,0);
	
	//PetscLogStagePop();
	
	PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "MatMultAdd_Schur"
PetscErrorCode MatMultAdd_Schur(Mat mat,Vec v1,Vec v2,Vec v3)
{
	PetscFunctionBegin;
	
	MatMult( mat,v1, v2 );
	VecAXPY( v3, 1.0, v2 );
	
	PetscFunctionReturn(0);
}


/* 5*/
#undef __FUNCT__
#define __FUNCT__ "MatMultTranspose_Schur_11"
PetscErrorCode MatMultTranspose_Schur_11(Mat mat,Vec x,Vec y)
{
	Mat_Schur s = (Mat_Schur)mat->data;
	PetscScalar final_scale;
	PetscTruth apply;
	
	PetscFunctionBegin; 
	
	//PetscLogStagePush( __PETScExt_MatSchur_MatMultTranspose_Stage );
	
	//PetscLogEventBegin(__PETScExt_MatSchur_MatMults_Event,0,0,0,0);
	
	/* this will save us doing one extra vector scale */
	apply = PETSC_FALSE;
	if (s->A22!=PETSC_NULL && s->scale_set == PETSC_TRUE ) {
		final_scale = - s->alpha;
		apply = PETSC_TRUE;
	}
	else if (s->A22!=PETSC_NULL && s->scale_set == PETSC_FALSE ) {
		final_scale = -1.0;
		apply = PETSC_TRUE;
	}
	else if (s->A22==PETSC_NULL && s->scale_set==PETSC_TRUE ) {
		final_scale = s->alpha;
		apply = PETSC_TRUE;
	}
	
	
	
	MatMultTranspose( s->A21, x, s->t1a );
	//PetscLogEventEnd(__PETScExt_MatSchur_MatMults_Event,0,0,0,0);
	
	//PetscLogEventBegin(__PETScExt_MatSchur_A11_inverse_Event,0,0,0,0);
	KSPSolveTranspose( s->ksp, s->t1a, s->t1 );
	//PetscLogEventEnd(__PETScExt_MatSchur_A11_inverse_Event,0,0,0,0);
	
	//PetscLogEventBegin(__PETScExt_MatSchur_MatMults_Event,0,0,0,0);
	MatMultTranspose( s->A12, s->t1, y );
	
	if( s->A22 != PETSC_NULL ) {
		VecScale( y, -1.0 );
		
		MatMultTransposeAdd( s->A22,x, y,y );
		
	}
	
	
	if (apply) {
		VecScale( y, final_scale );
	}
	//PetscLogEventEnd(__PETScExt_MatSchur_MatMults_Event,0,0,0,0);
	
	//PetscLogStagePop();
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatMultTranspose_Schur_22"
PetscErrorCode MatMultTranspose_Schur_22(Mat mat,Vec x,Vec y)
{
	Mat_Schur s = (Mat_Schur)mat->data;
	PetscScalar final_scale;
	PetscTruth apply;
	
	PetscFunctionBegin; 
	
	//PetscLogStagePush( __PETScExt_MatSchur_MatMultTranspose_Stage );
	
	//PetscLogEventBegin(__PETScExt_MatSchur_MatMults_Event,0,0,0,0);
	
	/* this will save us doing one extra vector scale */
	apply = PETSC_FALSE;
	if (s->A11!=PETSC_NULL && s->scale_set == PETSC_TRUE ) {
		final_scale = - s->alpha;
		apply = PETSC_TRUE;
	}
	else if (s->A11!=PETSC_NULL && s->scale_set == PETSC_FALSE ) {
		final_scale = -1.0;
		apply = PETSC_TRUE;
	}
	else if (s->A11==PETSC_NULL && s->scale_set==PETSC_TRUE ) {
		final_scale = s->alpha;
		apply = PETSC_TRUE;
	}
	
	MatMultTranspose( s->A12, x, s->t2a );
	//PetscLogEventEnd(__PETScExt_MatSchur_MatMults_Event,0,0,0,0);
	
	//PetscLogEventBegin(__PETScExt_MatSchur_A22_inverse_Event,0,0,0,0);
	KSPSolveTranspose( s->ksp, s->t2a, s->t2 );
	//PetscLogEventEnd(__PETScExt_MatSchur_A22_inverse_Event,0,0,0,0);
	
	//PetscLogEventBegin(__PETScExt_MatSchur_MatMults_Event,0,0,0,0);
	MatMultTranspose( s->A21, s->t2, y );
	
	if( s->A11 != PETSC_NULL ) {
		VecScale( y, -1.0 );
		
		MatMultTransposeAdd( s->A11,x, y,y );
		
	}
	
	if (apply) {
		VecScale( y, final_scale );
	}
	//PetscLogEventEnd(__PETScExt_MatSchur_MatMults_Event,0,0,0,0);
	
	//PetscLogStagePop();
	
	PetscFunctionReturn(0);
}





#undef __FUNCT__
#define __FUNCT__ "MatMultTransposeAdd_Schur"
PetscErrorCode MatMultTransposeAdd_Schur(Mat mat,Vec v1,Vec v2,Vec v3)
{
	PetscFunctionBegin; 
	
	MatMultTranspose( mat, v1, v2 );
	VecAXPY( v3, 1.0, v2 );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSolve_Schur"
PetscErrorCode MatSolve_Schur(Mat mat,Vec b,Vec x)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSolveAdd_Schur"
PetscErrorCode MatSolveAdd_Schur(Mat mat,Vec b,Vec y,Vec x)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatSolveTranspose_Schur"
PetscErrorCode MatSolveTranspose_Schur(Mat mat,Vec b,Vec x)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*10*/
#undef __FUNCT__
#define __FUNCT__ "MatSolveTransposeAdd_Schur"
PetscErrorCode MatSolveTransposeAdd_Schur(Mat mat,Vec b,Vec y,Vec x)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatLUFactor_Schur"
PetscErrorCode MatLUFactor_Schur(Mat mat,IS row,IS col,MatFactorInfo *info)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatCholeskyFactor_Schur"
PetscErrorCode MatCholeskyFactor_Schur(Mat mat,IS perm,MatFactorInfo *info)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatRelax_Schur"
PetscErrorCode MatRelax_Schur(Mat mat,Vec b,PetscReal omega,MatSORType flag,PetscReal shift,PetscInt its,PetscInt lits,Vec x)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatTranspose_Schur"
PetscErrorCode MatTranspose_Schur(Mat mat,MatReuse reuse,Mat *B)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*15*/

#undef __FUNCT__
#define __FUNCT__ "MatGetInfo_Schur"
PetscErrorCode MatGetInfo_Schur(Mat mat,MatInfoType flag,MatInfo *info)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatEqual_Schur"
PetscErrorCode MatEqual_Schur(Mat A,Mat B,PetscTruth *flg)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatGetDiagonal_Schur"
PetscErrorCode MatGetDiagonal_Schur(Mat mat,Vec v)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatDiagonalScale_Schur"
PetscErrorCode MatDiagonalScale_Schur(Mat mat,Vec l,Vec r)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatNorm_Schur"
PetscErrorCode MatNorm_Schur(Mat mat,NormType type,PetscReal *nrm)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*20*/
#undef __FUNCT__
#define __FUNCT__ "MatAssemblyBegin_Schur"
PetscErrorCode MatAssemblyBegin_Schur(Mat mat,MatAssemblyType type)
{
	Mat_Schur s = (Mat_Schur)mat->data;
	PetscTruth assembled;
	
	PetscFunctionBegin; 
	
	
	if( s->A11 != PETSC_NULL ) {
		MatAssembled( s->A11, &assembled );
		if( assembled == PETSC_FALSE ) {	MatAssemblyBegin( s->A11, MAT_FINAL_ASSEMBLY ); 	}
	}
	if( s->A12 != PETSC_NULL ) {
		MatAssembled( s->A12, &assembled );
		if( assembled == PETSC_FALSE ) {	MatAssemblyBegin( s->A12, MAT_FINAL_ASSEMBLY ); 	}
	}
	if( s->A21 != PETSC_NULL ) {
		MatAssembled( s->A21, &assembled );
		if( assembled == PETSC_FALSE ) {	MatAssemblyBegin( s->A21, MAT_FINAL_ASSEMBLY ); 	}
	}
	if( s->A22 != PETSC_NULL ) {
		MatAssembled( s->A22, &assembled );
		if( assembled == PETSC_FALSE ) {	MatAssemblyBegin( s->A22, MAT_FINAL_ASSEMBLY ); 	}
	}
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatAssemblyEnd_Schur"
PetscErrorCode MatAssemblyEnd_Schur(Mat mat,MatAssemblyType type)
{
	Mat_Schur s = (Mat_Schur)mat->data;
	PetscTruth assembled;
	
	PetscFunctionBegin; 
	
	
	if( s->A11 != PETSC_NULL ) {
		MatAssembled( s->A11, &assembled );
		if( assembled == PETSC_FALSE ) {	MatAssemblyEnd( s->A11, MAT_FINAL_ASSEMBLY ); 	}
	}
	if( s->A12 != PETSC_NULL ) {
		MatAssembled( s->A12, &assembled );
		if( assembled == PETSC_FALSE ) {	MatAssemblyEnd( s->A12, MAT_FINAL_ASSEMBLY ); 	}
	}
	if( s->A21 != PETSC_NULL ) {
		MatAssembled( s->A21, &assembled );
		if( assembled == PETSC_FALSE ) {	MatAssemblyEnd( s->A21, MAT_FINAL_ASSEMBLY ); 	}
	}
	if( s->A22 != PETSC_NULL ) {
		MatAssembled( s->A22, &assembled );
		if( assembled == PETSC_FALSE ) {	MatAssemblyEnd( s->A22, MAT_FINAL_ASSEMBLY ); 	}
	}
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatCompress_Schur"
PetscErrorCode MatCompress_Schur(Mat mat)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSetOption_Schur"
PetscErrorCode MatSetOption_Schur(Mat mat,MatOption op)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatZeroEntries_Schur"
PetscErrorCode MatZeroEntries_Schur(Mat mat)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*25*/
#undef __FUNCT__
#define __FUNCT__ "MatZeroRows_Schur"
PetscErrorCode MatZeroRows_Schur(Mat mat,PetscInt numRows,const PetscInt rows[],PetscScalar diag)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatLUFactorSymbolic_Schur"
PetscErrorCode MatLUFactorSymbolic_Schur(Mat mat,IS row,IS col,MatFactorInfo *info,Mat *fact)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatLUFactorNumeric_Schur"
PetscErrorCode MatLUFactorNumeric_Schur(Mat mat,MatFactorInfo *info,Mat *fact)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatCholeskyFactorSymbolic_Schur"
PetscErrorCode MatCholeskyFactorSymbolic_Schur(Mat mat,IS perm,MatFactorInfo *info,Mat *fact)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatCholeskyFactorNumeric_Schur"
PetscErrorCode MatCholeskyFactorNumeric_Schur(Mat mat,MatFactorInfo *info,Mat *fact)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*30*/
#undef __FUNCT__
#define __FUNCT__ "MatSetUpPreallocation_Schur"
PetscErrorCode MatSetUpPreallocation_Schur(Mat B)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatILUFactorSymbolic_Schur"
PetscErrorCode MatILUFactorSymbolic_Schur(Mat mat,IS row,IS col,MatFactorInfo *info,Mat *fact)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatICCFactorSymbolic_Schur"
PetscErrorCode MatICCFactorSymbolic_Schur(Mat mat,IS perm,MatFactorInfo *info,Mat *fact)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatGetArray_Schur"
PetscErrorCode MatGetArray_Schur(Mat mat,PetscScalar *v[])
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatRestoreArray_Schur"
PetscErrorCode MatRestoreArray_Schur(Mat mat,PetscScalar *v[])
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*35*/
#undef __FUNCT__
#define __FUNCT__ "MatDuplicate_Schur"
PetscErrorCode MatDuplicate_Schur(Mat mat,MatDuplicateOption op,Mat *M)
{
	
	PetscFunctionBegin; 
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatForwardSolve_Schur"
PetscErrorCode MatForwardSolve_Schur(Mat mat,Vec b,Vec x)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatBackwardSolve_Schur"
PetscErrorCode MatBackwardSolve_Schur(Mat mat,Vec b,Vec x)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatILUFactor_Schur"
PetscErrorCode MatILUFactor_Schur(Mat mat,IS row,IS col,MatFactorInfo *info)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatICCFactor_Schur"
PetscErrorCode MatICCFactor_Schur(Mat mat,IS row,MatFactorInfo* info)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*40*/
#undef __FUNCT__
#define __FUNCT__ "MatAXPY_Schur"
PetscErrorCode MatAXPY_Schur(Mat Y,PetscScalar a,Mat X,MatStructure str)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatGetSubMatrices_Schur"
PetscErrorCode MatGetSubMatrices_Schur(Mat mat,PetscInt n,const IS irow[],const IS icol[],MatReuse scall,Mat *submat[])
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatIncreaseOverlap_Schur"
PetscErrorCode MatIncreaseOverlap_Schur(Mat mat,PetscInt n,IS is[],PetscInt ov)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatGetValues_Schur"
PetscErrorCode MatGetValues_Schur(Mat mat,PetscInt m,const PetscInt idxm[],PetscInt n,const PetscInt idxn[],PetscScalar v[])
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatCopy_Schur"
PetscErrorCode MatCopy_Schur(Mat A,Mat B,MatStructure str)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*45*/
#undef __FUNCT__
#define __FUNCT__ "MatGetRowMax_Schur"
PetscErrorCode MatGetRowMax_Schur(Mat mat,Vec v,PetscInt idx[])
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatScale_Schur"
PetscErrorCode MatScale_Schur(Mat mat,PetscScalar a)
{
	PetscFunctionBegin;
	
	MatSchurSetScalar(mat,a);
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatShift_Schur"
PetscErrorCode MatShift_Schur( Mat mat, PetscScalar a )
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatDiagonalSet_Schur"
PetscErrorCode MatDiagonalSet_Schur( Mat Y, Vec D, InsertMode is )
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatILUDTFactor_Schur"
PetscErrorCode MatILUDTFactor_Schur(Mat mat,IS row,IS col,MatFactorInfo *info,Mat *fact)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*50*/
#undef __FUNCT__
#define __FUNCT__ "MatSetBlockSize_Schur"
PetscErrorCode MatSetBlockSize_Schur(Mat mat,PetscInt bs)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatGetRowIJ_Schur"
PetscErrorCode MatGetRowIJ_Schur(Mat mat,PetscInt shift,PetscTruth symmetric,PetscTruth blockcompressed,PetscInt *n,PetscInt *ia[],PetscInt* ja[],PetscTruth *done)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatRestorRowIJ_Schur"
PetscErrorCode MatRestorRowIJ_Schur(Mat mat,PetscInt shift,PetscTruth symmetric,PetscTruth blockcompressed,PetscInt *n,PetscInt *ia[],PetscInt* ja[],PetscTruth *done)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatGetColumnIJ_Schur"
PetscErrorCode MatGetColumnIJ_Schur(Mat mat,PetscInt shift,PetscTruth symmetric,PetscTruth blockcompressed,PetscInt *n,PetscInt *ia[],PetscInt* ja[],PetscTruth *done)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatRestoreColumnIJ_Schur"
PetscErrorCode MatRestoreColumnIJ_Schur(Mat mat,PetscInt shift,PetscTruth symmetric,PetscTruth blockcompressed,PetscInt *n,PetscInt *ia[],PetscInt* ja[],PetscTruth *done)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*55*/

#undef __FUNCT__
#define __FUNCT__ "MatFDColoringCreate_Schur"
PetscErrorCode MatFDColoringCreate_Schur(Mat mat,ISColoring iscoloring,MatFDColoring color)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatColoringPatch_Schur"
PetscErrorCode MatColoringPatch_Schur(Mat mat,PetscInt ncolors,PetscInt n,ISColoringValue colorarray[],ISColoring *iscoloring)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSetUnfactored_Schur"
PetscErrorCode MatSetUnfactored_Schur(Mat mat)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatPermute_Schur"
PetscErrorCode MatPermute_Schur( Mat mat, IS row, IS col, Mat *B )
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSetValuesBlocked_Schur"
PetscErrorCode MatSetValuesBlocked_Schur(Mat mat,PetscInt m,const PetscInt idxm[],PetscInt n,const PetscInt idxn[],const PetscScalar v[],InsertMode addv)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*60*/
#undef __FUNCT__
#define __FUNCT__ "MatGetSubMatrix_Schur"
PetscErrorCode MatGetSubMatrix_Schur(Mat mat,IS isrow,IS iscol,PetscInt csize,MatReuse cll,Mat *newmat)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
/*
#undef __FUNCT__
#define __FUNCT__ "MatDestroy_Schur"
PetscErrorCode MatDestroy_Schur( Mat mat )
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
*/
#undef __FUNCT__
#define __FUNCT__ "MatView_Schur"
PetscErrorCode MatView_Schur( Mat mat, PetscViewer viewer )
{
	Mat_Schur s = (Mat_Schur)mat->data;
	PetscTruth isascii;
	PetscViewerFormat format;
	
	Stg_PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&isascii);
	if (isascii) {
		
		PetscViewerGetFormat(viewer,&format);
		PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_INFO);
		
		PetscViewerASCIIPrintf(viewer,"complement type: %s\n", s->stype );
		if (s->scale_set) {
			PetscViewerASCIIPrintf(viewer,"scaling factor=%lf\n", s->alpha );
		}
		
		if (s->A11) {
			PetscViewerASCIIPrintf(viewer,"A11:---------------------\n");
			MatView(s->A11,viewer);
		}
		if (s->A12) {
			PetscViewerASCIIPrintf(viewer,"A12:---------------------\n");
			MatView(s->A12,viewer);
		}
		if (s->A21) {
			PetscViewerASCIIPrintf(viewer,"A21:---------------------\n");
			MatView(s->A21,viewer);
		}
		if (s->A22) {
			PetscViewerASCIIPrintf(viewer,"A22:---------------------\n");
			MatView(s->A22,viewer);
		}
		
		if (s->ksp) {
			PetscViewerASCIIPrintf(viewer,"Inner KSP:---------------------\n");
			KSPView(s->ksp,viewer);
		}
		
		PetscViewerSetFormat(viewer,format);
	}
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatConvertFrom_Schur"
PetscErrorCode MatConvertFrom_Schur(Mat mat, MatType newtype,MatReuse reuse,Mat *M)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatUseScaledForm_Schur"
PetscErrorCode MatUseScaledForm_Schur( Mat mat, PetscTruth scaled )
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*65*/
#undef __FUNCT__
#define __FUNCT__ "MatScaleSystem_Schur"
PetscErrorCode MatScaleSystem_Schur( Mat mat, Vec b, Vec x )
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatUnScaleSystem_Schur"
PetscErrorCode MatUnScaleSystem_Schur(Mat mat,Vec b,Vec x)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSetLocalToGlobalMapping_Schur"
PetscErrorCode MatSetLocalToGlobalMapping_Schur(Mat x,ISLocalToGlobalMapping mapping)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSetValuesLocal_Schur"
PetscErrorCode MatSetValuesLocal_Schur(Mat mat,PetscInt nrow,const PetscInt irow[],PetscInt ncol,const PetscInt icol[],const PetscScalar y[],InsertMode addv)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatZeroRowsLocal_Schur"
PetscErrorCode MatZeroRowsLocal_Schur(Mat mat,PetscInt numRows,const PetscInt rows[],PetscScalar diag)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*70*/
#undef __FUNCT__
#define __FUNCT__ "MatGetRowMaxAbs_Schur"
PetscErrorCode MatGetRowMaxAbs_Schur(Mat mat,Vec v,PetscInt idx[])
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatConvert_Schur"
PetscErrorCode MatConvert_Schur(Mat mat, MatType newtype,MatReuse reuse,Mat *M)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSetColoring_Schur"
PetscErrorCode MatSetColoring_Schur( Mat mat, ISColoring colour )
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSetValuesAdic_Schur"
PetscErrorCode MatSetValuesAdic_Schur(Mat mat,void *v)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSetValuesAdifor_Schur"
PetscErrorCode MatSetValuesAdifor_Schur(Mat mat,PetscInt nl,void *v)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}


/*75*/
#undef __FUNCT__
#define __FUNCT__ "MatFDColoringApply_Schur"
PetscErrorCode MatFDColoringApply_Schur(Mat J,MatFDColoring coloring,Vec x1,MatStructure *flag,void *sctx)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSetFromOptions_Schur"
PetscErrorCode MatSetFromOptions_Schur( Mat mat )
{
	PetscTruth     flg;
	PetscScalar    val;
	char           type[PETSC_MAX_PATH_LEN];
    PetscErrorCode ierr;
	
	PetscFunctionBegin; 
	ierr = PetscOptionsBegin(PETSC_COMM_WORLD, PETSC_NULL, "Schur Options", "PC");CHKERRQ(ierr);
	PetscOptionsScalar("-mat_schur_scalar","Scaling factor","MatSchurSetScalar",1.0,&val,&flg);
	if (flg) {			MatSchurSetScalar(mat,val);					}
	PetscOptionsString("-mat_schur_complement_type","Schur complemenet","MatSchurSetSchurComplementType",_MatSchurComplementType[0],type,PETSC_MAX_PATH_LEN-1,&flg);
	if (flg) {			MatSchurSetSchurComplementType(mat,type);								}
	ierr = PetscOptionsEnd();CHKERRQ(ierr);
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatMultConstrained_Schur"
PetscErrorCode MatMultConstrained_Schur(Mat mat,Vec x,Vec y)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatMultTransposeConstrained_Schur"
PetscErrorCode MatMultTransposeConstrained_Schur(Mat mat,Vec x,Vec y)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatILUFactorSymbolicConstrained_Schur"
PetscErrorCode MatILUFactorSymbolicConstrained_Schur(Mat mat,IS row,IS col,double d,PetscInt i,PetscInt j,Mat *fact)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}


/*80*/
#undef __FUNCT__
#define __FUNCT__ "MatPermuteSparsify_Schur"
PetscErrorCode MatPermuteSparsify_Schur(Mat A, PetscInt band, PetscReal frac, PetscReal tol, IS rowp, IS colp, Mat *B)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatMults_Schur"
PetscErrorCode MatMults_Schur(Mat A,Vecs x,Vecs y)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSolves_Schur"
PetscErrorCode MatSolves_Schur(Mat A,Vecs b,Vecs x)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatGetInertia_Schur"
PetscErrorCode MatGetInertia_Schur(Mat mat,PetscInt *nneg,PetscInt *nzero,PetscInt *npos)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatLoad_Schur"
PetscErrorCode MatLoad_Schur(PetscViewer viewer, MatType outtype,Mat *newmat)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}


/*85*/
#undef __FUNCT__
#define __FUNCT__ "MatIsSymmetric_Schur"
PetscErrorCode MatIsSymmetric_Schur(Mat A,PetscReal tol,PetscTruth *flg)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatIsHermitian_Schur"
PetscErrorCode MatIsHermitian_Schur(Mat A,PetscTruth *flg)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatIsStructurallySymmetric_Schur"
PetscErrorCode MatIsStructurallySymmetric_Schur(Mat A,PetscTruth *flg)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatPBRelax_Schur"
PetscErrorCode MatPBRelax_Schur(Mat mat,Vec b,PetscReal omega,MatSORType flag,PetscReal shift,PetscInt its,PetscInt lits,Vec x)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatGetVecs_Schur"
PetscErrorCode MatGetVecs_Schur(Mat mat,Vec *right,Vec *left)
{
	Mat_Schur s = (Mat_Schur)mat->data;
	PetscTruth flg;
	
	
	PetscFunctionBegin;
	
	PetscStrcmp( s->stype, "MatSchur_A11", &flg );
	if( flg == PETSC_TRUE ) {
		if( left != PETSC_NULL ) {		MatGetVecs( s->A21, PETSC_NULL, left );		}
		if( right != PETSC_NULL ){		MatGetVecs( s->A21, PETSC_NULL, right );	}
	}
	
	PetscStrcmp( s->stype, "MatSchur_A22", &flg );
	if( flg == PETSC_TRUE ) {
		
		if( left != PETSC_NULL ) {		MatGetVecs( s->A12, PETSC_NULL, left );		}
		if( right != PETSC_NULL ){		MatGetVecs( s->A12, PETSC_NULL, right );	}
		
	}
	
	PetscFunctionReturn(0);
}



/*90*/
#undef __FUNCT__
#define __FUNCT__ "MatMatMult_Schur"
PetscErrorCode MatMatMult_Schur(Mat A,Mat B,MatReuse scall,PetscReal fill,Mat *C) 
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatMatMultSymbolic_Schur"
PetscErrorCode MatMatMultSymbolic_Schur(Mat A,Mat B,PetscReal fill,Mat *C) 
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatMatMultNumeric_Schur"
PetscErrorCode MatMatMultNumeric_Schur(Mat A,Mat B,Mat C)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatPtAP_Schur"
PetscErrorCode MatPtAP_Schur(Mat A,Mat P,MatReuse scall,PetscReal fill,Mat *C) 
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatPtAPSymbolic_Schur"
PetscErrorCode MatPtAPSymbolic_Schur(Mat A,Mat P,PetscReal fill,Mat *C) 
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}


/*95*/
#undef __FUNCT__
#define __FUNCT__ "MatPtAPNumeric_Schur"
PetscErrorCode MatPtAPNumeric_Schur(Mat A,Mat P,Mat C) 
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatMatMultTranspose_Schur"
PetscErrorCode MatMatMultTranspose_Schur(Mat A,Mat B,MatReuse scall,PetscReal fill,Mat *C) 
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatMatMultTransposeSymbolic_Schur"
PetscErrorCode MatMatMultTransposeSymbolic_Schur(Mat A,Mat B,PetscReal fill,Mat *C) 
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatMatMultTransposeNumeric_Schur"
PetscErrorCode MatMatMultTransposeNumeric_Schur(Mat A,Mat B,Mat C)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatPtAPSymbolic__SEQAIJ_Schur"
PetscErrorCode MatPtAPSymbolic__SEQAIJ_Schur(Mat A,Mat P,PetscReal fill,Mat *C) 
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}


/*100*/
#undef __FUNCT__
#define __FUNCT__ "MatPtAPNumeric_SEQAIJ_Schur"
PetscErrorCode MatPtAPNumeric_SEQAIJ_Schur(Mat A,Mat P,Mat C) 
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatPtAPSymbolic_MPIAIJ_Schur"
PetscErrorCode MatPtAPSymbolic_MPIAIJ_Schur(Mat A,Mat P,PetscReal fill,Mat *C) 
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatPtAPNumeric_MPIAIJ_Schur"
PetscErrorCode MatPtAPNumeric_MPIAIJ_Schur(Mat A,Mat P,Mat C) 
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatConjugate_Schur"
PetscErrorCode MatConjugate_Schur(Mat mat)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSetSizes_Schur"
PetscErrorCode MatSetSizes_Schur(Mat mat,PetscInt m,PetscInt n,PetscInt M,PetscInt N)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}


/*105*/
#undef __FUNCT__
#define __FUNCT__ "MatSetValuesRow_Schur"
PetscErrorCode MatSetValuesRow_Schur(Mat mat,PetscInt row,const PetscScalar v[])
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatRealPart_Schur"
PetscErrorCode MatRealPart_Schur(Mat mat)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatImaginaryPart_Schur"
PetscErrorCode MatImaginaryPart_Schur(Mat mat)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatGetRowUpperTriangular_Schur"
PetscErrorCode MatGetRowUpperTriangular_Schur(Mat mat)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatRestoreRowUpperTriangular_Schur"
PetscErrorCode MatRestoreRowUpperTriangular_Schur(Mat mat)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}


/*110*/
#undef __FUNCT__
#define __FUNCT__ "MatMatSolve_Schur"
PetscErrorCode MatMatSolve_Schur(Mat mat,Mat b,Mat x)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatGetRedundantMatrix_Schur"
PetscErrorCode MatGetRedundantMatrix_Schur(Mat mat,PetscInt nsubcomm,MPI_Comm subcomm,PetscInt mlocal_red,MatReuse reuse,Mat *matredundant)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatGetRowMin_Schur"
PetscErrorCode MatGetRowMin_Schur(Mat mat,Vec v,PetscInt idx[])
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatGetColumnVector_Schur"
PetscErrorCode MatGetColumnVector_Schur(Mat A,Vec yy,PetscInt col)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

