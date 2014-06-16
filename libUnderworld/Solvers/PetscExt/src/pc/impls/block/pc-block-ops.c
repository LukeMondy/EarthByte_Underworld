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

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
#include <stdlib.h>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscpc.h>
#if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 3) )
  #include <petsc-private/kspimpl.h>
#else
  #include <private/kspimpl.h>
#endif
#if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 3) )
  #include <petsc-private/pcimpl.h>
#else
  #include <private/pcimpl.h>
#endif
#if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 3) )
  #include <petsc-private/matimpl.h>
#else
  #include <private/matimpl.h>
#endif

#include "private/vec/petscvec-block.h"
#include "private/mat/petscmat-block.h"

#include "private/pc/petscpc-block.h"
#include "pc-block-impl.h"

/* externs */
extern const char *PCBlockTypes[];



PetscErrorCode _PCSetUp_Block_DetermineBlockApplyType( PC pc, PCBlockType *type )
{
	PC_Block       s;
	PetscInt       i,j;
	PetscInt       mat_cnt, pc_cnt;
	PetscInt       upper_tri_cnt, lower_tri_cnt;
	
	s = (PC_Block)pc->data;
	
	/* Count number of entries */
	pc_cnt = 0;
	upper_tri_cnt = 0;
	lower_tri_cnt = 0;
	
	for( i=0; i<s->nr; i++ ) {
		/* count diagonals */
		if( s->diag[i] != PETSC_NULL ) {
			pc_cnt++;
		}
		for( j=0; j<s->nc; j++ ) {
			
			if( i == j ) {
				
			}
			else {
				
				if( (i < j) && (s->mat[i][j] != PETSC_NULL) ) {
					/* upper triangular */
					upper_tri_cnt++;
				}
				else if( (i > j) && (s->mat[i][j] != PETSC_NULL) ) {
					/* lower triangular */
					lower_tri_cnt++;
				}
				
			}
		}
	}
	mat_cnt = upper_tri_cnt + lower_tri_cnt;
	
	/* Look for an error */
	if( pc_cnt != s->nr ) {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "PCSetUP_Block requires all PC's to be not PETSC_NULL \n" );
	}
	
	
	/* Look for any special cases */
	
	/* 1) Full block inversion for 2x2 block */
	if( s->nr == 2 ) {
		mat_cnt = 0;
		pc_cnt = 0;
		
		if( pc_cnt == 2 && mat_cnt == 2 ) {
			*type = PC_BLOCK_FULL;
		}
		PetscFunctionReturn(0);
	}
	
	
	
	/* Look for standard types */
	
	/* Diagonal */
	if( upper_tri_cnt == 0 && lower_tri_cnt == 0 ) {
		*type = PC_BLOCK_DIAGONAL;
		PetscFunctionReturn(0);
	}
	
	/* Upper */
	if( upper_tri_cnt != 0 && lower_tri_cnt == 0 ) {
		*type = PC_BLOCK_UPPER;
		PetscFunctionReturn(0);
	}
	else if( lower_tri_cnt != 0 && upper_tri_cnt == 0 ) {
		*type = PC_BLOCK_LOWER;
		PetscFunctionReturn(0);
	}
	else if( upper_tri_cnt != 0 && lower_tri_cnt != 0 ) {
		/* Look for an error */
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "PCSetUP_Block requires one of the following configurations for general size blocks \n"
				"1) No off diagonal matrices, => PC_BLOCK_DIAGONAL \n"
				"2) Only matrices in the upper triangular block, => PC_BLOCK_UPPER \n"
				"3) Only matrices in the lower triangular block, => PC_BLOCK_LOWER \n"
				"User has provided matrices in both the upper and lower triangular blocks.\n"
				"Revise configuration of object. \n" );
	}
	
	PetscFunctionReturn(0);
}
/* 
Ensure all pc's have 
0) PC's are not NULL
1) had operators set on them
2) had some pc type specified

Test 1) is not appriopriate for something like stokes with a NULL (1,1) block.

*/
PetscErrorCode _PCSetUp_CheckIntegrityOfPCs( PC pc )
{
	PC_Block       s;
	PetscInt       i;
	
	s = (PC_Block)pc->data;
	
	/* 0) */
	for( i=0; i<s->nr; i++ ) {
		if( s->diag[i] == PETSC_NULL ) {
			Stg_SETERRQ1( PETSC_ERR_ARG_WRONG, "PCSetUP_Block requires KSP at index %d to be not PETSC_NULL \n", i );
		}
	}
	
	/* 1) */
	/*
	Mat Amat, Pmat;
	MatStructure flg;
	for( i=0; i<s->nr; i++ ) {
		KSPGetOperators( s->diag[i], &Amat, &Pmat, &flg );
		if( Pmat == PETSC_NULL ) {
			Stg_SETERRQ1( PETSC_ERR_ARG_WRONG, "PCSetUP_Block requires PC at index %d to have Pmat specified \n", i );
		}
	}
	*/
	/* 2) */
	/*
	for( i=0; i<s->nr; i++ ) {
		if( !s->diag[i]->type_name ) {
			Stg_SETERRQ1( PETSC_ERR_ARG_WRONG, "PCSetUP_Block requires PC at index %d to have type specified \n", i );
		}
	}
	*/
	PetscFunctionReturn(0);
}

/*
The setup routines PCSetUp_PCApply(), PCSetUp_PCApplyTranspose() will create
only the necessay temporary vectors they need to perform there respective applications,
i.e. PCApply() or PCApplyTranspose(). If we have already called PCApply() and then we
use PCApplyTranspose(), only the additional temporary vectors will be created. In this
case a temp. vector will have been created for every row in the block. But in cases where
we just use PCApply() or  just PCApplyTranspose(), we will save on allocating one vector.
*/
PetscErrorCode PCSetUp_PCBlockApply(PC pc)
{
	PC_Block       s = (PC_Block)pc->data;
	Mat Amat, Pmat;
	MatStructure flg;
	PetscInt i;
	
	
	PetscFunctionBegin; 
	
	
	if( s->application_type == PC_BLOCK_DIAGONAL ) {
		s->apply_setupcalled = PETSC_TRUE;
		PetscFunctionReturn(0);
	}
	
	
	/* Create some temporary vectors */
	if( s->t == PETSC_NULL ) {
		s->t = (Vec*)malloc( sizeof(Vec) * s->nr );
		for( i=0; i<s->nr; i++ ) {
			s->t[i] = PETSC_NULL;
		}
	}
	
	/* 
	There is little overhead for creating the array of vector which is
	the full size of the block, even though one of the s->t[] vectors will
	be null. By allocating the full size it also makse indexing though the 
	vector when we compute A[i][j] x[j] simple as well.
	*/
	if( s->application_type == PC_BLOCK_UPPER ) {
		/* If type == UPPER, leave t[nr-1] null and fetch the rest */
		for( i=0; i<s->nr-1; i++ ) {
			if( s->t[i] == PETSC_NULL ) {
				KSPGetOperators( s->diag[i], &Amat, &Pmat, &flg );
				MatGetVecs( Amat, PETSC_NULL, &s->t[i] );  // {t} = [Amat] {x}
			}
		}
	}
	else if ( s->application_type == PC_BLOCK_LOWER ) {
		/* If type == LOWER, leave t[0] null and fetch the rest */
		for( i=1; i<s->nr; i++ ) {
			if( s->t[i] == PETSC_NULL ) {
				KSPGetOperators( s->diag[i], &Amat, &Pmat, &flg );
				MatGetVecs( Amat, PETSC_NULL, &s->t[i] );
			}
		}
	}
	
	
	s->apply_setupcalled = PETSC_TRUE;
		
	
	PetscFunctionReturn(0);
	
}

PetscErrorCode PCSetUp_PCBlockApplyTranspose(PC pc)
{
	PC_Block       s = (PC_Block)pc->data;
	Mat Amat, Pmat;
	MatStructure flg;
	PetscInt i;
	
	
	PetscFunctionBegin; 
	
	
	if( s->application_type == PC_BLOCK_DIAGONAL ) {
		s->applytranspose_setupcalled = PETSC_TRUE;
		PetscFunctionReturn(0);
	}
	
	
	/* Create some temporary vectors */
	if( s->t_trans == PETSC_NULL ) {
		s->t_trans = (Vec*)malloc( sizeof(Vec) * s->nr );
		for( i=0; i<s->nr; i++ ) {
			s->t_trans[i] = PETSC_NULL;
		}
	}
	
	/*
	Note that the MatGetVecs has the args in the correct order as Amat must be square.
	The off diagonal matrice can be rectangular, but it is only the result which will be
	sotored in the tmp. vector t_trans. For the block system to have off diagonal blocks which
	are compoatible, ALL sub matrices (along any row) MUST have the same number of rows.
	
	If a non-square matrix is allowed, then we should swap the order of args in MatGetVecs below.
	*/
	if( s->application_type == PC_BLOCK_UPPER ) { /* need 1 <= i <= nr-1 */
		for( i=1; i<s->nr; i++ ) {
			if( s->t_trans[i] == PETSC_NULL ) {
				KSPGetOperators( s->diag[i], &Amat, &Pmat, &flg );
				MatGetVecs( Amat, PETSC_NULL, &s->t_trans[i] );	// {x} = [Amat]^T {t_trans}}
			}
		}
	}
	else if ( s->application_type == PC_BLOCK_LOWER ) {
		for( i=0; i<s->nr-1; i++ ) {
			if( s->t_trans[i] == PETSC_NULL ) {
				KSPGetOperators( s->diag[i], &Amat, &Pmat, &flg );
				MatGetVecs( Amat, PETSC_NULL, &s->t_trans[i] );
			}
		}
	}
	
	
	s->applytranspose_setupcalled = PETSC_TRUE;
		
	
	PetscFunctionReturn(0);
	
}





/* 0 */
PetscErrorCode PCSetUp_Block(PC pc)
{
	PC_Block       s = (PC_Block)pc->data;
	Mat Amat, Pmat;
	MatStructure flg;
	PetscTruth same;
	PetscInt i,j;//, size;
	MPI_Comm comm;
	Mat mat;
	PetscTruth wasSetup;
	char *sub_name;
	
	
	PetscFunctionBegin; 
	
	/* Now allocate space */
	PCGetOperators( pc, &Amat, &Pmat, &flg );
	/* check Pmat is block */
	Stg_PetscTypeCompare( (PetscObject)Pmat, "block", &same );
	if( same != PETSC_TRUE ) {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "PCBlock requires Pmat to be of type \"block\" \n" );
	}
		
	if( !pc->setupcalled ) {
		wasSetup = PETSC_FALSE;
		
		/* Get block size */
		MatGetSize( Pmat, &s->nr, &s->nc );
		if( s->nr != s->nc ) {
			Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "PCBlock requires Pmat must be square \n" );
		}
		
	#if 0

		/* Allocate space for matrices. The diagonals mat[i][i] will always remain empty */
		size = s->nr * s->nc;
		
		/*
		Message to self....
		The code below to allocate a contiguous block of mem works,
		but I'm a little nervous about how to free it. I thought
		free( s->mat ) should release the whole lot as this is the pointer
		to the array. It doesn't though, I have to do
		for( i=0; i<nr; i++ ) free( s->mat[i] );
		free( s->mat );
		So I'm not sure I understand this stuff well enough. I'm gonna just 
		keep it simple and do it the way I understand. 
		*/
		/* Allocate and init contiguous block */
		tmp_mat = (Mat*)malloc( sizeof(Mat) * size ); 
		for( k=0; k<size; k++ ) {
			tmp_mat[k] = PETSC_NULL;
		}
		
		/* Define as 2d array */
		s->mat = (Mat**)malloc( sizeof(Mat*) * s->nr );
		for( i=0; i<s->nr; i++ ) {
			s->mat[i] = &tmp_mat[ i*s->nc ];
		}
	#endif
		
		/* Create memory for matrices */
		s->mat = (Mat**)malloc( sizeof(Mat*) * s->nr );
		for( i=0; i<s->nr; i++ ) {
			s->mat[i] = (Mat*)malloc( sizeof(Mat) * s->nc );
			// init
			for( j=0; j<s->nc; j++ ) {
				s->mat[i][j] = PETSC_NULL;
			}
		}
		
		/* Check for matrices to be used in off diagonal regions*/
		for( i=0; i<s->nr; i++ ) {
			for( j=0; j<s->nc; j++ ) {
				MatBlockGetSubMatrix( Pmat, i,j, &mat );
				if( i == j ) continue;
				
				if( mat != PETSC_NULL ) {
					s->mat[i][j] = mat;
				}
				
			}
		}
		
		
		
		
		
		/* 
		Allocate space for pc's and init.
		These block pc's have same communication as Pmat.
		*/
		PetscObjectGetComm( (PetscObject)Pmat, &comm );
		s->diag = (KSP*)malloc( sizeof(KSP) * s->nr );
		for( i=0; i<s->nr; i++ ) {
			const char *prefix;

			KSPCreate( comm, &s->diag[i] );
			KSPSetType(s->diag[i],KSPPREONLY);
			PCGetOptionsPrefix(pc,&prefix);
			KSPSetOptionsPrefix(s->diag[i],prefix);
			asprintf( &sub_name, "pc_block_Q%d%d_", i+1,i+1 );
			KSPAppendOptionsPrefix( s->diag[i], sub_name );
			free( sub_name );
		}
		
		
#if 0		
		/* Find a vector */
		PCGetOperators( pc, &Amat, &Pmat, &flg );
		MatGetVecs( Amat, PETSC_NULL, &s->block_vec );
		VecBlockGetSubVectors( s->block_vec, &s->t );
#endif
	}
	else {
		wasSetup = PETSC_TRUE;
	}
	
#if 0	
	_PCSetUp_Block_DetermineBlockApplyType( pc, &s->application_type );
	/* Set the appropriate methods */
	if( s->application_type == PC_BLOCK_DIAGONAL ) {
		/* default is already set to be DIAGONAL so do nothing */
	}
	else if( s->application_type == PC_BLOCK_UPPER ) {
		pc->ops->apply = PCApply_Block_UPPER;
	}	
	else if( s->application_type == PC_BLOCK_LOWER ) {
		pc->ops->apply = PCApply_Block_LOWER;
	}
	else if( s->application_type == PC_BLOCK_FULL ) {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "PCSetUp_Block: PC_BLOCK_FULL not implemented \n" );
	}
#endif
	
	/* update the operators */
	for( i=0; i<s->nr; i++ ) {
		MatBlockGetSubMatrix( Pmat, i,i, &mat );
		if( mat != PETSC_NULL ) {
			KSPSetOperators( s->diag[i], mat, mat, pc->flag );
		}
	}
	
	if( !wasSetup && pc->setfromoptionscalled ) {
		for( i=0; i<s->nr; i++ ) {
			KSPSetFromOptions( s->diag[i] );
		}
	}
	
	_PCSetUp_CheckIntegrityOfPCs( pc );
	

	
	PetscFunctionReturn(0);
}


PetscErrorCode PCApply_Block_UPPER(PC pc,Vec x,Vec y)
{
	PC_Block       s = (PC_Block)pc->data;
	PetscInt       i,j,M,N;
	Vec            *sub_x, *sub_y;
	
	
	PetscFunctionBegin; 
	if( !s->apply_setupcalled ) {
		PCSetUp_PCBlockApply(pc);
	}
	
	
	/* lower triangular */
	M = s->nr;
	N = s->nc;
	
	VecBlockGetSubVectors( y, &sub_y );
	VecBlockGetSubVectors( x, &sub_x );
	
	
	/* Do the last row as it has no coupling */
	i = M-1;
	/*printf("PCApply[Q_%d, t_%d, y%d] \n", i,i,i ); */
	KSPSolve( s->diag[i], sub_x[i], sub_y[i] );
	
	/* Do the other rows from 0 <= i < M-1 */
	for( i=M-2; i>=0; i-- ) {
		VecZeroEntries( s->t[i] );
		for( j=N-1; j>=i+1; j-- ) {
			
			if( s->mat[i][j] != PETSC_NULL ){
				MatMultAdd( s->mat[i][j], sub_y[j], s->t[i], s->t[i] );
				/*printf("t_%d <= t_%d + A_{%d,%d} y_%d \n",	i,i, i,j, j ); */
			}
		}
		/*printf("t_%d <= x_%d - t_%d \n", i,i,i );*/
		VecAYPX( s->t[i], -1.0, sub_x[i] );
		
		/*printf("PCApply[Q_%d, t_%d, y%d] \n", i,i,i );*/
		KSPSolve( s->diag[i], s->t[i], sub_y[i] );
	} 
	
	VecBlockRestoreSubVectors( y );
	VecBlockRestoreSubVectors( x );
	
	PetscFunctionReturn(0);
}

PetscErrorCode PCApply_Block_LOWER(PC pc,Vec x,Vec y)
{
	PC_Block       s = (PC_Block)pc->data;
	PetscInt       i,j,M;
	Vec            *sub_x, *sub_y;
	
	
	PetscFunctionBegin; 
	if( !s->apply_setupcalled ) {
		PCSetUp_PCBlockApply(pc);
	}
	
	
	/* upper triangular */
	M = s->nr;
	//N = s->nc;
	
	VecBlockGetSubVectors( y, &sub_y );
	VecBlockGetSubVectors( x, &sub_x );
	
	
	
	/* Lower triangular */
	/* Do the first row as it has no coupling */
	i = 0;
	/*printf("PCApply[Q_%d, t_%d, y%d] \n", i,i,i );*/
	KSPSolve( s->diag[i], sub_x[i], sub_y[i] );
	
	/* Do the other rows from 1 <= i < M-1 */
	for( i=1; i<M; i++ ) {
		VecZeroEntries( s->t[i] );
		for( j=0; j<=i-1; j++ ) {
			
			if( s->mat[i][j] != PETSC_NULL ){
				MatMultAdd( s->mat[i][j], sub_y[j], s->t[i], s->t[i] );
				/*printf("t_%d <= t_%d + A_{%d,%d} y_%d \n",		i,i, i,j, j );*/
			}
		}
		/*printf("t_%d <= x_%d - t_%d \n", i,i,i );*/
		VecAYPX( s->t[i], -1.0, sub_x[i] );
		
		/*printf("PCApply[Q_%d, t_%d, y%d] \n", i,i,i );*/
		KSPSolve( s->diag[i], s->t[i], sub_y[i] );
	} 
	
	VecBlockRestoreSubVectors( y );
	VecBlockRestoreSubVectors( x );
	
	PetscFunctionReturn(0);
}





PetscErrorCode PCApply_Block_DIAGONAL(PC pc,Vec x,Vec y)
{
	PC_Block       s = (PC_Block)pc->data;
	PetscInt       i,M;
	Vec            *sub_x, *sub_y;
	
	
	PetscFunctionBegin;
	if( !s->apply_setupcalled ) {
		PCSetUp_PCBlockApply(pc);
	}
	
	
	M = s->nr;
	
	VecBlockGetSubVectors( y, &sub_y );
	VecBlockGetSubVectors( x, &sub_x );
	
	for( i=0; i<M; i++ ) {
		KSPSolve( s->diag[i], sub_x[i], sub_y[i] );
	} 
	
	VecBlockRestoreSubVectors( y );
	VecBlockRestoreSubVectors( x );
	
	PetscFunctionReturn(0);
}
PetscErrorCode PCApplyRichardson_Block(PC pc,Vec x,Vec y,Vec w,PetscReal rtol,PetscReal abstol, PetscReal dtol,PetscInt its)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	PetscFunctionReturn(0);
}
PetscErrorCode PCApplyBA_Block(PC pc,PCSide side,Vec x,Vec y,Vec work)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	PetscFunctionReturn(0);
}


PetscErrorCode PCApplyTranspose_Block_LOWER(PC pc,Vec x,Vec y)
{
	PC_Block       s = (PC_Block)pc->data;
	PetscInt       i,j,M,N;
	Vec            *sub_x, *sub_y;
	Vec            *TMP;
	
	PetscFunctionBegin; 
	if( !s->applytranspose_setupcalled ) {
		PCSetUp_PCBlockApplyTranspose(pc);
	}
	TMP = s->t_trans;
	
	/* upper triangular^T => lower triangular */
	M = s->nr;
	N = s->nc;
	
	VecBlockGetSubVectors( y, &sub_y );
	VecBlockGetSubVectors( x, &sub_x );
	
	
	/* Do the last row as it has no coupling */
	i = M-1;
	KSPSolveTranspose( s->diag[i], sub_x[i], sub_y[i] );
	
	/* Do the other rows from 0 <= i < M-1 */
	for( i=M-2; i>=0; i-- ) {
		VecZeroEntries( TMP[i] );
		for( j=N-1; j>=i+1; j-- ) {
			
			if( s->mat[j][i] != PETSC_NULL ){
				MatMultTransposeAdd( s->mat[j][i], sub_y[j], TMP[i], TMP[i] );
			}
		}
		VecAYPX( TMP[i], -1.0, sub_x[i] );
		
		KSPSolveTranspose( s->diag[i], TMP[i], sub_y[i] );
	} 
	
	VecBlockRestoreSubVectors( y );
	VecBlockRestoreSubVectors( x );
	
	PetscFunctionReturn(0);
}


PetscErrorCode PCApplyTranspose_Block_UPPER(PC pc,Vec x,Vec y)
{
	PC_Block       s = (PC_Block)pc->data;
	PetscInt       i,j,M;
	Vec            *sub_x, *sub_y;
	Vec            *TMP;
	
	PetscFunctionBegin; 
	if( !s->applytranspose_setupcalled ) {
		PCSetUp_PCBlockApplyTranspose(pc);
	}
	TMP = s->t_trans;
	
	/* upper triangular */
	M = s->nr;
	//N = s->nc;
	
	VecBlockGetSubVectors( y, &sub_y );
	VecBlockGetSubVectors( x, &sub_x );
	
	
	/* Lower triangular */
	/* Do the first row as it has no coupling */
	i = 0;
	KSPSolveTranspose( s->diag[i], sub_x[i], sub_y[i] );
	
	/* Do the other rows from 1 <= i < M-1 */
	for( i=1; i<M; i++ ) {
		VecZeroEntries( TMP[i] );
		for( j=0; j<=i-1; j++ ) {
			
			if( s->mat[j][i] != PETSC_NULL ){
				MatMultTransposeAdd( s->mat[j][i], sub_y[j], TMP[i], TMP[i] );
			}
		}
		VecAYPX( TMP[i], -1.0, sub_x[i] );
		
		KSPSolveTranspose( s->diag[i], TMP[i], sub_y[i] );
	} 
	
	VecBlockRestoreSubVectors( y );
	VecBlockRestoreSubVectors( x );
	
	PetscFunctionReturn(0);
}


PetscErrorCode PCApplyTranspose_Block_DIAGONAL(PC pc,Vec x,Vec y)
{
	PC_Block       s = (PC_Block)pc->data;
	PetscInt       i,M;
	Vec            *sub_x, *sub_y;
	
	
	PetscFunctionBegin;
	if( !s->applytranspose_setupcalled ) {
		PCSetUp_PCBlockApplyTranspose(pc);
	}
	
	M = s->nr;
	
	VecBlockGetSubVectors( y, &sub_y );
	VecBlockGetSubVectors( x, &sub_x );
	
	for( i=0; i<M; i++ ) {
		KSPSolveTranspose( s->diag[i], sub_x[i], sub_y[i] );
	} 
	
	VecBlockRestoreSubVectors( y );
	VecBlockRestoreSubVectors( x );
	
	PetscFunctionReturn(0);
}


/* 5 */
PetscErrorCode PCApplyBATranspose_Block(PC pc,PetscInt m,Vec x,Vec y,Vec work)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	PetscFunctionReturn(0);
}

/*
Activates the following options
-pc_block_type DIAGONAL,UPPER,LOWER

*/
PetscErrorCode PCSetFromOptions_Block(PC pc)
{
	PC_Block       s = (PC_Block)pc->data;
	PetscTruth     flg,same;
	PCBlockType    type;
	Mat            Amat, Pmat;
	MatStructure   str;
    PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	ierr = PetscOptionsBegin(PETSC_COMM_WORLD, PETSC_NULL, "PC Block Options", "PC");CHKERRQ(ierr);
	
	PCGetOperators( pc, &Amat, &Pmat, &str );
	/* check Pmat is block */
	Stg_PetscTypeCompare( (PetscObject)Pmat, "block", &same );
	if( same != PETSC_TRUE ) {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "PCBlock requires Pmat to be of type \"block\" \n" );
	}
	/* Get block size */
	MatGetSize( Pmat, &s->nr, &s->nc );
	if( s->nr != s->nc ) {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "PCBlock requires Pmat must be square \n" );
	}
	
	
	PetscOptionsEnum("-pc_block_type","Specifies the block structure of the preconditioner","PCBlockSetBlockType",
			PCBlockTypes, (PetscEnum)s->application_type,
			(PetscEnum*)&type, &flg );
	if (flg) {
		PCBlockSetBlockType( pc, type );
	}

	ierr = PetscOptionsEnd();CHKERRQ(ierr);
	PetscFunctionReturn(0);
}
PetscErrorCode PCPreSolve_Block(PC pc,KSP ksp,Vec rhs,Vec x)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	PetscFunctionReturn(0);
}
PetscErrorCode PCPostSolve_Block(PC pc,KSP ksp,Vec rhs,Vec x)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	PetscFunctionReturn(0);
}
PetscErrorCode PCGetFactoredMatrix_Block(PC pc,Mat *A)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	PetscFunctionReturn(0);
}
/* 10 */
PetscErrorCode PCApplySymmetricLeft_Block(PC pc,Vec x,Vec y)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	PetscFunctionReturn(0);
}
PetscErrorCode PCApplySymmetricRight_Block(PC pc,Vec x,Vec y)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	PetscFunctionReturn(0);
}
PetscErrorCode PCSetUpOnBlocks_Block(PC pc)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	PetscFunctionReturn(0);
}
/*
PetscErrorCode Stg_PCDestroy_Block(PC pc)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	PetscFunctionReturn(0);
}
*/


/*
Prints memory address or name if one exists
*/
#define _WRITE_NAME( i,j, pc, viewer ) \
if( (pc) == PETSC_NULL ) { \
	PetscViewerASCIIPrintf( (viewer),"(%d,%d)  %p ",(i),(j), (pc) ); \
} \
else if ( (pc)->prefix) { \
	PetscViewerASCIIPrintf( (viewer),"(%d,%d)  %s ", (i),(j), (pc)->prefix); \
} else { \
	PetscViewerASCIIPrintf( (viewer),"(%d,%d)  %p ",(i),(j), (pc) ); \
} \
		

#define PRINT_OBJ_IDNETIFER( type, i,j, obj, viewer ) \
if( (obj) != PETSC_NULL ) { \
	if( (obj)->prefix ) PetscViewerASCIIPrintf( (viewer),"[%s] (%d,%d)  %s\n", (type), (i),(j), (obj)->prefix); \
	else                PetscViewerASCIIPrintf( (viewer),"[%s] (%d,%d)  %p\n", (type), (i),(j), (obj) ); \
} else { \
    PetscViewerASCIIPrintf( (viewer),"[%s] (%d,%d)  %p\n", (type), (i),(j), (obj) ); \
} \


/* This is NOT recursive */
void print_pc_structure( PC pc, PetscViewer viewer )
{
	PC_Block	bA = (PC_Block)pc->data;
	PCType		pc_type;
	KSPType		ksp_type;
	MatType		mat_type;
	PetscInt	m,n;
	int i,j;
	const char *mat_name, *pc_name;
	Mat Aij;
	KSP Qii;
	PC _Qii;
	
	
	if(bA->diag==PETSC_NULL) {
		PetscViewerASCIIPrintf( viewer, "KSP list is not allocated <structure>\n" ); 
		return;
	}
	
	for( i=0; i<bA->nr; i++ ) {
		for( j=0; j<bA->nc; j++ ) {
			
			if( i== j ) { // PC
				Qii = bA->diag[i];
				
				if( Qii == PETSC_NULL ) {
					PetscViewerASCIIPrintf( viewer, "[KSP]  (%d,%d) - (PETSC_NULL) \n", i,j ); 
					continue;
				}
				pc_name = ((PetscObject)Qii)->prefix;
				KSPGetType( Qii, &ksp_type );
				KSPGetPC( Qii, &_Qii );
				PCGetType( _Qii, &pc_type );

				if( pc_name == PETSC_NULL ) { PetscViewerASCIIPrintf( viewer, "[KSP/PC] (%d,%d) - (%p), type=%s/%s \n", i,j,Qii, ksp_type, pc_type ); }
				else { PetscViewerASCIIPrintf( viewer, "[KSP/PC]  (%d,%d) - (%s), type=%s/%s \n", i,j, pc_name, ksp_type, pc_type ); }
				
				
			}
			else { // Mat
				Aij = bA->mat[i][j];
				
				if( Aij == PETSC_NULL ) {
					PetscViewerASCIIPrintf( viewer, "[Mat] (%d,%d) - (PETSC_NULL) \n", i,j ); 
					continue;
				}
				mat_name = ((PetscObject)Aij)->prefix;
				MatGetType( Aij, &mat_type );
				MatGetSize( Aij, &m, &n );
				
				if( mat_name == PETSC_NULL ) { PetscViewerASCIIPrintf( viewer, "[Mat] (%d,%d) - (%p), type=%s, rows=%d, cols=%d \n", i,j,Aij, mat_type, m,n ); }
				else { PetscViewerASCIIPrintf( viewer, "[Mat] (%d,%d) - (%s), type=%s, rows=%d, cols=%d \n", i,j, mat_name, mat_type, m,n ); }
			}
		}
	}
	
	
}


void print_pc_contents( PC pc, PetscViewer viewer )
{
	PC_Block	bA = (PC_Block)pc->data;
	PetscInt	m,n;
	int i,j;
	const char *my_pc_name;
	const char *mat_name, *pc_name;
	Mat Aij;
	KSP Qii;
	void *obj;
	
	if(bA->diag==PETSC_NULL) {
		PetscViewerASCIIPrintf( viewer, "KSP list is not allocated <contents>\n" ); 
		return;
	}
	
	
	my_pc_name = ((PetscObject)pc)->prefix;
	
	for( i=0; i<bA->nr; i++ ) {
		for( j=0; j<bA->nc; j++ ) {
			
			PetscViewerASCIIPushTab( viewer );
			
			if( i== j ) {     obj = (void*)bA->diag[i];       }   // PC
			else {            obj = (void*)bA->mat[i][j];     }   // Mat
			
			if( obj != PETSC_NULL ) {
				PetscViewerASCIIPrintf( viewer, "-------------------- Begin ( %p ) --------------------\n", obj );
			}
			
	//		PetscViewerASCIIPrintf( viewer, "(i,j) - [sub TYPE]: (Name / ptr),  [parent pc]: (Name / ptr) \n" );
			
			
			if( my_pc_name != PETSC_NULL ) {
			
				if( i== j ) { // PC
					Qii = bA->diag[i];
					
					if( !Qii ) {   PetscViewerASCIIPrintf( viewer, "(%d,%d) - [sub pc]: (PETSC_NULL), [parent pc]: %s \n", i,j, my_pc_name );    }
					else {
						pc_name = ((PetscObject)Qii)->prefix;
						
						if( !pc_name ) { PetscViewerASCIIPrintf( viewer, "(%d,%d) - [sub pc]: (%p), [parent pc]: %s \n", i,j,Qii, my_pc_name );       }
						else {           PetscViewerASCIIPrintf( viewer, "(%d,%d) - [sub pc]: (%s), [parent pc]: %s \n", i,j, pc_name, my_pc_name );  }
						KSPView( Qii, viewer );
					}
				}
				else { // Mat
					Aij = bA->mat[i][j];
					
					if( !Aij ) {  PetscViewerASCIIPrintf( viewer, "(%d,%d) - [sub mat]: (PETSC_NULL), [parent pc]: %s \n", i,j, my_pc_name );   }
					else {
						mat_name = ((PetscObject)Aij)->prefix;
						MatGetSize( Aij, &m, &n );
						
						if( !mat_name ) { PetscViewerASCIIPrintf( viewer, "(%d,%d) - [sub mat]: (%p), [parent pc]: %s \n", i,j,Aij, my_pc_name );       }
						else {            PetscViewerASCIIPrintf( viewer, "(%d,%d) - [sub mat]: (%s), [parent pc]: %s \n", i,j, mat_name, my_pc_name ); }
						
						PetscViewerSetFormat( viewer, PETSC_VIEWER_ASCII_INFO );
						MatView( Aij, viewer );
						PetscViewerSetFormat( viewer, PETSC_VIEWER_ASCII_COMMON );
					}
				}
				
			}
			else {
				
				if( i== j ) { // PC
					Qii = bA->diag[i];
					
					if( !Qii ) {   PetscViewerASCIIPrintf( viewer, "(%d,%d) - [sub pc]: (PETSC_NULL), [parent pc]: %p \n", i,j, pc );    }
					else {
						pc_name = ((PetscObject)Qii)->prefix;
						
						if( !pc_name ) { PetscViewerASCIIPrintf( viewer, "(%d,%d) - [sub pc]: (%p), [parent pc]: %p \n", i,j,Qii, pc );       }
						else {           PetscViewerASCIIPrintf( viewer, "(%d,%d) - [sub pc]: (%s), [parent pc]: %p \n", i,j, pc_name, pc );  }
						KSPView( Qii, viewer );
					}
				}
				else { // Mat
					Aij = bA->mat[i][j];
					
					if( !Aij ) {  PetscViewerASCIIPrintf( viewer, "(%d,%d) - [sub mat]: (PETSC_NULL), [parent pc]: %p \n", i,j, pc );   }
					else {
						mat_name = ((PetscObject)Aij)->prefix;
						MatGetSize( Aij, &m, &n );
						
						if( !mat_name ) { PetscViewerASCIIPrintf( viewer, "(%d,%d) - [sub mat]: (%p), [parent pc]: %p \n", i,j,Aij, pc );       }
						else {            PetscViewerASCIIPrintf( viewer, "(%d,%d) - [sub mat]: (%s), [parent pc]: %p \n", i,j, mat_name, pc ); }
						
						PetscViewerSetFormat( viewer, PETSC_VIEWER_ASCII_INFO );
						MatView( Aij, viewer );
						PetscViewerSetFormat( viewer, PETSC_VIEWER_ASCII_COMMON );
					}
				}
			
			}
			
			if( obj != PETSC_NULL ) {
				PetscViewerASCIIPrintf( viewer, "--------------------   End ( %p ) --------------------\n\n", obj );
			}
			else {
				PetscViewerASCIIPrintf( viewer, "\n" );
			}
			
			PetscViewerASCIIPopTab( viewer );
			
		}
	}
	
	
}



PetscErrorCode PCView_Block(PC pc,PetscViewer viewer)
{
	PC_Block       s = (PC_Block)pc->data;
	PetscTruth isascii;
		
	
	PetscFunctionBegin; 
	
	Stg_PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&isascii);
	if (isascii) {
	
		PetscViewerASCIIPushTab( viewer );
		PetscViewerASCIIPrintf(viewer,"Block size:   %d x %d \n", s->nr, s->nc );
		PetscViewerASCIIPrintf(viewer,"BlockPC type: %s \n", PCBlockTypes[ s->application_type ] );
		
		if (!pc->setupcalled) {
			PetscViewerASCIIPrintf(viewer,"PCSetUp_Block not yet called. \n" );
		}
		if (!s->apply_setupcalled) {
			PetscViewerASCIIPrintf(viewer,"PCSetUp_PCBlockApply not yet called. \n" );
		}
		if (!s->applytranspose_setupcalled) {
			PetscViewerASCIIPrintf(viewer,"PCSetUp_PCBlockApplyTranspose not yet called. \n" );
		}
	
		PetscViewerASCIIPrintf(viewer,"PCBlock structure: \n" );
		PetscViewerASCIIPrintf( viewer, "[TYPE] (i,j) - (Name / ptr), type=, rows=, cols=\n" );
	
		PetscViewerASCIIPushTab( viewer );
		print_pc_structure( pc, viewer );
		PetscViewerASCIIPopTab( viewer );
	
		PetscViewerASCIIPrintf(viewer,"PCBlock contents: \n" );
		PetscViewerASCIIPrintf( viewer, "(i,j) - [sub TYPE]: (Name / ptr),  [parent pc]: (Name / ptr) \n" );
		PetscViewerASCIIPushTab( viewer );		// push1
		print_pc_contents( pc, viewer );
		PetscViewerASCIIPopTab( viewer );		// pop1
	
	
		PetscViewerASCIIPopTab( viewer );
	}
	PetscFunctionReturn(0);
}




