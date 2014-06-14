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
#include <petscmat.h>
#include <petscvec.h>
#include <private/matimpl.h>


#include "src/vec/impls/block/vec_block_impl.h"
#include "private/vec/petscvec-block.h"

#include "src/mat/impls/block/mat_block_impl.h"
#include "private/mat/petscmat-block.h"



extern const char *MatBlockStructureName[];


PetscErrorCode _check_mat_mat_compatibility( Mat A, Mat B )
{
	PetscTruth isAblock, isBblock;
	
	PetscTypeCompare( (PetscObject)A, "block", &isAblock );
	PetscTypeCompare( (PetscObject)B, "block", &isBblock );
	
	if( isAblock == PETSC_TRUE && isBblock == PETSC_TRUE ) {
		PetscFunctionReturn(0);
	} else {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "Operation only valid for (MAT-MAT) BLOCK combinations.");
		PetscFunctionReturn(0);
	}
}

PetscErrorCode _check_mat_vec_compatibility( Mat A, Vec x )
{
	PetscTruth isAblock, isxblock;
	
	PetscTypeCompare( (PetscObject)A, "block", &isAblock );
	PetscTypeCompare( (PetscObject)x, "block", &isxblock );
	
	if( isAblock == PETSC_TRUE && isxblock == PETSC_TRUE ) {
		PetscFunctionReturn(0);
	} else {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "Operation only valid for (MAT-VEC) BLOCK combinations.");
		PetscFunctionReturn(0);
	}
}

PetscErrorCode _check_mat_vec_compatibility2( Mat A, Vec x, Vec y )
{
	PetscTruth isAblock, isxblock, isyblock;
	
	PetscTypeCompare( (PetscObject)A, "block", &isAblock );
	PetscTypeCompare( (PetscObject)x, "block", &isxblock );
	PetscTypeCompare( (PetscObject)y, "block", &isyblock );
	
	if( 		(isAblock == PETSC_TRUE)
			&& 	(isxblock == PETSC_TRUE)
			&&	(isyblock == PETSC_TRUE) ) {
		PetscFunctionReturn(0);
	} else {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "Operation only valid for (MAT-VEC) BLOCK combinations.");
		PetscFunctionReturn(0);
	}
}

PetscErrorCode _check_mat_vec_compatibility3( Mat A, Vec x, Vec y, Vec w )
{
	PetscTruth isAblock, isxblock, isyblock, iswblock;
	
	PetscTypeCompare( (PetscObject)A, "block", &isAblock );
	PetscTypeCompare( (PetscObject)x, "block", &isxblock );
	PetscTypeCompare( (PetscObject)y, "block", &isyblock );
	PetscTypeCompare( (PetscObject)w, "block", &iswblock );
	
	if( 		(isAblock == PETSC_TRUE)
			&& 	(isxblock == PETSC_TRUE)
			&&	(isyblock == PETSC_TRUE)
			&&	(iswblock == PETSC_TRUE) ) {
		PetscFunctionReturn(0);
	} else {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "Operation %s only valid for (MAT-VEC) BLOCK combinations.");
		PetscFunctionReturn(0);
	}
}





PetscErrorCode MatMult_Block( Mat A, Vec x, Vec y )
{
	Mat_Block *bA = (Mat_Block*)A->data;
	Vec_Block *bx = (Vec_Block*)x->data;
	Vec_Block *by = (Vec_Block*)y->data;
	PetscInt i,j;
	
	
	_check_mat_vec_compatibility2( A,x,y );
	
	VecSet( y, 0.0 );
	
	for( i=0; i<bA->nr; i++ ) {
		for( j=0; j<bA->nc; j++ ) {
			
			if( bA->m[i][j] == PETSC_NULL  || bx->v[j] == PETSC_NULL ) {
				continue;
			}
			
			/* y[i] <- y[i] + A[i][j] * x[j] */
			MatMultAdd( bA->m[i][j],bx->v[j], by->v[i],by->v[i] );
		}
	}
	
	
	PetscFunctionReturn(0);
}


PetscErrorCode MatMultTranspose_Block( Mat A, Vec x, Vec y )
{
	
	Mat_Block *bA = (Mat_Block*)A->data;
	Vec_Block *bx = (Vec_Block*)x->data;
	Vec_Block *by = (Vec_Block*)y->data;
	PetscInt i,j;
	
	
	_check_mat_vec_compatibility2( A,x,y );
	
	if( A->symmetric == PETSC_TRUE ) {
		MatMult_Block(A,x,y);
		PetscFunctionReturn(0);
	}
	
	VecSet( y, 0.0 );
	
	for( i=0; i<bA->nr; i++ ) {
		for( j=0; j<bA->nc; j++ ) {
			
			if( bA->m[j][i] == PETSC_NULL  || bx->v[j] == PETSC_NULL ) {
				continue;
			}
			
			/* y[i] <- y[i] + A^T[i][j] * x[j], so we swap i,j in mat[][] */
			MatMultTransposeAdd( bA->m[j][i], bx->v[j], by->v[i], by->v[i] );
		}
	}
	
	PetscFunctionReturn(0);
}

/*
A clear description of these vectors is;
  {left} = [A] {right}
*/
PetscErrorCode MatGetVecs_Block( Mat A, Vec *right, Vec *left )
{
  //Vec_Block *br;
  //Vec_Block *bl;
	Mat_Block *bA = (Mat_Block*)A->data;
	Vec *L, *R;
	MPI_Comm comm;
	PetscInt i,j;
	
	PetscObjectGetComm( (PetscObject)A, &comm );
	
	
	
	
	if( right != PETSC_NULL ) {
		/* allocate R */
		R = (Vec*)malloc( sizeof(Vec) * bA->nc );
		
		
		/* Create the right vectors */
	//	MatGetVecs( bm->A11, &R[0], PETSC_NULL );
	//	MatGetVecs( bm->A12, &R[1], PETSC_NULL );
		for( j=0; j<bA->nc; j++ ) {
			for( i=0; i<bA->nr; i++ ) {
				if( bA->m[i][j] != PETSC_NULL ) {
					MatGetVecs( bA->m[i][j], &R[j], PETSC_NULL );
					
					break;
				}
			}
			if( i==bA->nr ) {
				/* have an empty column */
				Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "MAT BLOCK contains a null column.");
			}
		}
		
		
		VecCreate( comm, right );
		VecSetSizes( *right, bA->nc, bA->nc );
		VecSetType( *right, "block" );
		//br = (Vec_Block*)( (*right)->data );
		
		for( j=0; j<bA->nc; j++ ) {
		//	VecBlock_SetBlock( *right, j, R[j] );
			VecBlockSetValue( *right, j, R[j], INSERT_VALUES );
			Stg_VecDestroy( &R[j] );
		}
		
		free( R );
	}
	
	if( left != PETSC_NULL ) {
		/* allocate L */
		L = (Vec*)malloc( sizeof(Vec) * bA->nr );
		
		
		/* Create the left vectors */
	//	MatGetVecs( bm->A11, PETSC_NULL, &L[0] );
	//	MatGetVecs( bm->A21, PETSC_NULL, &L[1] );
		
		for( i=0; i<bA->nr; i++ ) {
			for( j=0; j<bA->nc; j++ ) {
				if( bA->m[i][j] != PETSC_NULL ) {
					MatGetVecs( bA->m[i][j], PETSC_NULL, &L[i]  );
					
					break;
				}
			}
			if( j==bA->nc ) {
				/* have an empty row */
				Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "MAT BLOCK contains a null row.");
			}
		}
		
		VecCreate( comm, left );
		VecSetSizes( *left, bA->nr, bA->nr );
		VecSetType( *left, "block" );
		//bl = (Vec_Block*)( (*left)->data );
		
		for( i=0; i<bA->nr; i++ ) {
			VecBlockSetValue( *left, i, L[i], INSERT_VALUES );
			Stg_VecDestroy( &L[i] );
		}
		
		free( L );
	}
	
	
	PetscFunctionReturn(0);
}

/*
Prints memory address or name if one exists
*/
#define _WRITE_NAME( mat, viewer ) \
if( (mat) == PETSC_NULL ) { \
	PetscViewerASCIIPrintf( (viewer),"%p ",(mat) ); \
} \
else if ( (mat)->prefix) { \
	PetscViewerASCIIPrintf( (viewer),"%s ", (mat)->prefix); \
} else { \
	PetscViewerASCIIPrintf( (viewer),"%p ",(mat) ); \
} \
		
#define PRINT_OBJ_IDNETIFER( type, i,j, obj, viewer ) \
if( (obj) != PETSC_NULL ) { \
	if( (obj)->prefix ) PetscViewerASCIIPrintf( (viewer),"[%s] (%d,%d)  %s\n", (type), (i),(j), (obj)->prefix); \
	else                PetscViewerASCIIPrintf( (viewer),"[%s] (%d,%d)  %p\n", (type), (i),(j), (obj) ); \
} else { \
    PetscViewerASCIIPrintf( (viewer),"[%s] (%d,%d)  %p\n", (type), (i),(j), (obj) ); \
} \



void print_mat_structure( Mat A, int pi,PetscInt pj, PetscViewer viewer )
{
	Mat_Block	*bA;
	const MatType		type;
	PetscInt	m,n;
	int i,j;
	const char *name;
	PetscTruth is_block;
	
	
	if( A == PETSC_NULL ) {
		PetscViewerASCIIPushTab( viewer );
		PetscViewerASCIIPrintf( viewer, "(%d,%d) - (PETSC_NULL) \n", pi,pj ); 
		PetscViewerASCIIPopTab( viewer );
		return;
	}
	
	
	MatGetType( A, &type );
	name = ((PetscObject)A)->prefix;
	is_block = PETSC_FALSE;
	PetscTypeCompare( (PetscObject)A, "block", &is_block );
	
	if( is_block == PETSC_FALSE ) {
		MatGetSize( A, &m, &n );
		
		PetscViewerASCIIPushTab( viewer );
		
		if( name == PETSC_NULL ) { PetscViewerASCIIPrintf( viewer, "(%d,%d) - (%p), type=%s, rows=%d, cols=%d \n", pi,pj,A, type, m,n ); }
		else { PetscViewerASCIIPrintf( viewer, "(%d,%d) - (%s), type=%s, rows=%d, cols=%d \n", pi,pj, name, type, m,n ); }
		
		PetscViewerASCIIPopTab( viewer );
		return;
	}
	else {
		MatGetSize( A, &m, &n );
		
		PetscViewerASCIIPushTab( viewer );
		
		if( name == PETSC_NULL ) { PetscViewerASCIIPrintf( viewer, "(%d,%d) - (%p), type=%s, rows=%d, cols=%d \n", pi,pj,A, type, m,n ); }
		else { PetscViewerASCIIPrintf( viewer, "(%d,%d) - (%s), type=%s, rows=%d, cols=%d \n", pi,pj, name, type, m,n ); }
		
		bA = (Mat_Block*)A->data;
		for( i=0; i<bA->nr; i++ ) {
			for( j=0; j<bA->nc; j++ ) {
				print_mat_structure( bA->m[i][j],i,j, viewer );
			}
		}
		PetscViewerASCIIPopTab( viewer );
	}
}

void print_mat_contents( Mat A, Mat parent_A, PetscInt pi, PetscInt pj, PetscViewer viewer )
{
	Mat_Block	*bA;
	const MatType		type;
	PetscInt	m,n;
	int i,j;
	const char *name, *pname;
	PetscTruth is_block;
	
	
	name = PETSC_NULL;
	pname = PETSC_NULL;
	
	if( parent_A != PETSC_NULL ) {
		pname = ((PetscObject)parent_A)->prefix;
	}
	else pname = PETSC_NULL;
	
	if( A != PETSC_NULL ) {
		name = ((PetscObject)A)->prefix;
	}
	else {
		PetscViewerASCIIPushTab( viewer );
		
		if( pname == PETSC_NULL ) PetscViewerASCIIPrintf( viewer, "(%d,%d) - [sub mat]: (PETSC_NULL), [parent mat]: %p \n\n", pi,pj, parent_A ); 
		else PetscViewerASCIIPrintf( viewer, "(%d,%d) - [sub mat]: (PETSC_NULL), [parent mat]: %s \n\n", pi,pj, pname ); 
		
		PetscViewerASCIIPopTab( viewer );
		return;
	}
	
	
	MatGetType( A, &type );
	is_block = PETSC_FALSE;
	PetscTypeCompare( (PetscObject)A, "block", &is_block );
	
	if( is_block  == PETSC_FALSE ) {
		MatGetSize( A, &m, &n );
		
		
		PetscViewerASCIIPushTab( viewer );
		
		//PetscViewerASCIIPrintf( viewer, "(i,j) - [sub mat]: (Name / ptr) , [parent mat]: (Name / ptr) \n" );
		PetscViewerASCIIPrintf( viewer, "-------------------- Begin ( %p ) --------------------\n", A );
		if( name == PETSC_NULL ) { 
			if( pname == PETSC_NULL ) PetscViewerASCIIPrintf( viewer, "(%d,%d) - [sub mat]: (%p), [parent mat]: %p \n", pi,pj, A, parent_A ); 
			else PetscViewerASCIIPrintf( viewer, "(%d,%d) - [sub mat]: (%p), [parent mat]: %s \n", pi,pj, A, pname ); 
		}
		else { 
			if( pname == PETSC_NULL ) PetscViewerASCIIPrintf( viewer, "(%d,%d) - [sub mat]: (%s), [parent mat]: %p \n", pi,pj, name, parent_A ); 
			else PetscViewerASCIIPrintf( viewer, "(%d,%d) - [sub mat]: (%s), [parent mat]: %s \n", pi,pj, name, pname ); 
		}
		MatView( A, viewer);
		PetscViewerASCIIPrintf( viewer, "--------------------   End ( %p ) --------------------\n\n", A );
		
		PetscViewerASCIIPopTab( viewer );
		
		return;
	}
	else {
		MatGetSize( A, &m, &n );
		
		bA = (Mat_Block*)A->data;
		for( i=0; i<bA->nr; i++ ) {
			for( j=0; j<bA->nc; j++ ) {
				print_mat_contents( bA->m[i][j],A, i,j, viewer );
			}
		}
	}
}


PetscErrorCode MatView_Block( Mat A, PetscViewer viewer )
{
        Mat_Block *block = (Mat_Block*)A->data;
	PetscTruth isascii;
	
	PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&isascii);
	if (isascii) {

		PetscViewerASCIIPushTab( viewer );		// push0
	
		PetscViewerASCIIPrintf(viewer,"MatBlock structure: \n" );
		PetscViewerASCIIPrintf( viewer, "(i,j) - (Name / ptr), type=, rows=, cols=\n" );
		print_mat_structure( A, 0,0, viewer );

		PetscViewerASCIIPrintf(viewer,"Block structure: %s \n", MatBlockStructureName[block->bstruct] );

	
		PetscViewerASCIIPrintf(viewer,"MatBlock contents: \n" );
		PetscViewerASCIIPrintf( viewer, "(i,j) - [sub mat]: (Name / ptr) , [parent mat]: (Name / ptr) \n" );
		PetscViewerASCIIPushTab( viewer );		// push1
		print_mat_contents( A, PETSC_NULL, 0,0, viewer );
		PetscViewerASCIIPopTab( viewer );		// pop1
	
		PetscViewerASCIIPopTab( viewer );		// pop0 
	}
	PetscFunctionReturn(0);
}

/*


*a = [ [a1] [0]  [0] ]
     [  [0] [a2] [0] ]
     [  [0] [0] [a3] ]

PetscErrorCode  MatGetDiagonalBlock_Block( Mat A, PetscTruth *iscopy, MatReuse reuse, Mat *a )
{
	Mat_Block *bA = (Mat_Block*)A->data;
	Mat ba;
	int i;
	
	PetscFunctionBegin;
	
	for( i=0; i<bA->nr; i++ ) {
		ba = bA->m[i][i];
		MatGetDiagonalBlock( ba, iscopy, reuse, a );
	}
	
	PetscFunctionReturn(0);
}


MatGetSubMatrices( Mat mat,PetscInt n,const IS irow[],const IS icol[],MatReuse scall,Mat *submat[] )

*/



extern PetscErrorCode MatSetUp_Block( Mat A );

PetscErrorCode MatMatMult_Block( Mat A, Mat B, MatReuse scall, PetscReal max_fill, Mat *_C )
{
	Mat C;
	Mat_Block *bA = (Mat_Block*)A->data;
	Mat_Block *bB = (Mat_Block*)B->data;
	Mat_Block *bC;
	PetscInt i,j,k;
	PetscInt a_nr,a_nc,b_nr,b_nc;
	MPI_Comm comm;
	Mat Aik,Bkj, tmp;
	char *prefix;
	
	_check_mat_mat_compatibility( A, B );
	
	if( scall == MAT_INITIAL_MATRIX ) {
		/* create a new matrix */
		MatGetSize( A, &a_nr, &a_nc );
		MatGetSize( B, &b_nr, &b_nc );
		
		PetscObjectGetComm( (PetscObject)A, &comm );
		MatCreate( comm, &C );
		MatSetSizes( C, a_nr,b_nc, a_nr,b_nc );
		MatSetType( C, "block" );
		MatSetUp_Block( C );
	}
	else {
		C = *_C;
	}
	bC = (Mat_Block*)C->data;
	
	
	/* init any blocks */
	for( i=0; i<bC->nr; i++ ) {
		for( j=0; j<bC->nc; j++ ) {
			if( bC->m[i][j] != PETSC_NULL ) {
				MatZeroEntries( bC->m[i][j] );
			}
		}
	}
	
	
	
	/* Do multiplication C_ij = A_ik B_kj */
	for( i=0; i<bA->nr; i++ ) {
		for( j=0; j<bB->nc; j++ ) {
          //Cij = bC->m[i][j];
			
			for( k=0; k<bA->nc; k++ ) {
				Aik = bA->m[i][k];
				Bkj = bB->m[k][j];
				
				if( Aik == PETSC_NULL || Bkj == PETSC_NULL ) {
					continue;
				}
				
				if( bC->m[i][j] == PETSC_NULL ) { /* no matrix exist in ij */
				//	printf("C[%d][%d] = A[%d][%d] * B[%d][%d] \n", i,j, i,k, k,j );
					MatMatMult( Aik, Bkj, MAT_INITIAL_MATRIX, max_fill, &tmp );
					asprintf( &prefix, "matbk_%d%d", i,j ); 
					MatSetOptionsPrefix( tmp, prefix );
					free( prefix );
					
					bC->m[i][j] = tmp;
				}
				else {
				//	printf("C[%d][%d] = C[%d][%d] + A[%d][%d] * B[%d][%d] \n", i,j,i,j, i,k, k,j );
					MatMatMult( Aik, Bkj, MAT_INITIAL_MATRIX, max_fill, &tmp );
					MatAXPY( bC->m[i][j], 1.0, tmp, DIFFERENT_NONZERO_PATTERN );
					Stg_MatDestroy( &tmp );
				}
				
			}
			
		}
	}

        /* update block structure */
        MatBlockGetMatBlockStructure_Block( C, &bC->bstruct );
	
	MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY);

	*_C = C;

	
	PetscFunctionReturn(0);
}

PetscErrorCode MatTranspose_Block( Mat A, MatReuse reuse, Mat *_C ) 
{
	Mat_Block *bA = (Mat_Block*)A->data;
	Mat_Block *bAt;
	Mat At;
	MPI_Comm comm;
	PetscInt i,j;
	
	PetscFunctionBegin;
	
	
	PetscObjectGetComm( (PetscObject)A, &comm );
	MatCreate( comm, &At );
	MatSetSizes( At, bA->nc,bA->nr, bA->nc,bA->nr );
	MatSetType( At, "block" );
	MatSetUp_Block( At );
	bAt = (Mat_Block*)At->data;
	
	
	for( i=0; i<bA->nr; i++ ) {
		for( j=0; j<bA->nc; j++ ) {
			if( bA->m[i][j] == PETSC_NULL ) continue;
			
			MatTranspose( bA->m[i][j], reuse, &bAt->m[j][i] );
		}
	}

        /* update block structure */
        MatBlockGetMatBlockStructure_Block( At, &bAt->bstruct );
	
	*_C = At;

	PetscFunctionReturn(0);
}



PetscErrorCode MatPtAP_Block( Mat A, Mat P, MatReuse scall, PetscReal fill, Mat *_C ) 
{
	Mat C;
	Mat_Block *bA = (Mat_Block*)A->data;
	Mat_Block *bP = (Mat_Block*)P->data;
	Mat_Block *bC, *bPt;
	PetscInt i,j;
	MPI_Comm comm;
	PetscInt nd_cnt, nod_cnt;
	Mat Pt, AijPj;
	
	
	PetscFunctionBegin; 
	
	_check_mat_mat_compatibility( A, P );
	
	/* check P is diagonal */
	nd_cnt = nod_cnt = 0;
	for( i=0; i<bP->nr; i++ ) {
		for( j=0; j<bP->nc; j++ ) {
			if( bP->m[i][j] == PETSC_NULL ) continue;
			
			if( i == j ) nd_cnt++;
			else nod_cnt++;
		}
	}
	
	if( nod_cnt != 0 ) {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "Operation MatPtAP_Block only valid if P is diagonal.");
	}
	if( bP->nr != bP->nc ) {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "Operation MatPtAP_Block only valid if P is square.");
	}
	
	
	
	
	PetscObjectGetComm( (PetscObject)A, &comm );
	MatCreate( comm, &C );
	MatSetSizes( C, bP->nr,bP->nc, bP->nr,bP->nc );
	MatSetType( C, "block" );
	MatSetUp_Block( C );
	bC = (Mat_Block*)C->data;
	
	
	/* Do diagonals, Pt_i A_i P_i */
	for( i=0; i<bP->nr; i++ ) {
		if( bA->m[i][i] == PETSC_NULL || bP->m[i][i] == PETSC_NULL ) continue;
		
		MatPtAP( bA->m[i][i], bP->m[i][i], MAT_INITIAL_MATRIX, fill, &bC->m[i][i] );
	}
	
	/* Pt_i A_ij P_j */
	MatTranspose( P, scall, &Pt );
	bPt = (Mat_Block*)Pt->data;
	for( i=0; i<bA->nr; i++ ) {
		for( j=0; j<bA->nc; j++ ) {
			
			if( i==j ) continue; /* skip diagonal */
			if( bA->m[i][j] == PETSC_NULL || bP->m[i][i] == PETSC_NULL ) continue;
			
			MatMatMult( bA->m[i][j], bP->m[j][j], MAT_INITIAL_MATRIX, fill, &AijPj );
			MatMatMult( bPt->m[i][i], AijPj, MAT_INITIAL_MATRIX, fill, &bC->m[i][j] );
			Stg_MatDestroy( &AijPj );
			
		}
	}
	Stg_MatDestroy( &Pt );
	
        /* update block structure */
        MatBlockGetMatBlockStructure_Block( C, &bC->bstruct );

	*_C = C;
	
	PetscFunctionReturn(0);
}

PetscErrorCode MatZeroEntries_Block( Mat A )
{
	Mat_Block *bA = (Mat_Block*)A->data;
	PetscInt i,j;
	
	
	PetscFunctionBegin; 
	
	for( i=0; i<bA->nr; i++ ) {
		for( j=0; j<bA->nc; j++ ) {
			if( bA->m[i][j] == PETSC_NULL ) continue;
			MatZeroEntries( bA->m[i][j] );
	}}
	
	PetscFunctionReturn(0);
}



PetscErrorCode MatDuplicate_Block( Mat A, MatDuplicateOption op, Mat *B )
{
	Mat_Block *bA = (Mat_Block*)A->data;
	Mat_Block *bB;
	PetscInt i,j;
	Mat subB;
	MPI_Comm comm;
	
	
	PetscFunctionBegin; 
	
	PetscObjectGetComm( (PetscObject)A, &comm );
	MatCreate( comm, B );
	MatSetSizes( *B, bA->nr,bA->nc, bA->nr,bA->nc );
	MatSetType( *B, "block" );
	
	bB = (Mat_Block*)( (*B)->data );
	for( i=0; i<bA->nr; i++ ) {
		for( j=0; j<bA->nc; j++ ) {
			if( bA->m[i][j] == PETSC_NULL ) continue;
			MatDuplicate( bA->m[i][j], op, &subB );
			MatBlockSetValue( *B, i,j, subB, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
			/* Return control of subB back to B */
			Stg_MatDestroy( &subB );
	}}

        /* update block structure */
        MatBlockGetMatBlockStructure_Block( *B, &bB->bstruct );
	
	MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY);

	
	PetscFunctionReturn(0);
}

PetscErrorCode MatConvert_Block(Mat Ab, const MatType mtype,MatReuse reuse,Mat *A)
{

	if( reuse == MAT_INITIAL_MATRIX ) {
		*A = PETSC_NULL;
	}

	MatBlockMergeSubBlocks( Ab, INSERT_VALUES, mtype, A );

	PetscFunctionReturn(0);
}


