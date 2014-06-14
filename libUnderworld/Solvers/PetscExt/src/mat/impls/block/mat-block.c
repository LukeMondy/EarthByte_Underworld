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


/* Implements a basic block matrix implementation. */

#include <stdlib.h>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <private/vecimpl.h>
#include <private/matimpl.h>

#include "src/mat/impls/block/mat_block_impl.h"
#include "src/mat/impls/block/mat-block-ops.h"
#include "private/mat/petscmat-block.h"


const char *MatBlockStructureName[] = { "NULL_BLOCK", "DEFAULT_BLOCK", "DIAGONAL_BLOCK", "UPPER_BLOCK", "LOWER_BLOCK" };

#undef __FUNCT__  
#define __FUNCT__ "MatBlockGetMatBlockStructure_Block"
PetscErrorCode MatBlockGetMatBlockStructure_Block( Mat A, MatBlockStructure *bstruct )
{
  Mat_Block *ctx = (Mat_Block*)A->data;
  PetscInt i,j;
  PetscTruth has_diag, has_upper, has_lower;
  PetscInt M,N;
        
  PetscFunctionBegin;

  M = A->rmap->N;
  N = A->cmap->N;

  *bstruct = DEFAULT_BLOCK;
  
  /* if any are not NULL, then we can assume it may be a DEFAULT_BLOCK */
  for( i=0; i<M; i++ ) {
    for( j=0; j<N; j++ ) {
      if( ctx->m[i][j] != PETSC_NULL ) {
        *bstruct = DEFAULT_BLOCK;
      }
    }
  }


  /* Look for special cases */
  has_diag = has_upper =has_lower = PETSC_FALSE;
  for( i=0; i<M; i++ ) {
    for( j=0; j<N; j++ ) {

      if( ctx->m[i][j] == PETSC_NULL ) { continue; }

       /* diagonal */
      if( i == j ) {  has_diag  = PETSC_TRUE;  }
      if( j > i ) {   has_upper = PETSC_TRUE;  }
      if( i < j ) {   has_lower = PETSC_TRUE;  }
  
    }
  }

  if( has_diag == PETSC_TRUE ) {   *bstruct = DIAGONAL_BLOCK;  }
  if( (has_diag==PETSC_TRUE) && (has_upper==PETSC_TRUE) ) {
    *bstruct = UPPER_BLOCK;
  }
  if( (has_diag==PETSC_TRUE) && (has_lower==PETSC_TRUE) ) {
    *bstruct = LOWER_BLOCK;
  }
  if( (has_lower==PETSC_TRUE) && (has_upper==PETSC_TRUE) ) {
    *bstruct = DEFAULT_BLOCK;
  }


  PetscFunctionReturn(0);
}




#undef __FUNCT__  
#define __FUNCT__ "MatSetUp_Block"
PetscErrorCode MatSetUp_Block( Mat A )
{
	Mat_Block *ctx = (Mat_Block*)A->data;
	PetscInt i,j;
	
	
	PetscFunctionBegin;
	
	if( ctx->setup_called == PETSC_TRUE ) PetscFunctionReturn(0);
	
	
	ctx->nr = A->rmap->N;
	ctx->nc = A->cmap->N;
	
	if( ctx->nr < 0 ) {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "Cannot create MAT_BLOCK with < 0 row blocks." );
	}
	if( ctx->nc < 0 ) {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "Cannot create MAT_BLOCK with < 0 col blocks." );
	}
	
	/* Create space */
	ctx->m = (Mat**)malloc( sizeof(Mat*) * ctx->nr );
	for( i=0; i<ctx->nr; i++ ) {
		ctx->m[i] = (Mat*)malloc( sizeof(Mat) * ctx->nc );
	}
	
	for( i=0; i<ctx->nr; i++ ) {
		for( j=0; j<ctx->nc; j++ ) {
			ctx->m[i][j]       = PETSC_NULL;
	}}
	
	ctx->setup_called = PETSC_TRUE;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatSetSizes_Block"
PetscErrorCode  MatSetSizes_Block( Mat A, PetscInt m, PetscInt n, PetscInt M, PetscInt N )
{
	Mat_Block *ctx = (Mat_Block*)A->data;
	
	
	PetscFunctionBegin;
	
	if( M == PETSC_DETERMINE ) {
		Stg_SETERRQ( PETSC_ERR_ARG_INCOMP, "MatBlock: Must specify global row size of matrix." );
	}
	if( N == PETSC_DETERMINE ) {
		Stg_SETERRQ( PETSC_ERR_ARG_INCOMP, "MatBlock: Must specify global col size of matrix." );
	}
	
	ctx->nr = M;
	ctx->nc = N;
	
	MatSetUp_Block( A );
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatDestroy_Block"
PetscErrorCode MatDestroy_Block( Mat A )
{
	Mat_Block        *vs = (Mat_Block*)A->data;
	PetscErrorCode ierr;
	PetscInt i,j;
	
	
	PetscFunctionBegin;
	
	/* if memory was published with AMS then destroy it */
	ierr = PetscObjectDepublish(A);CHKERRQ(ierr);
	
	/* release the matrices and the place holders */
	if( vs->m != PETSC_NULL ) {
		for( i=0; i<vs->nr; i++ ) {
			for( j=0; j<vs->nc; j++ ) {
				
				if( vs->m[i][j] != PETSC_NULL ) {
                    Stg_MatDestroy( &(vs->m[i][j]) );
					vs->m[i][j] = PETSC_NULL;
				}
			}
			free( vs->m[i] );
			vs->m[i] = PETSC_NULL;
		}
		free( vs->m );
		vs->m = PETSC_NULL;
	}
	
	ierr=PetscFree( vs );CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}




/*
This function is called every time we insert a sub matrix the block.
It ensures that every block along row r, and col c has the same dimensions
as the submat being inserted.
*/
#undef __FUNCT__  
#define __FUNCT__ "MatBlock_CheckConsistency"
PetscErrorCode MatBlock_CheckConsistency( Mat A, Mat submat, PetscInt r, PetscInt c )
{
	Mat_Block *b = (Mat_Block*)A->data;
	PetscInt i,j;
	PetscInt nr,nc;
	PetscInt sM, sN, M,N;
	Mat mat;
	
	nr = b->nr;
	nc = b->nc;
	PetscFunctionBegin;
	
	
	MatGetSize( submat, &sM, &sN );
	for( i=0; i<nr; i++ ) {
		mat = b->m[i][c];
		if( mat != PETSC_NULL ) {
			MatGetSize( mat, &M, &N );
			/* Check columns match */
			
			if( sN != N ) {
				Stg_SETERRQ3(PETSC_ERR_SUP,"Inserted incompatible submatrix into block at (%D,%D). Submatrix must have %D rows",r,c,N );
			}
			
		}
	}
	
	for( j=0; j<nc; j++ ) {
		mat = b->m[r][j];
		if( mat != PETSC_NULL ) {
			MatGetSize( mat, &M, &N );
			/* Check rows match */
			
			if( sM != M ) {
				Stg_SETERRQ3(PETSC_ERR_SUP,"Inserted incompatible submatrix into block at (%D,%D). Submatrix must have %D cold",r,c,M );
			}
			
		}
	}
	
	
	PetscFunctionReturn(0);
}

/* ================================== Public functions ===================================== */

EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "MatBlockSetValues_Block"
PetscErrorCode MatBlockSetValues_Block( 
		Mat A, 
		PetscInt m, const PetscInt idxm[],
		PetscInt n ,const PetscInt idxn[],
		const Mat mat[], MatStructure str, InsertMode addv )
{
	Mat_Block *b = (Mat_Block*)A->data;
	PetscInt i,j;
	PetscInt row,col;
	
	
	PetscFunctionBegin;
	
	if (!m || !n) return(0);
	
	MatSetUp_Block( A );
	
	for( i=0; i<m; i++ ) {
		row = idxm[i];
		if( row >= A->rmap->N ) Stg_SETERRQ2(PETSC_ERR_ARG_OUTOFRANGE,"Row too large: row %D max %D",row,A->rmap->N-1);
		
		for( j=0; j<n; j++ ) {
			col = idxn[j];
			if( col >= A->cmap->N ) Stg_SETERRQ2(PETSC_ERR_ARG_OUTOFRANGE,"Col too large: col %D max %D",col,A->rmap->N-1);
			
			
			MatBlock_CheckConsistency( A, mat[j+i*n], row,col );
			
			if( b->m[row][col] == PETSC_NULL ) {
				PetscObjectReference( (PetscObject)mat[ j + i*n ] );
				b->m[row][col] = mat[ j + i*n ];
				
				continue;
			}
			else {
				
				if( addv == INSERT_VALUES ) {
                    Stg_MatDestroy( &(b->m[row][col]) );
					
					PetscObjectReference( (PetscObject)mat[ j + i*n ] );
					b->m[row][col] = mat[ j + i*n ]; // if( row_oriented )
				}
				else {
					MatAXPY( b->m[row][col], 1.0, mat[ j + i*n ], str );
				}
				
			}
		}
	}

	/* update block structure */
	MatBlockGetMatBlockStructure_Block( A, &b->bstruct );
	
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "MatBlockSetValues"
PetscErrorCode MatBlockSetValues( 
		Mat A, 
		PetscInt m, const PetscInt idxm[],
		PetscInt n ,const PetscInt idxn[],
		const Mat mat[], MatStructure str, InsertMode addv )
{
	PetscErrorCode ierr,(*f)(Mat,PetscInt,const PetscInt*,PetscInt,const PetscInt*,const Mat*,MatStructure,InsertMode);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)A,"MatBlockSetValues_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(A,m,idxm,n,idxn,mat,str,addv);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}


/* =========================================================================== */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "MatBlockSetValue_Block"
PetscErrorCode MatBlockSetValue_Block( Mat A, PetscInt idxm, PetscInt idxn, Mat mat, MatStructure str, InsertMode addv )
{
	Mat_Block *b = (Mat_Block*)A->data;
	PetscInt row,col;
	
	
	PetscFunctionBegin;
	
	MatSetUp_Block( A );
	
	row = idxm;
	col = idxn;
	if( row >= A->rmap->N ) Stg_SETERRQ2(PETSC_ERR_ARG_OUTOFRANGE,"Row too large: row %D max %D",row,A->rmap->N-1);
	if( col >= A->cmap->N ) Stg_SETERRQ2(PETSC_ERR_ARG_OUTOFRANGE,"Col too large: col %D max %D",col,A->rmap->N-1);
	
	MatBlock_CheckConsistency( A, mat, row,col );
	
	if( b->m[row][col] == PETSC_NULL ) {
		PetscObjectReference( (PetscObject)mat );
		b->m[row][col] = mat;
		
	        /* update block structure */
        	MatBlockGetMatBlockStructure_Block( A, &b->bstruct );

		PetscFunctionReturn(0);
	}
	else {
		if( addv == INSERT_VALUES ) {
            Stg_MatDestroy( &(b->m[row][col]) );
			
			PetscObjectReference( (PetscObject)mat );
			b->m[row][col] = mat; // if( row_oriented )
		}
		else {
			MatAXPY( b->m[row][col], 1.0, mat, str );
		}
	}

        /* update block structure */
        MatBlockGetMatBlockStructure_Block( A, &b->bstruct );
	
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "MatBlockSetValue"
PetscErrorCode MatBlockSetValue( Mat A, PetscInt idxm, PetscInt idxn, Mat mat, MatStructure str, InsertMode addv )
{
	PetscErrorCode ierr,(*f)(Mat,PetscInt,PetscInt,Mat,MatStructure,InsertMode);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)A,"MatBlockSetValue_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(A,idxm,idxn,mat,str,addv);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}



/* =========================================================================== */
#undef __FUNCT__  
#define __FUNCT__ "_MatBlockGetValue"
PetscErrorCode _MatBlockGetValue( Mat A, PetscInt idxm, PetscInt idxn, Mat *mat )
{
	Mat_Block *b = (Mat_Block*)A->data;
	PetscInt row,col;
	
	
	PetscFunctionBegin;
	
	row = idxm;
	col = idxn;
	if( row >= A->rmap->N ) Stg_SETERRQ2(PETSC_ERR_ARG_OUTOFRANGE,"Row too large: row %D max %D",row,A->rmap->N-1);
	if( col >= A->cmap->N ) Stg_SETERRQ2(PETSC_ERR_ARG_OUTOFRANGE,"Col too large: col %D max %D",col,A->rmap->N-1);
	
	
	*mat = b->m[ row ][ col ];
	
	PetscFunctionReturn(0);
}


/*
The functions MatBlockGetValues() and MatBlockGetValue() should possibly be renamed
to MatBlockGetSubMatrices() to be more consistent with VecGetArray() / VecRestoreArray()

If either
        i) MatBlockGetSubMatrix()
        ii) MatBlockGetSubMatrices()
are called, you MUST call
        MatBlockRestoreSubMatrices()
otherwise the block matrices state will not be increased, even through its 
submatrices (subordinates) have potentially been modified and have had their
state increased.
*/


/* =========================================================================== */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "MatBlockGetSubMatrix_Block"
PetscErrorCode MatBlockGetSubMatrix_Block( Mat A, PetscInt ridx, PetscInt cidx, Mat *sa )
{
	PetscFunctionBegin;
	_MatBlockGetValue( A, ridx, cidx, sa );
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "MatBlockGetSubMatrix"
PetscErrorCode MatBlockGetSubMatrix( Mat A, PetscInt ridx, PetscInt cidx, Mat *sa )
{
	PetscErrorCode ierr,(*f)(Mat,PetscInt,PetscInt,Mat*);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)A,"MatBlockGetSubMatrix_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(A,ridx,cidx,sa);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}

/* =========================================================================== */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "MatBlockGetSubMatrices_Block"
PetscErrorCode MatBlockGetSubMatrices_Block( Mat A, Mat ***sa )
{
	Mat_Block *b = (Mat_Block*)A->data;
	
	
	PetscFunctionBegin;
	*sa = b->m;
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "MatBlockGetSubMatrices"
PetscErrorCode MatBlockGetSubMatrices( Mat A, Mat ***sa )
{
	PetscErrorCode ierr,(*f)(Mat,Mat***);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)A,"MatBlockGetSubMatrices_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(A,sa);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}

/* =========================================================================== */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "MatBlockRestoreSubMatrices_Block"
PetscErrorCode MatBlockRestoreSubMatrices_Block( Mat A )
{
	PetscFunctionBegin;
	PetscObjectStateIncrease((PetscObject)A);
	PetscFunctionReturn(0);
}
EXTERN_C_END
		
#undef __FUNCT__  
#define __FUNCT__ "MatBlockRestoreSubMatrices"
PetscErrorCode MatBlockRestoreSubMatrices( Mat A )
{
	PetscErrorCode ierr,(*f)(Mat);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)A,"MatBlockRestoreSubMatrices_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(A);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}


/* =========================================================================== */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "MatBlockCreateSubMatrix_Block"
PetscErrorCode MatBlockCreateSubMatrix_Block( Mat A, const MatType mtype, PetscInt r, PetscInt c, Mat *sA )
{
	Mat_Block *bA = (Mat_Block*)A->data;
	PetscInt i,j;
	PetscInt m,n, M,N;
	Mat subA,s;
	MPI_Comm comm;
	
	
	PetscFunctionBegin;
	
	if( bA->m[r][c] != PETSC_NULL ) {
		Stg_SETERRQ2( PETSC_ERR_ARG_WRONG, "Cannot create new sub matrix at (%D,%D) as object already present.", r,c );
	}
	
	m = M = n = N = -1;
	for( i=0; i<bA->nr; i++ ) {
		if( bA->m[i][c] == PETSC_NULL ) continue;
		
		MatBlockGetSubMatrix( A, i,c, &s );
		MatGetSize( s, PETSC_NULL, &N );
		MatGetLocalSize( s, PETSC_NULL, &n );
	}
	if( n==-1 || N==-1 ) {
		Stg_SETERRQ1( PETSC_ERR_ARG_WRONG, "MAT BLOCK contains a null column at %D.", c );
	}
	
	for( j=0; j<bA->nc; j++ ) {
		if( bA->m[r][j] == PETSC_NULL ) continue;
		
		MatBlockGetSubMatrix( A, r,j, &s );
		MatGetSize( s, &M, PETSC_NULL );
		MatGetLocalSize( s, &m, PETSC_NULL );
	}
	if( m==-1 || M==-1 ) {
		Stg_SETERRQ1( PETSC_ERR_ARG_WRONG, "MAT BLOCK contains a null row at %D.", r );
	}
	
	MatBlockRestoreSubMatrices( A );
	
	
	PetscObjectGetComm( (PetscObject)A, &comm );
	MatCreate( comm, &subA );
	MatSetSizes( subA, m,n, M,N );
	if(mtype!=PETSC_NULL) {
		MatSetType( subA, mtype );
	}
	MatSetFromOptions(subA);	

	MatBlockSetValue( A, r,c, subA, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	
	*sA = subA;
	
#if 0
	/* Return control of releasing subA to the block matrix A */
	Stg_MatDestroy( &subA );
#endif
	
	
	PetscFunctionReturn(0);
}
EXTERN_C_END

/*
Given the block matrix A, we generates a new submatrix *sA at location (r,c).
If (r,c) has either a null row or column then this will fail. If you wish to 
return control of the sub matrix sA to the block matrix, call Stg_MatDestroy(sA) 
just after call MatBlockCreateSubMatrix(). This routine WILL automatically set
the submatrix sA into the location (r,c).
*/
#undef __FUNCT__  
#define __FUNCT__ "MatBlockCreateSubMatrix"
PetscErrorCode MatBlockCreateSubMatrix( Mat A, const MatType mtype, PetscInt r, PetscInt c, Mat *sA )
{
	PetscErrorCode ierr,(*f)(Mat,const MatType,PetscInt,PetscInt,Mat*);
	
        PetscFunctionBegin;
        ierr = PetscObjectQueryFunction((PetscObject)A,"MatBlockCreateSubMatrix_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(A,mtype,r,c,sA);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}

/* =========================================================================== */

/*
Create a matrix which we can merge a block matix into
*/

EXTERN_C_BEGIN

#undef __FUNCT__
#define __FUNCT__ "MatCreateFromMatBlock_Block"
PetscErrorCode MatCreateFromMatBlock_Block( Mat bA, const MatType type,Mat *mA )
{
  Mat_Block *block = (Mat_Block*)bA->data;
  PetscInt i,j, AM_sum,AN_sum, M,N;
  //MPI_Comm comm;

  PetscFunctionBegin;

  /* Determine size of new matrix */      
  AM_sum = AN_sum = 0;
        
        
  /* rows */      
  for( i=0; i<block->nr; i++ ) {
    for( j=0; j<block->nc; j++ ) {
      if( block->m[i][j] != PETSC_NULL ) {
        MatGetSize( block->m[i][j], &M, &N );
        AM_sum = AM_sum + M;
        break;
      }
    }
    if( j==(block->nc) ) {
      Stg_SETERRQ(PETSC_ERR_SUP, "Found an empty row" );
    }
  }
   
  /* cols */
  for( j=0; j<block->nc; j++ ) {
    for( i=0; i<block->nr; i++ ) {
      if( block->m[i][j] != PETSC_NULL ) {
        MatGetSize( block->m[i][j], &M, &N );
        AN_sum = AN_sum + N;
        break;
      }
    }
    if( i==(block->nr) ) {
      Stg_SETERRQ(PETSC_ERR_SUP, "Found an empty column" );
    }
  }


                
  MatCreate( PETSC_COMM_WORLD, mA );
  MatSetSizes( *mA, PETSC_DECIDE,PETSC_DECIDE, AM_sum,AN_sum );
  if( type!=PETSC_NULL) {
    MatSetType( *mA, type );
  }
  MatSetFromOptions( *mA );


  PetscFunctionReturn(0);
}
EXTERN_C_END



/* =========================================================================== */

#undef __FUNCT__  
#define __FUNCT__ "MatCreateFromMatBlock"
PetscErrorCode MatCreateFromMatBlock(  Mat Ab, const MatType mtype, Mat *A  )
{
        PetscErrorCode ierr,(*f)(Mat,const MatType,Mat*);
        
        PetscFunctionBegin;
        ierr = PetscObjectQueryFunction((PetscObject)Ab,"MatCreateFromMatBlock_C",(void (**)(void))&f);CHKERRQ(ierr);
        if (f) {
                ierr = (*f)(Ab,mtype,A);CHKERRQ(ierr);
        }
        PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MergeSubBlocksPreallocateMatAIJ"
PetscErrorCode MergeSubBlocksPreallocateMatAIJ( Mat bA, Mat A )
{
  Mat_Block *block = (Mat_Block*)bA->data;
  PetscTruth is_seqaij, is_mpiaij;
  PetscInt *dnnz;
  PetscInt I,J;
  PetscLogDouble t0,t1;
  PetscTruth has_info, flg;


  PetscTypeCompare( (PetscObject)A, MATSEQAIJ, &is_seqaij );
  PetscTypeCompare( (PetscObject)A, MATSEQAIJ, &is_mpiaij );

  /* check type is aij */
  if( (is_seqaij==PETSC_FALSE) && (is_mpiaij==PETSC_FALSE) ) {
    PetscFunctionReturn(0);
  }
  /* check we need to allocate */
  if( A->num_ass != 0 ) {
  /*  preallocated = PETSC_TRUE; */
    PetscFunctionReturn(0);
  }


  has_info = PETSC_FALSE;
  PetscOptionsGetTruth( PETSC_NULL, "-info", &has_info, &flg );



  PetscGetTime(&t0);
  if( is_seqaij==PETSC_TRUE ) {
    PetscInt MA,Mm,Nm;
    PetscInt row,row_offset,si,start_row,end_row,ncols;

    MatGetSize( A, &MA, PETSC_NULL );
    PetscMalloc( sizeof(PetscInt)*MA, &dnnz );
    for( si=0; si<MA; si++ ) {  dnnz[si] = 0;  }


    row_offset = 0;
    for( I=0; I<block->nr; I++ ) {
      for( J=0; J<block->nc; J++ ) {

        if( block->m[I][J]==PETSC_NULL ) {  continue;  }

        MatGetSize( block->m[I][J], &Mm, &Nm );

        MatGetOwnershipRange( block->m[I][J], &start_row, &end_row );
        for( si=start_row; si<end_row; si++ ) {
          MatGetRow( block->m[I][J], si, &ncols, PETSC_NULL, PETSC_NULL );

          row = si + row_offset;
	  dnnz[ row ] = dnnz[ row ] + ncols;
        }
        
      }
      row_offset = row_offset + Mm;
    }

    MatSeqAIJSetPreallocation( A, PETSC_NULL, dnnz );
  
    PetscFree( dnnz );
  }
  else { /* Must be MATMPIAIJ */
    PetscPrintf( PETSC_COMM_WORLD, "  WARNING: MatMergeSubBlocks is not performing preallocation for MATMPIAIJ \n");





  }

  PetscGetTime(&t1);

  if(has_info==PETSC_TRUE) {
    PetscPrintf( PETSC_COMM_WORLD, "  MatMergeSubBlocks: Preallocation (%lf sec) \n", t1-t0 );
  }


  PetscFunctionReturn(0);
}





/*
Merge all non NULL sub blocks into a single matrix
*/
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "MatBlockMergeSubBlocks_Block"
PetscErrorCode MatBlockMergeSubBlocks_Block( Mat Ab, InsertMode addv, const MatType mtype, Mat *A )
{
	Mat_Block *block = (Mat_Block*)Ab->data;
	PetscInt i,j;
	PetscInt m,n, M,N;
	PetscInt *AM, *AN;
	PetscInt AM_sum, AN_sum;
	MPI_Comm comm;
	PetscInt si,sj,start_row,end_row,row;
	PetscInt ncols, *shifted_cols;
	const PetscInt *cols;
	const PetscScalar *vals;
	/* for preallocation */
	const PetscInt *ranges_A;
	PetscInt rank, nprocs;
	PetscTruth is_MATAIJ;
        const MatType _mtype;	
	PetscLogDouble tt0,tt1,t0,t1,dt0,dt1,time_mgr,time_msv;
	PetscTruth has_info,flg;
	
	
	PetscFunctionBegin;
	

	has_info = PETSC_FALSE;
	PetscOptionsGetTruth( PETSC_NULL, "-info", &has_info, &flg );



	PetscGetTime(&tt0);
	PetscMalloc( sizeof(PetscInt)*block->nr, &AM );
	PetscMalloc( sizeof(PetscInt)*block->nc, &AN );
	
	/* Determine size of new matrix */
	AM_sum = AN_sum = 0;
	
	/* rows */
	for( i=0; i<block->nr; i++ ) {
		for( j=0; j<block->nc; j++ ) {
			if( block->m[i][j] != PETSC_NULL ) {
				MatGetSize( block->m[i][j], &M, &N );
				AM[i] = AM_sum;
				AM_sum = AM_sum + M;
				break;
			}
		}
		if( j==(block->nc) ) {
			Stg_SETERRQ(PETSC_ERR_SUP, "Found an empty row" );
		}
	}
	/* cols */
	for( j=0; j<block->nc; j++ ) {
		for( i=0; i<block->nr; i++ ) {
			if( block->m[i][j] != PETSC_NULL ) {
				MatGetSize( block->m[i][j], &M, &N );
				AN[j] = AN_sum;
				AN_sum = AN_sum + N;
				break;
			}
		}
		if( i==(block->nr) ) {
			Stg_SETERRQ(PETSC_ERR_SUP, "Found an empty column" );
		}
	}
	
	
	
	PetscObjectGetComm( (PetscObject)Ab, &comm );

	if( *A == PETSC_NULL ) {
		MatCreate( comm, A );
		MatSetSizes( *A, PETSC_DECIDE,PETSC_DECIDE, AM_sum,AN_sum );
		if( mtype!=PETSC_NULL) {
			MatSetType( *A, mtype );
		}
		MatSetFromOptions( *A );
	}
	else {
		PetscInt A_M, A_N;
		PetscTruth is_same;

		/* check type */
		is_same = PETSC_FALSE;

		/* need a special test for MATAIJ. Convert it to SEQAIJ or MPIAIJ */
		is_MATAIJ = PETSC_FALSE;
		PetscStrcmp(mtype,MATAIJ,&is_MATAIJ);
		if(is_MATAIJ==PETSC_TRUE) {
			PetscMPIInt size;

			MPI_Comm_size(comm,&size);
			_mtype = MATSEQAIJ;
			if(size>1) {  _mtype = MATMPIAIJ;  }
		}
		else {
			_mtype = mtype;
		}

		PetscTypeCompare( (PetscObject)(*A), _mtype, &is_same );
		if( is_same==PETSC_FALSE ) {
			Stg_SETERRQ( PETSC_ERR_ARG_NOTSAMETYPE, "Re-used matrix was wrong type" );
		}

		/* check size */
		MatGetSize( *A, &A_M, &A_N );
		if( (A_M!=AM_sum) || (A_N!=AN_sum) ) {	
			Stg_SETERRQ( PETSC_ERR_ARG_SIZ, "Re-used matrix has wrong global dimensions" );
		}
	}	


	/* < prealloacte space  > */
	MatGetOwnershipRanges( (*A), &ranges_A );
	if(ranges_A==PETSC_NULL){
		MatSetUp(*A);
		MatGetOwnershipRanges( (*A), &ranges_A );
	}
	MPI_Comm_size( comm, &nprocs );
	MPI_Comm_rank( comm, &rank );
	
#if 0

	MatGetVecs( (*A), &local_count, &off_proc_count );
	VecSet( local_count, 0.0 );
	VecSet( off_proc_count, 0.0 );
	
	PetscGetTime(&t0);
	one = 1.0;
	for( i=0; i<block->nr; i++ ) {
		for( j=0; j<block->nc; j++ ) {
			
			if( block->m[i][j] != PETSC_NULL ) {
				
				MatGetSize( block->m[i][j], &M, &N );
				MatGetLocalSize( block->m[i][j], &m, &n );
				MatGetOwnershipRange( block->m[i][j], &start_row, &end_row );
				for( si=start_row; si<end_row; si++ ) {
					row = si + AM[i];
					
					MatGetRow( block->m[i][j], si, &ncols, &cols, &vals );
					
					if( cols != PETSC_NULL ) {
						PetscMalloc( ncols*sizeof(PetscInt), &shifted_cols );
						for( sj=0; sj<ncols; sj++ ) {
							shifted_cols[sj] = cols[sj] + AN[j];
						}
						
						/* now scan and check each column */
						for( sj=0; sj<ncols; sj++ ) {
							for( p=0; p<nprocs; p++ ) {
								if( (shifted_cols[sj] >= ranges_A[p]) && (shifted_cols[sj] <ranges_A[p+1]) ) {
									break;
								}
							}
							if( p==rank ) {
								VecSetValue( local_count, row, one, ADD_VALUES );
							}
							else {
								VecSetValue( off_proc_count, row, one, ADD_VALUES );
							}
						}
						
						PetscFree( shifted_cols );
						MatRestoreRow( block->m[i][j], row, &ncols, &cols, &vals );
					}
				}
			}
			
		}
	}
	
	/* finish off building you integer lists for nnz and onnz */
	VecAssemblyBegin(local_count);
	VecAssemblyEnd(local_count);
	
	VecAssemblyBegin(off_proc_count);
	VecAssemblyEnd(off_proc_count);
        PetscGetTime(&t1);
        PetscPrintf(PETSC_COMM_WORLD,"  MergeSubBlocks: dnnz,onnz (%lf sec) \n", t1 -t0 );

	PetscGetTime(&t0);	
	VecGetArray( local_count, &s_nnz );
	VecGetArray( off_proc_count, &s_onnz );
	
	VecGetLocalSize( off_proc_count, &m );
	PetscMalloc( m*sizeof(PetscInt), &nnz );
	PetscMalloc( m*sizeof(PetscInt), &onnz );
	
	for( j=0; j<m; j++ ){
		nnz[j]  = (PetscInt)s_nnz[j];
		onnz[j] = (PetscInt)s_onnz[j];
	}
	
	VecRestoreArray( off_proc_count, &s_onnz );
	VecRestoreArray( local_count, &s_nnz );
	
	Stg_VecDestroy( &local_count );
	Stg_VecDestroy( &off_proc_count );

	
	/* call preallocate */
	preallocated = PETSC_FALSE;

	/* This seems to be PETSC_TRUE when matrix is created */
/*
	preallocated = (*A)->preallocated;
*/
	
	/* nz_allocated seems to get reset after each assembly */
/*
	if( (*A)!=PETSC_NULL ) {
		MatInfo info;

		MatGetInfo( (*A) , MAT_LOCAL, &info );
		if( info.nz_allocated == 0 ) {
			preallocated = PETSC_TRUE;
		}
	}
*/	

        if( (*A)!=PETSC_NULL ) {
                if( (*A)->num_ass != 0 ) {
                        preallocated = PETSC_TRUE;
                }
        }	


	if( preallocated == PETSC_FALSE ) {
		PetscTypeCompare( (PetscObject)(*A), MATSEQAIJ, &is_aij );
		if(is_aij==PETSC_TRUE) {
			MatSeqAIJSetPreallocation( *A, PETSC_NULL, nnz );
		}
		PetscTypeCompare( (PetscObject)(*A), MATMPIAIJ, &is_aij );
		if(is_aij==PETSC_TRUE) {
			MatMPIAIJSetPreallocation( *A, PETSC_NULL,nnz, PETSC_NULL,onnz );
		}
	}
	PetscFree(nnz);
	PetscFree(onnz);
        PetscGetTime(&t1);
        PetscPrintf(PETSC_COMM_WORLD,"  MergeSubBlocks: dnnz,onnz (II) (%lf sec) \n", t1 -t0 );
	/* < end preallocation > */
#endif	


	MergeSubBlocksPreallocateMatAIJ( Ab, *A );

	MatSetOption( *A, MAT_USE_HASH_TABLE, PETSC_TRUE );	


        time_msv = time_mgr = 0.0;
        PetscGetTime(&t0);

	for( i=0; i<block->nr; i++ ) {
		for( j=0; j<block->nc; j++ ) {
			
			if( block->m[i][j] != PETSC_NULL ) {
				
				MatGetSize( block->m[i][j], &M, &N );
				MatGetLocalSize( block->m[i][j], &m, &n );
				MatGetOwnershipRange( block->m[i][j], &start_row, &end_row );
				for( si=start_row; si<end_row; si++ ) {
					row = si + AM[i];
				
					PetscGetTime(&dt0);	
					MatGetRow( block->m[i][j], si, &ncols, &cols, &vals );
					PetscGetTime(&dt1);					
					time_mgr = time_mgr + (dt1-dt0);

					if( cols != PETSC_NULL ) {
						PetscMalloc( ncols*sizeof(PetscInt), &shifted_cols );
						for( sj=0; sj<ncols; sj++ ) {
							shifted_cols[sj] = cols[sj] + AN[j];
						}
						PetscGetTime(&dt0);
						MatSetValues( *A, 1, &row, ncols, shifted_cols, vals, addv );
						PetscGetTime(&dt1);
                                                time_msv = time_msv + (dt1-dt0);						

						PetscFree( shifted_cols );
					}
	                                MatRestoreRow( block->m[i][j], row, &ncols, &cols, &vals );
				}
			}
			
		}
	}

        PetscGetTime(&t1);
	if(has_info==PETSC_TRUE) {
          PetscPrintf(PETSC_COMM_WORLD,"  MatMergeSubBlocks: merge (%lf sec) \n", t1 -t0 );
          PetscPrintf(PETSC_COMM_WORLD,"  MatMergeSubBlocks:   merge->MatSetValues (%lf sec) \n", time_msv );
          PetscPrintf(PETSC_COMM_WORLD,"  MatMergeSubBlocks:   merge->MatGetRow (%lf sec) \n", time_mgr );
	}

#if 0
{
  PetscScalar *all_values;
  PetscInt *_start_row, *_end_row;
  PetscInt offset;

  PetscMalloc( sizeof(PetscInt)*AM_sum, &_start_row );
  PetscMalloc( sizeof(PetscInt)*AN_sum, &_end_row );


  for( i=0; i<block->nr; i++ ) {
    for( j=0; j<block->nc; j++ ) {
      if( block->m[i][j] == PETSC_NULL ) {  continue;  }
        MatGetOwnershipRange( block->m[i][j], &_start_row[i], &_end_row[i] );
        break;
    }
  }

        time_msv = time_mgr = 0.0;
        PetscGetTime(&t0);

	PetscMalloc( sizeof(PetscScalar)*AN_sum, &all_values );
        for( i=0; i<block->nr; i++ ) {

	        for( si=_start_row[i]; si<_end_row[i]; si++ ) {
        	        row = si + AM[i];        

			offset = 0;
	        	for( j=0; j<block->nc; j++ ) {
				if( block->m[i][j] == PETSC_NULL ) {  continue;  }

				PetscGetTime(&dt0);
				MatGetRow( block->m[i][j], si, &ncols, PETSC_NULL, &vals );			
				PetscGetTime(&dt1);                                     
                                time_mgr = time_mgr + (dt1-dt0);

				PetscMemcpy( &all_values[offset], vals, sizeof(PetscScalar)*ncols );
				offset = offset + ncols;
			}
			PetscGetTime(&dt0);
			MatSetValuesRow(*A,row,all_values);
			PetscGetTime(&dt1);
                        time_msv = time_msv + (dt1-dt0); 
		}

	}

	PetscGetTime(&t1);

  PetscFree( all_values );
  PetscFree( _start_row );
  PetscFree( _end_row );

        PetscPrintf(PETSC_COMM_WORLD,"  MatMergeSubBlocks: merge (%lf sec) \n", t1 -t0 );
        PetscPrintf(PETSC_COMM_WORLD,"  MatMergeSubBlocks:   merge->MatSetValues (%lf sec) \n", time_msv );
        PetscPrintf(PETSC_COMM_WORLD,"  MatMergeSubBlocks:   merge->MatGetRow (%lf sec) \n", time_mgr );

}
#endif

	MatAssemblyBegin( *A, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( *A, MAT_FINAL_ASSEMBLY );
	
	PetscFree( AM );
	PetscFree( AN );
	
	PetscGetTime(&tt1);
	if(has_info==PETSC_TRUE) {
	  PetscPrintf(PETSC_COMM_WORLD, "  MatMergeSubBlocks: total (%lf sec) \n", tt1-tt0 );	
	}

	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "MatBlockMergeSubBlocks"
PetscErrorCode MatBlockMergeSubBlocks(  Mat Ab, InsertMode addv, const MatType mtype, Mat *A  )
{
	PetscErrorCode ierr,(*f)(Mat,InsertMode,const MatType,Mat*);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)Ab,"MatBlockMergeSubBlocks_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(Ab,addv,mtype,A);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}


/* Constructors */
#undef __FUNCT__  
#define __FUNCT__ "MatSetOps_Block"
PetscErrorCode MatSetOps_Block( struct _MatOps* ops )
{
	PetscFunctionBegin;

	/* 0*/
	ops->setvalues  = 0; // MatSetValues_EMPTY;
	ops->getrow     = 0; // MatGetRow_EMPTY;
	ops->restorerow = 0; // MatRestoreRow_EMPTY;
	ops->mult       = MatMult_Block;
	ops->multadd    = 0;
	/* 5*/
	ops->multtranspose    = MatMultTranspose_Block;
	ops->multtransposeadd = 0;
	ops->solve            = 0; // MatSolve_EMPTY;
	ops->solveadd         = 0; // MatSolveAdd_EMPTY;
	ops->solvetranspose   = 0; // MatSolveTranspose_EMPTY;
	/*10*/
	ops->solvetransposeadd = 0; // MatSolveTransposeAdd_EMPTY;
	ops->lufactor          = 0; // MatLUFactor_EMPTY;
	ops->choleskyfactor    = 0; // MatCholeskyFactor_EMPTY;
	ops->sor               = 0; // MatSOR_EMPTY;
	ops->transpose         = 0; // MatTranspose_EMPTY;
	/*15*/
	ops->getinfo       = 0; // MatGetInfo_EMPTY;
	ops->equal         = 0; // MatEqual_EMPTY;
	ops->getdiagonal   = 0; // MatGetDiagonal_EMPTY;
	ops->diagonalscale = 0; // MatDiagonalScale_EMPTY;
	ops->norm          = 0; // MatNorm_EMPTY;
	/*20*/
	ops->assemblybegin = 0;
	ops->assemblyend   = 0;
	ops->setoption     = 0; // MatSetOption_EMPTY;
	ops->zeroentries   = MatZeroEntries_Block; // MatZeroEntries_EMPTY;
	/*24*/
	ops->zerorows               = 0; // MatZeroRows_EMPTY;
	ops->lufactorsymbolic       = 0; // MatLUFactorSymbolic_EMPTY;
	ops->lufactornumeric        = 0; // MatLUFactorNumeric_EMPTY;
	ops->choleskyfactorsymbolic = 0; // MatCholeskyFactorSymbolic_EMPTY;
	ops->choleskyfactornumeric  = 0; // MatCholeskyFactorNumeric_EMPTY;
	/*29*/
	ops->setuppreallocation = 0; // MatSetUpPreallocation_EMPTY;
	ops->ilufactorsymbolic  = 0; // MatILUFactorSymbolic_EMPTY;
	ops->iccfactorsymbolic  = 0; // MatICCFactorSymbolic_EMPTY;
	ops->getarray           = 0; // MatGetArray_EMPTY;
	ops->restorearray       = 0; // MatRestoreArray_EMPTY;
	/*34*/
	ops->duplicate     = MatDuplicate_Block; // MatDuplicate_EMPTY;
	ops->forwardsolve  = 0; // MatForwardSolve_EMPTY;
	ops->backwardsolve = 0; // MatBackwardSolve_EMPTY;
	ops->ilufactor     = 0; // MatILUFactor_EMPTY;
	ops->iccfactor     = 0; // MatICCFactor_EMPTY;
	/*39*/
	ops->axpy            = 0; // MatAXPY_EMPTY;
	ops->getsubmatrices  = 0; // MatGetSubMatrices_EMPTY;
	ops->increaseoverlap = 0; // MatIncreaseOverlap_EMPTY;
	ops->getvalues       = 0; // MatGetValues_EMPTY;
	ops->copy            = 0; // MatCopy_EMPTY;
	/*44*/
	ops->getrowmax   = 0; // MatGetRowMax_EMPTY;
	ops->scale       = 0;
	ops->shift       = 0; // MatShift_EMPTY;
	ops->diagonalset = 0; // MatDiagonalSet_EMPTY;
	////ops->dummy       = 0; // MatILUDTFactor_EMPTY;
	/*49*/
	ops->setblocksize    = 0; // MatSetBlockSize_EMPTY;
	ops->getrowij        = 0; // MatGetRowIJ_EMPTY;
	ops->restorerowij    = 0; // MatRestorRowIJ_EMPTY;
	ops->getcolumnij     = 0; // MatGetColumnIJ_EMPTY;
	ops->restorecolumnij = 0; // MatRestoreColumnIJ_EMPTY;
	/*54*/
	ops->fdcoloringcreate = 0; // MatFDColoringCreate_EMPTY;
	ops->coloringpatch    = 0; // MatColoringPatch_EMPTY;
	ops->setunfactored    = 0; // MatSetUnfactored_EMPTY;
	ops->permute          = 0; // MatPermute_EMPTY;
	ops->setvaluesblocked = 0; // MatSetValuesBlocked_EMPTY;
	/*59*/
	ops->getsubmatrix  = 0; // MatGetSubMatrix_EMPTY;
	ops->destroy       = MatDestroy_Block;
	ops->view          = MatView_Block;
	ops->convertfrom   = 0; // MatConvertFrom_EMPTY;
	ops->usescaledform = 0; // MatUseScaledForm_EMPTY;
	/*64*/
	ops->scalesystem             = 0; // MatScaleSystem_EMPTY;
	ops->unscalesystem           = 0; // MatUnScaleSystem_EMPTY;
	ops->setlocaltoglobalmapping = 0; // MatSetLocalToGlobalMapping_EMPTY;
	ops->setvalueslocal          = 0; // MatSetValuesLocal_EMPTY;
	ops->zerorowslocal           = 0; // MatZeroRowsLocal_EMPTY;
	/*69*/
	ops->getrowmaxabs    = 0; // MatGetRowMaxAbs_EMPTY;
	ops->getrowminabs    = 0; // 
	ops->convert         = MatConvert_Block;
	ops->setcoloring     = 0; // MatSetColoring_EMPTY;
	ops->setvaluesadic   = 0; // MatSetValuesAdic_EMPTY;
	/* 74 */
	ops->setvaluesadifor = 0; // MatSetValuesAdifor_EMPTY;
	ops->fdcoloringapply              = 0; // MatFDColoringApply_EMPTY;
	ops->setfromoptions               = 0;
	ops->multconstrained              = 0; // MatMultConstrained_EMPTY;
	ops->multtransposeconstrained     = 0; // MatMultTransposeConstrained_EMPTY;
	/*79*/
	//ops->permutesparsify = 0; // MatPermuteSparsify_EMPTY;
	ops->mults           = 0; // MatMults_EMPTY;
	ops->solves          = 0; // MatSolves_EMPTY;
	ops->getinertia      = 0; // MatGetInertia_EMPTY;
	ops->load            = 0; // Stg_MatLoad_EMPTY;
	/*84*/
	ops->issymmetric             = 0; // MatIsSymmetric_EMPTY;
	ops->ishermitian             = 0; // MatIsHermitian_EMPTY;
	ops->isstructurallysymmetric = 0; // MatIsStructurallySymmetric_EMPTY;
	//ops->dummy4                  = 0; // MatDummay_EMPTY;
	ops->getvecs                 = MatGetVecs_Block;
	/*89*/
	ops->matmult         = MatMatMult_Block; // MatMatMult_EMPTY;
	ops->matmultsymbolic = 0; // MatMatMultSymbolic_EMPTY;
	ops->matmultnumeric  = 0; // MatMatMultNumeric_EMPTY;
	ops->ptap            = 0; // MatPtAP_EMPTY;
	ops->ptapsymbolic    = 0; // MatPtAPSymbolic_EMPTY;
	/*94*/
	ops->ptapnumeric              = 0; // MatPtAPNumeric_EMPTY;
	ops->matmulttranspose         = 0; // MatMatMultTranspose_EMPTY;
	ops->matmulttransposesymbolic = 0; // MatMatMultTransposeSymbolic_EMPTY;
	ops->matmulttransposenumeric  = 0; // MatMatMultTransposeNumeric_EMPTY;
	ops->ptapsymbolic_seqaij      = 0; // MatPtAPSymbolic_SEQAIJ_EMPTY;
	/*99*/
	ops->ptapnumeric_seqaij  = 0; // MatPtAPNumeric_SEQAIJ_EMPTY;
	ops->ptapsymbolic_mpiaij = 0; // MatPtAPSymbolic_MPIAIJ_EMPTY;
	ops->ptapnumeric_mpiaij  = 0; // MatPtAPNumeric_MPIAIJ_EMPTY;
	ops->conjugate           = 0; // MatConjugate_EMPTY;
	ops->setsizes            = MatSetSizes_Block; // MatSetSizes_EMPTY;
	/*104*/
	ops->setvaluesrow              = 0; // MatSetValuesRow_EMPTY;
	ops->realpart                  = 0; // MatRealPart_EMPTY;
	ops->imaginarypart             = 0; // MatImaginaryPart_EMPTY;
	ops->getrowuppertriangular     = 0; // MatGetRowUpperTriangular_EMPTY;
	ops->restorerowuppertriangular = 0; // MatRestoreRowUpperTriangular_EMPTY;
	/*109*/
	ops->matsolve           = 0; // MatMatSolve_EMPTY;
	ops->getredundantmatrix = 0; // MatGetRedundantMatrix_EMPTY;
	ops->getrowmin          = 0; // MatGetRowMin_EMPTY;
	ops->getcolumnvector    = 0;  // MatGetColumnVector_EMPTY;
	ops->missingdiagonal    = 0;
	/* 114 */
	ops->getseqnonzerostructure = 0;
	ops->create                 = 0;
	ops->getghosts              = 0;
	//////ops->dummy2                 = 0;
	//////ops->dummy3                 = 0;
	/* 119 */
	ops->multdiagonalblock = 0;
	ops->hermitiantranspose = 0;
	ops->multhermitiantranspose = 0;
	ops->multhermitiantransposeadd = 0;
	
	PetscFunctionReturn(0);
}


/*
The following operations are not to used with MatBlock
 - MatGetOwnershipRange(A,&s,&e):  Result, s=0,e=0.
 - MatGetOwnershripRanges(A,&r):   Result, error.
*/
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "MatCreate_Block"
PetscErrorCode MatCreate_Block( Mat A )
{
	Mat_Block      *s;
	PetscErrorCode ierr;
	PetscMPIInt size;
	
	PetscFunctionBegin;
	
	ierr = PetscMalloc( sizeof(Mat_Block), &s );CHKERRQ(ierr);
	ierr = PetscMemzero( s, sizeof(Mat_Block) );CHKERRQ(ierr);
	
	
	/* define the operations */
	MatSetOps_Block( A->ops );
//	ierr = PetscMemcpy( A->ops, &MatBlockOps_Values,sizeof(struct _MatOps) );CHKERRQ(ierr);
	
//	A->precision       = PETSC_SCALAR;
	A->data            = (void*)s;
	
	//A->factor           = MAT_FACTOR_NONE;
	//A->mapping          = 0;
	A->spptr            = 0;
	A->same_nonzero     = PETSC_FALSE;
	A->assembled        = PETSC_FALSE;
	
	
	s->setup_called = PETSC_FALSE;
	s->nr = s->nc   = -1;
	s->m            = PETSC_NULL;
	s->bstruct      = NULL_BLOCK;
	
	PetscObjectChangeTypeName((PetscObject) A, "block" );
	
	/* 
	For whatever reason this breaks my parallel tests. I really only
	needed to setup the maps to make map->range not NULL. This requirment
	was needed as MatDuplicate copys the map and complains if the range is NULL.
	I think I'll added a little hack in MatDuplicate_Block to shut this up. 09 Jan, 2008
	*/
	
	/*
	M = A->rmap.N;		N = A->cmap.N;
	
	if (A->rmap.bs == -1) A->rmap.bs = 1;
	PetscMapSetUp( &A->rmap );
	PetscMapSetLocalSize( &A->rmap, M );
	PetscMapSetSize( &A->rmap, M );
	
	if (A->cmap.bs == -1) A->cmap.bs = 1;
	PetscMapSetUp( &A->cmap );
	PetscMapSetLocalSize( &A->cmap, N );
	PetscMapSetSize( &A->cmap, N );
	*/
	
	/* The interface MatDuplicate() will try to copy rmap->range, but the manner
	in which we constructed the block objects leaves that NULL. So here I'll
	allocate it just to keep the interface function happpy. */
	MPI_Comm_size( ((PetscObject)A)->comm, &size );
	PetscMalloc( (size+1)*sizeof(PetscInt), &A->rmap->range );
	PetscMalloc( (size+1)*sizeof(PetscInt), &A->cmap->range );
	
	
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)A,"MatBlockSetValues_C",
				"MatBlockSetValues_Block",
				MatBlockSetValues_Block);CHKERRQ(ierr);
		
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)A,"MatBlockSetValue_C",
				"MatBlockSetValue_Block",
				MatBlockSetValue_Block);CHKERRQ(ierr);
	
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)A,"MatBlockGetSubMatrix_C",
				"MatBlockGetSubMatrix_Block",
				MatBlockGetSubMatrix_Block);CHKERRQ(ierr);
	
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)A,"MatBlockGetSubMatrices_C",
				"MatBlockGetSubMatrices_Block",
				MatBlockGetSubMatrices_Block);CHKERRQ(ierr);
	
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)A,"MatBlockRestoreSubMatrices_C",
				"MatBlockRestoreSubMatrices_Block",
				MatBlockRestoreSubMatrices_Block);CHKERRQ(ierr);
	
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)A,"MatBlockCreateSubMatrix_C",
				"MatBlockCreateSubMatrix_Block",
				MatBlockCreateSubMatrix_Block);CHKERRQ(ierr);
	
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)A,"MatBlockMergeSubBlocks_C",
				"MatBlockMergeSubBlocks_Block",
				MatBlockMergeSubBlocks_Block);CHKERRQ(ierr);

        ierr = PetscObjectComposeFunctionDynamic((PetscObject)A,"MatCreateFromMatBlock_C",
                                "MatCreateFromMatBlock_Block",
                                MatCreateFromMatBlock_Block);CHKERRQ(ierr);

	
	PetscFunctionReturn(0);
}
EXTERN_C_END
