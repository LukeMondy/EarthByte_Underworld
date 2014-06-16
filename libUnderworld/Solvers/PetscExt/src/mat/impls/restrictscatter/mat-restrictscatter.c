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

#include "mat-restrictscatter-impl.h"
#include "mat-restrictscatter-ops.h"
#include "private/mat/mat-restrictscatter.h"
#include "petscext_utils.h"



PetscErrorCode Stg_MatDestroy_RestrictScatter( Mat A )
{
	PetscErrorCode ierr;
	Mat_RestrictScatter s = (Mat_RestrictScatter)A->data;
	
	PetscFunctionBegin;
	
	/* release implemenation internals */
	Stg_MatDestroy( & s->A );
	Stg_ISDestroy( & s->row );
	Stg_ISDestroy( & s->col );
	Stg_VecScatterDestroy( & s->scat_x );
	Stg_VecScatterDestroy( & s->scat_y );
	Stg_VecDestroy( & s->x_work );
	Stg_VecDestroy( & s->y_work );
	
	/* release implementation data pointer */
	ierr=PetscFree( s );CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


/*----------------------------------------------------------------------------*/

/*  We permit the row/col IS to change. */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "_MatRestrictScatterSetOperator_RestrictScatter"
PetscErrorCode _MatRestrictScatterSetMatIS_RestrictScatter( Mat rsA, Mat A, IS row, IS col )
{
	Mat_RestrictScatter   s = (Mat_RestrictScatter)rsA->data;
	MPI_Comm comm;
	PetscTruth is_updated = PETSC_FALSE;
	
	
	PetscFunctionBegin;
	
	if(A!=PETSC_NULL) {
		if( s->A != PETSC_NULL) {
			Stg_SETERRQ( PETSC_ERR_SUP, "You cannot change the original operator A");
		}
		s->A = A;
		PetscObjectReference( (PetscObject)A );
	}
	
	if(row!=PETSC_NULL) {
		PetscObjectReference( (PetscObject)row );
		if( s->row != PETSC_NULL ) {
			Stg_ISDestroy( & s->row );
			is_updated = PETSC_TRUE;
		}
		s->row = row;
	}
	if(col!=PETSC_NULL) {
		PetscObjectReference( (PetscObject)col );
		if( s->col != PETSC_NULL ) {
			Stg_ISDestroy( & s->col );
			is_updated = PETSC_TRUE;
		}
		s->col = col;
	}
	
	PetscObjectGetComm( (PetscObject)rsA, &comm );
	
	/* create internal objects if required */
	if(s->y_work==PETSC_NULL) {
		MatGetVecs( s->A, 0, &s->y_work );
	}
	if(s->x_work==PETSC_NULL) {
		MatGetVecs( s->A, &s->x_work, 0 );
	}
	if(s->scat_x==PETSC_NULL) {
		Vec xx;
		
		MatGetVecs( rsA, &xx, 0 );
		VecScatterCreate( xx,PETSC_NULL, s->x_work,s->col, &s->scat_x );
		Stg_VecDestroy( & xx );
	}
	if(s->scat_y==PETSC_NULL) {
		Vec yy;
		
		MatGetVecs( rsA, 0, &yy );
		VecScatterCreate( s->y_work,s->row , yy,PETSC_NULL, &s->scat_y );
		Stg_VecDestroy( & yy );
	}
	
	if( is_updated == PETSC_TRUE ) {
		Vec xx,yy;
		
		/* free everything */
		Stg_VecDestroy( & s->x_work );
		Stg_VecDestroy( & s->y_work );
		Stg_VecScatterDestroy( & s->scat_x );
		Stg_VecScatterDestroy( & s->scat_y );
		
		/* re-create everything */
		MatGetVecs( s->A, &s->x_work, &s->y_work );
		
		MatGetVecs( rsA, &xx, &yy );
		VecScatterCreate( xx,PETSC_NULL, s->x_work,s->col, &s->scat_x );
		VecScatterCreate( s->y_work,s->row , yy,PETSC_NULL, &s->scat_y );
		
		Stg_VecDestroy( & xx );
		Stg_VecDestroy( & yy );
	}
	
	
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "MatRestrictScatterSetOperator"
PetscErrorCode MatRestrictScatterSetMatIS( Mat rsA, Mat A, IS row, IS col )
{
	PetscErrorCode (*f)(Mat,Mat,IS,IS);
	
	
	PetscObjectQueryFunction( (PetscObject)rsA, "MatRestrictScatterSetMatIS_C", (void (**)(void))&f );
	if( f ) {
		(*f)(rsA, A,row,col );
	} else {
		Stg_SETERRQ(PETSC_ERR_ARG_WRONG,"Cannot set Mat and IS for this matrix");
	}
	
	PetscFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "MatCreateRestrictScatter"
PetscErrorCode MatCreateRestrictScatter( MPI_Comm comm, Mat A, IS row, IS col, Mat *rsA )
{	
	MPI_Comm commA;
	PetscInt M,N, m,n;
	Vec L,R;
	PetscMPIInt size;
	VecType vt;
	
	PetscObjectGetComm( (PetscObject)A, &commA );
	if( comm != commA ) {
//		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "Communicators on A must match comm passed in \n" );
	}
	MPI_Comm_size( commA, &size );
	vt = VECSEQ;
	if( size > 1 ) {
		vt = VECMPI;
	}
	
	ISGetSize( row, &M );
	ISGetLocalSize( row, &m );
	
	ISGetSize( col, &N );
	ISGetLocalSize( col, &n );
	
	VecCreate( commA, &L );
	VecSetSizes( L, PETSC_DECIDE, M );
	VecSetType(L,vt);
	VecGetLocalSize( L, &m );
	
	VecCreate( commA, &R );
	VecSetSizes( R, PETSC_DECIDE, N );
	VecSetType(R,vt);
	VecGetLocalSize( R, &n );
	
	Stg_VecDestroy( &L);
	Stg_VecDestroy( &R);
	
	MatCreate( comm, rsA );
	MatSetSizes( *rsA, m,n, M,N );
	MatSetType( *rsA, "restrictscatter" );
	MatRestrictScatterSetMatIS( *rsA, A, row, col );
	
	MatAssemblyBegin(*rsA,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*rsA,MAT_FINAL_ASSEMBLY);
	
	PetscFunctionReturn(0);
}


/* Constructors */
PetscErrorCode MatSetOps_RestrictScatter( struct _MatOps* ops )
{
	PetscFunctionBegin;
	
	/* 0*/
	ops->setvalues  = 0;
	ops->getrow     = 0;
	ops->restorerow = 0;
	ops->mult       = MatMult_RestrictScatter;
	ops->multadd    = 0;
	/* 5*/
	ops->multtranspose    = 0;
	ops->multtransposeadd = 0;
	ops->solve            = 0;
	ops->solveadd         = 0;
	ops->solvetranspose   = 0;
	/*10*/
	ops->solvetransposeadd = 0;
	ops->lufactor          = 0;
	ops->choleskyfactor    = 0; 
	ops->sor               = 0; 
	ops->transpose         = 0; 
	/*15*/
	ops->getinfo       = 0; 
	ops->equal         = 0; 
	ops->getdiagonal   = 0; 
	ops->diagonalscale = 0; 
	ops->norm          = 0; 
	/*20*/
	ops->assemblybegin = 0;
	ops->assemblyend   = 0;
	ops->setoption     = 0; 
	ops->zeroentries   = 0; 
	/*24*/
	ops->zerorows               = 0; 
	ops->lufactorsymbolic       = 0; 
	ops->lufactornumeric        = 0; 
	ops->choleskyfactorsymbolic = 0; 
	ops->choleskyfactornumeric  = 0; 
	/*29*/

	ops->ilufactorsymbolic  = 0; 
	ops->iccfactorsymbolic  = 0; 


	/*34*/
	ops->duplicate     = 0; 
	ops->forwardsolve  = 0; 
	ops->backwardsolve = 0; 
	ops->ilufactor     = 0; 
	ops->iccfactor     = 0; 
	/*39*/
	ops->axpy            = 0; 
	ops->getsubmatrices  = 0; 
	ops->increaseoverlap = 0; 
	ops->getvalues       = 0; 
	ops->copy            = 0; 
	/*44*/
	ops->getrowmax   = 0; 
	ops->scale       = 0;
	ops->shift       = 0; 
	ops->diagonalset = 0; 
	//ops->dummy       = 0; 
	/*49*/

	ops->getrowij        = 0; 
	ops->restorerowij    = 0; 
	ops->getcolumnij     = 0; 
	ops->restorecolumnij = 0; 
	/*54*/
	ops->fdcoloringcreate = 0; 
	ops->coloringpatch    = 0; 
	ops->setunfactored    = 0; 
	ops->permute          = 0; 
	ops->setvaluesblocked = 0; 
	/*59*/
	ops->getsubmatrix  = 0; 
	ops->destroy       = Stg_MatDestroy_RestrictScatter;
	ops->view          = 0;
	ops->convertfrom   = 0; 

	/*64*/


	ops->setlocaltoglobalmapping = 0; 
	ops->setvalueslocal          = 0; 
	ops->zerorowslocal           = 0; 
	/*69*/
	ops->getrowmaxabs    = 0; 
	ops->getrowminabs    = 0; // 
	ops->convert         = 0; 
	ops->setcoloring     = 0; 

	/* 74 */
	ops->setvaluesadifor = 0; 
	ops->fdcoloringapply              = 0; 
	ops->setfromoptions               = 0;
	ops->multconstrained              = 0; 
	ops->multtransposeconstrained     = 0; 
	/*79*/
	//ops->permutesparsify = 0; 
	ops->mults           = 0; 
	ops->solves          = 0; 
	ops->getinertia      = 0; 
	ops->load            = 0; 
	/*84*/
	ops->issymmetric             = 0; 
	ops->ishermitian             = 0; 
	ops->isstructurallysymmetric = 0; 
	//ops->dummy4                  = 0; 
	ops->getvecs                 = 0;
	/*89*/
	ops->matmult         = 0; 
	ops->matmultsymbolic = 0; 
	ops->matmultnumeric  = 0; 
	ops->ptap            = 0; 
	ops->ptapsymbolic    = 0; 
	/*94*/
	ops->ptapnumeric              = 0; 




	/*99*/



	ops->conjugate           = 0; 

	/*104*/
	ops->setvaluesrow              = 0; 
	ops->realpart                  = 0; 
	ops->imaginarypart             = 0; 
	ops->getrowuppertriangular     = 0; 
	ops->restorerowuppertriangular = 0; 
	/*109*/
	ops->matsolve           = 0; 
	ops->getredundantmatrix = 0; 
	ops->getrowmin          = 0; 
	ops->getcolumnvector    = 0;  
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

MatResrictScatter:
Given A, MatRestrictScatter defines the mat-vec product 
on y = B x, where B is a sub-matrix of A. It is similar to
MatScatter, except the dimensions of B are smaller than A.

The shape and values defining B are given by an index two
sets defining the row and colummns of A we wish to use.

Note that the cost of computing the mat-vec product of y = B x,
will cost approximately the same (actually, slightly more) than
the cost of computing the mat-vec prodcut y' = A x'.

This may be a problem, however for very specialised matrix-free
implementations, it might be a suitable option. In this case,
the matrix-free routine will (should be) highly optimised so
the extra overhead my not be too bad, especially if B is not
significantly smaller than A. In addition, this code also prevents
users from having to write another optimised matrix-free mat-vec
routine to define the particular restricted mat-vec product 
derived from the operator A.

                     y = B x
      scat_y                          scat_x
y <============ y_work = A x_work <============= x

*/
EXTERN_C_BEGIN
PetscErrorCode MatCreate_RestrictScatter( Mat A )
{
	Mat_RestrictScatter    s;
	PetscErrorCode  ierr;
	
	
	PetscFunctionBegin;
	
	/* define the operations */
	MatSetOps_RestrictScatter( A->ops );
	
//	A->precision       = PETSC_SCALAR;
	//A->factor           = MAT_FACTOR_NONE;
	//A->mapping          = 0;
	A->spptr            = 0;
	A->same_nonzero     = PETSC_FALSE;
	A->assembled        = PETSC_FALSE;
	
	/* allocate and set pointer for implememtation data */
	ierr = PetscMalloc( sizeof(struct _Mat_RestrictScatter), &s );CHKERRQ(ierr);
	ierr = PetscMemzero( s, sizeof(struct _Mat_RestrictScatter) );CHKERRQ(ierr);
	A->data            = (void*)s;
	
	/* Init the implementation data */
	s->A = PETSC_NULL;
	s->row = PETSC_NULL;
	s->col = PETSC_NULL;
	s->scat_x = PETSC_NULL;
	s->scat_y = PETSC_NULL;
	s->x_work = PETSC_NULL;
	s->y_work = PETSC_NULL;
	
	/* Set type */
	PetscObjectChangeTypeName((PetscObject) A, "restrictscatter" );
	
	
	/* Compose functions */
	PetscObjectComposeFunctionDynamic( (PetscObject)A, "MatRestrictScatterSetMatIS_C",
			"_MatRestrictScatterSetMatIS_RestrictScatter", _MatRestrictScatterSetMatIS_RestrictScatter );
	
	
	PetscFunctionReturn(0);
}
EXTERN_C_END
