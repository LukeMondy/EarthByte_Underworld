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

#include "mat-symbolic_transpose-impl.h"
#include "mat-symbolic_transpose-ops.h"
#include "private/mat/mat-symbolic_transpose.h"



#undef __FUNCT__
#define __FUNCT__ "Stg_MatDestroy_SymTrans"
PetscErrorCode Stg_MatDestroy_SymTrans( Mat A )
{
	PetscErrorCode ierr;
	Mat_SymTrans      s = (Mat_SymTrans)A->data;
	
	PetscFunctionBegin;
	
	/* release implemenation internals */
	Stg_MatDestroy( &s->A );
	
	/* release implementation data pointer */
	ierr=PetscFree( s );CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


/*----------------------------------------------------------------------------*/
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "_MatSymTransSetOperator_SymTrans"
PetscErrorCode _MatSymTransSetOperator_SymTrans( Mat symAt, Mat A )
{
	Mat_SymTrans   s = (Mat_SymTrans)symAt->data;
	PetscTruth is_sym;
	
	PetscFunctionBegin;
	
	is_sym = PETSC_FALSE;
	Stg_PetscTypeCompare( (PetscObject)A, "symtrans", &is_sym );
	if( is_sym == PETSC_TRUE ) {
		Stg_SETERRQ(PETSC_ERR_ARG_WRONG,"The operator for MatSymTrans must not be symbolic.");
	}

	s->A = A;
	PetscObjectReference( (PetscObject)A );
	
	MatAssemblyBegin( symAt, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( symAt, MAT_FINAL_ASSEMBLY );
	
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "MatSymTransSetOperator"
PetscErrorCode MatSymTransSetOperator( Mat symAt, Mat A )
{
	PetscErrorCode (*f)(Mat,Mat);
	
	
	PetscObjectQueryFunction( (PetscObject)symAt, "MatSymTransSetOperator_C", (void (**)(void))&f );
	if( f ) {
		(*f)(symAt, A );
	} else {
		Stg_SETERRQ(PETSC_ERR_ARG_WRONG,"Cannot set operator for this matrix");
	}
	
	PetscFunctionReturn(0);
}

EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "_MatSymTransGetOperator_SymTrans"
PetscErrorCode _MatSymTransGetOperator_SymTrans( Mat symAt, Mat *A )
{
	Mat_SymTrans   s = (Mat_SymTrans)symAt->data;
	
	PetscFunctionBegin;
	*A = s->A;
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "MatSymTransGetOperator"
PetscErrorCode MatSymTransGetOperator( Mat symAt, Mat *A )
{
	PetscErrorCode (*f)(Mat,Mat*);
	
	
	PetscObjectQueryFunction( (PetscObject)symAt, "MatSymTransGetOperator_C", (void (**)(void))&f );
	if( f ) {
		(*f)(symAt, A );
	} else {
		Stg_SETERRQ(PETSC_ERR_ARG_WRONG,"Cannot get operator for this matrix");
	}
	
	PetscFunctionReturn(0);
}




EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "_MatSymTransGetExplicitOperator_SymTrans"
PetscErrorCode _MatSymTransGetExplicitOperator_SymTrans( Mat symA, Mat *A )
{
	Mat_SymTrans   s;
	PetscTruth is_sym;
	
	PetscFunctionBegin;
	
	is_sym = PETSC_FALSE;
	Stg_PetscTypeCompare( (PetscObject)symA, "symtrans", &is_sym );
	if( is_sym == PETSC_FALSE ) {
		*A = symA;
		PetscObjectReference( (PetscObject)symA );
		PetscFunctionReturn(0);
	}
	
	s = (Mat_SymTrans)symA->data;
	MatTranspose( s->A, MAT_INITIAL_MATRIX, A );
	
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "MatSymTransGetExplicitOperator"
PetscErrorCode MatSymTransGetExplicitOperator( Mat symAt, Mat *A )
{
	PetscErrorCode (*f)(Mat,Mat*);
	
	
	PetscObjectQueryFunction( (PetscObject)symAt, "MatSymTransGetExplicitOperator_C", (void (**)(void))&f );
	if( f ) {
		(*f)(symAt, A );
	} else {
		Stg_SETERRQ(PETSC_ERR_ARG_WRONG,"Cannot get explicit operator for this matrix");
	}
	
	PetscFunctionReturn(0);
}





#undef __FUNCT__  
#define __FUNCT__ "MatCreateSymTrans"
PetscErrorCode MatCreateSymTrans( MPI_Comm comm, Mat A, Mat *symAt )
{	
	MPI_Comm commA;
	PetscInt M,N, m,n;
	
	PetscObjectGetComm( (PetscObject)A, &commA );
	if( comm != commA ) {
//		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "Communicators on A must match comm passed in \n" );
	}
	MatGetSize( A, &M, &N );
	MatGetLocalSize( A, &m, &n );
	
	MatCreate( comm, symAt );
	MatSetSizes( *symAt, n,m, N,M );
	MatSetType( *symAt, "symtrans" );
        MatSetUp( *symAt );
	MatSymTransSetOperator( *symAt, A );
	
	
	PetscFunctionReturn(0);
}


/* Constructors */
#undef __FUNCT__  
#define __FUNCT__ "MatSetOps_SymTrans"
PetscErrorCode MatSetOps_SymTrans( struct _MatOps* ops )
{
	PetscFunctionBegin;

	/* 0*/
	ops->setvalues  = 0; // MatSetValues_EMPTY;
	ops->getrow     = 0; // MatGetRow_EMPTY;
	ops->restorerow = 0; // MatRestoreRow_EMPTY;
	ops->mult       = MatMult_SymTrans;
	ops->multadd    = MatMultAdd_SymTrans;
	/* 5*/
	ops->multtranspose    = MatMultTranspose_SymTrans;
	ops->multtransposeadd = MatMultTransposeAdd_SymTrans;
	ops->solve            = MatSolve_SymTrans; // MatSolve_EMPTY;
	ops->solveadd         = MatSolveAdd_SymTrans; // MatSolveAdd_EMPTY;
	ops->solvetranspose   = MatSolveTranspose_SymTrans; // MatSolveTranspose_EMPTY;
	/*10*/
	ops->solvetransposeadd = MatSolveTransposeAdd_SymTrans; // MatSolveTransposeAdd_EMPTY;
	ops->lufactor          = 0; // MatLUFactor_EMPTY;
	ops->choleskyfactor    = 0; // MatCholeskyFactor_EMPTY;
	ops->sor               = 0; // MatSOR_EMPTY;
	ops->transpose         = MatTranspose_SymTrans; // MatTranspose_EMPTY;
	/*15*/
	ops->getinfo       = MatGetInfo_SymTrans; // MatGetInfo_EMPTY;
	ops->equal         = 0; // MatEqual_EMPTY;
	ops->getdiagonal   = 0; // MatGetDiagonal_EMPTY;
	ops->diagonalscale = 0; // MatDiagonalScale_EMPTY;
	ops->norm          = MatNorm_SymTrans; // MatNorm_EMPTY;
	/*20*/
	ops->assemblybegin = 0;
	ops->assemblyend   = 0;
	ops->setoption     = 0; // MatSetOption_EMPTY;
	ops->zeroentries   = 0; // MatZeroEntries_EMPTY;
	/*24*/
	ops->zerorows               = 0; // MatZeroRows_EMPTY;
	ops->lufactorsymbolic       = 0; // MatLUFactorSymbolic_EMPTY;
	ops->lufactornumeric        = 0; // MatLUFactorNumeric_EMPTY;
	ops->choleskyfactorsymbolic = 0; // MatCholeskyFactorSymbolic_EMPTY;
	ops->choleskyfactornumeric  = 0; // MatCholeskyFactorNumeric_EMPTY;
	/*29*/

	ops->ilufactorsymbolic  = 0; // MatILUFactorSymbolic_EMPTY;
	ops->iccfactorsymbolic  = 0; // MatICCFactorSymbolic_EMPTY;


	/*35*/
	ops->duplicate     = 0; // MatDuplicate_EMPTY;
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
	//ops->dummy       = 0; // MatILUDTFactor_EMPTY;
	/*49*/

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
	ops->destroy       = Stg_MatDestroy_SymTrans;
	ops->view          = 0;
	ops->convertfrom   = 0; // MatConvertFrom_EMPTY;

	/*64*/


	ops->setlocaltoglobalmapping = 0; // MatSetLocalToGlobalMapping_EMPTY;
	ops->setvalueslocal          = 0; // MatSetValuesLocal_EMPTY;
	ops->zerorowslocal           = 0; // MatZeroRowsLocal_EMPTY;
	/*69*/
	ops->getrowmaxabs    = 0; // MatGetRowMaxAbs_EMPTY;
	ops->getrowminabs    = 0; // 
	ops->convert         = 0; // MatConvert_EMPTY;
	ops->setcoloring     = 0; // MatSetColoring_EMPTY;

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
	//ops->dummy4                  = 0; // MatPBRelax_EMPTY;
	ops->getvecs                 = MatGetVecs_SymTrans;
	/*89*/
	ops->matmult         = MatMatMult_SymTrans; // MatMatMult_EMPTY;
	ops->matmultsymbolic = 0; // MatMatMultSymbolic_EMPTY;
	ops->matmultnumeric  = 0; // MatMatMultNumeric_EMPTY;
	ops->ptap            = MatPtAP_SymTrans; // MatPtAP_EMPTY;
	ops->ptapsymbolic    = 0; // MatPtAPSymbolic_EMPTY;
	/*94*/
	ops->ptapnumeric              = 0; // MatPtAPNumeric_EMPTY;




	/*99*/



	ops->conjugate           = 0; // MatConjugate_EMPTY;

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
	////ops->dummy2                 = 0;
	////ops->dummy3                 = 0;
	/* 119 */
	ops->multdiagonalblock      = 0;
	ops->hermitiantranspose = 0;
	ops->multhermitiantranspose = 0;
	ops->multhermitiantransposeadd = 0;
	
	PetscFunctionReturn(0);
}




/*
The following operations are not to used with MatSymTrans
 - MatGetOwnershipRange(A,&s,&e):  Result, s=0,e=0.
 - MatGetOwnershripRanges(A,&r):   Result, error.
*/
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "MatCreate_SymTrans"
PetscErrorCode MatCreate_SymTrans( Mat A )
{
	Mat_SymTrans    s;
	PetscErrorCode  ierr;
	
	
	PetscFunctionBegin;
	
	/* define the operations */
//	ierr = PetscMemcpy( A->ops, &MatSymTransOps_Values,sizeof(struct _MatOps) );CHKERRQ(ierr);
	MatSetOps_SymTrans( A->ops );
	
//	A->precision       = PETSC_SCALAR;
	//A->factor           = MAT_FACTOR_NONE;
	//A->mapping          = 0;
	A->spptr            = 0;
	A->same_nonzero     = PETSC_FALSE;
	A->assembled        = PETSC_FALSE;
	
	/* allocate and set pointer for implememtation data */
	ierr = PetscMalloc( sizeof(struct _Mat_SymTrans), &s );CHKERRQ(ierr);
	ierr = PetscMemzero( s, sizeof(struct _Mat_SymTrans) );CHKERRQ(ierr);
	A->data            = (void*)s;
	
	/* Init the implementation data */
	s->A = PETSC_NULL;
	
	/* Set type */
	PetscObjectChangeTypeName((PetscObject) A, "symtrans" );
	
	
	/* Compose functions */
	PetscObjectComposeFunctionDynamic( (PetscObject)A, "MatSymTransSetOperator_C",
			"_MatSymTransSetOperator_SymTrans", _MatSymTransSetOperator_SymTrans );
	
	PetscObjectComposeFunctionDynamic( (PetscObject)A, "MatSymTransGetOperator_C",
			"_MatSymTransGetOperator_SymTrans", _MatSymTransGetOperator_SymTrans );
	
	PetscObjectComposeFunctionDynamic( (PetscObject)A, "MatSymTransGetExplicitOperator_C",
			"_MatSymTransGetExplicitOperator_SymTrans", _MatSymTransGetExplicitOperator_SymTrans );
	
	
	
	PetscFunctionReturn(0);
}
EXTERN_C_END
