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
#include <petscversion.h>
#include <petscvec.h>
#include <petscmat.h>
#if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 3) )
  #include <petsc-private/matimpl.h>
#else
  #include <private/matimpl.h>
#endif

#include "mat-symbolic_transpose-impl.h"
#include "private/mat/mat-symbolic_transpose.h"

#include "petscext_helpers.h"

/* 0*/
#undef __FUNCT__
#define __FUNCT__ "MatSetValues_SymTrans"
PetscErrorCode MatSetValues_SymTrans(Mat mat,PetscInt m,const PetscInt idxm[],PetscInt n,const PetscInt idxn[],const PetscScalar v[],InsertMode addv)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );

	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatGetRow_SymTrans"
PetscErrorCode MatGetRow_SymTrans(Mat mat,PetscInt row,PetscInt *ncols,PetscInt *cols[],PetscScalar *vals[])
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatRestoreRow_SymTrans"
PetscErrorCode MatRestoreRow_SymTrans(Mat mat,PetscInt row,PetscInt *ncols,PetscInt *cols[],PetscScalar *vals[])
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatMult_SymTrans"
PetscErrorCode MatMult_SymTrans(Mat mat,Vec x,Vec y)
{
	Mat_SymTrans   s = (Mat_SymTrans)mat->data;	
    PetscErrorCode ierr;

	PetscFunctionBegin; 

    ierr = MatMultTranspose( s->A, x, y );CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatMultAdd_SymTrans"
PetscErrorCode MatMultAdd_SymTrans(Mat mat,Vec v1,Vec v2,Vec v3)
{
	Mat_SymTrans   s = (Mat_SymTrans)mat->data;	

	PetscFunctionBegin; 

	MatMultTransposeAdd( s->A, v1, v2, v3 );
	
	PetscFunctionReturn(0);
}


/* 5*/
#undef __FUNCT__
#define __FUNCT__ "MatMultTranspose_SymTrans"
PetscErrorCode MatMultTranspose_SymTrans(Mat mat,Vec x,Vec y)
{
	
	PetscFunctionBegin;
	
	Mat_SymTrans   s = (Mat_SymTrans)mat->data;
	MatMult( s->A, x, y );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatMultTransposeAdd_SymTrans"
PetscErrorCode MatMultTransposeAdd_SymTrans(Mat mat,Vec v1,Vec v2,Vec v3)
{
	Mat_SymTrans   s = (Mat_SymTrans)mat->data;
	
	PetscFunctionBegin; 
	
	MatMultAdd( s->A, v1, v2, v3 );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSolve_SymTrans"
PetscErrorCode MatSolve_SymTrans(Mat mat,Vec b,Vec x)
{
	Mat_SymTrans   s = (Mat_SymTrans)mat->data;
	
	PetscFunctionBegin; 
	
	MatSolveTranspose( s->A, b, x );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSolveAdd_SymTrans"
PetscErrorCode MatSolveAdd_SymTrans(Mat mat,Vec b,Vec y,Vec x)
{
	Mat_SymTrans   s = (Mat_SymTrans)mat->data;
	
	PetscFunctionBegin; 
	
	MatSolveTransposeAdd( s->A, b, y, x );
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatSolveTranspose_SymTrans"
PetscErrorCode MatSolveTranspose_SymTrans(Mat mat,Vec b,Vec x)
{
	Mat_SymTrans   s = (Mat_SymTrans)mat->data;
	
	PetscFunctionBegin; 
	
	MatSolve( s->A, b, x );
	
	PetscFunctionReturn(0);
}

/*10*/
#undef __FUNCT__
#define __FUNCT__ "MatSolveTransposeAdd_SymTrans"
PetscErrorCode MatSolveTransposeAdd_SymTrans(Mat mat,Vec b,Vec y,Vec x)
{
	Mat_SymTrans   s = (Mat_SymTrans)mat->data;
	
	PetscFunctionBegin; 
	
	MatSolveAdd( s->A, b, y, x );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatLUFactor_SymTrans"
PetscErrorCode MatLUFactor_SymTrans(Mat mat,IS row,IS col,MatFactorInfo *info)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatCholeskyFactor_SymTrans"
PetscErrorCode MatCholeskyFactor_SymTrans(Mat mat,IS perm,MatFactorInfo *info)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatRelax_SymTrans"
PetscErrorCode MatRelax_SymTrans(Mat mat,Vec b,PetscReal omega,MatSORType flag,PetscReal shift,PetscInt its,PetscInt lits,Vec x)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatTranspose_SymTrans"
PetscErrorCode MatTranspose_SymTrans(Mat mat,MatReuse reuse,Mat *B)
{
	Mat_SymTrans   s = (Mat_SymTrans)mat->data;
	
	PetscFunctionBegin; 
	
	*B = s->A;
	PetscObjectReference( (PetscObject)s->A );
	
	PetscFunctionReturn(0);
}

/*15*/
#undef __FUNCT__
#define __FUNCT__ "MatGetInfo_SymTrans"
PetscErrorCode MatGetInfo_SymTrans(Mat mat,MatInfoType flag,MatInfo *info)
{
	PetscErrorCode ierr;
	Mat_SymTrans   s = (Mat_SymTrans)mat->data;
	
	PetscFunctionBegin; 

	/* pull info out of the original operator */	
	ierr = MatGetInfo(s->A,flag,info);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatEqual_SymTrans"
PetscErrorCode MatEqual_SymTrans(Mat A,Mat B,PetscTruth *flg)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatGetDiagonal_SymTrans"
PetscErrorCode MatGetDiagonal_SymTrans(Mat mat,Vec v)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatDiagonalScale_SymTrans"
PetscErrorCode MatDiagonalScale_SymTrans(Mat mat,Vec l,Vec r)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatNorm_SymTrans"
PetscErrorCode MatNorm_SymTrans(Mat mat,NormType type,PetscReal *nrm)
{
	Mat_SymTrans   s = (Mat_SymTrans)mat->data;
	PetscErrorCode ierr;
	
	PetscFunctionBegin; 
	
	ierr = MatNorm(s->A,type,nrm);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/*20*/
#undef __FUNCT__
#define __FUNCT__ "MatAssemblyBegin_SymTrans"
PetscErrorCode MatAssemblyBegin_SymTrans(Mat mat,MatAssemblyType type)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatAssemblyEnd_SymTrans"
PetscErrorCode MatAssemblyEnd_SymTrans(Mat mat,MatAssemblyType type)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatCompress_SymTrans"
PetscErrorCode MatCompress_SymTrans(Mat mat)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSetOption_SymTrans"
PetscErrorCode MatSetOption_SymTrans(Mat mat,MatOption op)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatZeroEntries_SymTrans"
PetscErrorCode MatZeroEntries_SymTrans(Mat mat)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*25*/


#undef __FUNCT__
#define __FUNCT__ "MatZeroRows_SymTrans"
PetscErrorCode MatZeroRows_SymTrans(Mat mat,PetscInt numRows,const PetscInt rows[],PetscScalar diag)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatLUFactorSymbolic_SymTrans"
PetscErrorCode MatLUFactorSymbolic_SymTrans(Mat mat,IS row,IS col,MatFactorInfo *info,Mat *fact)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatLUFactorNumeric_SymTrans"
PetscErrorCode MatLUFactorNumeric_SymTrans(Mat mat,MatFactorInfo *info,Mat *fact)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatCholeskyFactorSymbolic_SymTrans"
PetscErrorCode MatCholeskyFactorSymbolic_SymTrans(Mat mat,IS perm,MatFactorInfo *info,Mat *fact)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatCholeskyFactorNumeric_SymTrans"
PetscErrorCode MatCholeskyFactorNumeric_SymTrans(Mat mat,MatFactorInfo *info,Mat *fact)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*30*/
#undef __FUNCT__
#define __FUNCT__ "MatSetUpPreallocation_SymTrans"
PetscErrorCode MatSetUpPreallocation_SymTrans(Mat B)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatILUFactorSymbolic_SymTrans"
PetscErrorCode MatILUFactorSymbolic_SymTrans(Mat mat,IS row,IS col,MatFactorInfo *info,Mat *fact)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatICCFactorSymbolic_SymTrans"
PetscErrorCode MatICCFactorSymbolic_SymTrans(Mat mat,IS perm,MatFactorInfo *info,Mat *fact)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatGetArray_SymTrans"
PetscErrorCode MatGetArray_SymTrans(Mat mat,PetscScalar *v[])
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatRestoreArray_SymTrans"
PetscErrorCode MatRestoreArray_SymTrans(Mat mat,PetscScalar *v[])
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*35*/
#undef __FUNCT__
#define __FUNCT__ "MatDuplicate_SymTrans"
PetscErrorCode MatDuplicate_SymTrans(Mat mat,MatDuplicateOption op,Mat *M)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatForwardSolve_SymTrans"
PetscErrorCode MatForwardSolve_SymTrans(Mat mat,Vec b,Vec x)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatBackwardSolve_SymTrans"
PetscErrorCode MatBackwardSolve_SymTrans(Mat mat,Vec b,Vec x)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatILUFactor_SymTrans"
PetscErrorCode MatILUFactor_SymTrans(Mat mat,IS row,IS col,MatFactorInfo *info)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatICCFactor_SymTrans"
PetscErrorCode MatICCFactor_SymTrans(Mat mat,IS row,MatFactorInfo* info)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*40*/
#undef __FUNCT__
#define __FUNCT__ "MatAXPY_SymTrans"
PetscErrorCode MatAXPY_SymTrans(Mat Y,PetscScalar a,Mat X,MatStructure str)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatGetSubMatrices_SymTrans"
PetscErrorCode MatGetSubMatrices_SymTrans(Mat mat,PetscInt n,const IS irow[],const IS icol[],MatReuse scall,Mat *submat[])
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatIncreaseOverlap_SymTrans"
PetscErrorCode MatIncreaseOverlap_SymTrans(Mat mat,PetscInt n,IS is[],PetscInt ov)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatGetValues_SymTrans"
PetscErrorCode MatGetValues_SymTrans(Mat mat,PetscInt m,const PetscInt idxm[],PetscInt n,const PetscInt idxn[],PetscScalar v[])
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatCopy_SymTrans"
PetscErrorCode MatCopy_SymTrans(Mat A,Mat B,MatStructure str)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*45*/
#undef __FUNCT__
#define __FUNCT__ "MatGetRowMax_SymTrans"
PetscErrorCode MatGetRowMax_SymTrans(Mat mat,Vec v,PetscInt idx[])
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatScale_SymTrans"
PetscErrorCode MatScale_SymTrans(Mat mat,PetscScalar a)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatShift_SymTrans"
PetscErrorCode MatShift_SymTrans( Mat mat, PetscScalar a )
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatDiagonalSet_SymTrans"
PetscErrorCode MatDiagonalSet_SymTrans( Mat Y, Vec D, InsertMode is )
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatILUDTFactor_SymTrans"
PetscErrorCode MatILUDTFactor_SymTrans(Mat mat,IS row,IS col,MatFactorInfo *info,Mat *fact)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*50*/
#undef __FUNCT__
#define __FUNCT__ "MatSetBlockSize_SymTrans"
PetscErrorCode MatSetBlockSize_SymTrans(Mat mat,PetscInt bs)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatGetRowIJ_SymTrans"
PetscErrorCode MatGetRowIJ_SymTrans(Mat mat,PetscInt shift,PetscTruth symmetric,PetscTruth blockcompressed,PetscInt *n,PetscInt *ia[],PetscInt* ja[],PetscTruth *done)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatRestorRowIJ_SymTrans"
PetscErrorCode MatRestorRowIJ_SymTrans(Mat mat,PetscInt shift,PetscTruth symmetric,PetscTruth blockcompressed,PetscInt *n,PetscInt *ia[],PetscInt* ja[],PetscTruth *done)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatGetColumnIJ_SymTrans"
PetscErrorCode MatGetColumnIJ_SymTrans(Mat mat,PetscInt shift,PetscTruth symmetric,PetscTruth blockcompressed,PetscInt *n,PetscInt *ia[],PetscInt* ja[],PetscTruth *done)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatRestorColumnIJ_SymTrans"
PetscErrorCode MatRestorColumnIJ_SymTrans(Mat mat,PetscInt shift,PetscTruth symmetric,PetscTruth blockcompressed,PetscInt *n,PetscInt *ia[],PetscInt* ja[],PetscTruth *done)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*55*/

#undef __FUNCT__
#define __FUNCT__ "MatFDColoringCreate_SymTrans"
PetscErrorCode MatFDColoringCreate_SymTrans(Mat mat,ISColoring iscoloring,MatFDColoring color)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatColoringPatch_SymTrans"
PetscErrorCode MatColoringPatch_SymTrans(Mat mat,PetscInt ncolors,PetscInt n,ISColoringValue colorarray[],ISColoring *iscoloring)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSetUnfactored_SymTrans"
PetscErrorCode MatSetUnfactored_SymTrans(Mat mat)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatPermute_SymTrans"
PetscErrorCode MatPermute_SymTrans( Mat mat, IS row, IS col, Mat *B )
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSetValuesBlocked_SymTrans"
PetscErrorCode MatSetValuesBlocked_SymTrans(Mat mat,PetscInt m,const PetscInt idxm[],PetscInt n,const PetscInt idxn[],const PetscScalar v[],InsertMode addv)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*60*/
#undef __FUNCT__
#define __FUNCT__ "MatGetSubMatrix_SymTrans"
PetscErrorCode MatGetSubMatrix_SymTrans(Mat mat,IS isrow,IS iscol,PetscInt csize,MatReuse cll,Mat *newmat)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
/*
#undef __FUNCT__
#define __FUNCT__ "Stg_MatDestroy_SymTrans"
PetscErrorCode Stg_MatDestroy_SymTrans( Mat mat )
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
*/
#undef __FUNCT__
#define __FUNCT__ "MatView_SymTrans"
PetscErrorCode MatView_SymTrans( Mat mat, PetscViewer viewer )
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatConvertFrom_SymTrans"
PetscErrorCode MatConvertFrom_SymTrans(Mat mat, MatType newtype,MatReuse reuse,Mat *M)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatUseScaledForm_SymTrans"
PetscErrorCode MatUseScaledForm_SymTrans( Mat mat, PetscTruth scaled )
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*65*/
#undef __FUNCT__
#define __FUNCT__ "MatScaleSystem_SymTrans"
PetscErrorCode MatScaleSystem_SymTrans( Mat mat, Vec b, Vec x )
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatUnScaleSystem_SymTrans"
PetscErrorCode MatUnScaleSystem_SymTrans(Mat mat,Vec b,Vec x)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSetLocalToGlobalMapping_SymTrans"
PetscErrorCode MatSetLocalToGlobalMapping_SymTrans(Mat x,ISLocalToGlobalMapping mapping)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSetValuesLocal_SymTrans"
PetscErrorCode MatSetValuesLocal_SymTrans(Mat mat,PetscInt nrow,const PetscInt irow[],PetscInt ncol,const PetscInt icol[],const PetscScalar y[],InsertMode addv)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatZeroRowsLocal_SymTrans"
PetscErrorCode MatZeroRowsLocal_SymTrans(Mat mat,PetscInt numRows,const PetscInt rows[],PetscScalar diag)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

/*70*/
#undef __FUNCT__
#define __FUNCT__ "MatGetRowMaxAbs_SymTrans"
PetscErrorCode MatGetRowMaxAbs_SymTrans(Mat mat,Vec v,PetscInt idx[])
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatConvert_SymTrans"
PetscErrorCode MatConvert_SymTrans(Mat mat, MatType newtype,MatReuse reuse,Mat *M)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSetColoring_SymTrans"
PetscErrorCode MatSetColoring_SymTrans( Mat mat, ISColoring colour )
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSetValuesAdic_SymTrans"
PetscErrorCode MatSetValuesAdic_SymTrans(Mat mat,void *v)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSetValuesAdifor_SymTrans"
PetscErrorCode MatSetValuesAdifor_SymTrans(Mat mat,PetscInt nl,void *v)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}


/*75*/
#undef __FUNCT__
#define __FUNCT__ "MatFDColoringApply_SymTrans"
PetscErrorCode MatFDColoringApply_SymTrans(Mat J,MatFDColoring coloring,Vec x1,MatStructure *flag,void *sctx)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSetFromOptions_SymTrans"
PetscErrorCode MatSetFromOptions_SymTrans( Mat mat )
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatMultConstrained_SymTrans"
PetscErrorCode MatMultConstrained_SymTrans(Mat mat,Vec x,Vec y)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatMultTransposeConstrained_SymTrans"
PetscErrorCode MatMultTransposeConstrained_SymTrans(Mat mat,Vec x,Vec y)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatILUFactorSymbolicConstrained_SymTrans"
PetscErrorCode MatILUFactorSymbolicConstrained_SymTrans(Mat mat,IS row,IS col,double d,PetscInt i,PetscInt j,Mat *fact)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}


/*80*/
#undef __FUNCT__
#define __FUNCT__ "MatPermuteSparsify_SymTrans"
PetscErrorCode MatPermuteSparsify_SymTrans(Mat A, PetscInt band, PetscReal frac, PetscReal tol, IS rowp, IS colp, Mat *B)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatMults_SymTrans"
PetscErrorCode MatMults_SymTrans(Mat A,Vecs x,Vecs y)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSolves_SymTrans"
PetscErrorCode MatSolves_SymTrans(Mat A,Vecs b,Vecs x)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatGetInertia_SymTrans"
PetscErrorCode MatGetInertia_SymTrans(Mat mat,PetscInt *nneg,PetscInt *nzero,PetscInt *npos)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "Stg_MatLoad_SymTrans"
PetscErrorCode Stg_MatLoad_SymTrans(PetscViewer viewer, MatType outtype,Mat *newmat)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}


/*85*/
#undef __FUNCT__
#define __FUNCT__ "MatIsSymmetric_SymTrans"
PetscErrorCode MatIsSymmetric_SymTrans(Mat A,PetscReal tol,PetscTruth *flg)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatIsHermitian_SymTrans"
PetscErrorCode MatIsHermitian_SymTrans(Mat A,PetscTruth *flg)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatIsStructurallySymmetric_SymTrans"
PetscErrorCode MatIsStructurallySymmetric_SymTrans(Mat A,PetscTruth *flg)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatPBRelax_SymTrans"
PetscErrorCode MatPBRelax_SymTrans(Mat mat,Vec b,PetscReal omega,MatSORType flag,PetscReal shift,PetscInt its,PetscInt lits,Vec x)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatGetVecs_SymTrans"
PetscErrorCode MatGetVecs_SymTrans(Mat mat,Vec *right,Vec *left)
{
	Mat_SymTrans   s = (Mat_SymTrans)mat->data;
	
	PetscFunctionBegin; 
	
	MatGetVecs( s->A, left, right );
	
	PetscFunctionReturn(0);
}



/*90*/
#undef __FUNCT__
#define __FUNCT__ "MatMatMult_SymTrans"
PetscErrorCode MatMatMult_SymTrans(Mat A,Mat B,MatReuse scall,PetscReal fill,Mat *C) 
{
	Mat_SymTrans   sA;
	Mat _A, _B;
	PetscTruth is_Asym;
	MPI_Comm comm;
	PetscInt M,N;	

	PetscFunctionBegin; 
	
	/* Special cases, 
		i)  A=symBt => use Pt I P, where P = B
	*/
	
	is_Asym = PETSC_FALSE;
	Stg_PetscTypeCompare( (PetscObject)A, "symtrans", &is_Asym );
	if( is_Asym == PETSC_TRUE ) {
		sA = (Mat_SymTrans)A->data;
		
		if( sA->A == B ) {
			/* Build identity */
			Mat Identity;
			
			PetscObjectGetComm( (PetscObject)B, &comm );
			MatGetSize( B, &M, &N );
			Identity = MATAIJIdentityNew( comm, M ); 
			MatPtAP( Identity, B, scall, fill, C );
			Stg_MatDestroy( &Identity );
			
			PetscFunctionReturn(0);
		}
	}
	
	MatSymTransGetExplicitOperator( A, &_A );
	MatSymTransGetExplicitOperator( B, &_B );
	
	MatMatMult( _A, _B, scall, fill, C );
	
	Stg_MatDestroy( &_B );
	Stg_MatDestroy( &_A );
	
	
	
	PetscFunctionReturn(0);
	
}
#undef __FUNCT__
#define __FUNCT__ "MatMatMultSymbolic_SymTrans"
PetscErrorCode MatMatMultSymbolic_SymTrans(Mat A,Mat B,PetscReal fill,Mat *C) 
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatMatMultNumeric_SymTrans"
PetscErrorCode MatMatMultNumeric_SymTrans(Mat A,Mat B,Mat C)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatPtAP_SymTrans"
PetscErrorCode MatPtAP_SymTrans(Mat A,Mat P,MatReuse scall,PetscReal fill,Mat *C) 
{
	Mat _A, _P;
	
	PetscFunctionBegin; 
	
	
	MatSymTransGetExplicitOperator( A, &_A );
	MatSymTransGetExplicitOperator( P, &_P );
	
	MatPtAP( _A, _P, scall, fill, C );
	
	Stg_MatDestroy( &_P );
	Stg_MatDestroy( &_A );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatPtAPSymbolic_SymTrans"
PetscErrorCode MatPtAPSymbolic_SymTrans(Mat A,Mat P,PetscReal fill,Mat *C) 
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}


/*95*/
#undef __FUNCT__
#define __FUNCT__ "MatPtAPNumeric_SymTrans"
PetscErrorCode MatPtAPNumeric_SymTrans(Mat A,Mat P,Mat C) 
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatMatMultTranspose_SymTrans"
PetscErrorCode MatMatMultTranspose_SymTrans(Mat A,Mat B,MatReuse scall,PetscReal fill,Mat *C) 
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatMatMultTransposeSymbolic_SymTrans"
PetscErrorCode MatMatMultTransposeSymbolic_SymTrans(Mat A,Mat B,PetscReal fill,Mat *C) 
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatMatMultTransposeNumeric_SymTrans"
PetscErrorCode MatMatMultTransposeNumeric_SymTrans(Mat A,Mat B,Mat C)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatPtAPSymbolic__SEQAIJ_SymTrans"
PetscErrorCode MatPtAPSymbolic__SEQAIJ_SymTrans(Mat A,Mat P,PetscReal fill,Mat *C) 
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}


/*100*/
#undef __FUNCT__
#define __FUNCT__ "MatPtAPNumeric_SEQAIJ_SymTrans"
PetscErrorCode MatPtAPNumeric_SEQAIJ_SymTrans(Mat A,Mat P,Mat C) 
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatPtAPSymbolic_MPIAIJ_SymTrans"
PetscErrorCode MatPtAPSymbolic_MPIAIJ_SymTrans(Mat A,Mat P,PetscReal fill,Mat *C) 
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatPtAPNumeric_MPIAIJ_SymTrans"
PetscErrorCode MatPtAPNumeric_MPIAIJ_SymTrans(Mat A,Mat P,Mat C) 
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatConjugate_SymTrans"
PetscErrorCode MatConjugate_SymTrans(Mat mat)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatSetSizes_SymTrans"
PetscErrorCode MatSetSizes_SymTrans(Mat mat,PetscInt m,PetscInt n,PetscInt M,PetscInt N)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}


/*105*/
#undef __FUNCT__
#define __FUNCT__ "MatSetValuesRow_SymTrans"
PetscErrorCode MatSetValuesRow_SymTrans(Mat mat,PetscInt row,const PetscScalar v[])
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatRealPart_SymTrans"
PetscErrorCode MatRealPart_SymTrans(Mat mat)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatImaginaryPart_SymTrans"
PetscErrorCode MatImaginaryPart_SymTrans(Mat mat)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatGetRowUpperTriangular_SymTrans"
PetscErrorCode MatGetRowUpperTriangular_SymTrans(Mat mat)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatRestroreRowUpperTriangular_SymTrans"
PetscErrorCode MatRestroreRowUpperTriangular_SymTrans(Mat mat)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}


/*110*/
#undef __FUNCT__
#define __FUNCT__ "MatMatSolve_SymTrans"
PetscErrorCode MatMatSolve_SymTrans(Mat mat,Mat b,Mat x)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatGetRedundantMatrix_SymTrans"
PetscErrorCode MatGetRedundantMatrix_SymTrans(Mat mat,PetscInt nsubcomm,MPI_Comm subcomm,PetscInt mlocal_red,MatReuse reuse,Mat *matredundant)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatGetRowMin_SymTrans"
PetscErrorCode MatGetRowMin_SymTrans(Mat mat,Vec v,PetscInt idx[])
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MatGetColumnVector_SymTrans"
PetscErrorCode MatGetColumnVector_SymTrans(Mat A,Vec yy,PetscInt col)
{
	PetscFunctionBegin; Stg_SETERRQ1( PETSC_ERR_SUP,"Operation: %s is not implemented", __func__ );
	
	PetscFunctionReturn(0);
}

