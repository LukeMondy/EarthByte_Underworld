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


#ifndef __PETSC_MAT_OPS_SYM_TRANS_H__
#define __PETSC_MAT_OPS_SYM_TRANS_H__
#include <petscversion.h>
#include "petscext_helpers.h"

/* 0*/
PetscErrorCode MatSetValues_SymTrans(Mat mat,PetscInt m,const PetscInt idxm[],PetscInt n,const PetscInt idxn[],const PetscScalar v[],InsertMode addv);
PetscErrorCode MatGetRow_SymTrans(Mat mat,PetscInt row,PetscInt *ncols,PetscInt *cols[],PetscScalar *vals[]);
PetscErrorCode MatRestoreRow_SymTrans(Mat mat,PetscInt row,PetscInt *ncols,PetscInt *cols[],PetscScalar *vals[]);
PetscErrorCode MatMult_SymTrans(Mat mat,Vec x,Vec y);
PetscErrorCode MatMultAdd_SymTrans(Mat mat,Vec v1,Vec v2,Vec v3);
/* 5*/
PetscErrorCode MatMultTranspose_SymTrans(Mat mat,Vec x,Vec y);
PetscErrorCode MatMultTransposeAdd_SymTrans(Mat mat,Vec v1,Vec v2,Vec v3);
PetscErrorCode MatSolve_SymTrans(Mat mat,Vec b,Vec x);
PetscErrorCode MatSolveAdd_SymTrans(Mat mat,Vec b,Vec y,Vec x);
PetscErrorCode MatSolveTranspose_SymTrans(Mat mat,Vec b,Vec x);
/*10*/
PetscErrorCode MatSolveTransposeAdd_SymTrans(Mat mat,Vec b,Vec y,Vec x);
PetscErrorCode MatLUFactor_SymTrans(Mat mat,IS row,IS col,MatFactorInfo *info);
PetscErrorCode MatCholeskyFactor_SymTrans(Mat mat,IS perm,MatFactorInfo *info);
PetscErrorCode MatRelax_SymTrans(Mat mat,Vec b,PetscReal omega,MatSORType flag,PetscReal shift,PetscInt its,PetscInt lits,Vec x);
PetscErrorCode MatTranspose_SymTrans(Mat mat,MatReuse reuse,Mat *B);
/*15*/
PetscErrorCode MatGetInfo_SymTrans(Mat mat,MatInfoType flag,MatInfo *info);
PetscErrorCode MatEqual_SymTrans(Mat A,Mat B,PetscTruth *flg);
PetscErrorCode MatGetDiagonal_SymTrans(Mat mat,Vec v);
PetscErrorCode MatDiagonalScale_SymTrans(Mat mat,Vec l,Vec r);
PetscErrorCode MatNorm_SymTrans(Mat mat,NormType type,PetscReal *nrm);
/*20*/
PetscErrorCode MatAssemblyBegin_SymTrans(Mat mat,MatAssemblyType type);
PetscErrorCode MatAssemblyEnd_SymTrans(Mat mat,MatAssemblyType type);
PetscErrorCode MatCompress_SymTrans(Mat mat);
PetscErrorCode MatSetOption_SymTrans(Mat mat,MatOption op);
PetscErrorCode MatZeroEntries_SymTrans(Mat mat);
/*25*/
PetscErrorCode MatZeroRows_SymTrans(Mat mat,PetscInt numRows,const PetscInt rows[],PetscScalar diag);
PetscErrorCode MatLUFactorSymbolic_SymTrans(Mat mat,IS row,IS col,MatFactorInfo *info,Mat *fact);
PetscErrorCode MatLUFactorNumeric_SymTrans(Mat mat,MatFactorInfo *info,Mat *fact);
PetscErrorCode MatCholeskyFactorSymbolic_SymTrans(Mat mat,IS perm,MatFactorInfo *info,Mat *fact);
PetscErrorCode MatCholeskyFactorNumeric_SymTrans(Mat mat,MatFactorInfo *info,Mat *fact);
/*30*/
PetscErrorCode MatSetUpPreallocation_SymTrans(Mat B);
PetscErrorCode MatILUFactorSymbolic_SymTrans(Mat mat,IS row,IS col,MatFactorInfo *info,Mat *fact);
PetscErrorCode MatICCFactorSymbolic_SymTrans(Mat mat,IS perm,MatFactorInfo *info,Mat *fact);
PetscErrorCode MatGetArray_SymTrans(Mat mat,PetscScalar *v[]);
PetscErrorCode MatRestoreArray_SymTrans(Mat mat,PetscScalar *v[]);
/*35*/
PetscErrorCode MatDuplicate_SymTrans(Mat mat,MatDuplicateOption op,Mat *M);
PetscErrorCode MatForwardSolve_SymTrans(Mat mat,Vec b,Vec x);
PetscErrorCode MatBackwardSolve_SymTrans(Mat mat,Vec b,Vec x);
PetscErrorCode MatILUFactor_SymTrans(Mat mat,IS row,IS col,MatFactorInfo *info);
PetscErrorCode MatICCFactor_SymTrans(Mat mat,IS row,MatFactorInfo* info);
/*40*/
PetscErrorCode MatAXPY_SymTrans(Mat Y,PetscScalar a,Mat X,MatStructure str);
PetscErrorCode MatGetSubMatrices_SymTrans(Mat mat,PetscInt n,const IS irow[],const IS icol[],MatReuse scall,Mat *submat[]);
PetscErrorCode MatIncreaseOverlap_SymTrans(Mat mat,PetscInt n,IS is[],PetscInt ov);
PetscErrorCode MatGetValues_SymTrans(Mat mat,PetscInt m,const PetscInt idxm[],PetscInt n,const PetscInt idxn[],PetscScalar v[]);
PetscErrorCode MatCopy_SymTrans(Mat A,Mat B,MatStructure str);
/*45*/
PetscErrorCode MatGetRowMax_SymTrans(Mat mat,Vec v,PetscInt idx[]);
PetscErrorCode MatScale_SymTrans(Mat mat,PetscScalar a);
PetscErrorCode MatShift_SymTrans( Mat mat, PetscScalar a );
PetscErrorCode MatDiagonalSet_SymTrans( Mat Y, Vec D, InsertMode is );
PetscErrorCode MatILUDTFactor_SymTrans(Mat mat,IS row,IS col,MatFactorInfo *info,Mat *fact);
/*50*/
PetscErrorCode MatSetBlockSize_SymTrans(Mat mat,PetscInt bs);
PetscErrorCode MatGetRowIJ_SymTrans(Mat mat,PetscInt shift,PetscTruth symmetric,PetscTruth blockcompressed,PetscInt *n,PetscInt *ia[],PetscInt* ja[],PetscTruth *done);
PetscErrorCode MatRestorRowIJ_SymTrans(Mat mat,PetscInt shift,PetscTruth symmetric,PetscTruth blockcompressed,PetscInt *n,PetscInt *ia[],PetscInt* ja[],PetscTruth *done);
PetscErrorCode MatGetColumnIJ_SymTrans(Mat mat,PetscInt shift,PetscTruth symmetric,PetscTruth blockcompressed,PetscInt *n,PetscInt *ia[],PetscInt* ja[],PetscTruth *done);
PetscErrorCode MatRestorColumnIJ_SymTrans(Mat mat,PetscInt shift,PetscTruth symmetric,PetscTruth blockcompressed,PetscInt *n,PetscInt *ia[],PetscInt* ja[],PetscTruth *done);
/*55*/
PetscErrorCode MatFDColoringCreate_SymTrans(Mat mat,ISColoring iscoloring,MatFDColoring color);
PetscErrorCode MatColoringPatch_SymTrans(Mat mat,PetscInt ncolors,PetscInt n,ISColoringValue colorarray[],ISColoring *iscoloring);
PetscErrorCode MatSetUnfactored_SymTrans(Mat mat);
PetscErrorCode MatPermute_SymTrans( Mat mat, IS row, IS col, Mat *B );
PetscErrorCode MatSetValuesBlocked_SymTrans(Mat mat,PetscInt m,const PetscInt idxm[],PetscInt n,const PetscInt idxn[],const PetscScalar v[],InsertMode addv);
/*60*/
PetscErrorCode MatGetSubMatrix_SymTrans(Mat mat,IS isrow,IS iscol,PetscInt csize,MatReuse cll,Mat *newmat);
/*PetscErrorCode MatDestroy_SymTrans( Mat mat );*/
PetscErrorCode MatView_SymTrans( Mat mat, PetscViewer viewer );
PetscErrorCode MatConvertFrom_SymTrans(Mat mat, MatType newtype,MatReuse reuse,Mat *M);
PetscErrorCode MatUseScaledForm_SymTrans( Mat mat, PetscTruth scaled );
/*65*/
PetscErrorCode MatScaleSystem_SymTrans( Mat mat, Vec b, Vec x );
PetscErrorCode MatUnScaleSystem_SymTrans(Mat mat,Vec b,Vec x);
PetscErrorCode MatSetLocalToGlobalMapping_SymTrans(Mat x,ISLocalToGlobalMapping mapping);
PetscErrorCode MatSetValuesLocal_SymTrans(Mat mat,PetscInt nrow,const PetscInt irow[],PetscInt ncol,const PetscInt icol[],const PetscScalar y[],InsertMode addv);
PetscErrorCode MatZeroRowsLocal_SymTrans(Mat mat,PetscInt numRows,const PetscInt rows[],PetscScalar diag);
/*70*/
PetscErrorCode MatGetRowMaxAbs_SymTrans(Mat mat,Vec v,PetscInt idx[]);
PetscErrorCode MatConvert_SymTrans(Mat mat, MatType newtype,MatReuse reuse,Mat *M);
PetscErrorCode MatSetColoring_SymTrans( Mat mat, ISColoring colour );
PetscErrorCode MatSetValuesAdic_SymTrans(Mat mat,void *v);
PetscErrorCode MatSetValuesAdifor_SymTrans(Mat mat,PetscInt nl,void *v);
/*75*/
PetscErrorCode MatFDColoringApply_SymTrans(Mat J,MatFDColoring coloring,Vec x1,MatStructure *flag,void *sctx);
PetscErrorCode MatSetFromOptions_SymTrans( Mat mat );
PetscErrorCode MatMultConstrained_SymTrans(Mat mat,Vec x,Vec y);
PetscErrorCode MatMultTransposeConstrained_SymTrans(Mat mat,Vec x,Vec y);
PetscErrorCode MatILUFactorSymbolicConstrained_SymTrans(Mat mat,IS row,IS col,double d,PetscInt i,PetscInt j,Mat *fact);
/*80*/
PetscErrorCode MatPermuteSparsify_SymTrans(Mat A, PetscInt band, PetscReal frac, PetscReal tol, IS rowp, IS colp, Mat *B);
PetscErrorCode MatMults_SymTrans(Mat A,Vecs x,Vecs y);
PetscErrorCode MatSolves_SymTrans(Mat A,Vecs b,Vecs x);
PetscErrorCode MatGetInertia_SymTrans(Mat mat,PetscInt *nneg,PetscInt *nzero,PetscInt *npos);
PetscErrorCode MatLoad_SymTrans(PetscViewer viewer, MatType outtype,Mat *newmat);
/*85*/
PetscErrorCode MatIsSymmetric_SymTrans(Mat A,PetscReal tol,PetscTruth *flg);
PetscErrorCode MatIsHermitian_SymTrans(Mat A,PetscTruth *flg);
PetscErrorCode MatIsStructurallySymmetric_SymTrans(Mat A,PetscTruth *flg);
PetscErrorCode MatPBRelax_SymTrans(Mat mat,Vec b,PetscReal omega,MatSORType flag,PetscReal shift,PetscInt its,PetscInt lits,Vec x);
PetscErrorCode MatGetVecs_SymTrans(Mat mat,Vec *right,Vec *left);
/*90*/
PetscErrorCode MatMatMult_SymTrans(Mat A,Mat B,MatReuse scall,PetscReal fill,Mat *C);
PetscErrorCode MatMatMultSymbolic_SymTrans(Mat A,Mat B,PetscReal fill,Mat *C);
PetscErrorCode MatMatMultNumeric_SymTrans(Mat A,Mat B,Mat C);
PetscErrorCode MatPtAP_SymTrans(Mat A,Mat P,MatReuse scall,PetscReal fill,Mat *C);
PetscErrorCode MatPtAPSymbolic_SymTrans(Mat A,Mat P,PetscReal fill,Mat *C);
/*95*/
PetscErrorCode MatPtAPNumeric_SymTrans(Mat A,Mat P,Mat C);
PetscErrorCode MatMatMultTranspose_SymTrans(Mat A,Mat B,MatReuse scall,PetscReal fill,Mat *C);
PetscErrorCode MatMatMultTransposeSymbolic_SymTrans(Mat A,Mat B,PetscReal fill,Mat *C);
PetscErrorCode MatMatMultTransposeNumeric_SymTrans(Mat A,Mat B,Mat C);
PetscErrorCode MatPtAPSymbolic_SEQAIJ_SymTrans(Mat A,Mat P,PetscReal fill,Mat *C);
/*100*/
PetscErrorCode MatPtAPNumeric_SEQAIJ_SymTrans(Mat A,Mat P,Mat C);
PetscErrorCode MatPtAPSymbolic_MPIAIJ_SymTrans(Mat A,Mat P,PetscReal fill,Mat *C);
PetscErrorCode MatPtAPNumeric_MPIAIJ_SymTrans(Mat A,Mat P,Mat C);
PetscErrorCode MatConjugate_SymTrans(Mat mat);
PetscErrorCode MatSetSizes_SymTrans(Mat mat,PetscInt m,PetscInt n,PetscInt M,PetscInt N);
/*105*/
PetscErrorCode MatSetValuesRow_SymTrans(Mat mat,PetscInt row,const PetscScalar v[]);
PetscErrorCode MatRealPart_SymTrans(Mat mat);
PetscErrorCode MatImaginaryPart_SymTrans(Mat mat);
PetscErrorCode MatGetRowUpperTriangular_SymTrans(Mat mat);
PetscErrorCode MatRestroreRowUpperTriangular_SymTrans(Mat mat);
/*110*/
PetscErrorCode MatMatSolve_SymTrans(Mat mat,Mat b,Mat x);
PetscErrorCode MatGetRedundantMatrix_SymTrans(Mat mat,PetscInt nsubcomm,MPI_Comm subcomm,PetscInt mlocal_red,MatReuse reuse,Mat *matredundant);
PetscErrorCode MatGetRowMin_SymTrans(Mat mat,Vec v,PetscInt idx[]);
PetscErrorCode MatGetColumnVector_SymTrans(Mat A,Vec yy,PetscInt col);


#endif


