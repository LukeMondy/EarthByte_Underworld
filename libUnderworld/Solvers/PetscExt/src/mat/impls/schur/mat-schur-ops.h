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


#ifndef __PETSC_MAT_OPS_SCHUR_H__
#define __PETSC_MAT_OPS_SCHUR_H__

/* 0*/
PetscErrorCode MatSetValues_Schur(Mat mat,PetscInt m,const PetscInt idxm[],PetscInt n,const PetscInt idxn[],const PetscScalar v[],InsertMode addv);
PetscErrorCode MatGetRow_Schur(Mat mat,PetscInt row,PetscInt *ncols,PetscInt *cols[],PetscScalar *vals[]);
PetscErrorCode MatRestoreRow_Schur(Mat mat,PetscInt row,PetscInt *ncols,PetscInt *cols[],PetscScalar *vals[]);
//PetscErrorCode MatMult_Schur(Mat mat,Vec x,Vec y);
PetscErrorCode MatMult_Schur_11(Mat mat,Vec x,Vec y);
PetscErrorCode MatMult_Schur_22(Mat mat,Vec x,Vec y);
PetscErrorCode MatMultAdd_Schur(Mat mat,Vec v1,Vec v2,Vec v3);
/* 5*/
//PetscErrorCode MatMultTranspose_Schur(Mat mat,Vec x,Vec y);
PetscErrorCode MatMultTranspose_Schur_11(Mat mat,Vec x,Vec y);
PetscErrorCode MatMultTranspose_Schur_22(Mat mat,Vec x,Vec y);

PetscErrorCode MatMultTransposeAdd_Schur(Mat mat,Vec v1,Vec v2,Vec v3);
PetscErrorCode MatSolve_Schur(Mat mat,Vec b,Vec x);
PetscErrorCode MatSolveAdd_Schur(Mat mat,Vec b,Vec y,Vec x);
PetscErrorCode MatSolveTranspose_Schur(Mat mat,Vec b,Vec x);
/*10*/
PetscErrorCode MatSolveTransposeAdd_Schur(Mat mat,Vec b,Vec y,Vec x);
PetscErrorCode MatLUFactor_Schur(Mat mat,IS row,IS col,MatFactorInfo *info);
PetscErrorCode MatCholeskyFactor_Schur(Mat mat,IS perm,MatFactorInfo *info);
PetscErrorCode MatRelax_Schur(Mat mat,Vec b,PetscReal omega,MatSORType flag,PetscReal shift,PetscInt its,PetscInt lits,Vec x);
PetscErrorCode MatTranspose_Schur(Mat mat,MatReuse reuse,Mat *B);
/*15*/
PetscErrorCode MatGetInfo_Schur(Mat mat,MatInfoType flag,MatInfo *info);
PetscErrorCode MatEqual_Schur(Mat A,Mat B,PetscTruth *flg);
PetscErrorCode MatGetDiagonal_Schur(Mat mat,Vec v);
PetscErrorCode MatDiagonalScale_Schur(Mat mat,Vec l,Vec r);
PetscErrorCode MatNorm_Schur(Mat mat,NormType type,PetscReal *nrm);
/*20*/
PetscErrorCode MatAssemblyBegin_Schur(Mat mat,MatAssemblyType type);
PetscErrorCode MatAssemblyEnd_Schur(Mat mat,MatAssemblyType type);
PetscErrorCode MatCompress_Schur(Mat mat);
PetscErrorCode MatSetOption_Schur(Mat mat,MatOption op);
PetscErrorCode MatZeroEntries_Schur(Mat mat);
/*25*/
PetscErrorCode MatZeroRows_Schur(Mat mat,PetscInt numRows,const PetscInt rows[],PetscScalar diag);
PetscErrorCode MatLUFactorSymbolic_Schur(Mat mat,IS row,IS col,MatFactorInfo *info,Mat *fact);
PetscErrorCode MatLUFactorNumeric_Schur(Mat mat,MatFactorInfo *info,Mat *fact);
PetscErrorCode MatCholeskyFactorSymbolic_Schur(Mat mat,IS perm,MatFactorInfo *info,Mat *fact);
PetscErrorCode MatCholeskyFactorNumeric_Schur(Mat mat,MatFactorInfo *info,Mat *fact);
/*30*/
PetscErrorCode MatSetUpPreallocation_Schur(Mat B);
PetscErrorCode MatILUFactorSymbolic_Schur(Mat mat,IS row,IS col,MatFactorInfo *info,Mat *fact);
PetscErrorCode MatICCFactorSymbolic_Schur(Mat mat,IS perm,MatFactorInfo *info,Mat *fact);
PetscErrorCode MatGetArray_Schur(Mat mat,PetscScalar *v[]);
PetscErrorCode MatRestoreArray_Schur(Mat mat,PetscScalar *v[]);
/*35*/
PetscErrorCode MatDuplicate_Schur(Mat mat,MatDuplicateOption op,Mat *M);
PetscErrorCode MatForwardSolve_Schur(Mat mat,Vec b,Vec x);
PetscErrorCode MatBackwardSolve_Schur(Mat mat,Vec b,Vec x);
PetscErrorCode MatILUFactor_Schur(Mat mat,IS row,IS col,MatFactorInfo *info);
PetscErrorCode MatICCFactor_Schur(Mat mat,IS row,MatFactorInfo* info);
/*40*/
PetscErrorCode MatAXPY_Schur(Mat Y,PetscScalar a,Mat X,MatStructure str);
PetscErrorCode MatGetSubMatrices_Schur(Mat mat,PetscInt n,const IS irow[],const IS icol[],MatReuse scall,Mat *submat[]);
PetscErrorCode MatIncreaseOverlap_Schur(Mat mat,PetscInt n,IS is[],PetscInt ov);
PetscErrorCode MatGetValues_Schur(Mat mat,PetscInt m,const PetscInt idxm[],PetscInt n,const PetscInt idxn[],PetscScalar v[]);
PetscErrorCode MatCopy_Schur(Mat A,Mat B,MatStructure str);
/*45*/
PetscErrorCode MatGetRowMax_Schur(Mat mat,Vec v,PetscInt idx[]);
PetscErrorCode MatScale_Schur(Mat mat,PetscScalar a);
PetscErrorCode MatShift_Schur( Mat mat, PetscScalar a );
PetscErrorCode MatDiagonalSet_Schur( Mat Y, Vec D, InsertMode is );
PetscErrorCode MatILUDTFactor_Schur(Mat mat,IS row,IS col,MatFactorInfo *info,Mat *fact);
/*50*/
PetscErrorCode MatSetBlockSize_Schur(Mat mat,PetscInt bs);
PetscErrorCode MatGetRowIJ_Schur(Mat mat,PetscInt shift,PetscTruth symmetric,PetscTruth blockcompressed,PetscInt *n,PetscInt *ia[],PetscInt* ja[],PetscTruth *done);
PetscErrorCode MatRestorRowIJ_Schur(Mat mat,PetscInt shift,PetscTruth symmetric,PetscTruth blockcompressed,PetscInt *n,PetscInt *ia[],PetscInt* ja[],PetscTruth *done);
PetscErrorCode MatGetColumnIJ_Schur(Mat mat,PetscInt shift,PetscTruth symmetric,PetscTruth blockcompressed,PetscInt *n,PetscInt *ia[],PetscInt* ja[],PetscTruth *done);
PetscErrorCode MatRestoreColumnIJ_Schur(Mat mat,PetscInt shift,PetscTruth symmetric,PetscTruth blockcompressed,PetscInt *n,PetscInt *ia[],PetscInt* ja[],PetscTruth *done);
/*55*/
PetscErrorCode MatFDColoringCreate_Schur(Mat mat,ISColoring iscoloring,MatFDColoring color);
PetscErrorCode MatColoringPatch_Schur(Mat mat,PetscInt ncolors,PetscInt n,ISColoringValue colorarray[],ISColoring *iscoloring);
PetscErrorCode MatSetUnfactored_Schur(Mat mat);
PetscErrorCode MatPermute_Schur( Mat mat, IS row, IS col, Mat *B );
PetscErrorCode MatSetValuesBlocked_Schur(Mat mat,PetscInt m,const PetscInt idxm[],PetscInt n,const PetscInt idxn[],const PetscScalar v[],InsertMode addv);
/*60*/
PetscErrorCode MatGetSubMatrix_Schur(Mat mat,IS isrow,IS iscol,PetscInt csize,MatReuse cll,Mat *newmat);
/*PetscErrorCode MatDestroy_Schur( Mat mat );*/
PetscErrorCode MatView_Schur( Mat mat, PetscViewer viewer );
PetscErrorCode MatConvertFrom_Schur(Mat mat, MatType newtype,MatReuse reuse,Mat *M);
PetscErrorCode MatUseScaledForm_Schur( Mat mat, PetscTruth scaled );
/*65*/
PetscErrorCode MatScaleSystem_Schur( Mat mat, Vec b, Vec x );
PetscErrorCode MatUnScaleSystem_Schur(Mat mat,Vec b,Vec x);
PetscErrorCode MatSetLocalToGlobalMapping_Schur(Mat x,ISLocalToGlobalMapping mapping);
PetscErrorCode MatSetValuesLocal_Schur(Mat mat,PetscInt nrow,const PetscInt irow[],PetscInt ncol,const PetscInt icol[],const PetscScalar y[],InsertMode addv);
PetscErrorCode MatZeroRowsLocal_Schur(Mat mat,PetscInt numRows,const PetscInt rows[],PetscScalar diag);
/*70*/
PetscErrorCode MatGetRowMaxAbs_Schur(Mat mat,Vec v,PetscInt idx[]);
PetscErrorCode MatConvert_Schur(Mat mat, MatType newtype,MatReuse reuse,Mat *M);
PetscErrorCode MatSetColoring_Schur( Mat mat, ISColoring colour );
PetscErrorCode MatSetValuesAdic_Schur(Mat mat,void *v);
PetscErrorCode MatSetValuesAdifor_Schur(Mat mat,PetscInt nl,void *v);
/*75*/
PetscErrorCode MatFDColoringApply_Schur(Mat J,MatFDColoring coloring,Vec x1,MatStructure *flag,void *sctx);
PetscErrorCode MatSetFromOptions_Schur( Mat mat );
PetscErrorCode MatMultConstrained_Schur(Mat mat,Vec x,Vec y);
PetscErrorCode MatMultTransposeConstrained_Schur(Mat mat,Vec x,Vec y);
PetscErrorCode MatILUFactorSymbolicConstrained_Schur(Mat mat,IS row,IS col,double d,PetscInt i,PetscInt j,Mat *fact);
/*80*/
PetscErrorCode MatPermuteSparsify_Schur(Mat A, PetscInt band, PetscReal frac, PetscReal tol, IS rowp, IS colp, Mat *B);
PetscErrorCode MatMults_Schur(Mat A,Vecs x,Vecs y);
PetscErrorCode MatSolves_Schur(Mat A,Vecs b,Vecs x);
PetscErrorCode MatGetInertia_Schur(Mat mat,PetscInt *nneg,PetscInt *nzero,PetscInt *npos);
PetscErrorCode MatLoad_Schur(PetscViewer viewer, MatType outtype,Mat *newmat);
/*85*/
PetscErrorCode MatIsSymmetric_Schur(Mat A,PetscReal tol,PetscTruth *flg);
PetscErrorCode MatIsHermitian_Schur(Mat A,PetscTruth *flg);
PetscErrorCode MatIsStructurallySymmetric_Schur(Mat A,PetscTruth *flg);
PetscErrorCode MatPBRelax_Schur(Mat mat,Vec b,PetscReal omega,MatSORType flag,PetscReal shift,PetscInt its,PetscInt lits,Vec x);
PetscErrorCode MatGetVecs_Schur(Mat mat,Vec *right,Vec *left);
/*90*/
PetscErrorCode MatMatMult_Schur(Mat A,Mat B,MatReuse scall,PetscReal fill,Mat *C);
PetscErrorCode MatMatMultSymbolic_Schur(Mat A,Mat B,PetscReal fill,Mat *C);
PetscErrorCode MatMatMultNumeric_Schur(Mat A,Mat B,Mat C);
PetscErrorCode MatPtAP_Schur(Mat A,Mat P,MatReuse scall,PetscReal fill,Mat *C);
PetscErrorCode MatPtAPSymbolic_Schur(Mat A,Mat P,PetscReal fill,Mat *C);
/*95*/
PetscErrorCode MatPtAPNumeric_Schur(Mat A,Mat P,Mat C);
PetscErrorCode MatMatMultTranspose_Schur(Mat A,Mat B,MatReuse scall,PetscReal fill,Mat *C);
PetscErrorCode MatMatMultTransposeSymbolic_Schur(Mat A,Mat B,PetscReal fill,Mat *C);
PetscErrorCode MatMatMultTransposeNumeric_Schur(Mat A,Mat B,Mat C);
PetscErrorCode MatPtAPSymbolic_SEQAIJ_Schur(Mat A,Mat P,PetscReal fill,Mat *C);
/*100*/
PetscErrorCode MatPtAPNumeric_SEQAIJ_Schur(Mat A,Mat P,Mat C);
PetscErrorCode MatPtAPSymbolic_MPIAIJ_Schur(Mat A,Mat P,PetscReal fill,Mat *C);
PetscErrorCode MatPtAPNumeric_MPIAIJ_Schur(Mat A,Mat P,Mat C);
PetscErrorCode MatConjugate_Schur(Mat mat);
PetscErrorCode MatSetSizes_Schur(Mat mat,PetscInt m,PetscInt n,PetscInt M,PetscInt N);
/*105*/
PetscErrorCode MatSetValuesRow_Schur(Mat mat,PetscInt row,const PetscScalar v[]);
PetscErrorCode MatRealPart_Schur(Mat mat);
PetscErrorCode MatImaginaryPart_Schur(Mat mat);
PetscErrorCode MatGetRowUpperTriangular_Schur(Mat mat);
PetscErrorCode MatRestoreRowUpperTriangular_Schur(Mat mat);
/*110*/
PetscErrorCode MatMatSolve_Schur(Mat mat,Mat b,Mat x);
PetscErrorCode MatGetRedundantMatrix_Schur(Mat mat,PetscInt nsubcomm,MPI_Comm subcomm,PetscInt mlocal_red,MatReuse reuse,Mat *matredundant);
PetscErrorCode MatGetRowMin_Schur(Mat mat,Vec v,PetscInt idx[]);
PetscErrorCode MatGetColumnVector_Schur(Mat A,Vec yy,PetscInt col);


#endif

