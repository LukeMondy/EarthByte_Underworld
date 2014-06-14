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


#ifndef __PETSCEXT_HELPERUTILS_H__
#define __PETSCEXT_HELPERUTILS_H__

PETSC_EXTERN_CXX_BEGIN

Mat MATAIJNew( MPI_Comm comm, PetscInt M, PetscInt N, const char prefix[] );

Mat MATAIJIdentityNew( MPI_Comm comm, PetscInt M );

PetscErrorCode MatFillStride( Mat A, PetscScalar start_val, PetscScalar stride );

PetscErrorCode KSPSetRelativeRhsConvergenceTest( KSP ksp );

PetscErrorCode VecLoad_MatrixMarket( MPI_Comm comm, const char fname[], VecType outtype, Vec *x );
PetscErrorCode VecView_MatrixMarket( Vec x, PetscViewer viewer );
PetscErrorCode MatAIJLoad_MatrixMarket( MPI_Comm comm, const char fname[], Mat *A );
PetscErrorCode MatView_MatrixMarket( Mat A, PetscViewer viewer );

PetscErrorCode __MatLoadMatrixMarket( MPI_Comm comm, const char fname[], const MatType type, Mat *A );
PetscErrorCode __MatViewMatrixMarket( Mat A, PetscViewer mm_viewer );

PetscErrorCode PETSCMAT_DLLEXPORT MatNullSpaceQuery(Mat mat,PetscReal *cnst_nullsp_val,PetscTruth *has_cnst_nullsp,PetscReal nullsp_val[],PetscTruth has_nullsp[]);

PetscErrorCode MatFreeComputeExplicitOperator( Mat mf, MatType type, PetscReal tol, Mat *_A);

PetscErrorCode MatViewContents(Mat A,const char name[]);
PetscErrorCode MatViewInfo(Mat A,const char name[]);
PetscErrorCode VecViewContents(Vec x,const char name[]);

PetscErrorCode PetscExtLoadMPIAIJOperations(Mat mat);

#define PETSCEXT_NULLSPACE_SIZE 1.0e-7
PetscErrorCode PetscExtMatComputeConstantNullSpace( Mat A, PetscReal *nrm );
PetscErrorCode PetscExtMatContainsConstantNullSpace( Mat A, PetscTruth *has_cnst_nullsp );
PetscErrorCode PetscExtVecRemoveConstantNullspace( Vec v );


PETSC_EXTERN_CXX_END
#endif


