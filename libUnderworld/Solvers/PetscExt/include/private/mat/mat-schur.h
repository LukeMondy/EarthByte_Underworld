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

#ifndef __PETSC_MAT_SCHUR_H__
#define __PETSC_MAT_SCHUR_H__

#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

#include "private/compat/petsccompat.h"

PETSC_EXTERN_CXX_BEGIN

typedef const char* MatSchurComplementType;


/* private constructor */
//PetscErrorCode MatCreate_Schur( Mat A );
// Declare this as extern when we register it.

/* public constructor */
PetscErrorCode MatCreateSchur( MPI_Comm comm, Mat A11,Mat A12,Mat A21,Mat A22, PetscScalar alpha, MatSchurComplementType type, Mat *A );
PetscErrorCode MatCreateSchurFromBlock( Mat bmat, PetscScalar alpha, MatSchurComplementType type, Mat *A );

/* implemenation specific interfaces */
PetscErrorCode MatSchurSetOperators( Mat A, Mat A11,Mat A12,Mat A21,Mat A22);
PetscErrorCode MatSchurSetOperatorsFromBlock( Mat A, Mat BlockA );

PetscErrorCode MatSchurApplyReductionToVec( Mat A, Vec f1, Vec f2, Vec subb );
PetscErrorCode MatSchurApplyReductionToVecFromBlock( Mat A, Vec F, Vec subb );

PetscErrorCode MatSchurSetSchurComplementType( Mat A, MatSchurComplementType type );
PetscErrorCode MatSchurGetSchurComplementType( Mat A, MatSchurComplementType *type );

PetscErrorCode MatSchurSetScalar( Mat A, PetscScalar alpha );
PetscErrorCode MatSchurGetScalar( Mat A, PetscScalar *alpha );

PetscErrorCode MatSchurSetKSP( Mat A, KSP ksp );
PetscErrorCode MatSchurGetKSP( Mat A, KSP *ksp );


PETSC_EXTERN_CXX_END
#endif

