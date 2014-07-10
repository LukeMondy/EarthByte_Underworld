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


#ifndef __vec_block_impl_h__
#define __vec_block_impl_h__

#include <petsc.h>
#include <petscmat.h>
#include <petscvec.h>
#include "private/compat/petsccompat.h"

typedef struct {
	PetscInt   nb; /* n blocks */
	Vec        *v;
	PetscTruth setup_called;
} Vec_Block;




/* operations */
PetscErrorCode VecBlockRestoreSubVectors_Block( Vec X );
PetscErrorCode PETSCMAT_DLLEXPORT VecBlockRestoreSubVectors( Vec x );
PetscErrorCode VecBlockSetValue( Vec x, PetscInt idxm, Vec vec, InsertMode addv );

PetscErrorCode VecCopy_Block( Vec x, Vec y );
PetscErrorCode VecDuplicate_Block( Vec x, Vec *y );

PetscErrorCode VecDot_Block( Vec x, Vec y, PetscScalar *val );
PetscErrorCode VecTDot_Block( Vec x, Vec y, PetscScalar *val );
PetscErrorCode VecAXPY_Block( Vec y, PetscScalar alpha, Vec x );
PetscErrorCode VecAYPX_Block( Vec y, PetscScalar alpha, Vec x );
PetscErrorCode VecAXPBY_Block( Vec y, PetscScalar alpha, PetscScalar beta, Vec x );

PetscErrorCode VecScale_Block( Vec x, PetscScalar alpha );
PetscErrorCode VecPointwiseMult_Block( Vec w, Vec x, Vec y );
PetscErrorCode VecPointwiseDivide_Block( Vec w, Vec x, Vec y );

PetscErrorCode VecReciprocal_Block( Vec x );
PetscErrorCode VecNorm_Block( Vec xin, NormType type, PetscReal* z );
PetscErrorCode VecMAXPY_Block( Vec y, PetscInt nv, const PetscScalar alpha[], Vec *x );
PetscErrorCode VecMDot_Block( Vec x, PetscInt nv, const Vec y[], PetscScalar *val );
PetscErrorCode VecMTDot_Block( Vec x, PetscInt nv, const Vec y[], PetscScalar *val );
PetscErrorCode VecSet_Block( Vec x, PetscScalar alpha );
PetscErrorCode VecConjugate_Block( Vec x );
PetscErrorCode VecSwap_Block( Vec x, Vec y );
PetscErrorCode VecWAXPY_Block( Vec w, PetscScalar alpha, Vec x, Vec y );

PetscErrorCode VecMax_Block( Vec x, PetscInt *p, PetscReal *val );
PetscErrorCode VecMin_Block( Vec x, PetscInt *p, PetscReal *val );
PetscErrorCode VecView_Block( Vec x, PetscViewer viewer );
PetscErrorCode VecGetSize_Block(Vec x,PetscInt *size);

PetscErrorCode VecMaxPointwiseDivide_Block(Vec x,Vec y,PetscReal *max);

/* NOT SUPPORTED */
PetscErrorCode __VecSetValues_Block( Vec x,PetscInt ni,const PetscInt ix[],const PetscScalar y[],InsertMode iora );
PetscErrorCode __VecGetArray_Block(Vec x,PetscScalar *a[]);
PetscErrorCode __VecRestoreArray_Block(Vec x,PetscScalar *a[]);
PetscErrorCode __VecSetRandom_Block(Vec x,PetscRandom rctx);

PetscErrorCode __VecSetValuesBlocked_Block(Vec x,PetscInt ni,const PetscInt ix[],const PetscScalar y[],InsertMode iora);
PetscErrorCode __VecPlaceArray_Block(Vec vec,const PetscScalar array[]);
PetscErrorCode __VecReplaceArray_Block(Vec vec,const PetscScalar array[]);
PetscErrorCode __VecLoadIntoVector_Block(PetscViewer viewer,Vec vec);
PetscErrorCode __VecResetArray_Block(Vec vec);
PetscErrorCode  __VecPointwiseMax_Block(Vec w,Vec x,Vec y);

PetscErrorCode  __VecPointwiseMaxAbs_Block(Vec w,Vec x,Vec y);
PetscErrorCode  __VecPointwiseMin_Block(Vec w,Vec x,Vec y);
PetscErrorCode  __VecGetValues_Block(Vec x,PetscInt ni,const PetscInt ix[],PetscScalar y[]);


/*
Debugging tips
To look at block data when we only have Vec b

p *(*(Vec_Block*)(b->data))->v[0]


*/


#endif
