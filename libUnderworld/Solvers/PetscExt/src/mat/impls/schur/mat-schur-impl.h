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
/*

[A11 A12]{x1}={f1} 
[A21 A22]{x2} {f2}

"MatSchur_A11" -----------------------------------------------------

Left multiply by eq (1) by A21 A11^{-1}

A21 {x1} + A21 A11^{-1} A12 {x2} = A21 A11^{-1} {f1}
{f2} - A22 {x2} + A21 A11^{-1} A12 {x2} = A21 A11^{-1} {f1}

( A21 A11^{-1} A12 - A22 ) {x2} = A21 A11^{-1} {f1} - {f2}


"MatSchur_A22" -----------------------------------------------------

Left multiply by eq (2) by A12 A22^{-1}

A12 A22^{-1} A21 {x1} + A12 {x2} = A12 A22^{-1} {f2}
A12 A22^{-1} A21 {x1} + {f1} - A11 {x1} = A12 A22^{-1} {f2}

( A12 A22^{-1} A21 - A11 ) {x1} = A12 A22^{-1} {f2} - {f1} 

*/

#ifndef __PETSC_MAT_IMPL_SCHUR_H__
#define __PETSC_MAT_IMPL_SCHUR_H__

#include "private/mat/mat-schur.h"
#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

typedef struct _Mat_Schur *Mat_Schur;

struct _Mat_Schur {
	MatSchurComplementType stype;
	PetscTruth operators_set, scale_set;
	PetscScalar alpha;
	KSP ksp;
	Mat A11, A12, A21, A22;
	Vec t1,t2 , t1a,t2a;
};



#endif

