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


#include <petscvec.h>
#include "petscext_vec.h"


/*
Does not support nested blocks.
*/
PetscErrorCode VecGetNormBlock( Vec x, NormType type, PetscInt *bs, PetscReal *val[] )
{
	PetscErrorCode ierr;
	PetscInt i, _bs;
	PetscReal *_val;
	PetscTruth is_block;
	Vec sub_x;
	
	PetscFunctionBegin;
	
	/* Get block size */
	Stg_PetscTypeCompare( (PetscObject)x, "block", &is_block );
	/* Catch case if Vec is not a block vector */
	if( is_block == PETSC_FALSE ) {
		ierr=PetscMalloc( sizeof(PetscReal), &_val );CHKERRQ(ierr);
		
		VecNorm( x, type, &_val[0] );
		
		*bs = 1;
		*val = _val;
		PetscFunctionReturn(0);
	}
	else {
		VecGetSize( x, &_bs );
	}
	
	/* compute norm for each block */
	PetscMalloc( _bs * sizeof(PetscReal), &_val );
	for( i=0; i<_bs; i++ ) {
		
		VecBlockGetSubVector( x, i, &sub_x );
		Stg_PetscTypeCompare( (PetscObject)sub_x, "block", &is_block );
		if( is_block == PETSC_TRUE ) {
			Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "VecGetNormBlock() does not support nested blocks." );
		}
		
		VecNorm( sub_x, type, &_val[i] );
		VecBlockRestoreSubVectors( x );
	}
	
	*bs = _bs;
	*val = _val;
	PetscFunctionReturn(0);
}

PetscErrorCode VecRestoreNormBlock( Vec x, NormType type, PetscInt *bs, PetscReal *val[] )
{
	PetscErrorCode ierr;

	PetscFunctionBegin;
	ierr=PetscFree( *val );CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

