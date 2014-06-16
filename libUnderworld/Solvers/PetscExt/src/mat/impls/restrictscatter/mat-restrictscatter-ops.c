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
#include "private/mat/mat-restrictscatter.h"

#include "petscext_utils.h"



PetscErrorCode MatMult_RestrictScatter(Mat mat,Vec x,Vec y)
{
	Mat_RestrictScatter s = (Mat_RestrictScatter)mat->data;
	
	PetscFunctionBegin;
	
	/* init x_work */
	/* We probably don't need to do this, but the index set may change.
	In this case, there may be junk lying around in x_work which should
	be killed off so the restricted mat-vec using B is correct */
	VecSet( s->x_work, 0.0 );
	
	VecSet( s->y_work, 0.0 );
	VecSet( y, 0.0 );
	
	/* x -> x_work */
	VecScatterBegin( s->scat_x, x,s->x_work, INSERT_VALUES, SCATTER_FORWARD );
	VecScatterEnd( s->scat_x, x,s->x_work,  INSERT_VALUES, SCATTER_FORWARD );
	
	/* do mat-vec on original operator */
	MatMult( s->A, s->x_work, s->y_work );
	
	/* y_work -> y_ */
	VecScatterBegin( s->scat_y, s->y_work,y,  INSERT_VALUES, SCATTER_FORWARD );
	VecScatterEnd( s->scat_y, s->y_work,y,  INSERT_VALUES, SCATTER_FORWARD );
	
	PetscFunctionReturn(0);
}

