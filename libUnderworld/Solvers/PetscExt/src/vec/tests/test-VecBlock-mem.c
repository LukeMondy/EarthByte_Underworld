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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>

#include "petscext_vec.h"



Vec gen_test_vector( MPI_Comm comm, PetscInt length, PetscInt start_value, PetscInt stride )
{
	int nproc;
	Vec v;
	PetscInt i;
	PetscScalar vx;
	
	MPI_Comm_size( comm, &nproc );
	
	VecCreate( comm, &v );
	VecSetSizes( v, PETSC_DECIDE, length );
	if( nproc == 1 ) VecSetType( v, VECSEQ );
	else VecSetType( v, VECMPI );
	
	for( i=0; i<length; i++ ) {
		vx = (PetscScalar)( start_value + i * stride );
		VecSetValue( v, i, vx, INSERT_VALUES );
	}
	VecAssemblyBegin( v );
	VecAssemblyEnd( v );
	
	return v;
	
}


void test_mem_usage( void )
{
	Vec x1,x2;
	Vec X;
	PetscReal val;
	PetscInt index;
	
	
	x1 = gen_test_vector( PETSC_COMM_WORLD, 4, 0, 1 );
	x2 = gen_test_vector( PETSC_COMM_WORLD, 5, 10, 2 );
	
	VecCreate( PETSC_COMM_WORLD, &X );
	VecSetSizes( X, 2,2 );
	VecSetType( X, "block" );
	VecBlockSetValue( X, 0, x1, INSERT_VALUES );
	VecBlockSetValue( X, 1, x2, INSERT_VALUES );
	VecAssemblyBegin(X);
	VecAssemblyEnd(X);
	Stg_VecDestroy( & x1 );			Stg_VecDestroy( & x2 );
	
	
	VecMax( X, &index, &val );
	
	Stg_VecDestroy( & X );
}







int main( int argc, char **args )
{
	int i;
	int BIG=100000;
	
	PetscInitialize( &argc, &args,(char *)0, PETSC_NULL );
	PetscExtVecRegisterAll();
	
	for( i=0; i<BIG; i++ ) {
		test_mem_usage();
	}
	
	
	
	PetscExtVecRegisterDestroyAll();
	PetscFinalize();
	return 0;
}
