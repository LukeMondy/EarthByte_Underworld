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



/* 
Usage;

VecCreate( );
VecSetSizes( PETSC_DECIDE, 2 ); // block size 
VecSetType( v, "block" ); // calls VecCreate_Block


*/


/*
X = ( a b )^T
a = ( c d )^T
b = ( e f )^T

c,d,e,f = seq


*/

void test_view( void )
{
	Vec X, a,b;
	Vec c,d,e,f;
	PetscInt index;
	PetscReal val;
	
	
	PetscPrintf( PETSC_COMM_WORLD, "\n\n============== %s ==============\n", __func__ );
	
	
	VecCreate( PETSC_COMM_WORLD, &X );
	VecSetSizes( X, 2, 2 );
	VecSetType( X, "block" );
	
	VecCreate( PETSC_COMM_WORLD, &a );
	VecSetSizes( a, 2, 2 );
	VecSetType( a, "block" );
	
	VecCreate( PETSC_COMM_WORLD, &b );
	VecSetSizes( b, 2, 2 );
	VecSetType( b, "block" );
	
	/* assemble X */
	VecBlockSetValue( X, 0, a, INSERT_VALUES );
	VecBlockSetValue( X, 1, b, INSERT_VALUES );
	VecAssemblyBegin(X);
	VecAssemblyEnd(X);
	
	
	VecCreate( PETSC_COMM_WORLD, &c );
	VecSetSizes( c, 3, 3 );
	VecSetType( c, VECSEQ );
	VecDuplicate( c, &d );
	VecDuplicate( c, &e );
	VecDuplicate( c, &f );
	
	VecSet( c, 1.0 );
	VecSet( d, 2.0 );
	VecSet( e, 3.0 );
	VecSetRandom( f, PETSC_NULL );
	VecScale( f, 10.0 );
	
	/* assemble a */
	VecBlockSetValue( a, 0, c, INSERT_VALUES );
	VecBlockSetValue( a, 1, d, INSERT_VALUES );
	VecAssemblyBegin(a);
	VecAssemblyEnd(a);
	
	
	/* assemble b */
	VecBlockSetValue( b, 0, e, INSERT_VALUES );
	VecBlockSetValue( b, 1, f, INSERT_VALUES );
	VecAssemblyBegin(b);
	VecAssemblyEnd(b);
	
//	PetscPrintf( PETSC_COMM_WORLD, "X \n");
//	VecView( X, PETSC_VIEWER_STDOUT_WORLD );
	
//	PetscPrintf( PETSC_COMM_WORLD, "b \n");
//	VecView( b, PETSC_VIEWER_STDOUT_WORLD );
	
	
	VecMax( b, &index, &val );
	PetscPrintf( PETSC_COMM_WORLD, "(max-b) = %f : index = %d \n", val, index );
	
	VecMin( b, &index, &val );
	PetscPrintf( PETSC_COMM_WORLD, "(min-b) = %f : index = %d \n", val, index );
	
	VecMax( X, &index, &val );
	PetscPrintf( PETSC_COMM_WORLD, "(max-X) = %f : index = %d \n", val, index );
	VecMin( X, &index, &val );
	PetscPrintf( PETSC_COMM_WORLD, "(min-X) = %f : index = %d \n", val, index );
	
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO_DETAIL );
	VecView( X, PETSC_VIEWER_STDOUT_WORLD );
}


void test_vec_ops( void )
{
	Vec X, a,b;
	Vec c,d,e,f;
	PetscInt index;
	PetscScalar val;
	
	
	PetscPrintf( PETSC_COMM_WORLD, "\n\n============== %s ==============\n", __func__ );
	
	VecCreate( PETSC_COMM_WORLD, &X );
	VecSetSizes( X, 2, 2 );
	VecSetType( X, "block" );
	
	VecCreate( PETSC_COMM_WORLD, &a );
	VecSetSizes( a, 2, 2 );
	VecSetType( a, "block" );
	
	VecCreate( PETSC_COMM_WORLD, &b );
	VecSetSizes( b, 2, 2 );
	VecSetType( b, "block" );
	
	/* assemble X */
	VecBlockSetValue( X, 0, a, INSERT_VALUES );
	VecBlockSetValue( X, 1, b, INSERT_VALUES );
	VecAssemblyBegin(X);
	VecAssemblyEnd(X);
	
	
	
	VecCreate( PETSC_COMM_WORLD, &c );
	VecSetSizes( c, 3, 3 );
	VecSetType( c, VECSEQ );
	VecDuplicate( c, &d );
	VecDuplicate( c, &e );
	VecDuplicate( c, &f );
	
	VecSet( c, 1.0 );
	VecSet( d, 2.0 );
	VecSet( e, 3.0 );
	VecSet( f, 4.0 );
	
	/* assemble a */
	VecBlockSetValue( a, 0, c, INSERT_VALUES );
	VecBlockSetValue( a, 1, d, INSERT_VALUES );
	VecAssemblyBegin(a);
	VecAssemblyEnd(a);
	
	/* assemble b */
	VecBlockSetValue( b, 0, e, INSERT_VALUES );
	VecBlockSetValue( b, 1, f, INSERT_VALUES );
	VecAssemblyBegin(b);
	VecAssemblyEnd(b);
	
//	PetscPrintf( PETSC_COMM_WORLD, "X \n");
//	VecView( X, PETSC_VIEWER_STDOUT_WORLD );
	
	VecDot( X,X, &val );
	PetscPrintf( PETSC_COMM_WORLD, "X.X = %f \n", val ); 
	
}



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


/*
X = ( [0,1,2,3], [10,12,14,16,18] )
Y = ( [4,7,10,13], [5,6,7,8,9] )

Y = aX + y = ( [4,8,12,16], (15,18,21,24,27] )
Y = aX + y = ( [4,9,14,19], (25,30,35,40,45] )

*/

void test_axpy_dot_max( void )
{
	Vec x1,y1, x2,y2;
	Vec X, Y;
	PetscReal real;
	PetscScalar scalar;
	PetscInt index;
	
	
	x1 = gen_test_vector( PETSC_COMM_WORLD, 4, 0, 1 );
	x2 = gen_test_vector( PETSC_COMM_WORLD, 5, 10, 2 );
	
	y1 = gen_test_vector( PETSC_COMM_WORLD, 4, 4, 3 );
	y2 = gen_test_vector( PETSC_COMM_WORLD, 5, 5, 1 );
	
	
	
	VecCreate( PETSC_COMM_WORLD, &X );
	VecSetSizes( X, 2,2 );
	VecSetType( X, "block" );
	VecBlockSetValue( X, 0, x1, INSERT_VALUES );
	VecBlockSetValue( X, 1, x2, INSERT_VALUES );
	VecAssemblyBegin(X);
	VecAssemblyEnd(X);
	Stg_VecDestroy( & x1 );			Stg_VecDestroy( & x2 );
	
	
	VecCreate( PETSC_COMM_WORLD, &Y );
	VecSetSizes( Y, 2,2 );
	VecSetType( Y, "block" );
	VecBlockSetValue( Y, 0, y1, INSERT_VALUES );
	VecBlockSetValue( Y, 1, y2, INSERT_VALUES );
	VecAssemblyBegin(Y);
	VecAssemblyEnd(Y);
	Stg_VecDestroy( & y1 );			Stg_VecDestroy( & y2 );
	
	
	PetscPrintf( PETSC_COMM_WORLD, "VecAXPY \n");
	VecAXPY( Y, 1.0, X ); /* Y <- a X + Y */
	VecBlockGetSubVector( Y, 0, &y1 );					VecBlockGetSubVector( Y, 1, &y2 );
	PetscPrintf( PETSC_COMM_WORLD, "(1) y1 = \n" );		VecView( y1, PETSC_VIEWER_STDOUT_WORLD );
	PetscPrintf( PETSC_COMM_WORLD, "(1) y2 = \n" );		VecView( y2, PETSC_VIEWER_STDOUT_WORLD );
	VecBlockRestoreSubVectors( Y );
	VecDot( X,Y, &scalar );
#ifdef PETSC_USE_COMPLEX
	PetscPrintf( PETSC_COMM_WORLD, "X.Y = %lf + %lfi \n", PetscRealPart(scalar), PetscImaginaryPart(scalar) ); 
#else
	PetscPrintf( PETSC_COMM_WORLD, "X.Y = %lf \n", scalar ); 
#endif 
	
	
	VecAXPY( Y, 1.0, X ); /* Y <- a X + Y */
	VecBlockGetSubVector( Y, 0, &y1 );					VecBlockGetSubVector( Y, 1, &y2 );
	PetscPrintf( PETSC_COMM_WORLD, "(2) y1 = \n" );		VecView( y1, PETSC_VIEWER_STDOUT_WORLD );
	PetscPrintf( PETSC_COMM_WORLD, "(2) y2 = \n" );		VecView( y2, PETSC_VIEWER_STDOUT_WORLD );
	VecBlockRestoreSubVectors( Y );
	VecDot( X,Y, &scalar );
#ifdef PETSC_USE_COMPLEX
	PetscPrintf( PETSC_COMM_WORLD, "X.Y = %lf + %lfi \n", PetscRealPart(scalar), PetscImaginaryPart(scalar) ); 
#else
	PetscPrintf( PETSC_COMM_WORLD, "X.Y = %lf \n", scalar ); 
#endif 
	
	
	
	VecMax( X, &index, &real );
	PetscPrintf( PETSC_COMM_WORLD, "(max-X) = %f : index = %d \n", real, index );
	VecMin( X, &index, &real );
	PetscPrintf( PETSC_COMM_WORLD, "(min-X) = %f : index = %d \n", real, index );
	
//	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO_DETAIL );
//	VecView( X, PETSC_VIEWER_STDOUT_WORLD );
	
	
	Stg_VecDestroy( & X );
	Stg_VecDestroy( & Y );
	
}



int main( int argc, char **args )
{
	int i;
	int BIG=100000;
	
	PetscInitialize( &argc, &args,(char *)0, PETSC_NULL );
	PetscExtVecRegisterAll();
	
//	test_view();	// this test will output pointers, thus we cannot diff result
	test_axpy_dot_max();
	
	PetscExtVecRegisterDestroyAll();
	PetscFinalize();
	return 0;
}
