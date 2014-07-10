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

#include "petsc.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscvec.h"

#include "petscext.h"
#include "petscext_ksp.h"




int test_cgr( void )
{
	Vec            x, b, b_perturbed, u;      /* approx solution, RHS, exact solution */
	Mat            A;            /* linear system matrix */
	PetscErrorCode ierr;
	PetscInt       i,n = 200,col[3],its;
	PetscScalar    one = 1.0,value[3];
	PetscReal      norm;         /* norm of solution error */
	PC             pc;
	double         epsilon, val, fac;
	KSP            ksp;
	
	
	
	ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
	ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
	ierr = VecSetType( x, VECSEQ ); CHKERRQ(ierr);
	
	ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&u);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&b_perturbed);CHKERRQ(ierr);
	
	
	ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
	ierr = MatSetType( A, MATSEQAIJ );CHKERRQ(ierr);
	
	/* 
	Assemble matrix
	*/
	value[0] = -1.0; value[1] = 2.0; value[2] = -1.0;
	for (i=1; i<n-1; i++) {
		col[0] = i-1; col[1] = i; col[2] = i+1;
		ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
	}
	i = n - 1; col[0] = n - 2; col[1] = n - 1;
	ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
	i = 0; col[0] = 0; col[1] = 1; value[0] = 2.0; value[1] = -1.0;
	ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	/* Set exact solution; then compute right-hand-side vector. */
	ierr = VecSet(u,one);CHKERRQ(ierr);
	ierr = MatMult(A,u,b);CHKERRQ(ierr);
	
	/* create pertubed rhs b +/- epsilon */
	fac = 1.0e-1;
	/*printf("||perturbation|| = %g \n", fac ); */
	for( i=0; i<n; i++ ) {
		epsilon = -1.0 + 2.0*( (double)random() ) / ( (double) RAND_MAX );
		epsilon = fac * epsilon;
		//val = 1.0 + epsilon;
		val = epsilon;
		VecSetValue( b_perturbed, i, val, INSERT_VALUES );
	}
	VecAssemblyBegin(b_perturbed);
	VecAssemblyEnd(b_perturbed);
	
	VecAXPY( b_perturbed, 1.0, b ); // b_perturbed <- 1.0 * b + b_perturbed 
	
	
	
	/* Create ctx for generalized cg solver */
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	ierr = Stg_KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
	
	KSPSetType( ksp, "cgr" );
	KSPGetPC( ksp, &pc );
	PCSetType( pc, "none" );
	
//	KSPCGRSetRestart( ksp, 200 );
//	KSPCGRSetPreAllocateVectors( ksp, PETSC_FALSE );
//	KSPCGRSetStoreAllDirections( ksp, PETSC_FALSE );
	
	KSPSetRelativeRhsConvergenceTest( ksp );
	
	KSPSetFromOptions( ksp );
	
	
	/* Initial solve */
	printf("Init solve \n");
	KSPSolve( ksp, b, x );
	
	MatMult( A, x, u );
	VecAYPX( u, -1.0, b );
	VecNorm( u, NORM_2, &norm );
	KSPGetIterationNumber( ksp, &its );
	/*printf("[%d] Norm of error %g \n", its, norm );*/
	printf("its %d \n", its );
	
	/* Same rhs, should solve in 0 iterations */
	printf("Same rhs \n");
	ierr = Stg_KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
	KSPSolve( ksp, b, x );
	MatMult( A, x, u );
	VecAYPX( u, -1.0, b ); /* u <- -u + b */
	VecNorm( u, NORM_2, &norm );
	KSPGetIterationNumber( ksp, &its );
	/*printf("[%d] Norm of error %g \n", its, norm );*/
	printf("its %d \n", its );
	
	/* Scale rhs, should converge in 0 iteration */
	printf("scale rhs solve \n");
	VecScale( b, 10.0 );
	KSPSolve( ksp, b, x );
	MatMult( A, x, u );
	VecAYPX( u, -1.0, b ); /* u <- -u + b */
	VecNorm( u, NORM_2, &norm );
	KSPGetIterationNumber( ksp, &its );
	/*printf("[%d] Norm of error %g \n", its, norm );*/
	printf("its %d \n", its );
	
	/* perturb rhs solve */
	printf("perturb rhs solve \n");
	KSPSolve( ksp, b_perturbed, x );
	MatMult( A, x, u );
	VecAYPX( u, -1.0, b_perturbed );
	VecNorm( u, NORM_2, &norm );
	KSPGetIterationNumber( ksp, &its );
	/*printf("[%d] Norm of error %g \n", its, norm );*/
	printf("its %d \n", its );
	
	printf("Reset search directions \n");
	KSPCGRResetSearchDirection( ksp );
	
	printf("Init perturbed solve \n");
/*	KSPSetFromOptions(ksp); */
	KSPSolve( ksp, b_perturbed, x );
	MatMult( A, x, u );
	VecAYPX( u, -1.0, b_perturbed );
	KSPGetIterationNumber( ksp, &its );
	/*printf("[%d] Norm of error %g \n", its, norm );*/
	printf("its %d \n", its );
	
	
	
	KSPDestroy( & ksp );
	
	
	ierr = Stg_VecDestroy( &x);CHKERRQ(ierr); 
	ierr = Stg_VecDestroy( &u);CHKERRQ(ierr);
	ierr = Stg_VecDestroy( &b);CHKERRQ(ierr); 
	ierr = Stg_MatDestroy(A);CHKERRQ(ierr);
	Stg_VecDestroy( & b_perturbed );
	
	return 0;
}



int main(int argc,char **args)
{
	int k;
	
	PetscInitialize(&argc,&args,(char *)0,PETSC_NULL);
	
	PetscExtKSPRegisterAll();
	
	for( k=0; k<1; k++ ) {
		test_cgr();
	}
	
	PetscExtKSPRegisterDestroyAll();
	
	PetscFinalize();
	
	return 0;
}

