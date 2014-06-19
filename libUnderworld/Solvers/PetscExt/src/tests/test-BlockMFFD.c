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


static char help[] = "Newton's method to solve a two-variable system, sequentially.\n"\
		"The same problem is solved twice - i) fully assembled system + ii) block system\n\n";
#if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >=3 ) )
  #include "petsc-private/snesimpl.h"
#else
  #include "private/snesimpl.h"
#endif

#include "petscsnes.h"
#include "petscext_vec.h"
#include "petscext_mat.h"

/*

- 
  4 SNES Function norm 1.936324324560e-03 
scale 3.0898464383e+06 
scale 3.0898464383e+06 
- 
  5 SNES Function norm 1.436805596222e-07 
scale 3.0901699197e+06 
scale 3.0901699197e+06 
- 
  6 SNES Function norm 1.098616031645e-14 
number of Newton iterations = 6

*/

/* User-defined routines */
extern PetscErrorCode FormFunction1_block(SNES,Vec,Vec,void*);
extern PetscErrorCode FormFunction2_block(SNES,Vec,Vec,void*);


#undef __FUNCT__
#define __FUNCT__ "block_system"
int block_system_force(int argc,char **argv)
{
	SNES           snes;         /* nonlinear solver context */
	KSP            ksp;         /* linear solver context */
	PC             pc;           /* preconditioner context */
	Vec            x,r;         /* solution, residual vectors */
	Mat            J;            /* Jacobian matrix */
	PetscErrorCode ierr;
	PetscInt       its;
	PetscScalar    pfive = .5;
	PetscTruth     flg;
	
	Vec x1, x2, r1, r2;
	Vec bv;
	
	
	PetscPrintf( PETSC_COMM_WORLD, "\n\n========================= Block system =========================\n\n" );
	
	ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Create matrix and vector data structures; set corresponding routines
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	
	/*
	Create sub vectors for solution and nonlinear function
	*/
	ierr = VecCreateSeq(PETSC_COMM_SELF,1,&x1);CHKERRQ(ierr);
	ierr = VecDuplicate(x1,&r1);CHKERRQ(ierr);
	
	ierr = VecCreateSeq(PETSC_COMM_SELF,1,&x2);CHKERRQ(ierr);
	ierr = VecDuplicate(x2,&r2);CHKERRQ(ierr);
	
	/*
	Create the block vectors
	*/
	VecCreate( PETSC_COMM_WORLD, &x );
	VecSetSizes( x, 2,2 );
	VecSetType( x, "block" );
	
	VecBlockSetValue( x, 0, x1, INSERT_VALUES );
	VecBlockSetValue( x, 1, x2, INSERT_VALUES );
	VecAssemblyBegin(x);
	VecAssemblyEnd(x);
	
	VecDuplicate( x, &r );
	VecBlockSetValue( r, 0, r1, INSERT_VALUES );
	VecBlockSetValue( r, 1, r2, INSERT_VALUES );
	VecAssemblyBegin(r);
	VecAssemblyEnd(r);
	
	
	ierr = PetscOptionsHasName(PETSC_NULL,"-hard",&flg);CHKERRQ(ierr);
	if (!flg) {
		/* Set function evaluation routine and vector. */
		ierr = SNESSetFunction(snes,r,FormFunction1_block,PETSC_NULL);CHKERRQ(ierr);
		MatCreateSNESBlockMFFD( snes, &J );
		SNESSetJacobian(snes,J,J,SNESDefaultComputeJacobian_Block, PETSC_NULL );
	} else {
		ierr = SNESSetFunction(snes,r,FormFunction2_block,PETSC_NULL);CHKERRQ(ierr);
		MatCreateSNESBlockMFFD( snes, &J );
		SNESSetJacobian(snes,J,J,SNESDefaultComputeJacobian_Block, PETSC_NULL );
	}
	
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Customize nonlinear solver; set runtime options
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	
	/* 
	Set linear solver defaults for this problem. By extracting the
	KSP, KSP, and PC contexts from the SNES context, we can then
	directly call any KSP, KSP, and PC routines to set various options.
	*/
	ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
	ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
	ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);
	ierr = KSPSetTolerances(ksp,1.e-4,PETSC_DEFAULT,PETSC_DEFAULT,20);CHKERRQ(ierr);
	
	/* 
	Set SNES/KSP/KSP/PC runtime options, e.g.,
	-snes_view -snes_monitor -ksp_type <ksp> -pc_type <pc>
	These options will override those specified above as long as
	SNESSetFromOptions() is called _after_ any other customization
	routines.
	*/
	ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
	ierr = SNESSetFromOptions_Block(snes);CHKERRQ(ierr);
	
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Evaluate initial guess; then solve nonlinear system
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	if (!flg) {
		ierr = VecSet(x,pfive);CHKERRQ(ierr);
	} else {
		VecBlockGetSubVector( x, 0, &bv );
		VecSetValue( bv, 0, 2.0, INSERT_VALUES );	//xx[0] = 2.0; 
		VecAssemblyBegin(bv);	VecAssemblyEnd(bv);
		
		VecBlockGetSubVector( x, 1, &bv );
		VecSetValue( bv, 0, 3.0, INSERT_VALUES );	//xx[1] = 3.0; 
		VecAssemblyBegin(bv);	VecAssemblyEnd(bv);
		
		VecBlockRestoreSubVectors( x );
	}
	/*
	Note: The user should initialize the vector, x, with the initial guess
	for the nonlinear solver prior to calling SNESSolve().  In particular,
	to employ an initial guess of zero, the user should explicitly set
	this vector to zero by calling VecSet().
	*/
	ierr = SNESSolve(snes,PETSC_NULL,x);CHKERRQ(ierr);
	ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
	if (flg) {
		Vec f;
		ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
		ierr = SNESGetFunction(snes,&f,0,0);CHKERRQ(ierr);
		ierr = VecView(r,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	}
	
	ierr = PetscPrintf(PETSC_COMM_SELF,"number of Newton iterations = %D\n\n",its);CHKERRQ(ierr);
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Free work space.  All PETSc objects should be destroyed when they
	are no longer needed.
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	
	ierr = Stg_VecDestroy( &x);CHKERRQ(ierr); ierr = Stg_VecDestroy( &r);CHKERRQ(ierr);
	ierr = Stg_MatDestroy( &J);CHKERRQ(ierr); ierr = Stg_SNESDestroy( &snes);CHKERRQ(ierr);
	
	return 0;
}
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormFunction1_block"
PetscErrorCode FormFunction1_block(SNES snes,Vec x,Vec f,void *dummy)
{
	Vec *xx, *ff, x1,x2, f1,f2;
	PetscScalar ff_0, ff_1;
	PetscScalar xx_0, xx_1;
	PetscInt index;
	
	/* get blocks for function */
	VecBlockGetSubVectors( f, &ff );
	f1 = ff[0];	f2 = ff[1];
//	VecBlockGetValue( f, 0, &f1 );
//	VecBlockGetValue( f, 1, &f2 );
	
	/* get blocks for solution */
	VecBlockGetSubVectors( x, &xx );
	x1 = xx[0];	x2 = xx[1];
//	VecBlockGetValue( x, 0, &x1 );
//	VecBlockGetValue( x, 1, &x2 );
	
	/* get solution values */
	index = 0;
	VecGetValues( x1,1, &index, &xx_0 );
	VecGetValues( x2,1, &index, &xx_1 );
	
	/* Compute function */
	ff_0 = xx_0*xx_0 + xx_0*xx_1 - 3.0;
	ff_1 = xx_0*xx_1 + xx_1*xx_1 - 6.0;
	
	/* set function values */
	VecSetValue( f1, index, ff_0, INSERT_VALUES ); 
	VecAssemblyBegin(f1);	VecAssemblyEnd(f1);
	
	VecSetValue( f2, index, ff_1, INSERT_VALUES );
	VecAssemblyBegin(f2);	VecAssemblyEnd(f2);
//	PetscPrintf( PETSC_COMM_WORLD, "\t%f %f \n", ff_0, ff_1 );


	VecBlockRestoreSubVectors( f );
	VecBlockRestoreSubVectors( x );

	
	return 0;
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormFunction2_block"
PetscErrorCode FormFunction2_block(SNES snes,Vec x,Vec f,void *dummy)
{
	Vec *xx, *ff, x1,x2, f1,f2;
	PetscScalar ff_0, ff_1;
	PetscScalar xx_0, xx_1;
	PetscInt index;
	
	/* get blocks for function */
	VecBlockGetSubVectors( f, &ff );
	f1 = ff[0];
	f2 = ff[1];
//	VecBlockGetValue( f, 0, &f1 );
//	VecBlockGetValue( f, 1, &f2 );
	
	/* get blocks for solution */
	VecBlockGetSubVectors( x, &xx );
	x1 = xx[0];	
	x2 = xx[1];
	
	/* get solution values */
	index = 0;
	VecGetValues( x1,1, &index, &xx_0 );
	VecGetValues( x2,1, &index, &xx_1 );
	
	/* Compute function */
	ff_0 = PetscSinScalar( 3.0*xx_0) + xx_0;
	ff_1 = xx_1;
	
	/* set function values */
	VecSetValue( f1, index, ff_0, INSERT_VALUES ); 
	VecAssemblyBegin(f1);	VecAssemblyEnd(f1);
	
	VecSetValue( f2, index, ff_1, INSERT_VALUES );
	VecAssemblyBegin(f2);	VecAssemblyEnd(f2);
//	PetscPrintf( PETSC_COMM_WORLD, "\t%f %f \n", ff_0, ff_1 );


	VecBlockRestoreSubVectors( f );
	VecBlockRestoreSubVectors( x );

	
	return 0;
}

/*

Notes,
  Using the options -mat_fd_type ds seems to produce closer answers than the default which is wp.


./test-BlockSNESMFFD.app -snes_block_fd -ass_snes_fd -mat_fd_type ds

*/


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	PetscMPIInt    size;
	PetscErrorCode ierr;
	
	PetscInitialize(&argc,&argv,(char *)0,help);
	
	PetscExtVecRegisterAll();
	PetscExtMatRegisterAll();
	
	
	
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	if (size != 1) Stg_SETERRQ(1,"This is a uniprocessor example only!");
	
	/* force snes monitor to be on */
	PetscOptionsInsertString( "-snes_monitor" );
	
	
	block_system_force(argc,argv);
	
	
	PetscExtVecRegisterDestroyAll();
	PetscExtMatRegisterDestroyAll();
	
	
	ierr = PetscFinalize();CHKERRQ(ierr);
	return 0;
}
