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

/*T
Concepts: SNES^basic uniprocessor example, block objects
Processors: 1
T*/

/* 
Include "petscsnes.h" so that we can use SNES solvers.  Note that this
file automatically includes:
petsc.h       - base PETSc routines   petscvec.h - vectors
petscsys.h    - system routines       petscmat.h - matrices
petscis.h     - index sets            petscksp.h - Krylov subspace methods
petscviewer.h - viewers               petscpc.h  - preconditioners
petscksp.h   - linear solvers
petscvec-block.h - block definition for vectors
petscmat-block.h - block definition for matrics
*/

/*
This example is block version of the test found at 
	${PETSC_DIR}/src/snes/examples/tutorials/ex1.c
In this test we replace the Jacobian systems
	[J]{x} = {F} 
where

[J] = ( j_00, j_01 ),	{x} = ( x_0, x_1 )^T, 	{F} = ( f_0, f_1 )^T
      ( j_10, j_11 )			
where [J] \in \mathbb^{2 \times 2}, {x},{F} \in \mathbb^{2 \times 1},

with a block system in which each block is of length 1.
i.e. The block system is thus

[J] = ( [j00], [j01] ),	{x} = ( {x0}, {x1} )^T, {F} = ( {f0}, {f1} )^T
      ( [j10], [j11] ) 
where 
[j00], [j01], [j10], [j11] \in \mathbb^{1 \times 1}
{x0}, {x1}, {f0}, {f1} \in \mathbb^{1 \times 1}

In practise we would not bother defing blocks of size one, and would instead assemble the 
full system. This is just a simple test to illustrate how to manipulate the blocks and
to confirm the implementation is correct.
*/



#include "petscsnes.h"
#include "petscext_vec.h"
#include "petscext_mat.h"

/* 
User-defined routines
*/
extern PetscErrorCode FormJacobian1(SNES,Vec,Mat*,Mat*,MatStructure*,void*);
extern PetscErrorCode FormFunction1(SNES,Vec,Vec,void*);
extern PetscErrorCode FormJacobian2(SNES,Vec,Mat*,Mat*,MatStructure*,void*);
extern PetscErrorCode FormFunction2(SNES,Vec,Vec,void*);

extern PetscErrorCode FormJacobian1_block(SNES,Vec,Mat*,Mat*,MatStructure*,void*);
extern PetscErrorCode FormFunction1_block(SNES,Vec,Vec,void*);
extern PetscErrorCode FormJacobian2_block(SNES,Vec,Mat*,Mat*,MatStructure*,void*);
extern PetscErrorCode FormFunction2_block(SNES,Vec,Vec,void*);


#if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >=3 ) )
  #include "petsc-private/snesimpl.h"
#else
  #include "private/snesimpl.h"
#endif

PetscErrorCode PETSCSNES_DLLEXPORT SNESBlockDefaultComputeJacobian(SNES snes,Vec x1,Mat *J,Mat *B,MatStructure *flag,void *ctx)
{
	Vec            j1a,j2a,x2;
	PetscErrorCode ierr;
	PetscInt       i,N,start,end,j,value, II,JJ, BM,BN,s,e;
	PetscScalar    dx,*y,scale,*xx,wscale;
	PetscReal      amax,epsilon = PETSC_SQRT_MACHINE_EPSILON;
	PetscReal      dx_min = 1.e-16,dx_par = 1.e-1,unorm;
	MPI_Comm       comm;
	PetscErrorCode (*eval_fct)(SNES,Vec,Vec)=0;
	PetscTruth     assembled,use_wp = PETSC_TRUE,flg;
	const char     *list[2] = {"ds","wp"};
	Vec sub_x1, sub_x2, sub_j2a;
	Mat **BB, **BkJac;
	Vec *sub_x1s, *sub_x2s, *sub_j2as;
	const char *prefix;
	
	PetscFunctionBegin;
	SNESGetOptionsPrefix( snes, &prefix );
	ierr = PetscOptionsGetReal(prefix,"-snes_test_err",&epsilon,0);CHKERRQ(ierr);
	eval_fct = SNESComputeFunction;
	
	ierr = PetscObjectGetComm((PetscObject)x1,&comm);CHKERRQ(ierr);
	ierr = MatAssembled(*B,&assembled);CHKERRQ(ierr);
	if (assembled) {
		ierr = MatZeroEntries(*B);CHKERRQ(ierr);
	}
	if (!snes->nvwork) {
		ierr = VecDuplicateVecs(x1,3,&snes->vwork);CHKERRQ(ierr);
		snes->nvwork = 3;
		ierr = PetscLogObjectParents(snes,3,snes->vwork);CHKERRQ(ierr);
	}
	j1a = snes->vwork[0]; j2a = snes->vwork[1]; x2 = snes->vwork[2];
	
	
	ierr = VecGetSize(x1,&BN);CHKERRQ(ierr);
	ierr = (*eval_fct)(snes,x1,j1a);CHKERRQ(ierr);
	
	ierr = PetscOptionsEList("-mat_fd_type","Algorithm to compute difference parameter","SNESDefaultComputeJacobian",list,2,"wp",&value,&flg);CHKERRQ(ierr);
	if (flg && !value) {
		use_wp = PETSC_FALSE;
	}
	if (use_wp) {
		ierr = VecNorm(x1,NORM_2,&unorm);CHKERRQ(ierr);
	}
	
	VecBlockGetSubVectors( x1, &sub_x1s );
	VecBlockGetSubVectors( x2, &sub_x2s );
	MatBlockGetSubMatrices( *B, &BB );
	MatBlockGetSubMatrices( *J, &BkJac );
	
	for( II=0; II<BN; II++ ) {
		sub_x1 = sub_x1s[II];
		sub_x2 = sub_x2s[II];
		
		ierr = VecGetSize(sub_x1,&N);CHKERRQ(ierr);
		ierr = VecGetOwnershipRange(sub_x1,&start,&end);CHKERRQ(ierr);
		
		/* Compute Jacobian approximation, 1 column at a time. 
		x1 = current iterate, j1a = F(x1)
		x2 = perturbed iterate, j2a = F(x2)
		*/
		for (i=0; i<N; i++) {
			ierr = VecCopy(x1,x2);CHKERRQ(ierr);
			if (i>= start && i<end) {
				ierr = VecGetArray(sub_x1,&xx);CHKERRQ(ierr);
				if (use_wp) {
					dx = 1.0 + unorm;
				} else {
					dx = xx[i-start];
				}
				ierr = VecRestoreArray(sub_x1,&xx);CHKERRQ(ierr);
				
				
#if !defined(PETSC_USE_COMPLEX)
				if (dx < dx_min && dx >= 0.0) dx = dx_par;
				else if (dx < 0.0 && dx > -dx_min) dx = -dx_par;
#else
				if (PetscAbsScalar(dx) < dx_min && PetscRealPart(dx) >= 0.0) dx = dx_par;
				else if (PetscRealPart(dx) < 0.0 && PetscAbsScalar(dx) < dx_min) dx = -dx_par;
#endif
				
				dx *= epsilon;
				wscale = 1.0/dx;
				ierr = VecSetValues(sub_x2,1,&i,&dx,ADD_VALUES);CHKERRQ(ierr);
			} else {
				wscale = 0.0;
			}
			
			ierr = (*eval_fct)(snes,x2,j2a);CHKERRQ(ierr);
			ierr = VecAXPY(j2a,-1.0,j1a);CHKERRQ(ierr);
			/* Communicate scale to all processors */
			ierr = MPI_Allreduce(&wscale,&scale,1,MPIU_SCALAR,MPI_SUM,comm);CHKERRQ(ierr);
			
			ierr = VecScale(j2a,scale);CHKERRQ(ierr);
			ierr = VecNorm(j2a,NORM_INFINITY,&amax);CHKERRQ(ierr); amax *= 1.e-14;
			
			VecBlockGetSubVectors( j2a, &sub_j2as );
			for( JJ=0; JJ<BN; JJ++ ) {
				sub_j2a = sub_j2as[JJ];
				
				ierr = VecGetOwnershipRange(sub_j2a,&s,&e);CHKERRQ(ierr);
				ierr = VecGetArray(sub_j2a,&y);CHKERRQ(ierr);
				for (j=s; j<e; j++) {
					if (PetscAbsScalar(y[j-s]) > amax) {
					//	ierr = MatSetValues(*B,1,&j,1,&i,y+j-start,INSERT_VALUES);CHKERRQ(ierr);
						
						ierr = MatSetValues(BB[JJ][II],1,&j,1,&i,y+j-s,INSERT_VALUES);CHKERRQ(ierr);
						
					}
				}
				ierr = VecRestoreArray(sub_j2a,&y);CHKERRQ(ierr);
			}
			VecBlockRestoreSubVectors( j2a );
			
		}
	
	}
	
	VecBlockRestoreSubVectors( x1 );
	VecBlockRestoreSubVectors( x2 );
	
	MatGetSize( *B, &BM, &BN );
	for( II=0; II<BM; II++ ) {
		for( JJ=0; JJ<BN; JJ++ ) {
			if( BB[II][JJ] != PETSC_NULL ) {
				ierr  = MatAssemblyBegin(BB[II][JJ],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
				ierr  = MatAssemblyEnd(BB[II][JJ],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
			}
			if (BB[II][JJ] != BkJac[II][JJ]) {
				if( BkJac[II][JJ] != PETSC_NULL ) {
					ierr  = MatAssemblyBegin(BkJac[II][JJ],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
					ierr  = MatAssemblyEnd(BkJac[II][JJ],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
				}
			}
		}
	}
	
	MatBlockRestoreSubMatrices( *B );
	MatBlockRestoreSubMatrices( *J );
	
	ierr  = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr  = MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	if (*B != *J) {
		ierr  = MatAssemblyBegin(*J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr  = MatAssemblyEnd(*J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}
	*flag =  DIFFERENT_NONZERO_PATTERN;
	
	
//	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
//	MatView( *B, PETSC_VIEWER_STDOUT_WORLD );
	
	
	PetscFunctionReturn(0);
}

PetscErrorCode SNESSetFromOptions_BlockFD( SNES snes )
{
	PetscTruth flg;
	PetscErrorCode ierr;
    PetscFunctionBegin;
	
	//PetscOptionsBegin( ((PetscObject)snes)->comm, ((PetscObject)snes)->prefix,"Nonlinear solver (SNES) options for block matrices", "SNES" );
	ierr = PetscOptionsHead("Nonlinear solver (SNES) options for block matrices");CHKERRQ(ierr);
	
	PetscOptionsName("-snes_block_fd","Use block finite difference (slow) to compute Jacobian", "SNESBlockDefaultComputeJacobian", &flg );
	if (flg) {
		SNESSetJacobian( snes, snes->jacobian, snes->jacobian_pre, SNESBlockDefaultComputeJacobian, PETSC_NULL) );
		PetscInfo( snes, "Setting block finite difference Jacobian matrix \n");
	}
	
    //PetscOptionsEnd();
	ierr = PetscOptionsTail();CHKERRQ(ierr);
	
    PetscFunctionReturn(0);
}



PetscErrorCode _SNESDefaultComputeJacobian(SNES snes,Vec x1,Mat *J,Mat *B,MatStructure *flag,void *ctx)
{
	Vec            j1a,j2a,x2;
	PetscErrorCode ierr;
	PetscInt       i,N,start,end,j,value;
	PetscScalar    dx,*y,scale,*xx,wscale;
	PetscReal      amax,epsilon = PETSC_SQRT_MACHINE_EPSILON;
	PetscReal      dx_min = 1.e-16,dx_par = 1.e-1,unorm;
	MPI_Comm       comm;
	PetscErrorCode (*eval_fct)(SNES,Vec,Vec)=0;
	PetscTruth     assembled,use_wp = PETSC_TRUE,flg;
	const char     *list[2] = {"ds","wp"};
	
	PetscFunctionBegin;
	ierr = PetscOptionsGetReal( ((PetscObject)snes)->prefix,"-snes_test_err",&epsilon,0);CHKERRQ(ierr);
	eval_fct = SNESComputeFunction;
	
	ierr = PetscObjectGetComm((PetscObject)x1,&comm);CHKERRQ(ierr);
	ierr = MatAssembled(*B,&assembled);CHKERRQ(ierr);
	if (assembled) {
		ierr = MatZeroEntries(*B);CHKERRQ(ierr);
	}
	if (!snes->nvwork) {
		ierr = VecDuplicateVecs(x1,3,&snes->vwork);CHKERRQ(ierr);
		snes->nvwork = 3;
		ierr = PetscLogObjectParents(snes,3,snes->vwork);CHKERRQ(ierr);
	}
	j1a = snes->vwork[0]; j2a = snes->vwork[1]; x2 = snes->vwork[2];
	
	ierr = VecGetSize(x1,&N);CHKERRQ(ierr);
	ierr = VecGetOwnershipRange(x1,&start,&end);CHKERRQ(ierr);
	ierr = (*eval_fct)(snes,x1,j1a);CHKERRQ(ierr);
	
	ierr = PetscOptionsEList("-mat_fd_type","Algorithm to compute difference parameter","SNESDefaultComputeJacobian",list,2,"wp",&value,&flg);CHKERRQ(ierr);
	if (flg && !value) {
		use_wp = PETSC_FALSE;
	}
	if (use_wp) {
		ierr = VecNorm(x1,NORM_2,&unorm);CHKERRQ(ierr);
	}
	
	/* Compute Jacobian approximation, 1 column at a time. 
	x1 = current iterate, j1a = F(x1)
	x2 = perturbed iterate, j2a = F(x2)
	*/
	for (i=0; i<N; i++) {
		ierr = VecCopy(x1,x2);CHKERRQ(ierr);
		if (i>= start && i<end) {
			ierr = VecGetArray(x1,&xx);CHKERRQ(ierr);
			if (use_wp) {
				dx = 1.0 + unorm;
			} else {
				dx = xx[i-start];
			}
			ierr = VecRestoreArray(x1,&xx);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
			if (dx < dx_min && dx >= 0.0) dx = dx_par;
			else if (dx < 0.0 && dx > -dx_min) dx = -dx_par;
#else
			if (PetscAbsScalar(dx) < dx_min && PetscRealPart(dx) >= 0.0) dx = dx_par;
			else if (PetscRealPart(dx) < 0.0 && PetscAbsScalar(dx) < dx_min) dx = -dx_par;
#endif
			dx *= epsilon;
			wscale = 1.0/dx;
			ierr = VecSetValues(x2,1,&i,&dx,ADD_VALUES);CHKERRQ(ierr);
		} else {
			wscale = 0.0;
		}
		ierr = (*eval_fct)(snes,x2,j2a);CHKERRQ(ierr);
		ierr = VecAXPY(j2a,-1.0,j1a);CHKERRQ(ierr);
		/* Communicate scale to all processors */
		ierr = MPI_Allreduce(&wscale,&scale,1,MPIU_SCALAR,MPI_SUM,comm);CHKERRQ(ierr);
		
		ierr = VecScale(j2a,scale);CHKERRQ(ierr);
		ierr = VecNorm(j2a,NORM_INFINITY,&amax);CHKERRQ(ierr); amax *= 1.e-14;
		
		
		ierr = VecGetArray(j2a,&y);CHKERRQ(ierr);
		for (j=start; j<end; j++) {
			if (PetscAbsScalar(y[j-start]) > amax) {
				ierr = MatSetValues(*B,1,&j,1,&i,y+j-start,INSERT_VALUES);CHKERRQ(ierr);
			}
		}
		ierr = VecRestoreArray(j2a,&y);CHKERRQ(ierr);
	}
	ierr  = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr  = MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	if (*B != *J) {
		ierr  = MatAssemblyBegin(*J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr  = MatAssemblyEnd(*J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}
	*flag =  DIFFERENT_NONZERO_PATTERN;
	
//	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
//	MatView( *B, PETSC_VIEWER_STDOUT_WORLD );
	
	
	
	PetscFunctionReturn(0);
}

PetscErrorCode _SNESSetFromOptions_FDDebug( SNES snes )
{
	PetscTruth flg;
	PetscErrorCode ierr;
    PetscFunctionBegin;
	
	//PetscOptionsBegin( ((PetscObject)snes)->comm, ((PetscObject)snes)->prefix,"Nonlinear solver (SNES) options for debug", "SNES" );
	ierr = PetscOptionsHead("Nonlinear solver (SNES) options for debug");CHKERRQ(ierr);
	
    PetscOptionsName("-snes_fd","Use block finite difference (slow) to compute Jacobian", "_SNESDefaultComputeJacobian", &flg );
	if (flg) {
		SNESSetJacobian( snes, snes->jacobian, snes->jacobian_pre, _SNESDefaultComputeJacobian, PETSC_NULL) );
		PetscInfo( snes, "Setting block finite difference Jacobian matrix \n");
	}
	
    
	//PetscOptionsEnd();
	ierr = PetscOptionsTail();CHKERRQ(ierr);

	PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "assembled_system"
int assembled_system(int argc,char **argv)
{
	SNES           snes;         /* nonlinear solver context */
	KSP            ksp;         /* linear solver context */
	PC             pc;           /* preconditioner context */
	Vec            x,r;         /* solution, residual vectors */
	Mat            J;            /* Jacobian matrix */
	PetscErrorCode ierr;
	PetscInt       its;
	PetscScalar    pfive = .5,*xx;
	PetscTruth     flg;
	
	
	PetscPrintf( PETSC_COMM_WORLD, "\n\n========================= Assembled system =========================\n\n" );
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Create nonlinear solver context
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	
	ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
	SNESAppendOptionsPrefix( snes, "ass_" );
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Create matrix and vector data structures; set corresponding routines
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	
	/*
	Create vectors for solution and nonlinear function
	*/
	ierr = VecCreateSeq(PETSC_COMM_SELF,2,&x);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&r);CHKERRQ(ierr);
	
	/*
	Create Jacobian matrix data structure
	*/
	ierr = MatCreate(PETSC_COMM_SELF,&J);CHKERRQ(ierr);
	ierr = MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,2,2);CHKERRQ(ierr);
	ierr = MatSetFromOptions(J);CHKERRQ(ierr);
	
	ierr = PetscOptionsHasName(PETSC_NULL,"-hard",&flg);CHKERRQ(ierr);
	if (!flg) {
		/* 
		Set function evaluation routine and vector.
		*/
		ierr = SNESSetFunction(snes,r,FormFunction1,PETSC_NULL);CHKERRQ(ierr);
		
		/* 
		Set Jacobian matrix data structure and Jacobian evaluation routine
		*/
		ierr = SNESSetJacobian(snes,J,J,FormJacobian1,PETSC_NULL);CHKERRQ(ierr);
	} else {
		ierr = SNESSetFunction(snes,r,FormFunction2,PETSC_NULL);CHKERRQ(ierr);
		ierr = SNESSetJacobian(snes,J,J,FormJacobian2,PETSC_NULL);CHKERRQ(ierr);
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
	_SNESSetFromOptions_FDDebug( snes );
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Evaluate initial guess; then solve nonlinear system
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	if (!flg) {
		ierr = VecSet(x,pfive);CHKERRQ(ierr);
	} else {
		ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
		xx[0] = 2.0; xx[1] = 3.0;
		ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);
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
#define __FUNCT__ "FormFunction1"
/* 
FormFunction1 - Evaluates nonlinear function, F(x).

Input Parameters:
.  snes - the SNES context
.  x - input vector
.  dummy - optional user-defined context (not used here)

Output Parameter:
.  f - function vector
*/
PetscErrorCode FormFunction1(SNES snes,Vec x,Vec f,void *dummy)
{
	PetscErrorCode ierr;
	PetscScalar    *xx,*ff;
	
	/*
	Get pointers to vector data.
	- For default PETSc vectors, VecGetArray() returns a pointer to
	the data array.  Otherwise, the routine is implementation dependent.
	- You MUST call VecRestoreArray() when you no longer need access to
	the array.
	*/
	ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
	ierr = VecGetArray(f,&ff);CHKERRQ(ierr);
	
	/*
	Compute function
	*/
	ff[0] = xx[0]*xx[0] + xx[0]*xx[1] - 3.0;
	ff[1] = xx[0]*xx[1] + xx[1]*xx[1] - 6.0;
	
	//PetscPrintf( PETSC_COMM_WORLD, "\t%f %f \n", ff[0], ff[1] );
	
	
	/*
	Restore vectors
	*/
	ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);
	ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr); 
	
	return 0;
}
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormJacobian1"
/*
FormJacobian1 - Evaluates Jacobian matrix.

Input Parameters:
.  snes - the SNES context
.  x - input vector
.  dummy - optional user-defined context (not used here)

Output Parameters:
.  jac - Jacobian matrix
.  B - optionally different preconditioning matrix
.  flag - flag indicating matrix structure
*/
PetscErrorCode FormJacobian1(SNES snes,Vec x,Mat *jac,Mat *B,MatStructure *flag,void *dummy)
{
	PetscScalar    *xx,A[4];
	PetscErrorCode ierr;
	PetscInt       idx[2] = {0,1};
	
	/*
	Get pointer to vector data
	*/
	ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
	
	/*
	Compute Jacobian entries and insert into matrix.
	- Since this is such a small problem, we set all entries for
	the matrix at once.
	*/
	A[0] = 2.0*xx[0] + xx[1]; A[1] = xx[0];
	A[2] = xx[1]; A[3] = xx[0] + 2.0*xx[1];
	ierr = MatSetValues(*jac,2,idx,2,idx,A,INSERT_VALUES);CHKERRQ(ierr);
	*flag = SAME_NONZERO_PATTERN;
	
	/*
	Restore vector
	*/
	ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);
	
	/*
	Assemble matrix
	*/
	ierr = MatAssemblyBegin(*jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	return 0;
}


/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormFunction2"
PetscErrorCode FormFunction2(SNES snes,Vec x,Vec f,void *dummy)
{
	PetscErrorCode ierr;
	PetscScalar    *xx,*ff;
	
	/*
	Get pointers to vector data.
	- For default PETSc vectors, VecGetArray() returns a pointer to
	the data array.  Otherwise, the routine is implementation dependent.
	- You MUST call VecRestoreArray() when you no longer need access to
	the array.
	*/
	ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
	ierr = VecGetArray(f,&ff);CHKERRQ(ierr);
	
	/*
	Compute function
	*/
	ff[0] = PetscSinScalar(3.0*xx[0]) + xx[0];
	ff[1] = xx[1];
	
	/*
	Restore vectors
	*/
	ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);
	ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr); 
	
	return 0;
}
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormJacobian2"
PetscErrorCode FormJacobian2(SNES snes,Vec x,Mat *jac,Mat *B,MatStructure *flag,void *dummy)
{
	PetscScalar    *xx,A[4];
	PetscErrorCode ierr;
	PetscInt       idx[2] = {0,1};
	
	/*
	Get pointer to vector data
	*/
	ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
	
	/*
	Compute Jacobian entries and insert into matrix.
	- Since this is such a small problem, we set all entries for
	the matrix at once.
	*/
	A[0] = 3.0*PetscCosScalar(3.0*xx[0]) + 1.0; A[1] = 0.0;
	A[2] = 0.0;                     A[3] = 1.0;
	ierr = MatSetValues(*jac,2,idx,2,idx,A,INSERT_VALUES);CHKERRQ(ierr);
	*flag = SAME_NONZERO_PATTERN;
	
	/*
	Restore vector
	*/
	ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);
	
	/*
	Assemble matrix
	*/
	ierr = MatAssemblyBegin(*jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	return 0;
}







#undef __FUNCT__
#define __FUNCT__ "block_system"
int block_system(int argc,char **argv)
{
	SNES           snes;         /* nonlinear solver context */
	KSP            ksp;         /* linear solver context */
	PC             pc;           /* preconditioner context */
	Vec            x,r;         /* solution, residual vectors */
	Mat            J,pJ;            /* Jacobian matrix */
	PetscErrorCode ierr;
	PetscInt       its;
	PetscScalar    pfive = .5;
	PetscTruth     flg;
	
	Mat pj11, pj12, pj21, pj22;
	Mat j11, j12, j21, j22;
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
	
	
	/*
	Create sub Jacobian matrix data structure
	*/
	MatCreate( PETSC_COMM_WORLD, &j11 );
	MatSetSizes( j11, PETSC_DECIDE, PETSC_DECIDE, 1, 1 );
	MatSetType( j11, MATSEQAIJ );
	
	MatCreate( PETSC_COMM_WORLD, &j12 );
	MatSetSizes( j12, PETSC_DECIDE, PETSC_DECIDE, 1, 1 );
	MatSetType( j12, MATSEQAIJ );
	
	MatCreate( PETSC_COMM_WORLD, &j21 );
	MatSetSizes( j21, PETSC_DECIDE, PETSC_DECIDE, 1, 1 );
	MatSetType( j21, MATSEQAIJ );
	
	MatCreate( PETSC_COMM_WORLD, &j22 );
	MatSetSizes( j22, PETSC_DECIDE, PETSC_DECIDE, 1, 1 );
	MatSetType( j22, MATSEQAIJ );

	
	
	/*
	Create block Jacobian matrix data structure
	*/
	MatCreate( PETSC_COMM_WORLD, &J );
	MatSetSizes_Block( J, 2,2, 2,2 );
	MatSetType( J, "block" );
	
	MatBlockSetValue( J, 0,0, j11, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( J, 0,1, j12, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( J, 1,0, j21, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( J, 1,1, j22, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	
	MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);
	
	//
	MatCreate( PETSC_COMM_WORLD, &pj11 );
	MatSetSizes( pj11, PETSC_DECIDE, PETSC_DECIDE, 1, 1 );
	MatSetType( pj11, MATSEQAIJ );
	
	MatCreate( PETSC_COMM_WORLD, &pj12 );
	MatSetSizes( pj12, PETSC_DECIDE, PETSC_DECIDE, 1, 1 );
	MatSetType( pj12, MATSEQAIJ );
	
	MatCreate( PETSC_COMM_WORLD, &pj21 );
	MatSetSizes( pj21, PETSC_DECIDE, PETSC_DECIDE, 1, 1 );
	MatSetType( pj21, MATSEQAIJ );
	
	MatCreate( PETSC_COMM_WORLD, &pj22 );
	MatSetSizes( pj22, PETSC_DECIDE, PETSC_DECIDE, 1, 1 );
	MatSetType( pj22, MATSEQAIJ );
	
	MatCreate( PETSC_COMM_WORLD, &pJ );
	MatSetSizes_Block( pJ, 2,2, 2,2 );
	MatSetType( pJ, "block" );
	
	MatBlockSetValue( pJ, 0,0, pj11, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( pJ, 0,1, pj12, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( pJ, 1,0, pj21, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( pJ, 1,1, pj22, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	
	MatAssemblyBegin(pJ,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(pJ,MAT_FINAL_ASSEMBLY);
	//
	
	
	ierr = PetscOptionsHasName(PETSC_NULL,"-hard",&flg);CHKERRQ(ierr);
	if (!flg) {
		/* Set function evaluation routine and vector. */
		ierr = SNESSetFunction(snes,r,FormFunction1_block,PETSC_NULL);CHKERRQ(ierr);
		
		/* Set Jacobian matrix data structure and Jacobian evaluation routine */
		ierr = SNESSetJacobian(snes,J,J,FormJacobian1_block,PETSC_NULL);CHKERRQ(ierr);
	} else {
		ierr = SNESSetFunction(snes,r,FormFunction2_block,PETSC_NULL);CHKERRQ(ierr);
		ierr = SNESSetJacobian(snes,J,J,FormJacobian2_block,PETSC_NULL);CHKERRQ(ierr);
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
	ierr = SNESSetFromOptions_BlockFD(snes);CHKERRQ(ierr);

	
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
#define __FUNCT__ "FormJacobian1_block"
PetscErrorCode FormJacobian1_block(SNES snes,Vec x,Mat *jac,Mat *B,MatStructure *flag,void *dummy)
{
	Vec *xx, x1,x2;
	PetscScalar xx_0, xx_1;
	PetscInt index;
	PetscScalar A_00, A_01, A_10, A_11;
	Mat j11, j12, j21, j22;
	
	/* get blocks for solution */
	VecBlockGetSubVectors( x, &xx );
	x1 = xx[0];	x2 = xx[1];
//	VecBlockGetValue( x, 0, &x1 );
//	VecBlockGetValue( x, 1, &x2 );
	
	/* get solution values */
	index = 0;
	VecGetValues( x1,1, &index, &xx_0 );
	VecGetValues( x2,1, &index, &xx_1 );
	
	/* get block matrices */
	MatBlockGetSubMatrix( *jac, 0,0, &j11 );
	MatBlockGetSubMatrix( *jac, 0,1, &j12 );
	MatBlockGetSubMatrix( *jac, 1,0, &j21 );
	MatBlockGetSubMatrix( *jac, 1,1, &j22 );
	
//	MatBlockGetValues( *jac, 2, idx, 2, idx, mats );
	
	/* compute jacobian entries */
	A_00 = 2.0*xx_0 + xx_1;
	A_01 = xx_0;
	A_10 = xx_1;
	A_11 = xx_0 + 2.0*xx_1;
	
	/* set jacobian values */
	MatSetValue( j11, 0,0, A_00, INSERT_VALUES);
	MatSetValue( j12, 0,0, A_01, INSERT_VALUES);
	MatSetValue( j21, 0,0, A_10, INSERT_VALUES);
	MatSetValue( j22, 0,0, A_11, INSERT_VALUES);
	
	
	
	*flag = SAME_NONZERO_PATTERN;
	
	/* Assemble sub matrix */
	MatAssemblyBegin(j11,MAT_FINAL_ASSEMBLY); MatAssemblyEnd(j11,MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(j12,MAT_FINAL_ASSEMBLY); MatAssemblyEnd(j12,MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(j21,MAT_FINAL_ASSEMBLY); MatAssemblyEnd(j21,MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(j22,MAT_FINAL_ASSEMBLY); MatAssemblyEnd(j22,MAT_FINAL_ASSEMBLY);

	VecBlockRestoreSubVectors( x );
	MatBlockRestoreSubMatrices( *jac );	

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
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormJacobian2_block"
PetscErrorCode FormJacobian2_block(SNES snes,Vec x,Mat *jac,Mat *B,MatStructure *flag,void *dummy)
{
	PetscScalar    *xx,A[4];
	PetscErrorCode ierr;
	PetscInt       idx[2] = {0,1};
	
	/*
	Get pointer to vector data
	*/
	ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
	
	/*
	Compute Jacobian entries and insert into matrix.
	- Since this is such a small problem, we set all entries for
	the matrix at once.
	*/
	A[0] = 3.0*PetscCosScalar(3.0*xx[0]) + 1.0; A[1] = 0.0;
	A[2] = 0.0;                     A[3] = 1.0;
	ierr = MatSetValues(*jac,2,idx,2,idx,A,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues(*B,2,idx,2,idx,A,INSERT_VALUES);CHKERRQ(ierr);
	*flag = SAME_NONZERO_PATTERN;
	
	/*
	Restore vector
	*/
	ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);
	
	/*
	Assemble matrix
	*/
	ierr = MatAssemblyBegin(*jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
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
	PetscOptionsInsertString( "-snes_monitor -ass_snes_monitor" );
	
	assembled_system(argc,argv);
	block_system(argc,argv);
	
	PetscExtVecRegisterDestroyAll();
	PetscExtMatRegisterDestroyAll();
	
	
	ierr = PetscFinalize();CHKERRQ(ierr);
	return 0;
}
