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


#include "petscsnes.h"
#include "petscext_vec.h"
#include "petscext_mat.h"

#include "private/snesimpl.h"

PetscErrorCode SNESDefaultComputeJacobian_Block(SNES snes,Vec x1,Mat *J,Mat *B,MatStructure *flag,void *ctx)
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

PetscErrorCode SNESSetFromOptions_Block( SNES snes )
{
	PetscTruth flg;
	PetscErrorCode ierr;
	
	ierr=PetscOptionsBegin( ((PetscObject)snes)->comm, ((PetscObject)snes)->prefix,"Nonlinear solver (SNES) options for block matrices", "SNES" );CHKERRQ(ierr);
	
	PetscOptionsName("-snes_fd","Use finite difference (slow) to compute Jacobian on block matrices", "SNESDefaultComputeJacobian_Block", &flg );
	if (flg) {
		SNESSetJacobian( snes, snes->jacobian, snes->jacobian_pre, SNESDefaultComputeJacobian_Block, snes->funP );
		PetscInfo( snes, "Setting block finite difference Jacobian matrix \n");
	}
	
	PetscOptionsEnd();
	
	PetscFunctionReturn(0);
}

PetscErrorCode MatCreateSNESBlockMFFD(SNES snes,Mat *J)
{
	PetscErrorCode ierr;
	PetscInt i,j;
	Mat subMat;
	PetscInt Mr,*M;
	PetscTruth is_block;
	Vec *sub_vecs;
	MPI_Comm comm;
	
	PetscFunctionBegin;
	
	
	if( !snes->vec_func) Stg_SETERRQ(PETSC_ERR_ARG_WRONGSTATE, "SNESSetFunction() must be called first");
	PetscTypeCompare( (PetscObject)snes->vec_func, "block", &is_block );
	if( !is_block ) Stg_SETERRQ(PETSC_ERR_ARG_WRONGSTATE, "MatCreateSNESBlockMFFD() only valid for block objects");
	
	PetscObjectGetComm( (PetscObject)snes, &comm );
	
	VecGetSize( snes->vec_func,&Mr );
	PetscMalloc( sizeof(PetscInt)*Mr, &M );
	
	VecBlockGetSubVectors( snes->vec_func, &sub_vecs );
	for( i=0; i<Mr; i++ ){
		VecGetSize( sub_vecs[i], &M[i] );
	}
	VecBlockRestoreSubVectors( snes->vec_func );
	
	
	MatCreate( comm, J );
	MatSetSizes( *J, Mr,Mr, Mr,Mr );
	MatSetType( *J, "block" );
	
	for( i=0; i<Mr; i++ ) {
		for( j=0; j<Mr; j++ ) {
			MatCreate( comm, &subMat );
			MatSetSizes( subMat, PETSC_DECIDE,PETSC_DECIDE, M[i],M[j] );
			MatSetType( subMat, MATAIJ );
			
			MatAssemblyBegin(subMat,MAT_FINAL_ASSEMBLY);
			MatAssemblyEnd(subMat,MAT_FINAL_ASSEMBLY);
			
			MatBlockSetValue( *J, i,j, subMat, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
			
			Stg_MatDestroy( &subMat );
		}
	}
	MatAssemblyBegin(*J,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*J,MAT_FINAL_ASSEMBLY);
	
	ierr=PetscFree( M );CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscErrorCode BlockMatMFFDComputeJacobian(SNES snes,Vec x,Mat *jac,Mat *B,MatStructure *flag,void *dummy)
{
	MatAssemblyBegin(*jac,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*jac,MAT_FINAL_ASSEMBLY);
	PetscFunctionReturn(0);
}

/* useage */
/*
SNESCreate(PETSC_COMM_WORLD,&snes);
SNESSetFunction(snes,r,FormFunction,NULL);
MatCreateSNESBlockMF( snes, &J );
MatBlockMFFDSetFromOptions(J);
SNESSetJacobian(snes,J,J,MatBlockMFFDComputeJacobian, FormFunction );
Stg_MatDestroy(J);


*/
