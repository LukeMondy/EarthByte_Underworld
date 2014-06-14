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
    This file implements preconditioned CGR (Conjugate Gradient with Reconjugasion) method.  
    References:
        Farhat C., Crivelli L. & Roux F. X., "Extending substructure based iterative solvers to multiple load and repeated analyses",
        Comput. Methods Appl. Mech. Engrg., 117, (1994), pp 195-209.
        
        Saad Y. "On the Lanczos method for solving symmetric linear systems with several right-hand sides",
        Mathematics of Computation, 48, (1987), pp 651-662.
        
        Rey C. & Lene F., "Reuse of krylov spaces in the solution of large-scale nonlinear elasticity problems",
        
*/



#include <stdlib.h>

#include "petsc.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscvec.h"

#include "private/kspimpl.h"
#include "private/ksp/cgr.h"


typedef struct _KSP_CGR* KSP_CGR;

struct _KSP_CGR {
	PetscInt   restart_its; /* max number of vectors we store */
	PetscTruth preallocate_vectors; /* preallocate all vectors at once (only is restart != -1) */
	PetscTruth store_all_vectors; 
	PetscInt   n_search_directions; /* number of search directions currently in use */
	PetscInt   n_allocated_vectors; /* number of vectors currently allocated */
	PetscInt   n_cnt;
	Vec         *S;
	PetscScalar *ST_iAS_i;
};


PetscErrorCode KSPSetFromOptions_CGR(KSP ksp)
{
	PetscErrorCode ierr;
	PetscInt       restart;
	KSP_CGR        cgr = (KSP_CGR)ksp->data;
	PetscTruth     flg;
	
	PetscFunctionBegin;
	ierr = PetscOptionsHead("KSP CGR Options");CHKERRQ(ierr);
	
	ierr = PetscOptionsInt("-ksp_cgr_restart","Number of search directions stored","KSPCGRSetRestart",cgr->restart_its,&restart,&flg);CHKERRQ(ierr);
	if (flg) { ierr = KSPCGRSetRestart(ksp,restart);CHKERRQ(ierr); }
	
	ierr = PetscOptionsName("-ksp_cgr_preallocate","Preallocate search directions","KSPCGRSetPreAllocateVectors",&flg);CHKERRQ(ierr);
	if (flg) {ierr = KSPCGRSetPreAllocateVectors(ksp,PETSC_TRUE);CHKERRQ(ierr);}
	
	ierr = PetscOptionsName("-ksp_cgr_store_all","Store ALL search directions","KSPCGRSetStoreAllDirections",&flg);CHKERRQ(ierr);
	if (flg) {ierr = KSPCGRSetStoreAllDirections(ksp,PETSC_TRUE);CHKERRQ(ierr);}
	
	/*
	ierr = PetscOptionsName("-ksp_gmres_krylov_monitor","Plot the Krylov directions","KSPSetMonitor",&flg);CHKERRQ(ierr);
	if (flg) {
		PetscViewers viewers;
		ierr = PetscViewersCreate(ksp->comm,&viewers);CHKERRQ(ierr);
		ierr = KSPSetMonitor(ksp,KSPGMRESKrylovMonitor,viewers,(PetscErrorCode (*)(void*))PetscViewersDestroy);CHKERRQ(ierr);
	}
	*/
	ierr = PetscOptionsTail();CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "KSPSetUp_CGR"
PetscErrorCode KSPSetup_CGR( KSP ksp )
{
	KSP_CGR cgr = (KSP_CGR)ksp->data;
	PetscInt k, min;
	Vec x;
	
	if( cgr->restart_its < 0 ) {
		Stg_SETERRQ( PETSC_ERR_USER, "KSPCGR: Need to set restart_its with KSPCGRSetRestart() \n" );
	}
	
	/* If we have changed the operator via KSPSetOperators(), we need to reset the history */	
	KSPCGRResetSearchDirection( ksp );
	
	/* 
	If we have already allocated vectors, leave function. 
	This situation will occur if KSPSetOperators() is called repeatedly.
	*/
	if(  cgr->S != PETSC_NULL && cgr->ST_iAS_i != PETSC_NULL ) {
		PetscFunctionReturn(0);
	}
	
	if( cgr->preallocate_vectors == PETSC_TRUE && cgr->store_all_vectors == PETSC_FALSE ) {
		cgr->n_allocated_vectors = cgr->restart_its;
	}
	else {
		min = PetscMin( (PetscInt)cgr->restart_its, (PetscInt)5 );
		cgr->n_allocated_vectors = min;
	}
	
	/* allocate space for the saved directions */
	cgr->ST_iAS_i = (PetscScalar*)malloc( sizeof(PetscScalar) * cgr->n_allocated_vectors );
	
	
	cgr->S = (Vec*)malloc( sizeof(Vec) * cgr->n_allocated_vectors );
	
	/* get the vector so we duplicate it */
	KSPGetSolution( ksp, &x );
	
	for( k=0; k<cgr->n_allocated_vectors; k++ ) {
		VecDuplicate( x, &cgr->S[k] );
	}
	
	PetscFunctionReturn(0);
}


PetscErrorCode KSPDestroy_CGR( KSP ksp )
{
	KSP_CGR        cgr = (KSP_CGR)ksp->data; 
	PetscInt k;
	
	
	/* Public functions declared here */
	PetscObjectComposeFunctionDynamic((PetscObject)ksp,"KSPCGRSetRestart_C", "", PETSC_NULL);
	PetscObjectComposeFunctionDynamic((PetscObject)ksp,"KSPCGRSetPreallocateVectors_C", "", PETSC_NULL);
	PetscObjectComposeFunctionDynamic((PetscObject)ksp,"KSPCGRSetStoreAllDirections_C", "", PETSC_NULL);
	PetscObjectComposeFunctionDynamic((PetscObject)ksp,"KSPCGRResetSearchDirection_C", "", PETSC_NULL);
	
	for( k=0; k<cgr->n_allocated_vectors; k++ ) {
		Stg_VecDestroy( & cgr->S[k] );
	}
	free( cgr->S );
	free( cgr->ST_iAS_i );
	
	KSPDefaultDestroy(ksp);
//	free( cgr );
	
	PetscFunctionReturn(0);
}

PetscErrorCode KSPView_CGR(KSP ksp,PetscViewer viewer)
{
	KSP_CGR        cgr = (KSP_CGR)ksp->data; 
	PetscTruth isascii;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&isascii);
	if (isascii) {
	
		PetscViewerASCIIPushTab(viewer); //1
		
		if( cgr->store_all_vectors == PETSC_TRUE ) {      
			ierr = PetscViewerASCIIPrintf(viewer,"CGR: store_all_vectors=TRUE \n" );CHKERRQ(ierr);
		}
		else {
			ierr = PetscViewerASCIIPrintf(viewer,"CGR: store_all_vectors=FALSE \n" );CHKERRQ(ierr); 
			ierr = PetscViewerASCIIPrintf(viewer,"CGR: restart=%D \n",cgr->restart_its );CHKERRQ(ierr);
		}
		
		if( cgr->preallocate_vectors == PETSC_TRUE ) {
			ierr = PetscViewerASCIIPrintf(viewer,"CGR: preallocate_vectors=TRUE \n" );CHKERRQ(ierr);
		}
		else { ierr = PetscViewerASCIIPrintf(viewer,"CGR: preallocate_vectors=FALSE \n" );CHKERRQ(ierr); }
	
		ierr = PetscViewerASCIIPrintf(viewer,"CGR: vectors_allocated=%d \n", cgr->n_allocated_vectors );CHKERRQ(ierr);
		PetscViewerASCIIPopTab(viewer); //1
	}	
	
	PetscFunctionReturn(0);
}












PetscErrorCode KSPCGRUpdateSearchDirectionHistory_CGR( KSP_CGR ctx, Vec p, PetscScalar p_dot_q )
{
	PetscInt i;
	Vec *tmp_v;
	PetscScalar *tmp_s;
	PetscInt new_size;
	
	
	/* If we have stored the max allowable number of vectors, reset the search direction counter */
	if( ctx->n_search_directions == ctx->restart_its && ctx->store_all_vectors == PETSC_FALSE) {
		ctx->n_search_directions = 0;
	}
	
	/* allocate more vectors if required */
	if( ctx->n_search_directions == ctx->n_allocated_vectors ) {
		/* increase by 10% */
		new_size = ctx->n_search_directions + 10;
		/* prevent allocated size from becoming larger than size of the restart number */
		if( new_size > ctx->restart_its && ctx->store_all_vectors == PETSC_FALSE ) {
			new_size = ctx->restart_its;
		}
		
		/* realloc lists */
		tmp_v = (Vec*)realloc( ctx->S, sizeof(Vec)*new_size );
		ctx->S = tmp_v;
		for( i=ctx->n_allocated_vectors; i<new_size; i++ ) {
			VecDuplicate( p, &ctx->S[i] );
		}
		
		tmp_s = (PetscScalar*)realloc( ctx->ST_iAS_i, sizeof(PetscScalar) * new_size );
		ctx->ST_iAS_i = tmp_s;
		
		ctx->n_allocated_vectors = new_size;
	}
	
	/* Store new direction */
	ctx->ST_iAS_i[ ctx->n_search_directions ] = p_dot_q;
	VecCopy( p, ctx->S[ ctx->n_search_directions ] );
	
	ctx->n_search_directions++;
	
	
	ctx->n_cnt++;
	
	PetscFunctionReturn(0);
}




PetscErrorCode CGR_iteration( KSP ksp )
{
	KSP_CGR        ctx = (KSP_CGR)ksp->data; 
	PetscErrorCode ierr;
	PetscInt i,k;
	PetscScalar rho_m1, rho_m2, p_dot_q, alpha, beta;
	PetscScalar alpha_q, skAs;
	Mat A, Pmat;
	Vec b, x;
	Vec sk, Asq;
	Vec r;
	Vec z,p,q;
	PetscInt n_old_vectors;
	MatStructure flag;
	PetscReal norm, dp = 0.0;
	
	n_old_vectors = ctx->n_search_directions;
	KSPGetOperators( ksp, &A, &Pmat, &flag );
	KSPGetSolution( ksp, &x );
	KSPGetRhs( ksp, &b );
	
	
	
	VecDuplicate( b, &r );
	VecDuplicate( b, &z );
	VecDuplicate( b, &p );
	VecDuplicate( b, &q );
	VecDuplicate( b, &sk );
	VecDuplicate( b, &Asq );
	
	
	
	MatMult( A, x, r );
	VecAYPX( r, -1.0, b ); /* r <- -r + b */
	
//	VecNorm( r, NORM_2, &norm );
//	printf("iteration[%d]- norm = %14.12e \n", 0, norm );
	if (ksp->normtype == KSP_NORM_PRECONDITIONED) {
		ierr = KSP_PCApply(ksp,r,z);CHKERRQ(ierr); /*     z <- Br         */
		ierr = VecNorm(z,NORM_2,&dp);CHKERRQ(ierr);                /*    dp <- z'*z = e'*A'*B'*B*A'*e'   */
	} else if (ksp->normtype == KSP_NORM_UNPRECONDITIONED) {
		ierr = VecNorm(r,NORM_2,&dp);CHKERRQ(ierr);                /*    dp <- r'*r = e'*A'*A*e         */
	} else dp = 0.0;
	KSPLogResidualHistory(ksp,dp);
	KSPMonitor(ksp,0,dp);                              /* call any registered monitor routines */
	ksp->rnorm = dp;
	
	ierr = (*ksp->converged)(ksp,0,dp,&ksp->reason,ksp->cnvP);CHKERRQ(ierr);      /* test for convergence */
	if (ksp->reason) PetscFunctionReturn(0);
	
	
	
	
	
	
	
	i = 1;
	do {
		n_old_vectors = ctx->n_search_directions;
		
		
		//VecCopy( r, z );
		KSP_PCApply(ksp,r,z); /*     z <- Br         */
		
		VecDot( r, z, &rho_m1 );
		
		if( i==1 ) {
			VecCopy( z, p );
		}
		else {
			beta = rho_m1/rho_m2;
			VecAYPX( p, beta, z );
		}
		
		/* Modify the search direction p */
		VecCopy( p, sk ); /* sk <- p */
	/*	
		for( k=0; k<n_old_vectors; k++ ) {
			MatMult( A, ctx->S[k], Asq );
			
			VecDot( sk, Asq, &skAs );
			alpha_q = skAs / ctx->ST_iAS_i[k];
			VecAXPY( p, -alpha_q, ctx->S[k] ); // p <- p -alpha_q s_q //
		}
	*/
		MatMult( A, sk, Asq );
		for( k=0; k<n_old_vectors; k++ ) {
			/* Use first expression below eq 21 in Farhat */
			VecDot( ctx->S[k], Asq, &skAs );
			alpha_q = skAs / ctx->ST_iAS_i[k];
			VecAXPY( p, -alpha_q, ctx->S[k] ); /* p <- p -alpha_q s_q */
		}

		/*********************************/
		
		MatMult( A, p, q );
		
		VecDot( p, q, &p_dot_q );
		alpha = rho_m1 / p_dot_q;
		
		
		VecAXPY( x, alpha, p ); /* x <- x + alpha p */
		VecAXPY( r, -alpha, q ); /* r <- r - alpha q */
		
		VecNorm( r, NORM_2, &norm );
		rho_m2 = rho_m1;
		
		KSPCGRUpdateSearchDirectionHistory_CGR( ctx, p, p_dot_q );
		
		
		
		if (ksp->normtype == KSP_NORM_PRECONDITIONED) {
			ierr = KSP_PCApply(ksp,r,z);CHKERRQ(ierr);        /*     z <- Br         */
			ierr = VecNorm(z,NORM_2,&dp);CHKERRQ(ierr);              /*    dp <- z'*z       */
		} else if (ksp->normtype == KSP_NORM_UNPRECONDITIONED) {
			ierr = VecNorm(r,NORM_2,&dp);CHKERRQ(ierr);              /*    dp <- r'*r       */
		} else {
			dp = 0.0;
		}
		ksp->rnorm = dp;
		KSPLogResidualHistory(ksp,dp);
		KSPMonitor(ksp,i,dp);
		ierr = (*ksp->converged)(ksp,i,dp,&ksp->reason,ksp->cnvP);CHKERRQ(ierr);
		if (ksp->reason) break;
		
		
		i++;
		ksp->its = i;
	} while (i<ksp->max_it);
	if (i >= ksp->max_it) {
		ksp->reason = KSP_DIVERGED_ITS;
	}
	
	
	Stg_VecDestroy( & r );
	Stg_VecDestroy( & q );
	Stg_VecDestroy( & p );
	Stg_VecDestroy( & z );
	Stg_VecDestroy( & sk );
	Stg_VecDestroy( & Asq );
	
	PetscFunctionReturn(0);

}




PetscErrorCode KSPSolve_CGR( KSP ksp )
{
	KSP_CGR     ctx = (KSP_CGR)ksp->data; 
	PetscScalar *STf, *yInit;
	PetscInt k;
	Vec x, b;
	
	
	KSPGetSolution( ksp, &x );
	KSPGetRhs( ksp, &b );
	ksp->its = 0;
	
	if( ctx->n_search_directions == 0 ) {
		VecSet( x, 0.0 );
	}
	else {
		/* compute initial guess using previous directions from the i^th right hand side */
		/* form S^T_i-1 f_i */
		STf = (PetscScalar*)malloc( sizeof(PetscScalar) * ctx->n_search_directions );
		yInit = (PetscScalar*)malloc( sizeof(PetscScalar) * ctx->n_search_directions );
		
		
//		VecMDot( ctx->n_search_directions, b, ctx->S, STf );
		VecMDot( b, ctx->n_search_directions, ctx->S, STf );
		
		/* solve S^T_i-1 A S_i-1 y^0_i = S^T_i-1 f_i */
		for( k=0; k<ctx->n_search_directions; k++ ) {
			yInit[k] = STf[k]/ctx->ST_iAS_i[k];
		}
		
		/* x^0_i = S_i-1 y^0_i */
		VecSet( x, 0.0 );
		VecMAXPY( x, ctx->n_search_directions, yInit, ctx->S ); /* X <- X + y[] * S[] */
		
		free( STf );
		free( yInit );
	}
	CGR_iteration( ksp );
	
	
	PetscFunctionReturn(0);
}




/* ----------------------- Public function for KSP_CGR ----------------------- */




EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "KSPCGRSetRestart_CGR"
PetscErrorCode KSPCGRSetRestart_CGR( KSP ksp, PetscInt i )
{
	KSP_CGR cgr = (KSP_CGR)ksp->data;
	
	cgr->restart_its = i;
	
	PetscFunctionReturn(0);
}
EXTERN_C_END


#undef __FUNCT__  
#define __FUNCT__ "KSPCGRSetRestart"
/*MC
   KSPCGRSetRestart - Sets the number of iterations at which CGR restarts.

   Collective on KSP

   Input Parameter:
+  ksp - the Krylov space context
-  restart - integer restart value

   Options Database:
.   -ksp_cgr_restart <positive integer>

   Level: intermediate

   Concepts: reconjugasion

.keywords: KSP, CGR, restart

.seealso: 
M*/
PetscErrorCode KSPCGRSetRestart( KSP ksp, PetscInt restart )
{
	PetscErrorCode ierr,(*f)(KSP,PetscInt);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)ksp,"KSPCGRSetRestart_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(ksp,restart);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}


EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "KSPCGRSetPreAllocateVectors_CGR"
PetscErrorCode KSPCGRSetPreAllocateVectors_CGR( KSP ksp, PetscTruth flg )
{	
	KSP_CGR cgr = (KSP_CGR)ksp->data;
	
	cgr->preallocate_vectors = flg;
	
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "KSPCGRSetPreAllocateVectors"
/*MC
   KSPCGRSetPreAllocateVectors - Indicates whether we wish to preallocate space for all the search directions.

   Collective on KSP

   Input Parameter:
+  ksp - the Krylov space context
-  flg - indicate whether we wish to preallocate stored vectors.

   Options Database:
.   -ksp_cgr_preallocate

   Level: beginner

   Concepts: reconjugasion

.seealso: 
M*/
PetscErrorCode KSPCGRSetPreAllocateVectors( KSP ksp, PetscTruth flg )
{	
	PetscErrorCode ierr,(*f)(KSP,PetscTruth);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)ksp,"KSPCGRSetPreAllocateVectors_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(ksp,flg);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}



EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "KSPCGRSetStoreAllDirections_CGR"
PetscErrorCode KSPCGRSetStoreAllDirections_CGR( KSP ksp, PetscTruth flg )
{
	KSP_CGR cgr = (KSP_CGR)ksp->data;
	
	cgr->store_all_vectors = flg;
	
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "KSPCGRSetStoreAllDirections"
/*MC
   KSPCGRSetStoreAllDirections - Indicates whether we wish to store all the search directions.

   Collective on KSP

   Input Parameter:
+  ksp - the Krylov space context
-  flg - indicate whether we wish to store the directions.

   Options Database:
.   -ksp_cgr_store_all

   Level: beginner

   Concepts: reconjugasion

.seealso: 
M*/
PetscErrorCode KSPCGRSetStoreAllDirections( KSP ksp, PetscTruth flg )
{
	PetscErrorCode ierr,(*f)(KSP,PetscTruth);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)ksp,"KSPCGRSetStoreAllDirections_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(ksp,flg);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}




EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "KSPCGRResetSearchDirection_CGR"
PetscErrorCode KSPCGRResetSearchDirection_CGR( KSP ksp )
{
	KSP_CGR cgr = (KSP_CGR)ksp->data;
	
	cgr->n_search_directions = 0;
	
	PetscFunctionReturn(0);
}
EXTERN_C_END


#undef __FUNCT__  
#define __FUNCT__ "KSPCGRResetSearchDirection"
/*MC
    KSPCGRResetSearchDirection - Resets the counter used to monitor the 
    number of stored directions.    

    Collective on KSP

    Input Parameters:
.   ksp - the Krylov space context

    Level: advanced
    
    Note:
    For basic usage of the KSPCGR the user need not explicitly call
    KSPCGRResetSearchDirection(), as this actions will happen automatically
    within KSPSolve().

.keywords: CGR, conjugate gradient reconjugasion, serch direction 
M*/
PetscErrorCode KSPCGRResetSearchDirection( KSP ksp )
{
	PetscErrorCode ierr,(*f)(KSP);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)ksp,"KSPCGRResetSearchDirection_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(ksp);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}




/*
KSPCreate_CGR - Creates the data structure for the Generalized Krylov Correction 
method CG and sets the function pointers for all the routines it needs to call (KSPSolve_CGR() etc)

It must be wrapped in EXTERN_C_BEGIN to be dynamically linkable in C++
*/
/*MC
KSPCGR - The preconditioned conjugate gradient with reconjugasion (PCGR) iterative method

   Options Database Keys:
+   -ksp_cgr_restart <restart> - the number of Krylov directions to orthogonalize against
.   -ksp_cgr_preallocate - preallocate all the Krylov search directions initially (otherwise groups of vectors are allocated as needed)
-   -ksp_cgr_store_all - stores all search directions

Level: beginner

Notes: The PCGR method requires both the matrix and preconditioner to 
be symmetric positive (semi) definite

.seealso:  KSPCreate(), KSPSetType(), KSPType (for list of available types), KSP,
KSPCGSetType()

M*/
EXTERN_C_BEGIN
PetscErrorCode KSPCreate_CGR( KSP ksp )
{
	KSP_CGR _cgr;
	PetscErrorCode ierr;
	
	/* Create ctx for generalized cg solver */
//	_cgr = (KSP_CGR)malloc( sizeof(struct _KSP_CGR) );
	
	ierr = PetscMalloc( sizeof(struct _KSP_CGR), &_cgr );CHKERRQ(ierr);
	ierr = PetscMemzero( _cgr, sizeof(struct _KSP_CGR) );CHKERRQ(ierr);
	PetscLogObjectMemory( ksp, sizeof(struct _KSP_CGR) );
	
//	PetscNew( struct _KSP_CGR,&_cgr );
//	PetscLogObjectMemory( ksp, sizeof(struct _KSP_CGR) );
	
	ksp->data                              = (void*)_cgr;
	ksp->pc_side                           = PC_LEFT;
	ksp->ops->buildsolution                = PETSC_NULL;
	ksp->ops->setup                        = KSPSetup_CGR;
	ksp->ops->solve                        = KSPSolve_CGR;
	ksp->ops->destroy                      = KSPDestroy_CGR;
	ksp->ops->view                         = KSPView_CGR;
	ksp->ops->setfromoptions               = KSPSetFromOptions_CGR;
	ksp->ops->computeextremesingularvalues = PETSC_NULL;
	ksp->ops->computeeigenvalues           = PETSC_NULL;
	
	
	/* Public functions declared here */
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)ksp,"KSPCGRSetRestart_C",
			"KSPCGRSetRestart_CGR",
			KSPCGRSetRestart_CGR);CHKERRQ(ierr);
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)ksp,"KSPCGRSetPreallocateVectors_C",
			"KSPCGRSetPreAllocateVectors_CGR",
			KSPCGRSetPreAllocateVectors_CGR);CHKERRQ(ierr);
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)ksp,"KSPCGRSetStoreAllDirections_C",
			"KSPCGRSetStoreAllDirections_CGR",
			KSPCGRSetStoreAllDirections_CGR);CHKERRQ(ierr);
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)ksp,"KSPCGRResetSearchDirection_C",
			"KSPCGRResetSearchDirection_CGR",
			KSPCGRResetSearchDirection_CGR);CHKERRQ(ierr);
	
	
	/* init variables */
	_cgr->restart_its           = 30;
	_cgr->preallocate_vectors   = PETSC_FALSE;
	_cgr->store_all_vectors     = PETSC_FALSE;
	_cgr->n_search_directions   = 0;
	_cgr->n_allocated_vectors   = 0;
	_cgr->n_cnt                 = 0;
	_cgr->S                     = PETSC_NULL;
	_cgr->ST_iAS_i              = PETSC_NULL;
	
	PetscFunctionReturn(0);
}
EXTERN_C_END

