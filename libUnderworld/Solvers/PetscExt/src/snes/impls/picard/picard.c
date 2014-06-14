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
#define PETSCSNES_DLL

#include "picard-impl.h" /* contains context definition */
#include "private/snes/snespicard.h"  /* contains prototypes for public functions */


PetscErrorCode _SNESPicard_FormJacobianNULL(SNES snes,Vec a,Mat *A,Mat *B,MatStructure *s,void *c) 
{
	PetscFunctionReturn(0);
}

/*
F(x) := Ax - b

Notes - This will get loaded and accessed via SNESComputeFunction().
Users should not call this function via SNESComputeFunction() unless the vector F
is the same as that specified via SNESSetFunction. The result F will not be consistent
if this is not followd.

If you really insist on calling this function, do the following.

		User_Context *ctx;
		Vec R;
		PetscErrorCode (*func)(SNES,Vec,Vec,void*);
		
		SNESGetFunction( snes, &R, &func, (void**)&ctx );
		func( snes,x,R,(void*)ctx );

[  F must be the same vector specified with SNESSetFunction()  ]
*/
PetscErrorCode _SNESPicard_ComputeFunction( SNES snes, Vec X, Vec F, void *_ctx )
{
	SNES_PICARD      *ctx = (SNES_PICARD*)_ctx;
	PetscErrorCode   ierr;
	MatStructure     flg = DIFFERENT_NONZERO_PATTERN;
	
	
	ierr = SNESPicardComputeOperator(snes,X,&ctx->Amat,&ctx->Pmat,&flg);CHKERRQ(ierr);
	ierr = SNESPicardComputeRhs(snes,X,ctx->rhs);CHKERRQ(ierr);
	
	MatMult( ctx->Amat, X, F );
	VecAXPY( F, -1.0, ctx->rhs );
	
	PetscFunctionReturn(0);
}




/* ----------------------- Public function for SNES_Picard ----------------------- */



EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "SNESPicardSetOperator_PICARD"
PetscErrorCode PETSCMAT_DLLEXPORT SNESPicardSetOperator_PICARD(
		SNES snes,Mat A,Mat B,
		PetscErrorCode (*compute_operator)(SNES,Vec,Mat*,Mat*,MatStructure*,void*),void *op_ctx )
{
	PetscErrorCode  ierr;
	SNES_PICARD      *ctx = (SNES_PICARD*)snes->data;
	
	PetscFunctionBegin;
	
	if (compute_operator) ctx->compute_operator = compute_operator;
	if (op_ctx) ctx->op_ctx = op_ctx;
	if (A) {
		PetscObjectReference((PetscObject)A);
		if (ctx->Amat) {Stg_MatDestroy( &ctx->Amat);}
		ctx->Amat = A;
	}
	if (B) {
		PetscObjectReference((PetscObject)B);
		if (ctx->Pmat) {Stg_MatDestroy( &ctx->Pmat);}
		ctx->Pmat = B;
	}
	
	/*
	We have to do this as
	1) This function must be called before SNESolve. Hence we can ensure that snes->compute function is
	called before SNESSetUp() is called.
	*/
	if( !ctx->F ) {
		ierr = MatGetVecs( A, &ctx->F, PETSC_NULL );CHKERRQ(ierr);
		/* define the ComputeFunction */
		SNESSetFunction( snes, ctx->F, _SNESPicard_ComputeFunction, (void*)ctx );
	}
	if( !snes->jacobian ) {
		ierr = SNESSetJacobian(snes,A,B,_SNESPicard_FormJacobianNULL,PETSC_NULL);CHKERRQ(ierr);
	}
	
	
	
	
	PetscFunctionReturn(0);
}
EXTERN_C_END
#undef __FUNCT__  
#define __FUNCT__ "SNESPicardSetOperator"
PetscErrorCode PETSCMAT_DLLEXPORT SNESPicardSetOperator( SNES snes,Mat A,Mat B,
		PetscErrorCode (*compute_operator)(SNES,Vec,Mat*,Mat*,MatStructure*,void*),void *op_ctx )
{
	PetscErrorCode ierr,(*f)(SNES,Mat,Mat, PetscErrorCode (*)(SNES,Vec,Mat*,Mat*,MatStructure*,void*), void*);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)snes,"SNESPicardSetOperator_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(snes,A,B,compute_operator,op_ctx);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}


EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "SNESPicardSetRhs_PICARD"
PetscErrorCode PETSCMAT_DLLEXPORT SNESPicardSetRhs_PICARD(
		SNES snes, Vec rhs,
		PetscErrorCode (*compute_rhs)(SNES,Vec,Vec,void*),void *rhs_ctx )
{
	PetscErrorCode  ierr;
	SNES_PICARD      *ctx = (SNES_PICARD*)snes->data;
	
	PetscFunctionBegin;
	
	if (compute_rhs) ctx->compute_rhs = compute_rhs;
	if (rhs_ctx) ctx->rhs_ctx = rhs_ctx;
	if (rhs) {
		PetscObjectReference((PetscObject)rhs);
		if (ctx->rhs) {Stg_VecDestroy( &ctx->rhs);}
		ctx->rhs = rhs;
	}
	
	/*
	We have to do this as
	1) This function must be called before SNESolve. Hence we can ensure that snes->compute function is
	called before SNESSetUp() is called.
	2) Note that if the rhs is NULL, we'll catch this when we set the operator
	*/
	if( !ctx->F && !rhs) {
		ierr = VecDuplicate(rhs, &ctx->F );CHKERRQ(ierr);
		/* define the ComputeFunction */
		SNESSetFunction( snes, ctx->F, _SNESPicard_ComputeFunction, (void*)ctx );
	}
	
	
	PetscFunctionReturn(0);
}
EXTERN_C_END
#undef __FUNCT__  
#define __FUNCT__ "SNESPicardSetRhs"
PetscErrorCode PETSCMAT_DLLEXPORT SNESPicardSetRhs( SNES snes,Vec rhs,
		PetscErrorCode (*compute_rhs)(SNES,Vec,Vec,void*),void *rhs_ctx )
{
	PetscErrorCode ierr,(*f)(SNES,Vec, PetscErrorCode (*)(SNES,Vec,Vec,void*), void*);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)snes,"SNESPicardSetRhs_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(snes,rhs,compute_rhs,rhs_ctx);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "SNESPicardComputeOperator"
PetscErrorCode SNESPicardComputeOperator(SNES snes,Vec X,Mat *A,Mat *B,MatStructure *flg)
{
	SNES_PICARD     *ctx = (SNES_PICARD*)snes->data;
	PetscTruth      picard;
	
	PetscFunctionBegin;
	
	picard = PETSC_TRUE;
	PetscTypeCompare( (PetscObject)snes, "picardext", &picard );
	if( !picard ) {
		Stg_SETERRQ(PETSC_ERR_ARG_WRONGSTATE, "SNESPicardComputeOperator() only valid if SNESType = picardext \n" );
	}
	
	if (!ctx->compute_operator) PetscFunctionReturn(0);
	
	ctx->compute_operator( snes,X,A,B,flg,ctx->op_ctx);
	
	PetscFunctionReturn(0);
}


/*
x   - contains the current solution
rhs - where will store the result
*/
#undef __FUNCT__  
#define __FUNCT__ "SNESPicardComputeRhs"
PetscErrorCode SNESPicardComputeRhs(SNES snes,Vec x,Vec rhs)
{
	SNES_PICARD     *ctx = (SNES_PICARD*)snes->data;
	PetscTruth      picard;
	
	PetscFunctionBegin;
	
	picard = PETSC_TRUE;
	PetscTypeCompare( (PetscObject)snes, "picardext", &picard );
	if( !picard ) {
		Stg_SETERRQ(PETSC_ERR_ARG_WRONGSTATE, "SNESPicardComputeRhs() only valid if SNESType = picardext \n" );
	}
	
	if( !ctx->compute_rhs ) PetscFunctionReturn(0);
	
	ctx->compute_rhs( snes,x,rhs,ctx->rhs_ctx);
	
	PetscFunctionReturn(0);
}

/* --------------------- Internal routines for SNES Package --------------------- */

#undef __FUNCT__  
#define __FUNCT__ "SNESSolve_PicardExt"
PetscErrorCode SNESSolve_PicardExt(SNES snes)
{ 
	SNES_PICARD     *ctx = (SNES_PICARD*)snes->data;
	PetscErrorCode  ierr;
	PetscInt maxits,i,lits;
	PetscReal fnorm,xnorm;
	Vec X, F;
	MatStructure flg = DIFFERENT_NONZERO_PATTERN;
	KSPConvergedReason kspreason;
	
	
	PetscFunctionBegin;
	
	snes->numFailures            = 0;
	snes->numLinearSolveFailures = 0;
	snes->reason                 = SNES_CONVERGED_ITERATING;
	snes->iter = 0;
	
	maxits  = snes->max_its;    /* maximum number of iterations */
	X       = snes->vec_sol;    /* solution vector */
	F       = ctx->F;           /* residual vector */
	
	
	/* Build operator [A(u)], and rhs {b(u)} */
	ierr = SNESPicardComputeOperator(snes,X,&ctx->Amat,&ctx->Pmat,&flg);
	ierr = SNESPicardComputeRhs(snes,X,ctx->rhs);
	
	
	/* Compute nonlinear residual */
	ierr = PetscObjectGrantAccess(snes);CHKERRQ(ierr);
	ierr = SNESComputeFunction(snes,X,F);CHKERRQ(ierr);
	ierr = VecNorm(F,NORM_2,&fnorm);CHKERRQ(ierr);	/* fnorm <- ||F||  */
	if (fnorm != fnorm) Stg_SETERRQ(PETSC_ERR_FP,"User provided compute function generated a Not-a-Number");
	ierr = PetscObjectTakeAccess(snes);CHKERRQ(ierr);
	snes->norm = fnorm;
	ierr = PetscObjectGrantAccess(snes);CHKERRQ(ierr);
	SNESLogConvHistory(snes,fnorm,0);
	SNESMonitor(snes,0,fnorm);
	
	i = 0;
	do {
		/* Copy current solution into F */
		VecCopy( X, F ); /* F <- X */
		
		
		/* Solve [A]{u} = {b} */
		ierr = KSPSetOperators(snes->ksp,ctx->Amat,ctx->Pmat,flg);CHKERRQ(ierr);
		ierr = SNES_KSPSolve(snes,snes->ksp,ctx->rhs,X);CHKERRQ(ierr);
		ierr = KSPGetConvergedReason(snes->ksp,&kspreason);CHKERRQ(ierr);
		if (kspreason < 0) {
			if (++snes->numLinearSolveFailures >= snes->maxLinearSolveFailures) {
				snes->reason = SNES_DIVERGED_LINEAR_SOLVE;
				PetscFunctionReturn(0);
			}
		}
		ierr = KSPGetIterationNumber(snes->ksp,&lits);CHKERRQ(ierr);
		snes->linear_its += lits;
		ierr = PetscInfo2(snes,"iter=%D, linear solve iterations=%D\n",snes->iter,lits);CHKERRQ(ierr);
		
		
		/* compute residuals */
		
		
		/* 1) Compute u_i - u_{i-1} */
		VecAXPY( F, -1.0, X ); /* F <- F - X */
		ierr = VecNorm(F,NORM_2,&xnorm);CHKERRQ(ierr);	/* xnorm = || X_i - X_{i-1} || */
		
		/* 2) Compute F */
		/* ReBuild operator [A(u)], and rhs {b(u)} */
		//SNESPicardComputeOperator(snes,X,&ctx->Amat,&ctx->Pmat,&flg);
		//SNESPicardComputeRhs(snes,X,ctx->rhs);
		
		ierr = SNESComputeFunction(snes,X,F);CHKERRQ(ierr);
		ierr = VecNorm(F,NORM_2,&fnorm);CHKERRQ(ierr);	/* fnorm = || F || */
		
		
		ierr = PetscObjectTakeAccess(snes);CHKERRQ(ierr);
		snes->iter = i+1;
		snes->norm = fnorm;
		ierr = PetscObjectGrantAccess(snes);CHKERRQ(ierr);
		SNESLogConvHistory(snes,fnorm,lits);
		SNESMonitor(snes,i+1,fnorm);
		
		
		/* check for convergence */
		if (snes->ops->converged) {
			ierr = (*snes->ops->converged)(snes,snes->iter,xnorm,1.0e6,fnorm,&snes->reason,snes->cnvP);CHKERRQ(ierr);
			if (snes->reason) break;
		}
		
		
		/* ReBuild operator [A(u)], and rhs {b(u)} */
		//SNESPicardComputeOperator(snes,X,&ctx->Amat,&ctx->Pmat,&flg);
		//SNESPicardComputeRhs(snes,X,ctx->rhs);
		
		i++;
		
	} while( i < snes->max_its );
	
	if (X != snes->vec_sol) {
		ierr = VecCopy(X,snes->vec_sol);CHKERRQ(ierr);
	}
	if (F != snes->vec_func) {
		ierr = VecCopy(F,snes->vec_func);CHKERRQ(ierr);
	}
//	snes->vec_sol_always  = snes->vec_sol;
//	snes->vec_func_always = snes->vec_func;
	
	
	if( i >= snes->max_its ) {
		ierr = PetscInfo1(snes,"Maximum number of iterations has been reached: %D\n",maxits);CHKERRQ(ierr);
		snes->reason = SNES_DIVERGED_MAX_IT;
	}
	
	
	PetscFunctionReturn(0);
}

/*
Inexact Newton constructs iterates,
  x_k+1 = x_k + d_k, where d_k is computed to satisfy
  || J(x_k) dx_k + F(x_k) || <= eta_k || F(x_k) ||
  
  where eta_k = tau || F(x_k) ||.
  The case eta_k = 0 corresponds to the exact Newton iteration.

For picard iterates, we only systems of the form
  A(x_k) x_k+1 = b(x_k).
When we solve such systems, we will use a stopping condition defined by eta_k above.
Not sure whether eta_k needs to be an absolute or relative tolerance.

*/
#undef __FUNCT__  
#define __FUNCT__ "SNESSolve_InexactPicard"
PetscErrorCode SNESSolve_InexactPicard(SNES snes)
{ 
	SNES_PICARD     *ctx = (SNES_PICARD*)snes->data;
	PetscErrorCode  ierr;
	PetscInt maxits,i,lits;
	PetscReal fnorm,xnorm;
	Vec X, F;
	MatStructure flg = DIFFERENT_NONZERO_PATTERN;
	KSPConvergedReason kspreason;
	PetscReal eta_k, tau;
	
	PetscFunctionBegin;
	
	snes->numFailures            = 0;
	snes->numLinearSolveFailures = 0;
	snes->reason                 = SNES_CONVERGED_ITERATING;
	snes->iter = 0;
	
	maxits  = snes->max_its;    /* maximum number of iterations */
	X       = snes->vec_sol;    /* solution vector */
	F       = ctx->F;           /* residual vector */
	
	
	tau = 1.0e-2;
	
	
	/* Build operator [A(u)], and rhs {b(u)} */
	ierr = SNESPicardComputeOperator(snes,X,&ctx->Amat,&ctx->Pmat,&flg);
	ierr = SNESPicardComputeRhs(snes,X,ctx->rhs);
	
	
	/* Compute nonlinear residual */
	ierr = PetscObjectGrantAccess(snes);CHKERRQ(ierr);
	ierr = SNESComputeFunction(snes,X,F);CHKERRQ(ierr);
	ierr = VecNorm(F,NORM_2,&fnorm);CHKERRQ(ierr);	/* fnorm <- ||F||  */
	if (fnorm != fnorm) Stg_SETERRQ(PETSC_ERR_FP,"User provided compute function generated a Not-a-Number");
	ierr = PetscObjectTakeAccess(snes);CHKERRQ(ierr);
	snes->norm = fnorm;
	ierr = PetscObjectGrantAccess(snes);CHKERRQ(ierr);
	SNESLogConvHistory(snes,fnorm,0);
	SNESMonitor(snes,0,fnorm);
	
	i = 0;
	do {
		/* Copy current solution into F */
		VecCopy( X, F ); /* F <- X */
		
		/* Compute the new tolerance */
		eta_k = tau * fnorm;
		
		/* Solve [A]{u} = {b} */
		ierr = KSPSetOperators(snes->ksp,ctx->Amat,ctx->Pmat,flg);CHKERRQ(ierr);
		
		ierr = KSPSetTolerances(snes->ksp, eta_k,eta_k, PETSC_DEFAULT,PETSC_DEFAULT );
		ierr = SNES_KSPSolve(snes,snes->ksp,ctx->rhs,X);CHKERRQ(ierr);
		
		ierr = KSPGetConvergedReason(snes->ksp,&kspreason);CHKERRQ(ierr);
		if (kspreason < 0) {
			if (++snes->numLinearSolveFailures >= snes->maxLinearSolveFailures) {
				snes->reason = SNES_DIVERGED_LINEAR_SOLVE;
				PetscFunctionReturn(0);
			}
		}
		ierr = KSPGetIterationNumber(snes->ksp,&lits);CHKERRQ(ierr);
		snes->linear_its += lits;
		ierr = PetscInfo2(snes,"iter=%D, linear solve iterations=%D\n",snes->iter,lits);CHKERRQ(ierr);
		
		
		/* compute residuals */
		
		
		/* 1) Compute u_i - u_{i-1} */
		VecAXPY( F, -1.0, X ); /* F <- F - X */
		ierr = VecNorm(F,NORM_2,&xnorm);CHKERRQ(ierr);	/* xnorm = || X_i - X_{i-1} || */
		
		/* 2) Compute F */
		/* ReBuild operator [A(u)], and rhs {b(u)} */
		ierr = SNESComputeFunction(snes,X,F);CHKERRQ(ierr);
		ierr = VecNorm(F,NORM_2,&fnorm);CHKERRQ(ierr);	/* fnorm = || F || */
		
		
		ierr = PetscObjectTakeAccess(snes);CHKERRQ(ierr);
		snes->iter = i+1;
		snes->norm = fnorm;
		ierr = PetscObjectGrantAccess(snes);CHKERRQ(ierr);
		SNESLogConvHistory(snes,fnorm,lits);
		SNESMonitor(snes,i+1,fnorm);
		
		
		/* check for convergence */
		if (snes->ops->converged) {
			ierr = (*snes->ops->converged)(snes,snes->iter,xnorm,1.0e6,fnorm,&snes->reason,snes->cnvP);CHKERRQ(ierr);
			if (snes->reason) break;
		}
		
		
		/* ReBuild operator [A(u)], and rhs {b(u)} */
		//SNESPicardComputeOperator(snes,X,&ctx->Amat,&ctx->Pmat,&flg);
		//SNESPicardComputeRhs(snes,X,ctx->rhs);
		
		i++;
		
	} while( i < snes->max_its );
	
	if (X != snes->vec_sol) {
		ierr = VecCopy(X,snes->vec_sol);CHKERRQ(ierr);
	}
	if (F != snes->vec_func) {
		ierr = VecCopy(F,snes->vec_func);CHKERRQ(ierr);
	}
//	snes->vec_sol_always  = snes->vec_sol;
//	snes->vec_func_always = snes->vec_func;
	
	
	if( i >= snes->max_its ) {
		ierr = PetscInfo1(snes,"Maximum number of iterations has been reached: %D\n",maxits);CHKERRQ(ierr);
		snes->reason = SNES_DIVERGED_MAX_IT;
	}
	
	
	PetscFunctionReturn(0);
}


/* -------------------------------------------------------------------------- */
/*
   SNESSetUp_Picard - Sets up the internal data structures for the later use
   of the SNESPicard nonlinear solver.

   Input Parameter:
.  snes - the SNES context

   Application Interface Routine: SNESSetUp()

   Notes:
   For basic use of the SNES solvers, the user need not explicitly call
   SNESSetUp(), since these actions will automatically occur during
   the call to SNESSolve().
 */
#undef __FUNCT__  
#define __FUNCT__ "SNESSetUp_PicardExt"
PetscErrorCode SNESSetUp_PicardExt(SNES snes)
{
	SNES_PICARD      *ctx = (SNES_PICARD*)snes->data;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	/* 
	We will accept a NULL function for either the compute_operator or compute_rhs 
	Either case indicates that this object is not non-linear.
	*/
	if (!ctx->Amat) {
		Stg_SETERRQ(PETSC_ERR_ARG_WRONGSTATE,"Must call SNESPicardSetOperator() before SNESSolve");
	}
	if (!ctx->rhs) {
		Stg_SETERRQ(PETSC_ERR_ARG_WRONGSTATE,"Must call SNESPicardSetRhs() before SNESSolve");
	}
	
	if( !ctx->F ) {
		ierr = VecDuplicate(snes->vec_sol, &ctx->F );CHKERRQ(ierr);
		
		/* define the ComputeFunction */
		SNESSetFunction( snes, ctx->F, _SNESPicard_ComputeFunction, (void*)ctx );
	}
	
	/* set picard solver to use last guess */
	KSPSetInitialGuessNonzero( snes->ksp, PETSC_TRUE );
	
	PetscFunctionReturn(0);
}

/* -------------------------------------------------------------------------- */
/*
   Stg_SNESDestroy_PicardExt - Destroys the private SNES_Picard context that was created
   with SNESCreate_PicardExt().

   Input Parameter:
.  snes - the SNES context

   Application Interface Routine: Stg_SNESDestroy( &)
 */
#undef __FUNCT__  
#define __FUNCT__ "Stg_SNESDestroy_PicardExt"
PetscErrorCode Stg_SNESDestroy_PicardExt(SNES snes)
{
	SNES_PICARD     *ctx = (SNES_PICARD*)snes->data;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	if( ctx->Amat )Stg_MatDestroy( & ctx->Amat );
	if( ctx->Pmat )Stg_MatDestroy( & ctx->Pmat );
	if( ctx->rhs )Stg_VecDestroy( & ctx->rhs );
	if( ctx->F )Stg_VecDestroy( & ctx->F );
	
	ierr = PetscFree(snes->data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/*
   SNESView_Picard - Prints info from the SNESPicard data structure.

   Input Parameters:
.  SNES - the SNES context
.  viewer - visualization context

   Application Interface Routine: SNESView()
*/
#undef __FUNCT__  
#define __FUNCT__ "SNESView_PicardExt"
static PetscErrorCode SNESView_PicardExt(SNES snes,PetscViewer viewer)
{
/*
	SNES_PICARD     *ctx = (SNES_PICARD*)snes->data;
	PetscErrorCode ierr;
*/	
	PetscFunctionBegin;
	
	PetscFunctionReturn(0);
}
/* -------------------------------------------------------------------------- */
/*
   SNESSetFromOptions_Picard - Sets various parameters for the SNESPicard method.

   Input Parameter:
.  snes - the SNES context

   Application Interface Routine: SNESSetFromOptions()
*/
#undef __FUNCT__  
#define __FUNCT__ "SNESSetFromOptions_Picard"
static PetscErrorCode SNESSetFromOptions_PicardExt(SNES snes)
{
/*
	SNES_PICARD     *ctx = (SNES_PICARD*)snes->data;
	PetscErrorCode ierr;
*/	
	PetscFunctionBegin;
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "SNESConverged_LS"
PetscErrorCode PETSCSNES_DLLEXPORT SNESConverged_PicardExt(SNES snes,PetscInt it,PetscReal xnorm,PetscReal pnorm,PetscReal fnorm,SNESConvergedReason *reason,void *dummy)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	PetscValidHeaderSpecific(snes,SNES_COOKIE,1);
	PetscValidType(snes,1);
	PetscValidPointer(reason,6);
	ierr = SNESDefaultConverged(snes,it,xnorm,pnorm,fnorm,reason,dummy);CHKERRQ(ierr);
	
	/* I don't like this diveregence test - it's too damm restrictive. */
	/* If you find convergence is not occuring, use -ksp_type richardson */
#if 0
	if (xnorm <= snes->abstol) {
		PetscInfo2(snes,"Converged due to norm(x_i - x_{i-1}) %G < %G (absolute tolerance)\n",xnorm,snes->abstol);
		*reason = 6699; // SNES_CONVERGED_DELTA_XNORM_ABS;
	}
#endif
	
	PetscFunctionReturn(0);
}

/* -------------------------------------------------------------------------- */
/*MC
      SNESPicard - Picard based nonlinear solver.

   Options Database:
+   -snes_picard_omega [positive number typically less than 1.0] - Relaxation parameter


   Level: beginner

.seealso:  SNESCreate(), SNES, SNESSetType()
M*/
EXTERN_C_BEGIN
PetscErrorCode SNESCreate_PicardExt(SNES snes)
{
	PetscErrorCode  ierr;
	SNES_PICARD      *neP;
	
	PetscFunctionBegin;
	snes->ops->setup           = SNESSetUp_PicardExt;
	snes->ops->solve           = SNESSolve_PicardExt;
	snes->ops->destroy         = Stg_SNESDestroy_PicardExt;
	snes->ops->converged       = SNESConverged_PicardExt;
	snes->ops->setfromoptions  = SNESSetFromOptions_PicardExt;
	snes->ops->view            = SNESView_PicardExt;
	snes->nwork                = 0;
	
	/* Create and set the Picard context */
	ierr       = PetscNew( SNES_PICARD, &neP );CHKERRQ(ierr);
	ierr       = PetscLogObjectMemory(snes,sizeof(SNES_PICARD));CHKERRQ(ierr);
	snes->data = (void*)neP;
	
	/* Init params on SNES_PICARD */
	neP->Amat     = PETSC_NULL;
	neP->Pmat     = PETSC_NULL;
	neP->rhs      = PETSC_NULL;
	neP->F        = PETSC_NULL;
	neP->op_ctx           = PETSC_NULL;
	neP->compute_operator = PETSC_NULL;
	neP->rhs_ctx          = PETSC_NULL;
	neP->compute_rhs      = PETSC_NULL;
	
	
	/* 
	Force a NULL function call for ComputeFunction and ComputeJacobian 
	so the SNESSetUp will not complain about it not being set.
	*/
//	SNESSetUp_Picard( snes );
//	ierr = SNESSetFunction(snes,PETSC_NULL,_SNESPicard_FormFunctionNULL,PETSC_NULL);CHKERRQ(ierr);
//	ierr = SNESSetJacobian(snes,PETSC_NULL,PETSC_NULL,_SNESPicard_FormJacobianNULL,PETSC_NULL);CHKERRQ(ierr);
	
	
	
	/* Public functions declared here */
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)snes,"SNESPicardSetOperator_C",
			"SNESPicardSetOperator_PICARD",
			SNESPicardSetOperator_PICARD);CHKERRQ(ierr);
	
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)snes,"SNESPicardSetRhs_C",
			"SNESPicardSetRhs_PICARD",
			SNESPicardSetRhs_PICARD);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}
EXTERN_C_END


