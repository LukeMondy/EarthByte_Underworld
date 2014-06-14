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
This test has been adopted from the petsc2.3.3 distribution.
File used was: src/snes/examples/tutorials/ex6.c
*/

static char help[] = "u`` + u^{2} = f. Different matrices for the Jacobian and the preconditioner.\n\
Demonstrates the use of matrix-free Newton-Krylov methods in conjunction\n\
with a user-provided preconditioner.  Input arguments are:\n\
   -snes_mf : Use matrix-free Newton methods\n\
   -user_precond : Employ a user-defined preconditioner.  Used only with\n\
                   matrix-free methods in this example.\n\n";

/*T
   Concepts: SNES^different matrices for the Jacobian and preconditioner;
   Concepts: SNES^matrix-free methods
   Concepts: SNES^user-provided preconditioner;
   Concepts: matrix-free methods
   Concepts: user-provided preconditioner;
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
*/
#include "petscsnes.h"
#include "petscext_snes.h"


/* 
   User-defined routines
*/
PetscErrorCode FormJacobian(SNES,Vec,Mat*,Mat*,MatStructure*,void*);
PetscErrorCode FormFunction(SNES,Vec,Vec,void*);
PetscErrorCode MatrixFreePreconditioner(PC,Vec,Vec);

int petsc_original( void )
{
  SNES           snes;                /* SNES context */
  KSP            ksp;                /* KSP context */
  PC             pc;                  /* PC context */
  Vec            x,r,F;               /* vectors */
  Mat            J,JPrec;             /* Jacobian,preconditioner matrices */
  PetscErrorCode ierr;
  PetscInt       it,n = 5,i;
  PetscMPIInt    size;
  PetscInt       *Shistit = 0,Khistl = 200,Shistl = 10;
  PetscReal      h,xp = 0.0,*Khist = 0,*Shist = 0;
  PetscScalar    v,pfive = .5;
  PetscTruth     flg;
	
	PetscPrintf( PETSC_COMM_WORLD, "In: %s \n", __func__ );
	
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size != 1) Stg_SETERRQ(1,"This is a uniprocessor example only!");
  ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&n,PETSC_NULL);CHKERRQ(ierr);
  h = 1.0/(n-1);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create nonlinear solver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create vector data structures; set function evaluation routine
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = VecCreate(PETSC_COMM_SELF,&x);CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&r);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&F);CHKERRQ(ierr);

  ierr = SNESSetFunction(snes,r,FormFunction,(void*)F);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create matrix data structures; set Jacobian evaluation routine
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,n,n,3,PETSC_NULL,&J);CHKERRQ(ierr);
  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,n,n,1,PETSC_NULL,&JPrec);CHKERRQ(ierr);

  /*
     Note that in this case we create separate matrices for the Jacobian
     and preconditioner matrix.  Both of these are computed in the
     routine FormJacobian()
  */
  ierr = SNESSetJacobian(snes,J,JPrec,FormJacobian,0);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Customize nonlinear solver; set runtime options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* Set preconditioner for matrix-free method */
  ierr = PetscOptionsHasName(PETSC_NULL,"-snes_mf",&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PetscOptionsHasName(PETSC_NULL,"-user_precond",&flg);CHKERRQ(ierr);
    if (flg) { /* user-defined precond */
      ierr = PCSetType(pc,PCSHELL);CHKERRQ(ierr);
      ierr = PCShellSetApply(pc,MatrixFreePreconditioner);CHKERRQ(ierr);
    } else {ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);}
  }

  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  /*
     Save all the linear residuals for all the Newton steps; this enables us
     to retain complete convergence history for printing after the conclusion
     of SNESSolve().  Alternatively, one could use the monitoring options
           -snes_monitor -ksp_monitor
     to see this information during the solver's execution; however, such
     output during the run distorts performance evaluation data.  So, the
     following is a good option when monitoring code performance, for example
     when using -log_summary.
  */
  ierr = PetscOptionsHasName(PETSC_NULL,"-rhistory",&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
    ierr = PetscMalloc(Khistl*sizeof(PetscReal),&Khist);CHKERRQ(ierr);
    ierr = KSPSetResidualHistory(ksp,Khist,Khistl,PETSC_FALSE);CHKERRQ(ierr);
    ierr = PetscMalloc(Shistl*sizeof(PetscReal),&Shist);CHKERRQ(ierr);
    ierr = PetscMalloc(Shistl*sizeof(PetscInt),&Shistit);CHKERRQ(ierr);
    ierr = SNESSetConvergenceHistory(snes,Shist,Shistit,Shistl,PETSC_FALSE);CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize application:
     Store right-hand-side of PDE and exact solution
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  xp = 0.0;
  for (i=0; i<n; i++) {
    v = 6.0*xp + pow(xp+1.e-12,6.0); /* +1.e-12 is to prevent 0^6 */
    ierr = VecSetValues(F,1,&i,&v,INSERT_VALUES);CHKERRQ(ierr);
    xp += h;
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Evaluate initial guess; then solve nonlinear system
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = VecSet(x,pfive);CHKERRQ(ierr);
  ierr = SNESSolve(snes,PETSC_NULL,x);CHKERRQ(ierr);
  ierr = SNESGetIterationNumber(snes,&it);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF,"Newton iterations = %D\n\n",it);CHKERRQ(ierr);
	
//	VecView( x, PETSC_VIEWER_STDOUT_SELF );
	FormFunction( snes, x, r, (void*)F );
	{ PetscReal Fnorm;
		VecNorm( r, NORM_2, &Fnorm );
		PetscPrintf( PETSC_COMM_WORLD, "Fnorm = %5.5e \n", Fnorm );
	}
	
	
  ierr = PetscOptionsHasName(PETSC_NULL,"-rhistory",&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = KSPGetResidualHistory(ksp,PETSC_NULL,&Khistl);CHKERRQ(ierr);
    ierr = PetscRealView(Khistl,Khist,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
    ierr = PetscFree(Khist);CHKERRQ(ierr);CHKERRQ(ierr);
    ierr = SNESGetConvergenceHistory(snes,PETSC_NULL,PETSC_NULL,&Shistl);CHKERRQ(ierr);
    ierr = PetscRealView(Shistl,Shist,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
    ierr = PetscIntView(Shistl,Shistit,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
    ierr = PetscFree(Shist);CHKERRQ(ierr);
    ierr = PetscFree(Shistit);CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = Stg_VecDestroy( &x);CHKERRQ(ierr);     ierr = Stg_VecDestroy( &r);CHKERRQ(ierr);
  ierr = Stg_VecDestroy( &F);CHKERRQ(ierr);     ierr = Stg_MatDestroy(J);CHKERRQ(ierr);
  ierr = Stg_MatDestroy(JPrec);CHKERRQ(ierr); ierr = Stg_SNESDestroy(snes);CHKERRQ(ierr);

  return 0;
}
/* ------------------------------------------------------------------- */
/* 
   FormInitialGuess - Forms initial approximation.

   Input Parameters:
   user - user-defined application context
   X - vector

   Output Parameter:
   X - vector
 */
PetscErrorCode FormFunction(SNES snes,Vec x,Vec f,void *dummy)
{
  PetscScalar    *xx,*ff,*FF,d;
  PetscErrorCode ierr;
  PetscInt       i,n;

  ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
  ierr = VecGetArray(f,&ff);CHKERRQ(ierr);
  ierr = VecGetArray((Vec)dummy,&FF);CHKERRQ(ierr);
  ierr = VecGetSize(x,&n);CHKERRQ(ierr);
  d = (PetscReal)(n - 1); d = d*d;
  ff[0]   = xx[0];
  for (i=1; i<n-1; i++) {
    ff[i] = d*(xx[i-1] - 2.0*xx[i] + xx[i+1]) + xx[i]*xx[i] - FF[i];
  }
  ff[n-1] = xx[n-1] - 1.0;
  ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);
  ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr);
  ierr = VecRestoreArray((Vec)dummy,&FF);CHKERRQ(ierr);
  return 0;
}
/* ------------------------------------------------------------------- */
/*
   FormJacobian - This routine demonstrates the use of different
   matrices for the Jacobian and preconditioner

   Input Parameters:
.  snes - the SNES context
.  x - input vector
.  ptr - optional user-defined context, as set by SNESSetJacobian()

   Output Parameters:
.  A - Jacobian matrix
.  B - different preconditioning matrix
.  flag - flag indicating matrix structure
*/
PetscErrorCode FormJacobian(SNES snes,Vec x,Mat *jac,Mat *prejac,MatStructure *flag,void *dummy)
{
  PetscScalar    *xx,A[3],d;
  PetscInt       i,n,j[3];
  PetscErrorCode ierr;

  ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
  ierr = VecGetSize(x,&n);CHKERRQ(ierr);
  d = (PetscReal)(n - 1); d = d*d;

  /* Form Jacobian.  Also form a different preconditioning matrix that 
     has only the diagonal elements. */
  i = 0; A[0] = 1.0; 
  ierr = MatSetValues(*jac,1,&i,1,&i,&A[0],INSERT_VALUES);CHKERRQ(ierr);
  ierr = MatSetValues(*prejac,1,&i,1,&i,&A[0],INSERT_VALUES);CHKERRQ(ierr);
  for (i=1; i<n-1; i++) {
    j[0] = i - 1; j[1] = i;                   j[2] = i + 1; 
    A[0] = d;     A[1] = -2.0*d + 2.0*xx[i];  A[2] = d; 
    ierr = MatSetValues(*jac,1,&i,3,j,A,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValues(*prejac,1,&i,1,&i,&A[1],INSERT_VALUES);CHKERRQ(ierr);
  }
  i = n-1; A[0] = 1.0; 
  ierr = MatSetValues(*jac,1,&i,1,&i,&A[0],INSERT_VALUES);CHKERRQ(ierr);
  ierr = MatSetValues(*prejac,1,&i,1,&i,&A[0],INSERT_VALUES);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(*jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(*prejac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*prejac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
//	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB_DENSE );
//	MatView( *jac, PETSC_VIEWER_STDOUT_SELF );
	
	
  ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);
  *flag = SAME_NONZERO_PATTERN;
  return 0;
}
/* ------------------------------------------------------------------- */
/*
   MatrixFreePreconditioner - This routine demonstrates the use of a
   user-provided preconditioner.  This code implements just the null
   preconditioner, which of course is not recommended for general use.

   Input Parameters:
.  ctx - optional user-defined context, as set by PCShellSetContext()
.  x - input vector

   Output Parameter:
.  y - preconditioned vector
*/
PetscErrorCode MatrixFreePreconditioner(PC pc,Vec x,Vec y)
{
  PetscErrorCode ierr;
  ierr = VecCopy(x,y);CHKERRQ(ierr);  
  return 0;
}







PetscErrorCode Picard_FormLinearRhs(SNES snes,Vec x,Vec f,void *dummy)
{
	PetscScalar    *xx,*ff,*FF,d;
	PetscErrorCode ierr;
	PetscInt       i,n;
	
	ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
	ierr = VecGetArray(f,&ff);CHKERRQ(ierr);
	ierr = VecGetArray((Vec)dummy,&FF);CHKERRQ(ierr);
	ierr = VecGetSize(x,&n);CHKERRQ(ierr);
	d = (PetscReal)(n - 1); d = d*d;
	ff[0]   = 0.0;
	for (i=1; i<n-1; i++) {
		// linear test
		ff[i] = FF[i];
	}
	ff[n-1] = 1.0;
	ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);
	ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr);
	ierr = VecRestoreArray((Vec)dummy,&FF);CHKERRQ(ierr);
	return 0;
}
PetscErrorCode Picard_FormRhs(SNES snes,Vec x,Vec f,void *dummy)
{
	PetscScalar    *xx,*ff,*FF,d;
	PetscErrorCode ierr;
	PetscInt       i,n;
	
	ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
	ierr = VecGetArray(f,&ff);CHKERRQ(ierr);
	ierr = VecGetArray((Vec)dummy,&FF);CHKERRQ(ierr);
	ierr = VecGetSize(x,&n);CHKERRQ(ierr);
	d = (PetscReal)(n - 1); d = d*d;
	ff[0]   = 0.0;
	for (i=1; i<n-1; i++) {
		// nonlinear test
		ff[i] = - xx[i]*xx[i] + FF[i];
	}
	ff[n-1] = 1.0;
	ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);
	ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr);
	ierr = VecRestoreArray((Vec)dummy,&FF);CHKERRQ(ierr);
	return 0;
}

PetscErrorCode FormLinearFunction(SNES snes,Vec x,Vec f,void *dummy)
{
	PetscScalar    *xx,*ff,*FF,d;
	PetscErrorCode ierr;
	PetscInt       i,n;
	
	ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
	ierr = VecGetArray(f,&ff);CHKERRQ(ierr);
	ierr = VecGetArray((Vec)dummy,&FF);CHKERRQ(ierr);
	ierr = VecGetSize(x,&n);CHKERRQ(ierr);
	d = (PetscReal)(n - 1); d = d*d;
	ff[0]   = xx[0];
	for (i=1; i<n-1; i++) {
		ff[i] = d*(xx[i-1] - 2.0*xx[i] + xx[i+1]) - FF[i];
	}
	ff[n-1] = xx[n-1] - 1.0;
	ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);
	ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr);
	ierr = VecRestoreArray((Vec)dummy,&FF);CHKERRQ(ierr);
	return 0;
}


PetscErrorCode Picard_FormOperator(SNES snes,Vec x,Mat *jac,Mat *prejac,MatStructure *flag,void *dummy)
{
	PetscScalar    *xx,A[3],d;
	PetscInt       i,n,j[3];
	PetscErrorCode ierr;
	PetscScalar    kappa;
	
	MatZeroEntries( *jac );
	MatZeroEntries( *prejac );
	
	ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
	ierr = VecGetSize(x,&n);CHKERRQ(ierr);
	d = (PetscReal)(n - 1); d = d*d;
	
	/* Form Jacobian.  Also form a different preconditioning matrix that 
	has only the diagonal elements. */
	i = 0; A[0] = 1.0; 
	ierr = MatSetValues(*jac,1,&i,1,&i,&A[0],INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues(*prejac,1,&i,1,&i,&A[0],INSERT_VALUES);CHKERRQ(ierr);
	for (i=1; i<n-1; i++) {
		j[0] = i - 1; j[1] = i;                   j[2] = i + 1; 
		kappa = xx[i]*xx[i];
		kappa = 1.0;
		
//		ff[i] = d*(xx[i-1] - 2.0*xx[i] + xx[i+1]) + xx[i]*xx[i] - FF[i];
		A[0] = d*kappa;     A[1] = -2.0*d*kappa;    A[2] = d*kappa; 
		ierr = MatSetValues(*jac,1,&i,3,j,A,INSERT_VALUES);CHKERRQ(ierr);
		ierr = MatSetValues(*prejac,1,&i,3,j,A,INSERT_VALUES);CHKERRQ(ierr);
	}
	i = n-1; A[0] = 1.0; 
	ierr = MatSetValues(*jac,1,&i,1,&i,&A[0],INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues(*prejac,1,&i,1,&i,&A[0],INSERT_VALUES);CHKERRQ(ierr);
	
	ierr = MatAssemblyBegin(*jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(*prejac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*prejac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	//	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB_DENSE );
	//	MatView( *jac, PETSC_VIEWER_STDOUT_SELF );
	
	
	ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);
	*flag = SAME_NONZERO_PATTERN;
	return 0;
}


PetscErrorCode HCPicard_FormOperator( Mat *jac,Mat *prejac,MatStructure *flag,void *dummy)
{
	PetscScalar    A[3],d;
	PetscInt       i,m,n,j[3];
	PetscErrorCode ierr;
	PetscScalar    kappa;
	
	MatZeroEntries( *jac );
	MatZeroEntries( *prejac );
	
	
	ierr = MatGetSize(*jac,&m,&n);CHKERRQ(ierr);
	d = (PetscReal)(n - 1); d = d*d;
	
	/* Form Jacobian.  Also form a different preconditioning matrix that 
	has only the diagonal elements. */
	i = 0; A[0] = 1.0; 
	ierr = MatSetValues(*jac,1,&i,1,&i,&A[0],INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues(*prejac,1,&i,1,&i,&A[0],INSERT_VALUES);CHKERRQ(ierr);
	for (i=1; i<n-1; i++) {
		j[0] = i - 1; j[1] = i;                   j[2] = i + 1; 
		kappa = 1.0;
		
//		ff[i] = d*(xx[i-1] - 2.0*xx[i] + xx[i+1]) + xx[i]*xx[i] - FF[i];
		A[0] = d*kappa;     A[1] = -2.0*d*kappa;    A[2] = d*kappa; 
		ierr = MatSetValues(*jac,1,&i,3,j,A,INSERT_VALUES);CHKERRQ(ierr);
		ierr = MatSetValues(*prejac,1,&i,1,&i,&A[1],INSERT_VALUES);CHKERRQ(ierr);
	}
	i = n-1; A[0] = 1.0; 
	ierr = MatSetValues(*jac,1,&i,1,&i,&A[0],INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues(*prejac,1,&i,1,&i,&A[0],INSERT_VALUES);CHKERRQ(ierr);
	
	ierr = MatAssemblyBegin(*jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(*prejac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*prejac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
//	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB_DENSE );
//	MatView( *jac, PETSC_VIEWER_STDOUT_SELF );
//	exit(0);
	
	*flag = SAME_NONZERO_PATTERN;
	return 0;
}



int picard_test( void )
{
	SNES           snes;                /* SNES context */
	KSP            ksp;                /* KSP context */
	PC             pc;                  /* PC context */
	Vec            x,rhs,FF;               /* vectors */
	Mat            Amat,Pmat;           /* Jacobian,preconditioner matrices */
	PetscErrorCode ierr;
	PetscInt       it,n = 5,i;
	PetscMPIInt    size;
	PetscReal      h,xp = 0.0;
	PetscScalar    v,pfive = .5;
	PetscTruth     flg;
	
	PetscPrintf( PETSC_COMM_WORLD, "In %s \n", __func__ );
	
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	if (size != 1) Stg_SETERRQ(1,"This is a uniprocessor example only!");
	ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&n,PETSC_NULL);CHKERRQ(ierr);
	h = 1.0/(n-1);
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Create nonlinear solver context
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	
	ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
	ierr = SNESSetType( snes, "picardext" );
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Create vector data structures; set function evaluation routine
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = VecCreate(PETSC_COMM_SELF,&x);CHKERRQ(ierr);
	ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
	ierr = VecSetFromOptions(x);CHKERRQ(ierr);
	
	ierr = VecDuplicate(x,&rhs);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&FF);CHKERRQ(ierr);
	
	
//	ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,n,n,3,PETSC_NULL,&Amat);CHKERRQ(ierr);
//	ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,n,n,3,PETSC_NULL,&Pmat);CHKERRQ(ierr);
	
	MatCreate( PETSC_COMM_SELF, &Amat );
	MatSetSizes( Amat, PETSC_DECIDE,PETSC_DECIDE, n,n );
	MatSetType( Amat, MATSEQAIJ );
	MatSeqAIJSetPreallocation( Amat, 3, PETSC_NULL );
	
	MatCreate( PETSC_COMM_SELF, &Pmat );
	MatSetSizes( Pmat, PETSC_DECIDE,PETSC_DECIDE, n,n );
	MatSetType( Pmat, MATSEQAIJ );
	MatSeqAIJSetPreallocation( Pmat, 3, PETSC_NULL );
	
	/*
	{
		MatStructure str;
		
		Picard_FormOperator( snes, x, &Amat,&Pmat,&str,PETSC_NULL );
		
		MatZeroEntries( Amat );
		MatZeroEntries( Pmat );
	}
	*/
	
	
//#if 0	
	SNESPicardSetOperator( snes,Amat,Pmat,		Picard_FormOperator, 	PETSC_NULL );
	
//	SNESPicardSetRhs( snes,FF,					Picard_FormLinearRhs, 	(void*)rhs );
	SNESPicardSetRhs( snes,FF,					Picard_FormRhs, 		(void*)rhs );
	
	
	
	ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
	
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Initialize application:
	Store right-hand-side of PDE and exact solution
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	
	xp = 0.0;
	for (i=0; i<n; i++) {
		v = 6.0*xp + pow(xp+1.e-12,6.0); /* +1.e-12 is to prevent 0^6 */
		ierr = VecSetValues(rhs,1,&i,&v,INSERT_VALUES);CHKERRQ(ierr);
		xp += h;
	}
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Evaluate initial guess; then solve nonlinear system
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	
	ierr = VecSet(x,pfive);CHKERRQ(ierr);
	/*
	{
		KSP _k;
		PC _p;
		MatStructure _flg;
		Vec _r;
		PetscReal nn;
		
		KSPCreate( PETSC_COMM_SELF, &_k );
		KSPSetOperators( _k, Amat,Pmat, DIFFERENT_NONZERO_PATTERN );
		_flg = DIFFERENT_NONZERO_PATTERN;
		SNESPicardComputeOperator( snes,x,&Amat,&Pmat,&_flg);
		SNESPicardComputeRhs( snes,x,FF);
		
		KSPSolve( _k, FF, x );
		//PetscPrintf( PETSC_COMM_WORLD, "linear solve \n");
		//VecView( x, PETSC_VIEWER_STDOUT_SELF );
		
		VecDuplicate( x, &_r );
		MatMult( Amat, x, _r );
		VecAXPY( _r, -1.0, FF );
		VecNorm( _r, NORM_2, &nn );
		PetscPrintf( PETSC_COMM_WORLD, "Linear norm = %5.5e \n", nn );
		
		ierr = VecSet(x,pfive);CHKERRQ(ierr);
	}
	*/
	
	
	
	ierr = SNESSolve(snes,PETSC_NULL,x);CHKERRQ(ierr);
//	VecView( x, PETSC_VIEWER_STDOUT_SELF );
	
	
	/*
	{ 
		Vec ffff;
		PetscReal Fnorm;
		//FormLinearFunction( snes, x, FF, (void*)rhs );
		//FormFunction( snes, x, FF, (void*)rhs );
		VecSet(FF,1.0);
		ffff = ((SNES_PICARD*)(snes->data))->F;
		
		SNESComputeFunction( snes, x, ffff );
		//SNESComputeFunction( snes, x, FF ); // this will not work. Must use the vector specified with SNESSetFunction()
		VecNorm( ffff, NORM_2, &Fnorm );
		PetscPrintf( PETSC_COMM_WORLD, "Fnorm(FF) = %e \n", Fnorm );
	}	
	*/
	{ 
		PetscReal Fnorm;
		Vec R, RHS;
		PetscErrorCode (*func)(SNES,Vec,Vec,void*);
		
		SNESGetFunction( snes, &R, &func, (void**)&RHS );
		func( snes,x,R,(void*)RHS );
		
		VecNorm( R, NORM_2, &Fnorm );
		PetscPrintf( PETSC_COMM_WORLD, "Fnorm(R) = %e \n", Fnorm );
	}
	
	
	ierr = SNESGetIterationNumber(snes,&it);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_SELF,"Newton iterations = %D\n\n",it);CHKERRQ(ierr);
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Free work space.  All PETSc objects should be destroyed when they
	are no longer needed.
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//#endif
	
	ierr = Stg_VecDestroy( &FF);CHKERRQ(ierr);
	ierr = Stg_VecDestroy( &x);CHKERRQ(ierr);
	ierr = Stg_VecDestroy( &rhs);CHKERRQ(ierr);
	ierr = Stg_MatDestroy(Amat);CHKERRQ(ierr);
	ierr = Stg_MatDestroy(Pmat);CHKERRQ(ierr);
	ierr = Stg_SNESDestroy(snes);CHKERRQ(ierr);
	
	return 0;
}
/* ------------------------------------------------------------------- */


int hand_coded_picard( void )
{
	KSP            ksp;                /* KSP context */
	PC             pc;                  /* PC context */
	Vec            x,rhs,S,r,xm,xdelta;               /* vectors */
	Mat            Amat,Pmat;           /* Jacobian,preconditioner matrices */
	PetscErrorCode ierr;
	PetscInt       it,n = 5,i;
	PetscMPIInt    size;
	PetscInt       *Shistit = 0,Khistl = 200,Shistl = 10;
	PetscReal      nn, h,xp = 0.0,*Khist = 0,*Shist = 0;
	PetscScalar    v,pfive = .5;
	PetscTruth     flg;
	MatStructure strflg;
	PetscScalar *_rhs, *_xx;
	
	PetscPrintf( PETSC_COMM_WORLD, "In %s \n", __func__ );
	
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	if (size != 1) Stg_SETERRQ(1,"This is a uniprocessor example only!");
	ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&n,PETSC_NULL);CHKERRQ(ierr);
	h = 1.0/(n-1);
	
	
	ierr = VecCreate(PETSC_COMM_SELF,&x);CHKERRQ(ierr);
	ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
	ierr = VecSetFromOptions(x);CHKERRQ(ierr);
	
	ierr = VecDuplicate(x,&rhs);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&S);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&r);CHKERRQ(ierr);
	
	ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,n,n,3,PETSC_NULL,&Amat);CHKERRQ(ierr);
	ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,n,n,3,PETSC_NULL,&Pmat);CHKERRQ(ierr);
	
	
	// setup source term
	xp = 0.0;
	for (i=0; i<n; i++) {
		v = 6.0*xp + pow(xp+1.e-12,6.0); /* +1.e-12 is to prevent 0^6 */
		ierr = VecSetValues(S,1,&i,&v,INSERT_VALUES);CHKERRQ(ierr);
		xp += h;
	}
	// setup 
	ierr = VecSet(x,pfive);CHKERRQ(ierr);
	
	KSPCreate( PETSC_COMM_SELF, &ksp );
	KSPSetFromOptions( ksp );
	KSPSetInitialGuessNonzero( ksp, PETSC_TRUE );
	
	strflg = DIFFERENT_NONZERO_PATTERN;
	
	for( it=0; it<20; it++ ) {
		/* build A,P, rhs */
		HCPicard_FormOperator( &Amat,&Pmat,&strflg,PETSC_NULL);
		
		VecCopy( S, rhs );
		
		ierr = VecGetArray(x,&_xx);CHKERRQ(ierr);
		ierr = VecGetArray(rhs,&_rhs);CHKERRQ(ierr);
		ierr = VecGetSize(x,&n);CHKERRQ(ierr);
		_rhs[0]   = 0.0;
		for (i=1; i<n-1; i++) {
			// nonlinear test
			_rhs[i] = _rhs[i] - _xx[i]*_xx[i];
		}
		_rhs[n-1] = 1.0;
		ierr = VecRestoreArray(x,&_xx);CHKERRQ(ierr);
		ierr = VecRestoreArray(rhs,&_rhs);CHKERRQ(ierr);
		
		
		
		
		KSPSetOperators( ksp, Amat, Amat, DIFFERENT_NONZERO_PATTERN );
		KSPSolve( ksp, rhs, x );
		
		// update F
		VecCopy( S, rhs );
		
		ierr = VecGetArray(x,&_xx);CHKERRQ(ierr);
		ierr = VecGetArray(rhs,&_rhs);CHKERRQ(ierr);
		ierr = VecGetSize(x,&n);CHKERRQ(ierr);
		_rhs[0]   = 0.0;
		for (i=1; i<n-1; i++) {
			// nonlinear test
			_rhs[i] = _rhs[i] - _xx[i]*_xx[i];
		}
		_rhs[n-1] = 1.0;
		ierr = VecRestoreArray(x,&_xx);CHKERRQ(ierr);
		ierr = VecRestoreArray(rhs,&_rhs);CHKERRQ(ierr);
		
		
		
		MatMult( Amat, x, r );
		VecAXPY( r, -1.0, rhs );
		VecNorm( r, NORM_2, &nn );
		PetscPrintf( PETSC_COMM_WORLD, "Linear norm = %5.5e \n", nn );
		
	}
//	VecView( x, PETSC_VIEWER_STDOUT_WORLD );
	
	
	
	PetscPrintf( PETSC_COMM_WORLD, "Applying corrections \n");
	VecDuplicate( x, &xm );
	VecDuplicate( x, &xdelta );
	VecCopy( x, xm );
	
	VecZeroEntries( xdelta );
	VecZeroEntries( x );
	
	VecCopy( xm, x );
	VecAXPY( x, 1.0, xdelta );
	
	
	for( it=0; it<20; it++ ) {
		
		/* build A,P, rhs */
		HCPicard_FormOperator( &Amat,&Pmat,&strflg,PETSC_NULL);
		
		VecCopy( S, rhs );
		
		ierr = VecGetArray(x,&_xx);CHKERRQ(ierr);
		ierr = VecGetArray(rhs,&_rhs);CHKERRQ(ierr);
		ierr = VecGetSize(x,&n);CHKERRQ(ierr);
		_rhs[0]   = 0.0;
		for (i=1; i<n-1; i++) {
			// nonlinear test
			_rhs[i] = _rhs[i] - _xx[i]*_xx[i];
		}
		_rhs[n-1] = 1.0;
		ierr = VecRestoreArray(x,&_xx);CHKERRQ(ierr);
		ierr = VecRestoreArray(rhs,&_rhs);CHKERRQ(ierr);
		
		
		MatMult( Amat, x, r );
		VecAYPX( r, -1.0, rhs ); /* r <- b - Ax */
	
		
		VecZeroEntries( xdelta );
		KSPSetOperators( ksp, Amat, Amat, DIFFERENT_NONZERO_PATTERN );
		KSPSolve( ksp, r, xdelta );
		VecCopy( xm, x );
		VecAXPY( x, 1.0, xdelta );
		
		// update F
		VecCopy( S, rhs );
		
		ierr = VecGetArray(x,&_xx);CHKERRQ(ierr);
		ierr = VecGetArray(rhs,&_rhs);CHKERRQ(ierr);
		ierr = VecGetSize(x,&n);CHKERRQ(ierr);
		_rhs[0]   = 0.0;
		for (i=1; i<n-1; i++) {
			// nonlinear test
			_rhs[i] = _rhs[i] - _xx[i]*_xx[i];
		}
		_rhs[n-1] = 1.0;
		ierr = VecRestoreArray(x,&_xx);CHKERRQ(ierr);
		ierr = VecRestoreArray(rhs,&_rhs);CHKERRQ(ierr);
		
		
		
		MatMult( Amat, x, r );
		VecAXPY( r, -1.0, rhs );
		VecNorm( r, NORM_2, &nn );
		PetscPrintf( PETSC_COMM_WORLD, "Linear norm = %5.5e \n", nn );
		
	}
//	VecView( x, PETSC_VIEWER_STDOUT_WORLD );
	
	
	
	
	
	
	
	
	
	
	
	ierr = Stg_VecDestroy( &x);CHKERRQ(ierr);
	ierr = Stg_VecDestroy( &rhs);CHKERRQ(ierr);
	ierr = Stg_MatDestroy(Amat);CHKERRQ(ierr);
	ierr = Stg_MatDestroy(Pmat);CHKERRQ(ierr);
	
	return 0;
}
/* ------------------------------------------------------------------- */





int main(int argc,char **argv)
{
	PetscErrorCode ierr;
	int i;
	
	
	PetscInitialize(&argc,&argv,(char *)0,help);
	PetscExtSNESRegisterAll();
	
//	petsc_original();
	
	
	for( i=0; i<1; i++ ) {
		picard_test(); 
	}
	
	
	/* 
	==== RECOMMEND YOU USE RICHARDSON ====
	./ex6.app -snes_monitor -n 100 -ksp_type richardson -snes_atol 1.0e-8 -snes_converged_reason 
	*/
	
//	hand_coded_picard();
	
	
	
	ierr = PetscExtSNESRegisterDestroyAll();CHKERRQ(ierr);
	ierr = PetscFinalize();CHKERRQ(ierr);
	
	return 0;
}

