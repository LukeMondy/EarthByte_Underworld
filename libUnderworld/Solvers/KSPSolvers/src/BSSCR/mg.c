#ifdef HAVE_PETSCEXT
#include "mg.h"

/* PetscErrorCode BSSCR_mgPCApply( void *_ctx, Vec x, Vec y ) { */
/*       MGContext *ctx = (MGContext*)_ctx; */
/*       PetscInt curIt; */
/*       PetscTruth found, flg; */
/*       static double wallTime = 0.0; */
/*       PetscInt numSmooths, numIts; */

/*       KSPGetIterationNumber( ctx->ksp, &curIt ); */

/*       /\* If we're on the first iteration, reset everything. *\/ */
/*       if( curIt == 0 ) { */
/* 	    ctx->curIt = -1; */
/* 	    ctx->myIts = 0; */
/* 	    ctx->remain = 0; */
/* 	    wallTime = MPI_Wtime();		 */
/*       } */

/*       /\* If we are still on the same iteration, don't check everything again. *\/ */

/*       if( curIt != ctx->curIt ) { */
/* 	    ctx->curIt = curIt; */

/* 	    if( ctx->remain ) */
/* 		  ctx->remain--; */



/* 	    else { */

/* 		  numSmooths = (PetscInt)(((PetscScalar)ctx->smoothsMin)*exp( ((PetscScalar)ctx->myIts)/ctx->smoothsFac )); */
/* 		  numIts = (PetscInt)(((PetscScalar)ctx->itsMax)*exp( ((PetscScalar)(-ctx->myIts))/ctx->itsFac )); */

/* 		  if( numSmooths > ctx->smoothsMax ) */
/* 			numSmooths = ctx->smoothsMax; */
/* 		  if( numIts < ctx->itsMin ) */
/* 			numIts = ctx->itsMin; */

/* 		  flg = PETSC_FALSE; */
/* 		  PetscOptionsGetTruth( PETSC_NULL, "-mg_fancy_view", &flg, &found ); */
/* 		  if( flg ) { */
/* 			PetscPrintf( PETSC_COMM_WORLD,   */
/* 		                     "%f seconds - update to use %4d smoothing cycles for the next %3d iterations\n", */
/* 				     MPI_Wtime() - wallTime, numSmooths, numIts ); */
/* 		  } */

/* 		  wallTime = MPI_Wtime(); */

/* 		  PCMGSetNumberSmoothDown( ctx->pc, numSmooths ); */
/* 		  PCMGSetNumberSmoothUp(   ctx->pc, numSmooths );  */
/* 		  ctx->remain = numIts; */
/* 		  ctx->myIts++; */
/* 	    } */
/*       } */

/*       PCApply( ctx->pc, x, y ); */

/*       PetscFunctionReturn( 0 ); */
/* } */

/* PetscErrorCode BSSCR_mgPCAccelerating( void *_ctx, Vec x, Vec y ) { */
/*       MGContext *ctx = (MGContext*)_ctx; */
/*       PetscErrorCode          ierr; */
/*       PetscInt curIt; */
/*       PetscInt numSmooths; */
/*       PetscTruth found, flg; */
/*       PetscReal Rnorm; */
/*       PetscReal logResidualRatio; */
/*       PetscTruth useNormInfStoppingConditions; */
	
/*       static double wallTime = 0.0; */
/*       Vec  R, w1, w2, work; */

/*       KSPGetIterationNumber( ctx->ksp, &curIt ); */
	
/*       /\* 	Why does KSPGetResidualNorm always return the 2 Norm ?  */
	
/* 	   	There's no choice but to over-ride this by hand, but, my guess is that, */
/* 	 	for some solvers this is going to be quite expensive - caveat emptor ! */
	
/*       *\/ */
	
	
/*       ierr = VecDuplicate(ctx->ksp->vec_rhs,&work);CHKERRQ(ierr); */
/*       ierr = VecDuplicate(ctx->ksp->vec_rhs,&w1);CHKERRQ(ierr); */
/*       ierr = VecDuplicate(ctx->ksp->vec_rhs,&w2);CHKERRQ(ierr);	 */
    
/*       KSPBuildResidual( ctx->ksp, w1,w2, &R ); */
	
/*       useNormInfStoppingConditions = PETSC_FALSE; */
/*       PetscOptionsGetTruth( PETSC_NULL ,"-A11_use_norm_inf_stopping_condition", &useNormInfStoppingConditions, &found ); */

/*       if(useNormInfStoppingConditions)  */
/* 	    VecNorm( R, NORM_INFINITY, &Rnorm ); */
/*       else */
/* 	    VecNorm( R, NORM_2, &Rnorm ); */


/*       //KSPGetResidualNorm( ctx->ksp, &Rnorm ); */

/*       /\* If we're on the first iteration, reset everything. *\/ */

/*       if( curIt == 0 ) { */
/* 	    ctx->curIt = -1; */
/* 	    ctx->myIts = 0; */
/* 	    ctx->remain = ctx->targetCyclesForTenfoldReduction; */
/* 	    ctx->residualToWatch = Rnorm; */
/* 	    wallTime = MPI_Wtime();	 */
		
/* 	    /\* Reset this each time if no tests are going to be done (smooths always increasing with iteration count) *\/ */
		
/* 	    if(!ctx->smoothingAccelerationTest)	 */
/* 		  ctx->currentNumberOfSmooths = ctx->smoothsToStartWith; */
		
	   
/* 	    flg = PETSC_FALSE;   */
/* 	    PetscOptionsGetTruth( PETSC_NULL, "-mg_accelerating_smoothing_view", &flg, &found ); */

/* 	    if( flg ) { */
/* 		  PetscPrintf( PETSC_COMM_WORLD,   */
/* 			       "MG:  Use %d smoothing cycles for the first %d iterations \n", */
/* 			       ctx->currentNumberOfSmooths, ctx->targetCyclesForTenfoldReduction); */
/* 	    } */
		
/*       } */

/*       /\* If we are still on the same iteration, don't need to check everything again. *\/ */

/*       if( curIt != ctx->curIt ) { */
/* 	    ctx->curIt = curIt; */

/* 	    /\* 	Check how we are going after a fixed number of iterations and */
/* 		make adjustments to the smoothing parameters accordingly *\/ */

/* 	    if( ctx->remain ) */
/* 		  ctx->remain--; */

/* 	    else { */
	  
/* 		  /\* How does the current residual match against the target ? *\/ */
	
/* 		  logResidualRatio = log10( ctx->residualToWatch / Rnorm );   /\* Should check if either of these is zero,  */
/* 										 but why would we be here if it was ? *\/ */
	
/* 		  numSmooths = ctx->currentNumberOfSmooths; */
			
/* 		  /\* 1) We hit the target very easily *\/  */
			
/* 		  if(ctx->smoothingAccelerationTest && logResidualRatio > 2.0 ) { */
/* 			numSmooths = (int) ((PetscScalar) ctx->currentNumberOfSmooths / ctx->smoothingAcceleration) - ctx->smoothingIncrement; */
						
/* 			if( numSmooths < ctx->smoothsMin ) */
/* 			      numSmooths = ctx->smoothsMin;	 */
/* 		  } */
			
/* 		  /\* 2) We scrape the target and therefore need more iterations, or we */
/* 		     don't bother testing and assume we are not making enough progress  */
/* 		     on the basis of having already done several iterations *\/  */
			
/* 		  if(!ctx->smoothingAccelerationTest || logResidualRatio < 1.0 ) { */
/* 			numSmooths = (int) ((PetscScalar) ctx->currentNumberOfSmooths * ctx->smoothingAcceleration) + ctx->smoothingIncrement; */
			
/* 			if( numSmooths > ctx->smoothsMax ) */
/* 			      numSmooths = ctx->smoothsMax;	 */
/* 		  } */
			
/* 		  flg = PETSC_FALSE;   */
/* 		  PetscOptionsGetTruth( PETSC_NULL, "-mg_accelerating_smoothing_view", &flg, &found ); */
/* 		  if( flg ) { */
/* 			PetscPrintf( PETSC_COMM_WORLD,   */
/* 		                     "MG: Update to use %d smoothing cycles for the next %d iterations ( %fs to reduce res by 10^%f [%f to %f])\n", */
/* 				     numSmooths, ctx->targetCyclesForTenfoldReduction, MPI_Wtime() - wallTime, */
/* 				     logResidualRatio, ctx->residualToWatch, Rnorm ); */
/* 		  } */

/* 		  wallTime = MPI_Wtime(); */
/* 		  PCMGSetNumberSmoothDown( ctx->pc, numSmooths ); */
/* 		  PCMGSetNumberSmoothUp(   ctx->pc, numSmooths );  */
/* 		  ctx->remain = ctx->targetCyclesForTenfoldReduction; */
/* 		  ctx->currentNumberOfSmooths = numSmooths; */
/* 		  ctx->myIts++; */
/* 		  ctx->residualToWatch = Rnorm; */
/* 	    } */
/*       } */

/*       PCApply( ctx->pc, x, y ); */

/*       Stg_VecDestroy(&work); */
/*       Stg_VecDestroy(&w1); */
/*       Stg_VecDestroy(&w2); */

/*       PetscFunctionReturn( 0 ); */
/* } */
/******************************************************************************************************/
PetscErrorCode MG_inner_solver_mgContext_initialise(MGContext *mgCtx) {

    PetscTruth found;

    /* 1) Default values into the mg context */

    mgCtx->useAcceleratingSmoothingMG   = PETSC_FALSE;
    mgCtx->acceleratingSmoothingMGView  = PETSC_FALSE;

    mgCtx->smoothsToStartWith = 5;
    mgCtx->smoothsMax = 250;
    mgCtx->smoothingIncrement = 2;
    mgCtx->targetCyclesForTenfoldReduction = 3;
    mgCtx->totalMgCycleCount=0;
    mgCtx->totalSmoothingCount=0;

    /* 2) Read options from petsc dictionary */
    	
    PetscOptionsGetTruth( PETSC_NULL, "-mg_accelerating_smoothing",               &mgCtx->useAcceleratingSmoothingMG,  &found );
    PetscOptionsGetTruth( PETSC_NULL, "-mg_accelerating_smoothing_view",          &mgCtx->acceleratingSmoothingMGView, &found );
    PetscOptionsGetInt(   PETSC_NULL, "-mg_smooths_to_start",                     &(mgCtx->smoothsToStartWith),         PETSC_NULL );
    PetscOptionsGetInt(   PETSC_NULL, "-mg_smooths_max",                          &(mgCtx->smoothsMax),                 PETSC_NULL );   
    PetscOptionsGetInt(   PETSC_NULL, "-mg_smoothing_increment",                  &(mgCtx->smoothingIncrement),              PETSC_NULL );
    PetscOptionsGetInt(   PETSC_NULL, "-mg_target_cycles_10fold_reduction",       &(mgCtx->targetCyclesForTenfoldReduction), PETSC_NULL );
    
    mgCtx->currentNumberOfSmooths = mgCtx->smoothsToStartWith;

    if(mgCtx->acceleratingSmoothingMGView && mgCtx->useAcceleratingSmoothingMG) {
	PetscPrintf( PETSC_COMM_WORLD,  "\nSCR MG OPTIONS\n" );
	PetscPrintf( PETSC_COMM_WORLD,  "  Smooths to start with:     %d\n", mgCtx->smoothsToStartWith );
	PetscPrintf( PETSC_COMM_WORLD,  "  Smooths maximum:           %d\n", mgCtx->smoothsMax );
	PetscPrintf( PETSC_COMM_WORLD,  "  Smoothing increment:       %d\n", mgCtx->smoothingIncrement );
	PetscPrintf( PETSC_COMM_WORLD,  "  Iterations between checks: %d\n\n", mgCtx->targetCyclesForTenfoldReduction );
    }	

    PetscFunctionReturn( 0 );
}
/* PetscErrorCode MG_inner_solver_pcmg_create( KSP ksp_inner, PC *pc_MG ) { */
/*     /\* Create PCs *\/ */
/*     PCCreate(PETSC_COMM_WORLD, pc_MG); */
/*     KSPSetPC(ksp_inner, *pc_MG); */
/*     PetscFunctionReturn( 0 ); */
/* } */
PetscErrorCode MG_inner_solver_pcmg_setup( KSP_BSSCR * bsscrp_self, MGContext *mgCtx, KSP ksp_inner, PC pc_MG, Mat K ) {

    PCSetType(pc_MG, PCMG);
    PCMGSetLevels(pc_MG, bsscrp_self->mg->nLevels, PETSC_NULL);
    PCMGSetType(pc_MG, PC_MG_MULTIPLICATIVE);
    #if ((PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR>=2) )
    PCMGSetGalerkin( pc_MG, PETSC_TRUE );
    #else
    PCMGSetGalerkin( pc_MG );
    #endif
    PCSetFromOptions(pc_MG);
    
    KSPSetOperators(ksp_inner, K, K, DIFFERENT_NONZERO_PATTERN);
    
    bsscrp_self->mg->mgData->ksp = ksp_inner;
	
    mgCtx->ksp = ksp_inner;
    mgCtx->pc = pc_MG;
 	
    PETScMGSolver_UpdateOps(bsscrp_self->mg);

    /* If we are using dynamically adjusting smoother settings then 
       this is implemented as a KSPMonitor (with side-effects)*/	

    if( mgCtx->useAcceleratingSmoothingMG  ) 
	KSPMonitorSet( ksp_inner, KSPCycleEffectivenessMonitorAndAdjust, mgCtx, PETSC_NULL );

    PetscFunctionReturn(0);

}
PetscErrorCode MG_inner_solver_pcmg_shutdown( PC pc_MG ) {
    PetscErrorCode ierr;
    ierr = Stg_PCDestroy(&pc_MG );CHKERRQ(ierr);
    PetscFunctionReturn(0);
}
#undef __FUNCT__  
#define __FUNCT__ "KSPCycleEffectivenessMonitorAndAdjust"
PetscErrorCode KSPCycleEffectivenessMonitorAndAdjust(KSP ksp, PetscInt n, PetscReal rnorm, MGContext *mgctx )
{
    PetscErrorCode          ierr;
    //PetscViewerASCIIMonitor viewer;
    PetscViewerASCIIMonitor viewer  = PETSC_VIEWER_STDOUT_(((PetscObject)ksp)->comm);			
    PetscFunctionBegin;
	
    if(n==0) {
	mgctx->smoothingCountThisSolve=0;
	mgctx->currentNumberOfSmooths = mgctx->smoothsToStartWith;
    }
    else if((n % mgctx->targetCyclesForTenfoldReduction) == 0 && mgctx->currentNumberOfSmooths < mgctx->smoothsMax) {
	mgctx->currentNumberOfSmooths += mgctx->smoothingIncrement; 		
    }
				
    ierr = PetscViewerASCIIMonitorCreate(((PetscObject)ksp)->comm,"stdout",0,&viewer); 
    CHKERRQ(ierr);			
	
    if(mgctx->acceleratingSmoothingMGView) {
	ierr = PetscViewerASCIIMonitorPrintf(viewer,"%3D MG smoothing cycles %d [%d] \n", n, 
					     mgctx->currentNumberOfSmooths,mgctx->smoothingCountThisSolve); 
	CHKERRQ(ierr);
    }
	
    PCMGSetNumberSmoothDown( mgctx->pc, mgctx->currentNumberOfSmooths );
    PCMGSetNumberSmoothUp(   mgctx->pc, mgctx->currentNumberOfSmooths ); 
	
    mgctx->totalMgCycleCount++;
    mgctx->smoothingCountThisSolve += mgctx->currentNumberOfSmooths;
    mgctx->totalSmoothingCount += mgctx->currentNumberOfSmooths;		
																		
    ierr = PetscViewerASCIIMonitorDestroy(viewer);
    CHKERRQ(ierr);
		
    PetscFunctionReturn(0);
}
double setupMG( KSP_BSSCR * bsscrp_self, KSP ksp_inner, PC pc_MG, Mat K, MGContext *mgCtx ){
    double mgSetupTime;

    mgSetupTime = MPI_Wtime();

    MG_inner_solver_mgContext_initialise( mgCtx ); 
    //MG_inner_solver_pcmg_create( ksp_inner, &pc_MG );
    KSPSetOptionsPrefix( ksp_inner, "A11_" );/* just in case.. for the moment */
    KSPSetFromOptions( ksp_inner );
    PetscPrintf( PETSC_COMM_WORLD,  "Populate MG PC \n");
    MG_inner_solver_pcmg_setup( bsscrp_self, mgCtx, ksp_inner, pc_MG, K );
    PetscPrintf( PETSC_COMM_WORLD,  "Populate MG PC ... done \n");

    mgSetupTime = MPI_Wtime() - mgSetupTime;

    return mgSetupTime;
}
/********************************************************************************************************/
/* double XXsetupMG( KSP_BSSCR * bsscrp_self, KSP ksp_inner, PC pc_MG, Mat K, PC *_shellPC){ */
/*     double mgSetupTime; */
/*     PetscTruth useFancySmoothingMG, useAcceleratingSmoothingMG, flg, found; */
/*     PC shellPC; */
/*     MGContext mgCtx; */
		
/*     mgSetupTime = MPI_Wtime(); */
		
/*     PCSetType(pc_MG, PCMG); */
/*     PCMGSetLevels(pc_MG, bsscrp_self->mg->nLevels, PETSC_NULL); */
/*     PCMGSetType(pc_MG, PC_MG_MULTIPLICATIVE); */
/*     PCMGSetGalerkin(pc_MG); */
/*     PCSetFromOptions(pc_MG); */

/*     KSPSetOperators(ksp_inner, K, K, DIFFERENT_NONZERO_PATTERN); */

/*     bsscrp_self->mg->mgData->ksp = ksp_inner; */
/*     PETScMGSolver_UpdateOps(bsscrp_self->mg); */
/*     /\* Do we want to use dynamically adjusting smoother settings ?  *\/ */
/*     useFancySmoothingMG = PETSC_FALSE; */
/*     PetscOptionsGetTruth( PETSC_NULL, "-mg_fancy", &useFancySmoothingMG, &found ); */

/*     useAcceleratingSmoothingMG = PETSC_FALSE; */
/*     PetscOptionsGetTruth( PETSC_NULL, "-mg_accelerating_smoothing", &useAcceleratingSmoothingMG, &found ); */

/*     if( useAcceleratingSmoothingMG || useFancySmoothingMG ) { */
/* 	/\* Create a shell preconditioner so we can influence the multigrid */
/* 	   options during iterations. *\/ */
/* 	PCCreate( PETSC_COMM_WORLD, &shellPC ); */
/* 	PCSetType( shellPC, PCSHELL ); */
/* 	if( useAcceleratingSmoothingMG) */
/* 	    PCShellSetApply( shellPC, BSSCR_mgPCAccelerating ); */
/* 	else */
/* 	    PCShellSetApply( shellPC, BSSCR_mgPCApply ); */
				
/* 	PCShellSetContext( shellPC, &mgCtx ); */
/* 	/\* Insert the shell PC. *\/ */
/* 	PetscObjectReference( (PetscObject)pc_MG ); */
/* 	KSPSetPC( ksp_inner, shellPC ); */
/* 	KSPSetOperators(ksp_inner, K, K, DIFFERENT_NONZERO_PATTERN); */
/* 	/\* Setup the MG context. */
/* 	   Default values, then values taken from the petsc options database *\/ */
/* 	mgCtx.ksp = ksp_inner; */
/* 	mgCtx.pc = pc_MG; */
/* 	mgCtx.smoothsMin = 1; */
/* 	mgCtx.smoothsToStartWith = 5; */
/* 	mgCtx.smoothsMax = 250; */
/* 	mgCtx.smoothsFac = 1.0; */
/* 	mgCtx.itsMin     = 3; */
/* 	mgCtx.itsMax     = 10; */
/* 	mgCtx.itsFac     = 2.0; */
/* 	mgCtx.smoothingIncrement = 2; */
/* 	mgCtx.smoothingAcceleration = 1.0; */
/* 	mgCtx.targetCyclesForTenfoldReduction = 3; */
/* 	mgCtx.smoothingAccelerationTest = PETSC_TRUE; */

/* 	PetscOptionsGetInt( PETSC_NULL,    "-mg_smooths_to_start", &mgCtx.smoothsToStartWith, PETSC_NULL ); */
/* 	PetscOptionsGetInt( PETSC_NULL,    "-mg_smooths_min", &mgCtx.smoothsMin, PETSC_NULL ); */
/* 	PetscOptionsGetInt( PETSC_NULL,    "-mg_smooths_max", &mgCtx.smoothsMax, PETSC_NULL ); */
/* 	PetscOptionsGetScalar( PETSC_NULL, "-mg_smooths_factor", &mgCtx.smoothsFac, PETSC_NULL ); */
/* 	PetscOptionsGetInt( PETSC_NULL,    "-mg_its_min", &mgCtx.itsMin, PETSC_NULL ); */
/* 	PetscOptionsGetInt( PETSC_NULL,    "-mg_its_max", &mgCtx.itsMax, PETSC_NULL ); */
/* 	PetscOptionsGetScalar( PETSC_NULL, "-mg_its_factor", &mgCtx.itsFac, PETSC_NULL ); */

/* 	PetscOptionsGetInt(    PETSC_NULL, "-mg_smoothing_increment",            &mgCtx.smoothingIncrement,              PETSC_NULL ); */
/* 	PetscOptionsGetInt(    PETSC_NULL, "-mg_target_cycles_10fold_reduction", &mgCtx.targetCyclesForTenfoldReduction, PETSC_NULL ); */
/* 	PetscOptionsGetScalar( PETSC_NULL, "-mg_smoothing_acceleration",         &mgCtx.smoothingAcceleration,           PETSC_NULL ); */
/* 	PetscOptionsGetTruth(  PETSC_NULL, "-mg_smoothing_adjust_on_convergence_rate",   &mgCtx.smoothingAccelerationTest,       &found ); */

/* 	mgCtx.currentNumberOfSmooths = mgCtx.smoothsToStartWith; */

/* 	flg = PETSC_FALSE; PetscOptionsGetTruth( PETSC_NULL, "-mg_fancy_view", &flg, &found ); */
/* 	if( flg ) { */
/* 	    PetscPrintf( PETSC_COMM_WORLD,  "\nSCR MG OPTIONS\n" ); */
/* 	    PetscPrintf( PETSC_COMM_WORLD,  "  Smooths minimum:    %d\n", mgCtx.smoothsMin ); */
/* 	    PetscPrintf( PETSC_COMM_WORLD,  "  Smooths maximum:    %d\n", mgCtx.smoothsMax ); */
/* 	    PetscPrintf( PETSC_COMM_WORLD,  "  Smooths factor:     %g\n", mgCtx.smoothsFac ); */
/* 	    PetscPrintf( PETSC_COMM_WORLD,  "  Iterations minimum: %d\n", mgCtx.itsMin ); */
/* 	    PetscPrintf( PETSC_COMM_WORLD,  "  Iterations maximum: %d\n", mgCtx.itsMax ); */
/* 	    PetscPrintf( PETSC_COMM_WORLD,  "  Iterations factor:  %g\n\n", mgCtx.itsFac ); */
/* 	}	 */

/* 	flg = PETSC_FALSE; */
/* 	PetscOptionsGetTruth( PETSC_NULL, "-mg_accelerating_smoothing_view", &flg, &found ); */
/* 	if( flg ) { */
/* 	    PetscPrintf( PETSC_COMM_WORLD,  "\nSCR MG OPTIONS\n" ); */
/* 	    PetscPrintf( PETSC_COMM_WORLD,  "  Smooths to start with:     %d\n", mgCtx.smoothsToStartWith ); */
/* 	    PetscPrintf( PETSC_COMM_WORLD,  "  Smooths minimum:           %d\n", mgCtx.smoothsMin ); */
/* 	    PetscPrintf( PETSC_COMM_WORLD,  "  Smooths maximum:           %d\n", mgCtx.smoothsMax ); */
/* 	    PetscPrintf( PETSC_COMM_WORLD,  "  Smoothing increment (or):  %d\n", mgCtx.smoothingIncrement ); */
/* 	    PetscPrintf( PETSC_COMM_WORLD,  "  Smoothing acceleration:    %f\n", mgCtx.smoothingAcceleration ); */
/* 	    PetscPrintf( PETSC_COMM_WORLD,  "  Iterations between checks: %d\n\n", mgCtx.targetCyclesForTenfoldReduction ); */
/* 	}	 */
/*     }  */
/*     mgSetupTime = MPI_Wtime() - mgSetupTime; */
/*     *_shellPC=shellPC;				 */
/*     return mgSetupTime; */
/* } */
#endif
