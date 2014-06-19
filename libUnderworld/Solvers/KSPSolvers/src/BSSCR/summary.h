#ifdef HAVE_PETSCEXT
#ifndef __BSSCR_SUM_h__
#define __BSSCR_SUM_h__

#include <petsc.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>
#include <petscpc.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include "Solvers/KSPSolvers/KSPSolvers.h"

#include "common-driver-utils.h"
#include "BSSCR.h" /* includes StokesBlockKSPInterface.h */

//void bsscr_summary(KSP_BSSCR * bsscrp_self, KSP ksp_S, KSP ksp_inner, Mat K,Mat Korig,Mat K2,Mat D,Mat G,Mat C,Vec u,Vec p,Vec f,Vec h,Vec t,

void bsscr_summary(KSP_BSSCR * bsscrp_self, KSP ksp_S, KSP ksp_inner, Mat K,Mat K2,Mat D,Mat G,Mat C,Vec u,Vec p,Vec f,Vec h,Vec t,
		   double penaltyNumber,PetscTruth KisJustK,double mgSetupTime,double scrSolveTime,double a11SingleSolveTime);

#endif
#endif
