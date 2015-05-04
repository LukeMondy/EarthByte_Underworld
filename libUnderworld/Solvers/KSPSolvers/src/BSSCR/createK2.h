#ifndef _createK2_h
#define _createK2_h

#include <petsc.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>

PetscErrorCode bsscr_buildK2(KSP ksp);
PetscErrorCode bsscr_DGMiGtD( Mat *_K2, Mat K, Mat G, Mat M);
PetscErrorCode bsscr_GMiGt( Mat *_K2, Mat K, Mat G, Mat M);
PetscErrorCode bsscr_GGt( Mat *_K2, Mat K, Mat G);
#endif
