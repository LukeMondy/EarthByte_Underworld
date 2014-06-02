#ifdef HAVE_PETSCEXT
#ifndef __KSP_PNS_h__
#define __KSP_PNS_h__
PetscErrorCode KSPBuildPressure_CB_Nullspace_BSSCR(KSP ksp);
PetscErrorCode KSPBuildPressure_Const_Nullspace_BSSCR(KSP ksp);
PetscErrorCode KSPRemovePressureNullspace_BSSCR(KSP ksp, Vec h);
#endif
#endif

