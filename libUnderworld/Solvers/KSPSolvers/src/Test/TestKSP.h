#ifdef HAVE_PETSCEXT
#ifndef __TestKSP_h__
#define __TestKSP_h__

typedef struct {
    Stokes_SLE* sle;
    Mat J, Jp, Gt;
} KSP_TEST;

PetscErrorCode PETSCKSP_DLLEXPORT KSPRegisterTEST(const char path[]);

#endif
#endif
