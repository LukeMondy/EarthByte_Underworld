#ifdef HAVE_PETSCEXT
#ifndef _writeMat_h
#define _writeMat_h
#include <petsc.h>
#include <petscmat.h>
#include <petscvec.h>

void bsscr_writeMat(Mat A, char name[], char message[]);
void bsscr_writeVec(Vec V, char name[], char message[]);
void bsscr_dirwriteMat(Mat A, char name[], char dir[], char message[]);
void bsscr_dirwriteVec(Vec V, char name[], char dir[], char message[]);
#endif
#endif
