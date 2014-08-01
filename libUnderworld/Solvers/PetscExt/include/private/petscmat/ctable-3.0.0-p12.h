/* Contributed by - Mark Adams */

#if !defined(__CTABLE300_H)
#define __CTABLE300_H

struct _n_PetscTable {
  PetscInt *keytable;
  PetscInt *table;
  PetscInt count;
  PetscInt tablesize;
  PetscInt head;
};

#include "petscctable.h"
#endif
