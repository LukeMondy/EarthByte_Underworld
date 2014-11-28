/* Contributed by - Mark Adams */

#if !defined(__CTABLE31_H)
#define __CTABLE31_H

struct _n_PetscTable {
  PetscInt *keytable;
  PetscInt *table;
  PetscInt count;
  PetscInt tablesize;
  PetscInt head;
};

#include "petscctable.h"
#endif
