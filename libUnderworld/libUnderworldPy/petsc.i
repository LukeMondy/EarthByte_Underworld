/* -*- C -*-  (not really, but good for syntax highlighting) */

%module petsc

%import "petsc.i"

%{
/* Includes the header in the wrapper code */
/* #include <petscext.h> */
#include <petscsys.h>
#include <petscviewer.h>
%}


%inline %{

void OptionsInsertString(char * string){
    PetscOptionsInsertString(string);
}

void OptionsPrint(){
//     PetscOptionsView(PETSC_VIEWER_STDOUT_WORLD);
#if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 2 ) )
    PetscOptionsView(PETSC_STDOUT);
#else
    PetscOptionsPrint(PETSC_STDOUT);
#endif
}

void OptionsClear(){
    PetscOptionsClear();
}

void OptionsSetValue(const char iname[],const char value[]){
    PetscOptionsSetValue(iname, value);
}

void OptionsClearValue(const char iname[]){
    /* option name must being with a - */
    PetscOptionsClearValue(iname);
}

void OptionsInsertFile(const char file[]){
    PetscOptionsInsertFile(PETSC_COMM_WORLD,file,PETSC_TRUE);
}

void EmptyCall(){
  return 1;
}
%}
