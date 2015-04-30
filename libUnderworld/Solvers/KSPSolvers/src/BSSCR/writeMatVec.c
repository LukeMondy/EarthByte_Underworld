#include "writeMatVec.h"
#include "../../../../StgDomain/Utils/src/PETScCompatibility.h"

void bsscr_writeMat(Mat K, char name[], char message[]){
    PetscViewer             mat_view_file;
    char str[100];
    PetscTruth dump=PETSC_FALSE;
    PetscOptionsGetTruth( PETSC_NULL ,"-dump_matvec", &dump, 0 );
    if(dump){
    if( K != NULL ) {
	PetscPrintf( PETSC_COMM_WORLD,"%s \n",message);
	PetscObjectSetName((PetscObject)K,name);
	sprintf(str,"%smatrixBin",name);
	PetscViewerBinaryOpen( PETSC_COMM_WORLD, str, FILE_MODE_WRITE, &mat_view_file );
	PetscViewerSetFormat( mat_view_file, PETSC_VIEWER_NATIVE );
	MatView( K, mat_view_file );
	Stg_PetscViewerDestroy(&mat_view_file );
	sprintf(str,"%smatrix.m",name);
	PetscViewerASCIIOpen( PETSC_COMM_WORLD, str, &mat_view_file );
	PetscViewerSetFormat( mat_view_file, PETSC_VIEWER_ASCII_MATLAB );
	MatView( K, mat_view_file );
	Stg_PetscViewerDestroy(&mat_view_file );
    }
    }
}
void bsscr_writeVec(Vec V, char name[], char message[]){
    PetscViewer             vec_view_file;
    char str[100];
    PetscTruth dump=PETSC_FALSE;
    PetscOptionsGetTruth( PETSC_NULL ,"-dump_matvec", &dump, 0 );
    if(dump){
    if( V != NULL ) {
	PetscPrintf( PETSC_COMM_WORLD,"%s \n",message);
	PetscObjectSetName((PetscObject)V,name);
	sprintf(str,"%svectorBin",name);
	PetscViewerBinaryOpen( PETSC_COMM_WORLD, str, FILE_MODE_WRITE, &vec_view_file );
	PetscViewerSetFormat( vec_view_file, PETSC_VIEWER_NATIVE );
	VecView( V, vec_view_file );
	Stg_PetscViewerDestroy(&vec_view_file );
	sprintf(str,"%svector.m",name);
	PetscViewerASCIIOpen( PETSC_COMM_WORLD, str, &vec_view_file );
	PetscViewerSetFormat( vec_view_file, PETSC_VIEWER_ASCII_MATLAB );
	VecView( V, vec_view_file );
	Stg_PetscViewerDestroy(&vec_view_file );
    }
    }
}
void bsscr_dirwriteMat(Mat K, char name[], char dir[], char message[]){
    PetscViewer             mat_view_file;
    char str[100];
    PetscTruth dump=PETSC_FALSE;
    PetscOptionsGetTruth( PETSC_NULL ,"-dump_matvec", &dump, 0 );
    if(dump){
    if( K != NULL ) {
	PetscPrintf( PETSC_COMM_WORLD,"%s \n",message);
	PetscObjectSetName((PetscObject)K,name);
	sprintf(str,"%s%smatrixBin",dir,name);
	PetscViewerBinaryOpen( PETSC_COMM_WORLD, str, FILE_MODE_WRITE, &mat_view_file );
	PetscViewerSetFormat( mat_view_file, PETSC_VIEWER_NATIVE );
	MatView( K, mat_view_file );
	Stg_PetscViewerDestroy(&mat_view_file );
	sprintf(str,"%s%smatrix.m",dir,name);
	PetscViewerASCIIOpen( PETSC_COMM_WORLD, str, &mat_view_file );
	PetscViewerSetFormat( mat_view_file, PETSC_VIEWER_ASCII_MATLAB );
	MatView( K, mat_view_file );
	Stg_PetscViewerDestroy(&mat_view_file );
    }
    }
}
void bsscr_dirwriteVec(Vec V, char name[], char dir[], char message[]){
    PetscViewer             vec_view_file;
    char str[100];
    PetscTruth dump=PETSC_FALSE;
    PetscOptionsGetTruth( PETSC_NULL ,"-dump_matvec", &dump, 0 );
    if(dump){
    if( V != NULL ) {
	PetscPrintf( PETSC_COMM_WORLD,"%s \n",message);
	PetscObjectSetName((PetscObject)V,name);
	sprintf(str,"%s%svectorBin",dir,name);
	PetscViewerBinaryOpen( PETSC_COMM_WORLD, str, FILE_MODE_WRITE, &vec_view_file );
	PetscViewerSetFormat( vec_view_file, PETSC_VIEWER_NATIVE );
	VecView( V, vec_view_file );
	Stg_PetscViewerDestroy(&vec_view_file );
	sprintf(str,"%s%svector.m",dir,name);
	PetscViewerASCIIOpen( PETSC_COMM_WORLD, str, &vec_view_file );
	PetscViewerSetFormat( vec_view_file, PETSC_VIEWER_ASCII_MATLAB );
	VecView( V, vec_view_file );
	Stg_PetscViewerDestroy(&vec_view_file );
    }
    }
}
