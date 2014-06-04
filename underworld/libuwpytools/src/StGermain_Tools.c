/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**   Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
**   Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**   AuScope - http://www.auscope.org
**   Monash University, AuScope SAM VIC - http://www.auscope.monash.edu.au
**   Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
** Contributors:
**   Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**   Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**   Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**   David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**   Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**   Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**   Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**   Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**   Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**   Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**   Kent Humphries, Software Engineer, VPAC. (kenth@vpac.org)
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Lesser General Public
**  License as published by the Free Software Foundation; either
**  version 2.1 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License along with this library; if not, write to the Free Software
**  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const Type StGermain_Type = "StGermain";

#include "StGermain_Tools.h"

static int hasBeenMPIInit = 0;

StgData* StgInit( int argc, char* argv[] ) {
   StgData* data = (StgData*) malloc(sizeof(StgData));
   *data = (StgData){.commworld = NULL, .rank=-1, .nProcs=-1, .dictionary=NULL, .sources=NULL, .cf=NULL, .argcCpy=NULL, .argvCpy=NULL, .ioHandler=NULL}; 

   //lets copy all this data for safety
   data->argcCpy = (int*) malloc(sizeof(int));
   *(data->argcCpy) = argc;
   data->argvCpy = (char***) malloc(sizeof(char**)); 
   *(data->argvCpy) = (char**) malloc((argc+1)*sizeof(char*));
   int ii;
   for(ii = 0; ii<argc; ii++){
      (*(data->argvCpy))[ii] = (char*)malloc(strlen(argv[ii])+1);
      strcpy((*(data->argvCpy))[ii],argv[ii]);
   }
   (*(data->argvCpy))[argc]=NULL;  //add sentinel

   Index                 i;
   Stg_ObjectList*       inputPaths = NULL;
   char*                 inputPath = NULL;
   /* Initialise PETSc, get world info */
   if(!hasBeenMPIInit){
	   MPI_Init( data->argcCpy, data->argvCpy );
	   hasBeenMPIInit = 1;
   }
   MPI_Comm_dup( MPI_COMM_WORLD, &data->commworld );
   MPI_Comm_size( data->commworld, &data->nProcs );
   MPI_Comm_rank( data->commworld, &data->rank );
   StGermain_Init( data->argcCpy, data->argvCpy );

   /* Ensures copyright info always come first in output */
   MPI_Barrier( data->commworld ); 

   /* 
    * Parse the input path command line argument... needed before we start parsing the input.
    * And add the path to the global xml path dictionary. 
    */
   inputPaths = stgParseInputPathCmdLineArg( data->argcCpy, data->argvCpy );
   for( i = 0; i < Stg_ObjectList_Count( inputPaths ); i++ ) {
      inputPath = (char*)Stg_ObjectAdaptor_Object( (Stg_ObjectAdaptor*)Stg_ObjectList_At( inputPaths, i ) );
      XML_IO_Handler_AddDirectory( (Name)"--inputPath", inputPath );
      File_AddPath( inputPath );
   }
   Stg_Class_Delete( inputPaths );

   /* Create the application's dictionary & read input. */
   data->dictionary = Dictionary_New();
   data->sources = Dictionary_New();
   data->ioHandler = XML_IO_Handler_New();
   IO_Handler_ReadAllFromCommandLine( data->ioHandler, *(data->argcCpy), *(data->argvCpy), data->dictionary, data->sources );

   return data;
}

char* StgDictAsXMLString(StgData* data){
   return _XML_IO_Handler_WriteAllMem(XML_IO_Handler_New(), data->dictionary, data->sources );

} 

void StgSetDictFromXMLString(StgData* data, const char* xmlString, const char* tag){
   Stg_Class_Delete( data->dictionary );
   Stg_Class_Delete( data->sources );
   Stg_Class_Delete( data->ioHandler );
   data->dictionary = Dictionary_New();
   data->sources = Dictionary_New();
   data->ioHandler = XML_IO_Handler_New();
   StgAddToDictFromXMLString(data, xmlString, tag);
}

void StgAddToDictFromXMLString(StgData* data, const char* xmlString, const char* tag){
   data->ioHandler->currSources = data->sources;
   IO_Handler_ReadAllFromBuffer( data->ioHandler, xmlString, data->dictionary, tag );
}

int StgConstruct(StgData* data){
   Journal_ReadFromDictionary( data->dictionary );
   /* now dereference aliases */
   DictionaryUtils_AliasDereferenceDictionary( data->dictionary );
   
   ModulesManager_Load( stgToolboxesManager, data->dictionary, (Name)"" );

   data->cf = stgMainConstruct( data->dictionary, data->sources, data->commworld, NULL );
   return 0;
}

int StgBuildAndInitialise(StgData* data){
   stgMainBuildAndInitialise( data->cf );
   return 0;
}

int StgMainLoop(StgData* data){
   stgMainLoop( data->cf );
   return 0;
}

int StgRun(StgData* data){
   StgConstruct( data );
   StgBuildAndInitialise( data );
   StgMainLoop( data );
   return 0;
}

int StgFinalise(StgData* data){
   /* Close off everything */
   Stg_Class_Delete( data->sources );
   Stg_Class_Delete( data->ioHandler );
   stgMainDestroy( data->cf );
   StGermain_Finalise();
   //MPI_Finalize();

   /* free up these guys created earlier */
   int ii;
   for(ii = 0; ii<*(data->argcCpy); ii++)
   	  free((*(data->argvCpy))[ii]);
   free(data->argvCpy);
   free(data->argcCpy);
   free(data);   
   return 0; /* success */
}


