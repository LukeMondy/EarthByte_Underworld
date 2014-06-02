/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "types.h"
#include "VoxelDataHandler_Abstract.h"
#include "VoxelDataHandler_GMT.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>

const Type VoxelDataHandler_GMT_Type = "VoxelDataHandler_GMT";

VoxelDataHandler_GMT* _VoxelDataHandler_GMT_New( VOXELDATAHANDLER_GMT_DEFARGS )
{
   VoxelDataHandler_GMT* self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof( VoxelDataHandler_GMT ) );
   self = (VoxelDataHandler_GMT*)_VoxelDataHandler_Abstract_New(  VOXELDATAHANDLER_ABSTRACT_PASSARGS  );

   self->dataType         = Variable_DataType_Double;
   strcpy( self->delim, " ,\t\r\n>");

   return self;
}

void _VoxelDataHandler_GMT_Delete( void* voxelDataHandler ) {
   VoxelDataHandler_GMT* self = (VoxelDataHandler_GMT*)voxelDataHandler;

   /* delete parent class */
   _VoxelDataHandler_Abstract_Delete( self );
}

void _VoxelDataHandler_GMT_Print( void* voxelDataHandler, Stream* stream ) {
   VoxelDataHandler_GMT* self = (VoxelDataHandler_GMT*)voxelDataHandler;

   /* General info */
   Journal_Printf( stream, "VoxelDataHandler_GMT (ptr): %p:\n", self );
   Stream_Indent( stream );

   /* Parent class info */
   _VoxelDataHandler_Abstract_Print( self, stream );

   /* VoxelDataHandler_GMT */
   Journal_Printf( stream, "filename: %s\n", self->filename );

   Stream_UnIndent( stream );
}

void* _VoxelDataHandler_GMT_DefaultNew( Name name ) {
   /* Variables set in this function */
   SizeT                                                                _sizeOfSelf = sizeof(VoxelDataHandler_GMT);
   Type                                                                        type = VoxelDataHandler_GMT_Type;
   Stg_Class_DeleteFunction*                                                _delete = _VoxelDataHandler_GMT_Delete;
   Stg_Class_PrintFunction*                                                  _print = _VoxelDataHandler_GMT_Print;
   Stg_Class_CopyFunction*                                                    _copy = _Stg_Component_Copy;
   AllocationType                                                nameAllocationType = NON_GLOBAL;
   Stg_Component_DefaultConstructorFunction*                    _defaultConstructor = _VoxelDataHandler_GMT_DefaultNew;
   Stg_Component_ConstructFunction*                                      _construct = _VoxelDataHandler_GMT_AssignFromXML;
   Stg_Component_BuildFunction*                                              _build = _VoxelDataHandler_GMT_Build;
   Stg_Component_InitialiseFunction*                                    _initialise = _VoxelDataHandler_GMT_Initialise;
   Stg_Component_ExecuteFunction*                                          _execute = _VoxelDataHandler_GMT_Execute;
   Stg_Component_DestroyFunction*                                          _destroy = _VoxelDataHandler_GMT_Destroy;
   VoxelDataHandler_GetTotalVoxelCount*                         _getTotalVoxelCount = _VoxelDataHandler_Abstract_GetTotalVoxelCount;
   VoxelDataHandler_GotoVoxelDataTop*                             _gotoVoxelDataTop = _VoxelDataHandler_GMT_GotoVoxelDataTop;
   VoxelDataHandler_IncrementVoxelIndex*                       _incrementVoxelIndex = _VoxelDataHandler_Abstract_IncrementVoxelIndex;
   VoxelDataHandler_GetCurrentVoxelIndex*                     _getCurrentVoxelIndex = _VoxelDataHandler_Abstract_GetCurrentVoxelIndex;
   VoxelDataHandler_GetCurrentVoxelCoord*                     _getCurrentVoxelCoord = _VoxelDataHandler_Abstract_GetCurrentVoxelCoord;
   VoxelDataHandler_GetCurrentVoxelDataFunc*                   _getCurrentVoxelData = _VoxelDataHandler_GMT_GetCurrentVoxelData;

   return (void*)_VoxelDataHandler_GMT_New( VOXELDATAHANDLER_GMT_PASSARGS  );
}
/*****************************************************************************************************
   _VoxelDataHandler_GMT_InitialFileScan assumes the file being read in has latitude/longitude values 
   in first two columns of data and that this data is evenly spaced and aligned along lat/lon grid
   lines. So the step sizes are always the same throughout the file. 
  
   The fastest moving coordinate is assumed to be Longitude in the range [-180:180]

   The indices of the coordinates are 0 1 2 from fastest to slowest.
   The functions here are AGNOSTIC about mappings between coordinates.
 *****************************************************************************************************/
void  _VoxelDataHandler_GMT_InitialFileScan( void* voxelDataHandler ){
    VoxelDataHandler_GMT*    self     = (VoxelDataHandler_GMT*) voxelDataHandler;
    double col1Prev,col1Curr, col2Prev,col2Curr, col1Step, col2Step;
    int    lineCounter        =  0;
    int    fastCoordFound     =  0;
    int    fastCoordFlag      =  0;
    int    fastCoord          = -1;
    int    slowCoord          = -1;
    int    bogus              =  0;
    //int    slowCrdCellsAcross =  0;
    int    fastCrdCellsAcross =  1;
    int    scanDone           =  0;
    double stepSizeFast       =  0.0;
    double stepSizeSlow       =  0.0;
    double firstFastCoord     = -DBL_MAX;
    double firstSlowCoord     = -DBL_MAX;
    double firstA, firstB;
    double temp;
    char* result;
    fpos_t strPos, strPosStart;
    int line;
    double colA[2], colB[2];
    /** rewind file and point to data */
    _VoxelDataHandler_GMT_GotoVoxelDataTop( self );
    fgetpos(CFile_Ptr( self->fileData ), &strPosStart);
    /** keep grabbing new lines until EOF or we have detected the fast moving coordinate, this should be quick and only scan a small part of the file */
    /** scan to find the fastCoord twice to handle case of repeated data on latitude +-90 from gplates */
    /** get first two lines */

//    while(!fastCoordFound){
    /** If data is good should find fast coord in first two lines */
    for(line=0;line<2;line++){
        //fgetpos(CFile_Ptr( self->fileData ), &strPos);
        result = fgets( (char *) &self->currentString, VOXELMAX_LINE_LENGTH_DEFINE, CFile_Ptr( self->fileData ) );
        if ( result ){
            self->currentToken = strtok_r( self->currentString, self->delim, &self->currentStrtokKey );
            sscanf( self->currentToken, "%lf", &col1Curr );
            self->currentToken = strtok_r( NULL, self->delim, &self->currentStrtokKey );
            sscanf( self->currentToken, "%lf", &col2Curr );
            colA[line]=col1Curr;
            colB[line]=col2Curr;
        }else{ Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s' - No result for file scan Something wrong with file.",__func__,self->type, self->name ); }       
    }
    //fsetpos(CFile_Ptr( self->fileData ), &strPos); /** step back one step */
    //lineCounter++;
    /** Handle case where we cross the +-180 line while trying to detect the fast step size */
    if(colA[1]*colA[0] < 0){ colA[1]+=360.0; }
    if(colB[1]*colB[0] < 0){ colB[1]+=360.0; }
    /***************************/
    col1Step=(colA[1]-colA[0]);
    col2Step=(colB[1]-colB[0]);

    if(fabs(col1Step) > 10*DBL_EPSILON){
        fastCoord      =0;
        slowCoord      =1;
        fastCoordFound++;
        fastCoordFlag++;
        stepSizeFast = col1Step;
        firstFastCoord = colA[0];
        firstSlowCoord = colB[0];
    }
    if(fabs(col2Step) > 10*DBL_EPSILON){
        fastCoord      =1;
        slowCoord      =0;
        fastCoordFound++;
        fastCoordFlag++;
        stepSizeFast = col2Step;
        firstFastCoord = colB[0];
        firstSlowCoord = colA[0];
    }
        
    if(!fastCoordFound){/** Then we have bogus data: All the longitudes (fast coord) are 0.0 for given latitude at start of array instead of incrementing like they should */
        bogus = 1;
    }
    /** if NOT bogus then we have found the fast step size but not the slow one yet */
    _VoxelDataHandler_GMT_GotoVoxelDataTop( self );
    if(bogus){
        result = fgets( (char *) &self->currentString, VOXELMAX_LINE_LENGTH_DEFINE, CFile_Ptr( self->fileData ) );
        if ( result ){
            self->currentToken = strtok_r( self->currentString, self->delim, &self->currentStrtokKey );
            sscanf( self->currentToken, "%lf", &col1Curr );
            self->currentToken = strtok_r( NULL, self->delim, &self->currentStrtokKey );
            sscanf( self->currentToken, "%lf", &col2Curr );
            colA[0]=col1Curr;
            colB[0]=col2Curr;
        }else{ Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s' - No result for file scan Something wrong with file.",__func__,self->type, self->name ); }
        col1Prev=col1Curr;
        col2Prev=col2Curr;      
    }
    fsetpos(CFile_Ptr( self->fileData ), &strPosStart); /** back to start of data. to be able to count to get cellnum for slow data? */
    while(bogus){ /** Skip ahead until one of the columns of data changes */            
        for(line=0;line<2;line++){
            fgetpos(CFile_Ptr( self->fileData ), &strPos);
            result = fgets( (char *) &self->currentString, VOXELMAX_LINE_LENGTH_DEFINE, CFile_Ptr( self->fileData ) );
            if ( result ){
                self->currentToken = strtok_r( self->currentString, self->delim, &self->currentStrtokKey );
                sscanf( self->currentToken, "%lf", &col1Curr );
                self->currentToken = strtok_r( NULL, self->delim, &self->currentStrtokKey );
                sscanf( self->currentToken, "%lf", &col2Curr );
                colA[line]=col1Curr;
                colB[line]=col2Curr;
            }else{ Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s' - No result for file scan Something wrong with file.",__func__,self->type, self->name ); }       
        }
        if(lineCounter == 0){
            firstA = colA[0];
            firstB = colB[0];
        }
        fsetpos(CFile_Ptr( self->fileData ), &strPos); /** step back one step */
        lineCounter++;
        /** Handle case where we cross the +-180 line while trying to detect the fast step size */
        if(colA[1]*colA[0] < 0){ colA[1]+=360.0; }
        if(colB[1]*colB[0] < 0){ colB[1]+=360.0; }
        /***************************/
        col1Step=(colA[1]-colA[0]);
        col2Step=(colB[1]-colB[0]);
        if(fabs(col1Step) > 10*DBL_EPSILON){
            bogus = 0; col1Step=0.0;
        }
        if(fabs(col2Step) > 10*DBL_EPSILON){
            bogus = 0; col2Step=0.0;
        }

        /** so now we are at start of good data */
        if(!bogus){
            for(line=0;line<2;line++){
                fgetpos(CFile_Ptr( self->fileData ), &strPos);
                result = fgets( (char *) &self->currentString, VOXELMAX_LINE_LENGTH_DEFINE, CFile_Ptr( self->fileData ) );
                if ( result ){
                    self->currentToken = strtok_r( self->currentString, self->delim, &self->currentStrtokKey );
                    sscanf( self->currentToken, "%lf", &col1Curr );
                    self->currentToken = strtok_r( NULL, self->delim, &self->currentStrtokKey );
                    sscanf( self->currentToken, "%lf", &col2Curr );
                    colA[line]=col1Curr;
                    colB[line]=col2Curr;
                }else{ Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s' - No result for file scan Something wrong with file.",__func__,self->type, self->name ); }       
            }
            fsetpos(CFile_Ptr( self->fileData ), &strPos); /** step back one step */
            lineCounter++;
            /** Handle case where we cross the +-180 line while trying to detect the fast step size */
            if(colA[1]*colA[0] < 0){ colA[1]+=360.0; }
            if(colB[1]*colB[0] < 0){ colB[1]+=360.0; }
            /***************************/
            col1Step=(colA[1]-colA[0]);
            col2Step=(colB[1]-colB[0]);
            if(fabs(col1Step) > 10*DBL_EPSILON){
                fastCoord      =0;
                slowCoord      =1;
                fastCoordFound++;
                fastCoordFlag++;
                stepSizeFast = col1Step;
                firstFastCoord = firstA;
                firstSlowCoord = firstB;
            }
            if(fabs(col2Step) > 10*DBL_EPSILON){
                fastCoord      =1;
                slowCoord      =0;
                fastCoordFound++;
                fastCoordFlag++;
                stepSizeFast = col2Step;
                firstFastCoord = firstB;
                firstSlowCoord = firstA;
            }
        }
    }//while(bogus)

    /** We should have the fast coord identified and have the fast step size when we are here */

    /** Now need start coords and total number of points to work out grid layout */
    /** Now start again at beginning of file to get the number of cells in the slow direction (should also be quick: as above)*/
    //double col1Min,col1Max, col2Min,col2Max;
    /** rewind file and point to data */
    lineCounter = 0;
    scanDone    = 0;
    _VoxelDataHandler_GMT_GotoVoxelDataTop( self );
    do{
        result = fgets( (char *)&self->currentString, VOXELMAX_LINE_LENGTH_DEFINE, CFile_Ptr( self->fileData ) );
        if ( result && (result[0]!='>')){
            if(lineCounter > 0){
                col1Prev=col1Curr;
                col2Prev=col2Curr;
            }
            self->currentToken = strtok_r( self->currentString, self->delim, &self->currentStrtokKey );
            sscanf( self->currentToken, "%lf", &col1Curr );
            self->currentToken = strtok_r( NULL, self->delim, &self->currentStrtokKey );
            sscanf( self->currentToken, "%lf", &col2Curr );
            if(lineCounter > 0){
                /** Need to handle case where we cross the +-180 line while trying to detect the fast step size */
                if(col1Curr*col1Prev < 0){ col1Curr+=360.0;  }
                if(col2Curr*col2Prev < 0){ col2Curr+=360.0;  }
                col1Step=(col1Curr-col1Prev);
                col2Step=(col2Curr-col2Prev);
                if(fastCoord == 1){/** then col1Step is Latitude and will jump to next level after row of Longitudes processed */
                    if(fabs(col1Step) < 10*DBL_EPSILON){
                        fastCrdCellsAcross++;           }
                    else{
                        scanDone=1;
                        stepSizeSlow=col1Step;    }
                }
                else{
                    if(fabs(col2Step) < 10*DBL_EPSILON){
                        fastCrdCellsAcross++;           }
                    else{
                        scanDone=1;
                        stepSizeSlow=col2Step;    } 
                }
            }
            lineCounter++;
        }/** if ( result ) */
    }while(!scanDone && result != NULL);/** or EOF */
    
    self->fastCoord = fastCoord;
    self->voxelSize[fastCoord] = stepSizeFast;
    //if(stepSizeSlow < 0.0){ stepSizeSlow = -stepSizeSlow; }
    self->voxelSize[slowCoord] = stepSizeSlow;
    self->voxelSize[2] = 1000.0; /** arbitrary value for K coord */
    self->startCrd[fastCoord] = firstFastCoord;
    self->startCrd[slowCoord] = firstSlowCoord;
    self->startCrd[2] = 0.0;
    self->numCells[fastCoord] = fastCrdCellsAcross;
    self->numCells[2] = 1;
    printf("fastCoord = %d slowCoord = %d\n",fastCoord, slowCoord);
    printf("firstFastCoord = %lf firstSlowCoord = %lf\n",firstFastCoord, firstSlowCoord);
    printf("stepSizeFast = %lf stepSizeSlow = %lf\n",stepSizeFast, stepSizeSlow);
    printf("fastCrdCellsAcross = %d\n",fastCrdCellsAcross);
    /** Now start again at beginning of file to get the number of cells in each direction */
    /** rewind file and point to data */
    lineCounter = 0;
    _VoxelDataHandler_GMT_GotoVoxelDataTop( self );
    do{
        result = fgets( (char *)&self->currentString, VOXELMAX_LINE_LENGTH_DEFINE, CFile_Ptr( self->fileData ) );
        if ( result && (result[0]!='>')){

            self->currentToken = strtok_r( self->currentString, self->delim, &self->currentStrtokKey );
            sscanf( self->currentToken, "%lf", &col1Curr );
            self->currentToken = strtok_r( NULL, self->delim, &self->currentStrtokKey );
            sscanf( self->currentToken, "%lf", &col2Curr );
            lineCounter++;
        }
    }while(result != NULL);
    self->numCells[slowCoord] = lineCounter/fastCrdCellsAcross;
    self->noDataValue = -DBL_MAX;
    printf("numCells[0] = %d numCells[1] = %d\n",self->numCells[0],self->numCells[1]);
    _VoxelDataHandler_GMT_GotoVoxelDataTop( self );
}
void _VoxelDataHandler_GMT_AssignFromXML( void* voxelDataHandler, Stg_ComponentFactory *cf, void* data ) {
   VoxelDataHandler_GMT*    self     = (VoxelDataHandler_GMT*) voxelDataHandler;
   char* dataType;

   /** config parent */
   _VoxelDataHandler_Abstract_AssignFromXML( voxelDataHandler, cf, data );
   /** deliminators seperating ascii row data */
   /** broken due to dictionary stripping or escaping characters from xml input.. disable for now */
   //sscanf( Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"fileDeliminator", " ,\t" ), "%[^\n]", self->delim );

   dataType     = Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"DataType", "double" );
   if(     !strcmp(dataType, "char"  )) self->dataType = Variable_DataType_Char;
   else if(!strcmp(dataType, "int"   )) self->dataType = Variable_DataType_Int;
   else if(!strcmp(dataType, "float" )) self->dataType = Variable_DataType_Float;
   else if(!strcmp(dataType, "double")) self->dataType = Variable_DataType_Double;
   else
      Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s' - Provided DataType %s not know or not supported.",__func__,self->type,self->name, dataType );
   /** stride for data reading.. so if "PosI1 PosJ1 PosK1 Data1 PosI2 PosJ2 PosK2 Data2...." stride is 4 */
   /** negative values indicate a stride with the allowance of newline as required */
   self->dataStride  = Stg_ComponentFactory_GetInt( cf, self->name, (Dictionary_Entry_Key)"DataStride", 1 );
   if( self->dataStride == 0 )
      Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s' - A DataStride of zero is not valid. Stride must be >0 or <0",__func__,self->type,self->name);
   /** Position of voxel data in datum. starts counting at 1... so 1,2,3,4..  */
   self->dataPos     = Stg_ComponentFactory_GetInt( cf, self->name, (Dictionary_Entry_Key)"DataPos", 1 );
   if(  self->dataStride < 0 ){
       if( (abs(self->dataStride)<self->dataPos) )
           Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s' - Provided DataStride (%i) cannot be less than DataPosition (%i).",__func__,self->type, self->name, self->dataStride, self->dataPos );
   }
   
   /** hard set so no coordinate switching atm */
   self->switchAxis[0] = 0;
   self->switchAxis[1] = 1;
   self->switchAxis[2] = 2;

   /** open file */
   self->fileData  = CFile_NewRead( self->filename );
   Journal_Firewall( self->fileData != NULL,
        self->errorStream,
        "Error in %s for %s '%s' - Cannot find or open file '%s'.",
        __func__,
        self->type,
        self->name,
        self->filename );

   _VoxelDataHandler_GMT_InitialFileScan( self );

   /** setup voxel data ready to read, incase gototop routine isn't called */
   _VoxelDataHandler_GMT_GotoVoxelDataTop( self );

   /** now go ahead and complete setup */
   _VoxelDataHandler_Abstract_CompleteSetup( voxelDataHandler );
}

void _VoxelDataHandler_GMT_Build( void* voxelDataHandler, void* data ) {
   VoxelDataHandler_GMT*    self     = (VoxelDataHandler_GMT*) voxelDataHandler;

   _VoxelDataHandler_Abstract_Build( voxelDataHandler, data);

   switch ( self->dataType )
   {
      case Variable_DataType_Char   : self->dataBufferSize = sizeof(signed char);  break;
      case Variable_DataType_Int    : self->dataBufferSize = sizeof(        int);  break;
      case Variable_DataType_Float  : self->dataBufferSize = sizeof(      float);  break;
      case Variable_DataType_Double : self->dataBufferSize = sizeof(     double);  break;
      default                       : Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s' - Datatype not know or not supported.",__func__,self->type,self->name ); break;
   }
   self->dataBuffer = Memory_Alloc_Bytes_Unnamed( self->dataBufferSize , "voxel data handler type" );
}

void _VoxelDataHandler_GMT_Initialise( void* voxelDataHandler, void* data ) { 

    _VoxelDataHandler_Abstract_Initialise( voxelDataHandler, data); 


}

void _VoxelDataHandler_GMT_Execute( void* voxelDataHandler, void* data )  { _VoxelDataHandler_Abstract_Execute( voxelDataHandler, data); }

void _VoxelDataHandler_GMT_Destroy( void* voxelDataHandler, void* data ) {
    //VoxelDataHandler_GMT*  self = (VoxelDataHandler_GMT*)voxelDataHandler;
   /** destroy parent */
   _VoxelDataHandler_Abstract_Destroy( voxelDataHandler, data );
}

int _VoxelDataHandler_GMT_GotoVoxelDataTop( void* voxelDataHandler ){
   VoxelDataHandler_GMT*  self = (VoxelDataHandler_GMT*)voxelDataHandler;

   /** if open, close, then reopen file */
   if( self->fileData) 
      File_Close( self->fileData );
   self->fileData  = CFile_NewRead( self->filename );

   rewind( CFile_Ptr( self->fileData ) );
   /** skip blank lines (if any) */
   _VoxelDataHandler_GMT_SkipBlank( self->fileData );

   /* now re-init cursor */
   _VoxelDataHandler_Abstract_GotoVoxelDataTop( voxelDataHandler );

   return 1;
}

int _VoxelDataHandler_GMT_GetCurrentVoxelData( void* voxelDataHandler ){
   VoxelDataHandler_GMT*  self = (VoxelDataHandler_GMT*)voxelDataHandler;
   double doubleData;

   _VoxelDataHandler_Abstract_GetASCIIVoxelData( voxelDataHandler );

   switch ( self->dataType )
   {
      case Variable_DataType_Char   :sscanf( self->currentToken,  "%c", (  char*)self->dataBuffer ); break;
      case Variable_DataType_Int    :sscanf( self->currentToken,  "%i", (   int*)self->dataBuffer ); break;
      case Variable_DataType_Float  :sscanf( self->currentToken,  "%f", ( float*)self->dataBuffer ); break;
      case Variable_DataType_Double :sscanf( self->currentToken, "%lf", (double*)self->dataBuffer ); break;
   }

   switch ( self->dataType )
   {
      case Variable_DataType_Char   : doubleData = (double) *(signed char*)self->dataBuffer; break;
      case Variable_DataType_Int    : doubleData = (double) *(        int*)self->dataBuffer; break;
      case Variable_DataType_Float  : doubleData = (double) *(      float*)self->dataBuffer; break;
      case Variable_DataType_Double : doubleData = (double) *(     double*)self->dataBuffer; break;
   }
   if(doubleData > self->noDataValue)
      return 1;
   else
      return 0;

}

void _VoxelDataHandler_GMT_SkipBlank( File* file ){
   int   testChar;
   Bool  blankLine = True;
   Bool  headerLine = True;
   char  *testStr;
   fpos_t strPos;
   int   ii;
   int N=1000;

   testStr=(char *)malloc(N*sizeof(char));
   /** first check if line is empty.  keep skipping until non-empty line is found */
   while( blankLine == True ){
      testChar = fgetc(CFile_Ptr( file ));
      if ( testChar != '\n' )
         blankLine = False;
   }
   /** we must now rewind a character */
   fseek( CFile_Ptr( file ), -1, SEEK_CUR );

   /** now skip potential header lines */
   while( headerLine == True ){
      fgetpos(CFile_Ptr( file ), &strPos);
      fgets( testStr, N, CFile_Ptr( file ));
      if(testStr[0] != '>') headerLine = False;
   }
   fsetpos(CFile_Ptr( file ), &strPos);

   free(testStr);
}

