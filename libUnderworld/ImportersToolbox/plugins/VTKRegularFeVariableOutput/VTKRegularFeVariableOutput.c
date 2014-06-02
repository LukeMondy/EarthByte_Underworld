/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ** Copyright (c) 2005-2010, Monash University
 ** All rights reserved.
 ** Redistribution and use in source and binary forms, with or without modification,
 ** are permitted provided that the following conditions are met:
 **
 **       * Redistributions of source code must retain the above copyright notice,
 **          this list of conditions and the following disclaimer.
 **       * Redistributions in binary form must reproduce the above copyright
 **         notice, this list of conditions and the following disclaimer in the
 **         documentation and/or other materials provided with the distribution.
 **       * Neither the name of the Monash University nor the names of its contributors
 **         may be used to endorse or promote products derived from this software
 **         without specific prior written permission.
 **
 ** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 ** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 ** THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 ** PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
 ** BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 ** CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 ** SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 ** HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 ** LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 ** OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 **
 **
 ** Contact:
 *%  Louis.Moresi - Louis.Moresi@monash.edu
 *%
 ** Contributors:
 **
 **~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/** \file 
\details
Generates a VTK Structured Point output from a regular fevariable
 **/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include <string.h>

#if defined( READ_HDF5) || defined( WRITE_HDF5)
#include <hdf5.h>
#endif


const Type VTKRegularFeVariableOutput_Type = "ImportersToolbox_VTKRegularFeVariableOutput";

void VTKRegularFeVariableOutput_OutputAll( void* _context ) ;

typedef struct {
    __Codelet
    Stg_ComponentFactory* cf;
} VTKRegularFeVariableOutput;

void _VTKRegularFeVariableOutput_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data ) {
    VTKRegularFeVariableOutput* self = (VTKRegularFeVariableOutput*) _self;

    self->context = (FiniteElementContext*)Stg_ComponentFactory_PluginConstructByKey( cf, self, (Dictionary_Entry_Key)"Context", FiniteElementContext, True, data );

    EntryPoint_Append_AlwaysLast(
            Context_GetEntryPoint( self->context, AbstractContext_EP_Save ),
            (Name)"VTKRegularFeVariableOutput_OutputAll",
            VTKRegularFeVariableOutput_OutputAll,
            VTKRegularFeVariableOutput_Type );
    EntryPoint_Append_AlwaysLast(
            Context_GetEntryPoint( self->context, AbstractContext_EP_DataSave ),
            (Name)"VTKRegularFeVariableOutput_OutputAll",
            VTKRegularFeVariableOutput_OutputAll,
            VTKRegularFeVariableOutput_Type );


}

void* _VTKRegularFeVariableOutput_DefaultNew( Name name ) {
    /* Variables set in this function */
    SizeT                                              _sizeOfSelf = sizeof(VTKRegularFeVariableOutput);
    Type                                                      type = VTKRegularFeVariableOutput_Type;
    Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
    Stg_Class_PrintFunction*                                _print = _Codelet_Print;
    Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
    Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _VTKRegularFeVariableOutput_DefaultNew;
    Stg_Component_ConstructFunction*                    _construct = _VTKRegularFeVariableOutput_AssignFromXML;
    Stg_Component_BuildFunction*                            _build = _Codelet_Build;
    Stg_Component_InitialiseFunction*                  _initialise = _Codelet_Initialise;
    Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
    Stg_Component_DestroyFunction*                        _destroy = _Codelet_Destroy;

    /*
     * Variables that are set to ZERO are variables that will be set either by the
     * current _New function or another parent _New function further up the hierachy.
     */
    AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;
#if defined( READ_HDF5) || defined( WRITE_HDF5)
#else
    Journal_Firewall(
            0,
            NULL,
            "Error:  VTKRegularFeVariableOutput plugin requires Underworld to be configured & compiled with HDF5.");
#endif

    return _Codelet_New(  CODELET_PASSARGS );
}

Index ImportersToolbox_VTKRegularFeVariableOutput_Register( PluginsManager* pluginsManager ) {
    return PluginsManager_Submit( pluginsManager, VTKRegularFeVariableOutput_Type, (Name)"0", _VTKRegularFeVariableOutput_DefaultNew );
}


void VTKRegularFeVariableOutput_OutputAll( void* _context ) {
#if defined( READ_HDF5) || defined( WRITE_HDF5)
    FiniteElementContext* context = (FiniteElementContext*) _context;
    FieldVariable*       fieldVar = NULL;
    FeVariable*          feVar    = NULL;
    FeMesh*              feMesh   = NULL;
    Index                var_I = 0;
    Stream*              errorStream = Journal_Register( Error_Type, (Name)"VTKRegularFeVariableOutput" );
    unsigned             componentCount = LiveComponentRegister_GetCount(stgLiveComponentRegister);
    unsigned             compI;
    Stg_Component*       stgComp;
    Bool                 fileOpened;
    Stream*              stream;
    char*                filename;
    char*                outputPathString;

    if( context->rank != 0 )
        return;

    /** search for entire live component register for feMesh types  **/
    for( compI = 0 ; compI < componentCount ; compI++ ){
        stgComp = LiveComponentRegister_At( stgLiveComponentRegister, compI );
        /* check that component is of type FeMesh, and that its element family is linear */
        if ( Stg_Class_IsInstance( stgComp, FeMesh_Type ) && (!strcmp( ((FeMesh*)stgComp)->feElFamily, "constant" ) || ((Mesh*)stgComp)->isCheckpointedAndReloaded ) ) {
            feMesh = (FeMesh*)stgComp;

            FeMesh* elFeMesh = feMesh;
            if ( !strcmp( ((FeMesh*)stgComp)->feElFamily, "constant" ) )
                elFeMesh =  ((C0Generator*)feMesh->generator)->elMesh;
            /* should we write this mesh? */
            if ( (elFeMesh->nElTypes!=1) || (!elFeMesh->isRegular) )
                continue;

            /* lets create the required file */
            stream = Journal_Register( InfoStream_Type, (Name)"VTKOutputFile"  );

            /** Set auto flush on stream **/
            Stream_SetAutoFlush( stream, True );

            outputPathString = Context_GetCheckPointWritePrefixString( context );

            /** Get name of VTK file **/
            Stg_asprintf( &filename, "%s%sFields.%05d.vtk", outputPathString, feMesh->name, context->timeStep );

            /** Init file, always overwriting any existing **/
            fileOpened = Stream_RedirectFile( stream, filename );
            Journal_Firewall( fileOpened, errorStream,
                    "Could not open file %s. Possibly directory %s does not exist or is not writable.\n"
                    "Check 'checkpointWritePath' in input file.", filename, outputPathString );

            /**----------------------- START GEOMETRY   ------------------------------------------------------------------------------------------------------------------- **/
            /* write header material */
            Journal_Printf( stream, "# vtk DataFile Version 2.0\n");
            Journal_Printf( stream, "Underworld VTKRegularFeVariableOutput: %s FeVariables, time=%g\n", feMesh->name, context->currentTime);
            Journal_Printf( stream, "ASCII\n");
            Journal_Printf( stream, "DATASET STRUCTURED_POINTS\n");

            /* write dimensions */
            /* check if we have a vert grid, and use if we do... else use the el grid */
            unsigned* sizes;
            Bool usingElGrid=False;
            if(feMesh->vertGridId != (unsigned)-1 )
                sizes = Grid_GetSizes( *(Grid**) Mesh_GetExtension( (Mesh*) feMesh, Grid*,  feMesh->vertGridId ) ); /** global no. of vertices in each dim */
            else
            {
                sizes = Grid_GetSizes( *(Grid**) Mesh_GetExtension( (Mesh*) feMesh, Grid*,  feMesh->elGridId ) ); /** global no. of elements in each dim */
                usingElGrid = True;
            }
            if( context->dim ==2 )
                Journal_Printf( stream, "DIMENSIONS %u %u 1\n", sizes[0], sizes[1]);
            else
                Journal_Printf( stream, "DIMENSIONS %u %u %u\n" , sizes[0], sizes[1], sizes[2]);

            /* write origin/spacing */
            double min[3], max[3];
            Mesh_GetGlobalCoordRange( feMesh, &min, &max );

            /* calculate origin & spacing */
            double origin[3];
            double spacing[3];

            /* if using vertgrid, reduce count by 1 */
            unsigned factor = 0;
            if( !usingElGrid )
                factor = 1;
            spacing[0] = (max[0]-min[0])/(sizes[0]-factor);
            spacing[1] = (max[1]-min[1])/(sizes[1]-factor);
            if( context->dim > 2 )
                spacing[2] = (max[2]-min[2])/(sizes[2]-factor);

            /* if on constant mesh (ie, usingElGrid), set origin to centre of element */
            double dfactor = 0;
            if( usingElGrid )
                dfactor = 1.;
            origin[0] = min[0] + dfactor*0.5*spacing[0];
            origin[1] = min[1] + dfactor*0.5*spacing[1];
            if( context->dim > 2 )
                origin[2] = min[2] + dfactor*0.5*spacing[2];



            if( context->dim == 2 )
                Journal_Printf( stream, "ORIGIN %g %g 0 \n", origin[0], origin[1]);
            else
                Journal_Printf( stream, "ORIGIN %g %g %g\n", origin[0], origin[1], origin[2]);

            if( context->dim ==2 )
                Journal_Printf( stream, "SPACING %g %g 1  \n", spacing[0], spacing[1]);
            else
                Journal_Printf( stream, "SPACING %g %g %g \n", spacing[0], spacing[1], spacing[2]);

            if( context->dim ==2 )
                Journal_Printf( stream, "POINT_DATA %u \n", sizes[0]*sizes[1]);
            else
                Journal_Printf( stream, "POINT_DATA %u \n", sizes[0]*sizes[1]*sizes[2]);
            /**----------------------- FINISH GEOMETRY  ------------------------------------------------------------------------------------------------------------------- **/

            /**----------------------- FIELDS  ---------------------------------------------------------------------------------------------------------------------------- **/
            for ( var_I = 0; var_I < context->fieldVariable_Register->objects->count; var_I++ ) {
                fieldVar = FieldVariable_Register_GetByIndex( context->fieldVariable_Register, var_I );

                if ( !Stg_Class_IsInstance( fieldVar, FeVariable_Type ) )
                    continue;
                feVar = (FeVariable*)fieldVar;
                if ( (feVar->isCheckpointedAndReloaded && (context->checkpointEvery != 0) && (context->timeStep % context->checkpointEvery == 0))  ||
                        (feVar->isCheckpointedAndReloaded && (context->checkpointAtTimeInc && (context->currentTime >= context->nextCheckpointTime))) ||
                        (feVar->isSavedData               && (context->saveDataEvery != 0) && (context->timeStep % context->saveDataEvery   == 0)) ){

                    /** make sure that the fevariable femesh is the same as that used above for the geometry definition, if so proceed **/
                    if( feVar->feMesh != feMesh )
                        continue;

                    hid_t                   file, fileSpace, fileData, error;
                    unsigned                totalNodes, ii, noffset;
                    hid_t                   memSpace;
                    hsize_t                 start[2], count[2], size[2], maxSize[2];
                    double*                 buf;

                    char* filenameFe;
                    Stg_asprintf( &filenameFe, "%s%s.%.5u.h5", outputPathString, feVar->name, context->timeStep );

                    /* Open the file and data set. */
                    file = H5Fopen( filenameFe, H5F_ACC_RDONLY, H5P_DEFAULT );

                    Journal_Firewall(
                            file >= 0,
                            errorStream,
                            "Error in %s for '%s' - Cannot open file %s.",
                            __func__,
                            "VTKRegularFeVariableOutput",
                            filenameFe );
#if(H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8) || H5Dopen_vers == 1
                    fileData = H5Dopen( file, "/data" );
#else
                    fileData = H5Dopen( file, "/data", H5P_DEFAULT );
#endif
                    fileSpace = H5Dget_space( fileData );
                    /* Get size of dataspace to determine if coords are in file */
                    H5Sget_simple_extent_dims( fileSpace, size, maxSize );

                    Bool saveCoords;
                    if( maxSize[1] > feVar->fieldComponentCount  )
                        saveCoords = True;
                    else
                        saveCoords = False;

                    start[1] = 0;
                    count[0] = 1;
                    count[1] = feVar->fieldComponentCount;

                    if( saveCoords )
                        count[1] += feVar->dim;

                    memSpace = H5Screate_simple( 2, count, NULL );
                    totalNodes = Mesh_GetGlobalSize( feVar->feMesh, 0 );

                    Journal_Firewall(
                            (maxSize[0] == totalNodes),
                            errorStream,
                            "Error in %s %s \n"
                            "Number of node values(%u) stored in %s does not correspond to total number of requested mesh nodes(%u).",
                            __func__,
                            "VTKRegularFeVariableOutput",
                            (unsigned int)maxSize[0],
                            filenameFe,
                            totalNodes );

                    /* print field header material */
                    if(feVar->fieldComponentCount == 1)
                    {
                        Journal_Printf( stream, "SCALARS %s double 1\n", feVar->name);
                        Journal_Printf( stream, "LOOKUP_TABLE default\n");
                    }
                    else if (feVar->fieldComponentCount <= 3 )
                        Journal_Printf( stream, "VECTORS %s double\n", feVar->name);
                    else
                        continue;
#if 0
                   Journal_Firewall(
                                0,
                                errorStream,
                                "Error in %s %s \n"
                                "FeVariable '%s' field component count is %d, while dimension is %d.  This combination is not supported.",
                                __func__,
                                "VTKRegularFeVariableOutput",
                                feVar->name,
                                feVar->fieldComponentCount, feVar->dim);
#endif

                    buf = Memory_Alloc_Array( double, count[1], "fileBuffer" );
                    /* Read data from HDF5 checkpoint file */
                    for( ii=0; ii<totalNodes; ii++ ) {
                        start[0] = ii;

                        H5Sselect_hyperslab( fileSpace, H5S_SELECT_SET, start, NULL, count, NULL );
                        H5Sselect_all( memSpace );

                        error = H5Dread( fileData, H5T_NATIVE_DOUBLE, memSpace, fileSpace, H5P_DEFAULT, buf );

                        Journal_Firewall(
                                error >= 0,
                                errorStream,
                                "Error in %s '%s' - Cannot read data in %s.",
                                __func__,
                                "VTKRegularFeVariableOutput",
                                filenameFe );

                        /* print data */
                        Index dof_I;
                        noffset = 0;
                        if( saveCoords )
                            noffset = feVar->dim;
                        for( dof_I = 0; dof_I < feVar->fieldComponentCount; dof_I++ ) {
                            Journal_Printf( stream, "%g ", buf[dof_I + noffset]);
                        }
                        if( feVar->fieldComponentCount == 2 )
                            Journal_Printf( stream, "0.");
                        Journal_Printf( stream, "\n");
                    }
                    Memory_Free( buf );

                    /* Close all handles */
                    H5Dclose( fileData );
                    H5Sclose( memSpace );
                    H5Sclose( fileSpace );
                    H5Fclose( file );

                }

            }
            Stream_CloseAndFreeFile(stream);
        }

    }

#endif
}   
