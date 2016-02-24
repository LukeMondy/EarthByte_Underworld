/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**      Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
**      Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**      AuScope - http://www.auscope.org
**      Monash University, AuScope SAM VIC - http://www.auscope.monash.edu.au
**      Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
*% Development Team :
*%  http://www.underworldproject.org/aboutus.html
**      Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**      Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**      Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**      David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**      Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**      Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**      Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**      Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**      Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**      Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**      Kent Humphries, Software Engineer, VPAC. (kenth@vpac.org)
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

#include <stdio.h>
#include <string.h>
#include <float.h>
#include <mpi.h>
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <Experimental/Experimental.h>

#include "Experimental/SurfaceProcess/grid.h"
#include "Experimental/SurfaceProcess/SurfaceProcess.h"

#include "Tellus.h"


const Type Experimental_Tellus_Type = "Experimental_Tellus";
Experimental_Tellus *Experimental_Tellus_self;
extern Experimental_SurfaceProcess *experimental_surfaceProcess_self;

void Experimental_Tellus_Process( Experimental_SurfaceProcess* sp ) {
    Experimental_Tellus *self = Experimental_Tellus_self;
    int rank;
   char buf[1000];

    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    // Export to Tellus
    if( rank == 0 ){
      Experimental_Tellus_Export( self, sp );
      sprintf( buf, "%s", self->curInputName );
    }
    MPI_Bcast( &buf, 1000, MPI_CHAR, 0, MPI_COMM_WORLD );

    // Call Tellus with the new XML input file
#ifdef HAVE_TELLUS
    tellus_call2f(buf) ;
#else
    Journal_Firewall( NULL,
         Journal_MyStream( Error_Type, self ),
         "\n\nError in %s for %s - Underworld appears to be unaware of any Tellus installation.\nPlease reconfigure Underworld and ensure Tellus is found.\n\n\n",
         __func__,
         self->type );
#endif
    MPI_Barrier( MPI_COMM_WORLD );

    // Import to Tellus
    if( rank == 0 )
        Experimental_Tellus_Import( self, sp );
    MPI_Bcast( sp->heights, sp->num_crds, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( sp->crds,  2*sp->num_crds, MPI_DOUBLE, 0, MPI_COMM_WORLD );

}

// Write the new vertical displacements field and create a new Tellus XML file
void Experimental_Tellus_Export( Experimental_Tellus* self, Experimental_SurfaceProcess* sp ) {
    FILE *inf, *outf, *outf2;
    int idx;

    // Create the vertical displacement file which is an ASCII file strutured as follow.
   char buf[1000];
   char buf2[1000];
   if( sp->timestep == 0 )
   {
	   sprintf( buf, "%s/%s-00000.vdisp",sp->context->outputPath, self->dispName);
       outf = fopen( buf, "w" );
   }
   else
   {
	  sprintf( buf, "%s/%s-%05d.vdisp",sp->context->outputPath,self->dispName, sp->timestep);
      outf = fopen( buf, "w" );
    }

   fprintf( outf, "* %lf %lf\n",
        sp->begin_time*sp->timeScaling/31536000.0 + sp->timeShift,
        sp->end_time*sp->timeScaling/31536000.0 + sp->timeShift );
   for( idx = 0; idx < sp->num_crds; idx++ ) {
       fprintf( outf, "%d %lf\n",idx+1, (sp->base_heights[idx] - sp->old_base_heights[idx])*sp->spaceScaling );
   }
   fprintf( outf, "**\n" );

    fclose(outf);

    /* We need to create a new Tellus XML input file. This is quite dodgy here, it will be much nicer
     * if I had used libxml2 for parsing and replacing the initial XML. */
   char out_fn[1000];
   char nbuf[1000];
   char inputF[1000];
   char obuf[1000];

   sprintf( out_fn, "%s/%s-%05d.xml", sp->context->outputPath, self->inputName, sp->timestep );
   sprintf( inputF, "%s.xml", self->inputName );
   outf = fopen( out_fn, "w" );

   if( sp->timestep > 0 ) {
      // New directory
      sprintf( nbuf, "%s-%05d", self->outdirName, sp->timestep );
      sprintf( obuf, "%s-%05d.vdisp",self->dispName, sp->timestep);
      sprintf( obuf, "%s-%05d", self->outdirName, sp->timestep - 1 );
   }
   else {
      sprintf( nbuf, "%s-00000", self->outdirName );
      sprintf( obuf, "%s-00000.vdisp",self->dispName);
   }

   // Change start and end time + directory name
   inf = fopen( inputF, "r" );

   while( !feof(inf) ) {
       char line[1000];
       if(!fgets( line, 1000, inf ))
      break;
       assert( line[strlen( line ) - 1] == '\n' );
       if( strstr( line, "startTime" ) ){
         fprintf( outf, "<startTime>%lf</startTime>\n",
             sp->begin_time*sp->timeScaling/31536000.0 + sp->timeShift );
       }
       else if( strstr( line, "endTime" ) ){
         fprintf( outf, "<endTime>%lf</endTime>\n",
             sp->end_time*sp->timeScaling/31536000.0 + sp->timeShift );
       }
       else if( strstr( line, "OutputDirectory" ) ){
         fprintf( outf, "<OutputDirectory>%s</OutputDirectory>\n",
             nbuf );
       }
       else if( strstr( line, "</tellus:tellusInput>" ) ){
         // Do nothing
       }
       else {
         fputs( line, outf );
       }
    }

   // Vertical displacement file
   fprintf( outf, "<VerticalDisplacement>\n" );
   fprintf( outf, "\t<nbDispInterval>%d</nbDispInterval>\n", 1); //sp->timestep+1 );
   fprintf( outf, "\t<VDisplacementsFile>%s/%s-%05d.vdisp</VDisplacementsFile>\n",
      sp->context->outputPath, self->dispName,sp->timestep );
   fprintf( outf, "</VerticalDisplacement>\n" );

   // Restart
   if( sp->timestep > 0 ) {

      fprintf( outf, "\n<Restart>\n" );
      fprintf( outf, "\t<RestartFolder>%s</RestartFolder>\n",obuf );
      fprintf( outf, "\t<RestartIt>%d</RestartIt>\n",self->checkfreq * sp->timestep );
      fprintf( outf, "</Restart>\n" );

   }

   // Vertical displacement file
   fprintf( outf, "\n</tellus:tellusInput>\n" );

   strcpy( self->curInputName, out_fn );

   fclose( inf );
   fclose( outf );

}


// Import new elevations due to erosion/deposition by surface processes simulated with Tellus.
// When running Tellus with the <underworldCoupling> flag on, it produces a serie of ascii files at each stepchecking interval
// called topsurf_XXXX.udw and containing the node index and top coordinates (X,Y,Z)
void Experimental_Tellus_Import( Experimental_Tellus* self, Experimental_SurfaceProcess* sp ) {
    FILE *in_f;
    char line[1000];

    char in_fn[100];

    sprintf( in_fn, "%s-%05d/outputs/topsurf_%05d.udw",self->outdirName, sp->timestep, sp->timestep );

   in_f = fopen( in_fn, "r" );
   assert( in_f );

     /* Don't need to initialise the current heights to the current basement.
    memcpy( sp->heights, sp->base_heights, sp->num_crds*sizeof(double) ); */

    /* Read each line, checking if we've hit another display time. */
    while( !feof( in_f ) ) {
      if(!fgets( line, 1000, in_f ))
         break;
      assert( line[strlen( line ) - 1] == '\n' );

      /* Scan the line. */
      {
         int idx;
         double CTval[3];

         sscanf( line, "%d %lf %lf %lf",
               &idx, CTval + 0, CTval + 1, CTval + 2);

         /* Get the surface elevations. */
         sp->heights[idx-1] = CTval[2]/sp->spaceScaling;
         /* Reset the crds */
         sp->crds[2*(idx-1) + 0] = CTval[0]/sp->spaceScaling;
         sp->crds[2*(idx-1) + 1] = CTval[1]/sp->spaceScaling;
      }
    }
}


void _Experimental_Tellus_AssignFromXML( void* component,
                Stg_ComponentFactory* cf,
                void* data )
{
    Experimental_Tellus* self = (Experimental_Tellus*)component;

    self->context = Stg_ComponentFactory_PluginConstructByKey( cf, self, (Dictionary_Entry_Key)"Context", AbstractContext, True, data  );
    self->initialSurfaceFile = Stg_ComponentFactory_PluginGetString( cf, self, "InitialSurfaceFile", "surface.nodes" );
    self->inputName = Stg_ComponentFactory_PluginGetString( cf, self, "InputName", "input.xml" );
    self->depositName = Stg_ComponentFactory_PluginGetString( cf, self, "DepositName", "input.dep" ) ;
    self->dispName = Stg_ComponentFactory_PluginGetString( cf, self, "DispName", "input.disp" ) ;
    self->outdirName = Stg_ComponentFactory_PluginGetString( cf, self, "OutdirName", "output_tellus" );
    self->checkfreq = Stg_ComponentFactory_PluginGetInt( cf, self, "CheckFreq", 10 );

    Experimental_Tellus_self = self;
}

void _Experimental_Tellus_Build( void* component, void* data ) {
    EntryPoint_Append( experimental_surfaceProcess_self->processEP,
                       "Experimental_Tellus_Plugin",
                       Experimental_Tellus_Process,
                       Experimental_Tellus_Type );
}

// Read the initial surface elevation
void _Experimental_Tellus_Initialise( void* component, void* data ) {
    Experimental_Tellus* self = (Experimental_Tellus*)component;
    Experimental_SurfaceProcess* sp = experimental_surfaceProcess_self;
    FILE* f;
    int numVertex;
    int idx ;
    char line[1000];

    /* Get the number of coordinates from the initial surface file. */
    f = fopen( self->initialSurfaceFile, "r" );
    assert( f );

    // NOT SURE WE NEED TO HAVE A STRUCTURED SURFACE
    // AND IF numRows numCols should be defined for udw
    // SO I JUST READ THE VERTEX NUMBER FOR THE UNSTRUCTURED MESH
    fscanf( f, "%d \n", &numVertex );
    /* step forward to next line */
    fgets( line, 1000, f );

    Experimental_SurfaceProcess_SetNumCoords( sp, numVertex );

    /* Generate the coordinates and heights. */
    /* Read each line */
    while( !feof( f ) ) {
        if(!fgets( line, 1000, f ))
                break;
        //assert( line[strlen( line ) - 1] == '\n' );

        /* Scan the line. */
        {
                double CTval[3];

                sscanf( line, "%d %lf %lf %lf",
                                &idx, CTval + 0, CTval + 1, CTval + 2);

                /* Set the surface elevations. */
      if( (idx > 0) && (idx <= numVertex) ) {
                  sp->crds[2*(idx-1) + 0] = CTval[0]/sp->spaceScaling;
                  sp->crds[2*(idx-1) + 1] = CTval[1]/sp->spaceScaling;
                  sp->heights[idx-1]      = CTval[2]/sp->spaceScaling;
                }
        }
    }

    fclose( f );
}


void* _Experimental_Tellus_DefaultNew( Name name ) {
    /* Variables set in this function */
    SizeT                                              _sizeOfSelf = sizeof(Experimental_Tellus);
    Type                                                      type = Experimental_Tellus_Type;
    Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
    Stg_Class_PrintFunction*                                _print = _Codelet_Print;
    Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
    Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Experimental_Tellus_DefaultNew;
    Stg_Component_ConstructFunction*                    _construct = _Experimental_Tellus_AssignFromXML;
    Stg_Component_BuildFunction*                            _build = _Experimental_Tellus_Build;
    Stg_Component_InitialiseFunction*                  _initialise = _Experimental_Tellus_Initialise;
    Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
    Stg_Component_DestroyFunction*                        _destroy = _Codelet_Destroy;

    /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
    AllocationType  nameAllocationType = ZERO;

    return _Codelet_New(  CODELET_PASSARGS  );
}

Index Experimental_Tellus_Register( PluginsManager* pluginsManager ) {
    return PluginsManager_Submit( pluginsManager, Experimental_Tellus_Type, (Name)"0", _Experimental_Tellus_DefaultNew  );
}
