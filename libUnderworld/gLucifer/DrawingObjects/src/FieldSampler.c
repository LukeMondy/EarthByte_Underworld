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
*%  Owen Kaluza - Owen.Kaluza(at)monash.edu
*%
*% Development Team :
*%  http://www.underworldproject.org/aboutus.html
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include <gLucifer/Base/Base.h>


#include "types.h"
#include <gLucifer/Base/DrawingObject.h>
#include "FieldSampler.h"

void lucFieldSampler_SampleLocal( void* drawingObject, void* _context );
void lucFieldSampler_SampleGlobal( void* drawingObject, void* _context );

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucFieldSampler_Type = "lucFieldSampler";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucFieldSampler* _lucFieldSampler_New(  LUCFIELDSAMPLER_DEFARGS  )
{
   lucFieldSampler*               self;

   /* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
   assert( _sizeOfSelf >= sizeof(lucFieldSampler) );
   self = (lucFieldSampler*) _lucDrawingObject_New(  LUCDRAWINGOBJECT_PASSARGS  );

   return self;
}

void _lucFieldSampler_Init(
   lucFieldSampler*           self,
   IJK                        resolution,
   float                      pointSize,
   Bool                       pointSmoothing,
   Bool                       sampleGlobal)
{
   self->scaling = pointSize;
   memcpy( self->resolution, resolution, sizeof(IJK) );
   self->sampleGlobal = sampleGlobal;

   self->elementRes[I_AXIS] = Dictionary_GetInt( self->context->CF->rootDict, (Dictionary_Entry_Key)"elementResI"  );
   self->elementRes[J_AXIS] = Dictionary_GetInt( self->context->CF->rootDict, (Dictionary_Entry_Key)"elementResJ"  );
   self->elementRes[K_AXIS] = Dictionary_GetInt( self->context->CF->rootDict, (Dictionary_Entry_Key)"elementResK"  );

   /* Setup sampling resolution */
   if (self->resolution[I_AXIS] == 0) self->resolution[I_AXIS] = self->elementRes[I_AXIS];
   if (self->resolution[J_AXIS] == 0) self->resolution[J_AXIS] = self->elementRes[J_AXIS];
   if (self->resolution[K_AXIS] == 0) self->resolution[K_AXIS] = self->elementRes[K_AXIS];

   /* Append to property string */
   lucDrawingObject_AppendProps(self, "pointSmooth=%d\npointSize=%g\n", pointSmoothing, pointSize); 
}

void _lucFieldSampler_Delete( void* drawingObject )
{
   lucFieldSampler*  self = (lucFieldSampler*)drawingObject;

   /* Free memory */
   Memory_Free( self->samples );

   _lucDrawingObject_Delete( self );
}

void _lucFieldSampler_Print( void* drawingObject, Stream* stream )
{
   lucFieldSampler*  self = (lucFieldSampler*)drawingObject;

   _lucDrawingObject_Print( self, stream );
}

void* _lucFieldSampler_DefaultNew( Name name )
{
   /* Variables set in this function */
   SizeT                                                     _sizeOfSelf = sizeof(lucFieldSampler);
   Type                                                             type = lucFieldSampler_Type;
   Stg_Class_DeleteFunction*                                     _delete = _lucFieldSampler_Delete;
   Stg_Class_PrintFunction*                                       _print = _lucFieldSampler_Print;
   Stg_Class_CopyFunction*                                         _copy = NULL;
   Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _lucFieldSampler_DefaultNew;
   Stg_Component_ConstructFunction*                           _construct = _lucFieldSampler_AssignFromXML;
   Stg_Component_BuildFunction*                                   _build = _lucFieldSampler_Build;
   Stg_Component_InitialiseFunction*                         _initialise = _lucFieldSampler_Initialise;
   Stg_Component_ExecuteFunction*                               _execute = _lucFieldSampler_Execute;
   Stg_Component_DestroyFunction*                               _destroy = _lucFieldSampler_Destroy;
   lucDrawingObject_SetupFunction*                                _setup = _lucFieldSampler_Setup;
   lucDrawingObject_DrawFunction*                                  _draw = _lucFieldSampler_Draw;
   lucDrawingObject_CleanUpFunction*                            _cleanUp = lucDrawingObject_CleanUp;

   /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
   AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

   return (void*) _lucFieldSampler_New(  LUCFIELDSAMPLER_PASSARGS  );
}

void _lucFieldSampler_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data )
{
   lucFieldSampler*         self               = (lucFieldSampler*)drawingObject;
   Index                  defaultRes;
   IJK                    resolution;

   /* Construct Parent */
   _lucDrawingObject_AssignFromXML( self, cf, data );

   defaultRes = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"resolution", 0);
   resolution[ I_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"resolutionX", defaultRes);
   resolution[ J_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"resolutionY", defaultRes);
   resolution[ K_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"resolutionZ", defaultRes);

   _lucFieldSampler_Init(
      self,
      resolution,
      (float) Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"pointSize", 1.0  ),
      Stg_ComponentFactory_GetBool( cf, self->name, (Dictionary_Entry_Key)"pointSmoothing", True  ),
      Stg_ComponentFactory_GetBool( cf, self->name, (Dictionary_Entry_Key)"sampleGlobal", False  ));

   /* No lighting */
   self->lit = False;
}

void _lucFieldSampler_Build( void* drawingObject, void* data ) 
{
   lucFieldSampler*  self = (lucFieldSampler*)drawingObject;
   AbstractContext* context = self->context;
   Stg_ComponentFactory* cf = context->CF;

   /* It seems fields need to be constructed in build phase or are sometimes incomplete */
   self->sampleField = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"sampleField", FieldVariable, True, data  );
   Stg_Component_Build( self->sampleField, data, False );
}

void _lucFieldSampler_Initialise( void* drawingObject, void* data ) 
{
   lucFieldSampler*  self = (lucFieldSampler*)drawingObject;
   FieldVariable*    sampleField = self->sampleField;
   Coord             min, max, globalMin, globalMax;
   int i;

   if (!self->sampleGlobal)
   {
      /* Ensure resolution is a multiple of element res... */
      if (self->resolution[I_AXIS] % self->elementRes[I_AXIS] > 0)
         self->resolution[I_AXIS] += (self->elementRes[I_AXIS] - self->resolution[I_AXIS] % self->elementRes[I_AXIS]);
      if (self->resolution[J_AXIS] % self->elementRes[J_AXIS] > 0)
         self->resolution[J_AXIS] += (self->elementRes[J_AXIS] - self->resolution[J_AXIS] % self->elementRes[J_AXIS]);
      if (self->resolution[K_AXIS] % self->elementRes[K_AXIS] > 0)
         self->resolution[K_AXIS] += (self->elementRes[K_AXIS] - self->resolution[K_AXIS] % self->elementRes[K_AXIS]);
      /* Get resolution per element */
      self->nx = self->resolution[I_AXIS] / self->elementRes[I_AXIS];
      self->ny = self->resolution[J_AXIS] / self->elementRes[J_AXIS];
      self->nz = self->resolution[K_AXIS] / self->elementRes[K_AXIS];
      if (sampleField->dim == 2) self->nz = 0;
      /*printf("(%s) Resolution %d,%d,%d (nx/y/z %d,%d,%d)\n", self->name, self->resolution[0], self->resolution[1], self->resolution[2], self->nx, self->ny, self->nz);*/
   }
   if (sampleField->dim == 2) self->resolution[K_AXIS] = 0;

   FieldVariable_GetMinAndMaxGlobalCoords( sampleField, globalMin, globalMax );
   FieldVariable_GetMinAndMaxLocalCoords( sampleField, min, max );

   /* Calculate a resolution for this local domain based on a sample size over the entire domain */
   /*fprintf(stderr, "Old size %d,%d,%d ", self->resolution[0], self->resolution[1], self->resolution[2]);*/
   for (i=0; i<sampleField->dim; i++)
   {
      self->cell[i] = 0;
      if (self->resolution[i] > 0)
      {
         self->cell[i] = (globalMax[i] - globalMin[i])/((double) self->resolution[i]);
         /* force round up */
         self->resolution[i] = ceil((max[i] - min[i]) / self->cell[i]);
      }
   }
   /*fprintf(stderr, "New size %d,%d,%d ", self->resolution[0], self->resolution[1], self->resolution[2]);*/

   /* Calculate number of samples and allocate memory */
   self->count = (self->resolution[I_AXIS]+1) * (self->resolution[J_AXIS]+1) * (self->resolution[K_AXIS]+1);
   self->samples = Memory_Alloc_Array( Sample, self->count, (Name)"sample array" );
   Journal_Firewall( self->samples != NULL, lucError,
                     "Error - in %s(): Allocate memory for %d samples failed\n", __func__, self->count);
}

void _lucFieldSampler_Execute( void* drawingObject, void* data ) {}
void _lucFieldSampler_Destroy( void* drawingObject, void* data ) {}

/* Comparison for an ascending X,Y,Z coordinate sort */
int compareSample(const void *a, const void *b)
{
   Sample *s1 = (Sample*)a;
   Sample *s2 = (Sample*)b;
   
   if (s1->pos[2] == s2->pos[2])
   {
      if (s1->pos[1] == s2->pos[1])
      {
         if (s1->pos[0] == s2->pos[0])
            return 0;
         else
            return s1->pos[0] > s2->pos[0] ? 1 : -1;
      }
      else
         return s1->pos[1] > s2->pos[1] ? 1 : -1;
   }
   return s1->pos[2] > s2->pos[2] ? 1 : -1;
}

void _lucFieldSampler_Setup( void* drawingObject, lucDatabase* database, void* _context )
{
   lucFieldSampler*             self = (lucFieldSampler*)drawingObject;

   lucDrawingObject_SyncShadowValues( self, self->sampleField );

   if (self->sampleGlobal)
      lucFieldSampler_SampleGlobal(drawingObject, _context );
   else
      lucFieldSampler_SampleLocal(drawingObject, _context );

   /* Sort by X,Y,Z coordinate */
   qsort(self->samples, self->count, sizeof(Sample), compareSample);
}


/* Sample by element in local coords, faster, handles deformed meshes */
void lucFieldSampler_SampleLocal( void* drawingObject, void* _context )
{
   lucFieldSampler*             self               = (lucFieldSampler*)drawingObject;
   FeVariable*                feVariable         = (FeVariable*) self->sampleField;
   FeMesh*    		            mesh               = feVariable->feMesh;
   Element_LocalIndex         lElement_I;
   Element_LocalIndex         elementLocalCount  = FeMesh_GetElementLocalSize( mesh );
   int                        i, j, k;
   Coord                      globalMin, globalMax;

   FieldVariable_GetMinAndMaxGlobalCoords( feVariable, globalMin, globalMax );

   /*printf("Sampling %d elements, %d samples...\n", elementLocalCount, self->count);*/
   self->count = 0;
   for ( lElement_I = 0 ; lElement_I < elementLocalCount ; lElement_I++ )
   {
      for (i = 0 ; i <= self->nx; i++)
      {
         for (j = 0 ; j <= self->ny; j++)
         {
            for (k = 0 ; k <= self->nz; k++)
            {
               double value;
               /* Calc position within element in local coords */
               Coord local = {-1.0 + (self->nx > 0 ? 2.0 * i / self->nx : 0.0), 
                              -1.0 + (self->ny > 0 ? 2.0 * j / self->ny : 0.0),
                              -1.0 + (self->nz > 0 ? 2.0 * k / self->nz : 0.0)};

               /* Get global coordinates */
               double pos[3] = {0,0,0};
               FeMesh_CoordLocalToGlobal( mesh, lElement_I, local, pos );

               /* Skip outer edges of element except on end boundaries */
               if (i > 0 && i == self->nx && pos[0] < globalMax[I_AXIS]) continue;
               if (j > 0 && j == self->ny && pos[1] < globalMax[J_AXIS]) continue;
               if (k > 0 && k == self->nz && pos[2] < globalMax[K_AXIS]) continue;

               /* Get value at coords (faster using element and local coords) */
               FeVariable_InterpolateWithinElement( feVariable, lElement_I, local, &value);
               /*if (self->count % 10 == 0) 
                 printf("%d) %d,%d,%d  %f,%f,%f val = %f\n", self->count, i, j, k, pos[0], pos[1], pos[2], value);*/

               self->samples[self->count].value = (float)value;
               self->samples[self->count].pos[0] = (float)pos[0];
               self->samples[self->count].pos[1] = (float)pos[1];
               self->samples[self->count].pos[2] = (float)pos[2];
               self->count++;
            }
         }
      }
   }
}

/* Sampling in global coords, slower, assumes regular grid */
void lucFieldSampler_SampleGlobal( void* drawingObject, void* _context )
{
   lucFieldSampler* self = (lucFieldSampler*)drawingObject;
   FieldVariable*   sampleField = self->sampleField;
   int              i, j, k;
   Coord            pos, min, max;

   FieldVariable_GetMinAndMaxLocalCoords( sampleField, min, max );

   /* Sample Field in in regular grid */
   self->count = 0;
   /*printf("Sampling globally...\n");*/
   for ( i = 0 ; i <= self->resolution[I_AXIS] ; i++ )
   {
      for ( j = 0 ; j <= self->resolution[J_AXIS] ; j++ )
      {
         for ( k = 0 ; k <= self->resolution[K_AXIS] ; k++ )
         {
            double value;
            pos[ I_AXIS ] = min[ I_AXIS ] + self->cell[I_AXIS] * (double) i;
            pos[ J_AXIS ] = min[ J_AXIS ] + self->cell[J_AXIS] * (double) j;
            if (sampleField->dim == 2)
               pos[ K_AXIS ] = 0.0;
            else
               pos[ K_AXIS ] = min[ K_AXIS ] + self->cell[K_AXIS] * (double) k;

            FieldVariable_InterpolateValueAt( sampleField, pos, &value );
            //printf("%d,%d,%d  %f,%f,%f val = %f\n", i, j, k, pos[0], pos[1], pos[2], value);

            self->samples[self->count].value = (float)value;
            self->samples[self->count].pos[0] = (float)pos[0];
            self->samples[self->count].pos[1] = (float)pos[1];
            self->samples[self->count].pos[2] = (float)pos[2];
            self->count++;
         }
      }
   }
}

void _lucFieldSampler_Draw( void* drawingObject, lucDatabase* database, void* _context )
{
   lucFieldSampler* self          = (lucFieldSampler*)drawingObject;
   FieldVariable*   sampleField = self->sampleField;
   lucColourMap*    colourMap     = self->colourMap;
   int              i;
   Coord            min, max;

   /*printf("(%s) Resolution %d,%d,%d (dx/y/z %f,%f,%f)\n", self->name, self->resolution[0], self->resolution[1], self->resolution[2], self->cell[I_AXIS], self->cell[J_AXIS], self->cell[K_AXIS]);*/
   /* Calibrate Colour Map using Colour Variable */
   lucColourMap_CalibrateFromFieldVariable( colourMap, sampleField );

   FieldVariable_GetMinAndMaxLocalCoords( sampleField, min, max );

   for (i = 0 ; i < self->count; i++)
   {
      /*printf("%d = %f,%f,%f\n", i, self->samples[i].pos[0], self->samples[i].pos[1], self->samples[i].pos[2]);*/

      /* Plot Vertex */
      lucColourMap_SetColourFromValue( self->colourMap, self->samples[i].value, self->opacity );

      if (database)
      {
         /* Dump vertex pos, value */
         lucDatabase_AddVertices(database, 1, lucPointType, self->samples[i].pos);
         lucDatabase_AddValues(database, 1, lucPointType, lucColourValueData, colourMap, &self->samples[i].value);
      }
   }
}

