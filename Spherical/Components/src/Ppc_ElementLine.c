#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Components/Components.h>

#include <string.h>
#include <math.h>
#include <float.h>

#include "types.h"
#include "Ppc_ElementLine.h"


/* Textual name of this class */
const Type Ppc_ElementLine_Type = "Ppc_ElementLine";


/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
Ppc_ElementLine* _Ppc_ElementLine_New(  PPC_CONSTANT_DEFARGS  )
{
   Ppc_ElementLine* self;

   assert( _sizeOfSelf >= sizeof(Ppc_ElementLine) );
   nameAllocationType = NON_GLOBAL;
   self = (Ppc_ElementLine*) _Ppc_New(  PPC_PASSARGS  );
   self->_get = _get;
   return self;
}


void* _Ppc_ElementLine_DefaultNew( Name name )
{
   /* Variables set in this function */
   SizeT                                              _sizeOfSelf = sizeof(Ppc_ElementLine);
   Type                                                      type = Ppc_ElementLine_Type;
   Stg_Class_DeleteFunction*                              _delete = _Ppc_ElementLine_Delete;
   Stg_Class_PrintFunction*                                _print = _Ppc_ElementLine_Print;
   Stg_Class_CopyFunction*                                  _copy = NULL;
   Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Ppc_ElementLine_DefaultNew;
   Stg_Component_ConstructFunction*                    _construct = _Ppc_ElementLine_AssignFromXML;
   Stg_Component_BuildFunction*                            _build = _Ppc_ElementLine_Build;
   Stg_Component_InitialiseFunction*                  _initialise = _Ppc_ElementLine_Initialise;
   Stg_Component_ExecuteFunction*                        _execute = _Ppc_ElementLine_Execute;
   Stg_Component_DestroyFunction*                        _destroy = _Ppc_ElementLine_Destroy;
   AllocationType                              nameAllocationType = NON_GLOBAL;
   Ppc_GetFunction*                                          _get = _Ppc_ElementLine_Get;

   return (void*)_Ppc_ElementLine_New(  PPC_CONSTANT_PASSARGS  );
}


void _Ppc_ElementLine_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data )
{
   Ppc_ElementLine* self = (Ppc_ElementLine*)_self;

   /* Construct parent */
   _Ppc_AssignFromXML( self, cf, data );

}


void _Ppc_ElementLine_Build( void* _self, void* data )
{
   Ppc_ElementLine* self = (Ppc_ElementLine*)_self;

   /* Build parent */
   _Ppc_Build( self, data );

}

void _Ppc_ElementLine_Initialise( void* _self, void* data )
{
   Ppc_ElementLine* self = (Ppc_ElementLine*)_self;

   /* Initialize parent */
   _Ppc_Initialise( self, data );
}

void _Ppc_ElementLine_Delete( void* _self )
{
   Ppc_ElementLine* self = (Ppc_ElementLine*)_self;

   /* Delete parent */
   _Ppc_Delete( self );
}

void _Ppc_ElementLine_Print( void* _self, Stream* stream )
{
   Ppc_ElementLine* self = (Ppc_ElementLine*)_self;

   /* Print parent */
   _Ppc_Print( self, stream );
}

void _Ppc_ElementLine_Execute( void* _self, void* data )
{
   Ppc_ElementLine* self = (Ppc_ElementLine*)_self;

   /* Execute parent */
   _Ppc_Execute( self, data );
}

void _Ppc_ElementLine_Destroy( void* _self, void* data )
{
   Ppc_ElementLine* self = (Ppc_ElementLine*)_self;

   /* Destroy parent */
   _Ppc_Destroy( self, data );
}

/*
 * Public functions
 *
 */

int _Ppc_ElementLine_Get( void* _self, unsigned lElement_I, IntegrationPoint* particle, double* result )
{
	Ppc_ElementLine* self = (Ppc_ElementLine*) _self;
	unsigned ind[3];
	FeMesh *mesh = self->manager->mesh;
	Grid **grid = (Grid** )Mesh_GetExtension(mesh,Grid**,mesh->elGridId);

	Grid_Lift( grid[0], lElement_I, ind );

	int endE = 19;
	int mid = endE / 2;

	if( MAX( abs(ind[0]-mid), abs(ind[1]-mid) ) == 4 ) {
		result[0] = 3;
		return 0;
	}

	if( ind[1] == 4 ) {
		result[0] = 2;
		return 0;
	}

	if( ind[0] == 2 ) {
		result[0] = 1;
		return 0;
	}

	result[0] = 0;
		
  return 0;
}
