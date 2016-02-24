#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include "Components.h"

#include <string.h>
#include <math.h>

#include "types.h"
#include "Ppc_Melt_Polynomial.h"


/* Textual name of this class */
const Type Ppc_Melt_Polynomial_Type = "Ppc_Melt_Polynomial";


/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
Ppc_Melt_Polynomial* _Ppc_Melt_Polynomial_New(  PPC_LINEARDENSITY_DEFARGS  )
{
   Ppc_Melt_Polynomial* self;

   assert( _sizeOfSelf >= sizeof(Ppc_Melt_Polynomial) );
   nameAllocationType = NON_GLOBAL;
   self = (Ppc_Melt_Polynomial*) _Ppc_New(  PPC_PASSARGS  );
   self->_get = _get;
   return self;
}


void* _Ppc_Melt_Polynomial_DefaultNew( Name name )
{
   /* Variables set in this function */
   SizeT                                              _sizeOfSelf = sizeof(Ppc_Melt_Polynomial);
   Type                                                      type = Ppc_Melt_Polynomial_Type;
   Stg_Class_DeleteFunction*                              _delete = _Ppc_Melt_Polynomial_Delete;
   Stg_Class_PrintFunction*                                _print = _Ppc_Melt_Polynomial_Print;
   Stg_Class_CopyFunction*                                  _copy = NULL;
   Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Ppc_Melt_Polynomial_DefaultNew;
   Stg_Component_ConstructFunction*                    _construct = _Ppc_Melt_Polynomial_AssignFromXML;
   Stg_Component_BuildFunction*                            _build = _Ppc_Melt_Polynomial_Build;
   Stg_Component_InitialiseFunction*                  _initialise = _Ppc_Melt_Polynomial_Initialise;
   Stg_Component_ExecuteFunction*                        _execute = _Ppc_Melt_Polynomial_Execute;
   Stg_Component_DestroyFunction*                        _destroy = _Ppc_Melt_Polynomial_Destroy;
   AllocationType                              nameAllocationType = NON_GLOBAL;
   Ppc_GetFunction*                                          _get = _Ppc_Melt_Polynomial_Get;

   return (void*)_Ppc_Melt_Polynomial_New(  PPC_LINEARDENSITY_PASSARGS  );
}


void _Ppc_Melt_Polynomial_Init( void* _self, int term0, int term1, int term2, int term3, int pressure )
{
   Ppc_Melt_Polynomial* self = (Ppc_Melt_Polynomial*)_self;

   /* Sanity check */
   self->t0Tag = term0;
   self->t1Tag = term1;
   self->t2Tag = term2;
   self->t3Tag = term3;
   self->pressureTag = pressure;
}


void _Ppc_Melt_Polynomial_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data )
{
   Ppc_Melt_Polynomial* self = (Ppc_Melt_Polynomial*)_self;

   /* Construct parent */
   _Ppc_AssignFromXML( self, cf, data );

   _Ppc_Melt_Polynomial_Init(
      self,
      PpcManager_GetPpcFromDict( self->manager, cf, (Name)self->name, "term0", "" ),
      PpcManager_GetPpcFromDict( self->manager, cf, (Name)self->name, "term1", "" ),
      PpcManager_GetPpcFromDict( self->manager, cf, (Name)self->name, "term2", "" ),
      PpcManager_GetPpcFromDict( self->manager, cf, (Name)self->name, "term3", "" ),
      PpcManager_GetPpcFromDict( self->manager, cf, (Name)self->name, "Pressure", "" ));
}


void _Ppc_Melt_Polynomial_Build( void* _self, void* data )
{
   Ppc_Melt_Polynomial* self = (Ppc_Melt_Polynomial*)_self;

   /* Build parent */
   _Ppc_Build( self, data );
}

void _Ppc_Melt_Polynomial_Initialise( void* _self, void* data )
{
   Ppc_Melt_Polynomial* self = (Ppc_Melt_Polynomial*)_self;

   /* Initialize parent */
   _Ppc_Initialise( self, data );
}

void _Ppc_Melt_Polynomial_Delete( void* _self )
{
   Ppc_Melt_Polynomial* self = (Ppc_Melt_Polynomial*)_self;

   /* Delete parent */
   _Ppc_Delete( self );
}

void _Ppc_Melt_Polynomial_Print( void* _self, Stream* stream )
{
   Ppc_Melt_Polynomial* self = (Ppc_Melt_Polynomial*)_self;

   /* Print parent */
   _Ppc_Print( self, stream );
}

void _Ppc_Melt_Polynomial_Execute( void* _self, void* data )
{
   Ppc_Melt_Polynomial* self = (Ppc_Melt_Polynomial*)_self;

   /* Execute parent */
   _Ppc_Execute( self, data );
}

void _Ppc_Melt_Polynomial_Destroy( void* _self, void* data )
{
   Ppc_Melt_Polynomial* self = (Ppc_Melt_Polynomial*)_self;

   /* Destroy parent */
   _Ppc_Destroy( self, data );
}


Ppc_Melt_Polynomial* Ppc_Melt_Polynomial_New( Name name, Stg_Component* _self )
{
   Ppc_Melt_Polynomial* self = (Ppc_Melt_Polynomial*) _self;
   return self;
}

/*
 * Public functions
 *
 */

int _Ppc_Melt_Polynomial_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result )
{
   Ppc_Melt_Polynomial* self = (Ppc_Melt_Polynomial*) _self;
   double t0, t1, t2, t3, p;
   int err;

   err = PpcManager_Get( self->manager, lElement_I, particle, self->t0Tag, &t0 ); assert( !err );
   err = PpcManager_Get( self->manager, lElement_I, particle, self->t1Tag, &t1 ); assert( !err );
   err = PpcManager_Get( self->manager, lElement_I, particle, self->t2Tag, &t2 ); assert( !err );
   err = PpcManager_Get( self->manager, lElement_I, particle, self->t3Tag, &t3 ); assert( !err );
   err = PpcManager_Get( self->manager, lElement_I, particle, self->pressureTag, &p ); assert( !err );
   

   result[0] = t0 + t1*p + t2*(p*p) + t3*(p*p*p);
   

   return 0;
}
