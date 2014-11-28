#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include <string.h>
#include <math.h>
#include <float.h>

#include "types.h"
#include "Ppc_VD.h"


/* Textual name of this class */
const Type Ppc_VD_Type = "Ppc_VD";


/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
Ppc_VD* _Ppc_VD_New(  PPC_CONSTANT_DEFARGS  )
{
   Ppc_VD* self;

   assert( _sizeOfSelf >= sizeof(Ppc_VD) );
   nameAllocationType = NON_GLOBAL;
   self = (Ppc_VD*) _Ppc_New(  PPC_PASSARGS  );
   self->_get = _get;
   return self;
}


void* _Ppc_VD_DefaultNew( Name name )
{
   /* Variables set in this function */
   SizeT                                              _sizeOfSelf = sizeof(Ppc_VD);
   Type                                                      type = Ppc_VD_Type;
   Stg_Class_DeleteFunction*                              _delete = _Ppc_VD_Delete;
   Stg_Class_PrintFunction*                                _print = _Ppc_VD_Print;
   Stg_Class_CopyFunction*                                  _copy = NULL;
   Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Ppc_VD_DefaultNew;
   Stg_Component_ConstructFunction*                    _construct = _Ppc_VD_AssignFromXML;
   Stg_Component_BuildFunction*                            _build = _Ppc_VD_Build;
   Stg_Component_InitialiseFunction*                  _initialise = _Ppc_VD_Initialise;
   Stg_Component_ExecuteFunction*                        _execute = _Ppc_VD_Execute;
   Stg_Component_DestroyFunction*                        _destroy = _Ppc_VD_Destroy;
   AllocationType                              nameAllocationType = NON_GLOBAL;
   Ppc_GetFunction*                                          _get = _Ppc_VD_Get;

   return (void*)_Ppc_VD_New(  PPC_CONSTANT_PASSARGS  );
}


void _Ppc_VD_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data )
{
   Ppc_VD* self = (Ppc_VD*)_self;

   /* The PpcManager */
	self->manager = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Manager", PpcManager, False, data );
	if( !self->manager  )
		self->manager = Stg_ComponentFactory_ConstructByName( cf, (Name)"default_ppcManager", PpcManager, True, data );
   /* Construct parent */
   _Ppc_AssignFromXML( self, cf, data );
	self->strainRateTag = PpcManager_GetPpcFromDict( self->manager, cf, self->name, "StrainRate", "" );
	self->viscosityTag = PpcManager_GetPpcFromDict( self->manager, cf, self->name, "Viscosity", "" );
}


void _Ppc_VD_Build( void* _self, void* data )
{
   Ppc_VD* self = (Ppc_VD*)_self;

   /* Build parent */
   _Ppc_Build( self, data );
}

void _Ppc_VD_Initialise( void* _self, void* data )
{
   Ppc_VD* self = (Ppc_VD*)_self;

   /* Initialize parent */
   _Ppc_Initialise( self, data );
}

void _Ppc_VD_Delete( void* _self )
{
   Ppc_VD* self = (Ppc_VD*)_self;

   /* Delete parent */
   _Ppc_Delete( self );
}

void _Ppc_VD_Print( void* _self, Stream* stream )
{
   Ppc_VD* self = (Ppc_VD*)_self;

   /* Print parent */
   _Ppc_Print( self, stream );
}

void _Ppc_VD_Execute( void* _self, void* data )
{
   Ppc_VD* self = (Ppc_VD*)_self;

   /* Execute parent */
   _Ppc_Execute( self, data );
}

void _Ppc_VD_Destroy( void* _self, void* data )
{
   Ppc_VD* self = (Ppc_VD*)_self;

   /* Destroy parent */
   _Ppc_Destroy( self, data );
}

/*
 * Public functions
 *
 */

int _Ppc_VD_Get( void* _self, unsigned lElement_I, IntegrationPoint* particle, double* result )
{
   Ppc_VD* self = (Ppc_VD*) _self;
   int err;
   double vd, tau[6], sr[6], viscosity;
   double tmp;

   err = PpcManager_Get( self->manager, lElement_I, particle, self->strainRateTag, sr );
   assert(!err);

   err = PpcManager_Get( self->manager, lElement_I, particle, self->viscosityTag, &viscosity );
   assert(!err);

   // build a stress tensor - assume isotropic rheology
   tau[0] = 2 * viscosity * sr[0];
   tau[1] = 2 * viscosity * sr[1];
   tau[2] = 2 * viscosity * sr[2];

   // take the tensor inner product of stress and strain-rate 
   vd = tau[0]*sr[0] + tau[1]*sr[1] + 2*tau[2]*sr[2];


   /* i was confused about a 2 and this was testing code
   tmp = 2 * viscosity * 0.5*(sr[0]*sr[0] + sr[1]*sr[1] + 2*sr[2]*sr[2]);

   if( fabs(tmp-vd) > 1e-3 ) {
     printf("VD issue t:e = %g while 2*eta*eII*eII = %g \n", vd, tmp );
   }
   */

   *result = vd;

   return 0;
}
