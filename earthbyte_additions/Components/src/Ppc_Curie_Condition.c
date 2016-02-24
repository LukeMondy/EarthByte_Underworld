#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include "Components.h"

#include <string.h>
#include <math.h>

#include "types.h"
#include "Ppc_Curie_Condition.h"


/* Textual name of this class */
const Type Ppc_Curie_Condition_Type = "Ppc_Curie_Condition";


/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
Ppc_Curie_Condition* _Ppc_Curie_Condition_New(  Ppc_Curie_Condition_DEFARGS  )
{
    Ppc_Curie_Condition* self;

    assert( _sizeOfSelf >= sizeof(Ppc_Curie_Condition) );
    nameAllocationType = NON_GLOBAL;
    self = (Ppc_Curie_Condition*) _Ppc_New(  PPC_PASSARGS  );    
    self->_get = _get;
    return self;
}


void* _Ppc_Curie_Condition_DefaultNew( Name name ) {
    /* Variables set in this function */
    SizeT                                              _sizeOfSelf = sizeof(Ppc_Curie_Condition);
    Type                                                      type = Ppc_Curie_Condition_Type;
    Stg_Class_DeleteFunction*                              _delete = _Ppc_Curie_Condition_Delete;
    Stg_Class_PrintFunction*                                _print = _Ppc_Curie_Condition_Print;
    Stg_Class_CopyFunction*                                  _copy = NULL;
    Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Ppc_Curie_Condition_DefaultNew;
    Stg_Component_ConstructFunction*                    _construct = _Ppc_Curie_Condition_AssignFromXML;
    Stg_Component_BuildFunction*                            _build = _Ppc_Curie_Condition_Build;
    Stg_Component_InitialiseFunction*                  _initialise = _Ppc_Curie_Condition_Initialise;
    Stg_Component_ExecuteFunction*                        _execute = _Ppc_Curie_Condition_Execute;
    Stg_Component_DestroyFunction*                        _destroy = _Ppc_Curie_Condition_Destroy;
    AllocationType                              nameAllocationType = NON_GLOBAL;
   Ppc_GetFunction*                                          _get = _Ppc_Curie_Condition_Get;

    return (void*)_Ppc_Curie_Condition_New(  Ppc_Curie_Condition_PASSARGS  );
}


void _Ppc_Curie_Condition_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data ) {
    Ppc_Curie_Condition* self = (Ppc_Curie_Condition*)_self;

    /* Construct parent */
    _Ppc_AssignFromXML( self, cf, data );    

    self->field = PpcManager_GetPpcFromDict( self->manager, cf, self->name, (Dictionary_Entry_Key)"TemperatureField", "" );
    self->value = PpcManager_GetPpcFromDict( self->manager, cf, self->name, (Dictionary_Entry_Key)"CuriePoint", "" );
    
    self->magorient = PpcManager_GetPpcFromDict( self->manager, cf, self->name, (Dictionary_Entry_Key)"MagneticOrientation", "" );

    self->storage = PpcManager_GetPpcFromDict( self->manager, cf, self->name, (Dictionary_Entry_Key)"StorageProperty", "" );
}   


void _Ppc_Curie_Condition_Build( void* _self, void* data ) {
    Ppc_Curie_Condition* self = (Ppc_Curie_Condition*)_self;

    /* Build parent */
    _Ppc_Build( self, data );
}

void _Ppc_Curie_Condition_Initialise( void* _self, void* data ) {
    Ppc_Curie_Condition* self = (Ppc_Curie_Condition*)_self;
    /* Initialize parent */
    _Ppc_Initialise( self, data );
}

void _Ppc_Curie_Condition_Delete( void* _self ) {
    Ppc_Curie_Condition* self = (Ppc_Curie_Condition*)_self;

    /* Delete parent */
    _Ppc_Delete( self );
}

void _Ppc_Curie_Condition_Print( void* _self, Stream* stream ) {
    Ppc_Curie_Condition* self = (Ppc_Curie_Condition*)_self;

    /* Print parent */
    _Ppc_Print( self, stream );
}

void _Ppc_Curie_Condition_Execute( void* _self, void* data ) {
    Ppc_Curie_Condition* self = (Ppc_Curie_Condition*)_self;

    /* Execute parent */
    _Ppc_Execute( self, data );
}

void _Ppc_Curie_Condition_Destroy( void* _self, void* data ) {
    Ppc_Curie_Condition* self = (Ppc_Curie_Condition*)_self;

    /* Destroy parent */
    _Ppc_Destroy( self, data );
}


Ppc_Curie_Condition* Ppc_Curie_Condition_New( Name name, Stg_Component* _self ) {
  Ppc_Curie_Condition* self = (Ppc_Curie_Condition*) _self;
  return self;
}

/* 
 * Public functions 
 *
 */

int _Ppc_Curie_Condition_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result ) {
    Ppc_Curie_Condition* self = (Ppc_Curie_Condition*) _self;
    double curiepoint, temperature, orientation, initialOrientation;
    double prevMagOrientation;
    Bool record_mag = False;
    int err = 0;

    /* get field value */
    err = PpcManager_Get( self->manager, lElement_I, particle, self->field, &temperature );
    assert( !err );

    /* get value to compare  */
    err = PpcManager_Get( self->manager, lElement_I, particle, self->value, &curiepoint );
    assert( !err );

    initialOrientation = -1;

    /* get previous value */
    err = PpcManager_GetPrevious( self->manager, lElement_I, particle, self->storage, &prevMagOrientation );
    assert( !err );


    record_mag = temperature <= curiepoint;

    
    if ( self->context->timeStep == 0) {
        result[0] = initialOrientation;  // Initialise the whole set to "no orientation"
    }
    else {
        if ( record_mag ) {     
            // we are below the curie point, we need to record it!
            err = PpcManager_Get( self->manager, lElement_I, particle, self->magorient, &orientation );
            assert( !err );

            if ( prevMagOrientation == initialOrientation )
                result[0] = orientation;         // Material has *just* dropped below the Curie point
            else
                result[0] = prevMagOrientation;  // Material has already dropped below the Curie point, 
                                                 // and does not need updating
        }
        else
        {            
            result[0] = initialOrientation;      // We have gone over the Curie point, so reset mag.
        }
    }
    return err;
}
