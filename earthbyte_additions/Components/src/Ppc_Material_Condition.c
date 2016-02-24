#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include "Components.h"

#include <string.h>
#include <math.h>

#include "types.h"
#include "Ppc_Material_Condition.h"


/* Textual name of this class */
const Type Ppc_Material_Condition_Type = "Ppc_Material_Condition";


/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
Ppc_Material_Condition* _Ppc_Material_Condition_New(  Ppc_Material_Condition_DEFARGS  )
{
    Ppc_Material_Condition* self;

    assert( _sizeOfSelf >= sizeof(Ppc_Material_Condition) );
    nameAllocationType = NON_GLOBAL;
    self = (Ppc_Material_Condition*) _Ppc_New(  PPC_PASSARGS  );    
    self->_get = _get;
    return self;
}


void* _Ppc_Material_Condition_DefaultNew( Name name ) {
    /* Variables set in this function */
    SizeT                                              _sizeOfSelf = sizeof(Ppc_Material_Condition);
    Type                                                      type = Ppc_Material_Condition_Type;
    Stg_Class_DeleteFunction*                              _delete = _Ppc_Material_Condition_Delete;
    Stg_Class_PrintFunction*                                _print = _Ppc_Material_Condition_Print;
    Stg_Class_CopyFunction*                                  _copy = NULL;
    Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Ppc_Material_Condition_DefaultNew;
    Stg_Component_ConstructFunction*                    _construct = _Ppc_Material_Condition_AssignFromXML;
    Stg_Component_BuildFunction*                            _build = _Ppc_Material_Condition_Build;
    Stg_Component_InitialiseFunction*                  _initialise = _Ppc_Material_Condition_Initialise;
    Stg_Component_ExecuteFunction*                        _execute = _Ppc_Material_Condition_Execute;
    Stg_Component_DestroyFunction*                        _destroy = _Ppc_Material_Condition_Destroy;
    AllocationType                              nameAllocationType = NON_GLOBAL;
   Ppc_GetFunction*                                          _get = _Ppc_Material_Condition_Get;

    return (void*)_Ppc_Material_Condition_New(  Ppc_Material_Condition_PASSARGS  );
}


void _Ppc_Material_Condition_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data ) {
    Ppc_Material_Condition* self = (Ppc_Material_Condition*)_self;

    /* Construct parent */
    _Ppc_AssignFromXML( self, cf, data );    
    self->condition =Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"Condition", "" );
    self->field = PpcManager_GetPpcFromDict( self->manager, cf, self->name, (Dictionary_Entry_Key)"Field", "" );
    self->value = PpcManager_GetPpcFromDict( self->manager, cf, self->name, (Dictionary_Entry_Key)"ValueToCompare", "" );
    
    self->materialTrue_name = Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"MaterialIfTrue", "" );
    self->materialFalse_name = Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"MaterialIfFalse", "" );
    
    self->currentTimeProperty = PpcManager_GetPpcFromDict( self->manager, cf, self->name, (Dictionary_Entry_Key)"CurrentTime", "" );
    self->timeStorage = PpcManager_GetPpcFromDict( self->manager, cf, self->name, (Dictionary_Entry_Key)"StorageProperty", "" );
}   


void _Ppc_Material_Condition_Build( void* _self, void* data ) {
    Ppc_Material_Condition* self = (Ppc_Material_Condition*)_self;

    /* Build parent */
    _Ppc_Build( self, data );
}

void _Ppc_Material_Condition_Initialise( void* _self, void* data ) {
    Ppc_Material_Condition* self = (Ppc_Material_Condition*)_self;
    Material* materialTrue;
    Material* materialFalse;

    /* Initialize parent */
    _Ppc_Initialise( self, data );
    
    materialTrue = Materials_Register_GetByName( self->manager->materialSwarm->materials_Register, self->materialTrue_name );
    Journal_Firewall( (Bool)materialTrue, 
                            self->error_stream,
                            "Error in configuration of Ppc_IsMaterial %s.\n"
                            "No material with name %s found.",
                            self->name, self->materialTrue_name);
    self->materialTrue = materialTrue->index;
    
    if( strcmp(self->materialFalse_name, "") != 0 ) {
        materialFalse = Materials_Register_GetByName( self->manager->materialSwarm->materials_Register, self->materialFalse_name );
        Journal_Firewall( (Bool)materialFalse, 
                                self->error_stream,
                                "Error in configuration of Ppc_IsMaterial %s.\n"
                                "No material with name %s found.",
                                self->name, self->materialFalse_name);
        self->materialFalse = materialFalse->index;
    }
    else {
        // No material has been given, so it means any material is OK.
        self->materialFalse = -1;
    }
        
    
}

void _Ppc_Material_Condition_Delete( void* _self ) {
    Ppc_Material_Condition* self = (Ppc_Material_Condition*)_self;

    /* Delete parent */
    _Ppc_Delete( self );
}

void _Ppc_Material_Condition_Print( void* _self, Stream* stream ) {
    Ppc_Material_Condition* self = (Ppc_Material_Condition*)_self;

    /* Print parent */
    _Ppc_Print( self, stream );
}

void _Ppc_Material_Condition_Execute( void* _self, void* data ) {
    Ppc_Material_Condition* self = (Ppc_Material_Condition*)_self;

    /* Execute parent */
    _Ppc_Execute( self, data );
}

void _Ppc_Material_Condition_Destroy( void* _self, void* data ) {
    Ppc_Material_Condition* self = (Ppc_Material_Condition*)_self;

    /* Destroy parent */
    _Ppc_Destroy( self, data );
}


Ppc_Material_Condition* Ppc_Material_Condition_New( Name name, Stg_Component* _self ) {
  Ppc_Material_Condition* self = (Ppc_Material_Condition*) _self;
  return self;
}

/* 
 * Public functions 
 *
 */

int _Ppc_Material_Condition_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result ) {
    Ppc_Material_Condition* self = (Ppc_Material_Condition*) _self;
    double comp_value, fieldValue;
    double currentTime = -1.0;
    double previousTime;
    char* op = self->condition;
    Bool conditionResult = False;
    Bool oneToMany;
    int err = 0;

    MaterialPoint* material = NULL;
    MaterialPoint* mat;

    MaterialPointsSwarm *ms = NULL;

    /* get field value */
    err = PpcManager_Get( self->manager, lElement_I, particle, self->field, &fieldValue );
    assert( !err );

    /* get value to compare  */
    err = PpcManager_Get( self->manager, lElement_I, particle, self->value, &comp_value );
    assert( !err );

    switch( op[0] ) {
        case 'e':
            conditionResult = fieldValue == comp_value;
            break;
        case 'n':
            conditionResult = fieldValue != comp_value;
            break;
        case 'l':
            if( op[1] == 'e' )
                conditionResult = fieldValue <= comp_value;
            else
                conditionResult = fieldValue < comp_value;
            break;
        case 'g':
            if( op[1] == 'e' )
                conditionResult = fieldValue >= comp_value;
            else
                conditionResult = fieldValue > comp_value;
            break;
    }
    
    if ( self->currentTimeProperty != -1 )
    {
        err = PpcManager_Get( self->manager, lElement_I, particle, self->currentTimeProperty, &currentTime);
        assert( !err );
        err = PpcManager_GetPrevious( self->manager, lElement_I, particle, self->timeStorage, &previousTime );
           assert( !err );
    }

    oneToMany = Stg_Class_IsInstance(self->manager->integrationSwarm->mapper, OneToManyMapper_Type);

    
    if ( conditionResult ) {
        if (!oneToMany)
        {
            //mat = IntegrationPointsSwarm_GetMaterialOn(self->manager->integrationSwarm->mapper, particle);
            mat = OneToOneMapper_GetMaterialPoint( self->manager->integrationSwarm->mapper, particle, &ms );
            mat->materialIndex = self->materialTrue;
        }
        else
        {
            OneToManyRef *ref;
            int kk;
            ref = OneToManyMapper_GetMaterialRef(self->manager->integrationSwarm->mapper, particle);

            for (kk = 0; kk < ref->numParticles; kk++)
            {
               //mat = MaterialPointsSwarm_GetMaterialAt(((OneToManyMapper*)self->manager->integrationSwarm->mapper)->materialSwarm, ref->particleInds[kk]);
               MaterialPointsSwarm_SetMaterialAt(((OneToManyMapper*)self->manager->integrationSwarm->mapper)->materialSwarm, ref->particleInds[kk], self->materialTrue);
               //mat->materialIndex = self->materialTrue;
               
            }
        }

        if ( self->currentTimeProperty != -1 )    // If currentTime has been supplied, output the current time.
        {
            result[0] = currentTime;
        }
        else
            result[0] = 1.0;
    }

    else {
        //if ( self->materialFalse != -1 ) // If a false material is defined
        //    material->materialIndex = self->materialFalse;    // NOTE: This stuff may never get called... maybe?
            
        if ( self->currentTimeProperty != -1 )
        {
            if ( currentTime == 0 )
                result[0] = 0;
            else
                result[0] = previousTime;    // If currentTime has been supplied, we need to make sure we preserve
        }
        else                            // the value that was found before.
            result[0] = 0.0; 
    }

    return err;
}
