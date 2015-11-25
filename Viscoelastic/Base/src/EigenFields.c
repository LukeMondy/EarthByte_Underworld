/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005-2010, Monash University
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
**       * Redistributions of source code must retain the above copyright notice,
**          this list of conditions and the following disclaimer.
**       * Redistributions in binary form must reproduce the above copyright
**       notice, this list of conditions and the following disclaimer in the
**       documentation and/or other materials provided with the distribution.
**       * Neither the name of the Monash University nor the names of its contributors
**       may be used to endorse or promote products derived from this software
**       without specific prior written permission.
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
*% Development Team :
*%  http://www.underworldproject.org/aboutus.html
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include "types.h"
#include "EigenFields.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type EigenFields_Type = "EigenFields";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

EigenFields* _EigenFields_New(  EIGENFIELDS_DEFARGS  )
{
   EigenFields* self;

   /* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
   assert( _sizeOfSelf >= sizeof(EigenFields) );
   /* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
   /* This means that any values of these parameters that are passed into this function are not passed onto the parent function
      and so should be set to ZERO in any children of this class. */

   self = (EigenFields*)_Stg_Component_New(  STG_COMPONENT_PASSARGS  );
   
   self->context = NULL;
   self->eigenValueField[0]  = NULL;
   self->eigenValueField[1]  = NULL;
   self->eigenValueField[2]  = NULL;
   self->eigenVectorField[0] = NULL;
   self->eigenVectorField[1] = NULL;
   self->eigenVectorField[2] = NULL;
   self->lastEvaluationStep  = -1;
   /* General info */

   /* Virtual Info */

   return self;
}

void _EigenFields_Init(
         void*             _self,
         DomainContext*    context,
         FeVariable*       tensorField,
         FeVariable**      eigenValueField,
         FeVariable**      eigenVectorField )
{
   EigenFields* self = (EigenFields*)_self;
   Stream* errorStream = Journal_Register( Error_Type, (Name)self->type );
   unsigned ii;

   self->context          = context;
   self->tensorField      = tensorField;
   memcpy( self->eigenValueField, eigenValueField, self->context->dim*sizeof(FeVariable*));
   memcpy( self->eigenVectorField, eigenVectorField, self->context->dim*sizeof(FeVariable*));

   /* grab femesh off first eigenfield, and ensure all fields use the same mesh */
   self->eigenMesh = eigenValueField[0]->feMesh;
   for(ii = 0; ii<self->context->dim; ii++){
      Journal_Firewall( self->eigenMesh == eigenValueField[ii]->feMesh, errorStream,
        "\n\nError:  All eigen fields should use the same mesh.\n\n" );
      Journal_Firewall( self->eigenMesh == eigenVectorField[ii]->feMesh, errorStream,
        "\n\nError:  All eigen fields should use the same mesh.\n\n" );
   }
   /* ensure that eigenMesh is identical to tensorField mesh */
   Journal_Firewall( self->eigenMesh == tensorField->feMesh, errorStream,
        "\n\nError:  eigen fields should use the same mesh as the tensorField.\n\n" );
     
   /* Add self as extensions to facilitate auto refresh of data */
   for(ii = 0; ii<self->context->dim; ii++){
      ExtensionManager_Add(  self->eigenValueField[ii]->extensionMgr, (Name)"OwningEigenFieldInstance", sizeof(EigenFields*) );
      ExtensionManager_Add( self->eigenVectorField[ii]->extensionMgr, (Name)"OwningEigenFieldInstance", sizeof(EigenFields*) );
      *(EigenFields** )ExtensionManager_Get(  self->eigenValueField[ii]->extensionMgr,  self->eigenValueField[ii], ExtensionManager_GetHandle(  self->eigenValueField[ii]->extensionMgr, (Name)"OwningEigenFieldInstance" ) ) = self;
      *(EigenFields** )ExtensionManager_Get( self->eigenVectorField[ii]->extensionMgr, self->eigenVectorField[ii], ExtensionManager_GetHandle( self->eigenVectorField[ii]->extensionMgr, (Name)"OwningEigenFieldInstance" ) ) = self;
   }
   /* modify redirect InterpolateAt and ValueAtNode functions to use custom functions */
   for(ii = 0; ii<self->context->dim; ii++){
      self->eigenValueField[ii]->_interpolateValueAt  = _EigenFields_InterpolateValueAt;
      self->eigenValueField[ii]->_getValueAtNode      = _EigenFields_GetValueAtNode;
      self->eigenVectorField[ii]->_interpolateValueAt = _EigenFields_InterpolateValueAt;
      self->eigenVectorField[ii]->_getValueAtNode     = _EigenFields_GetValueAtNode;
  }

  /* sync routines run on this guy envoked on an entry point.  this avoids barriers issues that can occur where processors are waiting at different points of the simulation */ 
  EP_AppendClassHook( Context_GetEntryPoint( self->context, AbstractContext_EP_PostSolvePreUpdateClass ), _EigenFields_SyncFields, self );

}


/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _EigenFields_Delete( void* _self ) {
   EigenFields* self = (EigenFields*)_self;
   _Stg_Component_Delete( self );
}

void  _EigenFields_Print( void* eigenFields, Stream* stream ) { assert(0); }
void* _EigenFields_Copy( void* eigenFields, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) { assert(0); }

void* _EigenFields_DefaultNew( Name name ) {
   /* Variables set in this function */
   SizeT                                                _sizeOfSelf = sizeof(EigenFields);
   Type                                                        type = EigenFields_Type;
   Stg_Class_DeleteFunction*                                _delete = _EigenFields_Delete;
   Stg_Class_PrintFunction*                                  _print = _EigenFields_Print;
   Stg_Class_CopyFunction*                                    _copy = _EigenFields_Copy;
   Stg_Component_DefaultConstructorFunction*    _defaultConstructor = _EigenFields_DefaultNew;
   Stg_Component_ConstructFunction*                      _construct = _EigenFields_AssignFromXML;
   Stg_Component_BuildFunction*                              _build = _EigenFields_Build;
   Stg_Component_InitialiseFunction*                    _initialise = _EigenFields_Initialise;
   Stg_Component_ExecuteFunction*                          _execute = _EigenFields_Execute;
   Stg_Component_DestroyFunction*                          _destroy = _EigenFields_Destroy;

   /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
   AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

   return (void*) _EigenFields_New(  EIGENFIELDS_PASSARGS  );
}


void _EigenFields_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data ) {
   EigenFields*         self     = (EigenFields*) _self;
   DomainContext*       context;
   FeVariable*          tensorField;
   FeVariable*          eigenValueField[3];
   FeVariable*          eigenVectorField[3];

   context = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Context", DomainContext, False, data );
   if( !context  )
      context = Stg_ComponentFactory_ConstructByName( cf, (Name)"context", DomainContext, True, data  );

   /** the tensor field from which will be extracted eigenvalues and eigenvectors */ 
   tensorField        = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"TensorField", FeVariable, True, data  );
   /** eigenvalue fields.  these are to be instantiated from xml, usually from the default EigenFields.xml file */
   eigenValueField[0] = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"EigenValueField1", FeVariable, True, data  );
   eigenValueField[1] = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"EigenValueField2", FeVariable, True, data  );
   if(context->dim == 3)
      eigenValueField[2] = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"EigenValueField3", FeVariable, True, data  );
   /** eigenvector fields.  these are to be instantiated from xml, usually from the default EigenFields.xml file */
   eigenVectorField[0] = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"EigenVectorField1", FeVariable, True, data  );
   eigenVectorField[1] = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"EigenVectorField2", FeVariable, True, data  );
   if(context->dim == 3)
      eigenVectorField[2] = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"EigenVectorField3", FeVariable, True, data  );
   
   _EigenFields_Init( _self, context, tensorField, eigenValueField, eigenVectorField );
}

void _EigenFields_Build( void* _self, void* data ) {
   EigenFields* self = (EigenFields*) _self;
   Stream* errorStream = Journal_Register( Error_Type, (Name)self->type );
   unsigned ii;

   Stg_Component_Build( self->context, data, False );
   Stg_Component_Build( self->tensorField, data, False );
   for(ii = 0; ii<self->context->dim; ii++){
      Stg_Component_Build( self->eigenValueField[ii] , data, False );
      Stg_Component_Build( self->eigenVectorField[ii], data, False );
   } 
  /* check tensorField component counts, ensure symmetric */
   if(self->context->dim == 2)
      Journal_Firewall( self->tensorField->fieldComponentCount == 3, errorStream, "\n\nEigenFields component only supports symmetric tensor fields.\nCurrent tensor field component count is %u.  Count of 3 is required.\nYou can use an OperatorFeVariable to take the symmetric part.  Check VelocityField.xml for an example.\n\n", self->tensorField->fieldComponentCount);
   else
      Journal_Firewall( self->tensorField->fieldComponentCount == 6, errorStream, "\n\nEigenFields component only supports symmetric tensor fields.\nCurrent tensor field component count is %u.  Count of 6 is required.\nYou can use an OperatorFeVariable to take the symmetric part.  Check VelocityField.xml for an example.\n\n", self->tensorField->fieldComponentCount);


}

void _EigenFields_Initialise( void* _self, void* data ) {
   EigenFields* self = (EigenFields*) _self;
   unsigned ii;
   
   Stg_Component_Initialise( self->tensorField, data, False );
   for(ii = 0; ii<self->context->dim; ii++){
      Stg_Component_Initialise( self->eigenValueField[ii] , data, False );
      Stg_Component_Initialise( self->eigenVectorField[ii], data, False );
   }
   /** now load all values onto fields */
   _EigenFields_UpdateFields( _self );
   _EigenFields_SyncFields( _self, self->context );

}

void _EigenFields_Execute( void* _self, void* data ) {}

void _EigenFields_Destroy( void* _self, void* data ) {
   EigenFields* self = (EigenFields*) _self;
   unsigned ii;

   Stg_Component_Destroy( self->tensorField, data, False );
   for(ii = 0; ii<self->context->dim; ii++){
      Stg_Component_Destroy( self->eigenValueField[ii] , data, False );
      Stg_Component_Destroy( self->eigenVectorField[ii], data, False );
   }

}

InterpolationResult _EigenFields_InterpolateValueAt( void* variable, double* globalCoord, double* value ) {
   FeVariable* feVariable = (FeVariable*) variable;
   EigenFields* self = *(EigenFields**)ExtensionManager_Get(  feVariable->extensionMgr,  feVariable, ExtensionManager_GetHandle(  feVariable->extensionMgr, (Name)"OwningEigenFieldInstance" ) );

   if( self->lastEvaluationStep < (int) self->context->timeStep ) {
      _EigenFields_UpdateFields( (void*) self );
      self->lastEvaluationStep = (int)self->context->timeStep;
   }

   return _FeVariable_InterpolateValueAt( variable, globalCoord, value ); 
}

void _EigenFields_GetValueAtNode( void* _feVariable, Node_DomainIndex dNode_I, double* value ) {
   FeVariable* feVariable = (FeVariable*) _feVariable;
   EigenFields* self = *(EigenFields**)ExtensionManager_Get(  feVariable->extensionMgr,  feVariable, ExtensionManager_GetHandle(  feVariable->extensionMgr, (Name)"OwningEigenFieldInstance" ) );

   if( self->lastEvaluationStep < (int) self->context->timeStep ) {
      _EigenFields_UpdateFields( (void*) self );
      self->lastEvaluationStep = (int)self->context->timeStep;
   }
  
   return _FeVariable_GetValueAtNode( _feVariable, dNode_I, value );
}

void _EigenFields_UpdateFields( void* _self ) {
   EigenFields* self = (EigenFields*)_self;
   Eigenvector eigenVectors[3]; 
   Node_LocalIndex  lNode_I, lNode_max;
   SymmetricTensor tensor;
   unsigned ii;

   lNode_max = Mesh_GetLocalSize( self->eigenMesh, MT_VERTEX );

   /* step through nodes on tensorField, calculate eigen data, set values onto field */
   for( lNode_I=0; lNode_I<lNode_max; lNode_I++ ) {
      FeVariable_GetValueAtNode( (void*)self->tensorField, lNode_I,(double*)&tensor );
      SymmetricTensor_CalcAllEigenvectors( tensor, self->context->dim, eigenVectors );
      for(ii=0; ii<self->context->dim; ii++){
         FeVariable_SetValueAtNode( self->eigenValueField[ii] , lNode_I, &(eigenVectors[ii].eigenvalue) );
         FeVariable_SetValueAtNode( self->eigenVectorField[ii], lNode_I, eigenVectors[ii].vector     );
      }
   }
   
}

void _EigenFields_SyncFields( void* _self, void* context ) {
   EigenFields* self = (EigenFields*)_self;
   unsigned ii;

   for(ii=0; ii<self->context->dim; ii++){
      FeVariable_SyncShadowValues( self->eigenValueField[ii] );
      FeVariable_SyncShadowValues( self->eigenVectorField[ii] );
   }

}

