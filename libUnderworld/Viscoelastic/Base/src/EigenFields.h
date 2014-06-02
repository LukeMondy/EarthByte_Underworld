/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**      Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**      Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**      Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**      Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**      Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**      Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**
** This file may be distributed under the terms of the VPAC Public License
** as defined by VPAC of Australia and appearing in the file
** LICENSE.VPL included in the packaging of this file.
**
** This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
** WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
**
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Experimental_Utils_EigenFields_h__
#define __Experimental_Utils_EigenFields_h__

   /* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
   extern const Type EigenFields_Type;

   /* EigenFields information */
   #define __EigenFields                                                                            \
      /* Macro defining parent goes here - This means you can cast this class as its parent */      \
      __Stg_Component                                                                               \
      /* Virtual Info */                                                                            \
      DomainContext*              context;                                                          \
      FeVariable*                 tensorField;                                                      \
      FeVariable*                 eigenValueField[3];                                               \
      FeVariable*                 eigenVectorField[3];                                              \
      FeMesh*                     eigenMesh;                                                        \
      int                         lastEvaluationStep;

   struct EigenFields { __EigenFields };
   
   /*---------------------------------------------------------------------------------------------------------------------
   ** Constructors
   */
   
   #ifndef ZERO
   #define ZERO 0
   #endif
   
   #define EIGENFIELDS_DEFARGS \
           STG_COMPONENT_DEFARGS
   
   #define EIGENFIELDS_PASSARGS \
           STG_COMPONENT_PASSARGS
   
   EigenFields* _EigenFields_New(  EIGENFIELDS_DEFARGS  );
   
   /* Stg_Class_Delete EigenFields implementation */
   void _EigenFields_Delete( void* eigenFields );
   void _EigenFields_Print( void* eigenFields, Stream* stream );
   #define EigenFields_Copy( self ) \
          (EigenFields*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
   #define EigenFields_DeepCopy( self ) \
          (EigenFields*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
   void* _EigenFields_Copy( void* eigenFields, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
   void* _EigenFields_DefaultNew( Name name ) ;
   void _EigenFields_AssignFromXML( void* shape, Stg_ComponentFactory* cf, void* data ) ;
   void _EigenFields_Build( void* eigenFields, void* data ) ;
   void _EigenFields_Initialise( void* eigenFields, void* data ) ;
   void _EigenFields_Execute( void* eigenFields, void* data );
   void _EigenFields_Destroy( void* eigenFields, void* data ) ;
   void _EigenFields_Init(
            void*             eigenFields,
            DomainContext*    context,
            FeVariable*       tensorField,
            FeVariable**      eigenValueField,
            FeVariable**       eigenVectorField ) ;

   void _EigenFields_GetValueAtNode( void* feVariable, Node_DomainIndex dNode_I, double* value );
   InterpolationResult _EigenFields_InterpolateValueAt( void* variable, double* coord, double* value );   
   void _EigenFields_UpdateFields( void* _self );
   void _EigenFields_SyncFields( void* _self, void* context );
#endif
