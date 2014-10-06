#ifndef __PICellerator_Common_Ppc_PointGravity_h__
#define __PICellerator_Common_Ppc_PointGravity_h__

/** Textual name of this class */
extern const Type Ppc_PointGravity_Type;

/** Ppc_PointGravity class contents */
#define __Ppc_PointGravity									 \
  /* Parent info */									 \
  __Ppc													 \
  /* General data */									 \
  double          point[3];				 \
  int             alphaTag;						 \

	struct Ppc_PointGravity { __Ppc_PointGravity };

	Ppc_PointGravity* Ppc_PointGravity_New( Name name, Stg_Component* _self );
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define PPC_POINTGRAVITY_DEFARGS \
                PPC_DEFARGS

	#define PPC_POINTGRAVITY_PASSARGS \
                PPC_PASSARGS

	Ppc_PointGravity* _Ppc_PointGravity_New(  PPC_POINTGRAVITY_DEFARGS  );
	
	void _Ppc_PointGravity_Delete( void* _self );
	void _Ppc_PointGravity_Print( void* _self, Stream* stream );
	void* _Ppc_PointGravity_DefaultNew( Name name );
   void _Ppc_PointGravity_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data );
	void _Ppc_PointGravity_Build( void* _self, void* data );
	void _Ppc_PointGravity_Initialise( void* _self, void* data );
	void _Ppc_PointGravity_Execute( void* _self, void* data );
	void _Ppc_PointGravity_Destroy( void* _self, void* data );

   /* Public functions */
   int _Ppc_PointGravity_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result );

#endif

