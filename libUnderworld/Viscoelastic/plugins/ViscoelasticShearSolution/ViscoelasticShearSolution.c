
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <Viscoelastic/Base/Base.h>

#include <assert.h>

typedef struct {
   __Codelet
	MaterialViscosity*               materialViscosity;
	ViscoelasticRheology*            viscoelasticRheology;
	FeMesh*				 feMesh;
	FeVariable* velocityField;
	FeVariable* strainRateField;
	FeVariable* stressField;
	double                           velocityTopOfBox;
	double                           maxTime;
	Bool                             solidBodyRotation;
	double                           sbrCentreX;
	double                           sbrCentreY;
	double                           sbrCentreZ;
	double                           omega;	
} ViscoelasticShearSolution;

const Type ViscoelasticShearSolution_Type = "Viscoelastic_ViscoelasticShearSolution";

void ViscoelasticShearSolution_VelocityBC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _data, void* _result ) {
	FiniteElementContext *	         context       = (FiniteElementContext*)_context;
	ViscoelasticShearSolution*   self          = (ViscoelasticShearSolution*) LiveComponentRegister_Get( context->CF->LCRegister, (Name)ViscoelasticShearSolution_Type );
	double*                          result        = (double*) _result;
	double                           shear;
	
	if ( context->dim == 3 ) 
		shear = sqrt( self->velocityTopOfBox * self->velocityTopOfBox / 2.0 );
	else
		shear = self->velocityTopOfBox;

	/* Stop shear once we have reached the 'maxTime' */
	if (context->currentTime >= self->maxTime) 
		shear = 0.0;
		
	*result = shear;
}

void ViscoelasticShearSolution_VelocityBCX( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _data, void* _result ) {
	FiniteElementContext *	            context            = (FiniteElementContext*)_context;
	ViscoelasticShearSolution*          self               = (ViscoelasticShearSolution*) LiveComponentRegister_Get( context->CF->LCRegister, ViscoelasticShearSolution_Type );
	FeVariable*                         feVariable         = NULL;
	FeMesh*                             mesh               = NULL;

	double*                             coord;
	double*                             result        = (double*) _result;
	double                              shearVelocity;
	double                              rotationVelocity = 0.0;
	double                              shearVelocityPrime[3], shearVelocityUnprimed[3];
	double                              rotation[3][3], rotationAxis[3];
	double                              theta;
	int                                 i,j;
	
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;	
	
	coord = Mesh_GetVertex( mesh, node_lI );
	
	/* Stop shear once we have reached the 'maxTime' */
	if (context->currentTime >= self->maxTime) 
		shearVelocity = 0.0;
	else {

	shearVelocityPrime[0] = self->velocityTopOfBox ;		shearVelocityPrime[1] = 0.0;		shearVelocityPrime[2] = 0.0 ;
	rotationAxis[0] = 0.0 ; 								rotationAxis[1] = 1.0 ;  			rotationAxis[2] = 0.0 ; 

	theta = -1.0 * self->omega * context->currentTime  ;
		
	rotation[0][0] = rotationAxis[0]*rotationAxis[0]*(1-cos(theta))+cos(theta);
	rotation[0][1] = rotationAxis[0]*rotationAxis[1]*(1-cos(theta))-rotationAxis[2]*sin(theta);
	rotation[0][2] = rotationAxis[0]*rotationAxis[2]*(1-cos(theta))+rotationAxis[1]*sin(theta);
	
	rotation[1][0] = rotationAxis[0]*rotationAxis[1]*(1-cos(theta))+rotationAxis[2]*sin(theta);
	rotation[1][1] = rotationAxis[1]*rotationAxis[1]*(1-cos(theta))+cos(theta);
	rotation[1][2] = rotationAxis[1]*rotationAxis[2]*(1-cos(theta))-rotationAxis[0]*sin(theta);
	
	rotation[2][0] = rotationAxis[0]*rotationAxis[2]*(1-cos(theta))-rotationAxis[1]*sin(theta);
	rotation[2][1] = rotationAxis[1]*rotationAxis[2]*(1-cos(theta))+rotationAxis[0]*sin(theta);
	rotation[2][2] = rotationAxis[2]*rotationAxis[2]*(1-cos(theta))+cos(theta);

	
	for ( i = 0 ; i < 3 ; i++ ) {
		shearVelocityUnprimed[i] = 0.0;
		for ( j = 0 ; j < 3 ; j++ ) {
				shearVelocityUnprimed[i] += rotation[i][j]*shearVelocityPrime[j] ;
			}
		}
	}

	shearVelocity = shearVelocityUnprimed[0] * coord[ J_AXIS ] ;  // TODO need to have this value rotated into the frame
	
	rotationVelocity = -1.0 * self->omega * ( coord[ K_AXIS ] - self->sbrCentreZ );

	*result = shearVelocity + rotationVelocity;
}

void ViscoelasticShearSolution_VelocityBCZ( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _data, void* _result ) {
	FiniteElementContext *	            context            = (FiniteElementContext*)_context;
	ViscoelasticShearSolution*          self               = (ViscoelasticShearSolution*) LiveComponentRegister_Get( context->CF->LCRegister, ViscoelasticShearSolution_Type );
	FeVariable*                         feVariable         = NULL;
	FeMesh*                             mesh               = NULL;

	double*                             coord;
	double*                             result        = (double*) _result;
	double                              shearVelocity;
	double                              rotationVelocity = 0.0;
	double                              shearVelocityPrime[3], shearVelocityUnprimed[3];
	double                              rotation[3][3], rotationAxis[3];
	double                              theta;
	int                                 i,j;
	
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;	
	
	coord = Mesh_GetVertex( mesh, node_lI );
	
	/* Stop shear once we have reached the 'maxTime' */
	if (context->currentTime >= self->maxTime) 
		shearVelocity = 0.0;
	else {

	shearVelocityPrime[0] = self->velocityTopOfBox ;		shearVelocityPrime[1] = 0.0;		shearVelocityPrime[2] = 0.0 ;
	rotationAxis[0] = 0.0 ; rotationAxis[1] = 1.0 ;  rotationAxis[2] = 0.0 ; 

	theta = -1.0 * self->omega * context->currentTime  ;
		
	rotation[0][0] = rotationAxis[0]*rotationAxis[0]*(1-cos(theta))+cos(theta);
	rotation[0][1] = rotationAxis[0]*rotationAxis[1]*(1-cos(theta))-rotationAxis[2]*sin(theta);
	rotation[0][2] = rotationAxis[0]*rotationAxis[2]*(1-cos(theta))+rotationAxis[1]*sin(theta);
	
	rotation[1][0] = rotationAxis[0]*rotationAxis[1]*(1-cos(theta))+rotationAxis[2]*sin(theta);
	rotation[1][1] = rotationAxis[1]*rotationAxis[1]*(1-cos(theta))+cos(theta);
	rotation[1][2] = rotationAxis[1]*rotationAxis[2]*(1-cos(theta))-rotationAxis[0]*sin(theta);
	
	rotation[2][0] = rotationAxis[0]*rotationAxis[2]*(1-cos(theta))-rotationAxis[1]*sin(theta);
	rotation[2][1] = rotationAxis[1]*rotationAxis[2]*(1-cos(theta))+rotationAxis[0]*sin(theta);
	rotation[2][2] = rotationAxis[2]*rotationAxis[2]*(1-cos(theta))+cos(theta);

	
	for ( i = 0 ; i < 3 ; i++ ) {
		shearVelocityUnprimed[i] = 0.0;
		for ( j = 0 ; j < 3 ; j++ ) {
				shearVelocityUnprimed[i] += rotation[i][j]*shearVelocityPrime[j] ;
			}
		}
	}

	shearVelocity = shearVelocityUnprimed[2] * coord[ J_AXIS ] ;  // TODO need to have this value rotated into the frame
	
	rotationVelocity = self->omega * ( coord[ I_AXIS ] - self->sbrCentreX );

	*result = shearVelocity + rotationVelocity;
}

double ViscoelasticShearSolution_ImposeTimestep( void* self, UnderworldContext* context ) {
	double prescribedTime = Dictionary_GetDouble_WithDefault( context->dictionary, (Dictionary_Entry_Key)"prescribedTime", -1.0  );
	return prescribedTime;
}

void _ViscoelasticShearSolution_AssignFromXML( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	ViscoelasticShearSolution* self = (ViscoelasticShearSolution*)analyticSolution;
	FiniteElementContext*          context;

   	Dictionary*	pluginDict	= Codelet_GetPluginDictionary( self, cf->rootDict );

   	self->context = Stg_ComponentFactory_ConstructByName( cf, Dictionary_GetString( pluginDict, (Dictionary_Entry_Key)"Context"  ), AbstractContext, True, data );
   	if( !self->context )
	   	self->context = Stg_ComponentFactory_ConstructByName( cf, (Name)"context", AbstractContext, True, data ); 
	
	context = Stg_CheckType( self->context, FiniteElementContext );
	
	self->maxTime 			= Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"maxTime", HUGE_VAL );
	self->velocityTopOfBox 	= Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"velocityTopOfBox", 0.0 );
	self->solidBodyRotation = Dictionary_GetBool_WithDefault( context->dictionary, (Dictionary_Entry_Key)"solidBodyRotation", False );
	self->sbrCentreX 		= Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"SolidBodyRotationCentreX", 0.0 );
	self->sbrCentreY 		= Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"SolidBodyRotationCentreY", 0.0 );
	self->sbrCentreZ 		= Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"SolidBodyRotationCentreZ", 0.0 );
	self->omega      		= Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"SolidBodyRotationOmega", 1.0 );
	
	ConditionFunction_Register_Add( condFunc_Register, 
			ConditionFunction_New( ViscoelasticShearSolution_VelocityBC, (Name)"ShearTrigger", NULL) );	
			
	ConditionFunction_Register_Add( condFunc_Register, 
			ConditionFunction_New( ViscoelasticShearSolution_VelocityBCX, (Name)"ShearTriggerX", NULL) );	
			
	ConditionFunction_Register_Add( condFunc_Register, 
			ConditionFunction_New( ViscoelasticShearSolution_VelocityBCZ, (Name)"ShearTriggerZ", NULL) );	

	EP_ReplaceAllClassHook( context->calcDtEP, ViscoelasticShearSolution_ImposeTimestep, context ); 
}

/* This function will provide StGermain the abilty to instantiate (create) this codelet on demand. */
void* _ViscoelasticShearSolution_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(ViscoelasticShearSolution);
	Type                                                      type = ViscoelasticShearSolution_Type;
	Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
	Stg_Class_PrintFunction*                                _print = _Codelet_Print;
	Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _ViscoelasticShearSolution_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _ViscoelasticShearSolution_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Codelet_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _Codelet_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _Codelet_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return _Codelet_New(  CODELET_PASSARGS  );
}
	
/* This function is automatically run by StGermain when this plugin is loaded. The name must be "<plugin-name>_Register". */
Index Viscoelastic_ViscoelasticShearSolution_Register( PluginsManager* pluginsManager ) {
	/* A plugin is only properly registered once it returns the handle provided when submitting a codelet to StGermain. */
	return PluginsManager_Submit( pluginsManager, ViscoelasticShearSolution_Type, (Name)"0", _ViscoelasticShearSolution_DefaultNew  );
}
