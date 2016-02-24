#ifndef __Underworld_Rheology_HKViscousCreep_h__
#define __Underworld_Rheology_HKViscousCreep_h__

    /** Textual name of this class */
    extern const Type HKViscousCreep_Type;
        
    /** Rheology class contents */
    #define __HKViscousCreep \
        /* Macro defining parent goes here */ \
        __Rheology \
        /* Virtual functions go here */ \
        /* Material Parameters */\
        int                                               strainRateInvTag; \
        int                                               temperatureLabel; \
        int                                               pressureLabel; \
        double                                              preExponentialFactor; \
        double                                              stressExponent; \
        int                                               grainSizeLabel; \
        double                                              grainSizeExponent; \
        double                                              waterFugacityExponent; \
        int                                               waterFugacityLabel; \
        int                                               meltFractionLabel; \
        double                                              alpha; \
        double                                              activationEnergy; \
        double                                              activationVolume; \
      double                                              defaultStrainRateInvariant; \

    struct HKViscousCreep { __HKViscousCreep };
    
    /** Private Constructor: This will accept all the virtual functions for this class as arguments. */
    
    #ifndef ZERO
    #define ZERO 0
    #endif

    #define HKViscousCreep_DEFARGS \
                RHEOLOGY_DEFARGS

    #define HKViscousCreep_PASSARGS \
                RHEOLOGY_PASSARGS

    HKViscousCreep* _HKViscousCreep_New(  HKViscousCreep_DEFARGS  );

    /* 'Stg_Component' implementations */
    void* _HKViscousCreep_DefaultNew( Name name ) ;
    void _HKViscousCreep_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data );
   void _HKViscousCreep_Build( void* _self, void* data );
   void _HKViscousCreep_Initialise( void* _self, void* data );
   void _HKViscousCreep_Destroy( void* _self, void* data );
   
    void _HKViscousCreep_ModifyConstitutiveMatrix( 
        void*                                              rheology, 
        ConstitutiveMatrix*                                constitutiveMatrix,
        MaterialPointsSwarm*                               swarm,
        Element_LocalIndex                                 lElement_I,
        MaterialPoint*                                     materialPoint,
        Coord                                              xi );
#endif

