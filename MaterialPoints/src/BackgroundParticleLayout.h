#ifndef __PICellerator_MaterialPoints_BackgroundParticleLayout_h__
#define __PICellerator_MaterialPoints_BackgroundParticleLayout_h__
	
/* Textual name of this class */
extern const Type BackgroundParticleLayout_Type;
	
/* ParticleLayout information */
#define __BackgroundParticleLayout                              \
    /* General info */                                          \
    __ParticleLayout                                            \

struct BackgroundParticleLayout { __BackgroundParticleLayout };
	
/** public constructor */
BackgroundParticleLayout* BackgroundParticleLayout_New( Name name,
   AbstractContext* context, 
   CoordSystem      coordSystem, 
   Bool             weightsInitialisedAtStartup );

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define BACKGROUNDPARTICLELAYOUT_DEFARGS \
                PARTICLELAYOUT_DEFARGS

	#define BACKGROUNDPARTICLELAYOUT_PASSARGS \
                PARTICLELAYOUT_PASSARGS

BackgroundParticleLayout* _BackgroundParticleLayout_New(  BACKGROUNDPARTICLELAYOUT_DEFARGS  );
	
/* Initialise implementation */
void _BackgroundParticleLayout_Init( 
    void*                  particleLayout );
	
void _BackgroundParticleLayout_Delete( void* particleLayout );
void _BackgroundParticleLayout_Print( void* particleLayout, Stream* stream );
#define BackgroundParticleLayout_Copy( self )                           \
    (BackgroundParticleLayout*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
#define BackgroundParticleLayout_DeepCopy( self )                       \
    (BackgroundParticleLayout*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
void* _BackgroundParticleLayout_Copy( void* particleLayout, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

void* _BackgroundParticleLayout_DefaultNew( Name name );
void  _BackgroundParticleLayout_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data );
void  _BackgroundParticleLayout_Build( void* component, void* data );
void  _BackgroundParticleLayout_Initialise( void* component, void* data );
void  _BackgroundParticleLayout_Execute( void* component, void* data );
void  _BackgroundParticleLayout_Destroy( void* component, void* data );
	
void _BackgroundParticleLayout_SetInitialCounts( void* particleLayout, void* _swarm );
void _BackgroundParticleLayout_InitialiseParticles( void* particleLayout, void* _swarm );

	
#endif /* __PICellerator_MaterialPoints_BackgroundParticleLayout_h__ */

