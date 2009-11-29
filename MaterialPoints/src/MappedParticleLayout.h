#ifndef __PICellerator_MaterialPoints_MappedParticleLayout_h__
#define __PICellerator_MaterialPoints_MappedParticleLayout_h__
	
/* Textual name of this class */
extern const Type MappedParticleLayout_Type;
	
/* ParticleLayout information */
#define __MappedParticleLayout                          \
    /* General info */                                  \
    __ParticleLayout                                    \

struct MappedParticleLayout { __MappedParticleLayout };
	
MappedParticleLayout* MappedParticleLayout_New( 
      Name name, 
      AbstractContext* context,
      CoordSystem      coordSystem,
      Bool             weightsInitialisedAtStartup);

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define MAPPEDPARTICLELAYOUT_DEFARGS \
                PARTICLELAYOUT_DEFARGS

	#define MAPPEDPARTICLELAYOUT_PASSARGS \
                PARTICLELAYOUT_PASSARGS

MappedParticleLayout* _MappedParticleLayout_New(  MAPPEDPARTICLELAYOUT_DEFARGS  );
	
/* Initialise implementation */
void _MappedParticleLayout_Init( 
    void*                  particleLayout );
	
void _MappedParticleLayout_Delete( void* particleLayout );
void _MappedParticleLayout_Print( void* particleLayout, Stream* stream );
#define MappedParticleLayout_Copy( self )                               \
    (MappedParticleLayout*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
#define MappedParticleLayout_DeepCopy( self )                           \
    (MappedParticleLayout*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
void* _MappedParticleLayout_Copy( void* particleLayout, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

void* _MappedParticleLayout_DefaultNew( Name name );
void  _MappedParticleLayout_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data );
void  _MappedParticleLayout_Build( void* component, void* data );
void  _MappedParticleLayout_Initialise( void* component, void* data );
void  _MappedParticleLayout_Execute( void* component, void* data );
void  _MappedParticleLayout_Destroy( void* component, void* data );
	
void _MappedParticleLayout_SetInitialCounts( void* particleLayout, void* _swarm );
void _MappedParticleLayout_InitialiseParticles( void* particleLayout, void* _swarm );

	
#endif /* __PICellerator_MaterialPoints_MappedParticleLayout_h__ */

