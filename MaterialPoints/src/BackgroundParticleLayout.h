#ifndef __PICellerator_MaterialPoints_BackgroundParticleLayout_h__
#define __PICellerator_MaterialPoints_BackgroundParticleLayout_h__
	
/* Textual name of this class */
extern const Type BackgroundParticleLayout_Type;
	
/* ParticleLayout information */
#define __BackgroundParticleLayout                              \
    /* General info */                                          \
    __ParticleLayout                                            \

struct BackgroundParticleLayout { __BackgroundParticleLayout };
	
BackgroundParticleLayout* _BackgroundParticleLayout_New( 
    SizeT                                                        _sizeOfSelf,
    Type                                                         type,
    Stg_Class_DeleteFunction*                                    _delete,
    Stg_Class_PrintFunction*                                     _print,
    Stg_Class_CopyFunction*                                      _copy,
    Stg_Component_DefaultConstructorFunction*                    _defaultConstructor,
    Stg_Component_ConstructFunction*                             _construct,
    Stg_Component_BuildFunction*                                 _build,
    Stg_Component_InitialiseFunction*                            _initialise,
    Stg_Component_ExecuteFunction*                               _execute,
    Stg_Component_DestroyFunction*                               _destroy,
    ParticleLayout_SetInitialCountsFunction*                     _setInitialCounts,
    ParticleLayout_InitialiseParticlesFunction*                  _initialiseParticles,
    Name                                                         name,
    Bool                                                         initFlag,
    CoordSystem                                                  coordSystem,
    Bool                                                         weightsInitialisedAtStartup );
	
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
void  _BackgroundParticleLayout_Construct( void* component, Stg_ComponentFactory* cf, void* data );
void  _BackgroundParticleLayout_Build( void* component, void* data );
void  _BackgroundParticleLayout_Initialise( void* component, void* data );
void  _BackgroundParticleLayout_Execute( void* component, void* data );
void  _BackgroundParticleLayout_Destroy( void* component, void* data );
	
void _BackgroundParticleLayout_SetInitialCounts( void* particleLayout, void* _swarm );
void _BackgroundParticleLayout_InitialiseParticles( void* particleLayout, void* _swarm );

	
#endif /* __PICellerator_MaterialPoints_BackgroundParticleLayout_h__ */
