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

MappedParticleLayout* _MappedParticleLayout_New( 
      SizeT                                            _sizeOfSelf,
      Type                                             type,
      Stg_Class_DeleteFunction*                        _delete,
      Stg_Class_PrintFunction*                         _print,
      Stg_Class_CopyFunction*                          _copy, 
      Stg_Component_DefaultConstructorFunction*        _defaultConstructor,
      Stg_Component_ConstructFunction*                 _construct,
      Stg_Component_BuildFunction*                     _build,
      Stg_Component_InitialiseFunction*                _initialise,
      Stg_Component_ExecuteFunction*                   _execute,
      Stg_Component_DestroyFunction*                   _destroy,
      Name                                             name,
      AllocationType                                   nameAllocationType,
      ParticleLayout_SetInitialCountsFunction*         _setInitialCounts,
      ParticleLayout_InitialiseParticlesFunction*      _initialiseParticles,
      CoordSystem                                      coordSystem,
      Bool                                             weightsInitialisedAtStartup );
	
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
