

#ifndef __Underworld_Utils_PorosityTerm_h__
#define __Underworld_Utils_PorosityTerm_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	typedef struct PorosityTerm PorosityTerm;

	extern const Type PorosityTerm_Type;

	typedef struct { 
		double          phi;                   /* value of porosity coefficient */
	} PorosityExt;

	/*} PorosityElement;

	typedef struct {
		PorosityElement* porosityElementList;
		Index            porosityElementCount;
	} PorosityTerm_ParticleExt;*/

	/** PorosityTerm class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __PorosityTerm \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__ForceTerm \
		\
		/* Virtual functions go here */ \
		\
		/* PorosityTerm info */ \
		AbstractContext*		context;   \
		/*ExtensionManager_Register*    particles_Register;*/   \
		ExtensionInfo_Index             particleExtHandle;   \
		double                    	referencePorosity;   \
		double                    	gravity;   \
		MaterialPointsSwarm*		materialSwarm;   \

	struct PorosityTerm { __PorosityTerm };

	PorosityTerm* PorosityTerm_New( 
		Name                                                name,
		ForceVector*                                        forceVector,
		Swarm*                                              integrationSwarm,
		AbstractContext*                                    context,
		Materials_Register*                                 materials_Register );
		/*ExtensionManager_Register*                          particle_Register );*/

	PorosityTerm* _PorosityTerm_New( 
		SizeT                                               sizeOfSelf,  
		Type                                                type,
		Stg_Class_DeleteFunction*                           _delete,
		Stg_Class_PrintFunction*                            _print,
		Stg_Class_CopyFunction*                             _copy, 
		Stg_Component_DefaultConstructorFunction*           _defaultConstructor,
		Stg_Component_ConstructFunction*                    _construct,
		Stg_Component_BuildFunction*                        _build,
		Stg_Component_InitialiseFunction*                   _initialise,
		Stg_Component_ExecuteFunction*                      _execute,
		Stg_Component_DestroyFunction*                      _destroy,
		ForceTerm_AssembleElementFunction*                  _assembleElement,		
		Name                                                name );
	
	void PorosityTerm_InitAll( 
		void*                                               buoyancyForceTerm,
		ForceVector*                                        forceVector,
		Swarm*                                              integrationSwarm,
		AbstractContext*                                    context,
		Materials_Register*                                 materials_Register );
		/*ExtensionManager_Register*                          particle_Register );*/

	void _PorosityTerm_Delete( void* buoyancyForceTerm );
	void _PorosityTerm_Print( void* buoyancyForceTerm, Stream* stream );

	void* _PorosityTerm_DefaultNew( Name name ) ;
	void _PorosityTerm_Construct( void* buoyancyForceTerm, Stg_ComponentFactory* cf, void* data ) ;
	void _PorosityTerm_Build( void* buoyancyForceTerm, void* data ) ;
	void _PorosityTerm_Initialise( void* buoyancyForceTerm, void* data ) ;
	void _PorosityTerm_Execute( void* buoyancyForceTerm, void* data ) ;
	void _PorosityTerm_Destroy( void* buoyancyForceTerm, void* data ) ;

	void _PorosityTerm_AssembleElement( void* buoyancyForceTerm, ForceVector* forceVector, Element_LocalIndex lElement_I, double* elForceVec ) ;

#endif
