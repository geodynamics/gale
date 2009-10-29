

#ifndef __Underworld_DepthDependentViscosity_h__
#define __Underworld_DepthDependentViscosity_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type DepthDependentViscosity_Type;

	/** Rheology class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __DepthDependentViscosity \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__Rheology \
		/* Virtual functions go here */ \
		/* Material Parameters */\
		FeMesh*    		                            feMesh;                             \
		double                                              eta0;                               \
		double                                              gamma;                              \
		Axis                                                variationAxis;                      \
		double                                              referencePoint;                     

	struct DepthDependentViscosity { __DepthDependentViscosity };
	
	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	DepthDependentViscosity* _DepthDependentViscosity_New( 
		SizeT                                              sizeOfSelf,
		Type                                               type,
		Stg_Class_DeleteFunction*                          _delete,
		Stg_Class_PrintFunction*                           _print,
		Stg_Class_CopyFunction*                            _copy, 
		Stg_Component_DefaultConstructorFunction*          _defaultConstructor,
		Stg_Component_ConstructFunction*                   _construct,
		Stg_Component_BuildFunction*                       _build,
		Stg_Component_InitialiseFunction*                  _initialise,
		Stg_Component_ExecuteFunction*                     _execute,
		Stg_Component_DestroyFunction*                     _destroy,
		Rheology_ModifyConstitutiveMatrixFunction*         _modifyConstitutiveMatrix,
		Name                                               name );

	
	/* 'Stg_Component' implementations */
	void* _DepthDependentViscosity_DefaultNew( Name name ) ;
	void _DepthDependentViscosity_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data );

	void _DepthDependentViscosity_ModifyConstitutiveMatrix( 
		void*                                              rheology, 
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               swarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi );
#endif
