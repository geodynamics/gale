

#ifndef __Underworld_Utils_RadiogenicHeatingTerm_h__
#define __Underworld_Utils_RadiogenicHeatingTerm_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type RadiogenicHeatingTerm_Type;

	typedef struct { 
		double          Q;                   /* value of heating coefficient */
		double          lambda;              /* decay rate for this element */
	} HeatingElement;

	typedef struct {
		HeatingElement* heatingElementList;
		Index           heatingElementCount;
	} RadiogenicHeatingTerm_MaterialExt;

	/** RadiogenicHeatingTerm class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __RadiogenicHeatingTerm \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__ForceTerm \
		\
		/* Virtual functions go here */ \
		\
		/* RadiogenicHeatingTerm info */ \
		AbstractContext*                                    context;                           \
		Materials_Register*                                 materials_Register;                \
		ExtensionInfo_Index                                 materialExtHandle;

	struct RadiogenicHeatingTerm { __RadiogenicHeatingTerm };

	RadiogenicHeatingTerm* RadiogenicHeatingTerm_New( 
		Name                                                name,
		ForceVector*                                        forceVector,
		Swarm*                                              integrationSwarm,
		AbstractContext*                                    context,
		Materials_Register*                                 materials_Register );

	RadiogenicHeatingTerm* _RadiogenicHeatingTerm_New( 
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
	
	void RadiogenicHeatingTerm_InitAll( 
		void*                                               forceTerm,
		ForceVector*                                        forceVector,
		Swarm*                                              integrationSwarm,
		AbstractContext*                                    context,
		Materials_Register*                                 materials_Register );

	void _RadiogenicHeatingTerm_Delete( void* forceTerm );
	void _RadiogenicHeatingTerm_Print( void* forceTerm, Stream* stream );

	void* _RadiogenicHeatingTerm_DefaultNew( Name name ) ;
void _RadiogenicHeatingTerm_Construct( void* forceTerm, Stg_ComponentFactory* cf, void* data ) ;
	void _RadiogenicHeatingTerm_Build( void* forceTerm, void* data ) ;
	void _RadiogenicHeatingTerm_Initialise( void* forceTerm, void* data ) ;
	void _RadiogenicHeatingTerm_Execute( void* forceTerm, void* data ) ;
	void _RadiogenicHeatingTerm_Destroy( void* forceTerm, void* data ) ;

	void _RadiogenicHeatingTerm_AssembleElement( void* forceTerm, ForceVector* forceVector, Element_LocalIndex lElement_I, double* elForceVec ) ;

#endif
