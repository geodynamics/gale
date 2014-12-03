

#ifndef __Underworld_Utils_RadiogenicHeatingTerm_h__
#define __Underworld_Utils_RadiogenicHeatingTerm_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type RadiogenicHeatingTerm_Type;

	typedef struct { 
		double	Q; /* value of heating coefficient */
		double	lambda; /* decay rate for this element */
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
		Materials_Register*	materials_Register; \
		ExtensionInfo_Index	materialExtHandle;

	struct RadiogenicHeatingTerm { __RadiogenicHeatingTerm };

	RadiogenicHeatingTerm* RadiogenicHeatingTerm_New( 
		Name							name,
		FiniteElementContext*	context,
		ForceVector*				forceVector,
		Swarm*						integrationSwarm,
		Materials_Register*		materials_Register );

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define RADIOGENICHEATINGTERM_DEFARGS \
                FORCETERM_DEFARGS

	#define RADIOGENICHEATINGTERM_PASSARGS \
                FORCETERM_PASSARGS

	RadiogenicHeatingTerm* _RadiogenicHeatingTerm_New(  RADIOGENICHEATINGTERM_DEFARGS  );

	void _RadiogenicHeatingTerm_Init( void* forceTerm, Materials_Register* materials_Register );
	
	void _RadiogenicHeatingTerm_Delete( void* forceTerm );

	void _RadiogenicHeatingTerm_Print( void* forceTerm, Stream* stream );

	void* _RadiogenicHeatingTerm_DefaultNew( Name name ) ;

	void _RadiogenicHeatingTerm_AssignFromXML( void* forceTerm, Stg_ComponentFactory* cf, void* data ) ;

	void _RadiogenicHeatingTerm_Build( void* forceTerm, void* data ) ;

	void _RadiogenicHeatingTerm_Initialise( void* forceTerm, void* data ) ;

	void _RadiogenicHeatingTerm_Execute( void* forceTerm, void* data ) ;

	void _RadiogenicHeatingTerm_Destroy( void* forceTerm, void* data ) ;

	void _RadiogenicHeatingTerm_AssembleElement( void* forceTerm, ForceVector* forceVector, Element_LocalIndex lElement_I, double* elForceVec ) ;

#endif

