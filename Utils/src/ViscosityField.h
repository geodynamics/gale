/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Monash Cluster Computing, Australia
** (C) 2003-2004 All Rights Reserved
**
** Primary Authors:
** Robert Turnbull, MCC
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Underworld_Utils_ViscosityField_h__
#define __Underworld_Utils_ViscosityField_h__
	
	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type ViscosityField_Type;

	/** ViscosityField class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __ViscosityField \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__ParticleFeVariable \
		\
		/* Virtual functions go here */ \
		\
		/* Passed in parameters */ \
		Variable_Register*	variable_Register; \
		ConstitutiveMatrix*	constitutiveMatrix;
		
	struct ViscosityField { __ViscosityField };
	
	/* --- Contstructors / Destructors --- */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define VISCOSITYFIELD_DEFARGS \
                PARTICLEFEVARIABLE_DEFARGS

	#define VISCOSITYFIELD_PASSARGS \
                PARTICLEFEVARIABLE_PASSARGS

	ViscosityField* _ViscosityField_New(  VISCOSITYFIELD_DEFARGS  );
	
	/** Print the contents of an ViscosityField construct */
	void* _ViscosityField_DefaultNew( Name name );

	void _ViscosityField_Delete( void* variable );

	void _ViscosityField_Print( void* variable, Stream* stream );

	void* _ViscosityField_Copy( void* feVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

	void _ViscosityField_AssignFromXML( void* variable, Stg_ComponentFactory* cf, void* data );

	void _ViscosityField_Build( void* variable, void* data );

	void _ViscosityField_Initialise( void* variable, void* data );

	void _ViscosityField_Execute( void* variable, void* data );

	void _ViscosityField_Destroy( void* variable, void* data );

	void _ViscosityField_ValueAtParticle( void* viscosityField, IntegrationPointsSwarm* swarm, Element_LocalIndex lElement_I, void* particle, double* viscosity );
	
#endif 

