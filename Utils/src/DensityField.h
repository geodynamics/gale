/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Monash Cluster Computing, Australia
** (C) 2003-2004 All Rights Reserved
**
** Primary Authors:
** Robert Turnbull, MCC
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Underworld_Utils_DensityField_h__
#define __Underworld_Utils_DensityField_h__
	
	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type DensityField_Type;

	/** DensityField class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __DensityField \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__ParticleFeVariable \
		\
		/* Virtual functions go here */ \
		\
		/* Passed in parameters */ \
		BuoyancyForceTerm*       buoyancyForceTerm; \
		Variable_Register*       variable_Register;
		
	struct DensityField { __DensityField };
	
	/* --- Contstructors / Destructors --- */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define DENSITYFIELD_DEFARGS \
                PARTICLEFEVARIABLE_DEFARGS

	#define DENSITYFIELD_PASSARGS \
                PARTICLEFEVARIABLE_PASSARGS

	DensityField* _DensityField_New(  DENSITYFIELD_DEFARGS  );
	
	/** Print the contents of an DensityField construct */
	void* _DensityField_DefaultNew( Name name );
	void _DensityField_Delete( void* variable );
	void _DensityField_Print( void* variable, Stream* stream );
	void* _DensityField_Copy( void* feVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

	void _DensityField_AssignFromXML( void* variable, Stg_ComponentFactory* cf, void* data ) ;
	void _DensityField_Build( void* variable, void* data ) ;
	void _DensityField_Initialise( void* variable, void* data ) ;
	void _DensityField_Execute( void* variable, void* data ) ;
	void _DensityField_Destroy( void* variable, void* data ) ;
	void _DensityField_ValueAtParticle( void* viscosityField, IntegrationPointsSwarm* swarm, Element_LocalIndex lElement_I, void* particle, double* viscosity ) ;
	

#endif 

