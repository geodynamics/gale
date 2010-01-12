/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Monash Cluster Computing, Australia
** (C) 2003-2004 All Rights Reserved
**
** Primary Authors:
** Luke Hodkinson
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Underworld_Utils_NodalPressureField_h__
#define __Underworld_Utils_NodalPressureField_h__

	extern const Type NodalPressureField_Type;

	#define __NodalPressureField \
	   __ParticleFeVariable \
		\
	   /* Virtual functions. */ \
		\
	   /* Members. */ \
	   Variable_Register* variable_Register; \
	   FeVariable* pressureField;

	struct NodalPressureField { __NodalPressureField };

	/*
	** Constructors/Destructors.
	*/

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define NODALPRESSUREFIELD_DEFARGS \
                PARTICLEFEVARIABLE_DEFARGS

	#define NODALPRESSUREFIELD_PASSARGS \
                PARTICLEFEVARIABLE_PASSARGS

	NodalPressureField* _NodalPressureField_New(  NODALPRESSUREFIELD_DEFARGS  );

	void _NodalPressureField_Init( NodalPressureField* self, Variable_Register* variable_Register, FeVariable* pressureField, SystemLinearEquations* sle );

	void* _NodalPressureField_DefaultNew( Name name );
	
	void _NodalPressureField_Delete( void* _self );

	/*
	** Methods.
	*/

	void _NodalPressureField_Print( void* _self, Stream* stream );

	void* _NodalPressureField_Copy( void* _self, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

	void _NodalPressureField_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data ) ;

	void _NodalPressureField_Build( void* _self, void* data ) ;

	void _NodalPressureField_Initialise( void* _self, void* data ) ;

	void _NodalPressureField_Execute( void* _self, void* data ) ;

	void _NodalPressureField_Destroy( void* _self, void* data ) ;

	void _NodalPressureField_ValueAtParticle(
		void*							_self,
		IntegrationPointsSwarm*	swarm,
		Element_LocalIndex		lElement_I,
		void*							particle,
		double*						pressure ) ;

   void NodalPressureField_NonLinearUpdate( void* _sle, void* _ctx );
#endif

