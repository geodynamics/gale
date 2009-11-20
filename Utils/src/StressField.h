/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Monash Cluster Computing, Australia
** (C) 2003-2004 All Rights Reserved
**
** Primary Authors:
** Robert Turnbull, MCC
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Underworld_MaterialPoints_StressField_h__
#define __Underworld_MaterialPoints_StressField_h__
	
	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type StressField_Type;

	/** StressField class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __StressField \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__ParticleFeVariable \
		\
		/* Virtual functions go here */ \
		\
		/* Passed in parameters */ \
		Variable_Register*	variable_Register; \
		FeVariable*				strainRateField; \
		ConstitutiveMatrix*	constitutiveMatrix; \
		Variable*				dataVariableList[6]; \
		Variable*				stressVariable;
		
	struct StressField { __StressField };
	
	/* --- Contstructors / Destructors --- */
	StressField* _StressField_New(
 		SizeT                                             _sizeOfSelf,
		Type                                              type,
		Stg_Class_DeleteFunction*                         _delete,
		Stg_Class_PrintFunction*                          _print,
		Stg_Class_CopyFunction*                           _copy, 
		Stg_Component_DefaultConstructorFunction*         _defaultConstructor,
		Stg_Component_ConstructFunction*                  _construct,
		Stg_Component_BuildFunction*                      _build,
		Stg_Component_InitialiseFunction*                 _initialise,
		Stg_Component_ExecuteFunction*                    _execute,
		Stg_Component_DestroyFunction*                    _destroy,
		FieldVariable_InterpolateValueAtFunction*         _interpolateValueAt,
		FieldVariable_GetValueFunction*	                  _getMinGlobalFeMagnitude,
		FieldVariable_GetValueFunction*                   _getMaxGlobalFeMagnitude,
		FieldVariable_GetCoordFunction*                   _getMinAndMaxLocalCoords,
		FieldVariable_GetCoordFunction*                   _getMinAndMaxGlobalCoords,		
		FeVariable_InterpolateWithinElementFunction*      _interpolateWithinElement,	
		FeVariable_GetValueAtNodeFunction*                _getValueAtNode,
		ParticleFeVariable_ValueAtParticleFunction*       _valueAtParticle,
		Name                                              name );
	
	/** Print the contents of an StressField construct */
	void* _StressField_DefaultNew( Name name );

	void _StressField_Delete( void* variable );

	void _StressField_Print( void* variable, Stream* stream );

	void* _StressField_Copy( void* feVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

	void _StressField_AssignFromXML( void* variable, Stg_ComponentFactory* cf, void* data );

	void _StressField_Init(	StressField* self, FeVariable* strainRateField,	ConstitutiveMatrix* constitutiveMatrix, Variable* stressVariable, Variable_Register* variable_Register, SystemLinearEquations* sle);

   void StressField_NonLinearUpdate( void* _sle, void* _ctx );

   void _StressField_Build( void* variable, void* data );

	void _StressField_Initialise( void* variable, void* data );

	void _StressField_Execute( void* variable, void* data );

	void _StressField_Destroy( void* variable, void* data );

	void _StressField_ValueAtParticle_Recalculate( void* stressField, IntegrationPointsSwarm* swarm, Element_LocalIndex lElement_I, void* particle, double* stress );

	void _StressField_ValueAtParticle_FromVariable( void* stressField, IntegrationPointsSwarm* swarm, Element_LocalIndex lElement_I, void* particle, double* stress );

#endif 
