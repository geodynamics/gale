/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Monash Cluster Computing, Australia
** (C) 2003-2004 All Rights Reserved
**
** Primary Authors:
** Luke Hodkinson
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Underworld_Utils_SmoothVelGradField_h__
#define __Underworld_Utils_SmoothVelGradField_h__

	/*
	** Class Definition.
	*/

	extern const Type SmoothVelGradField_Type;

	#define __SmoothVelGradField \
   	__ParticleFeVariable \
		\
   	/* Virtual functions. */ \
		\
   	/* Members. */ \
   	Variable_Register*	variable_Register; \
   	Variable*				dataVariableList[9]; \
   	FeVariable*				velField;

	struct SmoothVelGradField { __SmoothVelGradField };

	/*
	** Constructors/Destructors.
	*/

	SmoothVelGradField* _SmoothVelGradField_New(
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
      FieldVariable_GetValueFunction*                    _getMinGlobalFeMagnitude,
      FieldVariable_GetValueFunction*                   _getMaxGlobalFeMagnitude,
      FieldVariable_GetCoordFunction*                   _getMinAndMaxLocalCoords,
      FieldVariable_GetCoordFunction*                   _getMinAndMaxGlobalCoords,
      FeVariable_InterpolateWithinElementFunction*      _interpolateWithinElement,
      FeVariable_GetValueAtNodeFunction*                _getValueAtNode,
      ParticleFeVariable_ValueAtParticleFunction*       _valueAtParticle,
      Name                                              name );

	void _SmoothVelGradField_Init(
		SmoothVelGradField* self,
		Variable_Register* variable_Register,
		FeVariable* velField,
		SystemLinearEquations* sle );

	void* _SmoothVelGradField_DefaultNew( Name name );

	void _SmoothVelGradField_Delete( void* _self );

	void SmoothVelGradField_NonLinearUpdate( void* _sle, void* _ctx );

	/*
	** Methods.
	*/

	void _SmoothVelGradField_Print( void* _self, Stream* stream );

	void* _SmoothVelGradField_Copy( void* _self, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

	void _SmoothVelGradField_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data );

	void _SmoothVelGradField_Build( void* _self, void* data );

	void _SmoothVelGradField_Initialise( void* _self, void* data );

	void _SmoothVelGradField_Execute( void* _self, void* data );

	void _SmoothVelGradField_Destroy( void* _self, void* data );

	void _SmoothVelGradField_ValueAtParticle(
		void*							_self,
		IntegrationPointsSwarm*	swarm,
		Element_LocalIndex		lElement_I,
		void*							particle,
		double*						pressure );

#endif
