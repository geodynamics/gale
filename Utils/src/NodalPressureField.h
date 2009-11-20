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

	NodalPressureField* _NodalPressureField_New( 
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
