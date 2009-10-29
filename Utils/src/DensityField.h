/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Monash Cluster Computing, Australia
** (C) 2003-2004 All Rights Reserved
**
** Primary Authors:
** Robert Turnbull, MCC
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Underworld_MaterialPoints_DensityField_h__
#define __Underworld_MaterialPoints_DensityField_h__
	
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
	DensityField* _DensityField_New(
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
