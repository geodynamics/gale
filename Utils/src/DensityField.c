/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Monash Cluster Computing, Australia
** (C) 2003-2004 All Rights Reserved
**
** Primary Authors:
** Robert Turnbull, MCC
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include <Underworld/Rheology/Rheology.h>

#include "types.h"
#include "DensityField.h"
#include <assert.h>

const Type DensityField_Type = "DensityField";

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
		Name                                              name )
{
	DensityField*		self;
	
	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(DensityField) );
	self = (DensityField*)
		_ParticleFeVariable_New(
			_sizeOfSelf, 
			type, 
			_delete,
			_print,
			_copy,
			_defaultConstructor,
			_construct,
			_build,
			_initialise,
			_execute, 
			_destroy,
			_interpolateValueAt,
			_getMinGlobalFeMagnitude, 
			_getMaxGlobalFeMagnitude,
			_getMinAndMaxLocalCoords, 
			_getMinAndMaxGlobalCoords,
			_interpolateWithinElement,
			_getValueAtNode,
			_valueAtParticle,
			name );
	
	return self;
}

void _DensityField_Init( 
		DensityField*                                   self,
		BuoyancyForceTerm*                               buoyancyForceTerm,
		Variable_Register*                                variable_Register )
{

	/* Assign Pointers */
	self->buoyancyForceTerm = buoyancyForceTerm;
	self->variable_Register = variable_Register;
	
	/* Set pointers to swarm to be the same as the one on the constitutive matrix */
	self->assemblyTerm->integrationSwarm = self->buoyancyForceTerm->integrationSwarm;
	self->massMatrixForceTerm->integrationSwarm = self->buoyancyForceTerm->integrationSwarm;	
}

/* --- Virtual Function Implementations --- */
void _DensityField_Delete( void* densityField ) {
	DensityField* self = (DensityField*) densityField;

	_ParticleFeVariable_Delete( self );
}

void _DensityField_Print( void* densityField, Stream* stream ) {
	DensityField* self = (DensityField*) densityField;
	
	/* General info */
	Journal_Printf( stream, "DensityField (ptr): %p\n", self );
	
	/* Print parent */
	_FeVariable_Print( self, stream );
	
	/* DensityField info */
	Journal_PrintPointer( stream, self->buoyancyForceTerm );
}


void* _DensityField_Copy( void* feVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	DensityField*	self = (DensityField*)feVariable;
	DensityField*	newDensityField;
	
	newDensityField = (DensityField*) _FeVariable_Copy( feVariable, dest, deep, nameExt, ptrMap );

	newDensityField->buoyancyForceTerm = self->buoyancyForceTerm;
	
	return (void*)newDensityField;
}

void* _DensityField_DefaultNew( Name name ) {
	return (void*) _DensityField_New(
		sizeof(DensityField),
		DensityField_Type,
		_DensityField_Delete,
		_DensityField_Print,
		_DensityField_Copy,
		_DensityField_DefaultNew,
		_DensityField_AssignFromXML,
		_DensityField_Build, 
		_DensityField_Initialise,
		_DensityField_Execute,
		_DensityField_Destroy,
		_FeVariable_InterpolateValueAt,
		_FeVariable_GetMinGlobalFieldMagnitude,
		_FeVariable_GetMaxGlobalFieldMagnitude,
		_FeVariable_GetMinAndMaxLocalCoords,
		_FeVariable_GetMinAndMaxGlobalCoords,
		_FeVariable_InterpolateNodeValuesToElLocalCoord,
		_FeVariable_GetValueAtNode,
		_DensityField_ValueAtParticle,
		name );
}

void _DensityField_AssignFromXML( void* densityField, Stg_ComponentFactory* cf, void* data ){
	DensityField*			self = (DensityField*) densityField;
	BuoyancyForceTerm*	buoyancyForceTerm;
	Variable_Register*	variable_Register;

	_ParticleFeVariable_AssignFromXML( self, cf, data );

	buoyancyForceTerm = Stg_ComponentFactory_ConstructByKey( cf, self->name, "BuoyancyForceTerm", BuoyancyForceTerm, True, data );
	variable_Register = self->context->variable_Register; 
	assert( variable_Register );

	_DensityField_Init( self, buoyancyForceTerm, variable_Register );
}

void _DensityField_Build( void* densityField, void* data ) {
	DensityField* self = (DensityField*) densityField;
	Variable_Register* variable_Register = (Variable_Register*) self->variable_Register;
	Name              tmpName;

	Stg_Component_Build( self->buoyancyForceTerm, data, False );

  	/* Create Dof Layout */
	assert( Class_IsSuper( self->feMesh->topo, IGraph ) );
	tmpName = Stg_Object_AppendSuffix( self, "densityVariable" );
	self->dataVariable = Variable_NewScalar( 	
			tmpName,
			Variable_DataType_Double, 
			&((IGraph*)self->feMesh->topo)->remotes[MT_VERTEX]->nDomains, 
			NULL,
			(void**)&self->data, 
			variable_Register );
	Memory_Free( tmpName );
	self->fieldComponentCount = 1;
	
	tmpName = Stg_Object_AppendSuffix( self, "densityDOF" );
	self->dofLayout = DofLayout_New( tmpName, variable_Register, 0, self->feMesh );
	self->dofLayout->_numItemsInLayout = FeMesh_GetNodeDomainSize( self->feMesh );
	DofLayout_AddAllFromVariableArray( self->dofLayout, 1, &self->dataVariable );
	Memory_Free( tmpName );
	self->eqNum->dofLayout = self->dofLayout;

	_ParticleFeVariable_Build( self, data );
}
void _DensityField_Initialise( void* densityField, void* data ) {
	DensityField* self = (DensityField*) densityField;

	Stg_Component_Initialise( self->buoyancyForceTerm, data, False );

	_ParticleFeVariable_Initialise( self, data );
}
void _DensityField_Execute( void* densityField, void* data ) {
	DensityField* self = (DensityField*) densityField;

	_ParticleFeVariable_Execute( self, data );
}
void _DensityField_Destroy( void* densityField, void* data ) {
	DensityField* self = (DensityField*) densityField;

   Stg_Component_Destroy( self->buoyancyForceTerm, data, False );

	_ParticleFeVariable_Destroy( self, data );
}

void _DensityField_ValueAtParticle( void* densityField, IntegrationPointsSwarm* swarm, Element_LocalIndex lElement_I, void* _particle, double* density ) {
	DensityField*                    self         = (DensityField*) densityField;
	IntegrationPoint*                particle     = (IntegrationPoint*) _particle;
	Material*                        material = NULL;
	BuoyancyForceTerm_MaterialExt*   materialExt;
	
	/* Calculate density from particle material */
	material = IntegrationPointsSwarm_GetMaterialOn( (IntegrationPointsSwarm*)swarm, particle );
	materialExt = ExtensionManager_Get( material->extensionMgr, material, self->buoyancyForceTerm->materialExtHandle );
	*density = materialExt->density;
}

