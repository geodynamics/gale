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
#include "ViscosityField.h"
#include <assert.h>

const Type ViscosityField_Type = "ViscosityField";

ViscosityField* _ViscosityField_New(
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
	ViscosityField*		self;
	
	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(ViscosityField) );
	self = (ViscosityField*)
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

void _ViscosityField_Init( 
		ViscosityField*                                   self,
		ConstitutiveMatrix*                               constitutiveMatrix,
		Variable_Register*                                variable_Register )
{
	/* Assign Pointers */
	self->variable_Register = variable_Register;
	self->constitutiveMatrix = constitutiveMatrix;
	
	/* Set pointers to swarm to be the same as the one on the constitutive matrix */
	self->assemblyTerm->integrationSwarm = self->constitutiveMatrix->integrationSwarm;
	self->massMatrixForceTerm->integrationSwarm = self->constitutiveMatrix->integrationSwarm;	
}

/* --- Virtual Function Implementations --- */
void _ViscosityField_Delete( void* viscosityField ) {
	ViscosityField* self = (ViscosityField*) viscosityField;

	Stg_Class_Delete( self->assemblyVector );
	Memory_Free( self->assemblyVectorName );

	_FeVariable_Delete( self );
}

void _ViscosityField_Print( void* viscosityField, Stream* stream ) {
	ViscosityField* self = (ViscosityField*) viscosityField;
	
	/* General info */
	Journal_Printf( stream, "ViscosityField (ptr): %p\n", self );
	
	/* Print parent */
	_FeVariable_Print( self, stream );
	
	/* ViscosityField info */
	Journal_PrintPointer( stream, self->constitutiveMatrix );
}


void* _ViscosityField_Copy( void* feVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	ViscosityField*	self = (ViscosityField*)feVariable;
	ViscosityField*	newViscosityField;
	
	newViscosityField = (ViscosityField*) _FeVariable_Copy( feVariable, dest, deep, nameExt, ptrMap );

	newViscosityField->constitutiveMatrix = self->constitutiveMatrix;
	
	return (void*)newViscosityField;
}

void* _ViscosityField_DefaultNew( Name name ) {
	return (void*) _ViscosityField_New(
		sizeof(ViscosityField),
		ViscosityField_Type,
		_ViscosityField_Delete,
		_ViscosityField_Print,
		_ViscosityField_Copy,
		_ViscosityField_DefaultNew,
		_ViscosityField_Construct,
		_ViscosityField_Build, 
		_ViscosityField_Initialise,
		_ViscosityField_Execute,
		_ViscosityField_Destroy,
		_FeVariable_InterpolateValueAt,
		_FeVariable_GetMinGlobalFieldMagnitude,
		_FeVariable_GetMaxGlobalFieldMagnitude,
		_FeVariable_GetMinAndMaxLocalCoords,
		_FeVariable_GetMinAndMaxGlobalCoords,
		_FeVariable_InterpolateNodeValuesToElLocalCoord,
		_FeVariable_GetValueAtNode,
		_ViscosityField_ValueAtParticle,
		name );
}

void _ViscosityField_Construct( void* viscosityField, Stg_ComponentFactory* cf, void* data ){
	ViscosityField*          self              = (ViscosityField*) viscosityField;
	ConstitutiveMatrix*   constitutiveMatrix;
	Variable_Register*    variable_Register;

	/* Construct Parent */
	_ParticleFeVariable_Construct( self, cf, data );

	/*_FieldVariable_Construct( self, cf, data );*/

	constitutiveMatrix = Stg_ComponentFactory_ConstructByKey( cf, self->name, "ConstitutiveMatrix", ConstitutiveMatrix, True, data );
	variable_Register      = (Variable_Register*) Stg_ObjectList_Get( cf->registerRegister, "Variable_Register" );
	assert( variable_Register );

	_ViscosityField_Init( self, constitutiveMatrix, variable_Register );
}

void _ViscosityField_Build( void* viscosityField, void* data ) {
	ViscosityField* self = (ViscosityField*) viscosityField;
	Name              tmpName;

	Stg_Component_Build( self->feMesh, data, False );

	/* Create Dof Layout */
	assert( Class_IsSuper( self->feMesh->topo, IGraph ) );
	tmpName = Stg_Object_AppendSuffix( self, "DataVariable" );
	self->dataVariable = Variable_NewScalar( 	
			tmpName,
			Variable_DataType_Double, 
			&((IGraph*)self->feMesh->topo)->remotes[MT_VERTEX]->nDomains, 
			NULL,
			(void**)&self->data, 
			self->variable_Register );
	Memory_Free( tmpName );
	self->fieldComponentCount = 1;
	
	tmpName = Stg_Object_AppendSuffix( self, "DofLayout" );
	self->dofLayout = DofLayout_New( tmpName, self->variable_Register, 0, self->feMesh );
	self->dofLayout->_numItemsInLayout = FeMesh_GetNodeDomainSize( self->feMesh );
	DofLayout_AddAllFromVariableArray( self->dofLayout, 1, &self->dataVariable );
	Memory_Free( tmpName );
	self->eqNum->dofLayout = self->dofLayout;

	_ParticleFeVariable_Build( self, data );
}
void _ViscosityField_Initialise( void* viscosityField, void* data ) {
	ViscosityField* self = (ViscosityField*) viscosityField;

	_ParticleFeVariable_Initialise( self, data );
}
void _ViscosityField_Execute( void* viscosityField, void* data ) {
	ViscosityField* self = (ViscosityField*) viscosityField;

	_ParticleFeVariable_Execute( self, data );
}
void _ViscosityField_Destroy( void* viscosityField, void* data ) {
	ViscosityField* self = (ViscosityField*) viscosityField;

	_ParticleFeVariable_Destroy( self, data );
}
void _ViscosityField_ValueAtParticle( void* viscosityField, IntegrationPointsSwarm* swarm, Element_LocalIndex lElement_I, void* _particle, double* viscosity ) {
	ViscosityField*   self         = (ViscosityField*) viscosityField;
	IntegrationPoint* particle     = (IntegrationPoint*) _particle;
	
	/* Calculate viscosity from constitutive matrix */
	ConstitutiveMatrix_Assemble( self->constitutiveMatrix, lElement_I, particle );
	*viscosity = ConstitutiveMatrix_GetIsotropicViscosity( self->constitutiveMatrix );
}

