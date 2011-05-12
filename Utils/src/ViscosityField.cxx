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

ViscosityField* _ViscosityField_New(  VISCOSITYFIELD_DEFARGS  ) {
	ViscosityField* self;
	
	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree.
		At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(ViscosityField) );
	self = (ViscosityField*) _ParticleFeVariable_New(  PARTICLEFEVARIABLE_PASSARGS  );
	
	return self;
}

void _ViscosityField_Init( 
	ViscosityField*		self,
	ConstitutiveMatrix*	constitutiveMatrix,
	Variable_Register*	variable_Register )
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


void* _ViscosityField_Copy( const void* feVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	ViscosityField*	self = (ViscosityField*)feVariable;
	ViscosityField*	newViscosityField;
	
	newViscosityField = (ViscosityField*) _FeVariable_Copy( feVariable, dest, deep, nameExt, ptrMap );

	newViscosityField->constitutiveMatrix = self->constitutiveMatrix;
	
	return (void*)newViscosityField;
}

void* _ViscosityField_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                       _sizeOfSelf = sizeof(ViscosityField);
	Type                                                               type = ViscosityField_Type;
	Stg_Class_DeleteFunction*                                       _delete = _ViscosityField_Delete;
	Stg_Class_PrintFunction*                                         _print = _ViscosityField_Print;
	Stg_Class_CopyFunction*                                           _copy = _ViscosityField_Copy;
	Stg_Component_DefaultConstructorFunction*           _defaultConstructor = _ViscosityField_DefaultNew;
	Stg_Component_ConstructFunction*                             _construct = _ViscosityField_AssignFromXML;
	Stg_Component_BuildFunction*                                     _build = _ViscosityField_Build;
	Stg_Component_InitialiseFunction*                           _initialise = _ViscosityField_Initialise;
	Stg_Component_ExecuteFunction*                                 _execute = _ViscosityField_Execute;
	Stg_Component_DestroyFunction*                                 _destroy = _ViscosityField_Destroy;
	FieldVariable_InterpolateValueAtFunction*           _interpolateValueAt = _FeVariable_InterpolateValueAt;
	FieldVariable_GetCoordFunction*                _getMinAndMaxLocalCoords = _FeVariable_GetMinAndMaxLocalCoords;
	FieldVariable_GetCoordFunction*               _getMinAndMaxGlobalCoords = _FeVariable_GetMinAndMaxGlobalCoords;
	FeVariable_InterpolateWithinElementFunction*  _interpolateWithinElement = _FeVariable_InterpolateNodeValuesToElLocalCoord;
	FeVariable_GetValueAtNodeFunction*                      _getValueAtNode = _FeVariable_GetValueAtNode;
	ParticleFeVariable_ValueAtParticleFunction*            _valueAtParticle = _ViscosityField_ValueAtParticle;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType                             nameAllocationType = (AllocationType)ZERO;
	FieldVariable_GetValueFunction*   _getMinGlobalFieldMagnitude = ZERO;
	FieldVariable_GetValueFunction*   _getMaxGlobalFieldMagnitude = ZERO;
	FeVariable_SyncShadowValuesFunc*            _syncShadowValues = ZERO;

	return (void*) _ViscosityField_New(  VISCOSITYFIELD_PASSARGS  );
}

void _ViscosityField_AssignFromXML( void* viscosityField, Stg_ComponentFactory* cf, void* data ){
	ViscosityField*		self = (ViscosityField*) viscosityField;
	ConstitutiveMatrix*	constitutiveMatrix;
	Variable_Register*	variable_Register;

	/* Construct Parent */
	_ParticleFeVariable_AssignFromXML( self, cf, data );

	constitutiveMatrix = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"ConstitutiveMatrix", ConstitutiveMatrix, True, data );
	variable_Register = self->context->variable_Register; 
	assert( variable_Register  );

	_ViscosityField_Init( self, constitutiveMatrix, variable_Register );
}

void _ViscosityField_Build( void* viscosityField, void* data ) {
	ViscosityField* self = (ViscosityField*) viscosityField;
	char*              tmpName;

	Stg_Component_Build( self->feMesh, data, False );
	Stg_Component_Build( self->constitutiveMatrix, data, False );

	/* Create Dof Layout */
	assert( Class_IsSuper( self->feMesh->topo, IGraph ) );
	tmpName = Stg_Object_AppendSuffix( self, (Name)"viscosityVariable"  );
	self->dataVariable = Variable_NewScalar( tmpName, (AbstractContext*)self->context, Variable_DataType_Double, (Index*)&((IGraph*)self->feMesh->topo)->remotes[MT_VERTEX]->nDomains, NULL, (void**)&self->data, self->variable_Register );
	Memory_Free( tmpName  );
	self->fieldComponentCount = 1;
	
	tmpName = Stg_Object_AppendSuffix( self, (Name)"viscosityDOF"  );
	self->dofLayout = DofLayout_New( tmpName, self->context, self->variable_Register, 0, self->feMesh );
	self->dofLayout->_numItemsInLayout = FeMesh_GetNodeDomainSize( self->feMesh );
	DofLayout_AddAllFromVariableArray( self->dofLayout, 1, &self->dataVariable );
	Memory_Free( tmpName );
	self->eqNum->dofLayout = self->dofLayout;

	_ParticleFeVariable_Build( self, data );
}
void _ViscosityField_Initialise( void* viscosityField, void* data ) {
	ViscosityField* self = (ViscosityField*) viscosityField;

   Stg_Component_Initialise( self->constitutiveMatrix, data, False );
	_ParticleFeVariable_Initialise( self, data );
}
void _ViscosityField_Execute( void* viscosityField, void* data ) {
	ViscosityField* self = (ViscosityField*) viscosityField;

	_ParticleFeVariable_Execute( self, data );
}
void _ViscosityField_Destroy( void* viscosityField, void* data ) {
	ViscosityField* self = (ViscosityField*) viscosityField;

   Stg_Component_Destroy( self->constitutiveMatrix, data, False );

	_ParticleFeVariable_Destroy( self, data );
}
void _ViscosityField_ValueAtParticle( void* viscosityField, IntegrationPointsSwarm* swarm, Element_LocalIndex lElement_I, void* _particle, double* viscosity ) {
	ViscosityField*   self         = (ViscosityField*) viscosityField;
	IntegrationPoint* particle     = (IntegrationPoint*) _particle;
	
	/* Calculate viscosity from constitutive matrix */
	ConstitutiveMatrix_Assemble( self->constitutiveMatrix, lElement_I,
                                     self->currentParticleIndex, particle );
	*viscosity = ConstitutiveMatrix_GetIsotropicViscosity( self->constitutiveMatrix );
}



