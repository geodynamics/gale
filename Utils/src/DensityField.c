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

DensityField* _DensityField_New(  DENSITYFIELD_DEFARGS  )
{
	DensityField*		self;
	
	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(DensityField) );
	self = (DensityField*)
		_ParticleFeVariable_New(  PARTICLEFEVARIABLE_PASSARGS  );
	
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
	/* Variables set in this function */
	SizeT                                                       _sizeOfSelf = sizeof(DensityField);
	Type                                                               type = DensityField_Type;
	Stg_Class_DeleteFunction*                                       _delete = _DensityField_Delete;
	Stg_Class_PrintFunction*                                         _print = _DensityField_Print;
	Stg_Class_CopyFunction*                                           _copy = _DensityField_Copy;
	Stg_Component_DefaultConstructorFunction*           _defaultConstructor = _DensityField_DefaultNew;
	Stg_Component_ConstructFunction*                             _construct = _DensityField_AssignFromXML;
	Stg_Component_BuildFunction*                                     _build = _DensityField_Build;
	Stg_Component_InitialiseFunction*                           _initialise = _DensityField_Initialise;
	Stg_Component_ExecuteFunction*                                 _execute = _DensityField_Execute;
	Stg_Component_DestroyFunction*                                 _destroy = _DensityField_Destroy;
	FieldVariable_InterpolateValueAtFunction*           _interpolateValueAt = _FeVariable_InterpolateValueAt;
	FieldVariable_GetCoordFunction*                _getMinAndMaxLocalCoords = _FeVariable_GetMinAndMaxLocalCoords;
	FieldVariable_GetCoordFunction*               _getMinAndMaxGlobalCoords = _FeVariable_GetMinAndMaxGlobalCoords;
	FeVariable_InterpolateWithinElementFunction*  _interpolateWithinElement = _FeVariable_InterpolateNodeValuesToElLocalCoord;
	FeVariable_GetValueAtNodeFunction*                      _getValueAtNode = _FeVariable_GetValueAtNode;
	ParticleFeVariable_ValueAtParticleFunction*            _valueAtParticle = _DensityField_ValueAtParticle;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType                             nameAllocationType = ZERO;
	FieldVariable_GetValueFunction*   _getMinGlobalFieldMagnitude = ZERO;
	FieldVariable_GetValueFunction*   _getMaxGlobalFieldMagnitude = ZERO;
	FeVariable_SyncShadowValuesFunc*            _syncShadowValues = ZERO;

	return (void*) _DensityField_New(  DENSITYFIELD_PASSARGS  );
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
			(AbstractContext*)self->context,
			Variable_DataType_Double, 
			&((IGraph*)self->feMesh->topo)->remotes[MT_VERTEX]->nDomains, 
			NULL,
			(void**)&self->data, 
			variable_Register );
	Memory_Free( tmpName );
	self->fieldComponentCount = 1;
	
	tmpName = Stg_Object_AppendSuffix( self, "densityDOF" );
	self->dofLayout = DofLayout_New( tmpName, self->context, variable_Register, 0, self->feMesh );
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



