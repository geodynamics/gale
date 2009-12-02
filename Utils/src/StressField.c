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
#include "StressField.h"
#include <assert.h>
#include <string.h>

const Type StressField_Type = "StressField";

StressField* _StressField_New(  STRESSFIELD_DEFARGS  ) {
	StressField* self;
	
	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree.
		At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(StressField) );
	self = (StressField*) _ParticleFeVariable_New(  PARTICLEFEVARIABLE_PASSARGS  );
	
	return self;
}

void _StressField_Init( 
	StressField*				self,
	FeVariable*					strainRateField,
	ConstitutiveMatrix*		constitutiveMatrix,
	Variable*					stressVariable,
	Variable_Register*		variable_Register,
	SystemLinearEquations*	sle)
{
	Dimension_Index dim = constitutiveMatrix->dim;

	/* Assign Pointers */
	self->strainRateField = strainRateField;
	self->constitutiveMatrix = constitutiveMatrix;
	self->variable_Register = variable_Register;

	self->fieldComponentCount = StGermain_nSymmetricTensorVectorComponents( dim );

	if ( stressVariable ) {
		self->stressVariable = stressVariable;
		self->_valueAtParticle = _StressField_ValueAtParticle_FromVariable;
	}
		
	/* Set pointers to swarm to be the same as the one on the constitutive matrix */
	self->assemblyTerm->integrationSwarm = self->constitutiveMatrix->integrationSwarm;
	self->massMatrixForceTerm->integrationSwarm = self->constitutiveMatrix->integrationSwarm;

	/*
	** If we're using this field for non-linear feedback, we'll need to update it in between
	** non-linear iterations. */
	if( sle )
		SystemLinearEquations_AddPostNonLinearEP( sle, StressField_Type, StressField_NonLinearUpdate );
}

/* --- Virtual Function Implementations --- */
void _StressField_Delete( void* stressField ) {
	StressField* self = (StressField*) stressField;

	_FeVariable_Delete( self );
}

void _StressField_Print( void* stressField, Stream* stream ) {
	StressField* self = (StressField*) stressField;
	
	/* General info */
	Journal_Printf( stream, "StressField (ptr): %p\n", self );
	
	/* Print parent */
	_FeVariable_Print( self, stream );
	
	/* StressField info */
	Journal_PrintPointer( stream, self->strainRateField );
	Journal_PrintPointer( stream, self->constitutiveMatrix );
}

void* _StressField_Copy( void* feVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	StressField*	self = (StressField*)feVariable;
	StressField*	newStressField;
	
	newStressField = (StressField*) _FeVariable_Copy( feVariable, dest, deep, nameExt, ptrMap );

	newStressField->strainRateField = self->strainRateField;
	newStressField->constitutiveMatrix = self->constitutiveMatrix;
	
	return (void*)newStressField;
}

void* _StressField_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                       _sizeOfSelf = sizeof(StressField);
	Type                                                               type = StressField_Type;
	Stg_Class_DeleteFunction*                                       _delete = _StressField_Delete;
	Stg_Class_PrintFunction*                                         _print = _StressField_Print;
	Stg_Class_CopyFunction*                                           _copy = _StressField_Copy;
	Stg_Component_DefaultConstructorFunction*           _defaultConstructor = _StressField_DefaultNew;
	Stg_Component_ConstructFunction*                             _construct = _StressField_AssignFromXML;
	Stg_Component_BuildFunction*                                     _build = _StressField_Build;
	Stg_Component_InitialiseFunction*                           _initialise = _StressField_Initialise;
	Stg_Component_ExecuteFunction*                                 _execute = _StressField_Execute;
	Stg_Component_DestroyFunction*                                 _destroy = _StressField_Destroy;
	FieldVariable_InterpolateValueAtFunction*           _interpolateValueAt = _FeVariable_InterpolateValueAt;
	FieldVariable_GetCoordFunction*                _getMinAndMaxLocalCoords = _FeVariable_GetMinAndMaxLocalCoords;
	FieldVariable_GetCoordFunction*               _getMinAndMaxGlobalCoords = _FeVariable_GetMinAndMaxGlobalCoords;
	FeVariable_InterpolateWithinElementFunction*  _interpolateWithinElement = _FeVariable_InterpolateNodeValuesToElLocalCoord;
	FeVariable_GetValueAtNodeFunction*                      _getValueAtNode = _FeVariable_GetValueAtNode;
	ParticleFeVariable_ValueAtParticleFunction*            _valueAtParticle = _StressField_ValueAtParticle_Recalculate;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType                             nameAllocationType = ZERO;
	FieldVariable_GetValueFunction*   _getMinGlobalFieldMagnitude = ZERO;
	FieldVariable_GetValueFunction*   _getMaxGlobalFieldMagnitude = ZERO;
	FeVariable_SyncShadowValuesFunc*            _syncShadowValues = ZERO;

	return (void*) _StressField_New(  STRESSFIELD_PASSARGS  );
}

void StressField_NonLinearUpdate( void* _sle, void* _ctx ) {
   SystemLinearEquations* sle = (SystemLinearEquations*)_sle;
   DomainContext* ctx = (DomainContext*)_ctx;
   FieldVariable_Register* fieldVar_Register;
   StressField* stressVar;

   fieldVar_Register = ctx->fieldVariable_Register;
   stressVar = (StressField*)FieldVariable_Register_GetByName( fieldVar_Register, "StressField" );
   ParticleFeVariable_Update( stressVar );
}

void _StressField_AssignFromXML( void* stressField, Stg_ComponentFactory* cf, void* data ){
	StressField*          self = (StressField*) stressField;
	FeVariable*           strainRateField;
	ConstitutiveMatrix*   constitutiveMatrix;
	Variable_Register*    variable_Register;
	Name                  stressVariableName;
	Variable*             stressVariable;
	SystemLinearEquations* sle;

	/* Construct Parent */
	_ParticleFeVariable_AssignFromXML( self, cf, data );

	strainRateField =  Stg_ComponentFactory_ConstructByKey( cf,  self->name, "StrainRateField", FeVariable, True, data );
	constitutiveMatrix = Stg_ComponentFactory_ConstructByKey( cf, self->name, "ConstitutiveMatrix", ConstitutiveMatrix, True, data );
	variable_Register = self->context->variable_Register; 
	assert( variable_Register );

	stressVariableName = Stg_ComponentFactory_GetString( cf, self->name, "StressVariable", "Stress" );
	stressVariable = Variable_Register_GetByName( variable_Register, stressVariableName );

   sle = Stg_ComponentFactory_ConstructByKey( cf, self->name, "SLE", SystemLinearEquations, False, data );

	_StressField_Init( self, strainRateField, constitutiveMatrix, stressVariable, variable_Register, sle );
}

void _StressField_Build( void* stressField, void* data ) {
	StressField* self = (StressField*) stressField;
	Name              tmpName;
	Name              variableName[6];
	Dimension_Index   dim = self->constitutiveMatrix->dim;
	Variable_Index variable_I;
	Node_DomainIndex  node_I;

	Stg_Component_Build( self->feMesh, data, False );
	Stg_Component_Build( self->strainRateField, data, False );
	Stg_Component_Build( self->constitutiveMatrix, data, False );
	if ( self->stressVariable ) Stg_Component_Build( self->stressVariable, data, False );

	if ( dim == 2 ) {
		variableName[0] = StG_Strdup( "tau_xx" );
		variableName[1] = StG_Strdup( "tau_yy" );
		variableName[2] = StG_Strdup( "tau_xy" );
	}
	else {
		variableName[0] = StG_Strdup( "tau_xx" );
		variableName[1] = StG_Strdup( "tau_yy" );
		variableName[2] = StG_Strdup( "tau_zz" );
		variableName[3] = StG_Strdup( "tau_xy" );
		variableName[4] = StG_Strdup( "tau_xz" );
		variableName[5] = StG_Strdup( "tau_yz" );
	}

	/* Create Variable to store data */
	assert( Class_IsSuper( self->feMesh->topo, IGraph ) );
	tmpName = Stg_Object_AppendSuffix( self, "DataVariable" );
	self->dataVariable = Variable_NewVector(
		tmpName,
		(AbstractContext*)self->context,
		Variable_DataType_Double, 
		self->fieldComponentCount,
		&((IGraph*)self->feMesh->topo)->remotes[MT_VERTEX]->nDomains, 
		NULL,
		(void**)&self->data, 
		self->variable_Register,
		variableName[0],
		variableName[1],
		variableName[2],
		variableName[3],
		variableName[4],
		variableName[5] );
	Memory_Free( tmpName );
	
	/* Create Dof Layout */
	tmpName = Stg_Object_AppendSuffix( self, "DofLayout" );
	self->dofLayout = DofLayout_New( tmpName, self->context, self->variable_Register, 0, self->feMesh );
	self->dofLayout->_numItemsInLayout = FeMesh_GetNodeDomainSize( self->feMesh );
	for( variable_I = 0; variable_I < self->fieldComponentCount ; variable_I++ ) {
		self->dataVariableList[ variable_I ] = Variable_Register_GetByName( self->variable_Register, 
                             variableName[ variable_I ] );
		for( node_I = 0; node_I < FeMesh_GetNodeDomainSize( self->feMesh ); node_I++ )
			DofLayout_AddDof_ByVarName( self->dofLayout, variableName[variable_I], node_I );
		/* Free Name */
		Memory_Free( variableName[ variable_I ] );
	}
	Memory_Free( tmpName );
	self->eqNum->dofLayout = self->dofLayout;

	/* Build and Update all Variables that this component has created */
	Stg_Component_Build( self->dataVariable, data, False); Variable_Update( self->dataVariable );
	for( variable_I = 0; variable_I < self->fieldComponentCount ; variable_I++ ) {
		Stg_Component_Build( self->dataVariableList[ variable_I ], data, False); Variable_Update( self->dataVariableList[ variable_I ] );
	}

	_ParticleFeVariable_Build( self, data );
	/* Update again, just in case things were changed/reallocated when ICs loaded */
	Variable_Update( self->dataVariable );
	for( variable_I = 0; variable_I < self->fieldComponentCount ; variable_I++ ) {
		Variable_Update( self->dataVariableList[ variable_I ] );
	}
}

void _StressField_Initialise( void* stressField, void* data ) {
	StressField* self = (StressField*) stressField;
	Variable_Index variable_I;

	Stg_Component_Initialise( self->strainRateField, data, False );
	Stg_Component_Initialise( self->constitutiveMatrix, data, False );
	/* Initialise and Update all Variables that this component has created */
	Stg_Component_Initialise( self->dataVariable, data, False); Variable_Update( self->dataVariable );
	for( variable_I = 0; variable_I < self->fieldComponentCount ; variable_I++ ) {
		Stg_Component_Initialise( self->dataVariableList[ variable_I ], data, False); Variable_Update( self->dataVariableList[ variable_I ] );
	}
	
	_ParticleFeVariable_Initialise( self, data );

	/* Do a post-update just in case */
	Variable_Update( self->dataVariable );
	for( variable_I = 0; variable_I < self->fieldComponentCount ; variable_I++ ) {
		Variable_Update( self->dataVariableList[ variable_I ] );
	}

}

void _StressField_Execute( void* stressField, void* data ) {
	StressField* self = (StressField*) stressField;

	_ParticleFeVariable_Execute( self, data );
}

void _StressField_Destroy( void* stressField, void* data ) {
	StressField* self = (StressField*) stressField;

	Stg_Component_Destroy( self->strainRateField, data, False );
	Stg_Component_Destroy( self->constitutiveMatrix, data, False );
	Stg_Component_Destroy( self->dataVariable, data, False);

	_ParticleFeVariable_Destroy( self, data );
}

void _StressField_ValueAtParticle_Recalculate( void* stressField, IntegrationPointsSwarm* swarm, Element_LocalIndex lElement_I, void* _particle, double* stress ) {
	StressField*      self         = (StressField*) stressField;
	SymmetricTensor   strainRate;
	IntegrationPoint* particle     = (IntegrationPoint*) _particle;
	
	/* Calculate stress from strain rate and constitutive matrix */
	ConstitutiveMatrix_Assemble( self->constitutiveMatrix, lElement_I,
                                     self->currentParticleIndex, particle );
	FeVariable_InterpolateWithinElement( self->strainRateField, lElement_I, particle->xi, strainRate );
	ConstitutiveMatrix_CalculateStress( self->constitutiveMatrix, strainRate, stress );

	assert ( !isnan( stress[0] ) );
	assert ( !isnan( stress[1] ) );
	assert ( !isnan( stress[2] ) );
	
}

void _StressField_ValueAtParticle_FromVariable( void* stressField, IntegrationPointsSwarm* swarm, Element_LocalIndex lElement_I, void* particle, double* stress ) {
	StressField*      self         = (StressField*) stressField;
	double*           stressParticleExt;
	MaterialPointsSwarm* materialPointsSwarm;
	void* materialPoint;
	
	/* Ok, we break any fancy mapper stuff here and assume the materialParticles correspond exactly to the
	 * Integration ones: */
	materialPoint = OneToOneMapper_GetMaterialPoint( swarm->mapper, particle, &materialPointsSwarm );	

	/* Get pointer to stress using variable */
	stressParticleExt = (double*) ((ArithPointer) materialPoint + (ArithPointer)self->stressVariable->offsets[0]);
	memcpy( stress, stressParticleExt, sizeof(double)*self->fieldComponentCount );
}


