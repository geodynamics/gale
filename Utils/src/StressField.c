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
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include <Underworld/Rheology/Rheology.h>

#include "types.h"
#include "StressField.h"
#include <assert.h>
#include <string.h>

const Type StressField_Type = "StressField";

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
		Name                                              name )
{
	StressField*		self;
	
	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(StressField) );
	self = (StressField*)
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

void _StressField_Init( 
		StressField*                                      self,
		FeVariable*                                       strainRateField,
		ConstitutiveMatrix*                               constitutiveMatrix,
		Variable*                                         stressVariable,
		Variable_Register*                                variable_Register )
{
	Dimension_Index   dim = constitutiveMatrix->dim;

	/* Assign Pointers */
	self->strainRateField    = strainRateField;
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
}

/* --- Virtual Function Implementations --- */
void _StressField_Delete( void* stressField ) {
	StressField* self = (StressField*) stressField;

	Stg_Class_Delete( self->assemblyVector );
	Memory_Free( self->assemblyVectorName );

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
	return (void*) _StressField_New(
		sizeof(StressField),
		StressField_Type,
		_StressField_Delete,
		_StressField_Print,
		_StressField_Copy,
		_StressField_DefaultNew,
		_StressField_Construct,
		_StressField_Build, 
		_StressField_Initialise,
		_StressField_Execute,
		_StressField_Destroy,
		_FeVariable_InterpolateValueAt,
		_FeVariable_GetMinGlobalFieldMagnitude,
		_FeVariable_GetMaxGlobalFieldMagnitude,
		_FeVariable_GetMinAndMaxLocalCoords,
		_FeVariable_GetMinAndMaxGlobalCoords,
		_FeVariable_InterpolateNodeValuesToElLocalCoord,
		_FeVariable_GetValueAtNode,
		_StressField_ValueAtParticle_Recalculate,
		name );
}

void _StressField_Construct( void* stressField, Stg_ComponentFactory* cf, void* data ){
	StressField*          self              = (StressField*) stressField;
	FeVariable*           strainRateField;
	ConstitutiveMatrix*   constitutiveMatrix;
	Variable_Register*    variable_Register;
	Name                  stressVariableName;
	Variable*             stressVariable;

	/* Construct Parent */
	_ParticleFeVariable_Construct( self, cf, data );

	/* _FieldVariable_Construct( self, cf, data ); */

	strainRateField =  Stg_ComponentFactory_ConstructByKey( cf,  self->name, "StrainRateField", FeVariable, True, data );
	constitutiveMatrix = Stg_ComponentFactory_ConstructByKey( cf, self->name, "ConstitutiveMatrix", ConstitutiveMatrix, True, data );
	variable_Register      = (Variable_Register*) Stg_ObjectList_Get( cf->registerRegister, "Variable_Register" );
	assert( variable_Register );

	stressVariableName = Stg_ComponentFactory_GetString( cf, self->name, "StressVariable", "Stress" );
	stressVariable = Variable_Register_GetByName( variable_Register, stressVariableName );

	_StressField_Init( self, strainRateField, constitutiveMatrix, stressVariable, variable_Register );
}

void _StressField_Build( void* stressField, void* data ) {
	StressField* self = (StressField*) stressField;
	Name              tmpName;
	Name              variableName[6];
	Dimension_Index   dim = self->constitutiveMatrix->dim;
	Variable_Index variable_I;
	Node_DomainIndex  node_I;

	Stg_Component_Build( self->feMesh, data, False );

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
	tmpName = Stg_Object_AppendSuffix( self, "DataVariable" );
	self->dataVariable = Variable_NewVector( 	
			tmpName,
			Variable_DataType_Double, 
			self->fieldComponentCount,
			&self->feMesh->topo->remotes[MT_VERTEX]->nDomains, 
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
	self->dofLayout = DofLayout_New( tmpName, self->variable_Register, 0, self->feMesh );
	self->dofLayout->_numItemsInLayout = FeMesh_GetNodeDomainSize( self->feMesh );
	for( variable_I = 0; variable_I < self->fieldComponentCount ; variable_I++ ) {
		self->dataVariableList[ variable_I ] = Variable_Register_GetByName( self->variable_Register, 
										    variableName[ variable_I ] );
		for( node_I = 0; node_I < FeMesh_GetNodeDomainSize( self->feMesh ); node_I++ ) {
			DofLayout_AddDof_ByVarName( self->dofLayout, variableName[variable_I], node_I );
		}
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

	_ParticleFeVariable_Destroy( self, data );
}
void _StressField_ValueAtParticle_Recalculate( void* stressField, IntegrationPointsSwarm* swarm, Element_LocalIndex lElement_I, void* _particle, double* stress ) {
	StressField*      self         = (StressField*) stressField;
	SymmetricTensor   strainRate;
	IntegrationPoint* particle     = (IntegrationPoint*) _particle;
	
	/* Calculate stress from strain rate and constitutive matrix */
	ConstitutiveMatrix_Assemble( self->constitutiveMatrix, lElement_I, particle );
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
