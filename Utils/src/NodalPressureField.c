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
#include "NodalPressureField.h"
#include <assert.h>
#include <string.h>

const Type NodalPressureField_Type = "NodalPressureField";

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
		FieldVariable_GetValueFunction*	                  _getMinGlobalFeMagnitude,
		FieldVariable_GetValueFunction*                   _getMaxGlobalFeMagnitude,
		FieldVariable_GetCoordFunction*                   _getMinAndMaxLocalCoords,
		FieldVariable_GetCoordFunction*                   _getMinAndMaxGlobalCoords,		
		FeVariable_InterpolateWithinElementFunction*      _interpolateWithinElement,	
		FeVariable_GetValueAtNodeFunction*                _getValueAtNode,
		ParticleFeVariable_ValueAtParticleFunction*       _valueAtParticle,
		Name                                              name )
{
	NodalPressureField*		self;
	
	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(NodalPressureField) );
	self = (NodalPressureField*)
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

void _NodalPressureField_Init( 
		NodalPressureField*                                      self,
		Variable_Register*                                variable_Register )
{
	/* Assign Pointers */
	self->variable_Register = variable_Register;

	self->fieldComponentCount = 1;
}

/* --- Virtual Function Implementations --- */
void _NodalPressureField_Delete( void* nodalPressureField ) {
	NodalPressureField* self = (NodalPressureField*) nodalPressureField;

	_ParticleFeVariable_Delete( self );
}

void _NodalPressureField_Print( void* nodalPressureField, Stream* stream ) {
	NodalPressureField* self = (NodalPressureField*) nodalPressureField;
	
	/* General info */
	Journal_Printf( stream, "NodalPressureField (ptr): %p\n", self );
	
	/* Print parent */
	_FeVariable_Print( self, stream );
}


void* _NodalPressureField_Copy( void* feVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	NodalPressureField*	self = (NodalPressureField*)feVariable;
	NodalPressureField*	newNodalPressureField;
	
	newNodalPressureField = (NodalPressureField*) _FeVariable_Copy( feVariable, dest, deep, nameExt, ptrMap );
	
	return (void*)newNodalPressureField;
}

void* _NodalPressureField_DefaultNew( Name name ) {
	return (void*) _NodalPressureField_New(
		sizeof(NodalPressureField),
		NodalPressureField_Type,
		_NodalPressureField_Delete,
		_NodalPressureField_Print,
		_NodalPressureField_Copy,
		_NodalPressureField_DefaultNew,
		_NodalPressureField_Construct,
		_NodalPressureField_Build, 
		_NodalPressureField_Initialise,
		_NodalPressureField_Execute,
		_NodalPressureField_Destroy,
		_FeVariable_InterpolateValueAt,
		_FeVariable_GetMinGlobalFieldMagnitude,
		_FeVariable_GetMaxGlobalFieldMagnitude,
		_FeVariable_GetMinAndMaxLocalCoords,
		_FeVariable_GetMinAndMaxGlobalCoords,
		_FeVariable_InterpolateNodeValuesToElLocalCoord,
		_FeVariable_GetValueAtNode,
		_NodalPressureField_ValueAtParticle,
		name );
}

void NodalPressureField_NonLinearUpdate( void* _sle, void* _ctx ) {
   SystemLinearEquations* sle = (SystemLinearEquations*)_sle;
   AbstractContext* ctx = (AbstractContext*)_ctx;
   FieldVariable_Register* fieldVar_Register;
   NodalPressureField* preVar;

   fieldVar_Register = Stg_ObjectList_Get( ctx->register_Register, "FieldVariable_Register" );
   preVar = (NodalPressureField*)FieldVariable_Register_GetByName( fieldVar_Register, "NodalPressureField" );
   ParticleFeVariable_Update( preVar );
}

void _NodalPressureField_Construct( void* nodalPressureField, Stg_ComponentFactory* cf, void* data ){
	NodalPressureField*          self              = (NodalPressureField*) nodalPressureField;
	Variable_Register*    variable_Register;
        SystemLinearEquations* sle;

	/* Construct Parent */
	_ParticleFeVariable_Construct( self, cf, data );

	/* _FieldVariable_Construct( self, cf, data ); */

	variable_Register      = (Variable_Register*) Stg_ObjectList_Get( cf->registerRegister, "Variable_Register" );
	assert( variable_Register );

        self->pressureField = Stg_ComponentFactory_ConstructByKey( cf, self->name, "PressureField", FeVariable, True, data );

	_NodalPressureField_Init( self, variable_Register );

        /*
        ** If we're using this field for non-linear feedback, we'll need to update it in between
        ** non-linear iterations. */
        sle = Stg_ComponentFactory_ConstructByKey( cf, self->name, "SLE", SystemLinearEquations, False, data );
        if( sle )
           SystemLinearEquations_AddPostNonLinearEP( sle, NodalPressureField_Type, NodalPressureField_NonLinearUpdate );
}

void _NodalPressureField_Build( void* nodalPressureField, void* data ) {
	NodalPressureField* self = (NodalPressureField*) nodalPressureField;
	Name              tmpName, tmpName2;
	Variable_Index variable_I;
	Node_DomainIndex  node_I;

	Stg_Component_Build( self->feMesh, data, False );

	/* Create Variable to store data */
	assert( Class_IsSuper( self->feMesh->topo, IGraph ) );
	tmpName = Stg_Object_AppendSuffix( self, "DataVariable" );
	self->dataVariable = Variable_NewScalar(
			tmpName,
			Variable_DataType_Double, 
			&((IGraph*)self->feMesh->topo)->remotes[MT_VERTEX]->nDomains, 
			NULL,
			(void**)&self->data, 
			self->variable_Register );
	
	/* Create Dof Layout */
	tmpName2 = Stg_Object_AppendSuffix( self, "DofLayout" );
	self->dofLayout = DofLayout_New( tmpName, self->variable_Register, 0, self->feMesh );
	self->dofLayout->_numItemsInLayout = FeMesh_GetNodeDomainSize( self->feMesh );
        for( node_I = 0; node_I < FeMesh_GetNodeDomainSize( self->feMesh ); node_I++ ) {
           DofLayout_AddDof_ByVarName( self->dofLayout, tmpName, node_I );
        }
	Memory_Free( tmpName );
	Memory_Free( tmpName2 );
	self->eqNum->dofLayout = self->dofLayout;

	/* Build and Update all Variables that this component has created */
	Stg_Component_Build( self->dataVariable, data, False); Variable_Update( self->dataVariable );

	_ParticleFeVariable_Build( self, data );
	/* Update again, just in case things were changed/reallocated when ICs loaded */
	Variable_Update( self->dataVariable );

}

void _NodalPressureField_Initialise( void* nodalPressureField, void* data ) {
	NodalPressureField* self = (NodalPressureField*) nodalPressureField;
	Variable_Index variable_I;
	
	/* Initialise and Update all Variables that this component has created */
	Stg_Component_Initialise( self->dataVariable, data, False); Variable_Update( self->dataVariable );
	
	_ParticleFeVariable_Initialise( self, data );

	/* Do a post-update just in case */
	Variable_Update( self->dataVariable );

}
void _NodalPressureField_Execute( void* nodalPressureField, void* data ) {
	NodalPressureField* self = (NodalPressureField*) nodalPressureField;

	_ParticleFeVariable_Execute( self, data );
}
void _NodalPressureField_Destroy( void* nodalPressureField, void* data ) {
	NodalPressureField* self = (NodalPressureField*) nodalPressureField;

	_ParticleFeVariable_Destroy( self, data );
}

void _NodalPressureField_ValueAtParticle( void* nodalPressureField, IntegrationPointsSwarm* swarm, Element_LocalIndex lElement_I, void* _particle, double* pressure ) {
	NodalPressureField*      self         = (NodalPressureField*) nodalPressureField;
	SymmetricTensor   strainRate;
	IntegrationPoint* particle     = (IntegrationPoint*) _particle;

        FeVariable_InterpolateWithinElement( self->pressureField, lElement_I, particle->xi, pressure );
	
}

