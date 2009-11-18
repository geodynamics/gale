/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Monash Cluster Computing, Australia
** (C) 2003-2004 All Rights Reserved
**
** Primary Authors:
** Luke Hodkinson
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Rheology/Rheology.h>

#include "types.h"
#include "NodalPressureField.h"

const Type NodalPressureField_Type = "NodalPressureField";

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
	Name                                              name )
{
   NodalPressureField* self;

   assert( _sizeOfSelf >= sizeof(NodalPressureField) );
   self = (NodalPressureField*)_ParticleFeVariable_New(
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
		name
  );

   return self;
}

void _NodalPressureField_Init( NodalPressureField* self, Variable_Register* variable_Register, FeVariable* pressureField) {
   self->variable_Register = variable_Register;
   self->fieldComponentCount = 1;
   self->pressureField = pressureField;
}

void* _NodalPressureField_DefaultNew( Name name ) {
   return (void*) _NodalPressureField_New(
      sizeof(NodalPressureField),
      NodalPressureField_Type,
      _NodalPressureField_Delete,
      _NodalPressureField_Print,
      _NodalPressureField_Copy,
      _NodalPressureField_DefaultNew,
      _NodalPressureField_AssignFromXML,
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

void _NodalPressureField_Delete( void* _self ) {
   NodalPressureField* self = (NodalPressureField*) _self;

   _ParticleFeVariable_Delete( self );
}

/*
** Methods.
*/

void _NodalPressureField_Print( void* _self, Stream* stream ) {
   NodalPressureField* self = (NodalPressureField*)_self;

   /* General info */
   Journal_Printf( stream, "NodalPressureField (ptr): %p\n", self );

   /* Print parent */
   _ParticleFeVariable_Print( self, stream );
}

void* _NodalPressureField_Copy( void* _self, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
   abort();
}

void NodalPressureField_NonLinearUpdate( void* _sle, void* _ctx ) {
   SystemLinearEquations*	sle = (SystemLinearEquations*)_sle;
   DomainContext*				ctx = (DomainContext*)_ctx;
   FieldVariable_Register*	fieldVar_Register;
   NodalPressureField*		preVar;

   fieldVar_Register = ctx->fieldVariable_Register; 
   preVar = (NodalPressureField*)FieldVariable_Register_GetByName( fieldVar_Register, "NodalPressureField" );
   ParticleFeVariable_Update( preVar );
}

void _NodalPressureField_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data ){
   NodalPressureField*		self = (NodalPressureField*) _self;
   Variable_Register*		variable_Register;
   SystemLinearEquations*	sle;
   FeVariable*					pressureField;

   _ParticleFeVariable_AssignFromXML( self, cf, data );

   variable_Register = self->variable_Register; 
   assert( variable_Register );

   pressureField = Stg_ComponentFactory_ConstructByKey( cf, self->name, "PressureField", FeVariable, True, data );

   _NodalPressureField_Init( self, variable_Register, pressureField );

   /*
   ** If we're using this field for non-linear feedback, we'll need to update it in between
   ** non-linear iterations. */
   sle = Stg_ComponentFactory_ConstructByKey( cf, self->name, "SLE", SystemLinearEquations, False, data );
   if( sle )
      SystemLinearEquations_AddPostNonLinearEP( sle, NodalPressureField_Type, NodalPressureField_NonLinearUpdate );
}

void _NodalPressureField_Build( void* _self, void* data ) {
   NodalPressureField* self = (NodalPressureField*) _self;
   Name tmpName, tmpName2;
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
   self->dofLayout = DofLayout_New( tmpName2, self->variable_Register, 0, self->feMesh );
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

void _NodalPressureField_Initialise( void* _self, void* data ) {
   NodalPressureField* self = (NodalPressureField*) _self;
   Variable_Index variable_I;
	
   /* Initialise and Update all Variables that this component has created */
   Stg_Component_Initialise( self->dataVariable, data, False); Variable_Update( self->dataVariable );
	
   _ParticleFeVariable_Initialise( self, data );

   /* Do a post-update just in case */
   Variable_Update( self->dataVariable );

}

void _NodalPressureField_Execute( void* _self, void* data ) {
   NodalPressureField* self = (NodalPressureField*) _self;

   _ParticleFeVariable_Execute( self, data );
}

void _NodalPressureField_Destroy( void* _self, void* data ) {
   NodalPressureField* self = (NodalPressureField*) _self;

   _ParticleFeVariable_Destroy( self, data );
}

void _NodalPressureField_ValueAtParticle( void* _self, IntegrationPointsSwarm* swarm,
					  Element_LocalIndex lElement_I, void* _particle,
					  double* pressure )
{
   NodalPressureField* self = (NodalPressureField*)_self;
   IntegrationPoint* particle = (IntegrationPoint*)_particle;

   FeVariable_InterpolateWithinElement( self->pressureField, lElement_I, particle->xi, pressure );
}
