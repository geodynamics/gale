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

NodalPressureField* _NodalPressureField_New(  NODALPRESSUREFIELD_DEFARGS  )
{
   NodalPressureField* self;

   assert( _sizeOfSelf >= sizeof(NodalPressureField) );
   self = (NodalPressureField*)_ParticleFeVariable_New(  PARTICLEFEVARIABLE_PASSARGS  );

   return self;
}

void _NodalPressureField_Init( NodalPressureField* self, Variable_Register* variable_Register, FeVariable* pressureField, SystemLinearEquations* sle) {
   self->variable_Register = variable_Register;
   self->fieldComponentCount = 1;
   self->pressureField = pressureField;
   if( sle )
      SystemLinearEquations_AddPostNonLinearEP( sle, NodalPressureField_Type, NodalPressureField_NonLinearUpdate );
   
}

void* _NodalPressureField_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                       _sizeOfSelf = sizeof(NodalPressureField);
	Type                                                               type = NodalPressureField_Type;
	Stg_Class_DeleteFunction*                                       _delete = _NodalPressureField_Delete;
	Stg_Class_PrintFunction*                                         _print = _NodalPressureField_Print;
	Stg_Class_CopyFunction*                                           _copy = _NodalPressureField_Copy;
	Stg_Component_DefaultConstructorFunction*           _defaultConstructor = _NodalPressureField_DefaultNew;
	Stg_Component_ConstructFunction*                             _construct = _NodalPressureField_AssignFromXML;
	Stg_Component_BuildFunction*                                     _build = _NodalPressureField_Build;
	Stg_Component_InitialiseFunction*                           _initialise = _NodalPressureField_Initialise;
	Stg_Component_ExecuteFunction*                                 _execute = _NodalPressureField_Execute;
	Stg_Component_DestroyFunction*                                 _destroy = _NodalPressureField_Destroy;
	FieldVariable_InterpolateValueAtFunction*           _interpolateValueAt = _FeVariable_InterpolateValueAt;
	FieldVariable_GetCoordFunction*                _getMinAndMaxLocalCoords = _FeVariable_GetMinAndMaxLocalCoords;
	FieldVariable_GetCoordFunction*               _getMinAndMaxGlobalCoords = _FeVariable_GetMinAndMaxGlobalCoords;
	FeVariable_InterpolateWithinElementFunction*  _interpolateWithinElement = _FeVariable_InterpolateNodeValuesToElLocalCoord;
	FeVariable_GetValueAtNodeFunction*                      _getValueAtNode = _FeVariable_GetValueAtNode;
	ParticleFeVariable_ValueAtParticleFunction*            _valueAtParticle = _NodalPressureField_ValueAtParticle;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType                             nameAllocationType = ZERO;
	FieldVariable_GetValueFunction*   _getMinGlobalFieldMagnitude = ZERO;
	FieldVariable_GetValueFunction*   _getMaxGlobalFieldMagnitude = ZERO;
	FeVariable_SyncShadowValuesFunc*            _syncShadowValues = ZERO;

   return (void*) _NodalPressureField_New(  NODALPRESSUREFIELD_PASSARGS  );
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

   /*
   ** If we're using this field for non-linear feedback, we'll need to update it in between
   ** non-linear iterations. */
   sle = Stg_ComponentFactory_ConstructByKey( cf, self->name, "SLE", SystemLinearEquations, False, data );

   _NodalPressureField_Init( self, variable_Register, pressureField, sle );

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
		(AbstractContext*)self->context,
      Variable_DataType_Double, 
      &((IGraph*)self->feMesh->topo)->remotes[MT_VERTEX]->nDomains, 
      NULL,
      (void**)&self->data, 
      self->variable_Register );
	
   /* Create Dof Layout */
   tmpName2 = Stg_Object_AppendSuffix( self, "DofLayout" );
   self->dofLayout = DofLayout_New( tmpName2, self->context, self->variable_Register, 0, self->feMesh );
   self->dofLayout->_numItemsInLayout = FeMesh_GetNodeDomainSize( self->feMesh );
   for( node_I = 0; node_I < FeMesh_GetNodeDomainSize( self->feMesh ); node_I++ ) {
      DofLayout_AddDof_ByVarName( self->dofLayout, tmpName, node_I );
   }
   Memory_Free( tmpName );
   Memory_Free( tmpName2 );
   self->eqNum->dofLayout = self->dofLayout;

   /* Build and Update all Variables that this component has created */
   Stg_Component_Build( self->dataVariable, data, False); Variable_Update( self->dataVariable );
   Stg_Component_Build( self->pressureField, data, False);

   _ParticleFeVariable_Build( self, data );
   /* Update again, just in case things were changed/reallocated when ICs loaded */
   Variable_Update( self->dataVariable );

}

void _NodalPressureField_Initialise( void* _self, void* data ) {
   NodalPressureField* self = (NodalPressureField*) _self;
   Variable_Index variable_I;
	
   /* Initialise and Update all Variables that this component has created */
   Stg_Component_Initialise( self->dataVariable, data, False); Variable_Update( self->dataVariable );
   
   Stg_Component_Initialise( self->pressureField, data, False);
	
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

   Stg_Component_Destroy( self->pressureField, data, False);

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


