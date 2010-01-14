/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Lesser General Public
**  License as published by the Free Software Foundation; either
**  version 2.1 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License along with this library; if not, write to the Free Software
**  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**
** $Id: SwarmVariable.c 4139 2007-06-12 02:39:52Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>

#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>
#include <StgDomain/Utils/Utils.h>

#include "types.h"
#include "SwarmClass.h"
#include "SwarmVariable_Register.h"
#include "SwarmVariable.h"

#include <assert.h>
#include <string.h>

const Type SwarmVariable_Type = "SwarmVariable";

SwarmVariable* SwarmVariable_New(		
	Name					name,
	AbstractContext*	context,
	Swarm*				swarm,
	Variable*			variable,
	Index					dofCount )
{
	SwarmVariable* self = (SwarmVariable*) _SwarmVariable_DefaultNew( name );

	self->isConstructed = True;
	_SwarmVariable_Init( self, context, swarm, variable, dofCount );

	return self;
}

void* _SwarmVariable_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                 _sizeOfSelf = sizeof(SwarmVariable);
	Type                                                         type = SwarmVariable_Type;
	Stg_Class_DeleteFunction*                                 _delete = _SwarmVariable_Delete;
	Stg_Class_PrintFunction*                                   _print = _SwarmVariable_Print;
	Stg_Class_CopyFunction*                                     _copy = _SwarmVariable_Copy;
	Stg_Component_DefaultConstructorFunction*     _defaultConstructor = _SwarmVariable_DefaultNew;
	Stg_Component_ConstructFunction*                       _construct = _SwarmVariable_AssignFromXML;
	Stg_Component_BuildFunction*                               _build = _SwarmVariable_Build;
	Stg_Component_InitialiseFunction*                     _initialise = _SwarmVariable_Initialise;
	Stg_Component_ExecuteFunction*                           _execute = _SwarmVariable_Execute;
	Stg_Component_DestroyFunction*                           _destroy = _SwarmVariable_Destroy;
	AllocationType                                 nameAllocationType = NON_GLOBAL;
	SwarmVariable_ValueAtFunction*                           _valueAt = _SwarmVariable_ValueAt;
	SwarmVariable_GetGlobalValueFunction*      _getMinGlobalMagnitude = _SwarmVariable_GetMinGlobalMagnitude;
	SwarmVariable_GetGlobalValueFunction*      _getMaxGlobalMagnitude = _SwarmVariable_GetMaxGlobalMagnitude;

	return (void*) _SwarmVariable_New(  SWARMVARIABLE_PASSARGS  );
}

SwarmVariable* _SwarmVariable_New(  SWARMVARIABLE_DEFARGS  ) {
	SwarmVariable* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(SwarmVariable) );
	self = (SwarmVariable*)_Stg_Component_New(  STG_COMPONENT_PASSARGS  );
	
	/* Virtual functions */
	self->_valueAt						= _valueAt;
	self->_getMinGlobalMagnitude	= _getMinGlobalMagnitude;
	self->_getMaxGlobalMagnitude	= _getMaxGlobalMagnitude;

	return self;
}

void _SwarmVariable_Init( void* swarmVariable, AbstractContext* context, Swarm* swarm, Variable* variable, Index dofCount ) {
	SwarmVariable* self = (SwarmVariable*) swarmVariable;

	self->context                   = context;
	self->swarm                     = swarm;
	self->variable                  = variable;
	self->dofCount                  = dofCount;
	self->swarmVariable_Register    = swarm->swarmVariable_Register;
	self->dim                       = swarm->dim;
   self->isCheckpointedAndReloaded = True;
	
	if( variable )
		Swarm_AddVariable( swarm, self );
	if ( self->swarmVariable_Register != NULL )	
		SwarmVariable_Register_Add( self->swarmVariable_Register, self );
}

void _SwarmVariable_Delete( void* swarmVariable ) {
	SwarmVariable* self = (SwarmVariable*) swarmVariable;

	_Stg_Component_Delete( self );
}

void _SwarmVariable_Print( void* _swarmVariable, Stream* stream ) {
	SwarmVariable* self = (SwarmVariable*) _swarmVariable;

	Journal_Printf( stream, "SwarmVariable - '%s'\n", self->name );
	Stream_Indent( stream );
	_Stg_Component_Print( self, stream );

	Journal_PrintPointer( stream, self->_valueAt );
	Journal_PrintPointer( stream, self->_getMinGlobalMagnitude );
	Journal_PrintPointer( stream, self->_getMaxGlobalMagnitude );

	Journal_Printf( stream, "Swarm - '%s'\n", self->swarm->name );
	if ( self->variable != NULL )
		Journal_Printf( stream, "Variable - '%s'\n", self->variable->name );

	Journal_PrintValue( stream, self->dofCount );
	Journal_PrintPointer( stream, self->swarmVariable_Register );
	Stream_UnIndent( stream );
}

void* _SwarmVariable_Copy( void* swarmVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	SwarmVariable*	self = (SwarmVariable*)swarmVariable;
	SwarmVariable*	newSwarmVariable;
	PtrMap*			map = ptrMap;
	Bool				ownMap = False;
	
	if( !map ) {
		map = PtrMap_New( 10 );
		ownMap = True;
	}
	
	newSwarmVariable = _Stg_Component_Copy( self, dest, deep, nameExt, map );
	
	newSwarmVariable->_valueAt						= self->_valueAt;
	newSwarmVariable->_getMinGlobalMagnitude	= self->_getMinGlobalMagnitude;
	newSwarmVariable->_getMaxGlobalMagnitude	= self->_getMaxGlobalMagnitude;

	newSwarmVariable->swarm							= self->swarm;
	newSwarmVariable->variable						= self->variable;
	newSwarmVariable->dofCount						= self->dofCount;
	newSwarmVariable->swarmVariable_Register	= self->swarmVariable_Register;

	if( ownMap ) {
		Stg_Class_Delete( map );
	}
				
	return (void*)newSwarmVariable;
}

void _SwarmVariable_AssignFromXML( void* swarmVariable, Stg_ComponentFactory* cf, void* data ) {
	SwarmVariable*		self = (SwarmVariable*)swarmVariable;
	Swarm*				swarm;
	Variable*			variable;
	Index					dofCount;
	AbstractContext*	context;

	context = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Context", AbstractContext, False, data );
	if( !context  )
		context = Stg_ComponentFactory_ConstructByName( cf, (Name)"context", AbstractContext, True, data  );

	swarm = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Swarm", Swarm, True, data  );
	variable =  Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Variable", Variable, False, data  );
	dofCount = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"dofCount", 0  );

	_SwarmVariable_Init( self, context, swarm, variable, dofCount );
}

void _SwarmVariable_Build( void* swarmVariable, void* data ) {
	SwarmVariable* self = (SwarmVariable*)swarmVariable;

	if ( self->variable )
		Stg_Component_Build( self->variable, data, False );
}

void _SwarmVariable_Initialise( void* swarmVariable, void* data ) {
	SwarmVariable* self = (SwarmVariable*)swarmVariable;

	if ( self->variable ) {
		Variable_Update( self->variable );
		Stg_Component_Initialise( self->variable, data, False );
	}
}

void _SwarmVariable_Execute( void* swarmVariable, void* data ) {
}

void _SwarmVariable_Destroy( void* swarmVariable, void* data ) {
	SwarmVariable* self = (SwarmVariable*)swarmVariable;
   
   Stg_Component_Destroy( self->swarm, data, False );
   if( self->variable )
		Stg_Component_Destroy(self->variable, data, False);
}

double SwarmVariable_GetMinGlobalMagnitude( void* swarmVariable ) {
	SwarmVariable*	self = (SwarmVariable*)swarmVariable;
	return self->_getMinGlobalMagnitude( self );
}

double SwarmVariable_GetMaxGlobalMagnitude( void* swarmVariable ) {
	SwarmVariable*	self = (SwarmVariable*)swarmVariable;
	return self->_getMaxGlobalMagnitude( self );
}

/*** Default Implementations ***/

void _SwarmVariable_ValueAt( void* swarmVariable, Particle_Index lParticle_I, double* value ) {
	SwarmVariable*	self = (SwarmVariable*)swarmVariable;
	Variable*		variable = self->variable;

	switch( variable->dataTypes[0] ) {
		case Variable_DataType_Double:
			self->_valueAt = _SwarmVariable_ValueAtDouble;
			break;
		case Variable_DataType_Int:
			self->_valueAt = _SwarmVariable_ValueAtInt;
			break;
		case Variable_DataType_Float:
			self->_valueAt = _SwarmVariable_ValueAtFloat;
			break;
		case Variable_DataType_Char:
			self->_valueAt = _SwarmVariable_ValueAtChar;
			break;
		case Variable_DataType_Short:
			self->_valueAt = _SwarmVariable_ValueAtShort;
			break;
		default:
			assert(0);
	}
	SwarmVariable_ValueAt( self, lParticle_I, value );
}
	
void _SwarmVariable_ValueAtDouble( void* swarmVariable, Particle_Index lParticle_I, double* value ) {
	SwarmVariable*	self = (SwarmVariable*)swarmVariable;
	double*			dataPtr;

	dataPtr = Variable_GetPtrDouble( self->variable, lParticle_I );
	memcpy( value, dataPtr, sizeof(double) * self->dofCount );
}

void _SwarmVariable_ValueAtInt( void* swarmVariable, Particle_Index lParticle_I, double* value ) {
	SwarmVariable*	self = (SwarmVariable*)swarmVariable;
	Variable*		variable = self->variable;
	Dof_Index		dofCount = self->dofCount;
	Dof_Index		dof_I;

	for ( dof_I = 0 ; dof_I < dofCount ; dof_I++ ) {
		value[ dof_I ] = (double) Variable_GetValueAtInt( variable, lParticle_I, dof_I );
	}
}

void _SwarmVariable_ValueAtFloat( void* swarmVariable, Particle_Index lParticle_I, double* value ) {
	SwarmVariable*	self = (SwarmVariable*)swarmVariable;
	Variable*		variable = self->variable;
	Dof_Index		dofCount = self->dofCount;
	Dof_Index		dof_I;

	for ( dof_I = 0 ; dof_I < dofCount ; dof_I++ ) {
		value[ dof_I ] = (double) Variable_GetValueAtFloat( variable, lParticle_I, dof_I );
	}
}

void _SwarmVariable_ValueAtChar( void* swarmVariable, Particle_Index lParticle_I, double* value ) {
	SwarmVariable*	self = (SwarmVariable*)swarmVariable;
	Variable*		variable = self->variable;
	Dof_Index		dofCount = self->dofCount;
	Dof_Index		dof_I;

	for ( dof_I = 0 ; dof_I < dofCount ; dof_I++ ) {
		value[ dof_I ] = (double) Variable_GetValueAtChar( variable, lParticle_I, dof_I );
	}
}

void _SwarmVariable_ValueAtShort( void* swarmVariable, Particle_Index lParticle_I, double* value ) {
	SwarmVariable*	self = (SwarmVariable*)swarmVariable;
	Variable*		variable = self->variable;
	Dof_Index		dofCount = self->dofCount;
	Dof_Index		dof_I;

	for ( dof_I = 0 ; dof_I < dofCount ; dof_I++ ) {
		value[ dof_I ] = (double) Variable_GetValueAtShort( variable, lParticle_I, dof_I );
	}
}
	
double _SwarmVariable_GetMinGlobalMagnitude( void* swarmVariable ) {
	SwarmVariable*	self = (SwarmVariable*)swarmVariable;
	double*			value;
	Swarm*			swarm = self->swarm;
	Particle_Index	particleLocalCount = swarm->particleLocalCount;
	Particle_Index	lParticle_I;
	double			localMin = HUGE_VAL;
	double			globalMin;
	Index				dofCount = self->dofCount;
	double			magnitude;

	value = Memory_Alloc_Array( double, dofCount, "value" );

	/* Search through all local particles and find smallest value of variable */
	for ( lParticle_I = 0 ; lParticle_I < particleLocalCount ; lParticle_I++ ) {
		SwarmVariable_ValueAt( self, lParticle_I, value );

		if ( dofCount == 1 )
			magnitude = value[0];
		else 
			assert(0); /*TODO */

		if ( localMin > magnitude )
			localMin = magnitude;
	}

	Memory_Free( value );
	MPI_Allreduce( &localMin, &globalMin, dofCount, MPI_DOUBLE, MPI_MIN, swarm->comm );

	return globalMin;
}


double _SwarmVariable_GetMaxGlobalMagnitude( void* swarmVariable ) {
	SwarmVariable*	self = (SwarmVariable*)swarmVariable;
	double*			value;
	Swarm*			swarm = self->swarm;
	Particle_Index	particleLocalCount = swarm->particleLocalCount;
	Particle_Index	lParticle_I;
	double			localMax = -HUGE_VAL;
	double			globalMax;
	Index				dofCount = self->dofCount;
	double			magnitude;

	value = Memory_Alloc_Array( double, dofCount, "value" );

	/* Search through all local particles and find smallest value of variable */
	for ( lParticle_I = 0 ; lParticle_I < particleLocalCount ; lParticle_I++ ) {
		SwarmVariable_ValueAt( self, lParticle_I, value );

		if ( dofCount == 1 )
			magnitude = value[0];
		else 
			assert(0); /*TODO */

		if ( localMax < magnitude )
			localMax = magnitude;
	}

	Memory_Free( value );
	MPI_Allreduce( &localMax, &globalMax, dofCount, MPI_DOUBLE, MPI_MAX, swarm->comm );
	
	return globalMax;
}

	






