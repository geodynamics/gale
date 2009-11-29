/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
** Copyright (c) 2005-2006, Monash Cluster Computing, Building 28, Monash University Clayton Campus,
**	Victoria, 3800, Australia
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**
** Contributors:
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Patrick D. Sunter, Software Engineer, VPAC. (patrick@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	David Stegman, Postdoctoral Fellow, Monash University. (david.stegman@sci.monash.edu.au)
**	Wendy Sharples, PhD Student, Monash University (wendy.sharples@sci.monash.edu.au)
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
** $Id: MaterialSwarmVariable.c 518 2007-10-11 08:07:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>
#include <PICellerator/MaterialPoints/MaterialPoints.h>

#include "types.h"
#include "MaterialSwarmVariable.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class */
const Type MaterialSwarmVariable_Type = "MaterialSwarmVariable";

MaterialSwarmVariable* MaterialSwarmVariable_New( 
	Name						name,
	AbstractContext*		context,
	MaterialPointsSwarm*	swarm,
	Index						dofCount,
	Materials_Register*	materials_Register,
	ExtensionInfo_Index	materialExtensionHandle,
	SizeT						offset )
{
	MaterialSwarmVariable* self = (MaterialSwarmVariable*) _MaterialSwarmVariable_DefaultNew( name );

	self->isConstructed = True;
	_SwarmVariable_Init( (SwarmVariable*)self, context, (Swarm*)swarm, NULL, dofCount );
	_MaterialSwarmVariable_Init( self, materials_Register, materialExtensionHandle, offset );

	return self;
}

void* _MaterialSwarmVariable_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                 _sizeOfSelf = sizeof(MaterialSwarmVariable);
	Type                                                         type = MaterialSwarmVariable_Type;
	Stg_Class_DeleteFunction*                                 _delete = _MaterialSwarmVariable_Delete;
	Stg_Class_PrintFunction*                                   _print = _MaterialSwarmVariable_Print;
	Stg_Class_CopyFunction*                                     _copy = NULL;
	Stg_Component_DefaultConstructorFunction*     _defaultConstructor = _MaterialSwarmVariable_DefaultNew;
	Stg_Component_ConstructFunction*                       _construct = _MaterialSwarmVariable_AssignFromXML;
	Stg_Component_BuildFunction*                               _build = _MaterialSwarmVariable_Build;
	Stg_Component_InitialiseFunction*                     _initialise = _MaterialSwarmVariable_Initialise;
	Stg_Component_ExecuteFunction*                           _execute = _MaterialSwarmVariable_Execute;
	Stg_Component_DestroyFunction*                           _destroy = _MaterialSwarmVariable_Destroy;
	AllocationType                                 nameAllocationType = NON_GLOBAL;
	SwarmVariable_ValueAtFunction*                           _valueAt = _MaterialSwarmVariable_ValueAt;
	SwarmVariable_GetGlobalValueFunction*      _getMinGlobalMagnitude = _MaterialSwarmVariable_GetMinGlobalMagnitude;
	SwarmVariable_GetGlobalValueFunction*      _getMaxGlobalMagnitude = _MaterialSwarmVariable_GetMaxGlobalMagnitude;

	return (void*)_MaterialSwarmVariable_New(  MATERIALSWARMVARIABLE_PASSARGS  );
}

/* Creation implementation / Virtual constructor */
MaterialSwarmVariable* _MaterialSwarmVariable_New(  MATERIALSWARMVARIABLE_DEFARGS  ) {
	MaterialSwarmVariable* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(MaterialSwarmVariable) );
	self = (MaterialSwarmVariable*) _SwarmVariable_New(  SWARMVARIABLE_PASSARGS  );
	
	/* Virtual info */
	
	return self;
}

void _MaterialSwarmVariable_Init(
	void*						swarmVariable,
	Materials_Register*	materials_Register,
	ExtensionInfo_Index	materialExtensionHandle,
	SizeT						offset )
{
	MaterialSwarmVariable* self = (MaterialSwarmVariable*)swarmVariable;

	self->materials_Register        = materials_Register;
	self->materialExtensionHandle   = materialExtensionHandle;
	self->offset                    = offset;
   /* variable does not store data, so is not checkpointed */
   self->isCheckpointedAndReloaded = False;
}

void _MaterialSwarmVariable_Delete( void* swarmVariable ) {
	MaterialSwarmVariable* self = (MaterialSwarmVariable*)swarmVariable;

	_SwarmVariable_Delete( self );
}

void _MaterialSwarmVariable_Print( void* swarmVariable, Stream* stream ) {
	MaterialSwarmVariable* self = (MaterialSwarmVariable*)swarmVariable;
	
	_SwarmVariable_Print( self, stream );

	/* General info */
	Journal_PrintPointer( stream, self->materials_Register );
	Journal_PrintValue( stream, self->materialExtensionHandle );
	Journal_PrintValue( stream, self->offset );
}

void _MaterialSwarmVariable_AssignFromXML( void* swarmVariable, Stg_ComponentFactory* cf, void* data ) {
	MaterialSwarmVariable*	self = (MaterialSwarmVariable*)swarmVariable;
	Materials_Register*		materials_Register;
	PICelleratorContext*		context = (PICelleratorContext*)self->context;

	/* Construct Parent */
	_SwarmVariable_AssignFromXML( self, cf, data );

	assert( Stg_CheckType( context, PICelleratorContext ) );
	materials_Register = context->materials_Register;
	assert( materials_Register );

	abort( );
	_MaterialSwarmVariable_Init( self, materials_Register, 0, 0 );
}

void _MaterialSwarmVariable_Build( void* swarmVariable, void* data ) {
	MaterialSwarmVariable* self = (MaterialSwarmVariable*)swarmVariable;

	_SwarmVariable_Build( self, data );
}

void _MaterialSwarmVariable_Initialise( void* swarmVariable, void* data ) {
	_SwarmVariable_Initialise( swarmVariable, data );
}

void _MaterialSwarmVariable_Execute( void* swarmVariable, void* data ) {
	_SwarmVariable_Execute( swarmVariable, data );
}

void _MaterialSwarmVariable_Destroy( void* swarmVariable, void* data ) {
	MaterialSwarmVariable*  self = (MaterialSwarmVariable*) swarmVariable;

	_SwarmVariable_Destroy( self, data );
}

void _MaterialSwarmVariable_ValueAt( void* swarmVariable, Particle_Index lParticle_I, double* value ) {
	MaterialSwarmVariable*  self = (MaterialSwarmVariable*) swarmVariable;
	MaterialPointsSwarm*    swarm = (MaterialPointsSwarm*)self->swarm;
	Material*               material;
	ArithPointer            materialExt;

	/* Get Material */
	material = MaterialPointsSwarm_GetMaterialAt( swarm, lParticle_I );

	/* Get Material Extension */
	materialExt = (ArithPointer) ExtensionManager_Get( material->extensionMgr, NULL, self->materialExtensionHandle );

	/* Copy value */
	memcpy( value, (void*) (materialExt + self->offset), sizeof(double) * self->dofCount );
}

double _MaterialSwarmVariable_GetMinGlobalMagnitude( void* swarmVariable ) {
	return _SwarmVariable_GetMinGlobalMagnitude( swarmVariable );
}
double _MaterialSwarmVariable_GetMaxGlobalMagnitude( void* swarmVariable ) {
	return _SwarmVariable_GetMaxGlobalMagnitude( swarmVariable );
}
	
	


