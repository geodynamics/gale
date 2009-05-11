/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**	Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
** Contributors:
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
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
** $Id: FeSwarmVariable.c 1117 2008-05-02 02:50:02Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <string.h>
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "types.h"
#include "FeVariable.h"
#include "FeMesh.h"
#include "FeSwarmVariable.h"

#include <assert.h>

const Type FeSwarmVariable_Type = "FeSwarmVariable";

FeSwarmVariable* _FeSwarmVariable_New(
 		SizeT                                              _sizeOfSelf, 
		Type                                               type,
		Stg_Class_DeleteFunction*                          _delete,
		Stg_Class_PrintFunction*                           _print, 
		Stg_Class_CopyFunction*                            _copy, 
		Stg_Component_DefaultConstructorFunction*          _defaultConstructor,
		Stg_Component_ConstructFunction*                   _construct,
		Stg_Component_BuildFunction*                       _build,
		Stg_Component_InitialiseFunction*                  _initialise,
		Stg_Component_ExecuteFunction*                     _execute,
		Stg_Component_DestroyFunction*                     _destroy,
		SwarmVariable_ValueAtFunction*                     _valueAt,
		SwarmVariable_GetGlobalValueFunction*              _getMinGlobalMagnitude,
		SwarmVariable_GetGlobalValueFunction*              _getMaxGlobalMagnitude,
		Name                                               name )
{
	FeSwarmVariable*		    self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(FeSwarmVariable) );
	self = (FeSwarmVariable*) _SwarmVariable_New( 
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
			_valueAt,
			_getMinGlobalMagnitude,
			_getMaxGlobalMagnitude,
			name );

	return self;
}

void _FeSwarmVariable_Init( void* swarmVariable, FeVariable* feVariable ) {
	FeSwarmVariable*         self              = (FeSwarmVariable*) swarmVariable;

	self->feVariable                = feVariable;
   /* variable does not store data, so is not checkpointed */
   self->isCheckpointedAndReloaded = False;
}

void _FeSwarmVariable_Delete( void* _swarmVariable ) {
	FeSwarmVariable* self = (FeSwarmVariable*) _swarmVariable;

	_SwarmVariable_Delete( self );
}

void _FeSwarmVariable_Print( void* _swarmVariable, Stream* stream ) {
	FeSwarmVariable* self = (FeSwarmVariable*) _swarmVariable;

	_SwarmVariable_Print( self, stream );

	Journal_PrintPointer( stream, self->feVariable );
}


void* _FeSwarmVariable_Copy( void* swarmVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	FeSwarmVariable*	self = (FeSwarmVariable*)swarmVariable;
	FeSwarmVariable*	newFeSwarmVariable;
	
	newFeSwarmVariable = _SwarmVariable_Copy( self, dest, deep, nameExt, ptrMap );
	
	newFeSwarmVariable->feVariable              = self->feVariable;
	
	return (void*)newFeSwarmVariable;
}

void* _FeSwarmVariable_DefaultNew( Name name ) {
		return _FeSwarmVariable_New( 
			sizeof(FeSwarmVariable), 
			FeSwarmVariable_Type, 
			_FeSwarmVariable_Delete, 
			_FeSwarmVariable_Print,
			_FeSwarmVariable_Copy, 
			_FeSwarmVariable_DefaultNew,
			_FeSwarmVariable_Construct,
			_FeSwarmVariable_Build, 
			_FeSwarmVariable_Initialise, 
			_FeSwarmVariable_Execute,
			_FeSwarmVariable_Destroy,
			_FeSwarmVariable_ValueAt,
			_FeSwarmVariable_GetMinGlobalMagnitude,
			_FeSwarmVariable_GetMaxGlobalMagnitude,
			name );
}

void _FeSwarmVariable_Construct( void* swarmVariable, Stg_ComponentFactory* cf, void* data ) {
	FeSwarmVariable*        self       = (FeSwarmVariable*) swarmVariable;

	_SwarmVariable_Construct( self, cf, data );

	_FeSwarmVariable_Init( 
		self, 
		Stg_ComponentFactory_ConstructByKey( cf, self->name, "FeVariable", FeVariable, True, data ) ) ;
}

void _FeSwarmVariable_Build( void* swarmVariable, void* data ) {
	FeSwarmVariable* self = (FeSwarmVariable*) swarmVariable;

	Stg_Component_Build( self->feVariable, data, False );
}

void _FeSwarmVariable_Initialise( void* swarmVariable, void* data ) {
	FeSwarmVariable* self = (FeSwarmVariable*) swarmVariable;

	Stg_Component_Initialise( self->feVariable, data, False );
}

void _FeSwarmVariable_Execute( void* swarmVariable, void* data ) {}
void _FeSwarmVariable_Destroy( void* swarmVariable, void* data ) {}

void _FeSwarmVariable_ValueAt( void* swarmVariable, Particle_Index lParticle_I, double* value ) {
	FeSwarmVariable*  self            = (FeSwarmVariable*) swarmVariable;
	FeMesh*           mesh            = self->feVariable->feMesh;
	Swarm*            swarm           = self->swarm;
	GlobalParticle*   particle        = (GlobalParticle*) Swarm_ParticleAt( swarm, lParticle_I );
	double            xi[3];

	/* check if the swarm is using a GlobalCoordSystem, if not FAIL */
	Journal_Firewall( (swarm->particleLayout->coordSystem == GlobalCoordSystem), 
			Journal_Register( Error_Type, "FeSwarmVariable" ),
			"Error in function %s: FeSwarmVariables require swarms with GlobalCoordSystems\n"
		        "This FeSwarmVariable, %s, uses the swarm %s, which doesn't have a GlobalCoordSystem\n",
			__func__, self->name, swarm->name );
	FeMesh_CoordGlobalToLocal( mesh, particle->owningCell, particle->coord, xi );
	FeVariable_InterpolateWithinElement( self->feVariable, particle->owningCell, xi, value );
}


double _FeSwarmVariable_GetMinGlobalMagnitude( void* swarmVariable ) {
	FeSwarmVariable*  self            = (FeSwarmVariable*) swarmVariable;

	return FieldVariable_GetMinGlobalFieldMagnitude( self->feVariable );
}
double _FeSwarmVariable_GetMaxGlobalMagnitude( void* swarmVariable ) {
	FeSwarmVariable*  self            = (FeSwarmVariable*) swarmVariable;

	return FieldVariable_GetMaxGlobalFieldMagnitude( self->feVariable );
}
