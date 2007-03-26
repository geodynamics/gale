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
** $Id:  $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>

#include "types.h"
#include "EscapedRoutine.h"

#include <assert.h>
#include <string.h>
#include <math.h>

/* Textual name of this class */
const Type EscapedRoutine_Type = "EscapedRoutine";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

EscapedRoutine* _EscapedRoutine_New(
		SizeT                                      _sizeOfSelf, 
		Type                                       type,
		Stg_Class_DeleteFunction*                  _delete,
		Stg_Class_PrintFunction*                   _print,
		Stg_Class_CopyFunction*                    _copy, 
		Stg_Component_DefaultConstructorFunction*  _defaultConstructor,
		Stg_Component_ConstructFunction*           _construct,
		Stg_Component_BuildFunction*               _build,
		Stg_Component_InitialiseFunction*          _initialise,
		Stg_Component_ExecuteFunction*             _execute,
		Stg_Component_DestroyFunction*             _destroy,		
		EscapedRoutine_SelectFunction*     	   _select,
		Name                                       name )
{
	EscapedRoutine* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(EscapedRoutine) );
	self = (EscapedRoutine*)_Stg_Component_New( 
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
			name ,
			NON_GLOBAL );
	
	/* General info */

	/* Virtual Info */
	self->_select = _select;
	
	return self;
}

void* _EscapedRoutine_DefaultNew( Name name ) {
	return (void*) _EscapedRoutine_New(
			sizeof(EscapedRoutine),
			EscapedRoutine_Type,
			_EscapedRoutine_Delete,
			_EscapedRoutine_Print,
			_EscapedRoutine_Copy,
			_EscapedRoutine_DefaultNew,
			_EscapedRoutine_Construct,
			_EscapedRoutine_Build,
			_EscapedRoutine_Initialise,
			_EscapedRoutine_Execute,
			_EscapedRoutine_Destroy,
			_EscapedRoutine_Select, 
			name );
}

void _EscapedRoutine_Init( 
		EscapedRoutine*                            self, 
		Dimension_Index                            dim, 
		Particle_Index                             particlesToRemoveDelta )
{
	self->isConstructed          = True;
	self->dim                    = dim;
	self->particlesToRemoveDelta = particlesToRemoveDelta;

	self->debug = Journal_Register( Debug_Type, EscapedRoutine_Type ); /* TODO Register Child */
}

/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _EscapedRoutine_Delete( void* escapedRoutine ) {
	EscapedRoutine* self = (EscapedRoutine*)escapedRoutine;
	
	Memory_Free( self->particlesToRemoveList );

	/* Delete parent */
	_Stg_Component_Delete( self );
}


void _EscapedRoutine_Print( void* escapedRoutine, Stream* stream ) {
	EscapedRoutine* self = (EscapedRoutine*)escapedRoutine;
	
	/* Print parent */
	_Stg_Component_Print( self, stream );
}


void* _EscapedRoutine_Copy( void* escapedRoutine, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	EscapedRoutine*	self = (EscapedRoutine*)escapedRoutine;
	EscapedRoutine*	newEscapedRoutine;
	
	newEscapedRoutine = (EscapedRoutine*)_Stg_Component_Copy( self, dest, deep, nameExt, ptrMap );
	
	return (void*)newEscapedRoutine;
}

void _EscapedRoutine_Construct( void* escapedRoutine, Stg_ComponentFactory* cf, void* data ) {
	EscapedRoutine*	     self          = (EscapedRoutine*) escapedRoutine;
	Dimension_Index      dim;
	Particle_Index       particlesToRemoveDelta;

	dim = Stg_ComponentFactory_GetRootDictUnsignedInt( cf, "dim", 0 );
	particlesToRemoveDelta = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "particlesToRemoveDelta", 20 );

	_EscapedRoutine_Init( self, dim, particlesToRemoveDelta );
}

void _EscapedRoutine_Build( void* escapedRoutine, void* data ) {
	EscapedRoutine*	     self          = (EscapedRoutine*) escapedRoutine;

	self->particlesToRemoveList = Memory_Alloc_Array( unsigned, self->particlesToRemoveDelta * 10, "particlesToRemoveList" );
	
}
void _EscapedRoutine_Initialise( void* escapedRoutine, void* data ) {
}
void _EscapedRoutine_Execute( void* escapedRoutine, void* data ) {
	EscapedRoutine_RemoveFromSwarm( escapedRoutine, data );
}

void _EscapedRoutine_Destroy( void* escapedRoutine, void* data ) {
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void EscapedRoutine_Select( void* escapedRoutine, void* _swarm ) {
	EscapedRoutine*	 self          = (EscapedRoutine*) escapedRoutine;

	self->_select( self, _swarm );
}


void _EscapedRoutine_Select( void* escapedRoutine, void* _swarm ) {
	EscapedRoutine*	self = (EscapedRoutine*)escapedRoutine;
	Swarm*		swarm = (Swarm*)_swarm;
	unsigned	p_i;

	assert( self );
	assert( swarm );

	/* Check all particles for removal. */
	for( p_i = 0; p_i < swarm->particleLocalCount; p_i++ ) {
		StandardParticle*	particle = (StandardParticle*)Swarm_ParticleAt( swarm, p_i );

		if( particle->owningCell >= swarm->cellDomainCount ) {
			EscapedRoutine_SetParticleToRemove( self, swarm, p_i );
		}
	}
}


void EscapedRoutine_RemoveFromSwarm( void* escapedRoutine, void* _swarm ) {
	EscapedRoutine*	     self                = (EscapedRoutine*) escapedRoutine;
	Swarm*               swarm               = (Swarm*) _swarm;
	
	EscapedRoutine_InitialiseParticleList( self );
	
	/* Select particles to remove. */
	EscapedRoutine_Select( self, swarm );

	/* Actually remove particles */
	EscapedRoutine_RemoveParticles( self, swarm );
}

void EscapedRoutine_InitialiseParticleList( void* escapedRoutine ) {
	EscapedRoutine*	     self                = (EscapedRoutine*) escapedRoutine;

	self->particlesToRemoveCount = 0;
	memset( self->particlesToRemoveList, 0, sizeof(unsigned) * self->particlesToRemoveAlloced );
}

void EscapedRoutine_SetParticleToRemove( void* escapedRoutine, Swarm* swarm, Particle_Index lParticle_I ) {
	EscapedRoutine*	      self                = (EscapedRoutine*) escapedRoutine;

	/* Check memory */
	if ( self->particlesToRemoveCount >= self->particlesToRemoveAlloced ) {
		self->particlesToRemoveAlloced += self->particlesToRemoveDelta;
		self->particlesToRemoveList = 
			Memory_Realloc_Array( self->particlesToRemoveList, unsigned, self->particlesToRemoveAlloced );
	}

	self->particlesToRemoveList[ self->particlesToRemoveCount ] = lParticle_I;

	self->particlesToRemoveCount++;
}

int _EscapedRoutine_SortParticles( const void* _aParticleInfo, const void* _bParticleInfo ) {
	return (*(unsigned*)_aParticleInfo - *(unsigned*)_bParticleInfo );
}

void EscapedRoutine_SortParticleList( void* escapedRoutine ) {
	EscapedRoutine*	     self                = (EscapedRoutine*) escapedRoutine;

	qsort( self->particlesToRemoveList, self->particlesToRemoveCount, 
			sizeof(unsigned), _EscapedRoutine_SortParticles );
}

void EscapedRoutine_RemoveParticles( void* escapedRoutine, Swarm* swarm ) {
	EscapedRoutine*	      self                = (EscapedRoutine*) escapedRoutine;
	Index                 array_I;
	StandardParticle*     particleToRemove;
	Particle_Index        particleToRemove_I;
	
	StandardParticle*     lastParticle;
	Cell_Index            lastParticle_CellIndex;
	Particle_Index        lastParticle_I;
	Particle_InCellIndex  lastParticle_IndexWithinCell;
	SizeT                 particleSize        = swarm->particleExtensionMgr->finalSize;

	#if DEBUG
	if ( Stream_IsPrintableLevel( self->debug, 2 ) ) {
		Journal_Printf( self->debug, "Particles to remove:\n{ " );
		for ( array_I = 0 ; array_I < self->particlesToRemoveCount - 1 ; array_I++ ) {
			Journal_Printf( self->debug, "%u, ", self->particlesToRemoveList[ array_I ] );
		}
		Journal_Printf( self->debug, "%u }\n", self->particlesToRemoveList[ array_I ] );
	}
	#endif


	EscapedRoutine_SortParticleList( self );

	#if DEBUG
	if ( Stream_IsPrintableLevel( self->debug, 2 ) ) {
		Journal_Printf( self->debug, "Particles to remove:\n{ " );
		for ( array_I = 0 ; array_I < self->particlesToRemoveCount - 1 ; array_I++ ) {
			Journal_Printf( self->debug, "%u, ", self->particlesToRemoveList[ array_I ] );
		}
		Journal_Printf( self->debug, "%u }\n", self->particlesToRemoveList[ array_I ] );
	}
	#endif

	for ( array_I = self->particlesToRemoveCount - 1 ; array_I < self->particlesToRemoveCount ; array_I-- ) {
		particleToRemove_I               = self->particlesToRemoveList[ array_I ];
		particleToRemove                 = Swarm_ParticleAt( swarm, particleToRemove_I );

		Journal_DPrintfL( self->debug, 2, "Removing particle %u\n", particleToRemove_I );
		
		/* Copy over particle with last particle in array - as long as it isn't the last one */
		lastParticle_I = swarm->particleLocalCount - 1;
		lastParticle   = Swarm_ParticleAt( swarm, lastParticle_I );
		if ( particleToRemove_I != lastParticle_I ) {
			/* Get last Particle information */
			lastParticle_CellIndex       = lastParticle->owningCell;
			lastParticle_IndexWithinCell = Swarm_GetParticleIndexWithinCell( swarm, lastParticle_CellIndex, lastParticle_I);

			Journal_DPrintfL( self->debug, 2, 
					"Copying over particle %u using last particle %u from cell %u (cell particle index - %u)\n", 
					particleToRemove_I, lastParticle_I, lastParticle_CellIndex, lastParticle_IndexWithinCell );

			/* Copy over particle */
			memcpy( particleToRemove, lastParticle, particleSize );
			
			/* Change value in cell particle table to point to new index in array */
			swarm->cellParticleTbl[lastParticle_CellIndex][ lastParticle_IndexWithinCell ] = particleToRemove_I;
		}

		/* Initialise memory to zero so it is clear that it's been deleted */
		memset( lastParticle, 0, particleSize );
		swarm->particleLocalCount--;
	}

	Swarm_Realloc( swarm );
}
