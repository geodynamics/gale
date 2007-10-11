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
** $Id: RemovalRoutine.c 518 2007-10-11 08:07:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "types.h"
#include "RemovalRoutine.h"

#include <assert.h>
#include <string.h>
#include <math.h>

/* Textual name of this class */
const Type RemovalRoutine_Type = "RemovalRoutine";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

RemovalRoutine* _RemovalRoutine_New(
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
		RemovalRoutine_RemoveFromCellFunction*     _removeFromCell,
		Name                                       name )
{
	RemovalRoutine* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(RemovalRoutine) );
	self = (RemovalRoutine*)_Stg_Component_New( 
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
	self->_removeFromCell = _removeFromCell;
	
	return self;
}

void _RemovalRoutine_Init( 
		RemovalRoutine*                            self, 
		Dimension_Index                            dim, 
		Particle_InCellIndex                       idealParticleCount,
		Particle_InCellIndex                       maxParticlesPerCell,
		Particle_Index                             particlesToRemoveDelta )
{
	self->isConstructed          = True;
	self->dim                    = dim;
	self->idealParticleCount     = idealParticleCount;
	self->maxParticlesPerCell    = maxParticlesPerCell;
	self->particlesToRemoveDelta = particlesToRemoveDelta;

	self->debug = Journal_Register( Debug_Type, RemovalRoutine_Type ); /* TODO Register Child */
}

/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _RemovalRoutine_Delete( void* removalRoutine ) {
	RemovalRoutine* self = (RemovalRoutine*)removalRoutine;
	
	Memory_Free( self->particlesToRemoveList );

	/* Delete parent */
	_Stg_Component_Delete( self );
}


void _RemovalRoutine_Print( void* removalRoutine, Stream* stream ) {
	RemovalRoutine* self = (RemovalRoutine*)removalRoutine;
	
	/* Print parent */
	_Stg_Component_Print( self, stream );
}


void* _RemovalRoutine_Copy( void* removalRoutine, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	RemovalRoutine*	self = (RemovalRoutine*)removalRoutine;
	RemovalRoutine*	newRemovalRoutine;
	
	newRemovalRoutine = (RemovalRoutine*)_Stg_Component_Copy( self, dest, deep, nameExt, ptrMap );
	
	return (void*)newRemovalRoutine;
}

void _RemovalRoutine_Construct( void* removalRoutine, Stg_ComponentFactory* cf, void* data ) {
	RemovalRoutine*	     self          = (RemovalRoutine*) removalRoutine;
	Dimension_Index      dim;
	Particle_InCellIndex idealParticleCount;
	Particle_InCellIndex maxParticlesPerCell;
	Particle_Index       particlesToRemoveDelta;

	dim = Stg_ComponentFactory_GetRootDictUnsignedInt( cf, "dim", 0 );

	idealParticleCount  = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "idealParticleCount",  0 );
	maxParticlesPerCell = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "maxParticlesPerCell", idealParticleCount );
	
	particlesToRemoveDelta = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "particlesToRemoveDelta", 20 );

	_RemovalRoutine_Init( self, dim, idealParticleCount, maxParticlesPerCell, particlesToRemoveDelta );
}

void _RemovalRoutine_Build( void* removalRoutine, void* data ) {
	RemovalRoutine*	     self          = (RemovalRoutine*) removalRoutine;

	self->particlesToRemoveList = Memory_Alloc_Array( ParticleToRemoveInfo, self->particlesToRemoveDelta * 10, "particlesToRemoveList" );
	
}
void _RemovalRoutine_Initialise( void* removalRoutine, void* data ) {
}
void _RemovalRoutine_Execute( void* removalRoutine, void* data ) {
}
void _RemovalRoutine_Destroy( void* removalRoutine, void* data ) {
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void RemovalRoutine_RemoveFromCell( void* removalRoutine, void* _swarm, Cell_LocalIndex lCell_I ) {
	RemovalRoutine*	 self          = (RemovalRoutine*) removalRoutine;

	self->_removeFromCell( self, _swarm, lCell_I );
}

void RemovalRoutine_RemoveFromSwarm( void* removalRoutine, void* _swarm ) {
	RemovalRoutine*	     self                = (RemovalRoutine*) removalRoutine;
	Swarm*               swarm               = (Swarm*) _swarm;
	Particle_InCellIndex maxParticlesPerCell = self->maxParticlesPerCell;
	Cell_LocalIndex	     cellLocalCount      = swarm->cellLocalCount;
	Cell_LocalIndex	     lCell_I;
	
	RemovalRoutine_InitialiseParticleList( self );
	
	/* Loop over all local cells */
	for ( lCell_I = 0 ; lCell_I < cellLocalCount ; lCell_I++ ) {
		if ( swarm->cellParticleCountTbl[ lCell_I ] > maxParticlesPerCell ) 
			RemovalRoutine_RemoveFromCell( self, swarm, lCell_I );
	}

	/* Actually remove particles */
	RemovalRoutine_RemoveParticles( self, swarm );
}

void RemovalRoutine_InitialiseParticleList( void* removalRoutine ) {
	RemovalRoutine*	     self                = (RemovalRoutine*) removalRoutine;

	self->particlesToRemoveCount = 0;
	memset( self->particlesToRemoveList, 0, sizeof(Particle_Index) * self->particlesToRemoveAlloced );
}

void RemovalRoutine_SetParticleToRemove( void* removalRoutine, Swarm* swarm, Cell_Index cell_I, Particle_InCellIndex cParticle_I ) {
	RemovalRoutine*	      self                = (RemovalRoutine*) removalRoutine;
	ParticleToRemoveInfo* particleRemoveInfo;

	/* Check memory */
	if ( self->particlesToRemoveCount >= self->particlesToRemoveAlloced ) {
		self->particlesToRemoveAlloced += self->particlesToRemoveDelta;
		self->particlesToRemoveList = 
			Memory_Realloc_Array( self->particlesToRemoveList, ParticleToRemoveInfo, self->particlesToRemoveAlloced );
	}

	particleRemoveInfo = &self->particlesToRemoveList[ self->particlesToRemoveCount ];

	particleRemoveInfo->cell_I       = cell_I;
	particleRemoveInfo->lParticle_I  = swarm->cellParticleTbl[cell_I][cParticle_I];

	self->particlesToRemoveCount++;
}

int _RemovalRoutine_SortParticles( const void* _aParticleInfo, const void* _bParticleInfo ) {
	ParticleToRemoveInfo* aParticleInfo = (ParticleToRemoveInfo*) _aParticleInfo;
	ParticleToRemoveInfo* bParticleInfo = (ParticleToRemoveInfo*) _bParticleInfo;

	return ((int) aParticleInfo->lParticle_I - (int) bParticleInfo->lParticle_I );
}

void RemovalRoutine_SortParticleList( void* removalRoutine ) {
	RemovalRoutine*	     self                = (RemovalRoutine*) removalRoutine;

	qsort( self->particlesToRemoveList, self->particlesToRemoveCount, 
			sizeof(ParticleToRemoveInfo), _RemovalRoutine_SortParticles );
}

void RemovalRoutine_RemoveParticles( void* removalRoutine, Swarm* swarm ) {
	RemovalRoutine*	      self                = (RemovalRoutine*) removalRoutine;
	Index                 array_I;
	GlobalParticle*       particleToRemove;
	Particle_Index        particleToRemove_I;
	Particle_InCellIndex  particleToRemove_IndexWithinCell;
	Cell_Index            particleToRemove_CellIndex;
	ParticleToRemoveInfo* particleInfo;
	
	GlobalParticle*       lastParticle;
	Cell_Index            lastParticle_CellIndex;
	Particle_Index        lastParticle_I;
	Particle_InCellIndex  lastParticle_IndexWithinCell;
	SizeT                 particleSize        = swarm->particleExtensionMgr->finalSize;

	#if DEBUG
	if ( Stream_IsPrintableLevel( self->debug, 2 ) ) {
		Journal_Printf( self->debug, "Particles to remove:\n{ " );
		for ( array_I = 0 ; array_I < self->particlesToRemoveCount - 1 ; array_I++ ) {
			Journal_Printf( self->debug, "%u, ", self->particlesToRemoveList[ array_I ].lParticle_I );
		}
		Journal_Printf( self->debug, "%u }\n", self->particlesToRemoveList[ array_I ].lParticle_I );
	}
	#endif


	RemovalRoutine_SortParticleList( self );

	#if DEBUG
	if ( Stream_IsPrintableLevel( self->debug, 2 ) ) {
		Journal_Printf( self->debug, "Particles to remove:\n{ " );
		for ( array_I = 0 ; array_I < self->particlesToRemoveCount - 1 ; array_I++ ) {
			Journal_Printf( self->debug, "%u, ", self->particlesToRemoveList[ array_I ].lParticle_I );
		}
		Journal_Printf( self->debug, "%u }\n", self->particlesToRemoveList[ array_I ].lParticle_I );
	}
	#endif

	for ( array_I = self->particlesToRemoveCount - 1 ; array_I < self->particlesToRemoveCount ; array_I-- ) {
		particleInfo                     = &self->particlesToRemoveList[ array_I ];
		particleToRemove_I               = particleInfo->lParticle_I;
		particleToRemove_CellIndex       = particleInfo->cell_I;
		particleToRemove_IndexWithinCell = Swarm_GetParticleIndexWithinCell( 
				swarm, particleToRemove_CellIndex, particleToRemove_I );
		particleToRemove                 = (GlobalParticle*)Swarm_ParticleAt( swarm, particleToRemove_I );

		Journal_DPrintfL( self->debug, 2, "Removing particle %u from cell %u (cell particle index - %u)\n",
				particleToRemove_I, particleToRemove_CellIndex, particleToRemove_IndexWithinCell );
		
		Swarm_RemoveParticleFromCell( swarm, particleToRemove_CellIndex, particleToRemove_IndexWithinCell );
		
		/* Copy over particle with last particle in array - as long as it isn't the last one */
		lastParticle_I = swarm->particleLocalCount - 1;
		lastParticle   = (GlobalParticle*)Swarm_ParticleAt( swarm, lastParticle_I );
		if ( particleToRemove_I != lastParticle_I ) {
			/* Get last Particle information */
			lastParticle                 = (GlobalParticle*)Swarm_ParticleAt( swarm, lastParticle_I );
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
