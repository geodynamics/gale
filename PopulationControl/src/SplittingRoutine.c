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
** $Id: SplittingRoutine.c 376 2006-10-18 06:58:41Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>

#include "types.h"
#include "SplittingRoutine.h"

#include <assert.h>
#include <string.h>
#include <math.h>

/* Textual name of this class */
const Type SplittingRoutine_Type = "SplittingRoutine";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

SplittingRoutine* _SplittingRoutine_New(
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
		SplittingRoutine_SplitParticlesInCellFunction*     _splitParticlesInCell,
		Name                                               name )
{
	SplittingRoutine* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(SplittingRoutine) );
	self = (SplittingRoutine*)_Stg_Component_New( 
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
	self->_splitParticlesInCell = _splitParticlesInCell;
	
	return self;
}

void _SplittingRoutine_Init( 
		SplittingRoutine*                                  self, 
		Dimension_Index                                    dim, 
		Particle_InCellIndex                               idealParticleCount,
		Particle_InCellIndex                               minParticlesPerCell )
{
	self->isConstructed          = True;
	self->dim                    = dim;
	self->idealParticleCount     = idealParticleCount;
	self->minParticlesPerCell    = minParticlesPerCell;

	self->debug = Journal_Register( Debug_Type, SplittingRoutine_Type ); /* TODO Register Child */
}

/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _SplittingRoutine_Delete( void* splittingRoutine ) {
	SplittingRoutine* self = (SplittingRoutine*)splittingRoutine;
	
	/* Delete parent */
	_Stg_Component_Delete( self );
}


void _SplittingRoutine_Print( void* splittingRoutine, Stream* stream ) {
	SplittingRoutine* self = (SplittingRoutine*)splittingRoutine;
	
	/* Print parent */
	_Stg_Component_Print( self, stream );
}


void* _SplittingRoutine_Copy( void* splittingRoutine, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	SplittingRoutine*	self = (SplittingRoutine*)splittingRoutine;
	SplittingRoutine*	newSplittingRoutine;
	
	newSplittingRoutine = (SplittingRoutine*)_Stg_Component_Copy( self, dest, deep, nameExt, ptrMap );
	
	return (void*)newSplittingRoutine;
}

void _SplittingRoutine_Construct( void* splittingRoutine, Stg_ComponentFactory* cf, void* data ) {
	SplittingRoutine*	 self          = (SplittingRoutine*) splittingRoutine;
	Dimension_Index      dim;
	Particle_InCellIndex idealParticleCount;
	Particle_InCellIndex minParticlesPerCell;

	dim = Stg_ComponentFactory_GetRootDictUnsignedInt( cf, "dim", 0 );

	idealParticleCount  = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "idealParticleCount",  0 );
	minParticlesPerCell = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "minParticlesPerCell", idealParticleCount );
	
	_SplittingRoutine_Init( self, dim, idealParticleCount, minParticlesPerCell );
}

void _SplittingRoutine_Build( void* splittingRoutine, void* data ) {
}
void _SplittingRoutine_Initialise( void* splittingRoutine, void* data ) {
}
void _SplittingRoutine_Execute( void* splittingRoutine, void* _swarm ) {
	SplittingRoutine*    self                = (SplittingRoutine*) splittingRoutine;
	Swarm*               swarm               = (Swarm*) _swarm;
	Particle_InCellIndex minParticlesPerCell = self->minParticlesPerCell;
	Cell_LocalIndex	     cellLocalCount      = swarm->cellLocalCount;
	Cell_LocalIndex	     lCell_I;
	
	/* Loop over all local cells */
	for ( lCell_I = 0 ; lCell_I < cellLocalCount ; lCell_I++ ) {
		if ( swarm->cellParticleCountTbl[ lCell_I ] < minParticlesPerCell ) 
			SplittingRoutine_SplitParticlesInCell( self, swarm, lCell_I );
	}
}
void _SplittingRoutine_Destroy( void* splittingRoutine, void* data ) {
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void SplittingRoutine_SplitParticlesInCell( void* splittingRoutine, void* _swarm, Cell_LocalIndex lCell_I ) {
	SplittingRoutine*	 self          = (SplittingRoutine*) splittingRoutine;

	self->_splitParticlesInCell( self, _swarm, lCell_I );
}

	


void SplittingRoutine_SplitParticle( void* splittingRoutine, void* _swarm, Cell_LocalIndex lCell_I, Particle_InCellIndex cParticle_I, Coord newCoord ) 
{
	Swarm*               swarm            = (Swarm*)            _swarm;
	Particle_Index       newParticle_I;
	GlobalParticle*      newParticle;
	GlobalParticle*      particleToSplit;

	/* Add new particle to end of swarm */
	newParticle     = (GlobalParticle*) Swarm_CreateNewParticle( swarm, &newParticle_I );
	
	/* Copy particle information */
	particleToSplit = (GlobalParticle*) Swarm_ParticleInCellAt( swarm, lCell_I, cParticle_I );
	memcpy( newParticle, particleToSplit, swarm->particleExtensionMgr->finalSize );
	Swarm_AddParticleToCell( swarm, lCell_I, newParticle_I );

	/* Copy new position */
	memcpy( newParticle->coord, newCoord, sizeof(Coord) );
}
