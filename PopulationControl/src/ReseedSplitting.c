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
** $Id: ReseedSplitting.c 518 2007-10-11 08:07:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include "types.h"
#include "SplittingRoutine.h"
#include "ReseedSplitting.h"

#include <assert.h>
#include <string.h>
#include <math.h>

const Type ReseedSplitting_Type = "ReseedSplitting";

/*-------------------------------------------------------------------------------------------------------------------------
** Constructors
*/
ReseedSplitting* _ReseedSplitting_New(
		SizeT                                            _sizeOfSelf, 
		Type                                             type,
		Stg_Class_DeleteFunction*                        _delete,
		Stg_Class_PrintFunction*                         _print,
		Stg_Class_CopyFunction*                          _copy, 
		Stg_Component_DefaultConstructorFunction*        _defaultConstructor,
		Stg_Component_ConstructFunction*                 _construct,
		Stg_Component_BuildFunction*                     _build,
		Stg_Component_InitialiseFunction*                _initialise,
		Stg_Component_ExecuteFunction*                   _execute,
		Stg_Component_DestroyFunction*                   _destroy,		
		SplittingRoutine_SplitParticlesInCellFunction*   _splitInCell,
		Name                                             name )
{
	ReseedSplitting* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(ReseedSplitting) );
	self = (ReseedSplitting*)_SplittingRoutine_New( 
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
			_splitInCell,
			name );

	
	/* General info */

	/* Virtual Info */
	
	return self;
}

void _ReseedSplitting_Init( void* reseedSplitting, Index baseNum ) {
	ReseedSplitting* self = (ReseedSplitting*)reseedSplitting;
	
	self->isConstructed   = True;

	self->baseNum = baseNum;
	self->regionCount = baseNum * baseNum;
	if (self->dim == 3)
		self->regionCount *= baseNum;
	self->regionContainsParticleTbl = Memory_Alloc_Array( Bool, self->regionCount, "regionContainsParticleTbl" );
	memset( self->regionContainsParticleTbl, 0, self->regionCount * sizeof(Bool) );
}

void _ReseedSplitting_Delete( void* reseedSplitting ) {
	ReseedSplitting* self = (ReseedSplitting*)reseedSplitting;

	Memory_Free( self->regionContainsParticleTbl );

	_SplittingRoutine_Delete( self );
}

/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/
void* _ReseedSplitting_DefaultNew( Name name ) {
	return (void*) _ReseedSplitting_New(
			sizeof(ReseedSplitting),
			ReseedSplitting_Type,
			_ReseedSplitting_Delete,
			_SplittingRoutine_Print,
			_SplittingRoutine_Copy,
			_ReseedSplitting_DefaultNew,
			_ReseedSplitting_Construct,
			_SplittingRoutine_Build,
			_SplittingRoutine_Initialise,
			_SplittingRoutine_Execute,
			_SplittingRoutine_Destroy,
			_ReseedSplitting_SplitParticlesInCell,
			name );
}


void _ReseedSplitting_Construct( void* reseedSplitting, Stg_ComponentFactory* cf, void* data ) {
	ReseedSplitting*	     self            = (ReseedSplitting*) reseedSplitting;
	Index                    baseNum;

	_SplittingRoutine_Construct( self, cf, data );

	baseNum = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "baseNum", 
			(Index) pow( (double) self->idealParticleCount, 1.0/(double)(self->dim) ) );

	_ReseedSplitting_Init( self, baseNum );
}

/*-------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
void _ReseedSplitting_SplitParticlesInCell( void* reseedSplitting, void* _swarm, Cell_LocalIndex lCell_I ) {
	ReseedSplitting*        self              = (ReseedSplitting*) reseedSplitting;
	Swarm*                  swarm             = (Swarm*)  _swarm;
	Particle_InCellIndex    cellParticleCount = swarm->cellParticleCountTbl[ lCell_I ];
	Particle_InCellIndex    cParticle_I;
	GlobalParticle*         particle;
	Coord                   newCoord;
	Coord                   xi;
	IJK                     regionIJK         = { 0, 0, 0 };
	IJK                     regionSizes       = { 1, 1, 1 };
	Index                   region_I;
	Dimension_Index         dim_I;
	Dimension_Index         dim               = self->dim;
	double                  regionLength;
	FeMesh*     mesh              = (FeMesh*)((ElementCellLayout*)swarm->cellLayout)->mesh;

	/* Initialise all Bools to false */
	memset( self->regionContainsParticleTbl, 0, self->regionCount * sizeof(Bool) );
	regionSizes[ I_AXIS ] = regionSizes[ J_AXIS ] = self->baseNum;
	if ( dim == 3 )
		regionSizes[ K_AXIS ] = self->baseNum;
	regionLength = 2.0 / (double) self->baseNum;
		
	/* Fill table of regions with true if there is a particle in the region */
	for ( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {
		particle = (GlobalParticle*) Swarm_ParticleInCellAt( swarm, lCell_I, cParticle_I ) ;
		
		/* Calculate local coordinates */
		ElementType_ConvertGlobalCoordToElLocal(
				FeMesh_GetElementType( mesh, lCell_I ),
				mesh, 
				lCell_I, 
				particle->coord,
				xi );

		for ( dim_I = 0 ; dim_I < dim ; dim_I++ ) 
			regionIJK[ dim_I ] = (Index) ( (xi[ dim_I ] + 1.0) / regionLength );

		Dimension_3DTo1D( regionIJK, regionSizes, &region_I );

		self->regionContainsParticleTbl[ region_I ] = True;
	}

	/* Go through and add particles to empty regions */
	for ( region_I = 0 ; region_I < self->regionCount ; region_I++ ) {

		/* Don't do anything if region has particle */
		if ( self->regionContainsParticleTbl[ region_I ] )
			continue;

		/* Find local coordinate of empty region */
		Dimension_1DTo3D( region_I, regionSizes, regionIJK );
		for ( dim_I = 0 ; dim_I < dim ; dim_I++ ) {
			xi[ dim_I ] = (double) regionIJK[ dim_I ] * regionLength - 1.0 + 0.5 * regionLength;
		}

		/* Convert Local Coordinate to Global Coordinate */
		FeMesh_CoordLocalToGlobal( mesh, lCell_I, xi, newCoord );
			
		/* Work out particle to split by finding closest particle in this cell */
		cParticle_I = Swarm_FindClosestParticleInCell( swarm, lCell_I, dim, newCoord, NULL );

		/* Split particle */
		SplittingRoutine_SplitParticle( self, swarm, lCell_I, cParticle_I, newCoord );
	}
}
