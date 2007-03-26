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
** $Id: CellularAutomataVoronoi.c 376 2006-10-18 06:58:41Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>

#include "types.h"
#include "DiscreteVoronoi.h"
#include "CellularAutomataVoronoi.h"

#include <assert.h>
#include <string.h>
#include <math.h>

const unsigned int CELLULAR_AUTOMATA_UNCLAIMED = ((unsigned int) -1);
#define NO_CHECK                     ((unsigned int) -2)
#define VICTOR                       0
#define LOSER                        1

#define MEMORY_DELTA                 5

/* Textual name of this class */
const Type CellularAutomataVoronoi_Type = "CellularAutomataVoronoi";

/*-------------------------------------------------------------------------------------------------------------------------
** Constructors
*/
CellularAutomataVoronoi* _CellularAutomataVoronoi_New(
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
	DiscreteVoronoi_CalculateForCellFunction*          _calculate,
	DiscreteVoronoi_GetParticleIndexFunction*          _getParticleIndex,
	DiscreteVoronoi_GetVolumeFunction*                 _getVolume,
	DiscreteVoronoi_GetCentroidFunction*               _getCentroid,
	Name                                               name )
{
	CellularAutomataVoronoi* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(CellularAutomataVoronoi) );
	self = (CellularAutomataVoronoi*)_DiscreteVoronoi_New( 
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
		_calculate,
		_getParticleIndex,
		_getVolume,
		_getCentroid,
		name );

	
	/* General info */

	/* Virtual Info */
	
	return self;
}

void _CellularAutomataVoronoi_Init( void* cellularAutomataVoronoi, IJK resolution, Bool diagonalNeighbours ) {
	CellularAutomataVoronoi* self = (CellularAutomataVoronoi*)cellularAutomataVoronoi;

	self->diagonalNeighbours = diagonalNeighbours;

	if (self->dim == 2)
		resolution[ K_AXIS ] = 1;

	memcpy( self->resolution, resolution, sizeof(IJK) );

	/* Calculate Total number of cells */
	self->claimedCellCount = resolution[ I_AXIS ] * resolution[ J_AXIS ] * resolution[ K_AXIS ];

	/* Allocate for the volumes of each cell. */
	self->cellVolumes = Memory_Alloc_Array( double, self->claimedCellCount, "CellularAutomataVoronoi::cellVolumes" );
}


/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _CellularAutomataVoronoi_Delete( void* cellularAutomataVoronoi ) {
	CellularAutomataVoronoi*     self = (CellularAutomataVoronoi*)cellularAutomataVoronoi;
	Voronoi_CellIndex            vCell_I;
	CellularAutomataVoronoiCell* cell;
	
	/* Delete Neighbours and Battle History */
	for ( vCell_I = 0 ; vCell_I < self->claimedCellCount ; vCell_I++ ) {
		cell = &self->cellList[ vCell_I ];
		Memory_Free( cell->neighbourCellList );
		Memory_Free( cell->battleHistory.battlePair );
	}
	
	Memory_Free( self->cellList );
	Memory_Free( self->cellsToGrow.array );
	Memory_Free( self->cellsToCheck.array );
	FreeArray( self->cellVolumes );
	
	/* Delete parent */
	_DiscreteVoronoi_Delete( self );
}


void _CellularAutomataVoronoi_Print( void* cellularAutomataVoronoi, Stream* stream ) {
	CellularAutomataVoronoi* self = (CellularAutomataVoronoi*)cellularAutomataVoronoi;
	
	/* Print parent */
	_DiscreteVoronoi_Print( self, stream );
}

void* _CellularAutomataVoronoi_Copy( void* cellularAutomataVoronoi, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	CellularAutomataVoronoi*	self = (CellularAutomataVoronoi*)cellularAutomataVoronoi;
	CellularAutomataVoronoi*	newCellularAutomataVoronoi;
	
	newCellularAutomataVoronoi = (CellularAutomataVoronoi*)_DiscreteVoronoi_Copy( self, dest, deep, nameExt, ptrMap );
	
	return (void*)newCellularAutomataVoronoi;
}

void* _CellularAutomataVoronoi_DefaultNew( Name name ) {
	return (void*) _CellularAutomataVoronoi_New(
		sizeof(CellularAutomataVoronoi),
		CellularAutomataVoronoi_Type,
		_CellularAutomataVoronoi_Delete,
		_CellularAutomataVoronoi_Print,
		_CellularAutomataVoronoi_Copy,
		_CellularAutomataVoronoi_DefaultNew,
		_CellularAutomataVoronoi_Construct,
		_CellularAutomataVoronoi_Build,
		_CellularAutomataVoronoi_Initialise,
		_CellularAutomataVoronoi_Execute,
		_CellularAutomataVoronoi_Destroy,
		_CellularAutomataVoronoi_CalculateForCell,
		_CellularAutomataVoronoi_GetParticleIndex,
		_CellularAutomataVoronoi_GetVolume,
		_CellularAutomataVoronoi_GetCentroid,
		name );
}


void _CellularAutomataVoronoi_Construct( void* cellularAutomataVoronoi, Stg_ComponentFactory* cf, void* data ) {
	CellularAutomataVoronoi*	     self          = (CellularAutomataVoronoi*) cellularAutomataVoronoi;
	unsigned int                     defaultResolution;
	IJK                              resolution;

	_DiscreteVoronoi_Construct( self, cf, data );

	defaultResolution = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolution", 10 );
	resolution[ I_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionX", defaultResolution );
	resolution[ J_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionY", defaultResolution );
	resolution[ K_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionZ", defaultResolution );
	
	_CellularAutomataVoronoi_Init( 
		self, 
		resolution, 
		Stg_ComponentFactory_GetBool( cf, self->name, "diagonalNeighbours", True ) );
}

void _CellularAutomataVoronoi_Build( void* cellularAutomataVoronoi, void* data ) {
	CellularAutomataVoronoi*	self = (CellularAutomataVoronoi*)cellularAutomataVoronoi;

	_DiscreteVoronoi_Build( self, data );

	/* Allocate Memory for cells */
	self->cellList = Memory_Alloc_Array( CellularAutomataVoronoiCell, self->claimedCellCount, "cellList" );
	memset( self->cellList, 0, sizeof( CellularAutomataVoronoiCell ) * self->claimedCellCount );

	CellularAutomataVoronoiCellList_Build( &self->cellsToGrow );
	CellularAutomataVoronoiCellList_Build( &self->cellsToCheck );
}

void _CellularAutomataVoronoi_Initialise( void* cellularAutomataVoronoi, void* data ) {
	CellularAutomataVoronoi*     self = (CellularAutomataVoronoi*)cellularAutomataVoronoi;
	IJK                          vIJK;
	Voronoi_CellIndex            vCell_I;
	CellularAutomataVoronoiCell* cell;

	_DiscreteVoronoi_Initialise( self, data );

	for ( vIJK[ K_AXIS ] = 0 ; vIJK[ K_AXIS ] < self->resolution[ K_AXIS ] ; vIJK[ K_AXIS ]++ ) {
		for ( vIJK[ J_AXIS ] = 0 ; vIJK[ J_AXIS ] < self->resolution[ J_AXIS ] ; vIJK[ J_AXIS ]++ ) {
			for ( vIJK[ I_AXIS ] = 0 ; vIJK[ I_AXIS ] < self->resolution[ I_AXIS ] ; vIJK[ I_AXIS ]++ ) {
				Dimension_3DTo1D( vIJK, self->resolution, &vCell_I );
				cell = &self->cellList[ vCell_I ];

				/* Set up battle history */
				cell->battleHistory.battlesAllocated = 3;
				cell->battleHistory.battlePair = 
					Memory_Alloc_Array( BattlePair, cell->battleHistory.battlesAllocated, "battlePairArray" );

				/* Set up the neighbours - 
				 * We don't have to worry about cells being on side walls of the element here because that is 
				 * taken care of in CellularAutomataVoronoi_AddNeighbourToCell */
				if ( self->diagonalNeighbours ) {
					int add[3];
					/* Use loop to get all 26 neighbours (8 in 2D) */
					for ( add[ I_AXIS ] = -1 ; add[ I_AXIS ] <= 1 ; add[ I_AXIS ]++ ) {
						for ( add[ J_AXIS ] = -1 ; add[ J_AXIS ] <= 1 ; add[ J_AXIS ]++ ) {
							for ( add[ K_AXIS ] = -1 ; add[ K_AXIS ] <= 1 ; add[ K_AXIS ]++ ) {
								
								/* Make sure that this isn't the same cell as itself */
								if ( add[ I_AXIS ] == 0 && add[ J_AXIS ] == 0 && add[ K_AXIS ] == 0 )
									continue;

								CellularAutomataVoronoi_AddNeighbourToCell( self, cell, 
													    vIJK[I_AXIS]+add[I_AXIS], vIJK[J_AXIS]+add[J_AXIS], vIJK[K_AXIS]+add[K_AXIS]);
							}
						}
					}
				}
				else {
					/* Not using neighbours as diagonals - There are 6 in 3D and 4 in 2D */
					CellularAutomataVoronoi_AddNeighbourToCell( self, cell, vIJK[I_AXIS] - 1, vIJK[J_AXIS], vIJK[K_AXIS]);
					CellularAutomataVoronoi_AddNeighbourToCell( self, cell, vIJK[I_AXIS] + 1, vIJK[J_AXIS], vIJK[K_AXIS]);
						
					CellularAutomataVoronoi_AddNeighbourToCell( self, cell, vIJK[I_AXIS], vIJK[J_AXIS] - 1, vIJK[K_AXIS]);
					CellularAutomataVoronoi_AddNeighbourToCell( self, cell, vIJK[I_AXIS], vIJK[J_AXIS] + 1, vIJK[K_AXIS]);
						
					CellularAutomataVoronoi_AddNeighbourToCell( self, cell, vIJK[I_AXIS], vIJK[J_AXIS], vIJK[K_AXIS] - 1);
					CellularAutomataVoronoi_AddNeighbourToCell( self, cell, vIJK[I_AXIS], vIJK[J_AXIS], vIJK[K_AXIS] + 1);
				}
			}
		}
	}
}
void _CellularAutomataVoronoi_Execute( void* cellularAutomataVoronoi, void* data ) {
	CellularAutomataVoronoi*	self = (CellularAutomataVoronoi*)cellularAutomataVoronoi;
	
	_DiscreteVoronoi_Execute( self, data );
}
void _CellularAutomataVoronoi_Destroy( void* cellularAutomataVoronoi, void* data ) {
	CellularAutomataVoronoi*	self = (CellularAutomataVoronoi*)cellularAutomataVoronoi;
	
	_DiscreteVoronoi_Destroy( self, data );
}

void _CellularAutomataVoronoi_CalculateForCell( void* cellularAutomataVoronoi, void* swarm, Cell_LocalIndex lCell_I ) {
	CellularAutomataVoronoi*             self            = (CellularAutomataVoronoi*)  cellularAutomataVoronoi;

	/* Initialise all values on cells */
	CellularAutomataVoronoi_InitialiseCell( self, swarm, lCell_I );

	/* Seed Cellular Automata */
	CellularAutomataVoronoi_Seed( self );
	
	/* Do CellularAutomataVoronoi - Keep claiming cells until there are no more being claimed */
	do {
		CellularAutomataVoronoi_GrowCells( self );
		CellularAutomataVoronoi_CheckCells( self );
	} while ( self->cellsToGrow.cellCount != 0 );
}

Particle_InCellIndex _CellularAutomataVoronoi_GetParticleIndex( void* discreteVoronoi, Voronoi_CellIndex vCell_I ){
	CellularAutomataVoronoi*             self            = (CellularAutomataVoronoi*)  discreteVoronoi;

	return self->cellList[ vCell_I ].claimedParticle;
}

double _CellularAutomataVoronoi_GetVolume(void* discreteVoronoi, Voronoi_CellIndex vCell_I ){
	CellularAutomataVoronoi*             self            = (CellularAutomataVoronoi*)  discreteVoronoi;

	return self->cellVolumes[vCell_I];
}

void _CellularAutomataVoronoi_GetCentroid(void* discreteVoronoi, Voronoi_CellIndex vCell_I, Coord centroid ){
	CellularAutomataVoronoi*             self            = (CellularAutomataVoronoi*)  discreteVoronoi;

	memcpy( centroid, self->cellList[ vCell_I ].centroid, self->dim * sizeof(double) );
}

/* Other functions */
void CellularAutomataVoronoi_CheckCells( void* cellularAutomataVoronoi ) {
	CellularAutomataVoronoi*             self            = (CellularAutomataVoronoi*)  cellularAutomataVoronoi;
	Voronoi_CellIndex                    checkCellCount;
	Voronoi_CellIndex                    checkCell_I;
	CellularAutomataVoronoiCell*         cell;
	Particle_InCellIndex                 claimedParticle;
	Particle_InCellIndex                 particleToCheck;
	
	checkCellCount = self->cellsToCheck.cellCount;

	for ( checkCell_I = 0 ; checkCell_I < checkCellCount ; checkCell_I++ ) {
		cell            = self->cellsToCheck.array[ checkCell_I ];
		claimedParticle = cell->claimedParticle;
		particleToCheck = cell->particleToCheck;

		/* If this cell is unclaimed - then claim it for this particle */
		if ( claimedParticle == CELLULAR_AUTOMATA_UNCLAIMED ) {
			cell->claimedParticle = particleToCheck;
			CellularAutomataVoronoiCellList_AddCell( &self->cellsToGrow, cell );
		}
		else {
			/* This cell has already been claimed - therefore there needs to be battle between these particles
			 * to see who gets to own it */
			if ( particleToCheck == CellularAutomataVoronoi_Battle( self, cell, claimedParticle, particleToCheck ) ) {
				cell->claimedParticle = particleToCheck;
				CellularAutomataVoronoiCellList_AddCell( &self->cellsToGrow, cell );
			}
		}
		cell->particleToCheck = NO_CHECK;
	}

	/* Now that we've looked through all the cells to check if they should be claimed - 
	 * we can set the number to be checked back to zero */
	self->cellsToCheck.cellCount = 0;
}

void CellularAutomataVoronoi_GrowCells( void* cellularAutomataVoronoi ) {
	CellularAutomataVoronoi*             self            = (CellularAutomataVoronoi*)  cellularAutomataVoronoi;
	Voronoi_CellIndex                    growCellCount;
	Voronoi_CellIndex                    growCell_I;
	Voronoi_CellIndex                    neighbourCount;
	Voronoi_CellIndex                    neighbour_I;
	CellularAutomataVoronoiCell*         cell;
	CellularAutomataVoronoiCell*         neighbour;
	Particle_InCellIndex                 claimedParticle;
	
	growCellCount = self->cellsToGrow.cellCount;

	for ( growCell_I = 0 ; growCell_I < growCellCount ; growCell_I++ ) {
		cell            = self->cellsToGrow.array[ growCell_I ];
		neighbourCount  = cell->neighbourCount;
		claimedParticle = cell->claimedParticle;
	
		for ( neighbour_I = 0 ; neighbour_I < neighbourCount ; neighbour_I++ ) {
			neighbour = cell->neighbourCellList[ neighbour_I ];

			/* If this particle has already claimed this cell - then ignore it */
			if ( neighbour->claimedParticle == claimedParticle )
				continue;

			/* If this cell is already going to try this particle - then we don't need to do it twice */
			if ( neighbour->particleToCheck == claimedParticle )
				continue;
			
			/* If this cell is going to try a different particle - then do battle */
			if ( neighbour->particleToCheck != NO_CHECK ) {
				neighbour->particleToCheck = 
					CellularAutomataVoronoi_Battle( self, neighbour, neighbour->particleToCheck, claimedParticle );
				continue;
			}

			/* Otherwise - Add this cell as one of the ones to be checked */
			neighbour->particleToCheck = claimedParticle;
			CellularAutomataVoronoiCellList_AddCell( &self->cellsToCheck, neighbour );
		}
	}

	/* Now that we've looked through all the cells to grow 
	 * we can set the number to grow back to zero */
	self->cellsToGrow.cellCount = 0;
}


Particle_InCellIndex CellularAutomataVoronoi_Battle( 
	void*                                              cellularAutomataVoronoi, 
	CellularAutomataVoronoiCell*                       cell, 
	Particle_InCellIndex                               champion_I, 
	Particle_InCellIndex                               contender_I )
{
	CellularAutomataVoronoi*  self            = (CellularAutomataVoronoi*)  cellularAutomataVoronoi;
	double*                   cellCentroid;
	void*                     champion;
	void*                     contender;
	Coord                     championCoord;
	Coord                     contenderCoord;
	double                    distanceToChampion;
	double                    distanceToContender;
	Index                     battle_I;
	Index                     battleCount = cell->battleHistory.battleCount;
	Swarm*                    swarm;
	Cell_LocalIndex           lCell_I;
	Particle_InCellIndex*     battlePair;

	FiniteElement_Mesh*       mesh;

	/* Check if battle has already been fought */
	for ( battle_I = 0 ; battle_I < battleCount ; battle_I++ ) {
		battlePair = cell->battleHistory.battlePair[battle_I];
		if (       ( champion_I == battlePair[LOSER ] && contender_I == battlePair[VICTOR] )
			   || ( champion_I == battlePair[VICTOR] && contender_I == battlePair[LOSER ] ) )
			return battlePair[VICTOR];
	}
	
	swarm               = self->swarm;
	lCell_I             = self->lCell_I;
	cellCentroid        = cell->centroid;
	champion            = Swarm_ParticleInCellAt( swarm, lCell_I, champion_I );
	contender           = Swarm_ParticleInCellAt( swarm, lCell_I, contender_I );
	mesh                = (FiniteElement_Mesh*)(((ElementCellLayout*)swarm->cellLayout)->mesh); /* Assume ElementCellLayout */
	
	if ( swarm->particleLayout->coordSystem == GlobalCoordSystem ) {
		memcpy( championCoord, ((GlobalParticle*)champion)->coord, sizeof(Coord) );
		memcpy( contenderCoord, ((GlobalParticle*)contender)->coord, sizeof(Coord) );
	}
	else {
		/* LocalCoordSystem need to convert to global */
		FiniteElement_Mesh_CalcGlobalCoordFromLocalCoord(
			mesh,
			swarm->dim,
			lCell_I,
			((LocalParticle*)champion)->xi,
			championCoord );
		FiniteElement_Mesh_CalcGlobalCoordFromLocalCoord(
			mesh,
			swarm->dim,
			lCell_I,
			((LocalParticle*)contender)->xi,
			contenderCoord );
	}
						


	/* Calculate 'distances' from cell to each particle in the battle 
	 * ('distance' isn't actually the true distance - but the distance squared */
	distanceToChampion = 
		(championCoord[I_AXIS] - cellCentroid[I_AXIS])*(championCoord[I_AXIS] - cellCentroid[I_AXIS]) +
		(championCoord[J_AXIS] - cellCentroid[J_AXIS])*(championCoord[J_AXIS] - cellCentroid[J_AXIS]) ;
	
	distanceToContender = 
		(contenderCoord[I_AXIS] - cellCentroid[I_AXIS])*(contenderCoord[I_AXIS] - cellCentroid[I_AXIS]) +
		(contenderCoord[J_AXIS] - cellCentroid[J_AXIS])*(contenderCoord[J_AXIS] - cellCentroid[J_AXIS]) ;

	if ( self->dim == 3 ) {
		distanceToChampion += 
			(championCoord[K_AXIS] - cellCentroid[K_AXIS])*(championCoord[K_AXIS] - cellCentroid[K_AXIS]);
		distanceToContender += 
			(contenderCoord[K_AXIS] - cellCentroid[K_AXIS])*(contenderCoord[K_AXIS] - cellCentroid[K_AXIS]);
	}

	/* Check if we need to add memory to battle history */
	if ( cell->battleHistory.battleCount >= cell->battleHistory.battlesAllocated ) {
		cell->battleHistory.battlesAllocated = cell->battleHistory.battleCount + 2;
		cell->battleHistory.battlePair = 
			Memory_Realloc_Array( cell->battleHistory.battlePair, BattlePair, cell->battleHistory.battlesAllocated );
	}
	battlePair = cell->battleHistory.battlePair[ cell->battleHistory.battleCount ];
	cell->battleHistory.battleCount++;

	/* Record victory */
	if ( distanceToContender < distanceToChampion ) {
		battlePair[ VICTOR ] = contender_I;
		battlePair[ LOSER  ] = champion_I;
	}
	else {
		battlePair[ VICTOR ] = champion_I;
		battlePair[ LOSER  ] = contender_I;
	}

	return battlePair[ VICTOR ];
}

void CellularAutomataVoronoi_InitialiseCell( void* cellularAutomataVoronoi, Swarm* swarm, Cell_LocalIndex lCell_I ) {
	CellularAutomataVoronoi*             self            = (CellularAutomataVoronoi*)  cellularAutomataVoronoi;
	CellularAutomataVoronoiCell*         cell;
	Voronoi_CellIndex                    vCell_I;
	Dimension_Index                      dim_I;
	Dimension_Index                      dim = self->dim;
	Index*                               resolution      = self->resolution;
	IJK                                  ijk_I;

	/* Find spacing of sub-cells in the parent domain. */
	for ( dim_I = 0 ; dim_I < dim ; dim_I++ )
		self->dx[ dim_I ] = 2.0 / ((double) resolution[ dim_I ]);

	/* Update 'Current' info */
	self->swarm   = swarm;
	self->lCell_I = lCell_I;

	/* Reset Lists */
	self->cellsToGrow.cellCount  = 0;
	self->cellsToCheck.cellCount = 0;

	/* Reset Cells */
	for ( vCell_I = 0 ; vCell_I < self->claimedCellCount ; vCell_I++ ) {
		Dimension_1DTo3D( vCell_I, resolution, ijk_I );
		cell = &self->cellList[ vCell_I ];

		cell->claimedParticle = CELLULAR_AUTOMATA_UNCLAIMED;
		cell->particleToCheck = NO_CHECK;
		cell->battleHistory.battleCount = 0;
	}

	/* Calculate the volume and centroid of each sub-cell. */
	CellularAutomataVoronoi_CalcSubCells( self, swarm, lCell_I );
}

void CellularAutomataVoronoi_AddNeighbourToCell( void* cellularAutomataVoronoi, CellularAutomataVoronoiCell* cell, Index iIndex, Index jIndex, Index kIndex ) {
	CellularAutomataVoronoi*             self            = (CellularAutomataVoronoi*)  cellularAutomataVoronoi;
	Index*                               resolution      = self->resolution;
	Voronoi_CellIndex                    neighbourIndex;

	if ( iIndex < resolution[ I_AXIS ] && jIndex < resolution[ J_AXIS ] && kIndex < resolution[ K_AXIS ] ) {
		/* See if memory needs to be allocated */
		if ( cell->neighbourCount >= cell->neighboursAllocated ) {
			cell->neighboursAllocated = cell->neighbourCount + 1;
			cell->neighbourCellList = 
				Memory_Realloc_Array( cell->neighbourCellList, CellularAutomataVoronoiCell*, cell->neighboursAllocated );
		}

		/* Add neigbour to cell's list */
		Dimension_3DTo1D_3( iIndex, jIndex, kIndex, resolution[0], resolution[1], resolution[2], &neighbourIndex );
		cell->neighbourCellList[ cell->neighbourCount ] = &self->cellList[ neighbourIndex ];
		cell->neighbourCount++;
	}
}

void CellularAutomataVoronoi_Seed( void* cellularAutomataVoronoi ) {
	CellularAutomataVoronoi*             self              = (CellularAutomataVoronoi*)  cellularAutomataVoronoi;
	Swarm*                               swarm             = self->swarm;
	Cell_LocalIndex                      lCell_I           = self->lCell_I;
	void*                                particle;
	Particle_InCellIndex                 cParticle_I;
	Particle_InCellIndex                 cellParticleCount = swarm->cellParticleCountTbl[lCell_I];
	IJK                                  vCell_IJK         = {0,0,0};
	Voronoi_CellIndex                    vCell_I;
	CellularAutomataVoronoiCell*         cell;
	Dimension_Index                      dim               = self->dim;
	Dimension_Index                      dim_I;
	Index*                               resolution        = self->resolution;
	FiniteElement_Mesh*                  mesh;
	Coord                                localCoord;

	mesh                = (FiniteElement_Mesh*)(((ElementCellLayout*)swarm->cellLayout)->mesh); /* Assume ElementCellLayout */

	/* Loop over all the particles in the cell - assigning it to the a voronoi cell it is in */
	for ( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {
		particle = Swarm_ParticleInCellAt( swarm, lCell_I, cParticle_I );

		if ( swarm->particleLayout->coordSystem == GlobalCoordSystem ) {
			/* Must convert global to local. */
			FiniteElement_Mesh_CalcLocalCoordFromGlobalCoord( mesh, 
									  lCell_I, 
									  ((GlobalParticle*)particle)->coord, 
									  localCoord );
		}
		else {
			/* Now we need to coordinate in locals. */
			memcpy( localCoord, ((LocalParticle*)particle)->xi, sizeof(Coord) );
		}

		/* Find which discrete voronoi cell this particle is in */
		for ( dim_I = 0 ; dim_I < dim ; dim_I++ ) {
			vCell_IJK[ dim_I ] = (Index)((localCoord[dim_I] + 1.0) / self->dx[dim_I]);
			
			/* Check if particle is right on edge of element */
			if ( vCell_IJK[ dim_I ] == resolution[ dim_I ] )
				vCell_IJK[ dim_I ]--;
		}
		Dimension_3DTo1D( vCell_IJK, resolution, &vCell_I );

		/* Get pointer to this cell */
		cell = &self->cellList[ vCell_I ];

		/* Check if there is another particle in this cell - Then do growing and checking stage before anything else */
		if ( cell->claimedParticle != CELLULAR_AUTOMATA_UNCLAIMED ) {
			Particle_InCellIndex previousClaimedParticle = cell->claimedParticle;
			CellularAutomataVoronoi_GrowCells( self );
			CellularAutomataVoronoi_CheckCells( self );
			cell->claimedParticle = cParticle_I;
			CellularAutomataVoronoiCellList_AddCell( &self->cellsToGrow, cell );
			CellularAutomataVoronoi_GrowCells( self );
			CellularAutomataVoronoi_CheckCells( self );

			cell->claimedParticle = CellularAutomataVoronoi_Battle( self, cell, previousClaimedParticle, cParticle_I );
		}
		else {
			cell->claimedParticle = cParticle_I;
			CellularAutomataVoronoiCellList_AddCell( &self->cellsToGrow, cell );
		}
	}
}

void CellularAutomataVoronoi_AsciiArt( void* cellularAutomataVoronoi ) {
	CellularAutomataVoronoi*             self              = (CellularAutomataVoronoi*)  cellularAutomataVoronoi;
	IJK                                  vCell_IJK         = { 0,0,0 };
	Stream*                              stream            = Journal_MyStream( Info_Type, self );
	Voronoi_CellIndex                    vCell_I;
	CellularAutomataVoronoiCell*         cell;
	char                                 c;

	assert( self->dim == 2 );
		
	/* Draw Line at top of box */
	Journal_Printf( stream, " " );
	for ( vCell_IJK[ I_AXIS ] = 0 ; vCell_IJK[ I_AXIS ] < self->resolution[ I_AXIS ] ; vCell_IJK[ I_AXIS ]++ ) {
		Journal_Printf( stream, "_" );
	}
	Journal_Printf( stream, "\n" );

	for ( vCell_IJK[J_AXIS] = self->resolution[J_AXIS] - 1 ; vCell_IJK[ J_AXIS ] != (Index) -1 ; vCell_IJK[ J_AXIS ]-- ) {
		Journal_Printf( stream, "|" );
		for ( vCell_IJK[ I_AXIS ] = 0 ; vCell_IJK[ I_AXIS ] < self->resolution[ I_AXIS ] ; vCell_IJK[ I_AXIS ]++ ) {
			Dimension_3DTo1D( vCell_IJK, self->resolution, &vCell_I );
			cell = &self->cellList[ vCell_I ];

			if ( cell->claimedParticle == CELLULAR_AUTOMATA_UNCLAIMED ) {
				if ( cell->particleToCheck == NO_CHECK ) {
					c = ' ';
				}
				else {
					c = '*';
				}
			}
			else { /* Cell is claimed */
				if ( cell->particleToCheck == NO_CHECK ) {
					c = '0' + cell->claimedParticle;
				}
				else {
					c = '?';
				}
			}


			Journal_Printf( stream, "%c", c );
		}
		Journal_Printf( stream, "|\n" );
	}
	/* Draw Line at bottom of box */
	Journal_Printf( stream, " " );
	for ( vCell_IJK[ I_AXIS ] = 0 ; vCell_IJK[ I_AXIS ] < self->resolution[ I_AXIS ] ; vCell_IJK[ I_AXIS ]++ )
		Journal_Printf( stream, "-" );
	Journal_Printf( stream, "\n" );
}


void CellularAutomataVoronoiCellList_Build( CellularAutomataVoronoiCellList* cellList ) {
	cellList->cellsAllocated = MEMORY_DELTA;
	cellList->array = Memory_Alloc_Array( CellularAutomataVoronoiCell*, cellList->cellsAllocated, "CellList" );
}

void CellularAutomataVoronoiCellList_AddCell( CellularAutomataVoronoiCellList* cellList, CellularAutomataVoronoiCell* cell ) {
	/* See if memory needs to be allocated */
	if ( cellList->cellCount >= cellList->cellsAllocated ) {
		cellList->cellsAllocated = cellList->cellCount + MEMORY_DELTA;
		cellList->array = Memory_Realloc_Array( cellList->array, CellularAutomataVoronoiCell*, cellList->cellsAllocated );
	}

	/* Add cell to this list */
	cellList->array[ cellList->cellCount ] = cell;
			
	/* Increment count for list */
	cellList->cellCount++;
}

void CellularAutomataVoronoi_CalcSubCells( CellularAutomataVoronoi* self, Swarm* swarm, unsigned cellInd ) {
	/* Sanity check. */
	assert( self );
	assert( swarm );

	if( self->dim == 2 )
		CellularAutomataVoronoi_CalcSubCells2D( self, swarm, cellInd );
	else
		CellularAutomataVoronoi_CalcSubCells3D( self, swarm, cellInd );
}

void CellularAutomataVoronoi_CalcSubCells2D( CellularAutomataVoronoi* self, Swarm* swarm, unsigned cellInd ) {
	FiniteElement_Mesh*	mesh;
	unsigned		nSubCells;
	unsigned*		res;
	double**		gCrds;
	double**		lCrds;
	unsigned		d_i, d_j;

	/* Sanity check. */
	assert( self );
	assert( swarm );

	/* Shortcuts. */
	nSubCells = self->claimedCellCount;
	res = self->resolution;
	mesh = (FiniteElement_Mesh*)(((ElementCellLayout*)swarm->cellLayout)->mesh);

	/* NOTE: It is assumed in an earlier function that the cell layout is based on the
	   mesh's elements; I continue that assumption here. */

	/* Need space for the coordinates that comprise this sub-cell. */
	gCrds = Memory_Alloc_2DArray( double, 4, 3, "" );
	lCrds = Memory_Alloc_2DArray( double, 4, 3, "" );

	/* Calculate the volume of each cell. */
	for( d_j = 0; d_j < res[1]; d_j++ ) {
		lCrds[0][1] = (double)d_j * self->dx[1] - 1.0;
		lCrds[1][1] = (double)d_j * self->dx[1] - 1.0;
		lCrds[2][1] = (double)(d_j + 1) * self->dx[1] - 1.0;
		lCrds[3][1] = (double)(d_j + 1) * self->dx[1] - 1.0;

		for( d_i = 0; d_i < res[0]; d_i++ ) {
			unsigned	subCellInd = d_j * res[0] + d_i;
			unsigned	c_i;

			lCrds[0][0] = (double)d_i * self->dx[0] - 1.0;
			lCrds[1][0] = (double)(d_i + 1) * self->dx[0] - 1.0;
			lCrds[2][0] = (double)d_i * self->dx[0] - 1.0;
			lCrds[3][0] = (double)(d_i + 1) * self->dx[0] - 1.0;

			/* Map the local coordinates back to globals. */
			for( c_i = 0; c_i < 4; c_i++ ) {
				FiniteElement_Mesh_CalcGlobalCoordFromLocalCoord( mesh, 
										  self->dim, 
										  cellInd, 
										  lCrds[c_i], 
										  gCrds[c_i] );
			}

			/* Calculate the volume. */
			self->cellVolumes[subCellInd] = CellularAutomataVoronoi_QuadArea( self, gCrds );

			/* Calculate centroid. */
			CellularAutomataVoronoi_QuadCentroid( self, gCrds, self->cellList[subCellInd].centroid );
		}
	}

	/* Free memory. */
	FreeArray( gCrds );
	FreeArray( lCrds );
}

void CellularAutomataVoronoi_CalcSubCells3D( CellularAutomataVoronoi* self, Swarm* swarm, unsigned cellInd ) {
	const double		volWeight = 8.0;
	const double		sign[3][8] = {{-1, 1, -1, 1, -1, 1, -1, 1}, 
					      {-1, -1, 1, 1, -1, -1, 1, 1}, 
					      {-1, -1, -1, -1, 1, 1, 1, 1}};
	FiniteElement_Mesh*	mesh;
	unsigned		nSubCells;
	unsigned*		res;
	double**		gCrds;
	double**		lCrds;
	Coord*			gCrdPtrs[8];
	ElementType*		elType;
	unsigned		d_i, d_j, d_k;

	/* Sanity check. */
	assert( self );
	assert( swarm );

	/* Shortcuts. */
	nSubCells = self->claimedCellCount;
	res = self->resolution;
	mesh = (FiniteElement_Mesh*)(((ElementCellLayout*)swarm->cellLayout)->mesh);
	elType = FiniteElement_Mesh_ElementTypeAt( mesh, cellInd );

	/* NOTE: It is assumed in an earlier function that the cell layout is based on the
	   mesh's elements; I continue that assumption here. */

	/* Need space for the coordinates that comprise this sub-cell. */
	gCrds = Memory_Alloc_2DArray( double, 8, 3, "" );
	lCrds = Memory_Alloc_2DArray( double, 8, 3, "" );

	/* Calculate the volume of each cell. */
	for( d_k = 0; d_k < res[2]; d_k++ ) {
		lCrds[0][2] = (double)d_k * self->dx[2] - 1.0;
		lCrds[1][2] = (double)d_k * self->dx[2] - 1.0;
		lCrds[2][2] = (double)d_k * self->dx[2] - 1.0;
		lCrds[3][2] = (double)d_k * self->dx[2] - 1.0;
		lCrds[4][2] = (double)(d_k + 1) * self->dx[2] - 1.0;
		lCrds[5][2] = (double)(d_k + 1) * self->dx[2] - 1.0;
		lCrds[6][2] = (double)(d_k + 1) * self->dx[2] - 1.0;
		lCrds[7][2] = (double)(d_k + 1) * self->dx[2] - 1.0;

		for( d_j = 0; d_j < res[1]; d_j++ ) {
			lCrds[0][1] = (double)d_j * self->dx[1] - 1.0;
			lCrds[1][1] = (double)d_j * self->dx[1] - 1.0;
			lCrds[2][1] = (double)(d_j + 1) * self->dx[1] - 1.0;
			lCrds[3][1] = (double)(d_j + 1) * self->dx[1] - 1.0;
			lCrds[4][1] = (double)d_j * self->dx[1] - 1.0;
			lCrds[5][1] = (double)d_j * self->dx[1] - 1.0;
			lCrds[6][1] = (double)(d_j + 1) * self->dx[1] - 1.0;
			lCrds[7][1] = (double)(d_j + 1) * self->dx[1] - 1.0;

			for( d_i = 0; d_i < res[0]; d_i++ ) {
				unsigned	subCellInd = d_k * res[0] * res[1] + d_j * res[0] + d_i;
				double		jac[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
				double		jacDet;
				unsigned	c_i;

				/* Calculate the local coordinates. */
				lCrds[0][0] = (double)d_i * self->dx[0] - 1.0;
				lCrds[1][0] = (double)(d_i + 1) * self->dx[0] - 1.0;
				lCrds[2][0] = (double)d_i * self->dx[0] - 1.0;
				lCrds[3][0] = (double)(d_i + 1) * self->dx[0] - 1.0;
				lCrds[4][0] = (double)d_i * self->dx[0] - 1.0;
				lCrds[5][0] = (double)(d_i + 1) * self->dx[0] - 1.0;
				lCrds[6][0] = (double)d_i * self->dx[0] - 1.0;
				lCrds[7][0] = (double)(d_i + 1) * self->dx[0] - 1.0;

				/* Map the local coordinates back to globals. */
				for( c_i = 0; c_i < 8; c_i++ ) {
					/* Build a list of global coordinate pointers. */
					gCrdPtrs[c_i] = (Coord*)(gCrds + c_i);

					/* Convert to local coordinates. */
					FiniteElement_Mesh_CalcGlobalCoordFromLocalCoord( mesh, 
											  self->dim, 
											  cellInd, 
											  lCrds[c_i], 
											  gCrds[c_i] );

					/* Calculate the jacobian. */
					jac[0][0] += 0.125 * sign[0][c_i] * gCrds[c_i][0];
					jac[1][0] += 0.125 * sign[1][c_i] * gCrds[c_i][0];
					jac[2][0] += 0.125 * sign[2][c_i] * gCrds[c_i][0];
					jac[0][1] += 0.125 * sign[0][c_i] * gCrds[c_i][1];
					jac[1][1] += 0.125 * sign[1][c_i] * gCrds[c_i][1];
					jac[2][1] += 0.125 * sign[2][c_i] * gCrds[c_i][1];
					jac[0][2] += 0.125 * sign[0][c_i] * gCrds[c_i][2];
					jac[1][2] += 0.125 * sign[1][c_i] * gCrds[c_i][2];
					jac[2][2] += 0.125 * sign[2][c_i] * gCrds[c_i][2];
				}

				/* Calculate centroid. */
				CellularAutomataVoronoi_HexCentroid( self, gCrds, self->cellList[subCellInd].centroid );

				/* Calculate the volume. Doing this for irregular hexahedron is a little tricky. */
				jacDet = jac[0][0] * (jac[1][1] * jac[2][2] - jac[2][1] * jac[1][2]) - 
					jac[1][0] * (jac[0][1] * jac[2][2] - jac[2][1] * jac[0][2]) + 
					jac[2][0] * (jac[0][1] * jac[1][2] - jac[1][1] * jac[0][2]);
				self->cellVolumes[subCellInd] = volWeight * jacDet;
			}
		}
	}

	/* Free memory. */
	FreeArray( gCrds );
	FreeArray( lCrds );
}

double CellularAutomataVoronoi_QuadArea( CellularAutomataVoronoi* self, double** gCrds ) {
	double	area;
	double	vecs[4][2];
	double	diags[2][2];

	/* Sanity check. */
	assert( self );
	assert( gCrds );

	vecs[0][0] = gCrds[1][0] - gCrds[0][0];	vecs[0][1] = gCrds[1][1] - gCrds[0][1];
	vecs[1][0] = gCrds[3][0] - gCrds[1][0];	vecs[1][1] = gCrds[3][1] - gCrds[1][1];
	vecs[2][0] = gCrds[2][0] - gCrds[3][0];	vecs[2][1] = gCrds[2][1] - gCrds[3][1];
	vecs[3][0] = gCrds[0][0] - gCrds[2][0];	vecs[3][1] = gCrds[0][1] - gCrds[2][1];

	diags[0][0] = vecs[0][0] + vecs[1][0];	diags[0][1] = vecs[0][1] + vecs[1][1];
	diags[1][0] = vecs[1][0] + vecs[2][0];	diags[1][1] = vecs[1][1] + vecs[2][1];

	/* Calculate the area from the diagonals. */
	area = diags[0][0] * diags[1][1] - diags[0][1] * diags[1][0];
	area *= (area < 0.0) ? -0.5 : 0.5;

	return area;
}

void CellularAutomataVoronoi_QuadCentroid( CellularAutomataVoronoi* self, double** gCrds, Coord centroid ) {
	unsigned	c_i;

	/* Sanity check. */
	assert( self );
	assert( gCrds );

	/* The centroid of a convex quadrilateral is the average of it's corners. */
	centroid[0] = gCrds[0][0];
	centroid[1] = gCrds[0][1];
	for( c_i = 1; c_i < 4; c_i++ ) {
		centroid[0] += gCrds[c_i][0];
		centroid[1] += gCrds[c_i][1];
	}
	centroid[0] *= 0.25;
	centroid[1] *= 0.25;
}

void CellularAutomataVoronoi_HexCentroid( CellularAutomataVoronoi* self, double** gCrds, Coord centroid ) {
	unsigned	c_i;

	/* Sanity check. */
	assert( self );
	assert( gCrds );

	/* The centroid of a convex hexahedron is the average of it's corners. */
	centroid[0] = gCrds[0][0];
	centroid[1] = gCrds[0][1];
	centroid[2] = gCrds[0][2];
	for( c_i = 1; c_i < 8; c_i++ ) {
		centroid[0] += gCrds[c_i][0];
		centroid[1] += gCrds[c_i][1];
		centroid[2] += gCrds[c_i][2];
	}
	centroid[0] *= 0.125;
	centroid[1] *= 0.125;
	centroid[2] *= 0.125;
}
