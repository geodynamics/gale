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
** $Id: DiscreteVoronoi.c 374 2006-10-12 08:59:41Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>

#include "types.h"
#include "DiscreteVoronoi.h"

#include <assert.h>
#include <string.h>
#include <math.h>

/* Textual name of this class */
const Type DiscreteVoronoi_Type = "DiscreteVoronoi";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/


DiscreteVoronoi* _DiscreteVoronoi_New(
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
		DiscreteVoronoi_CalculateForCellFunction*          _calculateForCell,
		DiscreteVoronoi_GetParticleIndexFunction*          _getParticleIndex,
		DiscreteVoronoi_GetVolumeFunction*                 _getVolume,
		DiscreteVoronoi_GetCentroidFunction*               _getCentroid,
		Name                                               name )
{
	DiscreteVoronoi* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(DiscreteVoronoi) );
	self = (DiscreteVoronoi*)_Stg_Component_New( 
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
	self->_calculateForCell = _calculateForCell;
	self->_getParticleIndex = _getParticleIndex;
	self->_getVolume        = _getVolume;
	self->_getCentroid      = _getCentroid;
	
	return self;
}

void _DiscreteVoronoi_Init( void* discreteVoronoi, Dimension_Index dim ) {
	DiscreteVoronoi* self = (DiscreteVoronoi*)discreteVoronoi;
	
	self->isConstructed = True;
	self->dim = dim;
}


void DiscreteVoronoi_InitAll( void* discreteVoronoi, Dimension_Index dim ){
	DiscreteVoronoi* self = (DiscreteVoronoi*)discreteVoronoi;

	_DiscreteVoronoi_Init( self, dim );
}
	

/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _DiscreteVoronoi_Delete( void* discreteVoronoi ) {
	DiscreteVoronoi* self = (DiscreteVoronoi*)discreteVoronoi;
	
	/* Delete parent */
	_Stg_Component_Delete( self );
}


void _DiscreteVoronoi_Print( void* discreteVoronoi, Stream* stream ) {
	DiscreteVoronoi* self = (DiscreteVoronoi*)discreteVoronoi;
	
	/* Print parent */
	_Stg_Component_Print( self, stream );
}


void* _DiscreteVoronoi_Copy( void* discreteVoronoi, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	DiscreteVoronoi*	self = (DiscreteVoronoi*)discreteVoronoi;
	DiscreteVoronoi*	newDiscreteVoronoi;
	
	newDiscreteVoronoi = (DiscreteVoronoi*)_Stg_Component_Copy( self, dest, deep, nameExt, ptrMap );
	
	return (void*)newDiscreteVoronoi;
}



void _DiscreteVoronoi_Construct( void* discreteVoronoi, Stg_ComponentFactory* cf, void* data ) {
	DiscreteVoronoi*	 self          = (DiscreteVoronoi*) discreteVoronoi;
	Dimension_Index      dim;

	dim = Stg_ComponentFactory_GetRootDictUnsignedInt( cf, "dim", 0 );

	_DiscreteVoronoi_Init( self, dim );
}

void _DiscreteVoronoi_Build( void* discreteVoronoi, void* data ) {
}
void _DiscreteVoronoi_Initialise( void* discreteVoronoi, void* data ) {
}
void _DiscreteVoronoi_Execute( void* discreteVoronoi, void* data ) {
}
void _DiscreteVoronoi_Destroy( void* discreteVoronoi, void* data ) {
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void DiscreteVoronoi_CalculateForCell( void* discreteVoronoi, void* _swarm, Cell_LocalIndex lCell_I ) {
	DiscreteVoronoi*	 self          = (DiscreteVoronoi*) discreteVoronoi;

	self->_calculateForCell( self, _swarm, lCell_I );
}

DiscreteVoronoiParticleInfo* DiscreteVoronoi_CreateParticleInfo( void* discreteVoronoi, void* _swarm, Cell_LocalIndex lCell_I ) {
	DiscreteVoronoi*	         self          = (DiscreteVoronoi*) discreteVoronoi;
	Swarm*                       swarm         = (Swarm*)                 _swarm;
	DiscreteVoronoiParticleInfo* particleInfoArray;
	DiscreteVoronoiParticleInfo* currParticleInfo;
	Particle_InCellIndex         cellParticleCount = swarm->cellParticleCountTbl[lCell_I];
	Particle_InCellIndex         cParticle_I;
	Voronoi_CellIndex            vCell_I;
	double                       volume;
	Coord                        centroid;
	Dimension_Index              dim           = self->dim;

	particleInfoArray = Memory_Alloc_Array( DiscreteVoronoiParticleInfo, cellParticleCount, "particle info" );
	memset(particleInfoArray, 0, sizeof( DiscreteVoronoiParticleInfo ) * cellParticleCount );
	
	/* Add volumes and centroids */
	for ( vCell_I = 0 ; vCell_I < self->claimedCellCount ; vCell_I++ ) {
		/* Get Information for this voronoi cell */
		cParticle_I = DiscreteVoronoi_GetParticleIndex( self, vCell_I );
		volume = DiscreteVoronoi_GetVolume( self, vCell_I );
		DiscreteVoronoi_GetCentroid( self, vCell_I, centroid );

		currParticleInfo = &particleInfoArray[ cParticle_I ];
		
		currParticleInfo->volume += volume;

		currParticleInfo->centroid[ I_AXIS ] += volume * centroid[ I_AXIS ];
		currParticleInfo->centroid[ J_AXIS ] += volume * centroid[ J_AXIS ];
		if ( dim == 3 )
			currParticleInfo->centroid[ K_AXIS ] += volume * centroid[ K_AXIS ];

		currParticleInfo->voronoiCellCount++;
	}

	/* Set particle numbers and normalise centroids */
	for ( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {
		currParticleInfo = &particleInfoArray[ cParticle_I ];
	
		currParticleInfo->particle_I = cParticle_I;
		volume = currParticleInfo->volume;
		
		/* Sort out coordinate for centroid */
		if ( currParticleInfo->voronoiCellCount > 0 ) {
			currParticleInfo->centroid[ I_AXIS ] /= volume;
			currParticleInfo->centroid[ J_AXIS ] /= volume;
			if ( dim == 3 )
				currParticleInfo->centroid[ K_AXIS ] /= volume;
		}
		else {
			GlobalParticle* particle = (GlobalParticle*) Swarm_ParticleInCellAt( swarm, lCell_I, cParticle_I );
			/* If this particle hasn't claimed any voronoi cells - 
			 * then the centroid is co-incident with the particle */
			memcpy( currParticleInfo->centroid, particle->coord, dim*sizeof(double) );
		}
	}

	return particleInfoArray;
}

int _DiscreteVoronoiParticleInfo_CompareVolume(const void * _pInfoA, const void * _pInfoB) {
	DiscreteVoronoiParticleInfo* pInfoA = (DiscreteVoronoiParticleInfo*) _pInfoA;
	DiscreteVoronoiParticleInfo* pInfoB = (DiscreteVoronoiParticleInfo*) _pInfoB;

	if ( pInfoA->volume < pInfoB->volume )
		return -1;
	else 
		return 1;
}
	
/* Sorts the particle info list in order of volume - from smallest volume to largest */
void DiscreteVoronoiParticleInfo_SortByVolume( DiscreteVoronoiParticleInfo* particleInfo, Particle_InCellIndex cellParticleCount ) {
	qsort(particleInfo, cellParticleCount, sizeof(DiscreteVoronoiParticleInfo), _DiscreteVoronoiParticleInfo_CompareVolume);
}
