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
** $Id: DiscreteVoronoiRemove.c 518 2007-10-11 08:07:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include <PICellerator/Voronoi/Voronoi.h>

#include "types.h"
#include "RemovalRoutine.h"
#include "DiscreteVoronoiRemove.h"

#include <assert.h>
#include <string.h>
#include <math.h>


/* Textual name of this class */
const Type DiscreteVoronoiRemove_Type = "DiscreteVoronoiRemove";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/
DiscreteVoronoiRemove* _DiscreteVoronoiRemove_New(
		SizeT                                   _sizeOfSelf, 
		Type                                    type,
		Stg_Class_DeleteFunction*               _delete,
		Stg_Class_PrintFunction*                _print,
		Stg_Class_CopyFunction*                 _copy, 
		Stg_Component_DefaultConstructorFunction*   _defaultConstructor,
		Stg_Component_ConstructFunction*        _construct,
		Stg_Component_BuildFunction*            _build,
		Stg_Component_InitialiseFunction*       _initialise,
		Stg_Component_ExecuteFunction*          _execute,
		Stg_Component_DestroyFunction*          _destroy,		
		RemovalRoutine_RemoveFromCellFunction*  _removeFromCell,
		Name                                    name )
{
	DiscreteVoronoiRemove* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(DiscreteVoronoiRemove) );
	self = (DiscreteVoronoiRemove*)_RemovalRoutine_New( 
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
			_removeFromCell,
			name );

	
	/* General info */

	/* Virtual Info */
	
	return self;
}

void _DiscreteVoronoiRemove_Init( void* discreteVoronoiRemove, DiscreteVoronoi* discreteVoronoi ) {
	DiscreteVoronoiRemove* self = (DiscreteVoronoiRemove*)discreteVoronoiRemove;
	
	self->isConstructed = True;
	self->discreteVoronoi = discreteVoronoi;
}

/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/
void* _DiscreteVoronoiRemove_DefaultNew( Name name ) {
	return (void*) _DiscreteVoronoiRemove_New(
			sizeof(DiscreteVoronoiRemove),
			DiscreteVoronoiRemove_Type,
			_RemovalRoutine_Delete,
			_RemovalRoutine_Print,
			_RemovalRoutine_Copy,
			_DiscreteVoronoiRemove_DefaultNew,
			_DiscreteVoronoiRemove_Construct,
			_DiscreteVoronoiRemove_Build,
			_DiscreteVoronoiRemove_Initialise,
			_RemovalRoutine_Execute,
			_RemovalRoutine_Destroy,
			_DiscreteVoronoiRemove_RemoveFromCell,
			name );
}


void _DiscreteVoronoiRemove_Construct( void* discreteVoronoiRemove, Stg_ComponentFactory* cf, void* data ) {
	DiscreteVoronoiRemove*	     self            = (DiscreteVoronoiRemove*) discreteVoronoiRemove;
	DiscreteVoronoi*             discreteVoronoi;

	_RemovalRoutine_Construct( self, cf, data );

	discreteVoronoi =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "DiscreteVoronoi", DiscreteVoronoi, True, data );
	_DiscreteVoronoiRemove_Init( self, discreteVoronoi );
}

void _DiscreteVoronoiRemove_Build( void* discreteVoronoiRemove, void* data ) {
	DiscreteVoronoiRemove*	     self            = (DiscreteVoronoiRemove*) discreteVoronoiRemove;

	/* Ensure that the discrete voronoi object I'm using is built */
	Stg_Component_Build( self->discreteVoronoi, data, False );
	_RemovalRoutine_Build( self, data );
}

void _DiscreteVoronoiRemove_Initialise( void* discreteVoronoiRemove, void* data ) {
	DiscreteVoronoiRemove*       self            = (DiscreteVoronoiRemove*) discreteVoronoiRemove;

	/* Ensure that the discrete voronoi object I'm using is initialised */
	Stg_Component_Initialise( self->discreteVoronoi, data, False );
	_RemovalRoutine_Initialise( self, data );
}

/*-------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
void _DiscreteVoronoiRemove_RemoveFromCell( void* discreteVoronoiRemove, void* _swarm, Cell_LocalIndex lCell_I ) {
	DiscreteVoronoiRemove*	     self            = (DiscreteVoronoiRemove*) discreteVoronoiRemove;
	Swarm*                       swarm           = (Swarm*)                 _swarm;
	DiscreteVoronoiParticleInfo* voronoiParticleInfo;
	Particle_InCellIndex         cellParticleCount = swarm->cellParticleCountTbl[ lCell_I ];
	Particle_InCellIndex         cParticle_I;
	Particle_InCellIndex         particlesToRemoveCount = cellParticleCount - self->idealParticleCount;

	DiscreteVoronoi_CalculateForCell( self->discreteVoronoi, swarm, lCell_I );
	voronoiParticleInfo = DiscreteVoronoi_CreateParticleInfo( self->discreteVoronoi, swarm, lCell_I );
	DiscreteVoronoiParticleInfo_SortByVolume( voronoiParticleInfo, cellParticleCount );

	/* Go through particles with smallest area and set them for removal */
	for ( cParticle_I = 0 ; cParticle_I < particlesToRemoveCount ; cParticle_I++ ) 
		RemovalRoutine_SetParticleToRemove( self, swarm, lCell_I, voronoiParticleInfo[ cParticle_I ].particle_I );

	DiscreteVoronoiParticleInfo_Delete( voronoiParticleInfo );
}
