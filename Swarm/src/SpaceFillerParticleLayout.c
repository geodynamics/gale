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
** $Id: SpaceFillerParticleLayout.c 4081 2007-04-27 06:20:07Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>

#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>
#include <StgDomain/Utils/Utils.h>

#include "types.h"
#include "shortcuts.h"
#include "ParticleLayout.h"
#include "GlobalParticleLayout.h"
#include "SpaceFillerParticleLayout.h"

#include "SwarmClass.h"
#include "StandardParticle.h"
#include "ShadowInfo.h"
#include "CellLayout.h"
#include "ElementCellLayout.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

const Type SpaceFillerParticleLayout_Type = "SpaceFillerParticleLayout";

const Index SpaceFillerParticleLayout_Invalid = (Index) 0;

SpaceFillerParticleLayout* SpaceFillerParticleLayout_New( 
		Name             name,
      AbstractContext* context, 
      CoordSystem      coordSystem,
      Bool             weightsInitialisedAtStartup,
      unsigned int     totalInitialParticles, 
      double           averageInitialParticlesPerCell,
		Dimension_Index  dim )
{
	SpaceFillerParticleLayout* self = (SpaceFillerParticleLayout*) _SpaceFillerParticleLayout_DefaultNew( name );

   _ParticleLayout_Init( self, context, coordSystem, weightsInitialisedAtStartup );
   _GlobalParticleLayout_Init( self, totalInitialParticles, averageInitialParticlesPerCell );
	_SpaceFillerParticleLayout_Init( self, dim );
	return self;
}

SpaceFillerParticleLayout* _SpaceFillerParticleLayout_New( 
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
      Name                                             name,
      AllocationType                                   nameAllocationType,
      ParticleLayout_SetInitialCountsFunction*         _setInitialCounts,
      ParticleLayout_InitialiseParticlesFunction*      _initialiseParticles,
      CoordSystem                                      coordSystem,
      Bool                                             weightsInitialisedAtStartup,
      GlobalParticleLayout_InitialiseParticleFunction* _initialiseParticle,
      Particle_Index                                   totalInitialParticles,
      double                                           averageInitialParticlesPerCell,
      Dimension_Index                                  dim )
{
   SpaceFillerParticleLayout* self;

   /* Allocate memory */
   self = (SpaceFillerParticleLayout*)_GlobalParticleLayout_New( 
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
      name,
      nameAllocationType,
      _setInitialCounts,
      _initialiseParticles,
      coordSystem,
      weightsInitialisedAtStartup,
      _initialiseParticle,
      totalInitialParticles,
      averageInitialParticlesPerCell );

   self->dim = dim;

	return self;
}

void _SpaceFillerParticleLayout_Init( 
		void*                   spaceFillerParticleLayout,
		Dimension_Index         dim )
{
	SpaceFillerParticleLayout*	self = (SpaceFillerParticleLayout*) spaceFillerParticleLayout;

	self->isConstructed = True;
	self->dim           = dim;

	self->sobolGenerator[I_AXIS] = SobolGenerator_NewFromTable( "xSobolGenerator" );
	self->sobolGenerator[J_AXIS] = SobolGenerator_NewFromTable( "ySobolGenerator" );
	if ( dim == 3 )
		self->sobolGenerator[K_AXIS] = SobolGenerator_NewFromTable( "zSobolGenerator" );

	/* Must set one or the other. Fail if both set, or none */
	Journal_Firewall( 
		(self->totalInitialParticles == SpaceFillerParticleLayout_Invalid)
			^ ((Index)self->averageInitialParticlesPerCell == SpaceFillerParticleLayout_Invalid),
		Journal_MyStream( Error_Type, self ),
		"Error in func %s for %s '%s' - Both averageInitialParticlesPerCell and totalInitialParticles%sspecified.\n",
		__func__, 
		self->type, 
		self->name,
		(self->totalInitialParticles == SpaceFillerParticleLayout_Invalid) ? " not " : " " );
	
}


	
void _SpaceFillerParticleLayout_Delete( void* spaceFillerParticleLayout ) {
	SpaceFillerParticleLayout* self = (SpaceFillerParticleLayout*)spaceFillerParticleLayout;
	Dimension_Index            dim_I;

	_GlobalParticleLayout_Delete( self );
}

void _SpaceFillerParticleLayout_Print( void* spaceFillerParticleLayout, Stream* stream ) {
	SpaceFillerParticleLayout* self  = (SpaceFillerParticleLayout*)spaceFillerParticleLayout;
	
	/* General info */
	Journal_Printf( stream, "SpaceFillerParticleLayout (ptr): %p:\n", self );
	Stream_Indent( stream );
	
	/* Parent class info */
	_GlobalParticleLayout_Print( self, stream );
	
	/* SpaceFillerParticleLayout */
	Journal_PrintValue( stream, self->dim );
	
	Stream_UnIndent( stream );
}

void* _SpaceFillerParticleLayout_Copy( void* spaceFillerParticleLayout, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	SpaceFillerParticleLayout*		self                    = (SpaceFillerParticleLayout*)spaceFillerParticleLayout;
	SpaceFillerParticleLayout*		newSpaceFillerParticleLayout;
	
	newSpaceFillerParticleLayout = (SpaceFillerParticleLayout*)_GlobalParticleLayout_Copy( self, dest, deep, nameExt, ptrMap );
	
	newSpaceFillerParticleLayout->dim = self->dim;

	return (void*)newSpaceFillerParticleLayout;
}

void* _SpaceFillerParticleLayout_DefaultNew( Name name ) {
	return (void*)_SpaceFillerParticleLayout_New( 
			sizeof(SpaceFillerParticleLayout),
			SpaceFillerParticleLayout_Type,
			_SpaceFillerParticleLayout_Delete,
			_SpaceFillerParticleLayout_Print,
			_SpaceFillerParticleLayout_Copy,
			_SpaceFillerParticleLayout_DefaultNew,
			_SpaceFillerParticleLayout_AssignFromXML,
			_SpaceFillerParticleLayout_Build,
			_SpaceFillerParticleLayout_Initialise,
			_SpaceFillerParticleLayout_Execute,
			_SpaceFillerParticleLayout_Destroy,
         name, NON_GLOBAL, 
			_GlobalParticleLayout_SetInitialCounts,
			_SpaceFillerParticleLayout_InitialiseParticles,
         GlobalCoordSystem, False,
			_SpaceFillerParticleLayout_InitialiseParticle,
         0, 0.0,
			0  /* dim */ );
}


void _SpaceFillerParticleLayout_AssignFromXML( void* spaceFillerParticleLayout, Stg_ComponentFactory *cf, void* data ) {
   SpaceFillerParticleLayout* self = (SpaceFillerParticleLayout*) spaceFillerParticleLayout;
   Dimension_Index            dim;

   _GlobalParticleLayout_AssignFromXML( self, cf, data );

   dim = Stg_ComponentFactory_GetRootDictUnsignedInt( cf, "dim", 0 );

   _SpaceFillerParticleLayout_Init( self, dim );
}
	
void _SpaceFillerParticleLayout_Build( void* spaceFillerParticleLayout, void* data ) {
}
void _SpaceFillerParticleLayout_Initialise( void* spaceFillerParticleLayout, void* data ) {
}	
void _SpaceFillerParticleLayout_Execute( void* spaceFillerParticleLayout, void* data ) {	
}
void _SpaceFillerParticleLayout_Destroy( void* spaceFillerParticleLayout, void* data ) {	
   SpaceFillerParticleLayout* self = (SpaceFillerParticleLayout*)spaceFillerParticleLayout;
   int dim_I;
   
   for ( dim_I = 0 ; dim_I < self->dim ; dim_I++ ) 
      Stg_Class_Delete( self->sobolGenerator[ dim_I ] );

   _GlobalParticleLayout_Destroy( self, data );
}

void _SpaceFillerParticleLayout_InitialiseParticles( void* spaceFillerParticleLayout, void* swarm ) {
	SpaceFillerParticleLayout*        self     = (SpaceFillerParticleLayout*)spaceFillerParticleLayout;

	/* Initialise random number generator */
	_GlobalParticleLayout_InitialiseParticles( self, swarm );
}

void _SpaceFillerParticleLayout_InitialiseParticle( 
		void*             spaceFillerParticleLayout, 
		void*             _swarm, 
		Particle_Index    newParticle_I,
		void*             _particle )
{
	SpaceFillerParticleLayout*  self     = (SpaceFillerParticleLayout*)spaceFillerParticleLayout;
	Swarm*                      swarm    = (Swarm*)_swarm;
	double				minCrd[3];
	double				maxCrd[3];
	Dimension_Index             dim_I;
	double*                     coord;
	GlobalParticle*             particle = (GlobalParticle*)_particle;
	/* Note: we have to assume element cell layout to get the full box coords - maybe these
	should be available as part of all cell layouts */
	ElementCellLayout*          eCellLayout = Stg_CheckType( swarm->cellLayout, ElementCellLayout );
	Mesh*                       mesh = eCellLayout->mesh;

	Mesh_GetGlobalCoordRange( mesh, minCrd, maxCrd );

	coord = particle->coord;
	memset( coord, 0, sizeof(Coord) );

	for ( dim_I = 0; dim_I < self->dim; dim_I++ ) {
		coord[ dim_I ] = SobolGenerator_GetNextNumber_WithMinMax(
			self->sobolGenerator[ dim_I ], 
			minCrd[ dim_I ], 
			maxCrd[ dim_I ]);
	}
}
