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
** $Id: PlaneParticleLayout.c 3851 2006-10-12 08:57:22Z SteveQuenette $
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
#include "PlaneParticleLayout.h"
#include "CellLayout.h"
#include "SwarmClass.h"
#include "StandardParticle.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

const Type PlaneParticleLayout_Type = "PlaneParticleLayout";

const Index PlaneParticleLayout_Invalid = (Index) 0;

PlaneParticleLayout* PlaneParticleLayout_New( 
      Name             name,
      AbstractContext* context, 
      CoordSystem      coordSystem,
      Bool             weightsInitialisedAtStartup,
      unsigned int     totalInitialParticles, 
      double           averageInitialParticlesPerCell,
      Dimension_Index  dim,
      Axis             planeAxis, 
      double           planeCoord )
{
	PlaneParticleLayout* self = (PlaneParticleLayout*) _PlaneParticleLayout_DefaultNew( name );

   _ParticleLayout_Init( self, context, coordSystem, weightsInitialisedAtStartup );
   _GlobalParticleLayout_Init( self, totalInitialParticles, averageInitialParticlesPerCell );
	_SpaceFillerParticleLayout_Init( self, dim );
	_PlaneParticleLayout_Init( self, planeAxis, planeCoord );
	
	return self;
}

PlaneParticleLayout* _PlaneParticleLayout_New(  PLANEPARTICLELAYOUT_DEFARGS  )
{
	PlaneParticleLayout* self;
	
	/* Allocate memory */
	self = (PlaneParticleLayout*)_SpaceFillerParticleLayout_New(  SPACEFILLERPARTICLELAYOUT_PASSARGS  );  /* dim */

   self->planeAxis = planeAxis;
   self->planeCoord = planeCoord;
	
	return self;
}

void _PlaneParticleLayout_Init( 
		void*           particleLayout, 
		Axis            planeAxis, 
		double          planeCoord )
{
	PlaneParticleLayout* self = (PlaneParticleLayout*) particleLayout;

	self->planeAxis  = planeAxis;
	self->planeCoord = planeCoord;
}


void* _PlaneParticleLayout_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                                _sizeOfSelf = sizeof(PlaneParticleLayout);
	Type                                                                        type = PlaneParticleLayout_Type;
	Stg_Class_DeleteFunction*                                                _delete = _PlaneParticleLayout_Delete;
	Stg_Class_PrintFunction*                                                  _print = _PlaneParticleLayout_Print;
	Stg_Class_CopyFunction*                                                    _copy = _PlaneParticleLayout_Copy;
	Stg_Component_DefaultConstructorFunction*                    _defaultConstructor = _PlaneParticleLayout_DefaultNew;
	Stg_Component_ConstructFunction*                                      _construct = _PlaneParticleLayout_AssignFromXML;
	Stg_Component_BuildFunction*                                              _build = _PlaneParticleLayout_Build;
	Stg_Component_InitialiseFunction*                                    _initialise = _PlaneParticleLayout_Initialise;
	Stg_Component_ExecuteFunction*                                          _execute = _PlaneParticleLayout_Execute;
	Stg_Component_DestroyFunction*                                          _destroy = _PlaneParticleLayout_Destroy;
	AllocationType                                                nameAllocationType = NON_GLOBAL;
	ParticleLayout_SetInitialCountsFunction*                       _setInitialCounts = _GlobalParticleLayout_SetInitialCounts;
	ParticleLayout_InitialiseParticlesFunction*                 _initialiseParticles = _SpaceFillerParticleLayout_InitialiseParticles;
	CoordSystem                                                          coordSystem = GlobalCoordSystem;
	Bool                                                 weightsInitialisedAtStartup = False;
	GlobalParticleLayout_InitialiseParticleFunction*             _initialiseParticle = _PlaneParticleLayout_InitialiseParticle;
	Particle_Index                                             totalInitialParticles = 0;
	double                                            averageInitialParticlesPerCell = 0.0;
	Dimension_Index                                                              dim = 0;
	Axis                                                                   planeAxis = 0;
	double                                                                planeCoord = 0;

   return (void*)_PlaneParticleLayout_New(  PLANEPARTICLELAYOUT_PASSARGS  );
}

	
void _PlaneParticleLayout_Destroy( void* particleLayout, void* data ){
   PlaneParticleLayout* self = (PlaneParticleLayout*)particleLayout;
   _SpaceFillerParticleLayout_Destroy( self, data );
}

void _PlaneParticleLayout_Delete( void* particleLayout ) {
	PlaneParticleLayout* self = (PlaneParticleLayout*)particleLayout;

	/* Stg_Class_Delete parent class */
	_SpaceFillerParticleLayout_Delete( self );

}

void _PlaneParticleLayout_Print( void* particleLayout, Stream* stream ) {
	PlaneParticleLayout* self  = (PlaneParticleLayout*)particleLayout;
	
	/* General info */
	Journal_Printf( stream, "PlaneParticleLayout (ptr): %p:\n", self );
	Stream_Indent( stream );
	
	/* Parent class info */
	_SpaceFillerParticleLayout_Print( self, stream );
	
	/* PlaneParticleLayout */
	Journal_PrintValue( stream, self->planeAxis );
	Journal_PrintValue( stream, self->planeCoord );
	
	Stream_UnIndent( stream );
}


void* _PlaneParticleLayout_Copy( void* particleLayout, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	PlaneParticleLayout*		self                    = (PlaneParticleLayout*)particleLayout;
	PlaneParticleLayout*		newPlaneParticleLayout;
	
	newPlaneParticleLayout = _SpaceFillerParticleLayout_Copy( self, dest, deep, nameExt, ptrMap );
	
	newPlaneParticleLayout->planeAxis  = self->planeAxis;
	newPlaneParticleLayout->planeCoord = self->planeCoord;

	return (void*)newPlaneParticleLayout;
}


void _PlaneParticleLayout_AssignFromXML( void* particleLayout, Stg_ComponentFactory *cf, void* data ){
	PlaneParticleLayout* self = (PlaneParticleLayout*) particleLayout;
	Axis   planeAxis;
	double planeCoord;
	char*  planeAxisString;
	
	_SpaceFillerParticleLayout_AssignFromXML( self, cf, data );

	planeAxisString = Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"planeAxis", ""  );
	planeCoord = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"planeCoord", 0.0 );

	/* Check to make sure that some value is given for plane axis */
	Journal_Firewall( strlen( planeAxisString  ) > 0, Journal_MyStream( Error_Type, self ),
		"Error for %s '%s': No axis given in param 'planeAxis'.\n", self->type, self->name );

	/* Make axis case insensitive */
	planeAxisString[0] = toupper(planeAxisString[0]);
	Journal_Firewall( planeAxisString[0] >= 'X' && planeAxisString[0] <= 'Z', Journal_MyStream( Error_Type, self ),
		"Error for %s '%s': Incorrect axis '%c' given for param 'planeAxis'.\n", self->type, self->name,planeAxisString[0]);

	planeAxis = planeAxisString[0] - 'X';

	_PlaneParticleLayout_Init( self, planeAxis, planeCoord );
	
}
	
void _PlaneParticleLayout_Build( void* particleLayout, void* data ){
   PlaneParticleLayout* self = (PlaneParticleLayout*)particleLayout;
   _SpaceFillerParticleLayout_Build( self, data );
}
	
void _PlaneParticleLayout_Initialise( void* particleLayout, void* data ){
   PlaneParticleLayout* self = (PlaneParticleLayout*)particleLayout;
   _SpaceFillerParticleLayout_Initialise( self, data );
}
	
void _PlaneParticleLayout_Execute( void* particleLayout, void* data ){
	
}

void _PlaneParticleLayout_InitialiseParticle( 
		void* spaceFillerParticleLayout, 
		void* _swarm, 
		Particle_Index newParticle_I,
		void* _particle )
{
	PlaneParticleLayout*        self     = (PlaneParticleLayout*)spaceFillerParticleLayout;
	Swarm*                      swarm    = (Swarm*)_swarm;
	GlobalParticle*             particle = (GlobalParticle*)_particle;

	_SpaceFillerParticleLayout_InitialiseParticle( self, swarm, newParticle_I, particle );
	particle->coord[ self->planeAxis ] = self->planeCoord;
}



