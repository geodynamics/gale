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
** $Id: ManualParticleLayout.c 4081 2007-04-27 06:20:07Z LukeHodkinson $
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
#include "ManualParticleLayout.h"

#include "SwarmClass.h"
#include "StandardParticle.h"
#include "ShadowInfo.h"
#include "CellLayout.h"
#include "ElementCellLayout.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

const Type ManualParticleLayout_Type = "ManualParticleLayout";

const Index ManualParticleLayout_Invalid = (Index)0;

ManualParticleLayout* ManualParticleLayout_New( Name name,
      AbstractContext* context, 
      CoordSystem      coordSystem,
      Bool             weightsInitialisedAtStartup,
      unsigned int     totalInitialParticles, 
      double           averageInitialParticlesPerCell,
      Dictionary*      dictionary ) 
{
	ManualParticleLayout* self = (ManualParticleLayout*) _ManualParticleLayout_DefaultNew( name );

   _ParticleLayout_Init( self, context, coordSystem, weightsInitialisedAtStartup );
   _GlobalParticleLayout_Init( self, totalInitialParticles, averageInitialParticlesPerCell );
	_ManualParticleLayout_Init( self, dictionary );
	return self;
}

ManualParticleLayout* _ManualParticleLayout_New( 
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
      Dictionary*                                      dictionary )
{
	ManualParticleLayout* self;
	
	/* Allocate memory */
	self = (ManualParticleLayout*)_GlobalParticleLayout_New( 
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

   self->dictionary = dictionary;

	return self;
}

void _ManualParticleLayout_Init( void* manualParticleLayout, Dictionary* dictionary )
{
	ManualParticleLayout*	self                            = (ManualParticleLayout*)manualParticleLayout;
	Dictionary_Entry_Value* manualParticlePositions         = NULL;

	self->isConstructed = True;
	self->dictionary    = dictionary;

	manualParticlePositions = Dictionary_Get( self->dictionary, "manualParticlePositions" );
	self->totalInitialParticles = Dictionary_Entry_Value_GetCount( manualParticlePositions );
	self->averageInitialParticlesPerCell = 0;
}

void _ManualParticleLayout_Delete( void* manualParticleLayout ) {
	ManualParticleLayout* self = (ManualParticleLayout*)manualParticleLayout;

	_GlobalParticleLayout_Delete( self );
}

void _ManualParticleLayout_Print( void* manualParticleLayout, Stream* stream ) {
	ManualParticleLayout* self  = (ManualParticleLayout*)manualParticleLayout;
	
	/* General info */
	Journal_Printf( stream, "ManualParticleLayout (ptr): %p:\n", self );
	Stream_Indent( stream );
	
	/* Parent class info */
	_GlobalParticleLayout_Print( self, stream );
	
	/* ManualParticleLayout */
	
	Stream_UnIndent( stream );
}

void* _ManualParticleLayout_Copy( void* manualParticleLayout, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	ManualParticleLayout*		self                    = (ManualParticleLayout*)manualParticleLayout;
	ManualParticleLayout*		newManualParticleLayout;
	
	newManualParticleLayout = (ManualParticleLayout*)_GlobalParticleLayout_Copy( self, dest, deep, nameExt, ptrMap );

	return (void*)newManualParticleLayout;
}

void* _ManualParticleLayout_DefaultNew( Name name ) {
	return (void*)_ManualParticleLayout_New( 
			sizeof(ManualParticleLayout),
			ManualParticleLayout_Type,
			_ManualParticleLayout_Delete,
			_ManualParticleLayout_Print,
			_ManualParticleLayout_Copy,
			_ManualParticleLayout_DefaultNew,
			_ManualParticleLayout_AssignFromXML,
			_ManualParticleLayout_Build,
			_ManualParticleLayout_Initialise,
			_ManualParticleLayout_Execute,
			_ManualParticleLayout_Destroy,
         name, NON_GLOBAL, 
			_GlobalParticleLayout_SetInitialCounts,
			_GlobalParticleLayout_InitialiseParticles,
         GlobalCoordSystem, False,
			_ManualParticleLayout_InitialiseParticle,
         0, 0.0,
			NULL /* dictionary */ );
}

void _ManualParticleLayout_AssignFromXML( void* manualParticleLayout, Stg_ComponentFactory *cf, void* data ) {
	ManualParticleLayout*      self       = (ManualParticleLayout*) manualParticleLayout;
	Dictionary*                dictionary;
	
   _GlobalParticleLayout_AssignFromXML( self, cf, data );

   dictionary = Dictionary_GetDictionary( cf->componentDict, self->name );

   _ManualParticleLayout_Init( self, dictionary );
}
	
void _ManualParticleLayout_Build( void* manualParticleLayout, void* data ) {
}
void _ManualParticleLayout_Initialise( void* manualParticleLayout, void* data ) {
}
void _ManualParticleLayout_Execute( void* manualParticleLayout, void* data ) {
}
void _ManualParticleLayout_Destroy( void* manualParticleLayout, void* data ) {
	ManualParticleLayout* self = (ManualParticleLayout*)manualParticleLayout;

   _GlobalParticleLayout_Destroy( self, data );
}

void _ManualParticleLayout_InitialiseParticle( 
		void* manualParticleLayout, 
		void* _swarm, 
		Particle_Index newParticle_I,
		void* _particle )
{
	ManualParticleLayout*       self                    = (ManualParticleLayout*)manualParticleLayout;
	Dictionary_Entry_Value*     manualParticlePositions = NULL;
	Dictionary_Entry_Value*     particlePositionEntry   = NULL;
	Dictionary*                 particlePositionDict    = NULL;
	GlobalParticle*             particle                = (GlobalParticle*)_particle;

	manualParticlePositions = Dictionary_Get( self->dictionary, "manualParticlePositions" );
	particlePositionEntry = Dictionary_Entry_Value_GetElement( manualParticlePositions, newParticle_I );
	particlePositionDict = Dictionary_Entry_Value_AsDictionary( particlePositionEntry );
	particle->coord[I_AXIS] = Dictionary_GetDouble_WithDefault( particlePositionDict, "x", 0.0 );
	particle->coord[J_AXIS] = Dictionary_GetDouble_WithDefault( particlePositionDict, "y", 0.0 );
	particle->coord[K_AXIS] = Dictionary_GetDouble_WithDefault( particlePositionDict, "z", 0.0 );
}

