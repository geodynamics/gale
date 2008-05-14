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
** $Id: PeriodicBoundariesManager.c 569 2008-05-14 02:42:54Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include <PICellerator/Voronoi/Voronoi.h>
#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>

#include "types.h"
#include "PeriodicBoundariesManager.h"

#include <string.h>
#include <assert.h>

char IJKTopology_DimNumToDimLetter[3] = {'I', 'J', 'K'};

/* Textual name of this class */
const Type PeriodicBoundariesManager_Type = "PeriodicBoundariesManager";

/* Constructors ------------------------------------------------------------------------------------------------*/
void* _PeriodicBoundariesManager_DefaultNew( Name name ) {
	return (void*) _PeriodicBoundariesManager_New(
		sizeof(PeriodicBoundariesManager),
		PeriodicBoundariesManager_Type,
		_PeriodicBoundariesManager_Delete,
		_PeriodicBoundariesManager_Print,
		NULL, 
		_PeriodicBoundariesManager_DefaultNew,
		_PeriodicBoundariesManager_Construct,
		_PeriodicBoundariesManager_Build,
		_PeriodicBoundariesManager_Initialise,
		_PeriodicBoundariesManager_Execute,
		_PeriodicBoundariesManager_Destroy,
		name,
		False,
		NULL,
		NULL,
		NULL);
}

PeriodicBoundariesManager* PeriodicBoundariesManager_New( 
		Name                        name,
		Mesh*			    mesh, 
		Swarm*                      swarm,
		Dictionary*                 dictionary )
{
	return _PeriodicBoundariesManager_New(
		sizeof(PeriodicBoundariesManager),
		PeriodicBoundariesManager_Type,
		_PeriodicBoundariesManager_Delete,
		_PeriodicBoundariesManager_Print,
		_PeriodicBoundariesManager_Copy, 
		_PeriodicBoundariesManager_DefaultNew,
		_PeriodicBoundariesManager_Construct,
		_PeriodicBoundariesManager_Build,
		_PeriodicBoundariesManager_Initialise,
		_PeriodicBoundariesManager_Execute,
		_PeriodicBoundariesManager_Destroy,
		name,
		True,
		mesh, 
		swarm,
		dictionary );
}	


PeriodicBoundariesManager* _PeriodicBoundariesManager_New( 
		SizeT                                  sizeOfSelf,
		Type                                   type,
		Stg_Class_DeleteFunction*              _delete,
		Stg_Class_PrintFunction*               _print,
		Stg_Class_CopyFunction*                _copy, 
		Stg_Component_DefaultConstructorFunction*  _defaultConstructor,
		Stg_Component_ConstructFunction*       _construct,
		Stg_Component_BuildFunction*           _build,
		Stg_Component_InitialiseFunction*      _initialise,
		Stg_Component_ExecuteFunction*         _execute,
		Stg_Component_DestroyFunction*         _destroy,
		Name                                   name,
		Bool                                   initFlag,
		Mesh*				       mesh, 
		Swarm*                                 swarm,
		Dictionary*                            dictionary )		
{
	PeriodicBoundariesManager* self;
	
	/* Allocate memory */
	self = (PeriodicBoundariesManager*)_Stg_Component_New( 
		sizeOfSelf, 
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
		initFlag );
	
	/* General info */
	
	/* Virtual info */
	
	if( initFlag ){
		_PeriodicBoundariesManager_Init( self, mesh, swarm, dictionary );
	}
	
	return self;
}


void _PeriodicBoundariesManager_Init(
		void*             periodicBCsManager,
		Mesh*		  mesh, 
		Swarm*            swarm,
		Dictionary*       dictionary )
{
	PeriodicBoundariesManager* self = (PeriodicBoundariesManager*)periodicBCsManager;

	self->isConstructed = True;
	self->dictionary = dictionary;
	self->mesh = mesh;
	self->swarm = swarm;
	self->count = 0;
	self->delta = 0;
	self->size = 0;
	self->boundaries = NULL;
	self->debug = Journal_Register( Debug_Type, self->type );
}


void _PeriodicBoundariesManager_Construct( void* periodicBCsManager, Stg_ComponentFactory* cf, void* data ) {
	PeriodicBoundariesManager*	self = (PeriodicBoundariesManager*)periodicBCsManager;
	Dictionary*			dictionary = NULL;
	Mesh*				mesh = NULL;
	Swarm*                          swarm = NULL;

	dictionary = Dictionary_GetDictionary( cf->componentDict, self->name );
	mesh =  Stg_ComponentFactory_ConstructByKey(  cf,  self->name,  "mesh", Mesh,  True, data  ) ;
	swarm =  Stg_ComponentFactory_ConstructByKey(  cf,  self->name,  "Swarm", Swarm,  True, data  ) ;

	_PeriodicBoundariesManager_Init( self, mesh, swarm, dictionary );
}


void _PeriodicBoundariesManager_Delete( void* perBCsManager ) {
	PeriodicBoundariesManager* self = (PeriodicBoundariesManager*)perBCsManager;
	
	Memory_Free( self->boundaries );
	
	/* Stg_Class_Delete parent */
	_Stg_Component_Delete( self );
}

/* Virtual Functions -------------------------------------------------------------------------------------------------------------*/

void _PeriodicBoundariesManager_Print( void* perBCsManager, Stream* stream ) {
	PeriodicBoundariesManager* self = (PeriodicBoundariesManager*)perBCsManager;
	Index		perBoundary_I = 0;
	
	/* General info */
	Journal_Printf( stream, "PeriodicBoundariesManager (ptr): %p\n", self );
	
	/* Print parent */
	_Stg_Component_Print( self, stream );

	Journal_Printf( stream, "%d periodic boundaries registered: %p\n", self );
	Stream_Indent( stream );
	for ( perBoundary_I = 0; perBoundary_I < self->count; perBoundary_I++ ) {
		Journal_Printf( stream, "Boundary %d: Axis %d, Min=%f, Max=%f\n", perBoundary_I,
			self->boundaries[perBoundary_I].axis,
			self->boundaries[perBoundary_I].minWall,
			self->boundaries[perBoundary_I].maxWall );
	}
	Stream_UnIndent( stream );
}


void* _PeriodicBoundariesManager_Copy( void* periodicBCsManager, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	PeriodicBoundariesManager*   self = (PeriodicBoundariesManager*)periodicBCsManager;
	PeriodicBoundariesManager*   newPeriodicBCsManager;
	PtrMap*	                     map = ptrMap;
	Bool                         ownMap = False;
	
	if( !map ) {
		map = PtrMap_New( 10 );
		ownMap = True;
	}

	newPeriodicBCsManager = _Stg_Class_Copy( self, dest, deep, nameExt, map );

	newPeriodicBCsManager->count = self->count;
	newPeriodicBCsManager->size = self->size;
	newPeriodicBCsManager->delta = self->delta;

	if ( deep ) {
		newPeriodicBCsManager->dictionary = (Dictionary*)Stg_Class_Copy( self->dictionary, NULL, deep, nameExt, map );
		newPeriodicBCsManager->mesh = (Mesh*)Stg_Class_Copy( self->mesh, NULL, deep, nameExt, map );
		newPeriodicBCsManager->swarm = (Swarm*)Stg_Class_Copy( self->swarm, NULL, deep, nameExt, map );
		newPeriodicBCsManager->debug = self->debug;
		newPeriodicBCsManager->boundaries = Memory_Alloc_Array( PeriodicBoundary, self->size,
			"PeriodicBoundaries" );
		memcpy( newPeriodicBCsManager->boundaries, self->boundaries, sizeof(PeriodicBoundary)*self->count );	
	}
	else {
		newPeriodicBCsManager->dictionary = self->dictionary;
		newPeriodicBCsManager->mesh = self->mesh;
		newPeriodicBCsManager->swarm = self->swarm;
		newPeriodicBCsManager->boundaries = self->boundaries;
		newPeriodicBCsManager->debug = self->debug;
	}
	
	if( ownMap ) {
		Stg_Class_Delete( map );
	}
	
	return (void*)newPeriodicBCsManager;


	
	return NULL;
}


void _PeriodicBoundariesManager_Build( void* periodicBCsManager, void* data ) {	
	PeriodicBoundariesManager* self = (PeriodicBoundariesManager*)periodicBCsManager;
	Dictionary_Entry_Value*    periodicBCsList = NULL;
	
	self->size = 4;
	self->boundaries = Memory_Alloc_Array( PeriodicBoundary, self->size, "PeriodicBoundariesManager->boundaries" );

	if ( self->dictionary ) {
		periodicBCsList = Dictionary_Get( self->dictionary, "PeriodicBoundaries" );
		
		/* Dictionary entry is optional - users may prefer to enter in code */
		if ( periodicBCsList ) {
			Index                   numPeriodicBCs = 0;
			Index                   periodicBC_I = 0;
			Dictionary_Entry_Value* periodicBC = NULL;
			char*                   perBCAxis = NULL;
			
			numPeriodicBCs = Dictionary_Entry_Value_GetCount( periodicBCsList );

			for ( periodicBC_I = 0; periodicBC_I < numPeriodicBCs; periodicBC_I++ ) {
				periodicBC = Dictionary_Entry_Value_GetElement( periodicBCsList, periodicBC_I );
				perBCAxis = Dictionary_Entry_Value_AsString( periodicBC );

				if ( 0 == strcmp( perBCAxis, "I_AXIS" ) ) {
					PeriodicBoundariesManager_AddPeriodicBoundary( self, I_AXIS );
				}
				else if ( 0 == strcmp( perBCAxis, "J_AXIS" ) ) {
					PeriodicBoundariesManager_AddPeriodicBoundary( self, J_AXIS );
				}
				else if ( 0 == strcmp( perBCAxis, "K_AXIS" ) ) {
					PeriodicBoundariesManager_AddPeriodicBoundary( self, K_AXIS );
				}
			}
		}
	}
	/* Test if mesh is periodic */
	else if ( Stg_Class_IsInstance( self->mesh->generator, CartesianGenerator_Type ) ) {
		CartesianGenerator* cartesianGenerator = (CartesianGenerator*) self->mesh->generator;
		Dimension_Index dim_I;

		for ( dim_I = 0 ; dim_I < self->swarm->dim ; dim_I++ ) {
			/* Add boundaries straight from mesh generator */
			if ( cartesianGenerator->periodic[ dim_I ] ) 
				PeriodicBoundariesManager_AddPeriodicBoundary( self, dim_I );
		}		
	}



}


void _PeriodicBoundariesManager_Initialise( void* periodicBCsManager, void* data ) {	
}

void _PeriodicBoundariesManager_Execute( void* periodicBCsManager, void* data ) {	
}

void _PeriodicBoundariesManager_Destroy( void* periodicBCsManager, void* data ) {	
}

/* Public Functions -------------------------------------------------------------------------------------------------------------*/

void PeriodicBoundariesManager_AddPeriodicBoundary( void* periodicBCsManager, Axis axis ) {
	PeriodicBoundariesManager*	self = (PeriodicBoundariesManager*)periodicBCsManager;	
	PeriodicBoundary*		newPeriodicBoundary;
	double				min[3], max[3];

	Mesh_GetGlobalCoordRange( self->mesh, min, max );
	
	if ( self->count == self->size ) {
		self->size += self->delta;
		self->boundaries = Memory_Realloc_Array( self->boundaries, PeriodicBoundary, self->size );
	}
	newPeriodicBoundary = &self->boundaries[self->count];
	newPeriodicBoundary->axis = axis;
	newPeriodicBoundary->minWall = min[axis];
	newPeriodicBoundary->maxWall = max[axis];
	newPeriodicBoundary->particlesUpdatedMinEndCount = 0;	
	newPeriodicBoundary->particlesUpdatedMaxEndCount = 0;	
	self->count++;
}


void PeriodicBoundariesManager_UpdateParticle( void* periodicBCsManager, Particle_Index lParticle_I ) {
	Axis				boundaryAxis;	
	PeriodicBoundariesManager*	self = (PeriodicBoundariesManager*)periodicBCsManager;
	double				difference = 0.0;
	GlobalParticle*                 particle = NULL;
	Index				perBoundary_I = 0;
	PeriodicBoundary*		perBoundary = NULL;

	Journal_DPrintfL( self->debug, 2, "In %s:\n", __func__ );
	Stream_Indent( self->debug );

	particle = (GlobalParticle*)Swarm_ParticleAt( self->swarm, lParticle_I );

	Journal_DPrintfL( self->debug, 2, "Checking particle %d at (%.4g,%.4g,%.4g)\n", lParticle_I,
		particle->coord[0], particle->coord[1], particle->coord[2] );

	for ( perBoundary_I = 0; perBoundary_I < self->count; perBoundary_I++ ) {

		perBoundary = &self->boundaries[perBoundary_I];
		boundaryAxis = perBoundary->axis;

		Journal_DPrintfL( self->debug, 2, "Checking axis %d:\n", boundaryAxis );

			
		Stream_Indent( self->debug );
		if ( particle->coord[boundaryAxis] < perBoundary->minWall ) {
			Journal_DPrintfL( self->debug, 3, "coord is < min wall %.4f:\n", perBoundary->minWall );
			difference = perBoundary->minWall - particle->coord[boundaryAxis];
			particle->coord[boundaryAxis] = perBoundary->maxWall - difference;
			perBoundary->particlesUpdatedMinEndCount++;
			Journal_DPrintfL( self->debug, 3, "moving to (%.4f,%.4f,%.4f).\n",
				particle->coord[I_AXIS], particle->coord[J_AXIS],
				particle->coord[K_AXIS] );
		}
		else if ( particle->coord[perBoundary->axis] > perBoundary->maxWall ) {
			Journal_DPrintfL( self->debug, 3, "coord is > max wall %.4f:\n", perBoundary->maxWall );
			difference = particle->coord[boundaryAxis] - perBoundary->maxWall; 
			particle->coord[boundaryAxis] = perBoundary->minWall + difference;
			perBoundary->particlesUpdatedMaxEndCount++;
			Journal_DPrintfL( self->debug, 3, "moving to (%.4f,%.4f,%.4f).\n",
				particle->coord[I_AXIS], particle->coord[J_AXIS],
				particle->coord[K_AXIS] );
		}
		Stream_UnIndent( self->debug );
	}	

	Stream_UnIndent( self->debug );

	/* TODO: this is a bit of a hack to print this here using the lParticleI = swarm->total - 1, but its
	the only way I can see given this func is part of the SwarmAdvector intermediate. Should really be a 
	function on this class that updates all the particles. -- Main.PatrickSunter 15 May 2006 */
	if ( lParticle_I == (self->swarm->particleLocalCount-1) ) {
		PeriodicBoundary*      boundary = NULL;
		Index                  perB_I;
	
		Journal_DPrintfL( self->debug, 1, "PeriodicBoundariesManager total particles updated:\n" );
		Stream_Indent( self->debug );
		for ( perB_I = 0; perB_I < self->count; perB_I++ ) {
			boundary = &self->boundaries[perB_I];

			Journal_DPrintfL( self->debug, 1, "Periodic Boundary in %c Axis {%.2f,%.2f}: %d min end, %d max end\n",
				IJKTopology_DimNumToDimLetter[boundary->axis], boundary->minWall, boundary->maxWall,
				boundary->particlesUpdatedMinEndCount, boundary->particlesUpdatedMaxEndCount );
			/* Reset the counters for next time */
			boundary->particlesUpdatedMinEndCount = 0;	
			boundary->particlesUpdatedMaxEndCount = 0;	
		}
		Stream_UnIndent( self->debug );
	}
}	

