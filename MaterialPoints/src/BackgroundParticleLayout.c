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
** This file may be distributed under the terms of the VPAC Public License
** as defined by VPAC of Australia and appearing in the file
** LICENSE.VPL included in the packaging of this file.
**
** This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
** WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
**
** $Id: BackgroundParticleLayout.c 3310 2005-10-26 07:10:18Z RobertTurnbull $
**
*/
/** \file
**  Role:
**	A particle layout which creates only 1 particle, for use in creating backgroud layers where materials are "everywhere"
**
** Assumptions:
**
** Comments:
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>

#include "types.h"
#include "BackgroundParticleLayout.h"

#include <assert.h>
#include <string.h>

const Type BackgroundParticleLayout_Type = "BackgroundParticleLayout";

BackgroundParticleLayout* _BackgroundParticleLayout_New( 
                SizeT                                                       _sizeOfSelf,
                Type                                                        type,
                Stg_Class_DeleteFunction*                                   _delete,
                Stg_Class_PrintFunction*                                    _print,
                Stg_Class_CopyFunction*                                     _copy,
                Stg_Component_DefaultConstructorFunction*                   _defaultConstructor,
                Stg_Component_ConstructFunction*                            _construct,
                Stg_Component_BuildFunction*                                _build,
                Stg_Component_InitialiseFunction*                           _initialise,
                Stg_Component_ExecuteFunction*                              _execute,
                Stg_Component_DestroyFunction*                              _destroy,
                ParticleLayout_SetInitialCountsFunction*                    _setInitialCounts,
                ParticleLayout_InitialiseParticlesFunction*                 _initialiseParticles,
                Name                                                        name,
                Bool                                                        initFlag,
		CoordSystem                                                 coordSystem,
                Bool                                                        weightsInitialisedAtStartup )
{
	BackgroundParticleLayout*		self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(BackgroundParticleLayout) );
	self = (BackgroundParticleLayout*)_ParticleLayout_New( 
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
			_setInitialCounts, 
			_initialiseParticles, 
			name, 
			initFlag,
			coordSystem,
			weightsInitialisedAtStartup );
	
	if( initFlag ){
		_BackgroundParticleLayout_Init( self, coordSystem, weightsInitialisedAtStartup );
	}
	
	return self;
}


void _BackgroundParticleLayout_Init(
		void*                  particleLayout,
		CoordSystem            coordSystem,
		Bool                   weightsInitialisedAtStartup )
{
	BackgroundParticleLayout* self = (BackgroundParticleLayout*)particleLayout;
	
	self->isConstructed = True;

	_ParticleLayout_Init( particleLayout, coordSystem, weightsInitialisedAtStartup );
}

void _BackgroundParticleLayout_Delete( void* particleLayout ) {
	BackgroundParticleLayout* self = (BackgroundParticleLayout*)particleLayout;
	
	_ParticleLayout_Delete( self );
}

void _BackgroundParticleLayout_Print( void* particleLayout, Stream* stream ) {
	BackgroundParticleLayout* self = (BackgroundParticleLayout*)particleLayout;
	
	Journal_Printf( stream, "BackgroundParticleLayout (ptr): %p\n", self );
	
	/* Parent class info */
	_ParticleLayout_Print( self, stream );
	
}


void* _BackgroundParticleLayout_Copy( void* particleLayout, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	BackgroundParticleLayout*		self = (BackgroundParticleLayout*)particleLayout;
	BackgroundParticleLayout*		newBackgroundParticleLayout;
	
	newBackgroundParticleLayout = (BackgroundParticleLayout*)_ParticleLayout_Copy( self, dest, deep, nameExt, ptrMap );

	return (void*)newBackgroundParticleLayout;
}

void* _BackgroundParticleLayout_DefaultNew( Name name ) {
	return _BackgroundParticleLayout_New(
			sizeof(BackgroundParticleLayout),
			BackgroundParticleLayout_Type,
			_BackgroundParticleLayout_Delete,
			_BackgroundParticleLayout_Print,
			_BackgroundParticleLayout_Copy,
			_BackgroundParticleLayout_DefaultNew,
			_BackgroundParticleLayout_Construct,
			_BackgroundParticleLayout_Build,
			_BackgroundParticleLayout_Initialise,
			_BackgroundParticleLayout_Execute,
			_BackgroundParticleLayout_Destroy,
			_BackgroundParticleLayout_SetInitialCounts,
			_BackgroundParticleLayout_InitialiseParticles,
			name,
			False,
			GlobalCoordSystem,
			False );
}
void  _BackgroundParticleLayout_Construct( void* component, Stg_ComponentFactory* cf, void* data )  {
	BackgroundParticleLayout*	self = (BackgroundParticleLayout*)component;

	self->context = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Context", DomainContext, False, data );
	if( !self->context )
		self->context = Stg_ComponentFactory_ConstructByName( cf, "context", DomainContext, True, data );

	_BackgroundParticleLayout_Init( component, GlobalCoordSystem, False );
}

void  _BackgroundParticleLayout_Build( void* component, void* data ) {}
void  _BackgroundParticleLayout_Initialise( void* component, void* data ) {}
void  _BackgroundParticleLayout_Execute( void* component, void* data ) {}
void  _BackgroundParticleLayout_Destroy( void* component, void* data ) {}


void _BackgroundParticleLayout_SetInitialCounts( void* particleLayout, void* _swarm )
{
	Swarm*			  swarm         = (Swarm*)_swarm;
	Cell_DomainIndex	  cell_I        = 0;
	char			  tempStr[100];

	swarm->particleLocalCount = 1;

	for( cell_I = 0; cell_I < swarm->cellLocalCount; cell_I++ ) {
		/* Set initial counts to empty, till we add the particles */
		swarm->cellParticleCountTbl[cell_I] = 0; 
		swarm->cellParticleSizeTbl[cell_I] = 1; /* Just to create array */

		sprintf( tempStr, "Swarm->cellParticleTbl[%d]", cell_I );
		swarm->cellParticleTbl[cell_I] = Memory_Alloc_Array( Particle_Index, swarm->cellParticleSizeTbl[cell_I], tempStr );
										
	}

	/* Now initialise the shadow cell particle counts */
	for (; cell_I < swarm->cellDomainCount; cell_I++ ) {
		swarm->cellParticleCountTbl[cell_I] = 0;
		swarm->cellParticleSizeTbl[cell_I] = 0;
		swarm->cellParticleTbl[cell_I] = NULL;
	}
}

void _BackgroundParticleLayout_InitialiseParticles( void* particleLayout, void* _swarm )
{
	Swarm*                    swarm         = (Swarm*)_swarm;
	GlobalParticle*           particle;

	particle = (GlobalParticle*)Swarm_ParticleAt( swarm, 0 );

	particle->owningCell = 0;
	particle->coord[0] = 0.5; /* just some value, doesn't matter where */
	particle->coord[1] = 0.5;
	particle->coord[2] = 0.5;

	Swarm_AddParticleToCell( swarm, 0, 0 );
}

