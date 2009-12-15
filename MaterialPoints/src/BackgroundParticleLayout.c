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

BackgroundParticleLayout* BackgroundParticleLayout_New( Name name,
   AbstractContext* context, 
   CoordSystem      coordSystem, 
   Bool             weightsInitialisedAtStartup ) {

   BackgroundParticleLayout* self = _BackgroundParticleLayout_DefaultNew( name );

   _ParticleLayout_Init( self, context, coordSystem, weightsInitialisedAtStartup );
   _BackgroundParticleLayout_Init( self );
   self->isConstructed = True;
   return self;
}
BackgroundParticleLayout* _BackgroundParticleLayout_New(  BACKGROUNDPARTICLELAYOUT_DEFARGS  )
{
    BackgroundParticleLayout*		self;
	
    /* Allocate memory */
    assert( _sizeOfSelf >= sizeof(BackgroundParticleLayout) );
    self = (BackgroundParticleLayout*)_ParticleLayout_New(  PARTICLELAYOUT_PASSARGS  );
	
    return self;
}


void _BackgroundParticleLayout_Init(
    void*                  particleLayout )
{
    BackgroundParticleLayout* self = (BackgroundParticleLayout*)particleLayout;
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
	/* Variables set in this function */
	SizeT                                                        _sizeOfSelf = sizeof(BackgroundParticleLayout);
	Type                                                                type = BackgroundParticleLayout_Type;
	Stg_Class_DeleteFunction*                                        _delete = _BackgroundParticleLayout_Delete;
	Stg_Class_PrintFunction*                                          _print = _BackgroundParticleLayout_Print;
	Stg_Class_CopyFunction*                                            _copy = _BackgroundParticleLayout_Copy;
	Stg_Component_DefaultConstructorFunction*            _defaultConstructor = _BackgroundParticleLayout_DefaultNew;
	Stg_Component_ConstructFunction*                              _construct = _BackgroundParticleLayout_AssignFromXML;
	Stg_Component_BuildFunction*                                      _build = _BackgroundParticleLayout_Build;
	Stg_Component_InitialiseFunction*                            _initialise = _BackgroundParticleLayout_Initialise;
	Stg_Component_ExecuteFunction*                                  _execute = _BackgroundParticleLayout_Execute;
	Stg_Component_DestroyFunction*                                  _destroy = _BackgroundParticleLayout_Destroy;
	AllocationType                                        nameAllocationType = NON_GLOBAL;
	ParticleLayout_SetInitialCountsFunction*               _setInitialCounts = _BackgroundParticleLayout_SetInitialCounts;
	ParticleLayout_InitialiseParticlesFunction*         _initialiseParticles = _BackgroundParticleLayout_InitialiseParticles;
	CoordSystem                                                  coordSystem = GlobalCoordSystem;
	Bool                                         weightsInitialisedAtStartup = False;

    return _BackgroundParticleLayout_New(  BACKGROUNDPARTICLELAYOUT_PASSARGS  );
}
void  _BackgroundParticleLayout_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data )  {
   BackgroundParticleLayout*	self = (BackgroundParticleLayout*)component;

   _ParticleLayout_AssignFromXML( self, cf, data );

   _BackgroundParticleLayout_Init( self );
}

void  _BackgroundParticleLayout_Build( void* component, void* data ) {}
void  _BackgroundParticleLayout_Initialise( void* component, void* data ) {}
void  _BackgroundParticleLayout_Execute( void* component, void* data ) {}
void  _BackgroundParticleLayout_Destroy( void* component, void* data ) {
   BackgroundParticleLayout*	self = (BackgroundParticleLayout*)component;
   _ParticleLayout_Destroy( self, data );
}


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



