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
** $Id: MappedParticleLayout.c 3310 2005-10-26 07:10:18Z RobertTurnbull $
**
*/
/** \file
 **  Role:
 **      A particle layout designed for IntegrationPointsSwarms where particle positions and weights are mapped by
 **	another matieral points swarm and a weights calculator. Hence this particle layout does nothing except for
 **	creating the swarm's cell table arrays
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

#include "MaterialPoints.h"

#include <assert.h>
#include <string.h>

const Type MappedParticleLayout_Type = "MappedParticleLayout";

MappedParticleLayout* MappedParticleLayout_New( 
	Name					name, 
	AbstractContext*	context,
	CoordSystem			coordSystem,
	Bool					weightsInitialisedAtStartup )
{

   MappedParticleLayout* self = _MappedParticleLayout_DefaultNew( name ); 

   _ParticleLayout_Init( self, context, coordSystem, weightsInitialisedAtStartup );
   _MappedParticleLayout_Init( self );

	return self;
}
MappedParticleLayout* _MappedParticleLayout_New(  MAPPEDPARTICLELAYOUT_DEFARGS  )
{
   MappedParticleLayout*		self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof(MappedParticleLayout) );
   self = (MappedParticleLayout*)_ParticleLayout_New(  PARTICLELAYOUT_PASSARGS  );

   return self;
}


void _MappedParticleLayout_Init( void* particleLayout ) {
	MappedParticleLayout* self;

	self = (MappedParticleLayout*)particleLayout;
}

void _MappedParticleLayout_Delete( void* particleLayout ) {
    MappedParticleLayout* self = (MappedParticleLayout*)particleLayout;
	
    _ParticleLayout_Delete( self );
}

void _MappedParticleLayout_Print( void* particleLayout, Stream* stream ) {
    MappedParticleLayout* self = (MappedParticleLayout*)particleLayout;
	
    Journal_Printf( stream, "MappedParticleLayout (ptr): %p\n", self );
	
    /* Parent class info */
    _ParticleLayout_Print( self, stream );
	
}

void* _MappedParticleLayout_Copy( void* particleLayout, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
    MappedParticleLayout*		self = (MappedParticleLayout*)particleLayout;
    MappedParticleLayout*		newMappedParticleLayout;
	
    newMappedParticleLayout = (MappedParticleLayout*)_ParticleLayout_Copy( self, dest, deep, nameExt, ptrMap );

    return (void*)newMappedParticleLayout;
}

void* _MappedParticleLayout_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                        _sizeOfSelf = sizeof(MappedParticleLayout);
	Type                                                                type = MappedParticleLayout_Type;
	Stg_Class_DeleteFunction*                                        _delete = _MappedParticleLayout_Delete;
	Stg_Class_PrintFunction*                                          _print = _MappedParticleLayout_Print;
	Stg_Class_CopyFunction*                                            _copy = _MappedParticleLayout_Copy;
	Stg_Component_DefaultConstructorFunction*            _defaultConstructor = _MappedParticleLayout_DefaultNew;
	Stg_Component_ConstructFunction*                              _construct = _MappedParticleLayout_AssignFromXML;
	Stg_Component_BuildFunction*                                      _build = _MappedParticleLayout_Build;
	Stg_Component_InitialiseFunction*                            _initialise = _MappedParticleLayout_Initialise;
	Stg_Component_ExecuteFunction*                                  _execute = _MappedParticleLayout_Execute;
	Stg_Component_DestroyFunction*                                  _destroy = _MappedParticleLayout_Destroy;
	AllocationType                                        nameAllocationType = NON_GLOBAL;
	ParticleLayout_SetInitialCountsFunction*               _setInitialCounts = _MappedParticleLayout_SetInitialCounts;
	ParticleLayout_InitialiseParticlesFunction*         _initialiseParticles = _MappedParticleLayout_InitialiseParticles;
	CoordSystem                                                  coordSystem = LocalCoordSystem;
	Bool                                         weightsInitialisedAtStartup = False;

    return _MappedParticleLayout_New(  MAPPEDPARTICLELAYOUT_PASSARGS  );
}
void  _MappedParticleLayout_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
    MappedParticleLayout*		self = (MappedParticleLayout*)component;

    _ParticleLayout_AssignFromXML( self, cf, data );

    _MappedParticleLayout_Init( self );
}
void  _MappedParticleLayout_Build( void* component, void* data ) {}
void  _MappedParticleLayout_Initialise( void* component, void* data ) {}
void  _MappedParticleLayout_Execute( void* component, void* data ) {}
void  _MappedParticleLayout_Destroy( void* component, void* data ) {
    MappedParticleLayout*		self = (MappedParticleLayout*)component;

    _ParticleLayout_Destroy( self, data );
}

void _MappedParticleLayout_SetInitialCounts( void* particleLayout, void* _swarm )
{
    Swarm*			swarm = (Swarm*)_swarm;
    Cell_DomainIndex	cell_I = 0;
    char			tempStr[100];

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

void _MappedParticleLayout_InitialiseParticles( void* particleLayout, void* _swarm )
{
    /* Don't need to do anything */
}



