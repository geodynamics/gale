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
** $Id: TriGaussParticleLayout.c 3851 2006-10-12 08:57:22Z SteveQuenette $
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
#include "PerCellParticleLayout.h"
#include "TriGaussParticleLayout.h"

#include "SwarmClass.h"
#include "StandardParticle.h"
#include "IntegrationPoint.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>


const Type TriGaussParticleLayout_Type = "TriGaussParticleLayout";

TriGaussParticleLayout* TriGaussParticleLayout_New( 
   Name name, 
   AbstractContext* context,
   CoordSystem      coordSystem,
   Bool             weightsInitialisedAtStartup,
   unsigned int dim, unsigned int particlesPerCell ) 
{
	TriGaussParticleLayout* self = (TriGaussParticleLayout*)_TriGaussParticleLayout_DefaultNew( name );

   _ParticleLayout_Init( self, context, coordSystem, weightsInitialisedAtStartup );
   _PerCellParticleLayout_Init( self );
	_TriGaussParticleLayout_Init( self, dim, particlesPerCell );

	return self;
}

TriGaussParticleLayout* _TriGaussParticleLayout_New(  TRIGAUSSPARTICLELAYOUT_DEFARGS  )
{
	TriGaussParticleLayout* self;
	
   /* hard-wire here */
   coordSystem = LocalCoordSystem;
   weightsInitialisedAtStartup = True;
	/* Allocate memory */
	self = (TriGaussParticleLayout*)_PerCellParticleLayout_New(  PERCELLPARTICLELAYOUT_PASSARGS  );
	
	return self;
}


void _TriGaussParticleLayout_Init(
		TriGaussParticleLayout* self,
		unsigned int            dim,
		unsigned int            particlesPerCell )
{
	self->isConstructed    = True;
	self->dim              = dim;
	self->particlesPerCell = particlesPerCell;
}

void _TriGaussParticleLayout_Delete( void* triGaussParticleLayout ) {
	TriGaussParticleLayout* self = (TriGaussParticleLayout*)triGaussParticleLayout;
	
	_PerCellParticleLayout_Delete( self );
}

void _TriGaussParticleLayout_Print( void* triGaussParticleLayout, Stream* stream ) {
	TriGaussParticleLayout* self = (TriGaussParticleLayout*)triGaussParticleLayout;
	
	/* Set the Journal for printing informations */
	Stream* triGaussParticleLayoutStream = stream;
	
	/* General info */
	Journal_Printf( triGaussParticleLayoutStream, "TriGaussParticleLayout (ptr): %p:\n", self );
	
	/* Parent class info */
	_PerCellParticleLayout_Print( self, stream );
	
	/* Virtual info */
	
	/* TriGaussParticleLayout */
	Journal_Printf( triGaussParticleLayoutStream, "\tdim: %u\n", self->dim );
	Journal_Printf( triGaussParticleLayoutStream, "\tparticlesPerCell: %u\n", self->particlesPerCell  );
}


void* _TriGaussParticleLayout_Copy( void* triGaussParticleLayout, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	TriGaussParticleLayout*		self = (TriGaussParticleLayout*)triGaussParticleLayout;
	TriGaussParticleLayout*		newTriGaussParticleLayout;
	
	newTriGaussParticleLayout = (TriGaussParticleLayout*)_PerCellParticleLayout_Copy( self, dest, deep, nameExt, ptrMap );
	
	return (void*)newTriGaussParticleLayout;
}

void* _TriGaussParticleLayout_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                                     _sizeOfSelf = sizeof(TriGaussParticleLayout);
	Type                                                                             type = TriGaussParticleLayout_Type;
	Stg_Class_DeleteFunction*                                                     _delete = _TriGaussParticleLayout_Delete;
	Stg_Class_PrintFunction*                                                       _print = _TriGaussParticleLayout_Print;
	Stg_Class_CopyFunction*                                                         _copy = _TriGaussParticleLayout_Copy;
	Stg_Component_DefaultConstructorFunction*                         _defaultConstructor = _TriGaussParticleLayout_DefaultNew;
	Stg_Component_ConstructFunction*                                           _construct = _TriGaussParticleLayout_AssignFromXML;
	Stg_Component_BuildFunction*                                                   _build = _TriGaussParticleLayout_Build;
	Stg_Component_InitialiseFunction*                                         _initialise = _TriGaussParticleLayout_Initialise;
	Stg_Component_ExecuteFunction*                                               _execute = _TriGaussParticleLayout_Execute;
	Stg_Component_DestroyFunction*                                               _destroy = _TriGaussParticleLayout_Destroy;
	AllocationType                                                     nameAllocationType = NON_GLOBAL;
	ParticleLayout_SetInitialCountsFunction*                            _setInitialCounts = _PerCellParticleLayout_SetInitialCounts;
	ParticleLayout_InitialiseParticlesFunction*                      _initialiseParticles = _PerCellParticleLayout_InitialiseParticles;
	CoordSystem                                                               coordSystem = LocalCoordSystem;
	Bool                                                      weightsInitialisedAtStartup = True;
	PerCellParticleLayout_InitialCountFunction*                             _initialCount = _TriGaussParticleLayout_InitialCount;
	PerCellParticleLayout_InitialiseParticlesOfCellFunction*   _initialiseParticlesOfCell = _TriGaussParticleLayout_InitialiseParticlesOfCell;

	return _TriGaussParticleLayout_New(  TRIGAUSSPARTICLELAYOUT_PASSARGS  );
}

void _TriGaussParticleLayout_AssignFromXML( void* triGaussParticleLayout, Stg_ComponentFactory* cf, void* data ){
	TriGaussParticleLayout *self = (TriGaussParticleLayout*)triGaussParticleLayout;
	unsigned int dim;
	unsigned int particlesPerCell;

   _PerCellParticleLayout_AssignFromXML( self, cf, data );
	dim = Dictionary_Entry_Value_AsUnsignedInt(
		Dictionary_GetDefault( cf->rootDict, "dim", Dictionary_Entry_Value_FromUnsignedInt( 3 ) ) );

	particlesPerCell = Dictionary_Entry_Value_AsUnsignedInt(
		Dictionary_GetDefault( cf->rootDict, "particlesPerCell", Dictionary_Entry_Value_FromUnsignedInt( 1 ) ) );

	_TriGaussParticleLayout_Init( self, dim, particlesPerCell );
}
	
void _TriGaussParticleLayout_Build( void* triGaussParticleLayout, void* data ) {
}
	
void _TriGaussParticleLayout_Initialise( void* triGaussParticleLayout, void* data ) {
}
	
void _TriGaussParticleLayout_Execute( void* triGaussParticleLayout, void* data ) {
}
	
void _TriGaussParticleLayout_Destroy( void* triGaussParticleLayout, void* data ) {
   TriGaussParticleLayout* self = (TriGaussParticleLayout*)triGaussParticleLayout;
   _PerCellParticleLayout_Destroy( self, data );
}

Particle_InCellIndex _TriGaussParticleLayout_InitialCount( void* triGaussParticleLayout, void* celllayout, Cell_Index cell_I )
{
	TriGaussParticleLayout* self = (TriGaussParticleLayout*)triGaussParticleLayout;
	Particle_InCellIndex count;
	int dim;	
	
	dim = self->dim;
	count = (Particle_InCellIndex)( self->particlesPerCell );
	
	return count;
}


/* remember this only has to initialise one cell of particles at a time */
void _TriGaussParticleLayout_InitialiseParticlesOfCell( void* triGaussParticleLayout, void* _swarm, Cell_Index cell_I )
{
	#define MAX_DIM 3
	#define MAX_GAUSS_POINTS_2D 1
	#define MAX_GAUSS_POINTS_3D 1
	TriGaussParticleLayout*   self = (TriGaussParticleLayout*)triGaussParticleLayout;
	Swarm*                    swarm = (Swarm*)_swarm;
	IntegrationPoint*         integrationPoint = NULL;

	Particle_InCellIndex      cParticle_I = 0;

	Particle_InCellIndex ppc;
	int dim;
	static double weight[20];
	static double xi[20][3];
	static int beenHere = 0;
	int d;
	
	dim = self->dim;
	ppc = self->particlesPerCell;
	
	if( dim == 2 ) {
		if( ppc == 1 ) {
			weight[0] = 0.5;

			xi[0][0] = 0.333333333333;
			xi[0][1] = 0.333333333333;				
		}
		if( ppc > MAX_GAUSS_POINTS_2D ) {
			Journal_Firewall(
				ppc > MAX_GAUSS_POINTS_2D,
				Journal_MyStream( Error_Type, self ),
				"In %s(), error: particlesPerCell greater than implementated tabulated gauss values of %d\n",
				__func__,
				MAX_GAUSS_POINTS_2D );
		}
	}
	else if ( dim == 3 ) {
		if( ppc == 1 ) {
			weight[0] = 0.5;			

			xi[0][0] = 0.333333333333;
			xi[0][1] = 0.333333333333;				
			xi[0][2] = 0.333333333333;				
		}
		if( ppc > MAX_GAUSS_POINTS_3D ) {
			Journal_Firewall(
				ppc > MAX_GAUSS_POINTS_3D,
				Journal_MyStream( Error_Type, self ),
				"In %s(), error: particlesPerCell greater than implementated tabulated gauss values of %d\n",
				__func__,
				MAX_GAUSS_POINTS_3D );
		}
	}

	assert( ppc <= swarm->cellParticleCountTbl[cell_I] );
	
	for ( cParticle_I = 0; cParticle_I < ppc; cParticle_I++ ) {
		integrationPoint = (IntegrationPoint*)Swarm_ParticleInCellAt( swarm, cell_I, cParticle_I );
		integrationPoint->owningCell = cell_I;
			
		for( d = 0; d < dim; d++ ) {
			/*integrationPoint->coord[d] = xi[cParticle_I][d];*/
			integrationPoint->xi[d] = xi[cParticle_I][d];
		}
		
		integrationPoint->weight = weight[cParticle_I];
		
		beenHere = 1;
	}	
	
	
}


