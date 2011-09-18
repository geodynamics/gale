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
** $Id:  $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>

#include "types.h"
#include "IntegrationPointMapper.h"
#include "OneToOneMapper.h"
#include "GaussCoincidentMapper.h"

#include "MaterialPointsSwarm.h"
#include "MaterialPoint.h"
#include "IntegrationPointsSwarm.h"

#include <assert.h>
#include <string.h>
#include <math.h>

const Type GaussCoincidentMapper_Type = "GaussCoincidentMapper";

GaussCoincidentMapper* GaussCoincidentMapper_New(
	Name								name,
	PICelleratorContext*			context,
	IntegrationPointsSwarm*		integrationSwarm,
	MaterialPointsSwarm*			materialSwarm )
{
  GaussCoincidentMapper* self = (GaussCoincidentMapper*)_GaussCoincidentMapper_DefaultNew( name );

	self->isConstructed = True;
	_IntegrationPointMapper_Init( self, context, integrationSwarm );
   _OneToOneMapper_Init( self, materialSwarm );
	_GaussCoincidentMapper_Init( self );

	return self;
}

void* _GaussCoincidentMapper_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                                 _sizeOfSelf = sizeof(GaussCoincidentMapper);
	Type                                                                         type = GaussCoincidentMapper_Type;
	Stg_Class_DeleteFunction*                                                 _delete = _GaussCoincidentMapper_Delete;
	Stg_Class_PrintFunction*                                                   _print = _GaussCoincidentMapper_Print;
	Stg_Class_CopyFunction*                                                     _copy = _GaussCoincidentMapper_Copy;
	Stg_Component_DefaultConstructorFunction*                     _defaultConstructor = _GaussCoincidentMapper_DefaultNew;
	Stg_Component_ConstructFunction*                                       _construct = _GaussCoincidentMapper_AssignFromXML;
	Stg_Component_BuildFunction*                                               _build = _GaussCoincidentMapper_Build;
	Stg_Component_InitialiseFunction*                                     _initialise = _GaussCoincidentMapper_Initialise;
	Stg_Component_ExecuteFunction*                                           _execute = _GaussCoincidentMapper_Execute;
	Stg_Component_DestroyFunction*                                           _destroy = _GaussCoincidentMapper_Destroy;
	AllocationType                                                 nameAllocationType = NON_GLOBAL;
	IntegrationPointMapper_MapFunction*                                          _map = _GaussCoincidentMapper_Map;
	IntegrationPointMapper_GetMaterialPointsSwarmsFunction*  _getMaterialPointsSwarms = _OneToOneMapper_GetMaterialPointsSwarms;
	IntegrationPointMapper_GetMaterialIndexOnFunction*            _getMaterialIndexOn = _OneToOneMapper_GetMaterialIndexOn;
	IntegrationPointMapper_GetExtensionOnFunction*                    _getExtensionOn = _OneToOneMapper_GetExtensionOn;
        IntegrationPointMapper_GetDoubleFromExtension*                  _getDoubleFromExtension = _OneToOneMapper_GetDoubleFromExtension;
    IntegrationPointMapper_GetDoubleFromMaterial*                  _getDoubleFromMaterial = _OneToOneMapper_GetDoubleFromMaterial;

	return _GaussCoincidentMapper_New(  GAUSSCOINCIDENTMAPPER_PASSARGS  );
}

GaussCoincidentMapper* _GaussCoincidentMapper_New(  GAUSSCOINCIDENTMAPPER_DEFARGS  ) {
	GaussCoincidentMapper* result;

	result = (GaussCoincidentMapper*)_OneToOneMapper_New(  ONETOONEMAPPER_PASSARGS  );

	return result;
}

void _GaussCoincidentMapper_AssignFromXML( void* mapper, Stg_ComponentFactory* cf, void* data ) {
	GaussCoincidentMapper* self = (GaussCoincidentMapper*)mapper;

	_OneToOneMapper_AssignFromXML( self, cf, data );
}

void _GaussCoincidentMapper_Init( void* mapper ) {
}

void _GaussCoincidentMapper_Delete( void* mapper ) {
	GaussCoincidentMapper* self = (GaussCoincidentMapper*)mapper;

	_OneToOneMapper_Delete( self );
}

void _GaussCoincidentMapper_Print( void* mapper, Stream* stream ) {
	GaussCoincidentMapper* self = (GaussCoincidentMapper*)mapper;

	_IntegrationPointMapper_Print( self, stream );
}

void* _GaussCoincidentMapper_Copy( const void* mapper, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	return _IntegrationPointMapper_Copy( mapper, dest, deep, nameExt, ptrMap );
}

void _GaussCoincidentMapper_Build( void* mapper, void* cf ) {
	GaussCoincidentMapper* self = (GaussCoincidentMapper*)mapper;

	_OneToOneMapper_Build( self, cf );
}

void _GaussCoincidentMapper_Initialise( void* mapper, void* cf ) {
	GaussCoincidentMapper* self = (GaussCoincidentMapper*)mapper;

	_OneToOneMapper_Initialise( self, cf );
}

void _GaussCoincidentMapper_Execute( void* mapper, void* data ) {
}

void _GaussCoincidentMapper_Destroy( void* mapper, void* data ) {
	GaussCoincidentMapper* self = (GaussCoincidentMapper*)mapper;

	_OneToOneMapper_Destroy( self, data );
}

void _GaussCoincidentMapper_Map( void* mapper ) {
	GaussCoincidentMapper*	self = (GaussCoincidentMapper*)mapper;
	IntegrationPointsSwarm*	integrationSwarm = self->integrationSwarm;
	MaterialPointsSwarm*    materialSwarm = self->materialSwarm;
	IntegrationPoint*			integrationPoint;
	MaterialPointRef*			ref;
	Particle_Index				particle_lI;
#if 0
	MaterialPoint*				materialPoint;
	FeMesh*						mesh = materialSwarm->mesh;
	Cell_Index					cell_dI;
#endif

#ifdef CAUTIOUS
    Index						dim_I;
    Stream*						errorStream = Journal_Register( Error_Type, self->type  );
#endif
    Stream*						debugStream = Swarm_Debug;
	
#if 0
	Journal_DPrintfL( debugStream, 1, "In %s(): Re-creating a new set of integration points, exactly\n" 
		"\tmapping to the current material points & their positions.\n", __func__ ) ;
	Stream_IndentBranch( debugStream );

	Journal_DPrintfL( debugStream, 2, "Reallocating the integration points array from size of %u points "
		"to the \n\tcurrent material swarm particle count of %u\n", integrationSwarm->particleLocalCount,
		materialSwarm->particleLocalCount );
	materialSwarm->particleLocalCount = integrationSwarm->particleLocalCount;
	Swarm_Realloc( materialSwarm );

	Journal_DPrintfL( debugStream, 2, "Clearing all the cell->particle ownership tables, "
		"ready to add new cell->particle\n\trelationships as new integration points are set up.\n" );

	for( cell_dI = 0; cell_dI < materialSwarm->cellDomainCount; cell_dI++ ) {
		materialSwarm->cellParticleCountTbl[cell_dI] = 0;
		materialSwarm->cellParticleSizeTbl[cell_dI] = 0;

		if ( materialSwarm->cellParticleTbl[cell_dI] ) {
			Memory_Free( materialSwarm->cellParticleTbl[cell_dI] );
		}
		materialSwarm->cellParticleTbl[cell_dI] = NULL;
	}

	Journal_DPrintfL( debugStream, 2, "For each material particle, setting up a corresponding integration "
		"point, and\n\tcalculating its element-local coord based on the material's global coord:\n" );
		Stream_IndentBranch( debugStream );
#endif

	/* Map each point */
	for ( particle_lI = 0; particle_lI < integrationSwarm->particleLocalCount; particle_lI++ ) {
		integrationPoint = (IntegrationPoint*)Swarm_ParticleAt( integrationSwarm, particle_lI );

#if 0
		materialPoint = (MaterialPoint*)Swarm_ParticleAt( materialSwarm, particle_lI );

		cell_dI = integrationPoint->owningCell;

		Journal_DPrintfL( debugStream, 3, "Referring to local material point %u, from material swarm cell %u:\n", particle_lI, cell_dI );
		Stream_IndentBranch( debugStream );

		Journal_DPrintfL( debugStream, 3, "Adding new integration point %u to integration swarm cell %u\n", particle_lI, cell_dI );

		Swarm_AddParticleToCell( materialSwarm, cell_dI, particle_lI );

		/* Convert local to global coordinates */
		FeMesh_CoordLocalToGlobal(mesh, cell_dI, integrationPoint->xi, materialPoint->coord);

		Journal_DPrintfL( debugStream, 3, "Based on material point's coord of (%.2f,%.2f,%.2f):\n"
			"calculated and set new integration point's local coord as (%.2f,%.2f,%.2f)\n",
			materialPoint->coord[0], materialPoint->coord[1], materialPoint->coord[2],
			integrationPoint->xi[0], integrationPoint->xi[1], integrationPoint->xi[2] );
#endif

		ref = OneToOneMapper_GetMaterialRef( self, integrationPoint );
		ref->swarm_I = materialSwarm->swarmReg_I;
		ref->particle_I = particle_lI;
		Journal_DPrintfL( debugStream, 3, "updated the coincident mapper's material reference for "
			"this integration point to map back to the material point.\n" );
		Stream_UnIndentBranch( debugStream );	
	}

#if 0
	Stream_UnIndentBranch( debugStream );	
	Journal_DPrintfL( debugStream, 2, "...finished updating local positions.\n" );

	Stream_UnIndentBranch( debugStream );
	Journal_DPrintfL( debugStream, 1, "...%s(): Done.\n", __func__ );
#endif

}


