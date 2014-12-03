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
#include "CoincidentMapper.h"

#include "MaterialPointsSwarm.h"
#include "MaterialPoint.h"
#include "IntegrationPointsSwarm.h"

#include <assert.h>
#include <string.h>
#include <math.h>

const Type CoincidentMapper_Type = "CoincidentMapper";

CoincidentMapper* CoincidentMapper_New(
	Name							name,
	PICelleratorContext*		context,
	IntegrationPointsSwarm*	integrationSwarm,
	MaterialPointsSwarm*		materialSwarm )
{
  CoincidentMapper* self = (CoincidentMapper*)_CoincidentMapper_DefaultNew( name );

	self->isConstructed = True;
	_IntegrationPointMapper_Init( self, context, integrationSwarm );
	_OneToOneMapper_Init( self, materialSwarm );
	_CoincidentMapper_Init( self );

	return self;
}

void* _CoincidentMapper_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                                 _sizeOfSelf = sizeof(CoincidentMapper);
	Type                                                                         type = CoincidentMapper_Type;
	Stg_Class_DeleteFunction*                                                 _delete = _CoincidentMapper_Delete;
	Stg_Class_PrintFunction*                                                   _print = _CoincidentMapper_Print;
	Stg_Class_CopyFunction*                                                     _copy = _CoincidentMapper_Copy;
	Stg_Component_DefaultConstructorFunction*                     _defaultConstructor = _CoincidentMapper_DefaultNew;
	Stg_Component_ConstructFunction*                                       _construct = _CoincidentMapper_AssignFromXML;
	Stg_Component_BuildFunction*                                               _build = _CoincidentMapper_Build;
	Stg_Component_InitialiseFunction*                                     _initialise = _CoincidentMapper_Initialise;
	Stg_Component_ExecuteFunction*                                           _execute = _CoincidentMapper_Execute;
	Stg_Component_DestroyFunction*                                           _destroy = _CoincidentMapper_Destroy;
	AllocationType                                                 nameAllocationType = NON_GLOBAL;
	IntegrationPointMapper_MapFunction*                                          _map = _CoincidentMapper_Map;
	IntegrationPointMapper_GetMaterialPointsSwarmsFunction*  _getMaterialPointsSwarms = _OneToOneMapper_GetMaterialPointsSwarms;
	IntegrationPointMapper_GetMaterialIndexOnFunction*            _getMaterialIndexOn = _OneToOneMapper_GetMaterialIndexOn;
	IntegrationPointMapper_GetExtensionOnFunction*                    _getExtensionOn = _OneToOneMapper_GetExtensionOn;
    IntegrationPointMapper_GetDoubleFromExtension*                  _getDoubleFromExtension = _OneToOneMapper_GetDoubleFromExtension;
	IntegrationPointMapper_GetDoubleFromExtension*                  _getDoubleFromMaterial = _OneToOneMapper_GetDoubleFromMaterial;

	return _CoincidentMapper_New(  COINCIDENTMAPPER_PASSARGS  );
}

CoincidentMapper* _CoincidentMapper_New(  COINCIDENTMAPPER_DEFARGS  ) {
	CoincidentMapper* result;

	result = (CoincidentMapper*)_OneToOneMapper_New(  ONETOONEMAPPER_PASSARGS  );

	return result;
}

void _CoincidentMapper_AssignFromXML( void* mapper, Stg_ComponentFactory* cf, void* data ) {
	_OneToOneMapper_AssignFromXML( mapper, cf, data );
}

void _CoincidentMapper_Init( void* mapper ) {
}

void _CoincidentMapper_Delete( void* mapper ) {
	CoincidentMapper* self = (CoincidentMapper*)mapper;

	_OneToOneMapper_Delete( self );
}

void _CoincidentMapper_Print( void* mapper, Stream* stream ) {
	_IntegrationPointMapper_Print( mapper, stream );
}

void* _CoincidentMapper_Copy( const void* mapper, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	return _IntegrationPointMapper_Copy( mapper, dest, deep, nameExt, ptrMap );
}

void _CoincidentMapper_Build( void* mapper, void* cf ) {
	CoincidentMapper* self = (CoincidentMapper*)mapper;

	_OneToOneMapper_Build( self, cf );
}

void _CoincidentMapper_Initialise( void* mapper, void* cf ) {
	CoincidentMapper* self = (CoincidentMapper*)mapper;

	_OneToOneMapper_Initialise( self, cf );
}

void _CoincidentMapper_Execute( void* mapper, void* data ) {
}

void _CoincidentMapper_Destroy( void* mapper, void* data ) {
	CoincidentMapper* self = (CoincidentMapper*)mapper;

	_OneToOneMapper_Destroy( self, data );
}

void _CoincidentMapper_Map( void* mapper ) {
	CoincidentMapper*			self = (CoincidentMapper*)mapper;
	IntegrationPointsSwarm*	integrationSwarm = self->integrationSwarm;
	MaterialPointsSwarm*		materialSwarm = self->materialSwarm;
	IntegrationPoint*			integrationPoint;
	MaterialPoint*				materialPoint;
	MaterialPointRef*			ref;
	FeMesh*						mesh = materialSwarm->mesh;
	Particle_Index				particle_lI;
	Cell_Index					cell_dI;

#ifdef CAUTIOUS
	Index							dim_I;
	Stream*						errorStream = Journal_Register( Error_Type, self->type  );
#endif
	Stream*						debugStream = Swarm_Debug;
	
	Journal_DPrintfL( debugStream, 1, "In %s(): Re-creating a new set of integration points, exactly\n" 
		"\tmapping to the current material points & their positions.\n", __func__ );

	Stream_IndentBranch( debugStream );

	Journal_DPrintfL( debugStream, 2, "Reallocating the integration points array from size of %u points "
		"to the \n\tcurrent material swarm particle count of %u\n", integrationSwarm->particleLocalCount,
		materialSwarm->particleLocalCount );

	integrationSwarm->particleLocalCount = materialSwarm->particleLocalCount;
	Swarm_Realloc( integrationSwarm );

	Journal_DPrintfL( debugStream, 2, "Clearing all the cell->particle ownership tables, "
		"ready to add new cell->particle\n\trelationships as new integration points are set up.\n" );

	for( cell_dI = 0; cell_dI < integrationSwarm->cellDomainCount; cell_dI++ ) {
		integrationSwarm->cellParticleCountTbl[cell_dI] = 0;
		integrationSwarm->cellParticleSizeTbl[cell_dI] = 0;

		if ( integrationSwarm->cellParticleTbl[cell_dI] ) {
			Memory_Free( integrationSwarm->cellParticleTbl[cell_dI] );
		}
		integrationSwarm->cellParticleTbl[cell_dI] = NULL;
	}

	Journal_DPrintfL( debugStream, 2, "For each material particle, setting up a corresponding integration "
                      "point, and\n\tcalculating its element-local coord based on the material's global coord:\n" );
	Stream_IndentBranch( debugStream );	

	/* Map each point */
	for ( particle_lI = 0; particle_lI < materialSwarm->particleLocalCount; particle_lI++ ) {
		integrationPoint = (IntegrationPoint*)Swarm_ParticleAt( integrationSwarm, particle_lI );
		materialPoint = (MaterialPoint*)Swarm_ParticleAt( materialSwarm, particle_lI );

		cell_dI = materialPoint->owningCell;

		Journal_DPrintfL( debugStream, 3, "Referring to local material point %u, from material swarm cell %u:\n", particle_lI, cell_dI );
		Stream_IndentBranch( debugStream );	

		Journal_DPrintfL( debugStream, 3, "Adding new integration point %u to integration swarm cell %u\n", particle_lI, cell_dI );

		Swarm_AddParticleToCell( integrationSwarm, cell_dI, particle_lI );

		/* Convert global to local coordinates */
		ElementType_ConvertGlobalCoordToElLocal(
			FeMesh_GetElementType( mesh, cell_dI ),
			mesh, 
			cell_dI, 
			materialPoint->coord,
			integrationPoint->xi );

		Journal_DPrintfL( debugStream, 3, "Based on material point's coord of (%.2f,%.2f,%.2f):\n"
			"calculated and set new integration point's local coord as (%.2f,%.2f,%.2f)\n",
			materialPoint->coord[0], materialPoint->coord[1], materialPoint->coord[2],
			integrationPoint->xi[0], integrationPoint->xi[1], integrationPoint->xi[2] );

#ifdef CAUTIOUS
		/* Check the result is between -1 to 1 in all dimensions : if not, something is stuffed */		
		for ( dim_I= 0; dim_I < materialSwarm->dim; dim_I++ ) {
			Journal_Firewall(
				(integrationPoint->xi[dim_I] >= -1.001) && (integrationPoint->xi[dim_I] <= 1.001 ),
				errorStream,
				"Error - in %s(): unable to map material point %d in cell %d of swarm \"%s\" (type %s) "
				"coord to a valid \"local\" coordinate (xi). Coord was (%.3f,%.3f,%.3f), swarm's "
				"particle layout type was %s, Xi result was (%.4f,%.4f,%.4f).\n",
				__func__, particle_lI, cell_dI, materialSwarm->name, materialSwarm->type,
				materialPoint->coord[0], materialPoint->coord[1], materialPoint->coord[2],
				materialSwarm->particleLayout->type, integrationPoint->xi[0],
				integrationPoint->xi[1], integrationPoint->xi[2] );
		}
#endif

		ref = OneToOneMapper_GetMaterialRef( self, integrationPoint );
		ref->swarm_I = materialSwarm->swarmReg_I;
		ref->particle_I = particle_lI;
		Journal_DPrintfL( debugStream, 3, "updated the coincident mapper's material reference for "
			"this integration point to map back to the material point.\n" );
		Stream_UnIndentBranch( debugStream );	
	}

	Stream_UnIndentBranch( debugStream );	
	Journal_DPrintfL( debugStream, 2, "...finished updating local positions.\n" );

	Stream_UnIndentBranch( debugStream );
	Journal_DPrintfL( debugStream, 1, "...%s(): Done.\n", __func__ );
}


