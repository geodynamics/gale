/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org) ) {
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
** $Id: GaussMapper.c 189 2005-10-20 00:39:29Z RobertTurnbull $
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
#include "IntegrationPointMapper.h"
#include "OneToOneMapper.h"
#include "GaussMapper.h"
#include "IntegrationPointsSwarm.h"
#include "BackgroundParticleLayout.h"
#include "MaterialPoint.h"
#include "MaterialPointsSwarm.h"

#include <assert.h>
#include <string.h>
#include <math.h>

const Type GaussMapper_Type = "GaussMapper";

GaussMapper* _GaussMapper_New(
		SizeT                                                           _sizeOfSelf,
		Type                                                            type,
		Stg_Class_DeleteFunction*                                       _delete,
		Stg_Class_PrintFunction*                                        _print,
		Stg_Class_CopyFunction*                                         _copy,
		Stg_Component_DefaultConstructorFunction*                       _defaultConstructor,
		Stg_Component_ConstructFunction*                                _construct,
		Stg_Component_BuildFunction*                                    _build,
		Stg_Component_InitialiseFunction*                               _initialise,
		Stg_Component_ExecuteFunction*                                  _execute,
		Stg_Component_DestroyFunction*                                  _destroy,
		IntegrationPointMapper_MapFunction*                             _map,
		IntegrationPointMapper_GetMaterialPointsSwarmsFunction*         _getMaterialPointsSwarms,
		IntegrationPointMapper_GetMaterialIndexOnFunction*              _getMaterialIndexOn,
		IntegrationPointMapper_GetExtensionOnFunction*                  _getExtensionOn,
		Name                                                            name,
		Bool                                                            initFlag,

		IntegrationPointsSwarm*                                         integrationSwarm,
		MaterialPointsSwarm*                                            materialSwarm )
{
	GaussMapper* result;

	result = (GaussMapper*)_OneToOneMapper_New(
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
			_map,
			_getMaterialPointsSwarms,
			_getMaterialIndexOn,
			_getExtensionOn,
			name,
			initFlag,
			integrationSwarm,
			materialSwarm );

	if (initFlag) {
		_GaussMapper_Init( result, integrationSwarm, materialSwarm );
	}
	return result;
}

void _GaussMapper_Init(
		void*                   mapper,
		IntegrationPointsSwarm* integrationSwarm,
		MaterialPointsSwarm*    materialSwarm ) 
{
	_OneToOneMapper_Init( mapper, integrationSwarm, materialSwarm );
}

void _GaussMapper_Delete( void* mapper ) {
	_IntegrationPointMapper_Delete( mapper );
}
void _GaussMapper_Print( void* mapper, Stream* stream ) {
	_IntegrationPointMapper_Print( mapper, stream );
}
void* _GaussMapper_Copy( void* mapper, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	return _IntegrationPointMapper_Copy( mapper, dest, deep, nameExt, ptrMap );
}

void* _GaussMapper_DefaultNew( Name name ) {
	return _GaussMapper_New(
			sizeof(GaussMapper),
			GaussMapper_Type,
			_GaussMapper_Delete,
			_GaussMapper_Print,
			_GaussMapper_Copy,
			_GaussMapper_DefaultNew,
			_GaussMapper_Construct,
			_GaussMapper_Build,
			_GaussMapper_Initialise,
			_GaussMapper_Execute,
			_GaussMapper_Destroy,
			_GaussMapper_Map,
			_OneToOneMapper_GetMaterialPointsSwarms,
			_OneToOneMapper_GetMaterialIndexOn,
			_OneToOneMapper_GetExtensionOn,
			name,
			False,
			NULL,
			NULL );
}

void _GaussMapper_Construct( void* mapper, Stg_ComponentFactory* cf, void* data ) {
	GaussMapper* self = (GaussMapper*)mapper;
	
	_OneToOneMapper_Construct( mapper, cf, data );

	/* Validate assumptions on the layouts the swarms use */
	
	Journal_Firewall(
			Stg_Class_IsInstance( self->integrationSwarm->particleLayout, GaussParticleLayout_Type ) ||
			Stg_Class_IsInstance( self->integrationSwarm->particleLayout, TriGaussParticleLayout_Type ),
			Journal_MyStream( Error_Type, self ),
			"In func %s, To use Gauss mapper, the integration point swarm %s must use a %s\n",
			__func__,
			self->integrationSwarm->name,
			GaussParticleLayout_Type );
	Journal_Firewall(
			Stg_Class_IsInstance( self->materialSwarm->particleLayout, BackgroundParticleLayout_Type ),
			Journal_MyStream( Error_Type, self ),
			"In func %s, To use Gauss mapper, the material point swarm %s must use a %s\n",
			__func__,
			self->materialSwarm->particleLayout,
			BackgroundParticleLayout_Type );
}

void _GaussMapper_Build( void* mapper, void* cf ) {
	_OneToOneMapper_Build( mapper, cf );
}
void _GaussMapper_Initialise( void* mapper, void* cf ) {
	_OneToOneMapper_Initialise( mapper, cf );
}
void _GaussMapper_Execute( void* mapper, void* data ) {
}
void _GaussMapper_Destroy( void* mapper, void* data ) {
}

void _GaussMapper_Map( void* mapper ) {
	GaussMapper*            self                   = (GaussMapper*)mapper;

	IntegrationPointsSwarm* integrationSwarm       = self->integrationSwarm;
	MaterialPointsSwarm*    materialSwarm          = self->materialSwarm;

	IntegrationPoint*       integrationPoint;
	MaterialPoint*          materialPoint;
	MaterialPointRef*       ref;

	Particle_Index          point_I;

	materialPoint = (MaterialPoint*)Swarm_ParticleAt( materialSwarm, 0 ); /* Get the first and only point */
	
	/* Map each point in integration to the single material point for its properties */
	for ( point_I = 0; point_I < integrationSwarm->particleLocalCount; point_I++ ) {
		integrationPoint = (IntegrationPoint*)Swarm_ParticleAt( integrationSwarm, point_I );

		ref = OneToOneMapper_GetMaterialRef( self, integrationPoint );
		ref->swarm_I = materialSwarm->swarmReg_I;
		ref->particle_I = 0;
	}
}
