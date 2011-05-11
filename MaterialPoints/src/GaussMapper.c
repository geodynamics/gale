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

GaussMapper* _GaussMapper_New(  GAUSSMAPPER_DEFARGS  ) {
	GaussMapper* result;

	result = (GaussMapper*)_OneToOneMapper_New(  ONETOONEMAPPER_PASSARGS  );

	return result;
}

void* _GaussMapper_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                                 _sizeOfSelf = sizeof(GaussMapper);
	Type                                                                         type = GaussMapper_Type;
	Stg_Class_DeleteFunction*                                                 _delete = _GaussMapper_Delete;
	Stg_Class_PrintFunction*                                                   _print = _GaussMapper_Print;
	Stg_Class_CopyFunction*                                                     _copy = _GaussMapper_Copy;
	Stg_Component_DefaultConstructorFunction*                     _defaultConstructor = _GaussMapper_DefaultNew;
	Stg_Component_ConstructFunction*                                       _construct = _GaussMapper_AssignFromXML;
	Stg_Component_BuildFunction*                                               _build = _GaussMapper_Build;
	Stg_Component_InitialiseFunction*                                     _initialise = _GaussMapper_Initialise;
	Stg_Component_ExecuteFunction*                                           _execute = _GaussMapper_Execute;
	Stg_Component_DestroyFunction*                                           _destroy = _GaussMapper_Destroy;
	AllocationType                                                 nameAllocationType = NON_GLOBAL;
	IntegrationPointMapper_MapFunction*                                          _map = _GaussMapper_Map;
	IntegrationPointMapper_GetMaterialPointsSwarmsFunction*  _getMaterialPointsSwarms = _OneToOneMapper_GetMaterialPointsSwarms;
	IntegrationPointMapper_GetMaterialIndexOnFunction*            _getMaterialIndexOn = _OneToOneMapper_GetMaterialIndexOn;
	IntegrationPointMapper_GetExtensionOnFunction*                    _getExtensionOn = _OneToOneMapper_GetExtensionOn;
    IntegrationPointMapper_GetDoubleFromExtension*                  _getDoubleFromExtension = _OneToOneMapper_GetDoubleFromExtension;
	IntegrationPointMapper_GetDoubleFromMaterial*                  _getDoubleFromMaterial = _OneToOneMapper_GetDoubleFromMaterial;

	return _GaussMapper_New(  GAUSSMAPPER_PASSARGS  );
}

void _GaussMapper_Init( void* mapper ) {
	GaussMapper* self;

	self = (GaussMapper*)mapper;
}

void _GaussMapper_Delete( void* mapper ) {
	GaussMapper* self = (GaussMapper*)mapper;

	_OneToOneMapper_Delete( self );
}

void _GaussMapper_Print( void* mapper, Stream* stream ) {
	GaussMapper* self = (GaussMapper*)mapper;

	_OneToOneMapper_Print( self, stream );
}

void* _GaussMapper_Copy( const void* mapper, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	return _IntegrationPointMapper_Copy( mapper, dest, deep, nameExt, ptrMap );
}

void _GaussMapper_AssignFromXML( void* mapper, Stg_ComponentFactory* cf, void* data ) {
	GaussMapper* self = (GaussMapper*)mapper;
	
	_OneToOneMapper_AssignFromXML( mapper, cf, data );

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

	_GaussMapper_Init( self );
}

void _GaussMapper_Build( void* mapper, void* cf ) {
	GaussMapper* self = (GaussMapper*)mapper;

	_OneToOneMapper_Build( self, cf );
}

void _GaussMapper_Initialise( void* mapper, void* cf ) {
	GaussMapper* self = (GaussMapper*)mapper;

	_OneToOneMapper_Initialise( self, cf );
}

void _GaussMapper_Execute( void* mapper, void* data ) {
}

void _GaussMapper_Destroy( void* mapper, void* data ) {
	GaussMapper* self = (GaussMapper*)mapper;

	_OneToOneMapper_Destroy( self, data );
}

void _GaussMapper_Map( void* mapper ) {
	GaussMapper*				self = (GaussMapper*)mapper;
	IntegrationPointsSwarm*	integrationSwarm = self->integrationSwarm;
	MaterialPointsSwarm*		materialSwarm = self->materialSwarm;
	IntegrationPoint*			integrationPoint;
	MaterialPoint*				materialPoint;
	MaterialPointRef*			ref;
	Particle_Index				point_I;

	materialPoint = (MaterialPoint*)Swarm_ParticleAt( materialSwarm, 0 ); /* Get the first and only point */
	
	/* Map each point in integration to the single material point for its properties */
	for ( point_I = 0; point_I < integrationSwarm->particleLocalCount; point_I++ ) {
		integrationPoint = (IntegrationPoint*)Swarm_ParticleAt( integrationSwarm, point_I );

		ref = OneToOneMapper_GetMaterialRef( self, integrationPoint );
		ref->swarm_I = materialSwarm->swarmReg_I;
		ref->particle_I = 0;
	}
}


