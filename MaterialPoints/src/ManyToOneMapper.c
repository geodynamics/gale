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
#include <StgFEM/StgFEM.h>

#include <PICellerator/Voronoi/Voronoi.h>
#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>

#include "MaterialPoints.h"

#include <assert.h>
#include <string.h>
#include <math.h>

const Type ManyToOneMapper_Type = "ManyToOneMapper";

ManyToOneMapper* _ManyToOneMapper_New(
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
		Bool                                            initFlag,
		IntegrationPointsSwarm*                                         integrationSwarm,
		MaterialPointsSwarm**                                           materialSwarms,
		Index                                                           materialSwarmCount ) 
{
	ManyToOneMapper* result;

	result = (ManyToOneMapper*)_IntegrationPointMapper_New(
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
			integrationSwarm );

	if (initFlag) {
		_ManyToOneMapper_Init( result, integrationSwarm, materialSwarms, materialSwarmCount );
	}

	return result;
}

void _ManyToOneMapper_Init( 
		void*                   mapper,
		IntegrationPointsSwarm* integrationSwarm,
		MaterialPointsSwarm**   materialSwarms,
		Index                   materialSwarmCount )
{
	ManyToOneMapper* self = (ManyToOneMapper*)mapper;

	self->materialSwarms = materialSwarms;
	self->materialSwarmCount = materialSwarmCount;

	_IntegrationPointMapper_Init( self, integrationSwarm );
}

void _ManyToOneMapper_Delete( void* mapper ) {
	ManyToOneMapper* self = (ManyToOneMapper*)mapper;

	int i;

	if ( self->materialSwarms != NULL ) {
		for ( i = 0; i < self->materialSwarmCount; ++i )  {
			Stg_Class_Delete( self->materialSwarms[i] );
		}
		Memory_Free( self->materialSwarms );
	}
	
	_IntegrationPointMapper_Delete( self );
}
void _ManyToOneMapper_Print( void* mapper, Stream* stream ) {
	ManyToOneMapper* self = (ManyToOneMapper*)mapper;
	int i;

	_IntegrationPointMapper_Print( self, stream );
	if ( self->materialSwarms != NULL ) {
		Stream_Indent( stream );
		for ( i = 0; i < self->materialSwarmCount; ++i ) {
			Print( self->materialSwarms[i], stream );
		}
		Stream_UnIndent( stream );
	}
}
void* _ManyToOneMapper_Copy( void* mapper, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	ManyToOneMapper* self = (ManyToOneMapper*)mapper;
	ManyToOneMapper* newCopy;
	
	newCopy = (ManyToOneMapper*)_IntegrationPointMapper_Copy( self, dest, deep, nameExt, ptrMap );

	if ( self->materialSwarms != NULL ) {
		newCopy->materialSwarms = PtrMap_Find( ptrMap, self->materialSwarms );
		if ( newCopy->materialSwarms == NULL ) {
			newCopy->materialSwarms = Memory_Alloc_Array( MaterialPointsSwarm*, self->materialSwarmCount, "componentList" );
			PtrMap_Append( ptrMap, self->materialSwarms, newCopy->materialSwarms );
		}
		newCopy->materialSwarmCount = self->materialSwarmCount;
	}
	else {
		newCopy->materialSwarms = NULL;
		newCopy->materialSwarmCount = 0;
	}

	return newCopy;
}

void _ManyToOneMapper_Construct( void* mapper, Stg_ComponentFactory* cf, void* data ) {
	ManyToOneMapper* self = (ManyToOneMapper*)mapper;

	MaterialPointsSwarm** materialSwarms;
	int i;
	
	_IntegrationPointMapper_Construct( self, cf, data );

	materialSwarms = (MaterialPointsSwarm**)Stg_ComponentFactory_ConstructByList( 
					cf, 
					self->name, 
					IntegrationPointsSwarm_Type, 
					Stg_ComponentFactory_Unlimited, 
					IntegrationPointsSwarm,
					True,
					&(self->materialSwarmCount), data );

	Journal_Firewall( 
		self->materialSwarmCount < 1,
		Journal_Register( Error_Type, self->type ),
		"In func %s, there must be at least one swarm in the material swarm list!\n", __func__ );

	_ManyToOneMapper_Init( self, self->integrationSwarm, materialSwarms, self->materialSwarmCount );

	/* Each integration point will have a reference to a material particle (one for each swarm) */
	ExtensionManager_SetLockDown( self->integrationSwarm->particleExtensionMgr, False );
	for ( i = 0; i < self->materialSwarmCount; ++i ) {
		ExtensionManager_Add( 
			self->integrationSwarm->particleExtensionMgr,
			materialSwarms[i]->name, 
			sizeof(MaterialPointRef) );
	}
	ExtensionManager_SetLockDown( self->integrationSwarm->particleExtensionMgr, True );
}

void _ManyToOneMapper_Build( void* mapper, void* data ) {
	ManyToOneMapper* self = (ManyToOneMapper*)mapper;
	int i;

	_IntegrationPointMapper_Build( mapper, data );
	
	for ( i = 0 ; i < self->materialSwarmCount; ++i ) {
		Stg_Component_Build( self->materialSwarms[i], data, False );
	}
	
}
void _ManyToOneMapper_Initialise( void* mapper, void* data ) {
	ManyToOneMapper* self = (ManyToOneMapper*)mapper;
	int i;

	_IntegrationPointMapper_Initialise( mapper, data );
	for ( i = 0; i < self->materialSwarmCount; ++i ) {
		Stg_Component_Initialise( self->materialSwarms[i], data, False );
	}
}
void _ManyToOneMapper_Execute( void* mapper, void* data ) {

}
void _ManyToOneMapper_Destroy( void* mapper, void* data ) {

}

MaterialPointsSwarm** ManyToOneMapper_GetMaterialPointsSwarms( void* mapper, Index* count ) {
	ManyToOneMapper*       self   = (ManyToOneMapper*)mapper;
	MaterialPointsSwarm** result = Memory_Alloc_Array( MaterialPointsSwarm*, self->materialSwarmCount, "Swarms" );
	Index                  i;
	
	*count = self->materialSwarmCount;

	for ( i = 0; i < self->materialSwarmCount; ++i ) {
		result[i] = self->materialSwarms[i];
	}

	return result;
}
