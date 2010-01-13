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

#include "MaterialPoints.h"

#include <assert.h>
#include <string.h>
#include <math.h>

const Type OneToOneMapper_Type = "OneToOneMapper";

OneToOneMapper* _OneToOneMapper_New(  ONETOONEMAPPER_DEFARGS  ) {
	OneToOneMapper* result;

	result = (OneToOneMapper*)_IntegrationPointMapper_New(  INTEGRATIONPOINTMAPPER_PASSARGS  );

	return result;
}

void _OneToOneMapper_Init( void* mapper, MaterialPointsSwarm* materialSwarm ) {
	OneToOneMapper* self = (OneToOneMapper*)mapper;
	
	self->errorStream = Journal_MyStream( Error_Type, self );
	self->materialSwarm = materialSwarm;

	ExtensionManager_SetLockDown( self->integrationSwarm->particleExtensionMgr, False );
	self->materialRefHandle = ExtensionManager_Add( self->integrationSwarm->particleExtensionMgr, (Name)materialSwarm->name, sizeof(MaterialPointRef)  );
	ExtensionManager_SetLockDown( self->integrationSwarm->particleExtensionMgr, True );
}

void _OneToOneMapper_Delete( void* mapper ) {
	OneToOneMapper* self = (OneToOneMapper*)mapper;

	_IntegrationPointMapper_Delete( self );
}

void _OneToOneMapper_Print( void* mapper, Stream* stream ) {
	OneToOneMapper* self = (OneToOneMapper*)mapper;
	
	_IntegrationPointMapper_Print( self, stream );
	Stream_Indent( stream );
	Stg_Class_Print( self->materialSwarm, stream );
	Stream_UnIndent( stream );
}

void* _OneToOneMapper_Copy( void* mapper, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	OneToOneMapper* self = (OneToOneMapper*)mapper;
	OneToOneMapper* newCopy;
	
	newCopy = (OneToOneMapper*)_IntegrationPointMapper_Copy( self, dest, deep, nameExt, ptrMap );
	newCopy->materialSwarm = (MaterialPointsSwarm*)Stg_Class_Copy( self->materialSwarm, NULL, deep, nameExt, ptrMap );

	return newCopy;
}

void _OneToOneMapper_AssignFromXML( void* mapper, Stg_ComponentFactory* cf, void* data ) {
	OneToOneMapper*		self = (OneToOneMapper*)mapper;
	MaterialPointsSwarm*	materialSwarm;

	_IntegrationPointMapper_AssignFromXML( self, cf, data );

	materialSwarm = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)MaterialPointsSwarm_Type, MaterialPointsSwarm, True, data  );

	_OneToOneMapper_Init( self, materialSwarm );
}

void _OneToOneMapper_Build( void* mapper, void* data ) {
	OneToOneMapper* self = (OneToOneMapper*)mapper;

	_IntegrationPointMapper_Build( mapper, data );
	Stg_Component_Build( self->materialSwarm, data, False );
}

void _OneToOneMapper_Initialise( void* mapper, void* data ) {
    OneToOneMapper* self = (OneToOneMapper*)mapper;

    _IntegrationPointMapper_Initialise( mapper, data );
    Stg_Component_Initialise( self->materialSwarm, data, False );
}

void _OneToOneMapper_Execute( void* mapper, void* data ) {}

void _OneToOneMapper_Destroy( void* mapper, void* data ) {
	OneToOneMapper* self = (OneToOneMapper*)mapper;

	_IntegrationPointMapper_Destroy( self, data );
}

MaterialPointRef* OneToOneMapper_GetMaterialRef( void* mapper, void* integrationPoint ) {
	OneToOneMapper* self = (OneToOneMapper*)mapper;

	return (MaterialPointRef*)ExtensionManager_Get( self->integrationSwarm->particleExtensionMgr, integrationPoint, self->materialRefHandle );
}

MaterialPoint* OneToOneMapper_GetMaterialPoint( void* mapper, void* integrationPoint, MaterialPointsSwarm** materialSwarm ) {
	OneToOneMapper*		self = (OneToOneMapper*)mapper;
	MaterialPointRef*		ref;
	MaterialPointsSwarm*	swarm;
	MaterialPoint*			materialPoint; /* Assumes that material swarm holds Material particle or derivative */

	ref = OneToOneMapper_GetMaterialRef( self, integrationPoint );
	Journal_Firewall( ref != NULL, self->errorStream, "In func %s, no MaterialPointRef found on point\n", __func__ );

	swarm = (MaterialPointsSwarm*)Swarm_Register_At( Swarm_Register_GetSwarm_Register(), ref->swarm_I );
	Journal_Firewall( swarm != NULL, self->errorStream, "In func %s, no swarm found on for index %d\n", __func__, ref->swarm_I );

	if ( materialSwarm != NULL ) {
		*materialSwarm = swarm;
	}

	materialPoint = (MaterialPoint*)Swarm_ParticleAt( swarm, ref->particle_I );
	Journal_Firewall(
		materialPoint != NULL,
		self->errorStream, 
		"In func %s, no MaterialPoint found for swarm index %d, point index %d\n",
		__func__,
		ref->swarm_I,
		ref->particle_I );

	return materialPoint;
}

MaterialPointsSwarm** _OneToOneMapper_GetMaterialPointsSwarms( void* mapper, Index* count ) {
	OneToOneMapper*			self = (OneToOneMapper*)mapper;
	MaterialPointsSwarm**	result = Memory_Alloc_Array( MaterialPointsSwarm*, 1,  "Swarms" );

	result[0] = self->materialSwarm;
	*count = 1;

	return result;
}

Material_Index _OneToOneMapper_GetMaterialIndexOn( void* mapper, void* point ) {
	OneToOneMapper*	self = (OneToOneMapper*)mapper;
	MaterialPoint*		materialPoint; /* Assumes that material swarm holds Material particle or derivative */
	
	materialPoint = OneToOneMapper_GetMaterialPoint( self, point, NULL );

	return materialPoint->materialIndex;
}

void* _OneToOneMapper_GetExtensionOn( void* mapper, void* point, ExtensionInfo_Index extHandle ) {
	OneToOneMapper*		self = (OneToOneMapper*)mapper;
	MaterialPointsSwarm*	swarm;
	MaterialPoint*			materialPoint; /* Assumes that material swarm holds Material particle or derivative */

	materialPoint = OneToOneMapper_GetMaterialPoint( self, point, &swarm );

	return ExtensionManager_Get( swarm->particleExtensionMgr, materialPoint, extHandle );
}

double _OneToOneMapper_GetDoubleFromExtension(void* mapper, void* intPoint, ExtensionInfo_Index extHandle, int offs) {
    return *(double*)(IntegrationPointMapper_GetExtensionOn(mapper, intPoint, extHandle) + offs);
}

double _OneToOneMapper_GetDoubleFromMaterial(void* mapper, void* intPoint, ExtensionInfo_Index extHandle, int offs) {
    MaterialPointsSwarm *matSwarm;
    MaterialPoint *matPoint;

    matPoint = OneToOneMapper_GetMaterialPoint(mapper, intPoint, &matSwarm);
    return *(double*)(MaterialPointsSwarm_GetMaterialExtensionOn(matSwarm, matPoint, extHandle) + offs);
}
