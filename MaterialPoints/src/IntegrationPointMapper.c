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

const Type IntegrationPointMapper_Type = "IntegrationPointMapper";

IntegrationPointMapper* _IntegrationPointMapper_New(
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
    IntegrationPointsSwarm*                                         integrationSwarm )
{
    IntegrationPointMapper* self;

    self = (IntegrationPointMapper*)_Stg_Component_New(
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
        name,
        NON_GLOBAL );

    self->_map                           = _map;
    self->_getMaterialPointsSwarms       = _getMaterialPointsSwarms;
    self->_getMaterialIndexOn            = _getMaterialIndexOn;
    self->_getExtensionOn                = _getExtensionOn;

    if (initFlag) {
        self->isConstructed = True;
        _IntegrationPointMapper_Init( self, integrationSwarm );
    }

    return self;
}

void _IntegrationPointMapper_Init( void* mapper, IntegrationPointsSwarm* integrationSwarm ) {
    IntegrationPointMapper* self = (IntegrationPointMapper*)mapper;

    self->integrationSwarm = integrationSwarm;
}

void _IntegrationPointMapper_Delete( void* mapper ) {
    IntegrationPointMapper* self = (IntegrationPointMapper*)mapper;

    Stg_Class_Delete( self->integrationSwarm );

    _Stg_Component_Delete( self );
}
void _IntegrationPointMapper_Print( void* mapper, Stream* stream ) {
    IntegrationPointMapper* self = (IntegrationPointMapper*)mapper;
	
    Journal_Printf( stream, "IntegrationPointMapper (ptr): %s\n", (void*)self );
    _Stg_Component_Print( self, stream );

    Stream_Indent( stream );

    Stg_Class_Print( self->integrationSwarm, stream );

    Stream_UnIndent( stream );
}
void* _IntegrationPointMapper_Copy( void* mapper, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
    IntegrationPointMapper* self = (IntegrationPointMapper*)mapper;
    IntegrationPointMapper* newCopy = (IntegrationPointMapper*)_Stg_Component_Copy( self, dest, deep, nameExt, ptrMap );

    newCopy->integrationSwarm = Stg_Class_Copy( self->integrationSwarm, NULL, deep, nameExt, ptrMap );

    return newCopy;
}

void _IntegrationPointMapper_Construct( void* mapper, Stg_ComponentFactory* cf, void* data ) {
    IntegrationPointMapper* self = (IntegrationPointMapper*)mapper;

    IntegrationPointsSwarm* integrationSwarm;

    self->context = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Context", PICelleratorContext, False, data );
    if( !self->context ) 
        self->context = Stg_ComponentFactory_ConstructByName( cf, "context", PICelleratorContext, True, data );

    integrationSwarm = Stg_ComponentFactory_ConstructByKey( 
        cf, 
        self->name, 
        IntegrationPointsSwarm_Type, 
        IntegrationPointsSwarm,
        True,
        data  );
	
    _IntegrationPointMapper_Init( self, integrationSwarm );
}

void _IntegrationPointMapper_Build( void* mapper, void* data ) {
}

void _IntegrationPointMapper_Initialise( void* mapper, void* data ) {
}

void _IntegrationPointMapper_Execute( void* mapper, void* data ) {
}

void _IntegrationPointMapper_Destroy( void* mapper, void* data ) {
}

void IntegrationPointMapper_Map( void* mapper ) {
    IntegrationPointMapper* self = (IntegrationPointMapper*)mapper;

    self->_map( self );
}

MaterialPointsSwarm** IntegrationPointMapper_GetMaterialPointsSwarms( void* mapper, Index* count ) {
    IntegrationPointMapper* self = (IntegrationPointMapper*)mapper;

    return self->_getMaterialPointsSwarms( mapper, count );
}

Material_Index IntegrationPointMapper_GetMaterialIndexOnFunc( void* mapper, void* point ) {
    IntegrationPointMapper* self = (IntegrationPointMapper*)mapper;

    return self->_getMaterialIndexOn( mapper, point );
}
void* IntegrationPointMapper_GetExtensionOnFunc( void* mapper, void* point, ExtensionInfo_Index extHandle ) {
    IntegrationPointMapper* self = (IntegrationPointMapper*)mapper;

    return self->_getExtensionOn( mapper, point, extHandle );
}
Material_Index IntegrationPointMapper_GetMaterialIndexAtFunc( void* mapper, Index point_I ) {
    IntegrationPointMapper* self  = (IntegrationPointMapper*)mapper;
    void*                   point = Swarm_ParticleAt( self->integrationSwarm, point_I );

    Journal_Firewall(
        point != NULL,
        Journal_MyStream( Error_Type, self ),
        "In func %s, no point in swarm of index %d\n",
        __func__,
        point_I );

    return IntegrationPointMapper_GetMaterialIndexOn( mapper, point );
}
void* IntegrationPointMapper_GetExtensionAtFunc( void* mapper, Index point_I, ExtensionInfo_Index extHandle ) {
    IntegrationPointMapper* self  = (IntegrationPointMapper*)mapper;
    void*                   point = Swarm_ParticleAt( self->integrationSwarm, point_I );

    Journal_Firewall(
        point != NULL,
        Journal_MyStream( Error_Type, self ),
        "In func %s, no point in swarm of index %d\n",
        __func__,
        point_I );
	
    return IntegrationPointMapper_GetExtensionOn( mapper, point, extHandle );
}
