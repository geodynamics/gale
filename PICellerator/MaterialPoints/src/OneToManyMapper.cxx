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

#include <assert.h>
#include <string.h>
#include <math.h>

#include "MaterialPoints.h"

const Type OneToManyMapper_Type = "OneToManyMapper";

OneToManyMapper* _OneToManyMapper_New( ONETOMANYMAPPER_DEFARGS ) {
	OneToManyMapper* result;

	result = (OneToManyMapper*)_IntegrationPointMapper_New( INTEGRATIONPOINTMAPPER_PASSARGS );
        result->swarm=_swarm;
        result->harmonic_average=_harmonic_average;
		
	return result;
}

void* _OneToManyMapper_DefaultNew( Name name ) {
  /* Variables set in this function */
  SizeT _sizeOfSelf = sizeof(OneToManyMapper);
  Type type = OneToManyMapper_Type;
  Stg_Class_DeleteFunction* _delete = _OneToManyMapper_Delete;
  Stg_Class_PrintFunction*_print = _OneToManyMapper_Print;
  Stg_Class_CopyFunction*_copy = _OneToManyMapper_Copy;
  Stg_Component_DefaultConstructorFunction*
    _defaultConstructor = _OneToManyMapper_DefaultNew;
  Stg_Component_ConstructFunction*
    _construct = _OneToManyMapper_AssignFromXML;
  Stg_Component_BuildFunction*
    _build = _OneToManyMapper_Build;
  Stg_Component_InitialiseFunction*
    _initialise = _OneToManyMapper_Initialise;
  Stg_Component_ExecuteFunction*
    _execute = _OneToManyMapper_Execute;
  Stg_Component_DestroyFunction*
    _destroy = _OneToManyMapper_Destroy;
  AllocationType nameAllocationType = NON_GLOBAL;
  IntegrationPointMapper_MapFunction* _map = _OneToManyMapper_Map;
  IntegrationPointMapper_GetMaterialPointsSwarmsFunction*
    _getMaterialPointsSwarms = _OneToManyMapper_GetMaterialPointsSwarms;
  IntegrationPointMapper_GetMaterialIndexOnFunction*
    _getMaterialIndexOn = _OneToManyMapper_GetMaterialIndexOn;
  IntegrationPointMapper_GetExtensionOnFunction*
    _getExtensionOn = _OneToManyMapper_GetExtensionOn;
  IntegrationPointMapper_GetDoubleFromExtension*
    _getDoubleFromExtension = _OneToManyMapper_GetDoubleFromExtension;
  IntegrationPointMapper_GetDoubleFromExtension*
    _getDoubleFromMaterial = _OneToManyMapper_GetDoubleFromMaterial;
  IntegrationPointsSwarm *_swarm=NULL;
  Bool _harmonic_average=True;

  return _OneToManyMapper_New( ONETOMANYMAPPER_PASSARGS );
}

void _OneToManyMapper_Init( void* mapper,
                            IntegrationPointsSwarm* swarm,
                            Bool harmonic_average ) {
	OneToManyMapper* self = (OneToManyMapper*)mapper;

	self->swarm = swarm;
	self->harmonic_average = harmonic_average;
}

void _OneToManyMapper_Delete( void* mapper ) {
	OneToManyMapper* self = (OneToManyMapper*)mapper;
	
	_IntegrationPointMapper_Delete( self );
}
void _OneToManyMapper_Print( void* mapper, Stream* stream ) {
	OneToManyMapper* self = (OneToManyMapper*)mapper;
	
	_IntegrationPointMapper_Print( self, stream );
	Stream_Indent( stream );
	Stg_Class_Print( self->swarm, stream );
	Stream_UnIndent( stream );
}
void* _OneToManyMapper_Copy( const void* mapper, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	OneToManyMapper* self = (OneToManyMapper*)mapper;
	OneToManyMapper* newCopy;
	
	newCopy = (OneToManyMapper*)_IntegrationPointMapper_Copy( self, dest, deep, nameExt, ptrMap );
	newCopy->swarm = (IntegrationPointsSwarm*)Stg_Class_Copy( self->swarm, NULL, deep, nameExt, ptrMap );

	return newCopy;
}

void _OneToManyMapper_AssignFromXML( void* mapper, Stg_ComponentFactory* cf, void* data ) {
	OneToManyMapper* self = (OneToManyMapper*)mapper;
	IntegrationPointsSwarm* swarm;
        Bool harmonic_average;
	
	_IntegrationPointMapper_AssignFromXML( self, cf, data );

	swarm =
          Stg_ComponentFactory_ConstructByKey(cf,self->name,"MappedSwarm",
                                              IntegrationPointsSwarm,True,data);
        harmonic_average=Stg_ComponentFactory_GetBool(cf,self->name,"HarmonicAverage",True);
	_OneToManyMapper_Init( self, swarm, harmonic_average );

}

void _OneToManyMapper_Build( void* mapper, void* data ) {
	OneToManyMapper* self = (OneToManyMapper*)mapper;

	_IntegrationPointMapper_Build( mapper, data );
	Stg_Component_Build( self->swarm, data, False );
	
}
void _OneToManyMapper_Initialise( void* mapper, void* data ) {
	OneToManyMapper* self = (OneToManyMapper*)mapper;

	_IntegrationPointMapper_Initialise( mapper, data );
	Stg_Component_Initialise( self->swarm, data, False );
}
void _OneToManyMapper_Execute( void* mapper, void* data ) {}

void _OneToManyMapper_Destroy( void* mapper, void* data ) {
	OneToManyMapper* self = (OneToManyMapper*)mapper;

	/*Stg_Class_Delete( self->swarm );*/

	_IntegrationPointMapper_Destroy( self, data );
}

/* Just call the embedded swarm remapper */
void _OneToManyMapper_Map(void* mapper) {
  OneToManyMapper* self = (OneToManyMapper*)mapper;
  IntegrationPointMapper_Map(self->swarm->mapper);
}

MaterialPointsSwarm** _OneToManyMapper_GetMaterialPointsSwarms( void* mapper, Index* count ) {
  OneToManyMapper* self = (OneToManyMapper*)mapper;
  return IntegrationPointMapper_GetMaterialPointsSwarms(self->swarm->mapper,
                                                        count);
}

Material_Index _OneToManyMapper_GetMaterialIndexOn( void* mapper, void* point ) {
  abort();
  return 0;
}

void* _OneToManyMapper_GetExtensionOn( void* mapper, void* point, ExtensionInfo_Index extHandle ) {
  /* This method cannot work with a one-to-many mapping. */
  abort();
  return NULL;
}

double _OneToManyMapper_GetDoubleFromExtension(void* mapper, void* intPoint, ExtensionInfo_Index extHandle, int offs) {
  abort();
  return 0;
    // OneToManyMapper *self = (OneToManyMapper*)mapper;
    // OneToManyRef *ref;
    // double v = 0.0, c;
    // int ii;

    // ref = OneToManyMapper_GetMaterialRef(mapper, intPoint);
    // for(ii = 0; ii < ref->numParticles; ii++) {
    //   c = *(double*)((char*)MaterialPointsSwarm_GetExtensionAt(self->swarm, ref->particleInds[ii], extHandle) + offs);
    //   v += ((double)ref->weights[ii])*c;
    // }

    // return v;
}

double _OneToManyMapper_GetDoubleFromMaterial(void* mapper, void* intPoint, ExtensionInfo_Index extHandle, int offs) {
  abort();
  return 0;
    // OneToManyMapper *self = (OneToManyMapper*)mapper;
    // OneToManyRef *ref;
    // double v = 0.0, c;
    // int ii;

    // ref = OneToManyMapper_GetMaterialRef(mapper, intPoint);
    // for(ii = 0; ii < ref->numParticles; ii++) {
    //   c = *(double*)((char*)MaterialPointsSwarm_GetMaterialExtensionAt(self->swarm, ref->particleInds[ii], extHandle) + offs);
    //   v += ((double)ref->weights[ii])*c;
    // }

    // return v;
}

