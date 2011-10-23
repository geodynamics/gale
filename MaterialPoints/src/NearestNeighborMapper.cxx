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
#include <limits>

#include "MaterialPoints.h"

const Type NearestNeighborMapper_Type = "NearestNeighborMapper";

NearestNeighborMapper* _NearestNeighborMapper_New( NEARESTNEIGHBORMAPPER_DEFARGS ) {
	NearestNeighborMapper* result;

	result = (NearestNeighborMapper*)_IntegrationPointMapper_New( INTEGRATIONPOINTMAPPER_PASSARGS );
        result->swarm=_swarm;

	return result;
}

void* _NearestNeighborMapper_DefaultNew( Name name ) {
  /* Variables set in this function */
  SizeT _sizeOfSelf = sizeof(NearestNeighborMapper);
  Type type = NearestNeighborMapper_Type;
  Stg_Class_DeleteFunction* _delete = _NearestNeighborMapper_Delete;
  Stg_Class_PrintFunction*_print = _NearestNeighborMapper_Print;
  Stg_Class_CopyFunction*_copy = _NearestNeighborMapper_Copy;
  Stg_Component_DefaultConstructorFunction*
    _defaultConstructor = _NearestNeighborMapper_DefaultNew;
  Stg_Component_ConstructFunction*
    _construct = _NearestNeighborMapper_AssignFromXML;
  Stg_Component_BuildFunction*
    _build = _NearestNeighborMapper_Build;
  Stg_Component_InitialiseFunction*
    _initialise = _NearestNeighborMapper_Initialise;
  Stg_Component_ExecuteFunction*
    _execute = _NearestNeighborMapper_Execute;
  Stg_Component_DestroyFunction*
    _destroy = _NearestNeighborMapper_Destroy;
  AllocationType nameAllocationType = NON_GLOBAL;
  IntegrationPointMapper_MapFunction* _map = _NearestNeighborMapper_Map;
  IntegrationPointMapper_GetMaterialPointsSwarmsFunction*
    _getMaterialPointsSwarms = _NearestNeighborMapper_GetMaterialPointsSwarms;
  IntegrationPointMapper_GetMaterialIndexOnFunction*
    _getMaterialIndexOn = _NearestNeighborMapper_GetMaterialIndexOn;
  IntegrationPointMapper_GetExtensionOnFunction*
    _getExtensionOn = _NearestNeighborMapper_GetExtensionOn;
  IntegrationPointMapper_GetDoubleFromExtension*
    _getDoubleFromExtension = _NearestNeighborMapper_GetDoubleFromExtension;
  IntegrationPointMapper_GetDoubleFromExtension*
    _getDoubleFromMaterial = _NearestNeighborMapper_GetDoubleFromMaterial;
  IntegrationPointsSwarm *_swarm=NULL;

  return _NearestNeighborMapper_New( NEARESTNEIGHBORMAPPER_PASSARGS );
}

void _NearestNeighborMapper_Init( void* mapper, IntegrationPointsSwarm* swarm ) {
	NearestNeighborMapper* self = (NearestNeighborMapper*)mapper;

	self->errorStream = Journal_MyStream( Error_Type, self );
	self->swarm = swarm;
}

void _NearestNeighborMapper_Delete( void* mapper ) {
	NearestNeighborMapper* self = (NearestNeighborMapper*)mapper;
	
	_IntegrationPointMapper_Delete( self );
}
void _NearestNeighborMapper_Print( void* mapper, Stream* stream ) {
	NearestNeighborMapper* self = (NearestNeighborMapper*)mapper;
	
	_IntegrationPointMapper_Print( self, stream );
	Stream_Indent( stream );
	Stg_Class_Print( self->swarm, stream );
	Stream_UnIndent( stream );
}
void* _NearestNeighborMapper_Copy( const void* mapper, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	NearestNeighborMapper* self = (NearestNeighborMapper*)mapper;
	NearestNeighborMapper* newCopy;
	
	newCopy = (NearestNeighborMapper*)_IntegrationPointMapper_Copy( self, dest, deep, nameExt, ptrMap );
	newCopy->swarm = (IntegrationPointsSwarm*)Stg_Class_Copy( self->swarm, NULL, deep, nameExt, ptrMap );

	return newCopy;
}

void _NearestNeighborMapper_AssignFromXML( void* mapper, Stg_ComponentFactory* cf, void* data ) {
	NearestNeighborMapper* self = (NearestNeighborMapper*)mapper;
	IntegrationPointsSwarm* swarm;
	
	_IntegrationPointMapper_AssignFromXML( self, cf, data );

	swarm =
          Stg_ComponentFactory_ConstructByKey(cf,self->name,"MappedSwarm",
                                              IntegrationPointsSwarm,True,data);

	_NearestNeighborMapper_Init( self, swarm );
}

void _NearestNeighborMapper_Build( void* mapper, void* data ) {
	NearestNeighborMapper* self = (NearestNeighborMapper*)mapper;

	_IntegrationPointMapper_Build( mapper, data );
	Stg_Component_Build( self->swarm, data, False );
	
}
void _NearestNeighborMapper_Initialise( void* mapper, void* data ) {
	NearestNeighborMapper* self = (NearestNeighborMapper*)mapper;

	_IntegrationPointMapper_Initialise( mapper, data );
	Stg_Component_Initialise( self->swarm, data, False );
}
void _NearestNeighborMapper_Execute( void* mapper, void* data ) {}

void _NearestNeighborMapper_Destroy( void* mapper, void* data ) {
	NearestNeighborMapper* self = (NearestNeighborMapper*)mapper;

	_IntegrationPointMapper_Destroy( self, data );
}

/* Just call the embedded swarm remapper */
void _NearestNeighborMapper_Map(void* mapper) {
  NearestNeighborMapper* self = (NearestNeighborMapper*)mapper;
  IntegrationPointMapper_Map(self->swarm->mapper);
}

MaterialPointsSwarm** _NearestNeighborMapper_GetMaterialPointsSwarms
( void* mapper, Index* count ) {
  NearestNeighborMapper* self = (NearestNeighborMapper*)mapper;
  return IntegrationPointMapper_GetMaterialPointsSwarms(self->swarm->mapper,
                                                        count);
}

Material_Index _NearestNeighborMapper_GetMaterialIndexOn
( void* mapper, void* point ) {
  abort();
  return 0;
}

void* _NearestNeighborMapper_GetExtensionOn( void* mapper, void* point,
                                             ExtensionInfo_Index extHandle ) {
  abort();
  return NULL;
}

double _NearestNeighborMapper_GetDoubleFromExtension
(void* mapper,
 void* intPoint,
 ExtensionInfo_Index extHandle,
 int offs) {
  abort();
  return 0;
}

double _NearestNeighborMapper_GetDoubleFromMaterial
(void* mapper,
 void* intPoint,
 ExtensionInfo_Index extHandle,
 int offs) {
  abort();
  return 0;
}

/* This does a search over all of the particles in the swarm in the
   element to find the one that is closest to the gauss point.  There
   may be more efficient ways of doing this, but this works for
   now. */
int NearestNeighbor_FindNeighbor(void* mapper, const Element_LocalIndex &lElement_I,
                                 const int &cell_I, double *xi, const int &dim)
{
  NearestNeighborMapper* self = (NearestNeighborMapper*)mapper;
  IntegrationPointsSwarm* swarm=self->swarm;
    
  int cellParticleCount(swarm->cellParticleCountTbl[cell_I]);

  Journal_Firewall(cellParticleCount!=0,
                   Journal_Register(Error_Type,(Name)NearestNeighborMapper_Type),
                   "In func %s: cellParticleCount is 0.\n",__func__);

  double min_dist(std::numeric_limits<double>::max());
  int nearest_particle(-1);

  for(int cParticle_I=0; cParticle_I<cellParticleCount; cParticle_I++)
    {
      IntegrationPoint* particle=
        (IntegrationPoint*)Swarm_ParticleInCellAt(swarm,cell_I,cParticle_I);
                                                  
      double dist=
        StGermain_DistanceBetweenPoints(xi,particle->xi,dim);
      if(dist<min_dist)
        {
          nearest_particle=cParticle_I;
          min_dist=dist;
        }
    }
  return nearest_particle;
}

/* A convenience function to replace the gauss swarm and particle with
   the mapped swarm and nearest particle */

void NearestNeighbor_Replace
(IntegrationPointsSwarm **swarm, IntegrationPoint **particle,
 int *particle_index, const Element_LocalIndex &lElement_I,
 const int &dim)
{
  if(Stg_Class_IsInstance((*swarm)->mapper,NearestNeighborMapper_Type))
    {
      IntegrationPointsSwarm* NNswarm=
        ((NearestNeighborMapper*)((*swarm)->mapper))->swarm;
      int NNcell_I=
        CellLayout_MapElementIdToCellId(NNswarm->cellLayout,lElement_I);
      int nearest_particle=
        NearestNeighbor_FindNeighbor((*swarm)->mapper,lElement_I,
                                     NNcell_I,(*particle)->xi,dim);
      IntegrationPoint *NNparticle=
        (IntegrationPoint*)Swarm_ParticleInCellAt(NNswarm,
                                                  NNcell_I,
                                                  nearest_particle);
      *swarm=NNswarm;
      *particle=NNparticle;
      *particle_index=NNswarm->cellParticleTbl[NNcell_I][nearest_particle];
    }
}
