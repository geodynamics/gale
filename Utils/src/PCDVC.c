/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2006, Monash Cluster Computing 
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
** 		* Redistributions of source code must retain the above copyright notice, 
** 			this list of conditions and the following disclaimer.
** 		* Redistributions in binary form must reproduce the above copyright 
**			notice, this list of conditions and the following disclaimer in the 
**			documentation and/or other materials provided with the distribution.
** 		* Neither the name of the Monash University nor the names of its contributors 
**			may be used to endorse or promote products derived from this software 
**			without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
** THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
** PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS 
** BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
** CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
** SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
** HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
** LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
** OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
**
** Contact:
*%		Louis Moresi - Louis.Moresi@sci.monash.edu.au
*%
** Author:
**              Mirko Velic - Mirko.Velic@sci.monash.edu.au
**              Patrick Sunter - patrick@vpac.org
**              Julian Giordani - julian.giordani@sci.monash.edu.au
**          
**  Assumptions:
**  	 I am assuming that the xi's (local coords) on the IntegrationPoint particles
**       are precalculated. i.e. We assume the Coincident mapper is being used.
**
**  Notes:
**         The PCDVC class should really be a class the next level up here.
**	   We should be able to swap out the WeightsCalculator_CalculateAll function instead of just setting
**                 a pointer inside that function.
**
**         If the function  getIntParticleMaterialRef_PointingToMaterialParticle ever gets called
**         then someone has messed up the mapping between integration and material points. This function
**         is potentially slow as it traverses the whole swarm. This should be avoided.
**
**         We do not allow particle deletion in interface cells (cells that have more than one type of material
**         in them). Splitting is optional. This may be inadequate. We may need to do some handling of the neighbours
**         to interface cells as well, in order to preserve particle density about an interface.
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


/****************************************************************************************************************

  The algorithm here-in uses the DVCWeights module to compute a discrete voronoi diagram per FEM cell given a set of local
  particle positions, in 3D and 2D. The volumes of the Voronoi regions are used as integration weights for
  the integration point swarm and the integration points are the centroids of the same volumes.

  The volumes are also used as criteria for splitting and deleting particles.
  At the moment we are only deleting or adding *one* particle per cell per time-step: This
  will be changed shortlyas it may be inadequate.

  For a description of the Voronoi algorithm, see the article by Velic et.al.
     "A Fast Robust Algorithm for computing Discrete Voronoi Diagrams in N-dimensions"
  
*****************************************************************************************************************/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/Voronoi/Voronoi.h>
#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>
#include <PICellerator/MaterialPoints/MaterialPoints.h>

//#include <PICellerator/PICellerator.h>
#include "types.h"

#include "PCDVC.h"

#include <assert.h>
#include <string.h>
#include <math.h>


/* Textual name of this class */
const Type PCDVC_Type = "PCDVC";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/
PCDVC* PCDVC_New( Name name, Dimension_Index dim ) {
	PCDVC* self = (PCDVC*) _PCDVC_DefaultNew( name );

	PCDVC_InitAll( self, dim );
	return self;
}

PCDVC* _PCDVC_New(
		SizeT                                 _sizeOfSelf, 
		Type                                  type,
		Stg_Class_DeleteFunction*             _delete,
		Stg_Class_PrintFunction*              _print,
		Stg_Class_CopyFunction*               _copy, 
		Stg_Component_DefaultConstructorFunction* _defaultConstructor,
		Stg_Component_ConstructFunction*      _construct,
		Stg_Component_BuildFunction*          _build,
		Stg_Component_InitialiseFunction*     _initialise,
		Stg_Component_ExecuteFunction*        _execute,
		Stg_Component_DestroyFunction*        _destroy,		
		WeightsCalculator_CalculateFunction*  _calculate,
		Name                                  name )
{
	PCDVC* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(PCDVC) );
	/* Initialise the parent class. Every class has a parent, bar Stg_Component, which needs to be called */
	self = (PCDVC*)_DVCWeights_New( 
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
			_calculate,
			name );

	
	/* General info */

	/* Virtual Info */
	return self;
}

void _PCDVC_Init( void* pcdvc, MaterialPointsSwarm* mps,  double upT, double lowT, int maxDeletions, int maxSplits, Bool splitInInterfaceCells, int *res ) {
	PCDVC* self = (PCDVC*)pcdvc;
	
	self->materialPointsSwarm = mps;
	self->upperT = upT;
	self->lowerT = lowT;
	self->maxDeletions = maxDeletions;
	self->maxSplits = maxSplits;
	self->resX = res[I_AXIS];
	self->resY = res[J_AXIS];
	self->resZ = res[K_AXIS];
	self->splitInInterfaceCells = splitInInterfaceCells;

}

void PCDVC_InitAll( void* pcdvc, Dimension_Index dim ) {
	PCDVC* self = (PCDVC*)pcdvc;
	WeightsCalculator_InitAll( self, dim );
}

/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _PCDVC_Delete( void* pcdvc ) {
	PCDVC* self = (PCDVC*)pcdvc;
	/* Delete parent */
	_DVCWeights_Delete( self );
}


void _PCDVC_Print( void* pcdvc, Stream* stream ) {
	PCDVC* self = (PCDVC*)pcdvc;
	/* Print parent */
	_DVCWeights_Print( self, stream );
}



void* _PCDVC_Copy( void* pcdvc, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	PCDVC*	self = (PCDVC*)pcdvc;
	PCDVC*	newPCDVC;
	
	newPCDVC = (PCDVC*)_DVCWeights_Copy( self, dest, deep, nameExt, ptrMap );
	return (void*)newPCDVC;
}

void* _PCDVC_DefaultNew( Name name ) {
	return (void*) _PCDVC_New(
			sizeof(PCDVC),
			PCDVC_Type,
			_PCDVC_Delete,
			_PCDVC_Print,
			_PCDVC_Copy,
			_PCDVC_DefaultNew,
			_PCDVC_Construct,
			_PCDVC_Build,
			_PCDVC_Initialise,
			_PCDVC_Execute,
			_PCDVC_Destroy,
			_PCDVC_Calculate,
			name );
}


void _PCDVC_Construct( void* pcdvc, Stg_ComponentFactory* cf, void *data ) {

	PCDVC*	     self          = (PCDVC*) pcdvc;
	MaterialPointsSwarm*       materialPointsSwarm;
	double upT, lowT;
	int defaultResolution;
	int resolution[3];
	int maxD, maxS;
	Bool splitInInterfaceCells;
	_DVCWeights_Construct( self, cf, data );
	materialPointsSwarm = Stg_ComponentFactory_ConstructByKey( cf, self->name, "MaterialPointsSwarm", MaterialPointsSwarm, True, data );
	Stream*  stream = Journal_Register( Info_Type, materialPointsSwarm->type );
	

	upT = Stg_ComponentFactory_GetDouble( cf, self->name, "upperT", 25 );
	lowT = Stg_ComponentFactory_GetDouble( cf, self->name, "lowerT", 0.6 );
	maxD = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "maxDeletions", 2);
	maxS = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "maxSplits", 3);
	splitInInterfaceCells = Stg_ComponentFactory_GetBool( cf, self->name, "splitInInterfaceCells", False);
	defaultResolution = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolution", 10 );
	resolution[ I_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionX", defaultResolution );
	resolution[ J_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionY", defaultResolution );
	resolution[ K_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionZ", defaultResolution );

	if(upT < lowT){
	      lowT = 0.6;
	      upT = 25;
	      Journal_Printf( stream,"On Proc %d: In func %s(): WARNING!! lowT and upT have been reset to some more reasonable values. (lowT = 0.6, upT = 25) now!",materialPointsSwarm->myRank, __func__);
	}

	_PCDVC_Init( self, materialPointsSwarm,  upT, lowT, maxD, maxS, splitInInterfaceCells, resolution );
}

void _PCDVC_Build( void* pcdvc, void* data ) {
	PCDVC*	self = (PCDVC*)pcdvc;
	_DVCWeights_Build( self, data );
}
void _PCDVC_Initialise( void* pcdvc, void* data ) {
	PCDVC*	self = (PCDVC*)pcdvc;
	_DVCWeights_Initialise( self, data );
}
void _PCDVC_Execute( void* pcdvc, void* data ) {
	PCDVC*	self = (PCDVC*)pcdvc;
	_DVCWeights_Execute( self, data );
}
void _PCDVC_Destroy( void* pcdvc, void* data ) {
	PCDVC*	self = (PCDVC*)pcdvc;
	_DVCWeights_Destroy( self, data );
}

/*-------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
/* this function loops over the intSwarm in order to find the MaterialPointRef for the integration point particle
   that points to the matLastParticle_IndexOnCPU of a material particle. It is used here in the case that the last integration point
   particle DOES NOT point to the last material point particle, as it should.
   If this function is being called ever, then some other module/component somewhere has messed up the mapping between the integration Swarm and the material Swarm*/
MaterialPointRef* getIntParticleMaterialRef_PointingToMaterialParticle( IntegrationPointsSwarm*  intSwarm, Particle_Index matLastParticle_IndexOnCPU ){
      IntegrationPoint* intTestParticle;
      MaterialPointRef*       ref;
      int i;
      Stream*  stream = Journal_Register( Info_Type, intSwarm->type );
      Journal_Printf( stream,"\n\n\e[31m\nOn Proc %d: In func %s(): WARNING!! If this function is being called, then some other module/component, somewhere, has messed up the mapping between the integration Swarm and the material Swarm\n\n", intSwarm->myRank, __func__);
      Journal_Printf( stream,"This function is potentially slow. Someone should fix the offending module so that it doesn not mess up the ordering\n\n");
      Journal_Printf( stream,"\e[0;32m");
      for(i=0;i<intSwarm->particleLocalCount;i++){
	    intTestParticle =  (IntegrationPoint*)Swarm_ParticleAt( intSwarm, i);
	    ref = OneToOneMapper_GetMaterialRef(intSwarm->mapper, intTestParticle);
	    if(ref->particle_I == matLastParticle_IndexOnCPU){
		  break;
	    }
      }
      return ref;
}
/****************************************************************************************************
 This function deletes integration particle number intParticleToRemove_IndexWithinCell in lCell_I
 as well as the corresponding material point.
 It assumes the one to one coincidence mapper is being used. i.e. for every integration point 
 there is one material point and vice versa.
******************************************************************************************************/
void deleteIntParticleByIndexWithinCell( IntegrationPointsSwarm*  intSwarm, MaterialPointsSwarm* matSwarm, Cell_LocalIndex lCell_I, Particle_Index intParticleToRemove_IndexWithinCell ) {
      MaterialPointRef*       ref;
      MaterialPointRef*       refTestLastParticle;
      MaterialPointRef*       ref_moved_particle;
      MaterialPointRef*   refToLastMatParticle;
      int refToLastMatParticleFlag = 0;
      Particle_Index    intParticleToRemove_IndexOnCPU;/* the particle number within the swarm on the local CPU */
      Particle_Index    matParticleToRemove_IndexOnCPU;
      IntegrationPoint* intParticleToRemove;
      MaterialPoint*    matParticleToRemove;
      //Particle_Index    intParticleToRemove_IndexWithinCell;/* the number of the particle within the cell */
      Particle_Index    matParticleToRemove_IndexWithinCell;
      IntegrationPoint* intSwarmLastParticle; /* The last particle in the swarm on the local CPU */
      MaterialPoint*    matSwarmLastParticle;
      Particle_Index    intLastParticle_IndexOnCPU;
      Particle_Index    matLastParticle_IndexOnCPU;
      Particle_Index    intLastParticle_IndexWithinCell;
      Particle_Index    matLastParticle_IndexWithinCell;

      SizeT                 intparticleSize     = intSwarm->particleExtensionMgr->finalSize;
      SizeT                 matparticleSize     = matSwarm->particleExtensionMgr->finalSize;


      intParticleToRemove_IndexOnCPU  = intSwarm->cellParticleTbl[ lCell_I ][ intParticleToRemove_IndexWithinCell ];
      intParticleToRemove =  (IntegrationPoint*)Swarm_ParticleAt( intSwarm, intParticleToRemove_IndexOnCPU);
      ref = OneToOneMapper_GetMaterialRef(intSwarm->mapper, intParticleToRemove); /* so we can get the material point */
      matParticleToRemove_IndexOnCPU = ref->particle_I;
      matParticleToRemove =  (MaterialPoint*)Swarm_ParticleAt( matSwarm, matParticleToRemove_IndexOnCPU);
      matParticleToRemove_IndexWithinCell = Swarm_GetParticleIndexWithinCell( matSwarm, matParticleToRemove->owningCell, matParticleToRemove_IndexOnCPU);

      /* The particles are homeless after these function-calls and for this reason alone must be destroyed */
      Swarm_RemoveParticleFromCell( intSwarm, lCell_I, intParticleToRemove_IndexWithinCell );
      Swarm_RemoveParticleFromCell( matSwarm, matParticleToRemove->owningCell, matParticleToRemove_IndexWithinCell );

      /* Copy over particle to remove with the last particle in the particle array - as long as it isn't the last one  */
      intLastParticle_IndexOnCPU = intSwarm->particleLocalCount - 1; /* lastParticle_I is the last particle in the array of particles residing on the current CPU */
      matLastParticle_IndexOnCPU = matSwarm->particleLocalCount - 1; /* lastParticle_I is the last particle in the array of particles residing on the current CPU */

      /* we are switching places between the particle we are deleting and the last particle in the array
	 then actually deleting the last particle in the array */
      intSwarmLastParticle = (IntegrationPoint*)Swarm_ParticleAt( intSwarm, intLastParticle_IndexOnCPU );
      matSwarmLastParticle = (MaterialPoint*)Swarm_ParticleAt( matSwarm, matLastParticle_IndexOnCPU );
      /* we have the last particles of both swarms now..
	 if the intSwarmLastParticle material ref does NOT point to the matSwarmLastParticle
	 then we might be in trouble...so test for this here */
      refTestLastParticle = OneToOneMapper_GetMaterialRef(intSwarm->mapper, intSwarmLastParticle); /* so we can get the material point */
      //refToLastMatParticleFlag = 0;
      if(refTestLastParticle->particle_I != matLastParticle_IndexOnCPU){
	    printf("\e[31mThe last int particle does NOT point to the last mat particle..we need to handle this\n\n");
	    printf("\e[0;32m");
	    refToLastMatParticle = getIntParticleMaterialRef_PointingToMaterialParticle( intSwarm, matLastParticle_IndexOnCPU );
	    refToLastMatParticle->particle_I = matParticleToRemove_IndexOnCPU;
	    refToLastMatParticleFlag = 1;
	    exit(0);
	    /* if the test here is true then that means that some other int particle points to the last material particle
               and therefore *its* material ref needs to be updated when we move the last material particle.
	       It also means that we do not update the ref for the last int particle as that must be pointing somewhere
	       other than the last material particle and the mat particle that is somewhere else is not being moved.
	    */
      }
      if ( intParticleToRemove_IndexOnCPU != intLastParticle_IndexOnCPU ) {                                                                                                                       
	    //printf("Deleting int particle number %d from cell %d local cell particle num = %d\n\n", intParticleToRemove_IndexOnCPU ,lCell_I, intParticleToRemove_IndexWithinCell);                    
                                                                                                                                                                                                
	    intLastParticle_IndexWithinCell = Swarm_GetParticleIndexWithinCell( intSwarm, intSwarmLastParticle->owningCell, intLastParticle_IndexOnCPU);                                              
                                                                                                                                                                                                
	    /* particleToRemove gets over-written by lastParticle  void *memcpy(void *dest, const void *src, size_t n); */                                                              
	    memcpy( intParticleToRemove, intSwarmLastParticle, intparticleSize );                                                                                                       
                                                                                                                                                                                                
	    /* Change value in cell particle table to point to new index in localCPU particle array.                                                                                    
	       ok.. lastParticle_I is a global index of the last particle...we no longer have that many particles                                                                       
	       so give it a new index in the swarm...it has now assumed the identity of the particle that got removed...complete identity theft*/                                       
	    intSwarm->cellParticleTbl[ intSwarmLastParticle->owningCell ] [intLastParticle_IndexWithinCell ] = intParticleToRemove_IndexOnCPU;
	    if(!refToLastMatParticleFlag){                                                         
		  ref_moved_particle = OneToOneMapper_GetMaterialRef(intSwarm->mapper, intParticleToRemove); /* so we can get the material point */
		  ref_moved_particle->particle_I = matParticleToRemove_IndexOnCPU;
	    }
      }                                                                                                                                                                                 
      intSwarm->particleLocalCount--; 
      //Swarm_Realloc( intSwarm ); /* do I really want to do this EVERY time I delete a particle? No I don't hmmm*/
      if ( matParticleToRemove_IndexOnCPU != matLastParticle_IndexOnCPU ) {                                                                                                                       
	    //printf("Deleting mat particle number %d from cell %d local cell particle num = %d\n\n", matParticleToRemove_IndexOnCPU, matParticleToRemove->owningCell, matParticleToRemove_IndexWithinCell);                    
                                                                                                                                                                                                
	    matLastParticle_IndexWithinCell = Swarm_GetParticleIndexWithinCell( matSwarm, matSwarmLastParticle->owningCell, matLastParticle_IndexOnCPU);                                              
                                                                                                                                                                                                
	    /* particleToRemove gets over-written by lastParticle  void *memcpy(void *dest, const void *src, size_t n); */                                                              
	    memcpy( matParticleToRemove, matSwarmLastParticle, matparticleSize );                                                                                                       
                                                                                                                                                                                                
	    /* Change value in cell particle table to point to new index in localCPU particle array.                                                                                    
	       ok.. lastParticle_I is a global index of the last particle...we no longer have that many particles                                                                       
	       so give it a new index in the swarm...it has now assumed the identity of the particle that got removed...complete identity theft*/                                       
	    matSwarm->cellParticleTbl[ matSwarmLastParticle->owningCell][matLastParticle_IndexWithinCell ] = matParticleToRemove_IndexOnCPU;                                                          
                                                                                                                                                                                                
      }                                                                                                                                                                                 
      matSwarm->particleLocalCount--; 
      //Swarm_Realloc( intSwarm ); /* do I really want to do this EVERY time I delete a particle? No I don't hmmm*/
	      
}
/****************************************************************************************************
 This function deletes integration particle number intParticleToRemove_IndexOnCPU
 on integration swarm as well as the corresponding material point.
 It assumes the one to one coincidence mapper is being used. i.e. for every integration point 
 there is one material point and vice versa.
******************************************************************************************************/
void deleteIntParticleByIndexOnCPU( IntegrationPointsSwarm*  intSwarm, MaterialPointsSwarm* matSwarm, Particle_Index intParticleToRemove_IndexOnCPU ) {
      MaterialPointRef*       ref;
      MaterialPointRef*       refTestLastParticle;
      MaterialPointRef*       ref_moved_particle;
      MaterialPointRef*   refToLastMatParticle;
      int refToLastMatParticleFlag = 0;
      //Particle_Index    intParticleToRemove_IndexOnCPU;/* the particle number within the swarm on the local CPU */
      Particle_Index    matParticleToRemove_IndexOnCPU;
      IntegrationPoint* intParticleToRemove;
      MaterialPoint*    matParticleToRemove;
      Particle_Index    intParticleToRemove_IndexWithinCell;/* the number of the particle within the cell */
      Particle_Index    matParticleToRemove_IndexWithinCell;
      IntegrationPoint* intSwarmLastParticle; /* The last particle in the swarm on the local CPU */
      MaterialPoint*    matSwarmLastParticle;
      Particle_Index    intLastParticle_IndexOnCPU;
      Particle_Index    matLastParticle_IndexOnCPU;
      Particle_Index    intLastParticle_IndexWithinCell;
      Particle_Index    matLastParticle_IndexWithinCell;
      Cell_LocalIndex lCell_I;


      SizeT                 intparticleSize     = intSwarm->particleExtensionMgr->finalSize;
      SizeT                 matparticleSize     = matSwarm->particleExtensionMgr->finalSize;

      //intParticleToRemove_IndexOnCPU  = intSwarm->cellParticleTbl[ lCell_I ][ intParticleToRemove_IndexWithinCell ];
      intParticleToRemove =  (IntegrationPoint*)Swarm_ParticleAt( intSwarm, intParticleToRemove_IndexOnCPU);

      intParticleToRemove_IndexWithinCell = Swarm_GetParticleIndexWithinCell( intSwarm, intParticleToRemove->owningCell, intParticleToRemove_IndexOnCPU);

      lCell_I = intParticleToRemove->owningCell;

      ref = OneToOneMapper_GetMaterialRef(intSwarm->mapper, intParticleToRemove); /* so we can get the material point */
      matParticleToRemove_IndexOnCPU = ref->particle_I;
      matParticleToRemove =  (MaterialPoint*)Swarm_ParticleAt( matSwarm, matParticleToRemove_IndexOnCPU);
      matParticleToRemove_IndexWithinCell = Swarm_GetParticleIndexWithinCell( matSwarm, matParticleToRemove->owningCell, matParticleToRemove_IndexOnCPU);

      /* The particles are homeless after these function-calls and for this reason alone must be destroyed */
      Swarm_RemoveParticleFromCell( intSwarm, lCell_I, intParticleToRemove_IndexWithinCell );
      Swarm_RemoveParticleFromCell( matSwarm, matParticleToRemove->owningCell, matParticleToRemove_IndexWithinCell );

      /* Copy over particle to remove with the last particle in the particle array - as long as it isn't the last one  */
      intLastParticle_IndexOnCPU = intSwarm->particleLocalCount - 1; /* lastParticle_I is the last particle in the array of particles residing on the current CPU */
      matLastParticle_IndexOnCPU = matSwarm->particleLocalCount - 1; /* lastParticle_I is the last particle in the array of particles residing on the current CPU */

      /* we are switching places between the particle we are deleting and the last particle in the array
	 then actually deleting the last particle in the array */
      intSwarmLastParticle = (IntegrationPoint*)Swarm_ParticleAt( intSwarm, intLastParticle_IndexOnCPU );
      matSwarmLastParticle = (MaterialPoint*)Swarm_ParticleAt( matSwarm, matLastParticle_IndexOnCPU );
      /* we have the last particles of both swarms now..
	 if the intSwarmLastParticle material ref does NOT point to the matSwarmLastParticle
	 then we might be in trouble...so test for this here */
      refTestLastParticle = OneToOneMapper_GetMaterialRef(intSwarm->mapper, intSwarmLastParticle); /* so we can get the material point */
      //refToLastMatParticleFlag = 0;
      if(refTestLastParticle->particle_I != matLastParticle_IndexOnCPU){
	    printf("\e[31mThe last int particle does NOT point to the last mat particle..we need to handle this\n\n");
	    printf("\e[0;32m");
	    refToLastMatParticle = getIntParticleMaterialRef_PointingToMaterialParticle( intSwarm, matLastParticle_IndexOnCPU );
	    refToLastMatParticle->particle_I = matParticleToRemove_IndexOnCPU;
	    refToLastMatParticleFlag = 1;
	    exit(0);
	    /* if the test here is true then that means that some other int particle points to the last material particle
               and therefore *its* material ref needs to be updated when we move the last material particle.
	       It also means that we do not update the ref for the last int particle as that must be pointing somewhere
	       other than the last material particle and the mat particle that is somewhere else is not being moved.
	    */
      }

      if ( intParticleToRemove_IndexOnCPU != intLastParticle_IndexOnCPU ) {
	    //printf("Deleting int particle number %d from cell %d local cell particle num = %d\n\n", intParticleToRemove_IndexOnCPU ,lCell_I, intParticleToRemove_IndexWithinCell);

	    intLastParticle_IndexWithinCell = Swarm_GetParticleIndexWithinCell( intSwarm, intSwarmLastParticle->owningCell, intLastParticle_IndexOnCPU);
      
	    /* particleToRemove gets over-written by lastParticle  void *memcpy(void *dest, const void *src, size_t n); */                                                              
	    memcpy( intParticleToRemove, intSwarmLastParticle, intparticleSize );

	    /* Change value in cell particle table to point to new index in localCPU particle array.
	       ok.. lastParticle_I is a global index of the last particle...we no longer have that many particles                                                                       
	       so give it a new index in the swarm...it has now assumed the identity of the particle that got removed...complete identity theft*/                                       
	    intSwarm->cellParticleTbl[ intSwarmLastParticle->owningCell ] [intLastParticle_IndexWithinCell ] = intParticleToRemove_IndexOnCPU;
	    if(!refToLastMatParticleFlag){
		  ref_moved_particle = OneToOneMapper_GetMaterialRef(intSwarm->mapper, intParticleToRemove); /* so we can get the material point */
		  ref_moved_particle->particle_I = matParticleToRemove_IndexOnCPU;
	    }
      }

      intSwarm->particleLocalCount--; 
      //Swarm_Realloc( intSwarm ); /* do I really want to do this EVERY time I delete a particle? No I don't hmmm*/
      if ( matParticleToRemove_IndexOnCPU != matLastParticle_IndexOnCPU ) {                                                                                                                       
	    //printf("Deleting mat particle number %d from cell %d local cell particle num = %d\n\n", matParticleToRemove_IndexOnCPU, matParticleToRemove->owningCell, matParticleToRemove_IndexWithinCell);                    
                                                                                                                                                                                                
	    matLastParticle_IndexWithinCell = Swarm_GetParticleIndexWithinCell( matSwarm, matSwarmLastParticle->owningCell, matLastParticle_IndexOnCPU);                                              
                                                                                                                                                                                                
	    /* particleToRemove gets over-written by lastParticle  void *memcpy(void *dest, const void *src, size_t n); */                                                              
	    memcpy( matParticleToRemove, matSwarmLastParticle, matparticleSize );                                                                                                       
                                                                                                                                                                                                
	    /* Change value in cell particle table to point to new index in localCPU particle array.                                                                                    
	       ok.. lastParticle_I is a global index of the last particle...we no longer have that many particles                                                                       
	       so give it a new index in the swarm...it has now assumed the identity of the particle that got removed...complete identity theft*/                                       
	    matSwarm->cellParticleTbl[ matSwarmLastParticle->owningCell][matLastParticle_IndexWithinCell ] = matParticleToRemove_IndexOnCPU;                                                          
                                                                                                                                                                                                
      }                                                                                                                                                                                 
      matSwarm->particleLocalCount--; 
      //Swarm_Realloc( intSwarm ); /* do I really want to do this EVERY time I delete a particle? No I don't hmmm*/
	      
}
void splitIntParticleByIndexWithinCell( IntegrationPointsSwarm*  intSwarm, MaterialPointsSwarm* matSwarm, Cell_LocalIndex lCell_I, Particle_Index intParticleToSplit_IndexWithinCell, Coord xi ) {
      MaterialPointRef*       ref;
      Particle_Index    intNewParticle_IndexOnCPU;/* the particle number within the swarm on the local CPU */
      Particle_Index    matNewParticle_IndexOnCPU;
      IntegrationPoint* intNewParticle;
      MaterialPoint*    matNewParticle;
      Particle_Index    intNewParticle_IndexWithinCell;/* the number of the particle within the cell */
      //Particle_Index    matNewParticle_IndexWithinCell;
      IntegrationPoint* intParticleToSplit;
      MaterialPoint*    matParticleToSplit;
      //Particle_Index    intParticleToSplit_IndexWithinCell;
      //Particle_Index    matParticleToSplit_IndexWithinCell;
      Particle_Index    matParticleToSplit_IndexOnCPU;
//      SizeT                 intparticleSize     = intSwarm->particleExtensionMgr->finalSize;
//      SizeT                 matparticleSize     = matSwarm->particleExtensionMgr->finalSize;
      Coord                   newCoord;
//      Coord                   xi;

      FiniteElement_Mesh*     mesh              = (FiniteElement_Mesh*)((ElementCellLayout*)matSwarm->cellLayout)->mesh;

	      
      //intParticleToSplit_IndexWithinCell = maxI;

      /* Add a new particle to end the end of each swarm */
      /* this calls Swarm_Realloc -- don't like reallocing every time we create a particle
	 need to do this differently */
      intNewParticle     = (IntegrationPoint*) Swarm_CreateNewParticle( intSwarm, &intNewParticle_IndexOnCPU );
      matNewParticle     = (MaterialPoint*) Swarm_CreateNewParticle( matSwarm, &matNewParticle_IndexOnCPU );
	      
      /* Copy particle information */
      intParticleToSplit = (IntegrationPoint*) Swarm_ParticleInCellAt( intSwarm, lCell_I, intParticleToSplit_IndexWithinCell );
      memcpy( intNewParticle, intParticleToSplit, intSwarm->particleExtensionMgr->finalSize );
      Swarm_AddParticleToCell( intSwarm, lCell_I, intNewParticle_IndexOnCPU );
      ref = OneToOneMapper_GetMaterialRef(intSwarm->mapper, intNewParticle); /* so we can set the reference to the material point */
      ref->particle_I = matNewParticle_IndexOnCPU; /* now the ref for the new int particle points to the new material particle -- the swarm id should be correct because we memcpy'd from the original */

      /* Now get the material point corresponding to the int point being split so that
	 we can copy the material properties to the newly created material particle */
      ref = OneToOneMapper_GetMaterialRef(intSwarm->mapper, intParticleToSplit);
      matParticleToSplit_IndexOnCPU = ref->particle_I;
      matParticleToSplit  =  (MaterialPoint*)Swarm_ParticleAt( matSwarm, matParticleToSplit_IndexOnCPU );
      memcpy( matNewParticle, matParticleToSplit, matSwarm->particleExtensionMgr->finalSize );
      Swarm_AddParticleToCell( matSwarm, matNewParticle->owningCell, matNewParticle_IndexOnCPU );

      /* Copy new local position to xi on new int particle */
      memcpy( intNewParticle->xi, xi, sizeof(Coord) );
	      
      /* Get new Global Coordinates from the Local Coordinates */
      FiniteElement_Mesh_CalcGlobalCoordFromLocalCoord( mesh, matSwarm->dim, lCell_I, xi, newCoord );

      /* Copy new global position to coord on new mat particle */
      memcpy( matNewParticle->coord, newCoord, sizeof(Coord) );
	      
      intNewParticle_IndexWithinCell = Swarm_GetParticleIndexWithinCell( intSwarm, lCell_I, intNewParticle_IndexOnCPU);
      /*
	printf("\e[1;33m");
	printf("Creating int particle number %d from cell %d local within cell particle num = %d\n\n", intNewParticle_IndexOnCPU ,lCell_I, intNewParticle_IndexWithinCell);
	printf("matSwarm particle count is now : %d\n",matSwarm->cellParticleCountTbl[lCell_I]);
	printf("intSwarm particle count is now : %d\n",intSwarm->cellParticleCountTbl[lCell_I]);
	printf("Original Population was %d\n",nump-1);
	printf("\e[0;32m");
      */

}
void splitIntParticleByIndexOnCPU( IntegrationPointsSwarm*  intSwarm, MaterialPointsSwarm* matSwarm, Particle_Index intParticleToSplit_IndexOnCPU, Coord xi ) {
      MaterialPointRef*       ref;
      Particle_Index    intNewParticle_IndexOnCPU;/* the particle number within the swarm on the local CPU */
      Particle_Index    matNewParticle_IndexOnCPU;
      IntegrationPoint* intNewParticle;
      MaterialPoint*    matNewParticle;
      Particle_Index    intNewParticle_IndexWithinCell;/* the number of the particle within the cell */
      //Particle_Index    matNewParticle_IndexWithinCell;
      IntegrationPoint* intParticleToSplit;
      MaterialPoint*    matParticleToSplit;
      Particle_Index    intParticleToSplit_IndexWithinCell;
      //Particle_Index    matParticleToSplit_IndexWithinCell;
      Particle_Index    matParticleToSplit_IndexOnCPU;
//      SizeT                 intparticleSize     = intSwarm->particleExtensionMgr->finalSize;
//      SizeT                 matparticleSize     = matSwarm->particleExtensionMgr->finalSize;
      Coord                   newCoord;
//      Coord                   xi;

      Cell_LocalIndex lCell_I;

      FiniteElement_Mesh*     mesh              = (FiniteElement_Mesh*)((ElementCellLayout*)matSwarm->cellLayout)->mesh;

	      
      //intParticleToSplit_IndexWithinCell = maxI;

      /* Add a new particle to end the end of each swarm */
      /* this calls Swarm_Realloc -- don't like reallocing every time we create a particle
	 need to do this differently */
      intNewParticle     = (IntegrationPoint*) Swarm_CreateNewParticle( intSwarm, &intNewParticle_IndexOnCPU );
      matNewParticle     = (MaterialPoint*) Swarm_CreateNewParticle( matSwarm, &matNewParticle_IndexOnCPU );
	      


      /* Copy particle information */
      intParticleToSplit_IndexWithinCell = Swarm_GetParticleIndexWithinCell( intSwarm, intParticleToSplit->owningCell, intParticleToSplit_IndexOnCPU);

      intParticleToSplit = (IntegrationPoint*) Swarm_ParticleInCellAt( intSwarm, lCell_I, intParticleToSplit_IndexWithinCell );


      lCell_I = intParticleToSplit->owningCell;


      memcpy( intNewParticle, intParticleToSplit, intSwarm->particleExtensionMgr->finalSize );
      Swarm_AddParticleToCell( intSwarm, lCell_I, intNewParticle_IndexOnCPU );
      ref = OneToOneMapper_GetMaterialRef(intSwarm->mapper, intNewParticle); /* so we can set the reference to the material point */
      ref->particle_I = matNewParticle_IndexOnCPU; /* now the ref for the new int particle points to the new material particle -- the swarm id should be correct because we memcpy'd from the original */

      /* Now get the material point corresponding to the int point being split so that
	 we can copy the material properties to the newly created material particle */
      ref = OneToOneMapper_GetMaterialRef(intSwarm->mapper, intParticleToSplit);
      matParticleToSplit_IndexOnCPU = ref->particle_I;
      matParticleToSplit  =  (MaterialPoint*)Swarm_ParticleAt( matSwarm, matParticleToSplit_IndexOnCPU );
      memcpy( matNewParticle, matParticleToSplit, matSwarm->particleExtensionMgr->finalSize );
      Swarm_AddParticleToCell( matSwarm, matNewParticle->owningCell, matNewParticle_IndexOnCPU );

      /* Copy new local position to xi on new int particle */
      memcpy( intNewParticle->xi, xi, sizeof(Coord) );
	      
      /* Get new Global Coordinates from the Local Coordinates */
      FiniteElement_Mesh_CalcGlobalCoordFromLocalCoord( mesh, matSwarm->dim, lCell_I, xi, newCoord );

      /* Copy new global position to coord on new mat particle */
      memcpy( matNewParticle->coord, newCoord, sizeof(Coord) );
	      
      intNewParticle_IndexWithinCell = Swarm_GetParticleIndexWithinCell( intSwarm, lCell_I, intNewParticle_IndexOnCPU);
      /*
	printf("\e[1;33m");
	printf("Creating int particle number %d from cell %d local within cell particle num = %d\n\n", intNewParticle_IndexOnCPU ,lCell_I, intNewParticle_IndexWithinCell);
	printf("matSwarm particle count is now : %d\n",matSwarm->cellParticleCountTbl[lCell_I]);
	printf("intSwarm particle count is now : %d\n",intSwarm->cellParticleCountTbl[lCell_I]);
	printf("Original Population was %d\n",nump-1);
	printf("\e[0;32m");
      */

}
/* sort indexOnCPU in reverse order */
int compare_indexOnCPU(const void * _particleA, const void * _particleB){
      struct deleteParticle * particleA = (struct deleteParticle *) _particleA;
      struct deleteParticle * particleB = (struct deleteParticle *) _particleB;
      
      if(particleA->indexOnCPU < particleB->indexOnCPU)
	    return 1;
      else
	    return -1;

}
/* Calculate the integration weights for each particle by contructing
   a voronoi diagram in an element in 3D*/
void _PCDVC_Calculate3D( void* pcdvc, void* _swarm, Cell_LocalIndex lCell_I ) {
	PCDVC*             self            = (PCDVC*)  pcdvc;
	IntegrationPointsSwarm*  intSwarm  = (IntegrationPointsSwarm*) _swarm;
        MaterialPointsSwarm* matSwarm =	(MaterialPointsSwarm*) self->materialPointsSwarm;
	/* CoincidentMapper is a special case of the one to one mapper */
	//CoincidentMapper* mapper  = (CoincidentMapper*)(intSwarm->mapper); /* need the mapper after-all to update the material ref */
	Particle_InCellIndex         cParticleCount;
	IntegrationPoint**           particle;
	static int visited = 0 ;
	//static int deleted = 0 ;
	double dx,dy,dz,da;
	static struct cell *cells;// the connected grid
	struct particle *pList;// particle List
	struct chain *bchain;//boundary chain
	int nump_orig,nump,numx,numy,numz;
	double BBXMIN = -1.0; // the ranges of the local coordinates of a FEM cell.
	double BBXMAX = 1.0;
	double BBYMIN = -1.0;
	double BBYMAX = 1.0;
	double BBZMIN = -1.0;
	double BBZMAX = 1.0;
	int i,k;

	/*************************************/
	/* stuff for particle removal/adding */
	struct deleteParticle* deleteList;
	double maxW,minW;
	int maxI, minI;
	double lowT = self->lowerT;
	double upT = self->upperT;
	int delete_flag, split_flag;
	Particle_Index  *splitList;
//	int deleteListStackPtr = -1;/* use a number to tell me how many particles we are going to delete: saves doing it some other less efficient way */
//	int splitListStackPtr = -1;
	int maxDeletions = self->maxDeletions;/* hard setting this till I get stuff from xml file */
	int maxSplits = self->maxSplits;
	int splitCount;
	int deleteCount;
	int Count;
	int matTypeFirst;
	int matType;
	Bool splitInInterfaceCells = self->splitInInterfaceCells;

//	SizeT                 intparticleSize     = intSwarm->particleExtensionMgr->finalSize;
//	SizeT                 matparticleSize     = matSwarm->particleExtensionMgr->finalSize;
//	Coord                   newCoord;
	Coord                   xi;

//	FiniteElement_Mesh*     mesh              = (FiniteElement_Mesh*)((ElementCellLayout*)matSwarm->cellLayout)->mesh;

	/* end decs needed for particle control */
	/*************************************/
	


	numx = self->resX;
	numy = self->resY;
	numz = self->resZ;

	nump_orig = nump = cParticleCount = intSwarm->cellParticleCountTbl[lCell_I];

	/* need a struct for the deletList because we must sort it bu indexOnCPU and delete in reverse order
           so we don't have the potential problem of  deleting a particle from the list that points to the last particle on the swarm */
	deleteList = (struct deleteParticle*)malloc(nump*sizeof(struct deleteParticle));/* I don't think I am going to let you delete more than half the particles in a given cell */
	splitList  = (Particle_Index*)malloc(nump*sizeof(Particle_Index));


	Journal_Firewall( nump , Journal_Register(Error_Type, "PCDVC"), "Error in %s: Problem has an under resolved cell (Cell Id = %d), add more particles to your model\n", __func__, lCell_I );

	dx = (BBXMAX - BBXMIN)/numx;
	dy = (BBYMAX - BBYMIN)/numy;
	dz = (BBZMAX - BBZMIN)/numz;
	da = dx*dy*dz;
	
	// Construct the grid for the Voronoi cells only once.
	// If we wanted to call this function again during a job with a different resolution
	// then we should destroy the grid once we have looped through the whole mesh.
	// I am assuming we are not going to do that for now.
	// Easy to implement this anyway, if needed.
	if(!visited){
	      /* The PCDVC class should really be a class the next level up here */
	      /* We should be able to swap out the WeightsCalculator_CalculateAll instead of just setting
                 a pointer inside that function */
	      visited++;
	      _DVCWeights_ConstructGrid(&cells,numz,numy,numx,BBXMIN,BBYMIN,BBZMIN,BBXMAX,BBYMAX,BBZMAX);
	}
	
	// init the data structures
	_DVCWeights_InitialiseStructs( &bchain, &pList, nump);
	_DVCWeights_ResetGrid(&cells,numz*numy*numx);
	
	particle = (IntegrationPoint**)malloc( (nump)*sizeof(IntegrationPoint*));
	
        // initialize the particle positions to be the local coordinates of the material swarm particles
	// I am assuming the xi's (local coords) are precalculated somewhere and get reset based on material
	// positions each time step.
	for(i=0;i<nump;i++){
	      
	      particle[i] = (IntegrationPoint*) Swarm_ParticleInCellAt( intSwarm, lCell_I, i );
	      pList[i].x = particle[i]->xi[0];
	      pList[i].y = particle[i]->xi[1];
	      pList[i].z = particle[i]->xi[2];
	      
	}
	_DVCWeights_CreateVoronoi( &bchain, &pList, &cells, dx, dy, dz, nump, numx, numy, numz, BBXMIN, BBXMAX, BBYMIN, BBYMAX, BBZMIN, BBZMAX);
	_DVCWeights_GetCentroids( cells, pList,numz,numy,numx,nump,da);

	/****************************/
	/****************************/
	/* Start Population Control */
	/****************************/
	/****************************/
	/* todo: put print statements into Journal statements */
	split_flag = 0;
	delete_flag = 0;
//	maxW = 0.0;
//	minW = 8.0;
	splitCount = 0;
	deleteCount = 0;
	/* shouldn't need maxI and minI now */
	maxW = upT*8/100.0;
	minW = lowT*8/100.0;
	/* check to see if we are in an interface cell.
           We never want to delete particles in an interface cell */
	matTypeFirst = IntegrationPointMapper_GetMaterialIndexAt(intSwarm->mapper,intSwarm->cellParticleTbl[ lCell_I ][ 0 ]);
	for(i=0;i<nump;i++){
	      matType = IntegrationPointMapper_GetMaterialIndexAt(intSwarm->mapper,intSwarm->cellParticleTbl[ lCell_I ][ i ]);
	      if(matType != matTypeFirst){
		    maxDeletions = 0; /* no deletions in an interface cell */
		    /* this may be inadequate...we may need to do something in the neighbouring cells to interface cells as well */
		    if(!splitInInterfaceCells){
			  maxSplits = 0;
		    }
		    break;
	      }
	}
	for(i=0;i<nump;i++){
	      if(pList[i].w > maxW){ /* maxW = pList[i].w; maxI = i;*/ splitList[splitCount] = i; splitCount++;}
	      if(pList[i].w < minW){
		    /* minW = pList[i].w; minI = i; */
		    deleteList[deleteCount].indexWithinCell = i;
		    deleteList[deleteCount].indexOnCPU  = intSwarm->cellParticleTbl[ lCell_I ][ i ];		    
		    deleteCount++;
	      }
	}
	/* sort the deleteList by indexOnCPU so we can delete the list in reverse order */
	qsort(deleteList, (deleteCount), sizeof(struct deleteParticle),compare_indexOnCPU);
	//deleteCount--; /* is going to be one size too large after the loop */
	/*
	for(i=0;i<deleteCount;i++){
	      printf("deleteCount = %d\n",deleteCount);
	      printf("indices are indexWithinCell %d indexOnCPU %d\n",deleteList[i].indexWithinCell,deleteList[i].indexOnCPU);
	}
	*/
	if(maxDeletions > nump-4){ maxDeletions = nump/2;}

	/* we now have our lists of particles to delete and split */

//	if(pList[maxI].w > upT*8/100.0){
	Count = maxSplits > splitCount ? splitCount : maxSplits;
	for(i=0;i<Count;i++){
	      maxI = splitList[i];
	      /* now get local coords from centroid of the cell that particleToSplit lives in */
	      xi[0] = pList[maxI].cx;
	      xi[1] = pList[maxI].cy;
	      xi[2] = pList[maxI].cz;

	      split_flag = 1;
	      nump++;

	      splitIntParticleByIndexWithinCell( intSwarm, matSwarm, lCell_I, maxI, xi );

	}

//	if( (pList[minI].w < lowT*8/100.0) ){
	Count = maxDeletions > deleteCount ? deleteCount : maxDeletions;
	for(i=0;i<Count;i++){

	      minI = deleteList[i].indexOnCPU;

	      deleteIntParticleByIndexOnCPU( intSwarm,  matSwarm, minI );

	      delete_flag = 1;
	      nump--;

	} /* if(pList[minI].w < lowT*8/100.0) */
	
	//printf("pList[maxI].w = %lf particle num = %d : %d\n", pList[maxI].w, pList[maxI].index,maxI);
	//printf("pList[minI].w = %lf particle num = %d : %d\n", pList[minI].w, pList[minI].index,minI);
	if(delete_flag || split_flag ){/* then we need to redo the Voronoi diagram */
	      for(k=0;k<nump_orig;k++){
		    free(bchain[k].new_claimed_cells);
		    free(bchain[k].new_bound_cells);
	      }
	      free(particle);
	      free(bchain);
	      free(pList);
	      if(nump < 3){
		    Journal_Firewall( 0 , Journal_Register(Error_Type, "PCDVC"), "Something went horribly wrong in %s: Problem has an under resolved cell (Cell Id = %d), check or tune your population control parameters\n", __func__, lCell_I );
	      }
	      // init the data structures
	      _DVCWeights_InitialiseStructs( &bchain, &pList, nump);
	      //_DVCWeights_ResetGrid(&cells,numz*numy*numx);
	      
	      particle = (IntegrationPoint**)malloc( (nump)*sizeof(IntegrationPoint*));
	
	      // re-initialize the particle positions to be the local coordinates of the material swarm particles
	      for(i=0;i<nump;i++){
		    
		    particle[i] = (IntegrationPoint*) Swarm_ParticleInCellAt( intSwarm, lCell_I, i );
		    pList[i].x = particle[i]->xi[0];
		    pList[i].y = particle[i]->xi[1];
		    pList[i].z = particle[i]->xi[2];
		    //pList[i].index = i; // to track which particle numbers we have after sorting this list */
		    
	      }
	      //printf("Population of matSwarm is %d\n",matSwarm->particleLocalCount);
	      //printf("Population of intSwarm is %d\n",intSwarm->particleLocalCount);
	      _DVCWeights_ResetGrid(&cells,numz*numy*numx);
	      //reset_grid(&cells,numz*numy*numx);/* adding this line fixed memory leak probs */
	      //create_voronoi( &bchain, &pList, &cells, dx, dy, dz, nump, numx, numy, numz, BBXMIN, BBXMAX, BBYMIN, BBYMAX, BBZMIN, BBZMAX);
	      //get_centroids( cells, pList,numz,numy,numx,nump,da);
	      _DVCWeights_CreateVoronoi( &bchain, &pList, &cells, dx, dy, dz, nump, numx, numy, numz, BBXMIN, BBXMAX, BBYMIN, BBYMAX, BBZMIN, BBZMAX);
	      _DVCWeights_GetCentroids( cells, pList,numz,numy,numx,nump,da);

	}/* if delete_flag */
	/****************************/
	/****************************/
	/* End Population Control */
	/****************************/
	/****************************/

	// We are setting the integration points to be the centroids of the Voronoi regions here and
	// the weight is the volume of each Voronoi region.
	for(i=0;i<nump;i++){

	      particle[i]->xi[0] = pList[i].cx;
	      particle[i]->xi[1] = pList[i].cy;
	      particle[i]->xi[2] = pList[i].cz;
	      particle[i]->weight = pList[i].w;

	}	
	for(k=0;k<nump;k++){
	      free(bchain[k].new_claimed_cells);
	      free(bchain[k].new_bound_cells);
	}
	free(particle);
	free(bchain);
	free(pList);
	free(deleteList);
	free(splitList);

}

/* Calculate the integration weighting for each particle by contructing
a voronoi diagram in an element in 2D*/
void _PCDVC_Calculate2D( void* pcdvc, void* _swarm, Cell_LocalIndex lCell_I ) {
	PCDVC*             self            = (PCDVC*)  pcdvc;
	IntegrationPointsSwarm*  intSwarm  = (IntegrationPointsSwarm*) _swarm;
        MaterialPointsSwarm* matSwarm =	(MaterialPointsSwarm*) self->materialPointsSwarm;
	/* CoincidentMapper is a special case of the one to one mapper */
	//CoincidentMapper* mapper  = (CoincidentMapper*)(intSwarm->mapper); /* need the mapper after-all to update the material ref */
	Particle_InCellIndex         cParticleCount;
	IntegrationPoint**           particle;
	static int visited = 0 ;
	//static int deleted = 0 ;

	double dx,dy,da;
	static struct cell2d *cells;// the connected grid
	struct particle2d *pList;// particle List
	struct chain *bchain;//boundary chain
	int nump_orig,nump,numx,numy;
	double BBXMIN = -1.0; // the ranges of the local coordinates of a FEM cell.
	double BBXMAX = 1.0;
	double BBYMIN = -1.0;
	double BBYMAX = 1.0;
	int i,k;

	/*************************************/
	/* stuff for particle removal/adding */
	double maxW,minW;
	int maxI, minI;
	double lowT = self->lowerT;
	double upT = self->upperT;
	int delete_flag, split_flag;
	struct deleteParticle* deleteList;
	Particle_Index  *splitList;
	int splitCount;
	int deleteCount;
	int maxDeletions = self->maxDeletions;/* hard setting this till I get stuff from xml file */
	int maxSplits = self->maxSplits;
	int Count;
	int matTypeFirst;
	int matType;
	Bool splitInInterfaceCells = self->splitInInterfaceCells;
//	SizeT                 intparticleSize     = intSwarm->particleExtensionMgr->finalSize;
//	SizeT                 matparticleSize     = matSwarm->particleExtensionMgr->finalSize;

//	Coord                   newCoord;
	Coord                   xi;

//	FiniteElement_Mesh*     mesh              = (FiniteElement_Mesh*)((ElementCellLayout*)matSwarm->cellLayout)->mesh;

	/* end decs needed for particle control */
	/*************************************/
	

	numx = self->resX;
	numy = self->resY;

	nump_orig = nump = cParticleCount = intSwarm->cellParticleCountTbl[lCell_I];

	/* need a struct for the deletList because we must sort it bu indexOnCPU and delete in reverse order
           so we don't have the potential problem of  deleting a particle from the list that points to the last particle on the swarm */
	deleteList = (struct deleteParticle*)malloc(nump*sizeof(struct deleteParticle));/* I don't think I am going to let you delete more than half the particles in a given cell */
	splitList  = (Particle_Index*)malloc(nump*sizeof(Particle_Index));

	Journal_Firewall( nump , Journal_Register(Error_Type, "PCDVC"), "Error in %s: Problem has an under resolved cell (Cell Id = %d), add more particles to your model\n", __func__, lCell_I );

	dx = (BBXMAX - BBXMIN)/numx;
	dy = (BBYMAX - BBYMIN)/numy;
	da = dx*dy;
	
	// Construct the grid for the Voronoi cells only once.
	// If we wanted to call this function again during a job with a different resolution
	// then we should destroy the grid once we have looped through the whole mesh.
	// I am assuming we are not going to do that for now.
	// Easy to implement this anyway, if needed.
	if(!visited){
	      /* The PCDVC class should really be a class the next level up here */
	      /* We should be able to swap out the WeightsCalculator_CalculateAll instead of just setting
                 a pointer inside that function */
	      visited++;
	      _DVCWeights_ConstructGrid2D(&cells,numy,numx,BBXMIN,BBYMIN,BBXMAX,BBYMAX);
	}
	
	
	// init the data structures
	_DVCWeights_InitialiseStructs2D( &bchain, &pList, nump);
	_DVCWeights_ResetGrid2D(&cells,numy*numx);
	
	particle = (IntegrationPoint**)malloc((nump)*sizeof(IntegrationPoint*));
	
        // initialize the particle positions to be the local coordinates of the material swarm particles
	// I am assuming the xi's (local coords) are precalculated somewhere and get reset based on material
	// positions each time step.
	for(i=0;i<nump;i++){
	      
	      particle[i] = (IntegrationPoint*) Swarm_ParticleInCellAt( intSwarm, lCell_I, i );
	      pList[i].x = particle[i]->xi[0];
	      pList[i].y = particle[i]->xi[1];
	      
	}
	_DVCWeights_CreateVoronoi2D( &bchain, &pList, &cells, dx, dy, nump, numx, numy, BBXMIN, BBXMAX, BBYMIN, BBYMAX);
	_DVCWeights_GetCentroids2D( cells, pList,numy,numx,nump,da);

	/****************************/
	/****************************/
	/* Start Population Control */
	/****************************/
	/****************************/

	split_flag = 0;
	delete_flag = 0;
	splitCount = 0;
	deleteCount = 0;
	/* shouldn't need maxI and minI now */
	maxW = upT*4/100.0;
	minW = lowT*4/100.0;
	/* check to see if we are in an interface cell.
           We never want to delete particles in an interface cell */
	matTypeFirst = IntegrationPointMapper_GetMaterialIndexAt(intSwarm->mapper,intSwarm->cellParticleTbl[ lCell_I ][ 0 ]);
	for(i=0;i<nump;i++){
	      matType = IntegrationPointMapper_GetMaterialIndexAt(intSwarm->mapper,intSwarm->cellParticleTbl[ lCell_I ][ i ]);
	      if(matType != matTypeFirst){
		    maxDeletions = 0; /* no deletions in an interface cell */
		    /* this may be inadequate...we may need to do something in the neighbouring cells to interface cells as well */
		    if(!splitInInterfaceCells){
			  maxSplits = 0;
		    }
		    //printf("------- FOUND an Interface Cell!! --------------\n");
		    break;
	      }
	}
	for(i=0;i<nump;i++){
	      if(pList[i].w > maxW){ /* maxW = pList[i].w; maxI = i;*/ splitList[splitCount] = i; splitCount++;}
	      if(pList[i].w < minW){
		    /* minW = pList[i].w; minI = i; */
		    deleteList[deleteCount].indexWithinCell = i;
		    deleteList[deleteCount].indexOnCPU  = intSwarm->cellParticleTbl[ lCell_I ][ i ];		    
		    deleteCount++;
	      }
	}
	/* sort the deleteList by indexOnCPU so we can delete the list in reverse order */
	qsort(deleteList, (deleteCount), sizeof(struct deleteParticle),compare_indexOnCPU);
	//deleteCount--; /* is going to be one size too large after the loop */
	/*
	for(i=0;i<deleteCount;i++){
	      printf("deleteCount = %d\n",deleteCount);
	      printf("indices are indexWithinCell %d indexOnCPU %d\n",deleteList[i].indexWithinCell,deleteList[i].indexOnCPU);
	}
*/
	if(maxDeletions > nump-4){ maxDeletions = nump/2;}

	/* we now have our lists of particles to delete and split */
	Count = maxSplits > splitCount ? splitCount : maxSplits;
	for(i=0;i<Count;i++){
	      maxI = splitList[i];
	      /* now get local coords from centroid of the cell that particleToSplit lives in */
	      xi[0] = pList[maxI].cx;
	      xi[1] = pList[maxI].cy;
	      //xi[2] = pList[maxI].cz;

	      split_flag = 1;
	      nump++;

	      splitIntParticleByIndexWithinCell( intSwarm, matSwarm, lCell_I, maxI, xi );

	}
	Count = maxDeletions > deleteCount ? deleteCount : maxDeletions;
	for(i=0;i<Count;i++){

	      minI = deleteList[i].indexOnCPU;

	      deleteIntParticleByIndexOnCPU( intSwarm,  matSwarm, minI );


	      delete_flag = 1;
	      nump--;

	      
	} /* if(pList[minI].w < lowT*8/100.0) */
	
	//printf("pList[maxI].w = %lf particle num = %d : %d\n", pList[maxI].w, pList[maxI].index,maxI);
	//printf("pList[minI].w = %lf particle num = %d : %d\n", pList[minI].w, pList[minI].index,minI);
	if(delete_flag || split_flag ){/* then we need to redo the Voronoi diagram */
	      for(k=0;k<nump_orig;k++){
		    free(bchain[k].new_claimed_cells);
		    free(bchain[k].new_bound_cells);
	      }
	      free(particle);
	      free(bchain);
	      free(pList);
	      if(nump < 3){
		    Journal_Firewall( 0 , Journal_Register(Error_Type, "PCDVC"), "Something went horribly wrong in %s: Problem has an under resolved cell (Cell Id = %d), check or tune your population control parameters\n", __func__, lCell_I );
	      }
	      particle = (IntegrationPoint**)malloc((nump)*sizeof(IntegrationPoint*));
              // init the data structures
	      _DVCWeights_InitialiseStructs2D( &bchain, &pList, nump);
	      // re-initialize the particle positions to be the local coordinates of the material swarm particles
	      for(i=0;i<nump;i++){
		    
		    particle[i] = (IntegrationPoint*) Swarm_ParticleInCellAt( intSwarm, lCell_I, i );
		    pList[i].x = particle[i]->xi[0];
		    pList[i].y = particle[i]->xi[1];
		    //pList[i].z = particle[i]->xi[2];
		    //pList[i].index = i; // to track which particle numbers we have after sorting this list */
		    
	      }
	      //printf("Population of matSwarm is %d\n",matSwarm->particleLocalCount);
	      //printf("Population of intSwarm is %d\n",intSwarm->particleLocalCount);

	      _DVCWeights_ResetGrid2D(&cells,numy*numx);
	      //reset_grid(&cells,numz*numy*numx);/* adding this line fixed memory leak probs */
	      //create_voronoi( &bchain, &pList, &cells, dx, dy, dz, nump, numx, numy, numz, BBXMIN, BBXMAX, BBYMIN, BBYMAX, BBZMIN, BBZMAX);
	      //get_centroids( cells, pList,numz,numy,numx,nump,da);
	      _DVCWeights_CreateVoronoi2D( &bchain, &pList, &cells, dx, dy, nump, numx, numy, BBXMIN, BBXMAX, BBYMIN, BBYMAX);
	      _DVCWeights_GetCentroids2D( cells, pList,numy,numx,nump,da);

	}/* if delete_flag */
	/****************************/
	/****************************/
	/* End Population Control */
	/****************************/
	/****************************/


	// We are setting the integration points to be the centroids of the Voronoi regions here and
	// the weight is the volume of each Voronoi region.
	for(i=0;i<nump;i++){

	      particle[i]->xi[0] = pList[i].cx;
	      particle[i]->xi[1] = pList[i].cy;
	      particle[i]->weight = pList[i].w;

	}	
	for(k=0;k<nump;k++){
	      free(bchain[k].new_claimed_cells);
	      free(bchain[k].new_bound_cells);
	}
	free(particle);
	free(bchain);
	free(pList);
	free(deleteList);
	free(splitList);

}

void _PCDVC_Calculate( void* pcdvc, void* _swarm, Cell_LocalIndex lCell_I ){
      Swarm* swarm = (Swarm*) _swarm;
      Dimension_Index dim = swarm->dim;
      Stream*  stream = Journal_Register( Info_Type, swarm->type );
      PCDVC*             self            = (PCDVC*)  pcdvc;
      MaterialPointsSwarm* matSwarm =	(MaterialPointsSwarm*) self->materialPointsSwarm;
      /* it might be nice to report the total deletions and splits as well as the final population here */
      /* One could set the parameters to be too aggressive and cause "swarm thrashing" where many particles
         are being created and destroyed while maintaining some population that it has converged on */
      if(lCell_I == 0){

	   Journal_Printf( stream, "\nOn Proc %d: In func %s(): for swarm \"%s\" Population is %d\n", swarm->myRank, __func__, swarm->name, swarm->particleLocalCount );
	   Journal_Printf( stream, "On Proc %d: In func %s(): for swarm \"%s\" Population is %d\n\n", matSwarm->myRank,__func__, matSwarm->name, matSwarm->particleLocalCount );

      }
      if(dim == 3){
	    _PCDVC_Calculate3D( pcdvc, _swarm, lCell_I);
      }
      else {
	    _PCDVC_Calculate2D( pcdvc, _swarm, lCell_I);
      }

}
/*-------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

