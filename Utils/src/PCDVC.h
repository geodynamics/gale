/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2006, Monash Cluster Computing 
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
**              * Redistributions of source code must retain the above copyright notice, 
**                      this list of conditions and the following disclaimer.
**              * Redistributions in binary form must reproduce the above copyright 
**                      notice, this list of conditions and the following disclaimer in the 
**                      documentation and/or other materials provided with the distribution.
**              * Neither the name of the Monash University nor the names of its contributors 
**                      may be used to endorse or promote products derived from this software 
**                      without specific prior written permission.
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
*%              Louis Moresi - Louis.Moresi@sci.monash.edu.au
*%
** Author:
**              Mirko Velic - Mirko.Velic@sci.monash.edu.au
**              Patrick Sunter - patrick@vpac.org
**              Julian Giordani - julian.giordani@sci.monash.edu.au
**
**  Assumptions:
**       I am assuming that the xi's (local coords) on the IntegrationPoint particles
**       are precalculated somewhere and get reset based on material PIC positions each time step.
**
**  Notes:
**         The PCDVC class should really be a class the next level up here.
**         We should be able to swap out the WeightsCalculator_CalculateAll function instead of just setting
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

  The algorithm here-in computes a discrete voronoi diagram per FEM cell given a set of local
  particle positions, in 3D and 2D. The volumes of the Voronoi regions are used as integration weights for
  the integration point swarm and the integration points are the centroids of the same volumes.

  For a description of this algorithm, see the article by Velic et.al.
     "A Fast Robust Algorithm for computing Discrete Voronoi Diagrams in N-dimensions"



*****************************************************************************************************************/

#ifndef __PICellerator_Weights_PCDVCClass_h__
#define __PICellerator_Weights_PCDVCClass_h__

/* Textual name of this class */
extern const Type PCDVC_Type;

/* PCDVC information */
//   #define __PCDVC 
/* Who's my daddy */ 
//    __DVCWeights 
/* My Data structures */ 
//  MaterialPointsSwarm* materialPointsSwarm; 
//  int upT; 
//  int lowT;

#define __PCDVC                                                         \
    __DVCWeights                                                        \
                                                                        \
    MaterialPointsSwarm*  materialPointsSwarm;                          \
    double                upperT;                                       \
    double                lowerT;                                       \
    Bool                  splitInInterfaceCells;                        \
    Bool                  deleteInInterfaceCells;                       \
    int                   maxDeletions;                                 \
    int                   maxSplits;                                    \
    Bool                  Inflow;                                       \
    double                CentPosRatio;                                 \
    int                   ParticlesPerCell;                             \
    double                Threshold;                                    \
    /* we also need to store some parameters, as everything is turned up      */ \
    /* for interpolation restarts (which may require significant repopulation */ \
    /* and after the first timestep needs to be set back to correct values    */ \
    int                   maxDeletions_orig;                            \
    int                   maxSplits_orig;                               \
    Bool                  Inflow_orig;                                  \
    Bool                  splitInInterfaceCells_orig;                   \
    Bool                  deleteInInterfaceCells_orig;

struct PCDVC { __PCDVC };

struct deleteParticle{
    Particle_Index indexWithinCell;
    Particle_Index indexOnCPU;
};

/*---------------------------------------------------------------------------------------------------------------------
** Constructors
*/



PCDVC* PCDVC_New( Name name, Dimension_Index dim, int* res,
                  MaterialPointsSwarm* mps, double upT, double lowT,
                  int maxDeletions, int maxSplits, Bool splitInInterfaceCells,
                  Bool deleteInInterfaceCells, Bool Inflow, double CentPosRatio,
                  int ParticlesPerCell, double Threshold ) ;


	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define PCDVC_DEFARGS \
                DVCWEIGHTS_DEFARGS

	#define PCDVC_PASSARGS \
                DVCWEIGHTS_PASSARGS

PCDVC* _PCDVC_New(  PCDVC_DEFARGS  );

void _PCDVC_Init( void* pcdvc, MaterialPointsSwarm* mps, double upT, double lowT,
                  int maxDeletions, int maxSplits, Bool splitInInterfaceCells,
                  Bool deleteInInterfaceCells, Bool Inflow, double CentPosRatio,
                  int ParticlesPerCell, double Threshold ) ;


/* Stg_Class_Delete PCDVC implementation */
void _PCDVC_Delete( void* pcdvc );
void _PCDVC_Print( void* pcdvc, Stream* stream );
#define PCDVC_Copy( self )                                      \
    (PCDVC*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
#define PCDVC_DeepCopy( self )                                  \
    (PCDVC*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
void* _PCDVC_Copy( void* pcdvc, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
        
void* _PCDVC_DefaultNew( Name name ) ;

void _PCDVC_AssignFromXML( void* pcdvc, Stg_ComponentFactory* cf, void* data ) ;

void _PCDVC_Build( void* pcdvc, void* data ) ;
void _PCDVC_Initialise( void* pcdvc, void* data ) ;
void _PCDVC_Execute( void* pcdvc, void* data );
MaterialPointRef* getIntParticleMaterialRef_PointingToMaterialParticle( IntegrationPointsSwarm*  intSwarm, Particle_Index matLastParticle_IndexOnCPU );
void splitIntParticleByIndexWithinCell( IntegrationPointsSwarm* intSwarm,  MaterialPointsSwarm* matSwarm, Cell_LocalIndex lCell_I, Particle_Index intParticleToSplit_IndexOnCPU, Coord xi );
void deleteIntParticleByIndexWithinCell( IntegrationPointsSwarm* intSwarm,  MaterialPointsSwarm* matSwarm,  Cell_LocalIndex lCell_I, Particle_Index intParticleToSplit_IndexWithinCell );
void deleteIntParticleByIndexOnCPU( IntegrationPointsSwarm* intSwarm,  MaterialPointsSwarm* matSwarm, Particle_Index intParticleToSplit_IndexWithinCell );
void splitIntParticleByIndexOnCPU( IntegrationPointsSwarm* intSwarm,  MaterialPointsSwarm* matSwarm, Particle_Index intParticleToSplit_IndexOnCPU, Coord xi );
void _PCDVC_Calculate3D( void* pcdvc, void* _swarm, Cell_LocalIndex lCell_I );
void _PCDVC_Calculate2D( void* pcdvc, void* _swarm, Cell_LocalIndex lCell_I );  
void _PCDVC_Calculate( void* pcdvc, void* _swarm, Cell_LocalIndex lCell_I ) ;

#endif

