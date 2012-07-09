/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**	Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
** Contributors:
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
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
** $Id: UpwindParameter.c 985 2007-11-21 00:20:24Z MirkoVelic $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <math.h>
#include <assert.h>
#include <stddef.h>

#include "mpi.h"
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include "types.h"
#include "AdvectionDiffusionSLE.h"
#include "UpwindParameter.h"
#include "Residual.h"

#define ISQRT15 0.25819888974716112567

/** AdvectionDiffusion_UpwindDiffusivity - See Brooks, Hughes 1982 Section 3.3 
 * All equations refer to this paper if not otherwise indicated */
double AdvDiffResidualForceTerm_UpwindDiffusivity(AdvDiffResidualForceTerm* self, 
                                                  AdvectionDiffusionSLE* sle, 
                                                  FeMesh* mesh, 
                                                  Element_LocalIndex lElement_I, 
                                                  Dimension_Index dim)
{
  FeVariable*                velocityField   = self->velocityField;
  Coord                      xiElementCentre = {0.0,0.0,0.0};
  double                     xiUpwind;
  double                     velocityCentre[3];
  double                     pecletNumber;
  double                     lengthScale;
  double                     upwindDiffusivity;
  Dimension_Index            dim_I;
  double*                    leastCoord;
  double*                    greatestCoord;
  Node_LocalIndex            nodeIndex_LeastValues, nodeIndex_GreatestValues;
  unsigned                   *inc;
  IArray*		 incArray;
	
  Cell_Index                 cell_I;
  IntegrationPoint*          particle;
  Particle_Index             lParticle_I;
  Particle_InCellIndex       cParticle_I;
  Particle_InCellIndex       particleCount;
	
  IntegrationPointsSwarm* swarm((IntegrationPointsSwarm*)(self->integrationSwarm));
  if(self->picSwarm)
    swarm=(IntegrationPointsSwarm*)(self->picSwarm);

  bool one_to_many=Stg_Class_IsInstance(swarm->mapper,OneToManyMapper_Type);
  if(one_to_many)
    swarm=((OneToManyMapper*)(swarm->mapper))->swarm;

  /* Compute the average diffusivity */
  /* Find Number of Particles in Element */
  cell_I=CellLayout_MapElementIdToCellId(swarm->cellLayout,lElement_I);
  particleCount=swarm->cellParticleCountTbl[cell_I];

  double averageDiffusivity(0), weight(0);
  for(cParticle_I=0; cParticle_I<particleCount; cParticle_I++)
    {
      lParticle_I=swarm->cellParticleTbl[cell_I][cParticle_I];
      particle=(IntegrationPoint*) Swarm_ParticleAt(swarm,lParticle_I);
                                                    

      IntegrationPointsSwarm* NNswarm(swarm);
      IntegrationPoint* NNparticle(particle);
      NearestNeighbor_Replace(&NNswarm,&NNparticle,lElement_I,dim);

      double temp=IntegrationPointMapper_GetDoubleFromMaterial
        (NNswarm->mapper,NNparticle,self->materialExtHandle,
         offsetof(AdvDiffResidualForceTerm_MaterialExt,diffusivity));
      averageDiffusivity+=temp*particle->weight;
      weight+=particle->weight;
    }
  averageDiffusivity/=weight;

  /* Change Diffusivity if it is too small */
  if(averageDiffusivity<SUPG_MIN_DIFFUSIVITY)
    averageDiffusivity = SUPG_MIN_DIFFUSIVITY;
	
  /* Calculate Velocity At Middle of Element - See Eq. 3.3.6 */
  FeVariable_InterpolateFromMeshLocalCoord(velocityField,mesh,lElement_I,
                                           xiElementCentre,velocityCentre);
	
  /* Calculate Length Scales - See Fig 3.4 - ASSUMES BOX MESH TODO - fix */
  incArray = self->incarray;
  FeMesh_GetElementNodes( mesh, lElement_I, incArray );
  inc = (unsigned*)IArray_GetPtr( incArray );

  nodeIndex_LeastValues = inc[0];
  nodeIndex_GreatestValues = inc[IArray_GetSize(incArray)-1];
  leastCoord    = Mesh_GetVertex( mesh, nodeIndex_LeastValues );
  greatestCoord = Mesh_GetVertex( mesh, nodeIndex_GreatestValues );
	
  upwindDiffusivity = 0.0;
  for(dim_I=0; dim_I<dim; dim_I++)
    {
      lengthScale=greatestCoord[dim_I]-leastCoord[dim_I];
		
      /* Calculate Peclet Number (alpha) - See Eq. 3.3.5 */
      pecletNumber=velocityCentre[dim_I]*lengthScale/(2.0*averageDiffusivity);
		
      /* Calculate Upwind Local Coordinate - See Eq. 3.3.4 and (2.4.2,
         3.3.1 and 3.3.2) */
      xiUpwind=AdvDiffResidualForceTerm_UpwindParam(self,pecletNumber);
		
      /* Calculate Upwind Thermal Diffusivity - See Eq. 3.3.3  */
      upwindDiffusivity+=xiUpwind*velocityCentre[dim_I]*lengthScale;
    }
  upwindDiffusivity*=ISQRT15;         /* See Eq. 3.3.11 */

  upwindDiffusivity/=16; /* This factor is because we are using Q2
                           instead of Q1 elements, so the spacing
                           needs to be adjusted.  That gives a factor
                           of 2, but we make it a factor of 16 to
                           ensure stability even when the mesh is
                           distorted. */
	
  return upwindDiffusivity;
}


/** AdvectionDiffusion_UpwindXiExact - Brooks, Hughes 1982 equation 2.4.2
 *\f$ \bar \xi = coth( \alpha ) - \frac{1}{\alpha} \f$ */
double AdvDiffResidualForceTerm_UpwindXiExact( void* residual, double pecletNumber ) {
	if (fabs(pecletNumber) < 1.0e-8 )
		return 0.33333333333333 * pecletNumber;
	else if (pecletNumber < -20.0)
		return -1.0 - 1.0/pecletNumber;
	else if (pecletNumber > 20.0)
		return +1.0 - 1.0/pecletNumber;
		
	return cosh( pecletNumber )/sinh( pecletNumber ) - 1.0/pecletNumber;
}

/** AdvectionDiffusion_UpwindXiDoublyAsymptoticAssumption - Brooks, Hughes 1982 equation 3.3.1
 * Simplification of \f$ \bar \xi = coth( \alpha ) - \frac{1}{\alpha} \f$ from Brooks, Hughes 1982 equation 2.4.2
 * \f[
\bar \xi \sim \left\{ \begin{array}{rl} 
             -1                 &for \quad \alpha <= -3 \\
             \frac{\alpha}{3}   &for \quad -3 < \alpha <= 3 \\
             +1                 &for \quad \alpha > +3
             \end{array} \right.  
           
\f]*/
double AdvDiffResidualForceTerm_UpwindXiDoublyAsymptoticAssumption( void* residual, double pecletNumber ) {
	if (pecletNumber <= -3.0)
		return -1;
	else if (pecletNumber <= 3.0)
		return 0.33333333333333 * pecletNumber;
	else
		return 1.0;
}
	
/** AdvectionDiffusion_UpwindXiCriticalAssumption - Brooks, Hughes 1982 equation 3.3.2
 * Simplification of \f$ \bar \xi = coth( \alpha ) - \frac{1}{\alpha} \f$ from Brooks, Hughes 1982 equation 2.4.2
 * \f[
  \bar \xi \sim \left\{ \begin{array}{rl}
              -1 - \frac{1}{\alpha}   &for \quad \alpha <= -1 \\
               0                      &for \quad -1 < \alpha <= +1 \\
              +1 - \frac{1}{\alpha}   &for \quad \alpha > +1          
              \end{array} \right.              
\f]    */

double AdvDiffResidualForceTerm_UpwindXiCriticalAssumption( void* residual, double pecletNumber ) {
	if (pecletNumber <= -1.0)
		return -1.0 - 1.0/pecletNumber;
	else if (pecletNumber <= 1.0)
		return 0.0;
	else
		return 1.0 - 1.0/pecletNumber;
}
	


