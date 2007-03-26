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
** $Id: DiscreteVoronoiWeights.c 376 2006-10-18 06:58:41Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>

#include <PICellerator/Voronoi/Voronoi.h>

#include "types.h"
#include "WeightsCalculator.h"
#include "DiscreteVoronoiWeights.h"

#include <assert.h>
#include <string.h>
#include <math.h>


/* Textual name of this class */
const Type DiscreteVoronoiWeights_Type = "DiscreteVoronoiWeights";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/
DiscreteVoronoiWeights* _DiscreteVoronoiWeights_New(
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
	DiscreteVoronoiWeights* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(DiscreteVoronoiWeights) );
	self = (DiscreteVoronoiWeights*)_WeightsCalculator_New( 
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

void _DiscreteVoronoiWeights_Init( void* discreteVoronoiWeights, DiscreteVoronoi* discreteVoronoi ) {
	DiscreteVoronoiWeights* self = (DiscreteVoronoiWeights*)discreteVoronoiWeights;
	
	self->isConstructed = True;
	self->discreteVoronoi = discreteVoronoi;
}

/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/
void* _DiscreteVoronoiWeights_DefaultNew( Name name ) {
	return (void*) _DiscreteVoronoiWeights_New(
			sizeof(DiscreteVoronoiWeights),
			DiscreteVoronoiWeights_Type,
			_WeightsCalculator_Delete,
			_WeightsCalculator_Print,
			_WeightsCalculator_Copy,
			_DiscreteVoronoiWeights_DefaultNew,
			_DiscreteVoronoiWeights_Construct,
			_DiscreteVoronoiWeights_Build,
			_DiscreteVoronoiWeights_Initialise,
			_WeightsCalculator_Execute,
			_WeightsCalculator_Destroy,
			_DiscreteVoronoiWeights_Calculate,
			name );
}


void _DiscreteVoronoiWeights_Construct( void* discreteVoronoiWeights, Stg_ComponentFactory* cf, void* data ) {
	DiscreteVoronoiWeights*	     self            = (DiscreteVoronoiWeights*) discreteVoronoiWeights;
	DiscreteVoronoi*             discreteVoronoi;

	_WeightsCalculator_Construct( self, cf, data );

	discreteVoronoi =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "DiscreteVoronoi", DiscreteVoronoi, True, data ) ;
	_DiscreteVoronoiWeights_Init( self, discreteVoronoi );
}

void _DiscreteVoronoiWeights_Build( void* discreteVoronoiWeights, void* data ) {
	DiscreteVoronoiWeights*	     self            = (DiscreteVoronoiWeights*) discreteVoronoiWeights;

	/* Ensure that the discrete voronoi object I'm using is built */
	Stg_Component_Build( self->discreteVoronoi, data, False );
	_WeightsCalculator_Build( self, data );
}

void _DiscreteVoronoiWeights_Initialise( void* discreteVoronoiWeights, void* data ) {
	DiscreteVoronoiWeights*	     self            = (DiscreteVoronoiWeights*) discreteVoronoiWeights;

	/* Ensure that the discrete voronoi object I'm using is initialised */
	Stg_Component_Initialise( self->discreteVoronoi, data, False );
	_WeightsCalculator_Initialise( self, data );
}

/*-------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
void _DiscreteVoronoiWeights_Calculate( void* discreteVoronoiWeights, void* _swarm, Cell_LocalIndex lCell_I ) {
	DiscreteVoronoiWeights*     self            = (DiscreteVoronoiWeights*)  discreteVoronoiWeights;
	Swarm*                      swarm           = (Swarm*) _swarm;
	FiniteElement_Mesh*         feMesh          = (FiniteElement_Mesh*)((ElementCellLayout*)swarm->cellLayout)->mesh;
	Voronoi_CellIndex           claimedVoronoiCellsCount;
	Voronoi_CellIndex           voronoiCell_I;
	ElementType*                elementType;
	IntegrationPoint*           particle;
	double                      volume;
	double                      detJac;	
	Dimension_Index             dim             = swarm->dim;
	Particle_InCellIndex        cParticle_I;
	#ifdef CAUTIOUS
	Stream*                     errorStream = Journal_Register( Error_Type, self->type );
	#endif

	/* If the cell has no particles of the current swarm - just return */
	if ( swarm->cellParticleCountTbl[lCell_I] == 0) {
		return;
	}

	WeightsCalculator_ZeroWeightsInCell( self, swarm, lCell_I );
	elementType    = FeMesh_ElementTypeAt( feMesh, lCell_I );

	/* Do Discrete Voronoi for this Cell */
	DiscreteVoronoi_CalculateForCell( self->discreteVoronoi, swarm, lCell_I );
	claimedVoronoiCellsCount = self->discreteVoronoi->claimedCellCount;

	for ( voronoiCell_I = 0 ; voronoiCell_I < claimedVoronoiCellsCount ; voronoiCell_I++ ) {
		
		cParticle_I = DiscreteVoronoi_GetParticleIndex( self->discreteVoronoi, voronoiCell_I );

		#ifdef CAUTIOUS
		Journal_Firewall( cParticle_I != CELLULAR_AUTOMATA_UNCLAIMED, errorStream,
			"Error - in %s(): while calculating weights for swarm \"%s\", in local cell %u, "
			"DiscreteVoronoi_GetParticleIndex() returned a cell particle index of UNCLAIMED for "
			"voronoi cell %u, which suggestes Voronoi thinks this element is empty, however the "
			"swarm's cellParticleCount for that element is %u.\n",
			__func__, swarm->name, lCell_I, voronoiCell_I, 
			swarm->cellParticleCountTbl[lCell_I] );
			
		Journal_Firewall( cParticle_I < swarm->cellParticleCountTbl[lCell_I], errorStream,
			"Error - in %s(): while calculating weights for swarm \"%s\", in local cell %u, "
			"DiscreteVoronoi_GetParticleIndex() returned a cell particle index of %u for "
			"voronoi cell %u, which is >= the swarm's cellParticleCount for that cell of %u.\n",
			__func__, swarm->name, lCell_I, cParticle_I, voronoiCell_I, 
			swarm->cellParticleCountTbl[lCell_I] );
		#endif
		
		particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, lCell_I, cParticle_I );

		/* Calculate Volume for Voronoi Cell */
		volume = DiscreteVoronoi_GetVolume( self->discreteVoronoi, voronoiCell_I );
		
		/* Calculate Determinant of Jacobian */
		detJac = ElementType_JacobianDeterminant( elementType, feMesh, lCell_I, particle->xi, dim );

		/* Add this volume to particle's weight - with modification for local coordinate transformation */
		particle->weight += volume/detJac;
	}

}

