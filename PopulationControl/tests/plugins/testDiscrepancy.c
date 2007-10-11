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
** $Id: testDiscrepancy.c 518 2007-10-11 08:07:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/Voronoi/Voronoi.h>
#include <PICellerator/PopulationControl/PopulationControl.h>

/* Discrepancy is a measure of the spacing of the particles 
 * For information on 'Discrepancy' check out:
 * P. Bratley, B.L. Fox,  Implementing Sobol's Quasirandom Sequence Generator,  
 * ACM Transactions on Mathematical Software, March 1988, Volume 14, pp. 88 - 100
 * http://portal.acm.org/ft_gateway.cfm?id=214372&type=pdf&coll=GUIDE&dl=GUIDE&CFID=45681167&CFTOKEN=8969368
 *
 * also http://en.wikipedia.org/wiki/Quasi-random */

#define VOLUME_OF_CUBE 1.0

double Swarm_ParticleDiscrepancy( void* swarm, Cell_LocalIndex lCell_I, Particle_InCellIndex cParticle_I ) {
	Swarm*               self    = (Swarm*) swarm;
	Particle_InCellIndex countSn = 0;
	Dimension_Index      dim     = self->dim;
	double               volumeGx;
	Particle_InCellIndex smallerCellParticle_I;
	Particle_InCellIndex cellParticleCount = self->cellParticleCountTbl[ lCell_I ];
	GlobalParticle*      particle          = (GlobalParticle*)Swarm_ParticleInCellAt( self, lCell_I, cParticle_I );
	double*              coord             = particle->coord;
	GlobalParticle*      smallerCellParticle;
	double*              startCellCoord    = (double*) self->cellPointTbl[lCell_I][0];
	double*              endCellCoord      = (double*) self->cellPointTbl[lCell_I][ dim == 2 ? 2 : 6 ];

	/* Calculate volume of Gx */
	volumeGx = 
		(coord[ I_AXIS ] - startCellCoord[ I_AXIS ])/( endCellCoord[ I_AXIS ] - startCellCoord[ I_AXIS ]) *
		(coord[ J_AXIS ] - startCellCoord[ J_AXIS ])/( endCellCoord[ J_AXIS ] - startCellCoord[ J_AXIS ]) ;
	if ( dim == 3 )
		volumeGx *= (coord[ K_AXIS ] - startCellCoord[ K_AXIS ])/( endCellCoord[ K_AXIS ] - startCellCoord[ K_AXIS ]) ;
	
	for ( smallerCellParticle_I = 0 ; smallerCellParticle_I < cellParticleCount ; smallerCellParticle_I++ ) {
		smallerCellParticle = (GlobalParticle*)Swarm_ParticleInCellAt( self, lCell_I, smallerCellParticle_I );
	
		/* Only include particles with lower coords than this particle in counting Sn */
		if ( smallerCellParticle_I == cParticle_I )
			continue;

		if ( smallerCellParticle->coord[ I_AXIS ] >= coord[ I_AXIS ] )
			continue;
		
		if ( smallerCellParticle->coord[ J_AXIS ] >= coord[ J_AXIS ] )
			continue;
		
		if ( dim == 3 && smallerCellParticle->coord[ K_AXIS ] >= coord[ K_AXIS ] )
			continue;

		
		/* Otherwise, particle is within Gx - Therfore contributes to Sn */
		countSn++;
	}

	return fabs( ((double) countSn / (double) cellParticleCount - volumeGx / VOLUME_OF_CUBE ) );
}

double Swarm_CellDiscrepancy( void* swarm, Cell_LocalIndex lCell_I ) {
	Swarm*               self    = (Swarm*) swarm;
	double               discrepancy = 0.0;
	double               result;
	Particle_InCellIndex cParticle_I;
	Particle_InCellIndex cellParticleCount = self->cellParticleCountTbl[ lCell_I ];

	for ( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {
		result = Swarm_ParticleDiscrepancy( self, lCell_I, cParticle_I );

		/* Discrepancy is the supremum of all the particles */
		if ( discrepancy < result )
			discrepancy = result;
	}

	return discrepancy;
}

double Swarm_AverageDiscrepancy( void* swarm ) {
	Swarm*               self    = (Swarm*) swarm;
	double               localDiscrepancy  = 0.0;
	double               globalDiscrepancy = 0.0;
	Cell_LocalIndex      lCell_I;
	Cell_Index           globalCellCount;

	for ( lCell_I = 0 ; lCell_I < self->cellLocalCount ; lCell_I++ ) {
		localDiscrepancy += Swarm_CellDiscrepancy( self, lCell_I );
	}

	MPI_Allreduce( &localDiscrepancy, &globalDiscrepancy, 1, MPI_DOUBLE, MPI_MAX, self->comm );
	MPI_Allreduce( &self->cellLocalCount, &globalCellCount, 1, MPI_UNSIGNED, MPI_MAX, self->comm );

	return globalDiscrepancy/(double)globalCellCount;
}


void testDiscrepancy_Function( DomainContext* context ) {
	Stream*           stream = Journal_Register( Info_Type, CURR_MODULE_NAME );
	double            discrepancy;
	double            maxDiscrepancy;
	Particle_Index    globalParticleCount;
	Swarm*            swarm            = Stg_ComponentFactory_ConstructByName( context->CF, "swarm", Swarm, True, 0 );
	SplittingRoutine* splittingRoutine = Stg_ComponentFactory_ConstructByName( 
						context->CF, 
						"splittingRoutine", 
						SplittingRoutine, 
						False, 
						0 );
	RemovalRoutine*   removalRoutine   = Stg_ComponentFactory_ConstructByName( 
						context->CF, 
						"removalRoutine", 
						RemovalRoutine, 
						False, 
						0 );

	/* Output before population control */
	discrepancy = Swarm_AverageDiscrepancy( swarm );
	MPI_Allreduce( &swarm->particleLocalCount, &globalParticleCount, 1, MPI_UNSIGNED, MPI_MAX, swarm->comm );
	Journal_Printf( stream, "Before population control - global particle count = %u - discrepancy = %.4g\n",
			globalParticleCount, discrepancy );

	/* Do population control */
	if ( removalRoutine ) 
		RemovalRoutine_RemoveFromSwarm( removalRoutine, swarm );
	if ( splittingRoutine ) 
		Stg_Component_Execute( splittingRoutine, swarm, True );

	/* Output after population control */
	discrepancy = Swarm_AverageDiscrepancy( swarm );
	MPI_Allreduce( &swarm->particleLocalCount, &globalParticleCount, 1, MPI_UNSIGNED, MPI_MAX, swarm->comm );
	Journal_Printf( stream, "Before population control - global particle count = %u - discrepancy = %.4g\n",
			globalParticleCount, discrepancy );
	
	/* Output information to file */
	Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, "discrepancy.dat" );
	if ( removalRoutine )
		Journal_Printf( stream, "Testing %s\n", removalRoutine->type );
	if ( splittingRoutine )
		Journal_Printf( stream, "Testing %s\n", splittingRoutine->type );
	maxDiscrepancy = Dictionary_GetDouble( context->dictionary, "maxDiscrepancy" );
	if ( discrepancy <= maxDiscrepancy ) 
		Journal_Printf( stream, "Test passed with max discrepancy = %.4g\n", maxDiscrepancy );
	else 
		Journal_Printf( stream, "Test failed with discrepancy = %g and max discrepancy = %.4g\n", 
				discrepancy, maxDiscrepancy );

	exit( EXIT_SUCCESS );
}
	
const Type testDiscrepancy_Type = "testDiscrepancy";

typedef struct {
	__Codelet
} testDiscrepancy;

void _testDiscrepancy_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	DomainContext* context;
	context = (DomainContext*)Stg_ComponentFactory_ConstructByName( cf, "context", DomainContext, True, data ); 
	ContextEP_Append( context, AbstractContext_EP_Initialise, testDiscrepancy_Function );
}

void* _testDiscrepancy_Default_New( Name name ) {
	return _Codelet_New(
		sizeof( Codelet ),
		testDiscrepancy_Type,
		_Codelet_Delete,
		_Codelet_Print,
		_Codelet_Copy,
		_testDiscrepancy_Default_New, 
		_testDiscrepancy_Construct,
		_Codelet_Build,
		_Codelet_Initialise,
		_Codelet_Execute,
		_Codelet_Destroy,
		name );
}

Index testDiscrepancy_Register( PluginsManager* pluginsManager ) {
        Index result;

        result = PluginsManager_Submit( pluginsManager, testDiscrepancy_Type, "0", _testDiscrepancy_Default_New );

        return result;
}
	
