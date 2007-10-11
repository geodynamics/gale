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
** $Id: testDiscreteVoronoi.c 518 2007-10-11 08:07:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/Voronoi/Voronoi.h>

#include <math.h>
#include <string.h>
#include <assert.h>

const Type TestDiscreteVoronoi_Type = "TestDiscreteVoronoi";

void PICellerator_testDiscreteVoronoi( DomainContext* context ) {
	Swarm*              swarm;
	DiscreteVoronoi*    discreteVoronoi;
	Bool                particlesCorrectlyAssociated = True;
	double              totalAreaError               = 0.0;
	Cell_Index          lCell_I;
	Stream*             stream                       = Journal_Register( Info_Type, CURR_MODULE_NAME );
	Particle_InCellIndex assignedParticle, closestParticle;

	/* Get Stg_Components */
	swarm = (Swarm*)                     LiveComponentRegister_Get( context->CF->LCRegister, "swarm" );
	assert(swarm);
	discreteVoronoi = (DiscreteVoronoi*) LiveComponentRegister_Get( context->CF->LCRegister, "discreteVoronoi" );
	assert(discreteVoronoi);

	for ( lCell_I = 0 ; lCell_I < swarm->cellLocalCount ; lCell_I++ ) {
		double                      area = 0.0;
		Cell_Index                  vCell_I;
		Voronoi_CellIndex           claimedCellsCount;
		double                      cellArea;
		Coord                       centroid;

		/* Do Voronoi */
		DiscreteVoronoi_CalculateForCell( discreteVoronoi, swarm, lCell_I );
		claimedCellsCount = discreteVoronoi->claimedCellCount;
		
		/* Test Area */
		for ( vCell_I = 0 ; vCell_I < claimedCellsCount ; vCell_I++ ) {
			area += DiscreteVoronoi_GetVolume( discreteVoronoi, vCell_I );
			DiscreteVoronoi_GetCentroid( discreteVoronoi, vCell_I, centroid );

			assignedParticle = DiscreteVoronoi_GetParticleIndex( discreteVoronoi, vCell_I );
			closestParticle = Swarm_FindClosestParticleInCell( swarm, lCell_I, context->dim, centroid, NULL );
			if ( assignedParticle != closestParticle ) {
				particlesCorrectlyAssociated = False;
				Journal_Printf( stream, "Dodgy voronoi cell in element %u and cell %u -- assigned particle = %u closest particle = %u - centroid = %g %g %g\n", 
						lCell_I, vCell_I, assignedParticle, closestParticle, centroid[0], centroid[1], centroid[2] );
			}
					
		}
		if (context->dim == 2)
			cellArea = StGermain_ConvexQuadrilateralArea( 
					(double*)swarm->cellPointTbl[lCell_I][0],
					(double*)swarm->cellPointTbl[lCell_I][1],
					(double*)swarm->cellPointTbl[lCell_I][2],
					(double*)swarm->cellPointTbl[lCell_I][3],
					context->dim);
		else 
			cellArea = StGermain_ParallelepipedVolume( 
					(double*)swarm->cellPointTbl[lCell_I][0],
					(double*)swarm->cellPointTbl[lCell_I][1],
					(double*)swarm->cellPointTbl[lCell_I][3],
					(double*)swarm->cellPointTbl[lCell_I][4] );

		totalAreaError += fabs( area - cellArea );
	}

	/* Write to file */
	Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, "output.dat" );
	if ( totalAreaError < 1.0e-8 )
		Journal_Printf( stream, "totalAreaError = 0.0\n" );
	else 
		Journal_PrintValue( stream, totalAreaError );
	Journal_PrintBool( stream, particlesCorrectlyAssociated );
}
	

void _testDiscreteVoronoi_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	DomainContext* context;
	context = (DomainContext*)Stg_ComponentFactory_ConstructByName( cf, "context", DomainContext, True, data );
	ContextEP_ReplaceAll( context, AbstractContext_EP_Execute, PICellerator_testDiscreteVoronoi );
}


void* _testDiscreteVoronoi_DefaultNew( Name name ) {
	return _Codelet_New(
			sizeof( Codelet ),
			TestDiscreteVoronoi_Type,
			_Codelet_Delete,
			_Codelet_Print,
			_Codelet_Copy,
			_testDiscreteVoronoi_DefaultNew,
			_testDiscreteVoronoi_Construct,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}


Index testDiscreteVoronoi_Register( PluginsManager* pluginsManager ) {
	Index result;

	result = PluginsManager_Submit( pluginsManager, TestDiscreteVoronoi_Type, "0",
		_testDiscreteVoronoi_DefaultNew );

	return result;
}
