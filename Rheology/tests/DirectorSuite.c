/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**   Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**   Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**   Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**   Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**   Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**   Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
** Role:
**   Tests the DirectorSuite
**
** $Id: testDirector.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include "DirectorSuite.h"

#define TOLERANCE 0.01
#define RANDOM_DIRECTOR_ERROR 1.0
#define CURR_MODULE_NAME "DirectorSuite"

typedef struct {
} DirectorSuiteData;

double dt( UnderworldContext* context ) { return 0.05; }

void test( UnderworldContext* context ) {
	AlignmentSwarmVariable* alignment              = (AlignmentSwarmVariable*)LiveComponentRegister_Get( context->CF->LCRegister, "alignment" );
	FeVariable*             velocityField          = Stg_ComponentFactory_ConstructByName( context->CF, "VelocityField", FeVariable, True, 0 );
	Director*               director;
	Particle_Index          lParticle_I;
	GlobalParticle*         particle;
	double                  time                   = context->currentTime + context->dt;
	Swarm*                  swarm;
	double                  error                  = 0.0;
	double                  alignmentValue;
	XYZ                     normal;
	double                  angle                  = 0.5 * M_PI - atan(1.0/(2.0*time) );
	XYZ                     velocity;
	double                  analyticAlignmentValue = 1.0 - cos( angle );

	swarm = alignment->swarm;
	director = alignment->director;

	for ( lParticle_I = 0 ; lParticle_I < swarm->particleLocalCount ; lParticle_I++ ) {
		particle = (GlobalParticle*)Swarm_ParticleAt( swarm, lParticle_I );
		SwarmVariable_ValueAt( alignment, lParticle_I, &alignmentValue );
		SwarmVariable_ValueAt( director->directorSwarmVariable, lParticle_I, normal );

		FieldVariable_InterpolateValueAt( velocityField, particle->coord, velocity );

		error += fabs( alignmentValue - analyticAlignmentValue );
	}
	error /= (double) swarm->particleLocalCount;

	pcu_check_true( error < TOLERANCE );
}

void testRandom( UnderworldContext* context ) {
	AlignmentSwarmVariable* alignment              = (AlignmentSwarmVariable*) LiveComponentRegister_Get( context->CF->LCRegister, "alignment" );
	Director*               director;
	Particle_Index          lParticle_I;
	GlobalParticle*         particle;
	Swarm*                  swarm;
	XYZ                     normal;
	int                     ii;
	double                  dotProduct;
	double                  angleBetween;
	double                  circleAngle;
	int                     circleAngleCounts[36] 	= {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	Bool                    circleErrorFlag = False;
	int                     circleAngleSum;
	double                  circleAngleAverage;
	double                  circleAngleStdDev;
	int                     circleAngleLowerBound;
	int                     circleAngleUpperBound;
	double			gCircleAngleAverage;
	int                     gCircleAngleCounts[36] 	= {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	int			gCircleAngleSum;

int rank, rank_i;
MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
	swarm = alignment->swarm;
	director = alignment->director;

	for ( lParticle_I = 0 ; lParticle_I < swarm->particleLocalCount ; lParticle_I++ ) {
		particle = (GlobalParticle*)Swarm_ParticleAt( swarm, lParticle_I );
		SwarmVariable_ValueAt( director->directorSwarmVariable, lParticle_I, normal );
		/* Calculate dot product between normal and (0,1), then get an angle */

printf( "%d: particle[%d] normal: [%e,%e]\n", rank, lParticle_I, normal[0], normal[1] );

		dotProduct = normal[0] * 0 + normal[1] * 1;
		angleBetween = acos( dotProduct ) * 180 / M_PI; 
		if ( normal[0] > 0 ) {
			circleAngle = angleBetween;
		}
		else {
			circleAngle = 360 - angleBetween;
		}
		circleAngleCounts[ (int)(circleAngle+0.5) / 10 ] += 1;
	}
for( ii=0; ii<36; ii++ )
  MPI_Allreduce( &circleAngleCounts[ii], &gCircleAngleCounts[ii], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
	//MPI_Allreduce( circleAngleCounts, gCircleAngleCounts, 36, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

	circleAngleAverage = (double)swarm->particleLocalCount / 36;
	MPI_Allreduce( &circleAngleAverage, &gCircleAngleAverage, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

for( rank_i=0; rank_i<2; rank_i++ ) {
  if( rank_i == rank ) {
    printf( "%d: angle avg. - local:\t%e\tglobal:%e\n", rank, circleAngleAverage, gCircleAngleAverage );
    for( ii=0; ii<36; ii++ )
      printf( "%d: angle counts[%2d] - local:\t%d\tglobal:\t%d\n", rank, ii, circleAngleCounts[ii], gCircleAngleCounts[ii] );
  }
  MPI_Barrier( MPI_COMM_WORLD );
}

	/*NB. This definition is determined based on a set no. of particlesPerCell. Currently this value is = 20 */
	#define TheoreticalStandardDeviation 13.64
	circleAngleLowerBound = (int)(floor(/*circleAngleAverage*/gCircleAngleAverage - 3*TheoreticalStandardDeviation));
	if (circleAngleLowerBound < 0 ) {
		circleAngleLowerBound = 0;
	}
	circleAngleUpperBound = (int)(floor(/*circleAngleAverage*/gCircleAngleAverage + 3*TheoreticalStandardDeviation));
	
        circleAngleSum = 0;
	for ( ii =0; ii < 36; ii++ ) {	
		if (( gCircleAngleCounts[ii] < circleAngleLowerBound )|| (gCircleAngleCounts[ii] > circleAngleUpperBound)) {
			circleErrorFlag = True;			
		}
		//circleAngleSum = ((circleAngleAverage - circleAngleCounts[ii]) * 
		//				(circleAngleAverage - circleAngleCounts[ii])) + circleAngleSum;
		circleAngleSum = ((gCircleAngleAverage - gCircleAngleCounts[ii]) * 
						(gCircleAngleAverage - gCircleAngleCounts[ii])) + circleAngleSum;
	}
	//MPI_Allreduce( &circleAngleSum, &gCircleAngleSum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
	/* Calculate standard deviation */
	//circleAngleStdDev = sqrt(gCircleAngleSum / (36-1));
	circleAngleStdDev = sqrt(circleAngleSum / (36-1));

printf( "%d: std. dev: %e\n", rank, circleAngleStdDev );
	
	pcu_check_true( fabs( circleAngleStdDev - TheoreticalStandardDeviation ) < RANDOM_DIRECTOR_ERROR );
	pcu_check_true( !circleErrorFlag );	
}

#define DIR_TEST_NUM_MAT 3
#define DIR_X_AXIS 0
#define DIR_Y_AXIS 1
#define DIR_Z_AXIS 2
#define ANGLE_ERROR 1e-15

void testPerMaterial( UnderworldContext* context ) {
	AlignmentSwarmVariable* alignment              = (AlignmentSwarmVariable*) LiveComponentRegister_Get( context->CF->LCRegister, "alignment" );
	Materials_Register*     materials_Register     = context->materials_Register;
	Director*               director;
	Particle_Index          lParticle_I;
	Swarm*                  swarm;
	XYZ                     testVector;
	Material_Index          materialsCount;
	int                     material_I;
	XYZ*                    matDirectionVectors;
	double                  angle;
	Bool			angleFailure		= False;
	Material*		material;
	
	swarm = alignment->swarm;
	director = alignment->director;
	materialsCount = Materials_Register_GetCount( materials_Register);
	
	/*  construct test for testDirectorPerMaterial.xml  */
	/* assume a direction for each material and check that */
	/* each particle has the same direction for it's material type. */
	/* and that it's equal to the default value for each material. */
	/* get materials and their initial directions */
	
	/* check that there are at least 3 materials */
	/* check that the initialDirections are == certain defaults */
	matDirectionVectors = Memory_Alloc_Array(XYZ, DIR_TEST_NUM_MAT, "materialDirectionVectors");

	/* Define vectors */
	matDirectionVectors[0][DIR_X_AXIS] = 1.0; matDirectionVectors[0][DIR_Y_AXIS] = 0.0; matDirectionVectors[0][DIR_Z_AXIS] = 0.0;
	matDirectionVectors[1][DIR_X_AXIS] = 0.0; matDirectionVectors[1][DIR_Y_AXIS] = 1.0; matDirectionVectors[1][DIR_Z_AXIS] = 0.0;
	matDirectionVectors[2][DIR_X_AXIS] = 1.0; matDirectionVectors[2][DIR_Y_AXIS] = 1.0; matDirectionVectors[2][DIR_Z_AXIS] = 0.0;
	
	for ( lParticle_I = 0 ; lParticle_I < swarm->particleLocalCount ; lParticle_I++ ) {
		/* get particle's material type */
		material_I = MaterialPointsSwarm_GetMaterialIndexAt( swarm, lParticle_I );

material = Materials_Register_GetByIndex( materials_Register, material_I );
		/* get particle's initialDirection */
		SwarmVariable_ValueAt( director->directorSwarmVariable, lParticle_I, testVector );
		/* check that angle is within error bar of correct value */
		angle = StGermain_AngleBetweenVectors(matDirectionVectors[material_I],testVector, swarm->dim);
	
		if( fabs(angle) > ANGLE_ERROR )
			angleFailure = True;
	}

	pcu_check_true( !angleFailure );
}

void DirectorSuite_Setup( DirectorSuiteData* data ) {
}

void DirectorSuite_Teardown( DirectorSuiteData* data ) {
}

void DirectorSuite_Test( DirectorSuiteData* data ) {
	UnderworldContext* 	context;
	Stg_ComponentFactory*	cf;
	char			xml_input[PCU_PATH_MAX];

	pcu_filename_input( "testDirector.xml", xml_input );
	context = _UnderworldContext_DefaultNew( "context" );
	cf = stgMainInitFromXML( xml_input, MPI_COMM_WORLD, context );
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, test );
	EP_AppendClassHook( Context_GetEntryPoint( context, FiniteElementContext_EP_CalcDt ), dt, context );
	stgMainBuildAndInitialise( cf );
	stgMainLoop( cf );
	stgMainDestroy( cf );
}

void DirectorSuite_TestRandom( DirectorSuiteData* data ) {
	UnderworldContext* 	context;
	Stg_ComponentFactory*	cf;
	char			xml_input[PCU_PATH_MAX];

	pcu_filename_input( "testDirectorRandom.xml", xml_input );
	context = _UnderworldContext_DefaultNew( "context" );
	cf = stgMainInitFromXML( xml_input, MPI_COMM_WORLD, context );
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, testRandom );
	stgMainBuildAndInitialise( cf );
	stgMainLoop( cf );
	stgMainDestroy( cf );
}

void DirectorSuite_TestPerMaterial( DirectorSuiteData* data ) {
	UnderworldContext* 	context;
	Stg_ComponentFactory*	cf;
	char			xml_input[PCU_PATH_MAX];

	pcu_filename_input( "testDirectorPerMaterial.xml", xml_input );
	context = _UnderworldContext_DefaultNew( "context" );
	cf = stgMainInitFromXML( xml_input, MPI_COMM_WORLD, context );
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, testPerMaterial );
	stgMainBuildAndInitialise( cf );
	stgMainLoop( cf );
	stgMainDestroy( cf );
}

void DirectorSuite_TestPerMaterial2( DirectorSuiteData* data ) {
	UnderworldContext* 	context;
	Stg_ComponentFactory*	cf;
	char			xml_input[PCU_PATH_MAX];

	pcu_filename_input( "testDirectorPerMaterial2.xml", xml_input );
	context = _UnderworldContext_DefaultNew( "context" );
	cf = stgMainInitFromXML( xml_input, MPI_COMM_WORLD, context );
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, testPerMaterial );
	stgMainBuildAndInitialise( cf );
	stgMainLoop( cf );
	stgMainDestroy( cf );
}

void DirectorSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, DirectorSuiteData );
   pcu_suite_setFixtures( suite, DirectorSuite_Setup, DirectorSuite_Teardown );
   pcu_suite_addTest( suite, DirectorSuite_Test );
   //pcu_suite_addTest( suite, DirectorSuite_TestRandom );
   pcu_suite_addTest( suite, DirectorSuite_TestPerMaterial );
   pcu_suite_addTest( suite, DirectorSuite_TestPerMaterial2 );
}
