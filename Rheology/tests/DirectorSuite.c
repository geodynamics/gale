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

#define TOLERANCE 1e-10
#define RANDOM_DIRECTOR_ERROR 1.0
#define CURR_MODULE_NAME "DirectorSuite"

typedef struct {
} DirectorSuiteData;

double dt( UnderworldContext* context ) { return 0.05; }

void _Director_Intermediate_Replace( void* director, Index lParticle_I ) {}

void test( UnderworldContext* context ) {
	Director*			director = (Director*)LiveComponentRegister_Get( context->CF->LCRegister, "director" );
	Particle_Index		lParticle_I;
	GlobalParticle*	particle;
	double				time = context->currentTime + context->dt;
	Swarm*				swarm	= (Swarm*)LiveComponentRegister_Get( context->CF->LCRegister, "materialSwarm" );
	XYZ					normal;
	double				error = 0.0;
   double				angle = 0.5 * M_PI - atan(1.0/(2.0*time) );
   double				angleDirector;
	double				gError;
	int					particleGlobalCount;
   int ierr;

	for ( lParticle_I = 0 ; lParticle_I < swarm->particleLocalCount ; lParticle_I++ ) {
		particle = (GlobalParticle*)Swarm_ParticleAt( swarm, lParticle_I );
		SwarmVariable_ValueAt( director->directorSwarmVariable, lParticle_I, normal );
      angleDirector = atan(-normal[1]/normal[0]);
		error += fabs( angleDirector - angle );
	}
	ierr = MPI_Allreduce( &error, &gError, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	ierr = MPI_Allreduce( &swarm->particleLocalCount, &particleGlobalCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

	//error /= (double) swarm->particleLocalCount;
	//pcu_check_true( error < TOLERANCE );
	gError /= particleGlobalCount;
	pcu_check_true( gError < TOLERANCE );
}

void testRandom( UnderworldContext* context ) {
	AlignmentSwarmVariable* alignment = (AlignmentSwarmVariable*) LiveComponentRegister_Get( context->CF->LCRegister, "alignment" );
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
	double						gCircleAngleAverage;
	int                     gCircleAngleCounts[36] 	= {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   int ierr;

	swarm = alignment->swarm;
	director = alignment->director;

	for ( lParticle_I = 0 ; lParticle_I < swarm->particleLocalCount ; lParticle_I++ ) {
		particle = (GlobalParticle*)Swarm_ParticleAt( swarm, lParticle_I );
		SwarmVariable_ValueAt( director->directorSwarmVariable, lParticle_I, normal );
		/* Calculate dot product between normal and (0,1), then get an angle */

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
	ierr = MPI_Allreduce( circleAngleCounts, gCircleAngleCounts, 36, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

	circleAngleAverage = (double)swarm->particleLocalCount / 36;
	ierr = MPI_Allreduce( &circleAngleAverage, &gCircleAngleAverage, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

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
	/* Calculate standard deviation */
	circleAngleStdDev = sqrt(circleAngleSum / (36-1));

	pcu_check_true( fabs( circleAngleStdDev - TheoreticalStandardDeviation ) < RANDOM_DIRECTOR_ERROR );
	pcu_check_true( !circleErrorFlag );	
}

#define DIR_TEST_NUM_MAT 3
#define DIR_X_AXIS 0
#define DIR_Y_AXIS 1
#define DIR_Z_AXIS 2
#define ANGLE_ERROR 1e-7

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

		/* get particle's initialDirection */
		SwarmVariable_ValueAt( director->directorSwarmVariable, lParticle_I, testVector );
		/* check that angle is within error bar of correct value */
		angle = StGermain_AngleBetweenVectors(matDirectionVectors[material_I],testVector, swarm->dim);
	
		if( fabs(angle) > ANGLE_ERROR )
			angleFailure = True;
	}

	pcu_check_true( !angleFailure );
}

void testPerMaterial2( UnderworldContext* context ) {
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
	
	swarm = alignment->swarm;
	director = alignment->director;
	materialsCount = Materials_Register_GetCount( materials_Register);
	
	matDirectionVectors = Memory_Alloc_Array(XYZ, DIR_TEST_NUM_MAT, "materialDirectionVectors");

	/* Define vectors */
	matDirectionVectors[0][DIR_X_AXIS] = 1.0; matDirectionVectors[0][DIR_Y_AXIS] = 0.0; matDirectionVectors[0][DIR_Z_AXIS] = 0.0;
	matDirectionVectors[1][DIR_X_AXIS] = 0.0; matDirectionVectors[1][DIR_Y_AXIS] = 1.0; matDirectionVectors[1][DIR_Z_AXIS] = 0.0;
	matDirectionVectors[2][DIR_X_AXIS] = 1.0; matDirectionVectors[2][DIR_Y_AXIS] = 1.0; matDirectionVectors[2][DIR_Z_AXIS] = 0.0;
	/* normalise the vectors */
	matDirectionVectors[2][DIR_X_AXIS] /= pow( 2.0, 0.5 );
	matDirectionVectors[2][DIR_Y_AXIS] /= pow( 2.0, 0.5 );
	
	for ( lParticle_I = 0 ; lParticle_I < swarm->particleLocalCount ; lParticle_I++ ) {
		/* get particle's material type */
		material_I = MaterialPointsSwarm_GetMaterialIndexAt( swarm, lParticle_I );

		/* get particle's initialDirection */
		SwarmVariable_ValueAt( director->directorSwarmVariable, lParticle_I, testVector );
		/* check that angle is within error bar of correct value */
		angle = StGermain_AngleBetweenVectors(matDirectionVectors[material_I],testVector, swarm->dim);
	
		/* the '1' material has a random normal vector, so the angle between this and the test vector should be non-zero */
		if( material_I == 1 && fabs(angle) < ANGLE_ERROR )
			angleFailure = True;

		if( material_I != 1 && fabs(angle) > ANGLE_ERROR )
			angleFailure = True;
	}

	pcu_check_true( !angleFailure );
}

void DirectorSuite_Setup( DirectorSuiteData* data ) {
	Journal_Enable_AllTypedStream( False );
}

void DirectorSuite_Teardown( DirectorSuiteData* data ) {
	Journal_Enable_AllTypedStream( True );
}

void DirectorSuite_Test( DirectorSuiteData* data ) {
	UnderworldContext*		context;
	Stg_ComponentFactory*	cf;
	char							xml_input[PCU_PATH_MAX];
   Director*					director;
   /* this test checks that the director will eventually spin to point in the direction of maximum shear.
      initially the director is set to point 90 degrees to this direction, so that the director routines should result in it 
      spinning 90 degrees eventually.   the test checks how the angle evolves in time against the analytic result, where we have disabled 
      the intermediate re-normalisation of the director vectors */
	pcu_filename_input( "testDirector.xml", xml_input );
	cf = stgMainInitFromXML( xml_input, MPI_COMM_WORLD, NULL );
	context = (UnderworldContext*)LiveComponentRegister_Get( cf->LCRegister, "context" );
	Stream_Enable( context->info, False );
	Stream_Enable( context->verbose, False );
	Stream_Enable( context->debug, False );

	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, test );
	EP_AppendClassHook( Context_GetEntryPoint( context, FiniteElementContext_EP_CalcDt ), dt, context );
	stgMainBuildAndInitialise( cf );
	director = (Director*) LiveComponentRegister_Get( context->CF->LCRegister, "director" );
	/* we disable the intermediate director step, which normalises the director.  this is done so that we can have an analytic solution, but a more inclusive test 
	  would be a good idea */
	director->_intermediate = _Director_Intermediate_Replace;
	stgMainLoop( cf );
	stgMainDestroy( cf );
}

void DirectorSuite_TestRandom( DirectorSuiteData* data ) {
	UnderworldContext*		context;
	Stg_ComponentFactory*	cf;
	char							xml_input[PCU_PATH_MAX];

	pcu_filename_input( "testDirectorRandom.xml", xml_input );
	cf = stgMainInitFromXML( xml_input, MPI_COMM_WORLD, NULL );
	context = (UnderworldContext*)LiveComponentRegister_Get( cf->LCRegister, "context" );
	Stream_Enable( context->info, False );
	Stream_Enable( context->verbose, False );
	Stream_Enable( context->debug, False );
	
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, testRandom );
	stgMainBuildAndInitialise( cf );
	stgMainLoop( cf );
	stgMainDestroy( cf );
}

void DirectorSuite_TestPerMaterial( DirectorSuiteData* data ) {
	UnderworldContext*		context;
	Stg_ComponentFactory*	cf;
	char							xml_input[PCU_PATH_MAX];

	pcu_filename_input( "testDirectorPerMaterial.xml", xml_input );
	cf = stgMainInitFromXML( xml_input, MPI_COMM_WORLD, NULL );
	context = (UnderworldContext*)LiveComponentRegister_Get( cf->LCRegister, "context" );
	Stream_Enable( context->info, False );
	Stream_Enable( context->verbose, False );
	Stream_Enable( context->debug, False );

	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, testPerMaterial );
	stgMainBuildAndInitialise( cf );
	stgMainLoop( cf );
	stgMainDestroy( cf );
}

void DirectorSuite_TestPerMaterial2( DirectorSuiteData* data ) {
	UnderworldContext*		context;
	Stg_ComponentFactory*	cf;
	char							xml_input[PCU_PATH_MAX];

	pcu_filename_input( "testDirectorPerMaterial2.xml", xml_input );
	cf = stgMainInitFromXML( xml_input, MPI_COMM_WORLD, NULL );
	context = (UnderworldContext*)LiveComponentRegister_Get( cf->LCRegister, "context" );
	Stream_Enable( context->info, False );
	Stream_Enable( context->verbose, False );
	Stream_Enable( context->debug, False );

	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, testPerMaterial2 );
	stgMainBuildAndInitialise( cf );
	stgMainLoop( cf );
	stgMainDestroy( cf );
}

void DirectorSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, DirectorSuiteData );
   pcu_suite_setFixtures( suite, DirectorSuite_Setup, DirectorSuite_Teardown );
   pcu_suite_addTest( suite, DirectorSuite_Test ); /* parallel bug */
   pcu_suite_addTest( suite, DirectorSuite_TestRandom ); /* parallel bug */
   pcu_suite_addTest( suite, DirectorSuite_TestPerMaterial );
   pcu_suite_addTest( suite, DirectorSuite_TestPerMaterial2 );
}


