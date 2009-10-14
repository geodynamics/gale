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

	test( context );

	stgMainDestroy( cf );
}


void DirectorSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, DirectorSuiteData );
   pcu_suite_setFixtures( suite, DirectorSuite_Setup, DirectorSuite_Teardown );
   pcu_suite_addTest( suite, DirectorSuite_Test );
}
