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
**   Tests the TimeIntegrationSuite
**
** $Id: testTemplate.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h> 
#include "StgDomain/Geometry/Geometry.h"
#include "StgDomain/Shape/Shape.h"
#include "StgDomain/Mesh/Mesh.h" 
#include "StgDomain/Utils/Utils.h"
#include "StgDomain/Swarm/Swarm.h"

#include "TimeIntegrationSuite.h"

#define  CURR_MODULE_NAME "TimeIntegrationSuite"

typedef struct {
	MPI_Comm	comm;
	unsigned	rank;
	unsigned	nProcs;
} TimeIntegrationSuiteData;

double TimeIntegrationSuite_GetDt( void* context ) {
   return 0.1;
}

Bool TimeIntegrationSuite_ConstantTimeDeriv( void* timeIntegratee, Index array_I, double* timeDeriv ) {
	timeDeriv[0] = 2.0 * array_I;
   timeDeriv[1] = -1.0 * array_I;

   return True;
}
Bool TimeIntegrationSuite_ConstantTimeDeriv2( void* timeIntegratee, Index array_I, double* timeDeriv ) {
   timeDeriv[0] = -0.5 * array_I;
   timeDeriv[1] = 3.0 * array_I;

   return True;
}
Bool TimeIntegrationSuite_LinearTimeDeriv( void* timeIntegratee, Index array_I, double* timeDeriv ) {
   double time = TimeIntegratee_GetTime( timeIntegratee );

   timeDeriv[0] = 2.0 * array_I * time;
   timeDeriv[1] = -1.0 * array_I * time;

   return True;
}
Bool TimeIntegrationSuite_LinearTimeDeriv2( void* timeIntegratee, Index array_I, double* timeDeriv ) {
   double time = TimeIntegratee_GetTime( timeIntegratee );

   timeDeriv[0] = -0.5 * array_I * time;
   timeDeriv[1] = 3.0 * array_I * time;

   return True;
}
Bool TimeIntegrationSuite_CubicTimeDeriv( void* timeIntegratee, Index array_I, double* timeDeriv ) {
   double time = TimeIntegratee_GetTime( timeIntegratee );

   timeDeriv[0] = 2.0 * array_I * ( time * time * time - time*time );
   timeDeriv[1] = -1.0 * array_I * ( time * time * time - time*time );

   return True;
}
Bool TimeIntegrationSuite_CubicTimeDeriv2( void* timeIntegratee, Index array_I, double* timeDeriv ) {
   double time = TimeIntegratee_GetTime( timeIntegratee );

   timeDeriv[0] = -0.5 * array_I * ( time * time * time - time*time );
   timeDeriv[1] = 3.0 * array_I * ( time * time * time - time*time );

   return True;
}

TimeIntegratee_CalculateTimeDerivFunction* TimeIntegrationSuite_GetFunctionPtr( Name derivName ) {
   if ( strcasecmp( derivName, "Linear" ) == 0 )
      return TimeIntegrationSuite_LinearTimeDeriv;
   else if ( strcasecmp( derivName, "Linear2" ) == 0 )
      return TimeIntegrationSuite_LinearTimeDeriv2;
   else if ( strcasecmp( derivName, "Cubic" ) == 0 )
      return TimeIntegrationSuite_CubicTimeDeriv;
   else if ( strcasecmp( derivName, "Cubic2" ) == 0 )
      return TimeIntegrationSuite_CubicTimeDeriv2;
   else if ( strcasecmp( derivName, "Constant" ) == 0 )
      return TimeIntegrationSuite_ConstantTimeDeriv;
   else if ( strcasecmp( derivName, "Constant2" ) == 0 )
      return TimeIntegrationSuite_ConstantTimeDeriv2;
   else
      Journal_Firewall( 0 , Journal_Register( Error_Type, CURR_MODULE_NAME ), "Don't understand DerivName '%s'\n", derivName  );

   return NULL;
}

void TimeIntegrationSuite_TestContextType( void* timeIntegratee, Stg_Class* data ) {
   Stream* stream = Journal_Register (Info_Type, "myStream");

   Journal_Printf( stream, "In func %s\n", __func__ );
   assert( data->type == DomainContext_Type );
}
void TimeIntegrationSuite_TestVariableType( void* timeIntegratee, Stg_Class* data ) {
   Stream* stream = Journal_Register (Info_Type, "myStream");

   Journal_Printf( stream, "In func %s\n", __func__ );
   assert( data->type == Variable_Type );
}


void TimeIntegrationSuite_Setup( TimeIntegrationSuiteData* data ) {
	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );
}

void TimeIntegrationSuite_Teardown( TimeIntegrationSuiteData* data ) {
}

void TimeIntegrationSuite_TestDriver( TimeIntegrationSuiteData* data, char *_name, char *_DerivName0, char *_DerivName1, int _order   ) {
	int					procToWatch;
	Stream*				stream;
	XML_IO_Handler*	ioHandler;
   Dictionary*			dictionary;
   TimeIntegrator*	timeIntegrator;
   TimeIntegratee*	timeIntegratee;
   TimeIntegratee*	timeIntegrateeList[2];
   DomainContext*		context;
   Variable*			variable;
   Variable*			variableList[2];
   double*				array;
   double*				array2;
   Index					size0	= 11;
   Index					size1	= 7;
	Index					array_I;
   Index					timestep	= 0;
   Index					maxTimesteps	= 10;
	Bool					simultaneous;
   unsigned				order;
   double				error	= 0.0;
   Name					derivName;
   double				tolerance	= 0.001;
   Index					integratee_I;
   Index					integrateeCount	= 2;
	char					expected_file[PCU_PATH_MAX];

	stream = Journal_Register( Info_Type, "EulerStream" );
	Stream_RedirectFile( stream, _name );

	if( data->nProcs >= 2 ) {
		procToWatch = 1;
	}
	else {
		procToWatch = 0;
	}

	if( data->rank == procToWatch ) {
		/* Create Context */
		dictionary = Dictionary_New();
		Dictionary_Add(dictionary, "outputPath", Dictionary_Entry_Value_FromString("./output"));
		Dictionary_Add(dictionary, "DerivName0", Dictionary_Entry_Value_FromString(_DerivName0));
		Dictionary_Add(dictionary, "DerivName1", Dictionary_Entry_Value_FromString(_DerivName1));
		context = _DomainContext_New(
			sizeof(DomainContext),
			DomainContext_Type,
			_DomainContext_Delete,
			_DomainContext_Print,
			NULL,
			NULL,
			NULL,
			NULL,
			NULL,
			NULL,
			NULL,
			"discretisationContext",
			True,
			NULL,
			0,0,
			data->comm, dictionary );
		ContextEP_Append( context, AbstractContext_EP_Dt, TimeIntegrationSuite_GetDt );

		Journal_Printf( stream, "!!! info %d\n", Stream_IsEnable( Journal_Register( Info_Type, "TimeIntegrator" ) ) );

		/* Create Stuff */
		order							= _order;
		simultaneous				= False;
		variableList[0]			= Variable_NewVector( "testVariable",  Variable_DataType_Double, 2, &size0, NULL, (void**)&array, NULL );
		variableList[1]			= Variable_NewVector( "testVariable2", Variable_DataType_Double, 2, &size1, NULL, (void**)&array2, NULL );
		timeIntegrator				= TimeIntegrator_New( "testTimeIntegrator", order, simultaneous, NULL, NULL );
		timeIntegrateeList[0]	= TimeIntegratee_New( "testTimeIntegratee0", timeIntegrator, variableList[0], 0, NULL, True );
		timeIntegrateeList[1]	= TimeIntegratee_New( "testTimeIntegratee1", timeIntegrator, variableList[1], 0, NULL, True );

		derivName = Dictionary_GetString( dictionary, "DerivName0" );
		timeIntegrateeList[0]->_calculateTimeDeriv = TimeIntegrationSuite_GetFunctionPtr( derivName );
		Journal_Printf( stream, "DerivName0 - %s\n", derivName );
		derivName = Dictionary_GetString( dictionary, "DerivName1" );
		timeIntegrateeList[1]->_calculateTimeDeriv = TimeIntegrationSuite_GetFunctionPtr( derivName );
		Journal_Printf( stream, "DerivName1 - %s\n", derivName );

		/* Print Stuff to file */
		Journal_PrintValue( stream, order );
		Journal_PrintBool( stream, simultaneous );

		/* Add stuff to EPs */
		TimeIntegrator_AppendSetupEP( timeIntegrator, "start1", TimeIntegrationSuite_TestContextType, CURR_MODULE_NAME, context );
		TimeIntegrator_AppendFinishEP( timeIntegrator, "finish1", TimeIntegrationSuite_TestVariableType, CURR_MODULE_NAME, variableList[0] );
		TimeIntegrator_PrependSetupEP( timeIntegrator, "start0", TimeIntegrationSuite_TestVariableType, CURR_MODULE_NAME, variableList[0] );
		TimeIntegrator_PrependFinishEP( timeIntegrator, "finish0", TimeIntegrationSuite_TestContextType, CURR_MODULE_NAME, context );

		/* Build */
		Stg_Component_Build( variableList[0], context, False );
		Stg_Component_Build( variableList[1], context, False );
		Stg_Component_Build( timeIntegrator, context, False );
		Stg_Component_Build( timeIntegrateeList[0], context, False );
		Stg_Component_Build( timeIntegrateeList[1], context, False );
		array = Memory_Alloc_Array( double, 2 * size0, "name" );
		array2 = Memory_Alloc_Array( double, 2 * size1, "name" );

		/* Initialise */
		memset( array, 0, sizeof(double) * 2 * size0 );
		memset( array2, 0, sizeof(double) * 2 * size1 );
		Stg_Component_Initialise( timeIntegrator, context, False );
		Stg_Component_Initialise( variableList[0], context, False );
		Stg_Component_Initialise( variableList[1], context, False );
		Stg_Component_Initialise( timeIntegrateeList[0], context, False );
		Stg_Component_Initialise( timeIntegrateeList[1], context, False );

		for ( timestep = 0.0 ; timestep < maxTimesteps ; timestep ++ ) {
			Journal_Printf( stream, "Step %u - Time = %.3g\n", timestep, context->currentTime );

			Stg_Component_Execute( timeIntegrator, context, True );
			context->currentTime += AbstractContext_Dt( context );

			for ( integratee_I = 0 ; integratee_I < integrateeCount ; integratee_I++ ) {
				timeIntegratee	= timeIntegrateeList[ integratee_I ];
				variable			= variableList[ integratee_I ];
				for ( array_I = 0 ; array_I < variable->arraySize ; array_I++ ) {
					if ( timeIntegratee->_calculateTimeDeriv == TimeIntegrationSuite_ConstantTimeDeriv ) {
						error += fabs( Variable_GetValueAtDouble( variable, array_I, 0 ) - 2.0 * array_I * context->currentTime );
						error += fabs( Variable_GetValueAtDouble( variable, array_I, 1 ) + array_I * context->currentTime );
					}
					else if ( timeIntegratee->_calculateTimeDeriv == TimeIntegrationSuite_ConstantTimeDeriv2 ) {
						error += fabs( Variable_GetValueAtDouble( variable, array_I, 0 ) + 0.5 * array_I * context->currentTime );
						error += fabs( Variable_GetValueAtDouble( variable, array_I, 1 ) - 3 * array_I * context->currentTime );
					}
					else if ( timeIntegratee->_calculateTimeDeriv == TimeIntegrationSuite_LinearTimeDeriv ) {
						error += fabs( Variable_GetValueAtDouble( variable, array_I, 0 ) - array_I * context->currentTime * context->currentTime );
						error += fabs( Variable_GetValueAtDouble( variable, array_I, 1 ) + 0.5 * array_I * context->currentTime * context->currentTime );
					}
					else if ( timeIntegratee->_calculateTimeDeriv == TimeIntegrationSuite_LinearTimeDeriv2 ) {
							error += fabs( Variable_GetValueAtDouble( variable, array_I, 0 ) + 0.25 * array_I * context->currentTime * context->currentTime );
							error += fabs( Variable_GetValueAtDouble( variable, array_I, 1 ) - 1.5 * array_I * context->currentTime * context->currentTime );
					}
					else if ( timeIntegratee->_calculateTimeDeriv == TimeIntegrationSuite_CubicTimeDeriv ) {
							error += fabs( Variable_GetValueAtDouble( variable, array_I, 0 ) - 2.0 * array_I * ( 0.25 * pow( context->currentTime, 4.0 ) - pow( context->currentTime, 3.0)/3.0));
							error += fabs( Variable_GetValueAtDouble( variable, array_I, 1 ) + array_I * ( 0.25 * pow( context->currentTime, 4.0 ) - pow( context->currentTime, 3.0 )/3.0));
					}
					else if ( timeIntegratee->_calculateTimeDeriv == TimeIntegrationSuite_CubicTimeDeriv2 ) {
							error += fabs( Variable_GetValueAtDouble( variable, array_I, 0 ) + 0.5 * array_I * ( 0.25 * pow( context->currentTime, 4.0 ) - pow( context->currentTime, 3.0)/3.0));
							error += fabs( Variable_GetValueAtDouble( variable, array_I, 1 ) - 3.0 * array_I * ( 0.25 * pow( context->currentTime, 4.0 ) - pow( context->currentTime, 3.0 )/3.0));
					}
					else
						Journal_Firewall( 0 , Journal_Register( Error_Type, CURR_MODULE_NAME ), "Don't understand _calculateTimeDeriv = %p\n", timeIntegratee->_calculateTimeDeriv );
				}
			}
		}
		pcu_check_lt( error, tolerance );

		if ( error < tolerance )
			Journal_Printf( stream, "Passed\n" );
		else
			Journal_Printf( stream, "Failed - Error = %lf\n", error );

		if ( timeIntegratee->_calculateTimeDeriv == TimeIntegrationSuite_ConstantTimeDeriv
			|| timeIntegratee->_calculateTimeDeriv == TimeIntegrationSuite_ConstantTimeDeriv2 ) {
			pcu_filename_expected( "testTimeIntegrationEulerOutput.expected", expected_file );
		}
		else if ( timeIntegratee->_calculateTimeDeriv == TimeIntegrationSuite_LinearTimeDeriv
			|| timeIntegratee->_calculateTimeDeriv == TimeIntegrationSuite_LinearTimeDeriv2 ) {
			pcu_filename_expected( "testTimeIntegrationRK2Output.expected", expected_file );
		}
		else if ( timeIntegratee->_calculateTimeDeriv == TimeIntegrationSuite_CubicTimeDeriv
			|| timeIntegratee->_calculateTimeDeriv == TimeIntegrationSuite_CubicTimeDeriv2 ) {
			pcu_filename_expected( "testTimeIntegrationRK4Output.expected", expected_file );
		}

		pcu_check_fileEq( _name, expected_file );

		/* Destroy stuff */
		Memory_Free( array );
		Memory_Free( array2 );
		Stg_Class_Delete( variable );
		Stg_Class_Delete( timeIntegrator );
		Stg_Class_Delete( timeIntegrateeList[0] );
		Stg_Class_Delete( timeIntegrateeList[1] );
		remove( _name );
	}
}
	
void TimeIntegrationSuite_TestEuler( TimeIntegrationSuiteData* data ) {
	TimeIntegrationSuite_TestDriver( data, "testIntegrationEuler", "Constant", "Constant2", 1 );
}
void TimeIntegrationSuite_TestRK2( TimeIntegrationSuiteData* data ) {
	TimeIntegrationSuite_TestDriver( data, "testIntegrationRK2", "Linear", "Linear2", 2 );
}
void TimeIntegrationSuite_TestRK4( TimeIntegrationSuiteData* data ) {
	TimeIntegrationSuite_TestDriver( data, "testIntegrationRK4", "Cubic", "Cubic2", 4 );
}

void TimeIntegrationSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, TimeIntegrationSuiteData );
	pcu_suite_setFixtures( suite, TimeIntegrationSuite_Setup, TimeIntegrationSuite_Teardown );
	pcu_suite_addTest( suite, TimeIntegrationSuite_TestEuler );
	pcu_suite_addTest( suite, TimeIntegrationSuite_TestRK2 );
	pcu_suite_addTest( suite, TimeIntegrationSuite_TestRK4 );
}
