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
**   Tests the OperatorFieldVariableSuite
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

#include "OperatorFieldVariableSuite.h"

typedef struct {
	MPI_Comm	comm;
	unsigned	rank;
	unsigned	nProcs;
} OperatorFieldVariableSuiteData;

/* Simulate Solid Body Rotation */
InterpolationResult OperatorFieldVariableSuite_dummyInterpolateValueAt( void* sdVariable, Coord coord, double* value ) {
	value[0] = -coord[1];
	value[1] =  coord[0];
	value[2] =  coord[2];
	return LOCAL;
}

void OperatorFieldVariableSuite_dummyWithinElement( void* sdVariable, Element_DomainIndex dEl_I, Coord coord, double* value ) { return ; }
double OperatorFieldVariableSuite_dummyGetMinGlobalValue( void* sdVariable ) { return 0; }
double OperatorFieldVariableSuite_dummyGetMaxGlobalValue( void* sdVariable ) { return 1; }
void OperatorFieldVariableSuite_dummyGetMinAndMaxLocalCoords( void* sdVariable, Coord min, Coord max ) {
	min[0] = -1.0;
   min[1] = -2.3;
   min[2] = -4.0;

   max[0] = 11.0;
   max[1] = 12.3;
   max[2] = 14.0;
   return ;
}
void OperatorFieldVariableSuite_dummyGetMinAndMaxGlobalCoords( void* sdVariable, Coord min, Coord max ) {   min[0] = -1.0;
   min[0] = -111.0;
   min[1] = -112.3;
   min[2] = -114.0;

   max[0] = 111.0;
   max[1] = 112.3;
   max[2] = 114.0;
   return ;
}


void OperatorFieldVariableSuite_Setup( OperatorFieldVariableSuiteData* data ) {
	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );
}

void OperatorFieldVariableSuite_Teardown( OperatorFieldVariableSuiteData* data ) {
}

void OperatorFieldVariableSuite_TestOperator( OperatorFieldVariableSuiteData* data ) {
	int							procToWatch;
	Stream*						stream; 	
	FieldVariable_Register*	fV_Register;
   FieldVariable*				velocityField;
   OperatorFieldVariable*	velSquared2D;
   OperatorFieldVariable*	velSquared3D;
   OperatorFieldVariable*	velMag2D;
   OperatorFieldVariable*	velMag3D;
   double						coord[3][3] = {{ 0.4 , 2.0 , 7.0 }, { -0.2 , 6.0 , 2.0 },{ 0.3 , -2.0 , -13.0 }} ;
   double						value[3];
   Index							index;
   Coord                 	min, max;
	char							expected_file[PCU_PATH_MAX];

	Journal_Enable_NamedStream( Info_Type, CartesianGenerator_Type, False );
	stream = Journal_Register( Info_Type, "OperatorFieldVariableStream" );
	Stream_RedirectFile( stream, "operatorFieldVariable.dat" );

	if( data->nProcs >= 2 ) {
		procToWatch = 1;
	}
	else {
		procToWatch = 0;
	}

	if( data->rank == procToWatch ) {
		fV_Register = FieldVariable_Register_New();

		/* constructor  */
		velocityField = _FieldVariable_New(
			sizeof(FieldVariable),
			FieldVariable_Type,
  			_FieldVariable_Delete,
			_FieldVariable_Print,
			_FieldVariable_Copy,
			(Stg_Component_DefaultConstructorFunction*)FieldVariable_DefaultNew,
			_FieldVariable_Construct,
			_FieldVariable_Build,
			_FieldVariable_Initialise,
			_FieldVariable_Execute,
			_FieldVariable_Destroy,
			"Velocity",
			True,
			OperatorFieldVariableSuite_dummyInterpolateValueAt,
			OperatorFieldVariableSuite_dummyGetMinGlobalValue,
			OperatorFieldVariableSuite_dummyGetMaxGlobalValue,
			OperatorFieldVariableSuite_dummyGetMinAndMaxLocalCoords,
			OperatorFieldVariableSuite_dummyGetMinAndMaxGlobalCoords,
			3,
			3,
			False,
         data->comm,
			fV_Register );

		/* Construct Unary Operator Field Variables */
		velSquared2D = OperatorFieldVariable_NewUnary( "VelocitySquaredField2D", velocityField, "VectorSquare" );
 		velSquared2D->_operator->operandDofs = 2;

		velSquared3D = OperatorFieldVariable_NewUnary( "VelocitySquaredField3D", velocityField, "VectorSquare" );

		velMag2D = OperatorFieldVariable_NewUnary( "VelocityMagnitudeField2D", velocityField, "Magnitude" );
		velMag2D->_operator->operandDofs = 2;

		velMag3D = OperatorFieldVariable_NewUnary( "VelocityMagnitudeField3D", velocityField, "Magnitude" );

		/* Do Interpolations */
		for ( index = 0 ; index < 3 ; index++ ) {
			Journal_Printf( stream, "coord = ");
			StGermain_PrintVector( stream, coord[ index ], 3 );

			Journal_Printf( stream, "velocity = ");
			FieldVariable_InterpolateValueAt( velocityField, coord[ index ], value );	
			StGermain_PrintVector( stream, value, 3 );

			Journal_Printf( stream, "velocity squared 2d = ");
			FieldVariable_InterpolateValueAt( velSquared2D, coord[ index ], value );
			StGermain_PrintVector( stream, value, 1 );
			Journal_Printf( stream, "velocity squared 3d = ");
			FieldVariable_InterpolateValueAt( velSquared3D, coord[ index ], value );
			StGermain_PrintVector( stream, value, 1 );

			Journal_Printf( stream, "velocity magnitude 2d = ");
			FieldVariable_InterpolateValueAt( velMag2D, coord[ index ], value );
			StGermain_PrintVector( stream, value, 1 );
			Journal_Printf( stream, "velocity magnitude 3d = ");
			FieldVariable_InterpolateValueAt( velMag3D, coord[ index ], value );
			StGermain_PrintVector( stream, value, 1 );
			Journal_Printf(stream, "\n");
		}

		Journal_Printf( stream , "testing min max local coords:\n" );
		Journal_Printf( stream, "velocity:\n");
		FieldVariable_GetMinAndMaxLocalCoords( velocityField, min, max );
		StGermain_PrintNamedVector( stream, min, 3 );
		StGermain_PrintNamedVector( stream, max, 3 );

		Journal_Printf( stream, "velocity squared 2d = \n");
		FieldVariable_GetMinAndMaxLocalCoords( velSquared2D, min, max );
		StGermain_PrintNamedVector( stream, min, 3 );
		StGermain_PrintNamedVector( stream, max, 3 );

		Journal_Printf( stream , "testing min max global coords:\n" );
		Journal_Printf( stream, "velocity:\n");
		FieldVariable_GetMinAndMaxGlobalCoords( velocityField, min, max );
		StGermain_PrintNamedVector( stream, min, 3 );
		StGermain_PrintNamedVector( stream, max, 3 );

		Journal_Printf( stream, "velocity squared 2d = \n");
		FieldVariable_GetMinAndMaxGlobalCoords( velSquared2D, min, max );
		StGermain_PrintNamedVector( stream, min, 3 );
		StGermain_PrintNamedVector( stream, max, 3 );

		pcu_filename_expected( "testOperatorFieldVariableOutput.expected", expected_file );
		pcu_check_fileEq( "operatorFieldVariable.dat", expected_file );

		/* Clean Up */
		Stg_Class_Delete(fV_Register);
	}
}
	
void OperatorFieldVariableSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, OperatorFieldVariableSuiteData );
	pcu_suite_setFixtures( suite, OperatorFieldVariableSuite_Setup, OperatorFieldVariableSuite_Teardown );
	pcu_suite_addTest( suite, OperatorFieldVariableSuite_TestOperator );
}
