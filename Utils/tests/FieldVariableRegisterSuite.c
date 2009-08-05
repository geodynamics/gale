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
**   Tests the FieldVariableRegisterSuite
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

#include "FieldVariableRegisterSuite.h"

typedef struct {
	MPI_Comm	comm;
	unsigned	rank;
	unsigned	nProcs;
} FieldVariableRegisterSuiteData;

InterpolationResult FieldVariableRegisterSuite_dummyInterpolateValueAt( void* sdVariable, Coord coord, double* value ) { return OUTSIDE_GLOBAL; }
void FieldVariableRegisterSuite_dummyWithinElement( void* sdVariable, Element_DomainIndex dEl_I, Coord coord, double* value ) { return ; }
double FieldVariableRegisterSuite_dummyGetMinGlobalValue( void* sdVariable ) { return 0; }
double FieldVariableRegisterSuite_dummyGetMaxGlobalValue( void* sdVariable ) { return 1; }
void FieldVariableRegisterSuite_dummyGetMinAndMaxLocalCoords( void* sdVariable, Coord min, Coord max ) { return ; }
void FieldVariableRegisterSuite_dummyGetMinAndMaxGlobalCoords( void* sdVariable, Coord min, Coord max ) { return ; }

void FieldVariableRegisterSuite_Setup( FieldVariableRegisterSuiteData* data ) {
	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );
}

void FieldVariableRegisterSuite_Teardown( FieldVariableRegisterSuiteData* data ) {
}

void FieldVariableRegisterSuite_TestRegister( FieldVariableRegisterSuiteData* data ) {
	int							procToWatch;
	Stream*						stream;
   FieldVariable_Register*	fV_Register;
   FieldVariable*				testFVs[] = { NULL, NULL, NULL };
   Name							fvNames[] = { "testFV1", "testFV2", "testFV3" };
   Index							ii;
   Index							fV_Index;
	char							expected_file[PCU_PATH_MAX];

	Journal_Enable_NamedStream( Info_Type, CartesianGenerator_Type, False );
	stream = Journal_Register( Info_Type, "FieldVariableRegisterStream" );

	if( data->nProcs >= 2 ) {
		procToWatch = 1;
	}
	else {
		procToWatch = 0;
	}

	if( data->rank == procToWatch ) {
		fV_Register = FieldVariable_Register_New();

		/* abstract class, so have to hack constructor  */
		for( ii=0; ii < 3; ii++ ) {
			testFVs[ii] = _FieldVariable_New( sizeof(FieldVariable), FieldVariable_Type, _Stg_Component_Delete,
								_Stg_Component_Print, NULL, (Stg_Component_DefaultConstructorFunction*)FieldVariable_DefaultNew, _Variable_Construct, _Variable_Build,
								_Variable_Initialise, _Variable_Execute, _Variable_Destroy, fvNames[ii], True,
								FieldVariableRegisterSuite_dummyInterpolateValueAt,
								FieldVariableRegisterSuite_dummyGetMinGlobalValue, FieldVariableRegisterSuite_dummyGetMaxGlobalValue,
								FieldVariableRegisterSuite_dummyGetMinAndMaxLocalCoords, FieldVariableRegisterSuite_dummyGetMinAndMaxGlobalCoords,
								0, /* field components */
								3, /* dim */
								False,
								data->comm,
								fV_Register );
		}

		/* test the get functions */
		for( ii=0; ii < 3; ii++ ) {
			fV_Index = FieldVariable_Register_GetIndex( fV_Register, fvNames[ii] );
			pcu_check_true( fV_Index == ii );
			pcu_check_streq( (FieldVariable_Register_GetByIndex( fV_Register, fV_Index ) )->name, fvNames[ii] );
			pcu_check_streq( (FieldVariable_Register_GetByName( fV_Register, fvNames[ii] ) )->name, fvNames[ii] );
		}

		Stg_Class_Delete(fV_Register);
	}
}
	
void FieldVariableRegisterSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, FieldVariableRegisterSuiteData );
	pcu_suite_setFixtures( suite, FieldVariableRegisterSuite_Setup, FieldVariableRegisterSuite_Teardown );
	pcu_suite_addTest( suite, FieldVariableRegisterSuite_TestRegister );
}
