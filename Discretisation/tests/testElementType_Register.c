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
** $Id: testElementType_Register.c 832 2007-05-16 01:11:18Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include "StgFEM/Discretisation/Discretisation.h"

#include <stdio.h>
#include <stdlib.h>
#include <petsc.h>

int main( int argc, char* argv[] ) {
	MPI_Comm			CommWorld;
	int				rank;
	int				numProcessors;
	int				procToWatch;
	ElementType_Register*		elementType_Register;
	Stream*				stream;
	
	/* Initialise MPI, get world info */
	MPI_Init( &argc, &argv );
	MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
	MPI_Comm_size( CommWorld, &numProcessors );
	MPI_Comm_rank( CommWorld, &rank );
	
	StGermain_Init( &argc, &argv );
	StgFEM_Discretisation_Init( &argc, &argv );
	MPI_Barrier( CommWorld ); /* Ensures copyright info always come first in output */
	
	stream = Journal_Register (Info_Type, "myStream");

	if( argc >= 2 ) {
		procToWatch = atoi( argv[1] );
	}
	else {
		procToWatch = 0;
	}
	if( rank == procToWatch ) printf( "Watching rank: %i\n", rank );
	
	/* Build the element type register */
	elementType_Register = ElementType_Register_New("elementTypeRegister");
	
	/* Add the prebuilt types */
	ElementType_Register_Add( elementType_Register, (ElementType*)ConstantElementType_New("constant") );
	ElementType_Register_Add( elementType_Register, (ElementType*)BilinearElementType_New("bilinear") );
	ElementType_Register_Add( elementType_Register, (ElementType*)TrilinearElementType_New("trilinear") );
	
	/* Manually create extra types to force test the list re-sizing */
	ElementType_Register_Add( elementType_Register, _ElementType_New( sizeof(ConstantElementType), "TestElementType0",
		_ConstantElementType_Delete, _ConstantElementType_Print, NULL, ConstantElementType_DefaultNew,
		_ConstantElementType_Construct, _ConstantElementType_Build, _ConstantElementType_Initialise,
		_ConstantElementType_Execute, _ConstantElementType_Destroy, "TestElementType0Name", True,
		_ConstantElementType_SF_allNodes, _ConstantElementType_SF_allLocalDerivs_allNodes,
		_ConstantElementType_ConvertGlobalCoordToElLocal,
		1 ) );
	
	ElementType_Register_Add( elementType_Register, _ElementType_New( sizeof(ConstantElementType), "TestElementType1",
		_ConstantElementType_Delete, _ConstantElementType_Print, NULL, ConstantElementType_DefaultNew,
		_ConstantElementType_Construct, _ConstantElementType_Build, _ConstantElementType_Initialise,
		_ConstantElementType_Execute, _ConstantElementType_Destroy, "TestElementType1Name", True,
		_ConstantElementType_SF_allNodes, _ConstantElementType_SF_allLocalDerivs_allNodes,
		_ConstantElementType_ConvertGlobalCoordToElLocal,
		1 ) );
	
	ElementType_Register_Add( elementType_Register, _ElementType_New( sizeof(ConstantElementType), "TestElementType2",
		_ConstantElementType_Delete, _ConstantElementType_Print, NULL, ConstantElementType_DefaultNew,
		_ConstantElementType_Construct, _ConstantElementType_Build, _ConstantElementType_Initialise,
		_ConstantElementType_Execute, _ConstantElementType_Destroy, "TestElementType2Name", True,
		_ConstantElementType_SF_allNodes, _ConstantElementType_SF_allLocalDerivs_allNodes,
		_ConstantElementType_ConvertGlobalCoordToElLocal,
		1 ) );
	
	ElementType_Register_Add( elementType_Register, _ElementType_New( sizeof(ConstantElementType), "TestElementType3",
		_ConstantElementType_Delete, _ConstantElementType_Print, NULL, ConstantElementType_DefaultNew,
		_ConstantElementType_Construct, _ConstantElementType_Build, _ConstantElementType_Initialise,
		_ConstantElementType_Execute, _ConstantElementType_Destroy, "TestElementType3Name", True,
		_ConstantElementType_SF_allNodes, _ConstantElementType_SF_allLocalDerivs_allNodes,
		_ConstantElementType_ConvertGlobalCoordToElLocal,
		1 ) );
	
	ElementType_Register_Add( elementType_Register, _ElementType_New( sizeof(ConstantElementType), "TestElementType4",
		_ConstantElementType_Delete, _ConstantElementType_Print, NULL, ConstantElementType_DefaultNew,
		_ConstantElementType_Construct, _ConstantElementType_Build, _ConstantElementType_Initialise,
		_ConstantElementType_Execute, _ConstantElementType_Destroy, "TestElementType4Name", True,
		_ConstantElementType_SF_allNodes, _ConstantElementType_SF_allLocalDerivs_allNodes,
		_ConstantElementType_ConvertGlobalCoordToElLocal,
		1 ) );
	
	ElementType_Register_Add( elementType_Register, _ElementType_New( sizeof(ConstantElementType), "TestElementType5",
		_ConstantElementType_Delete, _ConstantElementType_Print, NULL, ConstantElementType_DefaultNew,
		_ConstantElementType_Construct, _ConstantElementType_Build, _ConstantElementType_Initialise,
		_ConstantElementType_Execute, _ConstantElementType_Destroy, "TestElementType5Name", True,
		_ConstantElementType_SF_allNodes, _ConstantElementType_SF_allLocalDerivs_allNodes,
		_ConstantElementType_ConvertGlobalCoordToElLocal,
		1 ) );
	
	
	Stg_Class_Print( elementType_Register, stream );
	
	/* Destroy stuff */
	Stg_Class_Delete( elementType_Register );
	
	StgFEM_Discretisation_Finalise();
	StGermain_Finalise();
	
	/* Close off MPI */
	MPI_Finalize();
	
	return 0; /* success */
}

