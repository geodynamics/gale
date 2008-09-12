/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
** $Id: Init.c 3901 2006-12-05 22:10:22Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <assert.h>

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "Discretisation/Discretisation.h"
#include "LinearAlgebra.h"


Bool StgFEM_SLE_LinearAlgebra_Init( int* argc, char** argv[] ) {
	int tmp;

	Journal_Printf( Journal_Register( DebugStream_Type, "Context" ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */
	tmp = Stream_GetPrintingRank( Journal_Register( InfoStream_Type, "Context" ) );
	Stream_SetPrintingRank( Journal_Register( InfoStream_Type, "Context" ), 0 );
	Journal_Printf( /* DO NOT CHANGE OR REMOVE */
		Journal_Register( InfoStream_Type, "Context" ), 
		"StGermain PETSc-LinearAlgebra Interface revision %s. Copyright (C) 2003-2005 VPAC.\n", VERSION );
	Stream_Flush( Journal_Register( InfoStream_Type, "Context" ) );
	Stream_SetPrintingRank( Journal_Register( InfoStream_Type, "Context" ), tmp );

	RegisterParent( Vector_Type, Stg_Component_Type );
	RegisterParent( Matrix_Type, Stg_Component_Type );
	RegisterParent( MatrixSolver_Type, Stg_Component_Type );

#ifdef HAVE_PETSC
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), 
				   PETScVector_Type, "0", (Stg_Component_DefaultConstructorFunction*)PETScVector_New );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), 
				   PETScMatrix_Type, "0", (Stg_Component_DefaultConstructorFunction*)PETScMatrix_New );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), 
				   PETScMatrixSolver_Type, "0", (Stg_Component_DefaultConstructorFunction*)PETScMatrixSolver_New );

	RegisterParent( PETScVector_Type, Vector_Type );
	RegisterParent( PETScMatrix_Type, Matrix_Type );
	RegisterParent( PETScMatrixSolver_Type, MatrixSolver_Type );

	{
		PetscErrorCode	ec;

		ec = PetscInitialize( argc, argv, (char*)0, NULL );
		CheckPETScError( ec );
	}
#endif

	return True;
}
