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
** $Id: Init.c 964 2007-10-11 08:03:06Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "StgFEM/Discretisation/Discretisation.h"
#include "SystemSetup.h"


Stream* StgFEM_SLE_Debug = NULL;
Stream* StgFEM_SLE_SystemSetup_Debug = NULL;

/** Initialises the Linear Algebra package, then any init for this package
such as streams etc */
Bool StgFEM_SLE_SystemSetup_Init( int* argc, char** argv[] ) {
	Journal_Printf( Journal_Register( DebugStream_Type, "Context" ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */
	
	/* initialise this level's streams */
	StgFEM_SLE_Debug = Stream_RegisterChild( StgFEM_Debug, "SLE" );
	StgFEM_SLE_SystemSetup_Debug = Stream_RegisterChild( StgFEM_SLE_Debug, "SystemSetup" );
	
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), FiniteElementContext_Type, "0", FiniteElementContext_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), ForceVector_Type, "0", _ForceVector_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), SolutionVector_Type, "0", SolutionVector_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), StiffnessMatrix_Type, "0", StiffnessMatrix_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), StiffnessMatrixTerm_Type, "0", _StiffnessMatrixTerm_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), SystemLinearEquations_Type, "0", _SystemLinearEquations_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), ForceTerm_Type, "0", _ForceTerm_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), 
				   MultigridSolver_Type, "0", 
				   (Stg_Component_DefaultConstructorFunction*)MultigridSolver_New );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), 
				   SROpGenerator_Type, "0", 
				   (Stg_Component_DefaultConstructorFunction*)SROpGenerator_New );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), 
				   PETScMGSolver_Type, "0", 
				   (Stg_Component_DefaultConstructorFunction*)PETScMGSolver_New );

	RegisterParent( SystemLinearEquations_Type,    Stg_Component_Type );
	RegisterParent( SLE_Solver_Type,               Stg_Component_Type );
	RegisterParent( StiffnessMatrix_Type,          Stg_Component_Type );
	RegisterParent( StiffnessMatrixTerm_Type,      Stg_Component_Type );
	RegisterParent( SolutionVector_Type,           Stg_Component_Type );
	RegisterParent( ForceVector_Type,              SolutionVector_Type );
	RegisterParent( ForceTerm_Type,                Stg_Component_Type );
	RegisterParent( Assembler_Type, Stg_Class_Type );
	RegisterParent( FiniteElementContext_Type,     DomainContext_Type );
	RegisterParent( PETScMGSolver_Type, Stg_Component_Type );
	RegisterParent( MultigridSolver_Type, Stg_Component_Type );
	RegisterParent( MGOpGenerator_Type, Stg_Component_Type );
	RegisterParent( SROpGenerator_Type, MGOpGenerator_Type);

	return True;
}
