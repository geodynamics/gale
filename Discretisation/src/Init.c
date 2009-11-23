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
** $Id: Init.c 1218 2008-09-04 06:18:44Z DavidLee $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <stdio.h>

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "Discretisation.h"


Stream* StgFEM_Debug = NULL;
Stream* StgFEM_Warning = NULL;
Stream* StgFEM_Discretisation_Debug = NULL;


/** Initialises the Linear Algebra package, then any init for this package
such as streams etc */
Bool StgFEM_Discretisation_Init( int* argc, char** argv[] ) {
	Stg_ComponentRegister*          componentRegister = Stg_ComponentRegister_Get_ComponentRegister();
	int tmp;
	
	Journal_Printf( Journal_Register( DebugStream_Type, "Context" ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */
	tmp = Stream_GetPrintingRank( Journal_Register( InfoStream_Type, "Context" ) );
	Stream_SetPrintingRank( Journal_Register( InfoStream_Type, "Context" ), 0 );
	Journal_Printf( /* DO NOT CHANGE OR REMOVE */
		Journal_Register( InfoStream_Type, "Context" ), 
		"StGermain FEM Discretisation Framework revision %s. Copyright (C) 2003-2005 VPAC.\n", VERSION );
	Stream_Flush( Journal_Register( InfoStream_Type, "Context" ) );
	Stream_SetPrintingRank( Journal_Register( InfoStream_Type, "Context" ), tmp );
	
	/* initialise this level's streams */
	StgFEM_Debug = Journal_Register( DebugStream_Type, "StgFEM" );
	StgFEM_Discretisation_Debug = Stream_RegisterChild( StgFEM_Debug, "Discretisation" );
	StgFEM_Warning = Journal_Register( ErrorStream_Type, "StgFEM" );
	
	Stg_ComponentRegister_Add( componentRegister, FeVariable_Type,         "0", FeVariable_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, LinkedDofInfo_Type,      "0", LinkedDofInfo_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, OperatorFeVariable_Type, "0", OperatorFeVariable_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, ShapeFeVariable_Type,    "0", ShapeFeVariable_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, FeSwarmVariable_Type,    "0", _FeSwarmVariable_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, FeMesh_Type, "0", (Stg_Component_DefaultConstructorFunction *)FeMesh_New );
	Stg_ComponentRegister_Add( componentRegister, C0Generator_Type, "0", (Stg_Component_DefaultConstructorFunction *)C0Generator_New );
	Stg_ComponentRegister_Add( componentRegister, C2Generator_Type, "0", (Stg_Component_DefaultConstructorFunction *)C2Generator_New );
/*
	Stg_ComponentRegister_Add( componentRegister, P1Generator_Type, "0", P1Generator_New );
*/
	Stg_ComponentRegister_Add( componentRegister, Inner2DGenerator_Type, "0", (Stg_Component_DefaultConstructorFunction *)Inner2DGenerator_New );
	Stg_ComponentRegister_Add( componentRegister, FieldTest_Type, "0", _FieldTest_DefaultNew );
	
	/** Register Parents for type checking */
	RegisterParent( FeMesh_Algorithms_Type, Mesh_Algorithms_Type );
	RegisterParent( FeMesh_ElementType_Type, Mesh_HexType_Type );
	RegisterParent( ElementType_Type,                  Stg_Component_Type );
	RegisterParent( BilinearElementType_Type,          ElementType_Type );
	RegisterParent( TrilinearElementType_Type,         ElementType_Type );
	RegisterParent( Biquadratic_Type, 		Biquadratic_Type );
	/* i'm assuming this is ok - doesn't seem to complain without it though - dave, 29.05.07 */
	RegisterParent( Triquadratic_Type,		Triquadratic_Type );
	RegisterParent( P1_Type, 			P1_Type );
	RegisterParent( RegularTrilinear_Type,			TrilinearElementType_Type );
	RegisterParent( ConstantElementType_Type,          ElementType_Type );
	RegisterParent( LinearTriangleElementType_Type,    ElementType_Type );
	RegisterParent( ElementType_Register_Type,         Stg_Component_Type );

	RegisterParent( FeEquationNumber_Type,             Stg_Component_Type );
	RegisterParent( LinkedDofInfo_Type,                Stg_Component_Type );
	RegisterParent( FeMesh_Type, Mesh_Type );
	RegisterParent( C0Generator_Type, MeshGenerator_Type );
	RegisterParent( C2Generator_Type, CartesianGenerator_Type );
/*
	RegisterParent( P1Generator_Type, MeshGenerator_Type );
*/
	RegisterParent( Inner2DGenerator_Type, MeshGenerator_Type );
	
	RegisterParent( FeVariable_Type,                   FieldVariable_Type );
	RegisterParent( OperatorFeVariable_Type,           FeVariable_Type );
	RegisterParent( ShapeFeVariable_Type,              FeVariable_Type );
	RegisterParent( FeSwarmVariable_Type,              SwarmVariable_Type );

	RegisterParent( FieldTest_Type,			   Stg_Component_Type );

	{
		PetscErrorCode	ec;
		ec = PetscInitialize( argc, argv, (char*)0, NULL );
		CheckPETScError( ec );
	}

	return True;
}
