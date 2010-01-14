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
** $Id: Init.c 4160 2007-07-30 06:17:06Z DavidLee $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <mpi.h>

#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>

#include "Utils.h"


Bool StgDomainUtils_Init( int* argc, char** argv[] ) {
	Journal_Printf( Journal_Register( DebugStream_Type, (Name)"Context"  ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */
	
	VariableCondition_Register_Add( variableCondition_Register, AllElementsVC_Type, AllElementsVC_Factory );
	VariableCondition_Register_Add( variableCondition_Register, AllNodesVC_Type, AllNodesVC_Factory );
	VariableCondition_Register_Add( variableCondition_Register, WallVC_Type, WallVC_Factory );
	VariableCondition_Register_Add( variableCondition_Register, ContactVC_Type, ContactVC_Factory );
	VariableCondition_Register_Add( variableCondition_Register, CornerVC_Type, CornerVC_Factory );
	VariableCondition_Register_Add( variableCondition_Register, InnerWallVC_Type, InnerWallVC_Factory );
	VariableCondition_Register_Add( variableCondition_Register, MeshShapeVC_Type, MeshShapeVC_Factory );
	
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), DomainContext_Type, (Name)"0", (void* (*)(Name))_DomainContext_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister( ), AllElementsVC_Type, "0", (void* (*)(Name))_AllElementsVC_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), AllNodesVC_Type, (Name)"0", (void* (*)(Name))_AllNodesVC_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister( ), DofLayout_Type, "0", (void* (*)(Name))_DofLayout_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), FieldVariable_Type, (Name)"0", (void* (*)(Name))_FieldVariable_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister( ), OperatorFieldVariable_Type, "0", (void* (*)(Name))_OperatorFieldVariable_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), WallVC_Type, (Name)"0", (void* (*)(Name))_WallVC_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister( ), ContactVC_Type, "0", (void* (*)(Name))_ContactVC_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), CornerVC_Type, (Name)"0", (void* (*)(Name))_CornerVC_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister( ), InnerWallVC_Type, "0", (void* (*)(Name))_InnerWallVC_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), MeshShapeVC_Type, (Name)"0", (void*  (*)(Name))_MeshShapeVC_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister( ), MeshBoundaryShape_Type, "0", (void*  (*)(Name))MeshBoundaryShape_New );
#ifdef HAVE_PETSC
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), RegularRemesherCmpt_Type, (Name)"0", (void* (*)(Name))_RegularRemesherCmpt_DefaultNew );
#endif
/*
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister( ), Remesher_Type, "0", (void* (*)(Name))_Remesher_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), StripRemesher_Type, (Name)"0", (void* (*)(Name))_StripRemesher_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister( ), CellRemesher_Type, "0", (void* (*)(Name))_CellRemesher_DefaultNew );
*/

	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), TimeIntegrator_Type, (Name)"0", (void*  (*)(Name))_TimeIntegrator_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister( ), TimeIntegrand_Type, "0", (void*  (*)(Name))_TimeIntegrand_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), ShapeAdvector_Type, (Name)"0", (void*  (*)(Name))_ShapeAdvector_DefaultNew  );

	RegisterParent( DomainContext_Type, AbstractContext_Type );
	RegisterParent( Operator_Type, Stg_Object_Type );
	RegisterParent( AllElementsVC_Type, VariableCondition_Type );
	RegisterParent( AllNodesVC_Type, VariableCondition_Type );
	RegisterParent( WallVC_Type, VariableCondition_Type );
	RegisterParent( ContactVC_Type, WallVC_Type );
	RegisterParent( CornerVC_Type, VariableCondition_Type );
	RegisterParent( InnerWallVC_Type, VariableCondition_Type );
	RegisterParent( MeshShapeVC_Type, VariableCondition_Type );
	RegisterParent( MeshBoundaryShape_Type, Stg_Shape_Type );
	RegisterParent( DofLayout_Type, Stg_Component_Type );
#ifdef HAVE_PETSC
	RegisterParent( RegularRemesherCmpt_Type, Remesher_Type );
#endif
/*
	RegisterParent( Remesher_Type, Stg_Component_Type );
	RegisterParent( StripRemesher_Type, Remesher_Type );
	RegisterParent( CellRemesher_Type, Remesher_Type );
*/

	RegisterParent( FieldVariable_Type, Stg_Component_Type );
	RegisterParent( OperatorFieldVariable_Type, FieldVariable_Type );
	RegisterParent( FieldVariable_Register_Type, NamedObject_Register_Type );
	RegisterParent( LinearRegression_Type, Stg_Class_Type );
	RegisterParent( TimeIntegrand_Type, Stg_Component_Type );
	RegisterParent( TimeIntegrator_Type, Stg_Component_Type );
	RegisterParent( ShapeAdvector_Type, Stg_Component_Type );
	
	return True;
}


