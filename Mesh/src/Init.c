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
** $Id: Init.c 4081 2007-04-27 06:20:07Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <mpi.h>
#include <StGermain/StGermain.h>

#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include "Mesh.h"


Stream* Mesh_VerboseConfig = NULL;
Stream* Mesh_Debug = NULL;
Stream* Mesh_Warning = NULL;
Stream* Mesh_Error = NULL;


Bool StgDomainMesh_Init( int* argc, char** argv[] ) {
	Mesh_VerboseConfig = Journal_Register( Info_Type, "Mesh_VerboseConfig" );
	Mesh_Debug = Journal_Register( Debug_Type, "Mesh" );
	Mesh_Warning = Journal_Register( Error_Type, "Mesh" );
	Mesh_Error = Journal_Register( Error_Type, "Mesh" );

	Journal_Printf( Journal_Register( DebugStream_Type, "Context" ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */

	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), 
				   Mesh_Algorithms_Type, "0", (Stg_Component_DefaultConstructorFunction*)Mesh_Algorithms_New );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), 
				   Mesh_HexAlgorithms_Type, "0", 
				   (Stg_Component_DefaultConstructorFunction*)Mesh_HexAlgorithms_New );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), 
				   Mesh_CentroidAlgorithms_Type, "0", 
				   (Stg_Component_DefaultConstructorFunction*)Mesh_CentroidAlgorithms_New );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), 
				   Mesh_RegularAlgorithms_Type, "0", 
				   (Stg_Component_DefaultConstructorFunction*)Mesh_RegularAlgorithms_New );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), 
				   MeshTopology_Type, "0", (Stg_Component_DefaultConstructorFunction*)MeshTopology_New );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), 
				   CartesianGenerator_Type, "0", (Stg_Component_DefaultConstructorFunction*)CartesianGenerator_New );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), 
				   Mesh_Type, "0", (Stg_Component_DefaultConstructorFunction*)Mesh_New );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), 
				   SurfaceAdaptor_Type, "0", (Stg_Component_DefaultConstructorFunction*)SurfaceAdaptor_New );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), 
				   CompressionAdaptor_Type, "0", (Stg_Component_DefaultConstructorFunction*)CompressionAdaptor_New );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), 
				   LinearSpaceAdaptor_Type, "0", (Stg_Component_DefaultConstructorFunction*)LinearSpaceAdaptor_New );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), 
				   MeshVariable_Type, "0", (Stg_Component_DefaultConstructorFunction*)MeshVariable_New );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), 
				   Remesher_Type, "0", (Stg_Component_DefaultConstructorFunction*)_Remesher_DefaultNew );

	RegisterParent( Mesh_ElementType_Type, Stg_Class_Type );
	RegisterParent( Mesh_HexType_Type, Mesh_ElementType_Type );
	RegisterParent( Mesh_CentroidType_Type, Mesh_ElementType_Type );
	RegisterParent( Mesh_Algorithms_Type, Stg_Component_Type );
	RegisterParent( Mesh_HexAlgorithms_Type, Mesh_Algorithms_Type );
	RegisterParent( Mesh_CentroidAlgorithms_Type, Mesh_Algorithms_Type );
	RegisterParent( Mesh_RegularAlgorithms_Type, Mesh_Algorithms_Type );
	RegisterParent( MeshTopology_Type, Stg_Component_Type );
	RegisterParent( Mesh_Type, Stg_Component_Type );
	RegisterParent( MeshGenerator_Type, Stg_Component_Type );
	RegisterParent( CartesianGenerator_Type, MeshGenerator_Type );
	RegisterParent( MeshAdaptor_Type, MeshGenerator_Type );
	RegisterParent( SurfaceAdaptor_Type, MeshAdaptor_Type );
	RegisterParent( CompressionAdaptor_Type, MeshAdaptor_Type );
	RegisterParent( LinearSpaceAdaptor_Type, MeshAdaptor_Type );
	RegisterParent( MeshVariable_Type, Variable_Type );
	RegisterParent( Remesher_Type, Stg_Component_Type );

	return True;
}


