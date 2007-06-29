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
*/
/** \file
**  Role:
**
** Assumptions:
**
** Invariants:
**
** Comments:
**
** $Id: MeshVariable.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Discretisaton_Mesh_MeshVariable_h__
#define __Discretisaton_Mesh_MeshVariable_h__

	/** Textual name of this class */
	extern const Type MeshVariable_Type;

	/** Virtual function types */

	/** Class contents */
	#define __MeshVariable				\
		/* General info */			\
		__Variable				\
							\
		/* Virtual info */			\
							\
		/* MeshVariable info */			\
		Mesh*			mesh;		\
		MeshTopology_Dim	topoDim;	\
		unsigned		meshArraySize;

	struct MeshVariable { __MeshVariable };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define MESHVARIABLE_DEFARGS	\
		VARIABLE_DEFARGS

	#define MESHVARIABLE_PASSARGS	\
		VARIABLE_PASSARGS

	MeshVariable* MeshVariable_New( Name name );
	MeshVariable* _MeshVariable_New( MESHVARIABLE_DEFARGS );
	void _MeshVariable_Init( MeshVariable* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _MeshVariable_Delete( void* meshVariable );
	void _MeshVariable_Print( void* meshVariable, Stream* stream );
	void _MeshVariable_Construct( void* meshVariable, Stg_ComponentFactory* cf, void* data );
	void _MeshVariable_Build( void* meshVariable, void* data );
	void _MeshVariable_Initialise( void* meshVariable, void* data );
	void _MeshVariable_Execute( void* meshVariable, void* data );
	void _MeshVariable_Destroy( void* meshVariable, void* data );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void MeshVariable_SetMesh( void* meshVariable, void* mesh );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	void MeshVariable_Destruct( MeshVariable* self );

	Index _MeshVariable_GetMeshArraySize( void* meshVariable );

#endif /* __Discretisaton_Mesh_MeshVariable_h__ */
