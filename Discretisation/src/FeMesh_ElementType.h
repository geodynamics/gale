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
** $Id: FeMesh_ElementType.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Domain_Mesh_FeMesh_ElementType_h__
#define __Domain_Mesh_FeMesh_ElementType_h__

	/** Textual name of this class */
	extern const Type FeMesh_ElementType_Type;

	/** Virtual function types */

	/** Class contents */
	#define __FeMesh_ElementType						\
		/* General info */						\
		__Mesh_HexType							\
										\
		/* Virtual info */						\
										\
		/* FeMesh_ElementType info */					\
		double	local[3];

	struct FeMesh_ElementType { __FeMesh_ElementType };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/



	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define FEMESH_ELEMENTTYPE_DEFARGS \
                MESH_HEXTYPE_DEFARGS

	#define FEMESH_ELEMENTTYPE_PASSARGS \
                MESH_HEXTYPE_PASSARGS

	FeMesh_ElementType* FeMesh_ElementType_New();
	FeMesh_ElementType* _FeMesh_ElementType_New(  FEMESH_ELEMENTTYPE_DEFARGS  );
	void _FeMesh_ElementType_Init( FeMesh_ElementType* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _FeMesh_ElementType_Delete( void* hexType );
	void _FeMesh_ElementType_Print( void* hexType, Stream* stream );

	void FeMesh_ElementType_Update( void* hexType );
	Bool FeMesh_ElementType_ElementHasPoint( void* hexType, unsigned elInd, double* point, 
					   MeshTopology_Dim* dim, unsigned* ind );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

#endif /* __Domain_Mesh_FeMesh_ElementType_h__ */

