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
** $Id: FeMesh_Algorithms.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_Discretisation_FeMesh_Algorithms_h__
#define __StgFEM_Discretisation_FeMesh_Algorithms_h__

	/** Textual name of this class */
	extern const Type FeMesh_Algorithms_Type;

	/** Virtual function types */

	/** Class contents */
	#define __FeMesh_Algorithms								\
		/* General info */								\
		__Mesh_Algorithms								\
												\
		/* Virtual info */								\
												\
		/* FeMesh_Algorithms info */

	struct FeMesh_Algorithms { __FeMesh_Algorithms };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define FEMESH_ALGORITHMS_DEFARGS	\
		MESH_ALGORITHMS_DEFARGS

	#define FEMESH_ALGORITHMS_PASSARGS	\
		MESH_ALGORITHMS_PASSARGS

	FeMesh_Algorithms* FeMesh_Algorithms_New( Name name, AbstractContext* context );
	FeMesh_Algorithms* _FeMesh_Algorithms_New( FEMESH_ALGORITHMS_DEFARGS );
	void _FeMesh_Algorithms_Init( FeMesh_Algorithms* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _FeMesh_Algorithms_Delete( void* algorithms );
	void _FeMesh_Algorithms_Print( void* algorithms, Stream* stream );
	void _FeMesh_Algorithms_AssignFromXML( void* algorithms, Stg_ComponentFactory* cf, void* data );
	void _FeMesh_Algorithms_Build( void* algorithms, void* data );
	void _FeMesh_Algorithms_Initialise( void* algorithms, void* data );
	void _FeMesh_Algorithms_Execute( void* algorithms, void* data );
	void _FeMesh_Algorithms_Destroy( void* algorithms, void* data );

	Bool _FeMesh_Algorithms_Search( void* algorithms, double* point, 
					MeshTopology_Dim* dim, unsigned* ind );
	Bool _FeMesh_Algorithms_SearchElements( void* algorithms, double* point, 
						unsigned* elInd );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	Bool FeMesh_Algorithms_SearchWithTree( void* _self, double* pnt, unsigned* dim, unsigned* el );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

#endif /* __StgFEM_Discretisation_FeMesh_Algorithms_h__ */
