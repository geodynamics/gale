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
** $Id: Mesh_ElementType.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgDomain_Mesh_ElementType_h__
#define __StgDomain_Mesh_ElementType_h__

	/** Textual name of this class */
	extern const Type Mesh_ElementType_Type;

	/** Virtual function types */
	typedef void (Mesh_ElementType_UpdateFunc)( void* elementType );
	typedef Bool (Mesh_ElementType_ElementHasPointFunc)( void* elementType, unsigned element, double* point, 
							     MeshTopology_Dim* dim, unsigned* ind );
	typedef double (Mesh_ElementType_GetMinimumSeparationFunc)( void* elementType, unsigned element, double* perDim );
	typedef void (Mesh_ElementType_GetCentroidFunc)( void* elementType, unsigned element, double* centroid );

	/** Class contents */
	#define __Mesh_ElementType								\
		/* General info */								\
		__Stg_Class									\
												\
		/* Virtual info */								\
		Mesh_ElementType_UpdateFunc*			updateFunc;			\
		Mesh_ElementType_ElementHasPointFunc*		elementHasPointFunc;		\
		Mesh_ElementType_GetMinimumSeparationFunc*	getMinimumSeparationFunc;	\
		Mesh_ElementType_GetCentroidFunc*		getCentroidFunc;		\
												\
		/* Mesh_ElementType info */							\
		Mesh*			mesh;

	struct Mesh_ElementType { __Mesh_ElementType };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/



	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define MESH_ELEMENTTYPE_DEFARGS \
                STG_CLASS_DEFARGS, \
                Mesh_ElementType_UpdateFunc*                              updateFunc, \
                Mesh_ElementType_ElementHasPointFunc*            elementHasPointFunc, \
                Mesh_ElementType_GetMinimumSeparationFunc*  getMinimumSeparationFunc, \
                Mesh_ElementType_GetCentroidFunc*                    getCentroidFunc

	#define MESH_ELEMENTTYPE_PASSARGS \
                STG_CLASS_PASSARGS, \
	        updateFunc,               \
	        elementHasPointFunc,      \
	        getMinimumSeparationFunc, \
	        getCentroidFunc         

	Mesh_ElementType* _Mesh_ElementType_New(  MESH_ELEMENTTYPE_DEFARGS  );
	void _Mesh_ElementType_Init( Mesh_ElementType* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _Mesh_ElementType_Delete( void* elementType );
	void _Mesh_ElementType_Print( void* elementType, Stream* stream );
	void _Mesh_ElementType_GetCentroid( void* elementType, unsigned element, double* centroid );

	#define Mesh_ElementType_Update( self )							\
		VirtualCall( self, updateFunc, self )

	#define Mesh_ElementType_ElementHasPoint( self, element, point, dim, ind )		\
		VirtualCall( self, elementHasPointFunc, self, element, point, dim, ind )

	#define Mesh_ElementType_GetMinimumSeparation( self, element, perDim )			\
		VirtualCall( self, getMinimumSeparationFunc, self, element, perDim )

	#define Mesh_ElementType_GetCentroid( self, element, centroid )				\
		VirtualCall( self, getCentroidFunc, self, element, centroid )

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void Mesh_ElementType_SetMesh( void* elementType, void* mesh );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

#endif /* __StgDomain_Mesh_ElementType_h__ */

