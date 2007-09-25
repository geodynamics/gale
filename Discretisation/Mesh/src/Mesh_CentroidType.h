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
** $Id: Mesh_CentroidType.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Discretisaton_Mesh_Mesh_CentroidType_h__
#define __Discretisaton_Mesh_Mesh_CentroidType_h__

	/** Textual name of this class */
	extern const Type Mesh_CentroidType_Type;

	/** Virtual function types */

	/** Class contents */
	#define __Mesh_CentroidType			\
		/* General info */			\
		__Mesh_ElementType			\
							\
		/* Virtual info */			\
							\
		/* Mesh_CentroidType info */		\
		Mesh*			elMesh;		\
		IArray*			incArray;

	struct Mesh_CentroidType { __Mesh_CentroidType };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define MESH_CENTROIDTYPE_DEFARGS		\
		MESH_ELEMENTTYPE_DEFARGS

	#define MESH_CENTROIDTYPE_PASSARGS		\
		MESH_ELEMENTTYPE_PASSARGS

	Mesh_CentroidType* Mesh_CentroidType_New();
	Mesh_CentroidType* _Mesh_CentroidType_New( MESH_CENTROIDTYPE_DEFARGS );
	void _Mesh_CentroidType_Init( Mesh_CentroidType* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _Mesh_CentroidType_Delete( void* centroidType );
	void _Mesh_CentroidType_Print( void* centroidType, Stream* stream );

	void Mesh_CentroidType_Update( void* centroidType );
	Bool Mesh_CentroidType_ElementHasPoint( void* centroidType, unsigned elInd, double* point, 
					   MeshTopology_Dim* dim, unsigned* ind );
	double Mesh_CentroidType_GetMinimumSeparation( void* centroidType, unsigned elInd, double* perDim );
	void Mesh_CentroidType_GetCentroid( void* centroidType, unsigned element, double* centroid );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void Mesh_CentroidType_SetElementMesh( void* centroidType, void* mesh );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

#endif /* __Discretisaton_Mesh_Mesh_CentroidType_h__ */
