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
** $Id: Mesh_HexType.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgDomain_Mesh_HexType_h__
#define __StgDomain_Mesh_HexType_h__

	/** Textual name of this class */
	extern const Type Mesh_HexType_Type;

	/** Virtual function types */

	/** Class contents */
	#define __Mesh_HexType							\
		/* General info */						\
		__Mesh_ElementType						\
										\
		/* Virtual info */						\
										\
		/* Mesh_HexType info */						\
		unsigned				mapSize;		\
		unsigned*				vertMap;		\
		unsigned*				inc;			\
		Mesh_ElementType_ElementHasPointFunc*	elementHasPoint;	\
		unsigned**				triInds;		\
		unsigned**				tetInds;		\
                unsigned                                num_simplexes[2];       \
		IArray*					incArray;

	struct Mesh_HexType { __Mesh_HexType };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/



	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define MESH_HEXTYPE_DEFARGS \
                MESH_ELEMENTTYPE_DEFARGS

	#define MESH_HEXTYPE_PASSARGS \
                MESH_ELEMENTTYPE_PASSARGS

	Mesh_HexType* Mesh_HexType_New();
	Mesh_HexType* _Mesh_HexType_New(  MESH_HEXTYPE_DEFARGS  );
	void _Mesh_HexType_Init( Mesh_HexType* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _Mesh_HexType_Delete( void* hexType );
	void _Mesh_HexType_Print( void* hexType, Stream* stream );

	void Mesh_HexType_Update( void* hexType );
	Bool Mesh_HexType_ElementHasPoint( void* hexType, unsigned elInd, double* point, 
					   MeshTopology_Dim* dim, unsigned* ind );
	double Mesh_HexType_GetMinimumSeparation( void* hexType, unsigned elInd, double* perDim );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void Mesh_HexType_SetVertexMap( void* hexType, unsigned* map );
        void Mesh_HexType_SetQ2Inds( void* hexType);
	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	Bool Mesh_HexType_ElementHasPoint3DGeneral( Mesh_HexType* self, unsigned elInd, double* point, 
						    MeshTopology_Dim* dim, unsigned* ind );
	Bool Mesh_HexType_ElementHasPoint3DWithIncidence( Mesh_HexType* self, unsigned elInd, double* point, 
							  MeshTopology_Dim* dim, unsigned* ind );
	Bool Mesh_HexType_ElementHasPoint2DGeneral( Mesh_HexType* self, unsigned elInd, double* point, 
						    MeshTopology_Dim* dim, unsigned* ind );
	Bool Mesh_HexType_ElementHasPoint2DWithIncidence( Mesh_HexType* self, unsigned elInd, double* point, 
							  MeshTopology_Dim* dim, unsigned* ind );
	Bool Mesh_HexType_ElementHasPoint1DGeneral( Mesh_HexType* self, unsigned elInd, double* point, 
						    MeshTopology_Dim* dim, unsigned* ind );
	Bool Mesh_HexType_ElementHasPoint1DWithIncidence( Mesh_HexType* self, unsigned elInd, double* point, 
							  MeshTopology_Dim* dim, unsigned* ind );
	void Mesh_HexType_TetBarycenter( double** verts, unsigned* inc, unsigned* inds, double* point, double* bc );
	void Mesh_HexType_TriBarycenter( double** verts, unsigned* inc, unsigned* inds, double* point, double* bc );

#endif /* __StgDomain_Mesh_HexType_h__ */

