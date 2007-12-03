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
** $Id: Mesh_Algorithms.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Domain_Mesh_Mesh_Algorithms_h__
#define __Domain_Mesh_Mesh_Algorithms_h__

	/** Textual name of this class */
	extern const Type Mesh_Algorithms_Type;

	/** Virtual function types */
	typedef void (Mesh_Algorithms_SetMeshFunc)( void* algorithms, void* mesh );
	typedef void (Mesh_Algorithms_UpdateFunc)( void* algorithms );
	typedef unsigned (Mesh_Algorithms_NearestVertexFunc)( void* algorithms, double* point );
	typedef Bool (Mesh_Algorithms_SearchFunc)( void* algorithms, double* point, 
						   MeshTopology_Dim* dim, unsigned* ind );
	typedef Bool (Mesh_Algorithms_SearchElementsFunc)( void* algorithms, double* point, 
							   unsigned* elInd );
	typedef double (Mesh_Algorithms_GetMinimumSeparationFunc)( void* algorithms, double* perDim );
	typedef void (Mesh_Algorithms_GetLocalCoordRangeFunc)( void* algorithms, double* min, double* max );
	typedef void (Mesh_Algorithms_GetDomainCoordRangeFunc)( void* algorithms, double* min, double* max );
	typedef void (Mesh_Algorithms_GetGlobalCoordRangeFunc)( void* algorithms, double* min, double* max );

	/** Class contents */
	#define __Mesh_Algorithms								\
		/* General info */								\
		__Stg_Component									\
												\
		/* Virtual info */								\
		Mesh_Algorithms_SetMeshFunc*			setMeshFunc;			\
		Mesh_Algorithms_UpdateFunc*			updateFunc;			\
		Mesh_Algorithms_NearestVertexFunc*		nearestVertexFunc;		\
		Mesh_Algorithms_SearchFunc*			searchFunc;			\
		Mesh_Algorithms_SearchElementsFunc*		searchElementsFunc;		\
		Mesh_Algorithms_GetMinimumSeparationFunc*	getMinimumSeparationFunc;	\
		Mesh_Algorithms_GetLocalCoordRangeFunc*		getLocalCoordRangeFunc;		\
		Mesh_Algorithms_GetDomainCoordRangeFunc*	getDomainCoordRangeFunc;	\
		Mesh_Algorithms_GetGlobalCoordRangeFunc*	getGlobalCoordRangeFunc;	\
												\
		/* Mesh_Algorithms info */							\
		Mesh_Algorithms_NearestVertexFunc*	nearestVertex;				\
		Mesh_Algorithms_SearchFunc*		search;					\
		Mesh*					mesh;					\
		SpatialTree*				tree;					\
		int					rank;					\
		IArray*					incArray;

	struct Mesh_Algorithms { __Mesh_Algorithms };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define MESH_ALGORITHMS_DEFARGS								\
		STG_COMPONENT_DEFARGS,								\
		Mesh_Algorithms_SetMeshFunc*			setMeshFunc,			\
		Mesh_Algorithms_UpdateFunc*			updateFunc,			\
		Mesh_Algorithms_NearestVertexFunc*		nearestVertexFunc, 		\
		Mesh_Algorithms_SearchFunc*			searchFunc, 			\
		Mesh_Algorithms_SearchElementsFunc*		searchElementsFunc, 		\
		Mesh_Algorithms_GetMinimumSeparationFunc*	getMinimumSeparationFunc, 	\
		Mesh_Algorithms_GetLocalCoordRangeFunc*		getLocalCoordRangeFunc,		\
		Mesh_Algorithms_GetDomainCoordRangeFunc*	getDomainCoordRangeFunc,	\
		Mesh_Algorithms_GetGlobalCoordRangeFunc*	getGlobalCoordRangeFunc

	#define MESH_ALGORITHMS_PASSARGS	\
		STG_COMPONENT_PASSARGS, 	\
		setMeshFunc, 			\
		updateFunc, 			\
		nearestVertexFunc,		\
		searchFunc, 			\
		searchElementsFunc,		\
		getMinimumSeparationFunc,	\
		getLocalCoordRangeFunc, 	\
		getDomainCoordRangeFunc, 	\
		getGlobalCoordRangeFunc

	Mesh_Algorithms* Mesh_Algorithms_New( Name name );
	Mesh_Algorithms* _Mesh_Algorithms_New( MESH_ALGORITHMS_DEFARGS );
	void _Mesh_Algorithms_Init( Mesh_Algorithms* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _Mesh_Algorithms_Delete( void* algorithms );
	void _Mesh_Algorithms_Print( void* algorithms, Stream* stream );
	void _Mesh_Algorithms_Construct( void* algorithms, Stg_ComponentFactory* cf, void* data );
	void _Mesh_Algorithms_Build( void* algorithms, void* data );
	void _Mesh_Algorithms_Initialise( void* algorithms, void* data );
	void _Mesh_Algorithms_Execute( void* algorithms, void* data );
	void _Mesh_Algorithms_Destroy( void* algorithms, void* data );

	void _Mesh_Algorithms_SetMesh( void* algorithms, void* mesh );
	void _Mesh_Algorithms_Update( void* algorithms );
	unsigned _Mesh_Algorithms_NearestVertex( void* algorithms, double* point );
	Bool _Mesh_Algorithms_Search( void* algorithms, double* point, 
				      MeshTopology_Dim* dim, unsigned* ind );
	Bool _Mesh_Algorithms_SearchElements( void* algorithms, double* point, 
					      unsigned* elInd );
	double _Mesh_Algorithms_GetMinimumSeparation( void* algorithms, double* perDim );
	void _Mesh_Algorithms_GetLocalCoordRange( void* algorithms, double* min, double* max );
	void _Mesh_Algorithms_GetDomainCoordRange( void* algorithms, double* min, double* max );
	void _Mesh_Algorithms_GetGlobalCoordRange( void* algorithms, double* min, double* max );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	#define Mesh_Algorithms_SetMesh( self, mesh )							\
		VirtualCall( self, setMeshFunc, self, mesh )

	#define Mesh_Algorithms_Update( self )								\
		VirtualCall( self, updateFunc, self )

	#define Mesh_Algorithms_NearestVertex( self, point )					\
		VirtualCall( self, nearestVertexFunc, self, point )

	#define Mesh_Algorithms_Search( algorithms, point, dim, ind )				\
		(assert( (algorithms) && ((Mesh_Algorithms*)algorithms)->searchFunc ),			\
		 ((Mesh_Algorithms*)algorithms)->searchFunc( algorithms, point, dim, ind ))

	#define Mesh_Algorithms_SearchElements( algorithms, point, elInd )			\
		(assert( (algorithms) && ((Mesh_Algorithms*)algorithms)->searchElementsFunc ),		\
		 ((Mesh_Algorithms*)algorithms)->searchElementsFunc( algorithms, point, elInd ))

	#define Mesh_Algorithms_GetMinimumSeparation( algorithms, perDim )			\
		(assert( (algorithms) && ((Mesh_Algorithms*)algorithms)->getMinimumSeparationFunc ),	\
		 ((Mesh_Algorithms*)algorithms)->getMinimumSeparationFunc( algorithms, perDim ))

	#define Mesh_Algorithms_GetLocalCoordRange( algorithms, min, max )			\
		(assert( (algorithms) && ((Mesh_Algorithms*)algorithms)->getLocalCoordRangeFunc ),	\
		 ((Mesh_Algorithms*)algorithms)->getLocalCoordRangeFunc( algorithms, min, max ))

	#define Mesh_Algorithms_GetDomainCoordRange( algorithms, min, max )			\
		(assert( (algorithms) && ((Mesh_Algorithms*)algorithms)->getDomainCoordRangeFunc ),	\
		 ((Mesh_Algorithms*)algorithms)->getDomainCoordRangeFunc( algorithms, min, max ))

	#define Mesh_Algorithms_GetGlobalCoordRange( algorithms, min, max )			\
		(assert( (algorithms) && ((Mesh_Algorithms*)algorithms)->getGlobalCoordRangeFunc ),	\
		 ((Mesh_Algorithms*)algorithms)->getGlobalCoordRangeFunc( algorithms, min, max ))

	unsigned Mesh_Algorithms_NearestVertexWithNeighbours( void* algorithms, double* point );
	unsigned Mesh_Algorithms_NearestVertexGeneral( void* algorithms, double* point );
	Bool Mesh_Algorithms_SearchWithFullIncidence( void* algorithms, double* point, 
						      MeshTopology_Dim* dim, unsigned* ind );
	Bool Mesh_Algorithms_SearchWithMinIncidence( void* algorithms, double* point, 
						     MeshTopology_Dim* dim, unsigned* ind );
	Bool Mesh_Algorithms_SearchGeneral( void* algorithms, double* point, 
					    MeshTopology_Dim* dim, unsigned* ind );

Bool Mesh_Algorithms_SearchWithTree( void* self, double* pnt, unsigned* dim, unsigned* el );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

#endif /* __Domain_Mesh_Mesh_Algorithms_h__ */
