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
** $Id: CartesianGenerator.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgDomain_Mesh_CartesianGenerator_h__
#define __StgDomain_Mesh_CartesianGenerator_h__

	/** Textual name of this class */
	extern const Type CartesianGenerator_Type;

	/** Virtual function types */
	typedef void (CartesianGenerator_SetTopologyParamsFunc)( void* meshGenerator, unsigned* sizes, 
								 unsigned maxDecompDims, unsigned* minDecomp, unsigned* maxDecomp );
	typedef void (CartesianGenerator_GenElementsFunc)( void* meshGenerator, IGraph* topo, Grid*** grids );
	typedef void (CartesianGenerator_GenFacesFunc)( void* meshGenerator, IGraph* topo, Grid*** grids );
	typedef void (CartesianGenerator_GenEdgesFunc)( void* meshGenerator, IGraph* topo, Grid*** grids );
	typedef void (CartesianGenerator_GenVerticesFunc)( void* meshGenerator, IGraph* topo, Grid*** grids );
	typedef void (CartesianGenerator_GenElementVertexIncFunc)( void* meshGenerator, IGraph* topo, Grid*** grids );
	typedef void (CartesianGenerator_GenVolumeEdgeIncFunc)( void* meshGenerator, IGraph* topo, Grid*** grids );
	typedef void (CartesianGenerator_GenVolumeFaceIncFunc)( void* meshGenerator, IGraph* topo, Grid*** grids );
	typedef void (CartesianGenerator_GenFaceVertexIncFunc)( void* meshGenerator, IGraph* topo, Grid*** grids );
	typedef void (CartesianGenerator_GenFaceEdgeIncFunc)( void* meshGenerator, IGraph* topo, Grid*** grids );
	typedef void (CartesianGenerator_GenEdgeVertexIncFunc)( void* meshGenerator, IGraph* topo, Grid*** grids );
	typedef void (CartesianGenerator_GenElementTypesFunc)( void* meshGenerator, Mesh* mesh );

	/** CartesianGenerator class contents */
	#define __CartesianGenerator								\
		/* General info */								\
		__MeshGenerator									\
												\
		/* Virtual info */								\
		CartesianGenerator_SetTopologyParamsFunc*	setTopologyParamsFunc;		\
		CartesianGenerator_GenElementsFunc*		genElementsFunc;		\
		CartesianGenerator_GenFacesFunc*		genFacesFunc;			\
		CartesianGenerator_GenEdgesFunc*		genEdgesFunc;			\
		CartesianGenerator_GenVerticesFunc*		genVerticesFunc;		\
		CartesianGenerator_GenElementVertexIncFunc*	genElementVertexIncFunc;	\
		CartesianGenerator_GenVolumeEdgeIncFunc*	genVolumeEdgeIncFunc;		\
		CartesianGenerator_GenVolumeFaceIncFunc*	genVolumeFaceIncFunc;		\
		CartesianGenerator_GenFaceVertexIncFunc*	genFaceVertexIncFunc;		\
		CartesianGenerator_GenFaceEdgeIncFunc*		genFaceEdgeIncFunc;		\
		CartesianGenerator_GenEdgeVertexIncFunc*	genEdgeVertexIncFunc;		\
		CartesianGenerator_GenElementTypesFunc*		genElementTypesFunc;		\
												\
		/* CartesianGenerator info */							\
		Comm*		comm;								\
		Bool		regular;							\
		Bool		periodic[3];							\
		/* read cartesian mesh data from checkpoint file? */                            \
		Bool		readFromFile;							\
		unsigned	maxDecompDims;							\
		unsigned*	minDecomp;							\
		unsigned*	maxDecomp;							\
		unsigned	shadowDepth;							\
		double*		crdMin;								\
		double*		crdMax;								\
												\
		Grid*		vertGrid;							\
		Grid*		elGrid;								\
		Grid*		procGrid;							\
		unsigned*	origin;								\
		unsigned*	range;								\
		unsigned*	vertOrigin;							\
		unsigned*	vertRange;                              \
                int             contactDepth[3][2];                     \
                double          contactGeom[3];

	struct CartesianGenerator { __CartesianGenerator };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/



	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define CARTESIANGENERATOR_DEFARGS \
                MESHGENERATOR_DEFARGS, \
                CartesianGenerator_SetTopologyParamsFunc*      setTopologyParamsFunc, \
                CartesianGenerator_GenElementsFunc*                  genElementsFunc, \
                CartesianGenerator_GenFacesFunc*                        genFacesFunc, \
                CartesianGenerator_GenEdgesFunc*                        genEdgesFunc, \
                CartesianGenerator_GenVerticesFunc*                  genVerticesFunc, \
                CartesianGenerator_GenElementVertexIncFunc*  genElementVertexIncFunc, \
                CartesianGenerator_GenVolumeEdgeIncFunc*        genVolumeEdgeIncFunc, \
                CartesianGenerator_GenVolumeFaceIncFunc*        genVolumeFaceIncFunc, \
                CartesianGenerator_GenFaceVertexIncFunc*        genFaceVertexIncFunc, \
                CartesianGenerator_GenFaceEdgeIncFunc*            genFaceEdgeIncFunc, \
                CartesianGenerator_GenEdgeVertexIncFunc*        genEdgeVertexIncFunc, \
                CartesianGenerator_GenElementTypesFunc*          genElementTypesFunc

	#define CARTESIANGENERATOR_PASSARGS \
                MESHGENERATOR_PASSARGS, \
	        setTopologyParamsFunc,   \
	        genElementsFunc,         \
	        genFacesFunc,            \
	        genEdgesFunc,            \
	        genVerticesFunc,         \
	        genElementVertexIncFunc, \
	        genVolumeEdgeIncFunc,    \
	        genVolumeFaceIncFunc,    \
	        genFaceVertexIncFunc,    \
	        genFaceEdgeIncFunc,      \
	        genEdgeVertexIncFunc,    \
	        genElementTypesFunc    

	CartesianGenerator* CartesianGenerator_New( Name name, AbstractContext* context );
	CartesianGenerator* _CartesianGenerator_New(  CARTESIANGENERATOR_DEFARGS  );
	void _CartesianGenerator_Init( CartesianGenerator* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _CartesianGenerator_Delete( void* meshGenerator );
	void _CartesianGenerator_Print( void* meshGenerator, Stream* stream );
	void _CartesianGenerator_AssignFromXML( void* meshGenerator, Stg_ComponentFactory* cf, void* data );
	void _CartesianGenerator_Build( void* meshGenerator, void* data );
	void _CartesianGenerator_Initialise( void* meshGenerator, void* data );
	void _CartesianGenerator_Execute( void* meshGenerator, void* data );
	void _CartesianGenerator_Destroy( void* meshGenerator, void* data );

	void CartesianGenerator_SetDimSize( void* meshGenerator, unsigned nDims );
	void CartesianGenerator_Generate( void* meshGenerator, void* mesh, void* data );
	void _CartesianGenerator_SetTopologyParams( void* meshGenerator, unsigned* sizes, 
						    unsigned maxDecompDims, unsigned* minDecomp, unsigned* maxDecomp );
	void _CartesianGenerator_GenElements( void* meshGenerator, IGraph* topo, Grid*** grids );
	void _CartesianGenerator_GenFaces( void* meshGenerator, IGraph* topo, Grid*** grids );
	void _CartesianGenerator_GenEdges( void* meshGenerator, IGraph* topo, Grid*** grids );
	void _CartesianGenerator_GenVertices( void* meshGenerator, IGraph* topo, Grid*** grids );
	void _CartesianGenerator_GenElementVertexInc( void* meshGenerator, IGraph* topo, Grid*** grids );
	void _CartesianGenerator_GenVolumeEdgeInc( void* meshGenerator, IGraph* topo, Grid*** grids );
	void _CartesianGenerator_GenVolumeFaceInc( void* meshGenerator, IGraph* topo, Grid*** grids );
	void _CartesianGenerator_GenFaceVertexInc( void* meshGenerator, IGraph* topo, Grid*** grids );
	void _CartesianGenerator_GenFaceEdgeInc( void* meshGenerator, IGraph* topo, Grid*** grids );
	void _CartesianGenerator_GenEdgeVertexInc( void* meshGenerator, IGraph* topo, Grid*** grids );
	void _CartesianGenerator_GenElementTypes( void* meshGenerator, Mesh* mesh );

	#define CartesianGenerator_SetTopologyParams( self, sizes, maxDecompDims, minDecomp, maxDecomp )	\
		VirtualCall( self, setTopologyParamsFunc, self, sizes, maxDecompDims, minDecomp, maxDecomp )
	#define CartesianGenerator_GenElements( self, topo, grids )						\
		VirtualCall( self, genElementsFunc, self, topo, grids )
	#define CartesianGenerator_GenFaces( self, topo, grids )						\
		VirtualCall( self, genFacesFunc, self, topo, grids )
	#define CartesianGenerator_GenEdges( self, topo, grids )						\
		VirtualCall( self, genEdgesFunc, self, topo, grids )
	#define CartesianGenerator_GenVertices( self, topo, grids )						\
		VirtualCall( self, genVerticesFunc, self, topo, grids )
	#define CartesianGenerator_GenElementVertexInc( self, topo, grids )					\
		VirtualCall( self, genElementVertexIncFunc, self, topo, grids )
	#define CartesianGenerator_GenVolumeEdgeInc( self, topo, grids )					\
		VirtualCall( self, genVolumeEdgeIncFunc, self, topo, grids )
	#define CartesianGenerator_GenVolumeFaceInc( self, topo, grids )					\
		VirtualCall( self, genVolumeFaceIncFunc, self, topo, grids )
	#define CartesianGenerator_GenFaceVertexInc( self, topo, grids )					\
		VirtualCall( self, genFaceVertexIncFunc, self, topo, grids )
	#define CartesianGenerator_GenFaceEdgeInc( self, topo, grids )						\
		VirtualCall( self, genFaceEdgeIncFunc, self, topo, grids )
	#define CartesianGenerator_GenEdgeVertexInc( self, topo, grids )					\
		VirtualCall( self, genEdgeVertexIncFunc, self, topo, grids )
	#define CartesianGenerator_GenElementTypes( self, mesh )						\
		VirtualCall( self, genElementTypesFunc, self, mesh )

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void CartesianGenerator_SetGeometryParams( void* meshGenerator, double* min, double* max );
	void CartesianGenerator_SetShadowDepth( void* meshGenerator, unsigned depth );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	void CartesianGenerator_BuildDecomp( CartesianGenerator* self );
	void CartesianGenerator_RecurseDecomps( CartesianGenerator* self, 
						unsigned dim, unsigned max, 
						unsigned* nSubDomains, 
						unsigned* nPos, unsigned*** posNSubDomains );
	void CartesianGenerator_GenTopo( CartesianGenerator* self, IGraph* topo );
	void CartesianGenerator_GenEdges2D( CartesianGenerator* self, IGraph* topo, Grid*** grids );
	void CartesianGenerator_GenEdges3D( CartesianGenerator* self, IGraph* topo, Grid*** grids );
	void CartesianGenerator_GenBndVerts( CartesianGenerator* self, IGraph* topo, Grid*** grids );
	void CartesianGenerator_CompleteVertexNeighbours( CartesianGenerator* self, IGraph* topo, Grid*** grids );
	void CartesianGenerator_MapToDomain( CartesianGenerator* self, 
                                             const Sync* sync, 
					     unsigned nIncEls, unsigned* incEls );
	void CartesianGenerator_GenGeom( CartesianGenerator* self, Mesh* mesh, void* data );
	void CartesianGenerator_CalcGeom( CartesianGenerator* self, Mesh* mesh, Sync* sync, Grid* grid, unsigned* inds, double* steps );
	void CartesianGenerator_Destruct( CartesianGenerator* self );
	void CartesianGenerator_DestructTopology( CartesianGenerator* self );
	void CartesianGenerator_DestructGeometry( CartesianGenerator* self );
	void CartesianGenerator_ReadFromHDF5(  CartesianGenerator* self, Mesh* mesh, const char* filename );
	void CartesianGenerator_ReadFromASCII( CartesianGenerator* self, Mesh* mesh, const char* filename );

#endif /* __StgDomain_Mesh_CartesianGenerator_h__ */

