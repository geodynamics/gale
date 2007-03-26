/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
** Copyright (c) 2005-2006, Monash Cluster Computing, Building 28, Monash University Clayton Campus,
**	Victoria, 3800, Australia
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**
** Contributors:
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Patrick D. Sunter, Software Engineer, VPAC. (patrick@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	David Stegman, Postdoctoral Fellow, Monash University. (david.stegman@sci.monash.edu.au)
**	Wendy Sharples, PhD Student, Monash University (wendy.sharples@sci.monash.edu.au)
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
** Comments:
**
** $Id: AdaptiveVoronoi.h 374 2006-10-12 08:59:41Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __PICellerator_Voronoi_AdaptiveVoronoi_h__
#define __PICellerator_Voronoi_AdaptiveVoronoi_h__

	/* Textual name of this class */
	extern const Type AdaptiveVoronoi_Type;

	typedef struct {
		Particle_InCellIndex particle_I;
		Coord                centroid;
		Coord                vertex[8];
	} AdaptiveVoronoiClaimedCell;

	typedef struct {
		Iteration_Index    level;
		Vertex_Index       vertex_I[8];
		Voronoi_CellIndex  nextCell_I;
	} AdaptiveVoronoiUnclaimedCell;

	typedef struct {
		Particle_Index particle_I;
		Coord          coord;
	} AdaptiveVoronoiVertex;

	/* AdaptiveVoronoi information */
	#define __AdaptiveVoronoi \
		/* General info */ \
		__DiscreteVoronoi \
		/* Virtual Info */\
		\
		/* General Info */\
		Iteration_Index                  maxIterations;          \
		Vertex_Index                     cellVertexCount;        \
		Vertex_Index                     newVertexCountForSplit; \
		/* Claimed Cell Info */ \
		AdaptiveVoronoiClaimedCell*      cellList;               \
		Voronoi_CellIndex                claimedCellsAllocated;  \
		/* Stored Cell Info */\
		Voronoi_CellIndex                cellsCount;             \
		Voronoi_CellIndex                cellsAllocatedCount;    \
		AdaptiveVoronoiUnclaimedCell*    cells;                  \
		/* Stored Vertex Info */\
		Vertex_Index                     verticiesCount;         \
		Vertex_Index                     verticiesAllocatedCount;\
		AdaptiveVoronoiVertex*           verticies;              \
		/* Current Info */ \
		Swarm*                           swarm;                  \
		Cell_LocalIndex                  cell_I;                 \

	struct AdaptiveVoronoi { __AdaptiveVoronoi };
	
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/
	AdaptiveVoronoi* _AdaptiveVoronoi_New(
		SizeT                                              _sizeOfSelf, 
		Type                                               type,
		Stg_Class_DeleteFunction*                          _delete,
		Stg_Class_PrintFunction*                           _print,
		Stg_Class_CopyFunction*                            _copy, 
		Stg_Component_DefaultConstructorFunction*          _defaultConstructor,
		Stg_Component_ConstructFunction*                   _construct,
		Stg_Component_BuildFunction*                       _build,
		Stg_Component_InitialiseFunction*                  _initialise,
		Stg_Component_ExecuteFunction*                     _execute,
		Stg_Component_DestroyFunction*                     _destroy,		
		DiscreteVoronoi_CalculateForCellFunction*          _calculate,		
		DiscreteVoronoi_GetParticleIndexFunction*          _getParticleIndex,
		DiscreteVoronoi_GetVolumeFunction*                 _getVolume,
		DiscreteVoronoi_GetCentroidFunction*               _getCentroid,	
		Name                                               name );

	/* Stg_Class_Delete AdaptiveVoronoi implementation */
	void _AdaptiveVoronoi_Delete( void* adaptiveVoronoi );
	void _AdaptiveVoronoi_Print( void* adaptiveVoronoi, Stream* stream );
	#define AdaptiveVoronoi_Copy( self ) \
		(AdaptiveVoronoi*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define AdaptiveVoronoi_DeepCopy( self ) \
		(AdaptiveVoronoi*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _AdaptiveVoronoi_Copy( void* adaptiveVoronoi, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	void* _AdaptiveVoronoi_DefaultNew( Name name ) ;
	void _AdaptiveVoronoi_Construct( void* adaptiveVoronoi, Stg_ComponentFactory* cf, void* data ) ;
	void _AdaptiveVoronoi_Build( void* adaptiveVoronoi, void* data ) ;
	void _AdaptiveVoronoi_Initialise( void* adaptiveVoronoi, void* data ) ;
	void _AdaptiveVoronoi_Execute( void* adaptiveVoronoi, void* data );
	void _AdaptiveVoronoi_Destroy( void* adaptiveVoronoi, void* data ) ;
	
	void _AdaptiveVoronoi_CalculateForCell( void* adaptiveVoronoi, void* _swarm, Cell_LocalIndex lCell_I ) ;
	Particle_InCellIndex _AdaptiveVoronoi_GetParticleIndex( void* discreteVoronoi, Voronoi_CellIndex vCell_I );
	double _AdaptiveVoronoi_GetVolume(void* discreteVoronoi, Voronoi_CellIndex vCell_I );
	void _AdaptiveVoronoi_GetCentroid(void* discreteVoronoi, Voronoi_CellIndex vCell_I, Coord centroid );

	/*---------------------------------------------------------------------------------------------------------------------
	** Private functions
	*/
	void AdaptiveVoronoi_InitialiseFirstCell( AdaptiveVoronoi* self ) ;
	Bool AdaptiveVoronoi_NeedSplitCell( AdaptiveVoronoi* self, AdaptiveVoronoiUnclaimedCell* cell ) ;
	void AdaptiveVoronoi_SplitCell( AdaptiveVoronoi* self, Index parentCell_I ) ;
	void AdaptiveVoronoi_ClaimCell( 
		AdaptiveVoronoi*              self,
		AdaptiveVoronoiUnclaimedCell* cell, 
		Particle_InCellIndex          particle_I );
	Voronoi_CellIndex AdaptiveVoronoi_Loop( AdaptiveVoronoi* self, Index cell_I ) ;
	void AdaptiveVoronoi_PrintLinkedList( AdaptiveVoronoi* self, Stream* stream );
	AdaptiveVoronoiClaimedCell* AdaptiveVoronoi_NewClaimedCell( AdaptiveVoronoi* self ) ;
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/
	
	
#endif 
