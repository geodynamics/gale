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
** $Id: CellularAutomataVoronoi.h 374 2006-10-12 08:59:41Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __PICellerator_Voronoi_CellularAutomataVoronoi_h__
#define __PICellerator_Voronoi_CellularAutomataVoronoi_h__

	/* Textual name of this class */
	extern const Type CellularAutomataVoronoi_Type;

	extern const unsigned int CELLULAR_AUTOMATA_UNCLAIMED;

	typedef Particle_InCellIndex                           BattlePair[2];

	typedef struct {
		BattlePair*                                        battlePair;
		Index                                              battleCount;
		Index                                              battlesAllocated;
	} CellularAutomataVoronoiBattleHistory;
	
	typedef struct CellularAutomataVoronoiCell             CellularAutomataVoronoiCell;
	struct CellularAutomataVoronoiCell {
		Coord                                              centroid;
		CellularAutomataVoronoiBattleHistory               battleHistory;
		CellularAutomataVoronoiCell**                      neighbourCellList;
		Voronoi_CellIndex                                  neighbourCount;
		Voronoi_CellIndex                                  neighboursAllocated;
		Particle_InCellIndex                               claimedParticle;
		Particle_InCellIndex                               particleToCheck;
	};

	typedef struct {
		CellularAutomataVoronoiCell**                      array;
		Voronoi_CellIndex                                  cellCount;
		Voronoi_CellIndex                                  cellsAllocated;
	} CellularAutomataVoronoiCellList;

	/* CellularAutomataVoronoi information */
	#define __CellularAutomataVoronoi \
		/* General info */ \
		__DiscreteVoronoi \
		/* Virtual Info */\
		\
		/* Other Info */\
		IJK                                                resolution;             \
		Bool                                               diagonalNeighbours;     \
		CellularAutomataVoronoiCell*                       cellList;               \
		XYZ                                                dx;                     \
		double*                                            cellVolumes;            \
		CellularAutomataVoronoiCellList                    cellsToGrow;            \
		CellularAutomataVoronoiCellList                    cellsToCheck;           \
		/* Current Info */ \
		Swarm*                                             swarm;                  \
		Cell_LocalIndex                                    lCell_I;                \
		
	struct CellularAutomataVoronoi { __CellularAutomataVoronoi };
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/
	CellularAutomataVoronoi* _CellularAutomataVoronoi_New(
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

	/* Stg_Class implementations */
	void _CellularAutomataVoronoi_Delete( void* cellularAutomataVoronoi );
	void _CellularAutomataVoronoi_Print( void* cellularAutomataVoronoi, Stream* stream );
	#define CellularAutomataVoronoi_Copy( self ) \
		(CellularAutomataVoronoi*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define CellularAutomataVoronoi_DeepCopy( self ) \
		(CellularAutomataVoronoi*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _CellularAutomataVoronoi_Copy( void* cellularAutomataVoronoi, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	void* _CellularAutomataVoronoi_DefaultNew( Name name ) ;
void _CellularAutomataVoronoi_Construct( void* cellularAutomataVoronoi, Stg_ComponentFactory* cf, void* data ) ;
	void _CellularAutomataVoronoi_Build( void* cellularAutomataVoronoi, void* data ) ;
	void _CellularAutomataVoronoi_Initialise( void* cellularAutomataVoronoi, void* data ) ;
	void _CellularAutomataVoronoi_Execute( void* cellularAutomataVoronoi, void* data );
	void _CellularAutomataVoronoi_Destroy( void* cellularAutomataVoronoi, void* data ) ;
	
	/* Discrete Voronoi implementations */
	void _CellularAutomataVoronoi_CalculateForCell( void* cellularAutomataVoronoi, void* _swarm, Cell_LocalIndex lCell_I ) ;
	Particle_InCellIndex _CellularAutomataVoronoi_GetParticleIndex( void* discreteVoronoi, Voronoi_CellIndex vCell_I );
	double _CellularAutomataVoronoi_GetVolume(void* discreteVoronoi, Voronoi_CellIndex vCell_I );
	void _CellularAutomataVoronoi_GetCentroid(void* discreteVoronoi, Voronoi_CellIndex vCell_I, Coord centroid );

	/* Other functions */
	void CellularAutomataVoronoi_CheckCells( void* cellularAutomataVoronoi ) ;
	void CellularAutomataVoronoi_GrowCells( void* cellularAutomataVoronoi ) ;
	Particle_InCellIndex CellularAutomataVoronoi_Battle(
		void*                                              cellularAutomataVoronoi, 
		CellularAutomataVoronoiCell*                       cell, 
		Particle_InCellIndex                               champion_I, 
		Particle_InCellIndex                               contender_I );

	void CellularAutomataVoronoi_InitialiseCell( void* cellularAutomataVoronoi, Swarm* swarm, Cell_LocalIndex lCell_I ) ;
	void CellularAutomataVoronoi_AddNeighbourToCell( void* cellularAutomataVoronoi, CellularAutomataVoronoiCell* cell, Index iIndex, Index jIndex, Index kIndex ) ;
	void CellularAutomataVoronoi_Seed( void* cellularAutomataVoronoi ) ;

	void CellularAutomataVoronoi_AsciiArt( void* cellularAutomataVoronoi ) ;

	void CellularAutomataVoronoiCellList_Build( CellularAutomataVoronoiCellList* cellList ) ;
	void CellularAutomataVoronoiCellList_AddCell( CellularAutomataVoronoiCellList* cellList, CellularAutomataVoronoiCell* cell ) ;

	void CellularAutomataVoronoi_CalcSubCells( CellularAutomataVoronoi* self, Swarm* swarm, unsigned cellInd );
	void CellularAutomataVoronoi_CalcSubCells2D( CellularAutomataVoronoi* self, Swarm* swarm, unsigned cellInd );
	void CellularAutomataVoronoi_CalcSubCells3D( CellularAutomataVoronoi* self, Swarm* swarm, unsigned cellInd );
	double CellularAutomataVoronoi_QuadArea( CellularAutomataVoronoi* self, double** gCrds );
	void CellularAutomataVoronoi_QuadCentroid( CellularAutomataVoronoi* self, double** gCrds, Coord centroid );
	void CellularAutomataVoronoi_HexCentroid( CellularAutomataVoronoi* self, double** gCrds, Coord centroid );
	
#endif 
