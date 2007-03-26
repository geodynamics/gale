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
** $Id: DiscreteVoronoi.h 374 2006-10-12 08:59:41Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __PICellerator_Voronoi_DiscreteVoronoi_h__
#define __PICellerator_Voronoi_DiscreteVoronoi_h__

	typedef void (DiscreteVoronoi_CalculateForCellFunction)
		( void* discreteVoronoi, void* _swarm, Cell_LocalIndex lCell_I );
	typedef Particle_InCellIndex (DiscreteVoronoi_GetParticleIndexFunction)    
		( void* discreteVoronoi, Voronoi_CellIndex vCell_I );
	typedef double (DiscreteVoronoi_GetVolumeFunction)     
		( void* discreteVoronoi, Voronoi_CellIndex vCell_I );
	typedef void (DiscreteVoronoi_GetCentroidFunction)     
		( void* discreteVoronoi, Voronoi_CellIndex vCell_I, Coord centroid );

	/* Textual name of this class */
	extern const Type DiscreteVoronoi_Type;

	/* DiscreteVoronoi information */
	#define __DiscreteVoronoi \
		/* General info */ \
		__Stg_Component \
		/* Virtual Info */\
		DiscreteVoronoi_CalculateForCellFunction*  _calculateForCell;        \
		DiscreteVoronoi_GetParticleIndexFunction*  _getParticleIndex;        \
		DiscreteVoronoi_GetVolumeFunction*         _getVolume;               \
		DiscreteVoronoi_GetCentroidFunction*       _getCentroid;             \
		/* Other Info */\
		Dimension_Index                            dim;                      \
		Voronoi_CellIndex                          claimedCellCount; 

	struct DiscreteVoronoi { __DiscreteVoronoi };
	
	struct DiscreteVoronoiParticleInfo {
		Coord                centroid;
		double               volume;
		Particle_InCellIndex particle_I;
		Voronoi_CellIndex    voronoiCellCount;
	};
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/
	DiscreteVoronoi* _DiscreteVoronoi_New(
		SizeT                                      _sizeOfSelf, 
		Type                                       type,
		Stg_Class_DeleteFunction*                  _delete,
		Stg_Class_PrintFunction*                   _print,
		Stg_Class_CopyFunction*                    _copy, 
		Stg_Component_DefaultConstructorFunction*  _defaultConstructor,
		Stg_Component_ConstructFunction*           _construct,
		Stg_Component_BuildFunction*               _build,
		Stg_Component_InitialiseFunction*          _initialise,
		Stg_Component_ExecuteFunction*             _execute,
		Stg_Component_DestroyFunction*             _destroy,		
		DiscreteVoronoi_CalculateForCellFunction*  _calculate,
		DiscreteVoronoi_GetParticleIndexFunction*  _getParticleIndex,
		DiscreteVoronoi_GetVolumeFunction*         _getVolume,
		DiscreteVoronoi_GetCentroidFunction*       _getCentroid,
		Name                                       name );
	
	void _DiscreteVoronoi_Init( void* discreteVoronoi , Dimension_Index dim) ;
	void DiscreteVoronoi_InitAll( void* discreteVoronoi, Dimension_Index dim );

	/* Stg_Class_Delete DiscreteVoronoi implementation */
	void _DiscreteVoronoi_Delete( void* discreteVoronoi );
	void _DiscreteVoronoi_Print( void* discreteVoronoi, Stream* stream );
	#define DiscreteVoronoi_Copy( self ) \
		(DiscreteVoronoi*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define DiscreteVoronoi_DeepCopy( self ) \
		(DiscreteVoronoi*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _DiscreteVoronoi_Copy( void* discreteVoronoi, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
void _DiscreteVoronoi_Construct( void* shape, Stg_ComponentFactory* cf, void* data ) ;
	void _DiscreteVoronoi_Build( void* discreteVoronoi, void* data ) ;
	void _DiscreteVoronoi_Initialise( void* discreteVoronoi, void* data ) ;
	void _DiscreteVoronoi_Execute( void* discreteVoronoi, void* data );
	void _DiscreteVoronoi_Destroy( void* discreteVoronoi, void* data ) ;
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Public member functions
	*/

	/* Wrappers for virtual functions */
	void DiscreteVoronoi_CalculateForCell( void* discreteVoronoi, void* _swarm, Cell_LocalIndex lCell_I );
	#define DiscreteVoronoi_GetParticleIndex( discreteVoronoi, vCell_I ) \
		( ( (DiscreteVoronoi*) discreteVoronoi )->_getParticleIndex( discreteVoronoi, vCell_I ) )
	#define DiscreteVoronoi_GetVolume( discreteVoronoi, vCell_I ) \
		( ( (DiscreteVoronoi*) discreteVoronoi )->_getVolume( discreteVoronoi, vCell_I ) )
	#define DiscreteVoronoi_GetCentroid( discreteVoronoi, vCell_I, centroid ) \
		( ( (DiscreteVoronoi*) discreteVoronoi )->_getCentroid( discreteVoronoi, vCell_I, centroid ) )

	
	DiscreteVoronoiParticleInfo* DiscreteVoronoi_CreateParticleInfo( void* discreteVoronoi, void* _swarm, Cell_LocalIndex lCell_I ) ;
	void DiscreteVoronoiParticleInfo_SortByVolume( DiscreteVoronoiParticleInfo* particleInfo, Particle_InCellIndex cellParticleCount ) ;
	
	#define DiscreteVoronoiParticleInfo_Delete( particleInfo ) \
		Memory_Free( (particleInfo) )

#endif 
