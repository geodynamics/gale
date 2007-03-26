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
** $Id: RemovalRoutine.h 374 2006-10-12 08:59:41Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __PICellerator_PopulationControl_RemovalRoutine_h__
#define __PICellerator_PopulationControl_RemovalRoutine_h__

	typedef void (RemovalRoutine_RemoveFromCellFunction)( void* removalRoutine, void* _swarm, Cell_LocalIndex lCell_I );

	/* Textual name of this class */
	extern const Type RemovalRoutine_Type;

	/* RemovalRoutine information */
	#define __RemovalRoutine \
		/* General info */ \
		__Stg_Component \
		/* Virtual Info */\
		RemovalRoutine_RemoveFromCellFunction*     _removeFromCell;          \
		/* Other Info */\
		Stream*                                    debug;                    \
		Dimension_Index                            dim;                      \
		Particle_InCellIndex                       idealParticleCount;       \
		Particle_InCellIndex                       maxParticlesPerCell;      \
		/* Removal Info */  \
		Particle_Index                             particlesToRemoveCount;   \
		Particle_Index                             particlesToRemoveAlloced; \
		Particle_Index                             particlesToRemoveDelta;   \
		ParticleToRemoveInfo*                      particlesToRemoveList;    


	typedef struct {
		Particle_Index       lParticle_I;
		Cell_Index           cell_I;
	} ParticleToRemoveInfo;

	struct RemovalRoutine { __RemovalRoutine };

	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/
	RemovalRoutine* _RemovalRoutine_New(
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
		RemovalRoutine_RemoveFromCellFunction*     _removeFromCell,
		Name                                       name );
	
	/* Stg_Class_Delete RemovalRoutine implementation */
	void _RemovalRoutine_Delete( void* removalRoutine );
	void _RemovalRoutine_Print( void* removalRoutine, Stream* stream );
	#define RemovalRoutine_Copy( self ) \
		(RemovalRoutine*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define RemovalRoutine_DeepCopy( self ) \
		(RemovalRoutine*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _RemovalRoutine_Copy( void* removalRoutine, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
void _RemovalRoutine_Construct( void* shape, Stg_ComponentFactory* cf, void* data ) ;
	void _RemovalRoutine_Build( void* removalRoutine, void* data ) ;
	void _RemovalRoutine_Initialise( void* removalRoutine, void* data ) ;
	void _RemovalRoutine_Execute( void* removalRoutine, void* data );
	void _RemovalRoutine_Destroy( void* removalRoutine, void* data ) ;
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Public member functions
	*/

	void RemovalRoutine_RemoveFromCell( void* removalRoutine, void* _swarm, Cell_LocalIndex lCell_I ) ;
	void RemovalRoutine_RemoveFromSwarm( void* removalRoutine, void* _swarm ) ;

	void RemovalRoutine_InitialiseParticleList( void* removalRoutine ) ;
	void RemovalRoutine_SetParticleToRemove( void* removalRoutine, Swarm* swarm, Cell_Index cell_I, Particle_InCellIndex cParticle_I ) ;

	void RemovalRoutine_SortParticleList( void* removalRoutine ) ;
	void RemovalRoutine_RemoveParticles( void* removalRoutine, Swarm* swarm ) ;

#endif 
