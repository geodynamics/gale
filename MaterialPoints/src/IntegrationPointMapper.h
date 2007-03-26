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
**     Abstract way to map integration points to a material properties.
**     Under this scheme, IntegrationPointsSwarms are essentailly 'virtual' swarms that are generated/mapped 
**     (and then reused/updated) at everytime time from the information obtained from MaterialPointsSwarm(s). 
**     MaterialPointsSwarm can undergo advection which makes it very dynamic.
**
**     This provides an abstract interface to access integration and material point properties.
**
** Assumptions:
**
** Comments:
**
** $Id: IntegrationPointMapper.h 189 2005-10-20 00:39:29Z RobertTurnbull $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __PICellerator_MaterialPoints_IntegrationPointMapper_h__
#define __PICellerator_MaterialPoints_IntegrationPointMapper_h__

	/* Textual name of this class */
	extern const Type IntegrationPointMapper_Type;

	/** @see IntegrtaionPointMapper_Map */
	typedef void (IntegrationPointMapper_MapFunction) ( void* mapper );

	/** @see IntegrationPointMapper_GetMaterialPointsSwarmsFunction */
	typedef MaterialPointsSwarm** (IntegrationPointMapper_GetMaterialPointsSwarmsFunction) ( void* mapper, Index* count );

	/** @see IntegrationPointMapper_GetMaterialIndexOn */
	typedef Material_Index (IntegrationPointMapper_GetMaterialIndexOnFunction) ( void* mapper, void* point );

	/** @see IntegrationPointMapper_GetExtensionOn */
	typedef void* (IntegrationPointMapper_GetExtensionOnFunction) ( 
				void*                   mapper, 
				void*                   points, 
				ExtensionInfo_Index     extHandle );
	
	/* IntegrationPointMapper information */
	#define __IntegrationPointMapper \
		__Stg_Component \
		\
		/* Virtual functions */ \
		IntegrationPointMapper_MapFunction*                             _map; \
		IntegrationPointMapper_GetMaterialPointsSwarmsFunction*         _getMaterialPointsSwarms; \
		IntegrationPointMapper_GetMaterialIndexOnFunction*              _getMaterialIndexOn; \
		IntegrationPointMapper_GetExtensionOnFunction*                  _getExtensionOn; \
		\
		/* General info */ \
		IntegrationPointsSwarm*                                     integrationSwarm;

	struct IntegrationPointMapper { __IntegrationPointMapper };
	
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/
	IntegrationPointMapper* _IntegrationPointMapper_New(
		SizeT                                                           _sizeOfSelf, 
		Type                                                            type,
		Stg_Class_DeleteFunction*                                       _delete,
		Stg_Class_PrintFunction*                                        _print,
		Stg_Class_CopyFunction*                                         _copy, 
		Stg_Component_DefaultConstructorFunction*                       _defaultConstructor,
		Stg_Component_ConstructFunction*                                _construct,
		Stg_Component_BuildFunction*                                    _build,
		Stg_Component_InitialiseFunction*                               _initialise,
		Stg_Component_ExecuteFunction*                                  _execute,
		Stg_Component_DestroyFunction*                                  _destroy,
		IntegrationPointMapper_MapFunction*                             _map, 
		IntegrationPointMapper_GetMaterialPointsSwarmsFunction*         _getMaterialPointsSwarms,
		IntegrationPointMapper_GetMaterialIndexOnFunction*              _getMaterialIndexOn, 
		IntegrationPointMapper_GetExtensionOnFunction*                  _getExtensionOn, 
		Name                                                            name,
		Bool                                                            initFlag,
		IntegrationPointsSwarm*                                         integrationSwarm );

	void _IntegrationPointMapper_Init( void* mapper, IntegrationPointsSwarm* integrationSwarm );

	void _IntegrationPointMapper_Delete( void* mapper );
	void _IntegrationPointMapper_Print( void* mapper, Stream* stream );
	#define IntegrationPointMapper_Copy( self ) \
		(IntegrationPointMapper*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define IntegrationPointMapper_DeepCopy( self ) \
		(IntegrationPointMapper*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _IntegrationPointMapper_Copy( void* mapper, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
void _IntegrationPointMapper_Construct( void* mapper, Stg_ComponentFactory* cf, void* data ) ;
	void _IntegrationPointMapper_Build( void* mapper, void* data ) ;
	void _IntegrationPointMapper_Initialise( void* mapper, void* data ) ;
	void _IntegrationPointMapper_Execute( void* mapper, void* data );
	void _IntegrationPointMapper_Destroy( void* mapper, void* data ) ;
	
	/** Performs a mapping between the MaterialPointSwarms and IntegrationPointsSwarm */
	void IntegrationPointMapper_Map( void* mapper );

	/** Returns the MaterialPointsSwarm(s) involved in this mapping. Allocates memory so the result must be deleted.
	 *     @param count Output parameter which stores the number of swarms returned */
	MaterialPointsSwarm** IntegrationPointMapper_GetMaterialPointsSwarms( void* mapper, Index* count );
	
	/** Returns the material index associated with this integration point by mapping to physical material swarm(s) */
	#ifdef MACRO_AS_FUNC
		#define IntegrationPointMapper_GetMaterialIndexOn IntegrationPointMapper_GetMaterialIndexOnFunc
	#else
		#define IntegrationPointMapper_GetMaterialIndexOn IntegrationPointMapper_GetMaterialIndexOnMacro
	#endif
	#define IntegrationPointMapper_GetMaterialIndexOnMacro( mapper, point ) \
			( (IntegrationPointMapper*)(mapper) )->_getMaterialIndexOn( (mapper), (point) ) 
	Material_Index IntegrationPointMapper_GetMaterialIndexOnFunc( void* mapper, void* point );

	
	/** Returns an extension associated with this integration point by mapping to physical material swarm(s) */
	#ifdef MACRO_AS_FUNC
		#define IntegrationPointMapper_GetExtensionOn IntegrationPointMapper_GetExtensionOnFunc
	#else
		#define IntegrationPointMapper_GetExtensionOn IntegrationPointMapper_GetExtensionOnMacro
	#endif
	#define IntegrationPointMapper_GetExtensionOnMacro( mapper, point, extHandle ) \
			( (IntegrationPointMapper*)(mapper) )->_getExtensionOn( (mapper), (point), (extHandle) )
	void* IntegrationPointMapper_GetExtensionOnFunc( void* mapper, void* point, ExtensionInfo_Index extHandle );


	/** Returns the material index associated with this integration point index by mapping to physical material swarm(s) */
	#ifdef MACRO_AS_FUNC
		#define IntegrationPointMapper_GetMaterialIndexAt IntegrationPointMapper_GetMaterialIndexAtFunc
	#else
		#define IntegrationPointMapper_GetMaterialIndexAt IntegrationPointMapper_GetMaterialIndexAtMacro
	#endif
	#define IntegrationPointMapper_GetMaterialIndexAtMacro( mapper, point_I ) \
			IntegrationPointMapper_GetMaterialIndexOn( \
					(mapper), \
					Swarm_ParticleAt( ((IntegrationPointMapper*)mapper)->integrationSwarm, point_I ) )
	Material_Index IntegrationPointMapper_GetMaterialIndexAtFunc( void* mapper, Index point_I );
	
	/** Returns an extension associated with this integration point index by mapping to physical material swarm(s) */
	#ifdef MACRO_AS_FUNC
		#define IntegrationPointMapper_GetExtensionAt
	#else
		#define IntegrationPointMapper_GetExtensionAt
	#endif
	#define IntegrationPointMapper_GetExtensionAtMacro( mapper, point_I, extHandle ) \
			IntegrationPointMapper_GetExtensionOn( \
					(mapper), \
					Swarm_ParticleAt( ((IntegrationPointMapper*)mapper)->integrationSwarm, point_I ), \
					extHandle )
	void* IntegrationPointMapper_GetExtensionAtFunc( void* mapper, Index point_I, ExtensionInfo_Index extHandle );

#endif
