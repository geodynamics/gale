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
**	A specific swarm used for integration which is "aware" of the materials used.
**	A mapper relationship links it with the physical material swarm(s) of which can actually be
**	advected.
**
** Assumptions:
**
** Comments:
**
** $Id: IntegrationPointsSwarm.h 435 2007-03-04 10:42:10Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __PICellerator_IntegrationPointsSwarm_IntegrationPointsSwarm_h__
#define __PICellerator_IntegrationPointsSwarm_IntegrationPointsSwarm_h__

	/* Textual name of this class */
	extern const Type IntegrationPointsSwarm_Type;

	/* IntegrationPointsSwarm information */
	#define __IntegrationPointsSwarm \
		__Swarm \
		\
		FiniteElement_Mesh*                   mesh;                 \
		TimeIntegrator*                       timeIntegrator;       \
		WeightsCalculator*                    weights;              \
		IntegrationPointMapper*               mapper;               \
		SwarmVariable*                        localCoordVariable;    /** Set only if a local coord system swarm. */ \
		SwarmVariable*                        weightVariable;       \
		Materials_Register*                   materials_Register;   \
		Bool                                  recalculateWeights;

	struct IntegrationPointsSwarm { __IntegrationPointsSwarm };
	
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	/** Classic C++-style constructor */
	IntegrationPointsSwarm* IntegrationPointsSwarm_New(
		Name                                  name,
		void*                                 cellLayout,
		void*                                 particleLayout,
		Dimension_Index                       dim,
		SizeT                                 particleSize,
		Particle_InCellIndex                  cellParticleTblDelta,
		double                                extraParticlesFactor,
		FiniteElement_Mesh*                   mesh,
		TimeIntegrator*                       timeIntegrator,
		WeightsCalculator*                    weights,
		IntegrationPointMapper*               mapper,
		Bool                                  recalculateWeights,
		ExtensionManager_Register*            extensionMgr_Register,
		Variable_Register*                    swarmVariable_Register,
		Materials_Register*                   materials_Register,
		MPI_Comm                              comm);
	
	void* _IntegrationPointsSwarm_DefaultNew( Name name ) ;

	/** Private New */
	IntegrationPointsSwarm* _IntegrationPointsSwarm_New(
		SizeT                                           _sizeOfSelf, 
		Type                                            type,
		Stg_Class_DeleteFunction*                       _delete,
		Stg_Class_PrintFunction*                        _print,
		Stg_Class_CopyFunction*                         _copy, 
		Stg_Component_DefaultConstructorFunction*       _defaultConstructor,
		Stg_Component_ConstructFunction*                _construct,
		Stg_Component_BuildFunction*                    _build,
		Stg_Component_InitialiseFunction*               _initialise,
		Stg_Component_ExecuteFunction*                  _execute,
		Stg_Component_DestroyFunction*                  _destroy,
		Name                                            name,
		Bool                                            initFlag,
		CellLayout*                                     cellLayout,
		ParticleLayout*                                 particleLayout,
		Dimension_Index                                 dim,
		SizeT                                           particleSize, 
		Particle_InCellIndex                            cellParticleTblDelta, 
		double                                          extraParticlesFactor,
		FiniteElement_Mesh*                             mesh, 
		TimeIntegrator*                                 timeIntegrator,
		WeightsCalculator*                              weights,
		IntegrationPointMapper*                         mapper,
		Bool                                            recalculateWeights,
		ExtensionManager_Register*                      extensionMgr_Register,
		Variable_Register*                              swarmVariable_Register,
		Materials_Register*                             materials_Register,
		MPI_Comm                                        comm
		);

	void _IntegrationPointsSwarm_Construct( void* shape, Stg_ComponentFactory* cf, void* data ) ;

	void _IntegrationPointsSwarm_Init(
		void*                                           swarm,
		FiniteElement_Mesh*                             mesh, 
		TimeIntegrator*                                 timeIntegrator,
		WeightsCalculator*                              weights,
		IntegrationPointMapper*                         mapper,
		Materials_Register*                             materials_Register,
		Bool                                            recalculateWeights );

	/* Stg_Class_Delete IntegrationPointsSwarm implementation */
	void _IntegrationPointsSwarm_Delete( void* integrationPoints );
	void _IntegrationPointsSwarm_Print( void* integrationPoints, Stream* stream );
	#define IntegrationPointsSwarm_Copy( self ) \
		(IntegrationPointsSwarm*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define IntegrationPointsSwarm_DeepCopy( self ) \
		(IntegrationPointsSwarm*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _IntegrationPointsSwarm_Copy( void* integrationPoints, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	void _IntegrationPointsSwarm_Build( void* integrationPoints, void* data ) ;
	void _IntegrationPointsSwarm_Initialise( void* integrationPoints, void* data ) ;
	void _IntegrationPointsSwarm_Execute( void* integrationPoints, void* data );
	void _IntegrationPointsSwarm_Destroy( void* integrationPoints, void* data ) ;

	void _IntegrationPointsSwarm_UpdateHook( void* timeIntegrator, void* swarm );
	
	/** Update the weights and local positions of all integration points, based on the 
		current mapped material particles */
	void IntegrationPointsSwarm_RemapIntegrationPointsAndRecalculateWeights( void* swarm );
	
	/** Get the material index associated with this integration point */
	Material_Index IntegrationPointsSwarm_GetMaterialIndexOn(
		IntegrationPointsSwarm* swarm,
		IntegrationPoint* point ); 
	
	/** Get the material associated with this integration point */
	Material* IntegrationPointsSwarm_GetMaterialOn(
		IntegrationPointsSwarm* swarm,
		IntegrationPoint* point );

	/** Get the extension associated with this integration point */
	#define IntegrationPointsSwarm_GetExtensionOn( swarm, point, extHandle ) \
		IntegrationPointMapper_GetExtensionOn( ((IntegrationPointsSwarm*)(swarm))->mapper, (point), (extHandle) )


	/** Get the material index associated with this integration point index */
	#define IntegrationPointsSwarm_GetMaterialIndexAt( swarm, point_I ) \
		IntegrationPointMapper_GetMaterialIndexAt( ((IntegrationPointsSwarm*)(swarm))->mapper, (point_I) )

	/** Get the material associated with this integration point index */
	#define IntegrationPointsSwarm_GetMaterialAt( swarm, point_I ) \
		Materials_Register_GetByIndex( \
				 ((IntegrationPointsSwarm*)(swarm))->materials_Register, \
				IntegrationPointsSwarm_GetMaterialIndexAt( (swarm), (point_I) ) )

	/** Get the extension associated with this integration point index */
	#define IntegrationPointsSwarm_GetExtensionAt( swarm, point_I, extHandle ) \
		IntegrationPointMapper_GetExtensionAt( ((IntegrationPointsSwarm*)(swarm))->mapper, (point_I), (extHandle) )

#endif
