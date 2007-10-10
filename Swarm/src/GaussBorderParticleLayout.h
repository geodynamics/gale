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
**	Gaussian particle layout with particles located on borders - intended for evaluation of boundary integral terms	
**
** Assumptions:
**	Cell is a right-angled cuboid.
**
** Comments:
**	
**
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Domain_Swarm_GaussBorderParticleLayout_h__
#define __Domain_Swarm_GaussBorderParticleLayout_h__
	

	/* Textual name of this class */
	extern const Type GaussBorderParticleLayout_Type;
	
	/* GaussBorderParticleLayout information */
	#define __GaussBorderParticleLayout \
                __GaussParticleLayout 	\
		Particle_InCellIndex*	particlesPerFace; // determined by particlesPerDim info, but calculated so many times it's worth making a member
	       
	struct GaussBorderParticleLayout { __GaussBorderParticleLayout };
	
	/* Create a new GaussBorderParticleLayout and initialise */
	GaussBorderParticleLayout* GaussBorderParticleLayout_New( Name name, Dimension_Index dim, Particle_InCellIndex* particlesPerDim ) ;
	
	/* Creation implementation / Virtual constructor */
	GaussBorderParticleLayout* _GaussBorderParticleLayout_New( 
		SizeT                                                       _sizeOfSelf,
		Type                                                        type,
		Stg_Class_DeleteFunction*                                   _delete,
		Stg_Class_PrintFunction*                                    _print,
		Stg_Class_CopyFunction*                                     _copy, 
		Stg_Component_DefaultConstructorFunction*                   _defaultConstructor,
		Stg_Component_ConstructFunction*                            _construct,
		Stg_Component_BuildFunction*                                _build,
		Stg_Component_InitialiseFunction*                           _initialise,
		Stg_Component_ExecuteFunction*                              _execute,
		Stg_Component_DestroyFunction*                              _destroy,
		ParticleLayout_SetInitialCountsFunction*                    _setInitialCounts,
		ParticleLayout_InitialiseParticlesFunction*                 _initialiseParticles,
		PerCellParticleLayout_InitialCountFunction*                 _initialCount,
		PerCellParticleLayout_InitialiseParticlesOfCellFunction*    _initialiseParticlesOfCell,
		Name                                                        name,
		Bool                                                        initFlag,
		Dimension_Index                                             dim,
		Particle_InCellIndex*                                       particlesPerDim );
		
	/* Initialise implementation */
	void _GaussBorderParticleLayout_Init( void* gaussBorderParticleLayout, Dimension_Index dim, Particle_InCellIndex* particlesPerDim );
	
	/* Stg_Class_Delete implementation */
	void _GaussBorderParticleLayout_Delete( void* gaussBorderParticleLayout );
	
	/* Print implementation */
	void _GaussBorderParticleLayout_Print( void* gaussBorderParticleLayout, Stream* stream );
	
	/* Copy */
	#define GaussBorderParticleLayout_Copy( self ) \
		(GaussBorderParticleLayout*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define GaussBorderParticleLayout_DeepCopy( self ) \
		(GaussBorderParticleLayout*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
	
	void* _GaussBorderParticleLayout_Copy( void* gaussBorderParticleLayout, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	void* _GaussBorderParticleLayout_DefaultNew( Name name );
	void  _GaussBorderParticleLayout_Construct( void* gaussBorderParticleLayout, Stg_ComponentFactory* cf, void* data );
	void  _GaussBorderParticleLayout_Build( void* gaussBorderParticleLayout, void* data );
	void  _GaussBorderParticleLayout_Initialise( void* gaussBorderParticleLayout, void* data );
	void  _GaussBorderParticleLayout_Execute( void* gaussBorderParticleLayout, void* data );
	void  _GaussBorderParticleLayout_Destroy( void* gaussBorderParticleLayout, void* data );
	
	Particle_InCellIndex _GaussBorderParticleLayout_InitialCount( void* gaussBorderParticleLayout, void* celllayout, Cell_Index cell_I );
	void _GaussBorderParticleLayout_InitialiseParticlesOfCell( void* gaussBorderParticleLayout, void* swarm, Cell_Index cell_I );

        Dimension_Index GaussBorderParticleLayout_GetFaceAxis( void* gaussBorderParticleLayout, Index face_I, Dimension_Index axis);
        Index GaussBorderParticleLayout_ParticleInCellIndexToFaceIndex( void* gaussBorderParticleLayout, Particle_InCellIndex cParticle_I );

        void _GaussBorderParticleLayout_InitialiseParticlesPerFace( GaussBorderParticleLayout* self );
	
#endif /* __Domain_Swarm_GaussBorderParticleLayout_h__ */
