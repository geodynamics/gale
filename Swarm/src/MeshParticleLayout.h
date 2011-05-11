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
**	Instantiates the ParticleLayout abstract class to a randomly distributed particle layout within each cell.
**
** Assumptions:
**	Cell is a right-angled cuboid.
**
** Comments:
**
** $Id: MeshParticleLayout.h 3851 2006-10-12 08:57:22Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgDomain_Swarm_MeshParticleLayout_h__
#define __StgDomain_Swarm_MeshParticleLayout_h__
	

	/* Textual name of this class */
	extern const Type MeshParticleLayout_Type;
	
	/* MeshParticleLayout information */
	#define __MeshParticleLayout \
		__PerCellParticleLayout \
		\
		Mesh*			     mesh; \
		Particle_InCellIndex         cellParticleCount; /**< Number of particles in each cell to create initially */ \
		unsigned int                 seed;              /**< The random seed used to generate positions. */

	struct MeshParticleLayout { __MeshParticleLayout };
	
	/* Create a new MeshParticleLayout and initialise */
   MeshParticleLayout* MeshParticleLayout_New( 
      Name                 name, 
      AbstractContext* context,
      CoordSystem      coordSystem,
      Bool             weightsInitialisedAtStartup,
      Mesh*            mesh,
      Particle_InCellIndex cellParticleCount, 
      unsigned int         seed );
	
	/* Creation implementation / Virtual constructor */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define MESHPARTICLELAYOUT_DEFARGS \
                PERCELLPARTICLELAYOUT_DEFARGS, \
                Mesh*                              mesh, \
                Particle_InCellIndex  cellParticleCount, \
                unsigned int                       seed

	#define MESHPARTICLELAYOUT_PASSARGS \
                PERCELLPARTICLELAYOUT_PASSARGS, \
	        mesh,              \
	        cellParticleCount, \
	        seed             

   MeshParticleLayout* _MeshParticleLayout_New(  MESHPARTICLELAYOUT_DEFARGS  );

	void _MeshParticleLayout_Init( 
			void*                meshParticleLayout, 
         Mesh*                mesh,
			Particle_InCellIndex cellParticleCount, 
			unsigned int         seed);

	/* 'Class' Stuff */
	void _MeshParticleLayout_Delete( void* meshParticleLayout );
	void _MeshParticleLayout_Print( void* meshParticleLayout, Stream* stream );
	
	#define MeshParticleLayout_Copy( self ) \
		(MeshParticleLayout*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define MeshParticleLayout_DeepCopy( self ) \
		(MeshParticleLayout*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _MeshParticleLayout_Copy( const void* meshParticleLayout, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	/* 'Stg_Component' Stuff */
	void* _MeshParticleLayout_DefaultNew( Name name ) ;
	void _MeshParticleLayout_AssignFromXML( void* meshParticleLayout, Stg_ComponentFactory* cf, void* data );
	void _MeshParticleLayout_Build( void* meshParticleLayout, void* data );
	void _MeshParticleLayout_Initialise( void* meshParticleLayout, void* data );
	void _MeshParticleLayout_Execute( void* meshParticleLayout, void* data );
	void _MeshParticleLayout_Destroy( void* meshParticleLayout, void* data );
	
	Particle_InCellIndex _MeshParticleLayout_InitialCount( void* meshParticleLayout, void* celllayout,  Cell_Index cell_I );
	
	void _MeshParticleLayout_InitialiseParticlesOfCell( void* meshParticleLayout, void* _swarm, Cell_Index cell_I );


	
#endif /* __StgDomain_Swarm_MeshParticleLayout_h__ */

