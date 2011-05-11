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
**	Instantiates the ParticleLayout abstract class to a manually distributed particle layout.
**
** Assumptions:
**	Cell is a right-angled cuboid.
**
** Comments:
**
** $Id: ManualParticleLayout.h 3851 2006-10-12 08:57:22Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgDomain_Swarm_ManualParticleLayout_h__
#define __StgDomain_Swarm_ManualParticleLayout_h__
	

	/* Textual name of this class */
	extern const Type ManualParticleLayout_Type;

	extern const Index ManualParticleLayout_Invalid;

	/* ManualParticleLayout information */
	#define __ManualParticleLayout \
		__GlobalParticleLayout \
		\
		Dictionary*             dictionary; /* TODO: don't think this dictionary needs to be defined on the component, JG, 7Dec09 */

	struct ManualParticleLayout { __ManualParticleLayout };
	
   /* Create a new ManualParticleLayout and initialise */
   ManualParticleLayout* ManualParticleLayout_New( Name name,
      AbstractContext* context, 
      CoordSystem      coordSystem,
      Bool             weightsInitialisedAtStartup,
      unsigned int     totalInitialParticles, 
      double           averageInitialParticlesPerCell,
      Dictionary*      dictionary );
	
	/* Creation implementation / Virtual constructor */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define MANUALPARTICLELAYOUT_DEFARGS \
                GLOBALPARTICLELAYOUT_DEFARGS, \
                Dictionary*  dictionary

	#define MANUALPARTICLELAYOUT_PASSARGS \
                GLOBALPARTICLELAYOUT_PASSARGS, \
	        dictionary

   ManualParticleLayout* _ManualParticleLayout_New(  MANUALPARTICLELAYOUT_DEFARGS  );
	
	void _ManualParticleLayout_Init( void* manualParticleLayout, Dictionary* dictionary );
	
	/* 'Stg_Class' Stuff */
	void _ManualParticleLayout_Delete( void* manualParticleLayout );
	void _ManualParticleLayout_Print( void* manualParticleLayout, Stream* stream );
	#define ManualParticleLayout_Copy( self ) \
		(ManualParticleLayout*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define ManualParticleLayout_DeepCopy( self ) \
		(ManualParticleLayout*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _ManualParticleLayout_Copy( const void* manualParticleLayout, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	/* 'Stg_Component' Stuff */
	void* _ManualParticleLayout_DefaultNew( Name name ) ;
	void _ManualParticleLayout_AssignFromXML( void* manualParticleLayout, Stg_ComponentFactory *cf, void* data );
	void _ManualParticleLayout_Build( void* manualParticleLayout, void* data );
	void _ManualParticleLayout_Initialise( void* manualParticleLayout, void* data );
	void _ManualParticleLayout_Execute( void* manualParticleLayout, void* data );
	void _ManualParticleLayout_Destroy( void* manualParticleLayout, void* data );
	
	void _ManualParticleLayout_InitialiseParticle( 
			void* manualParticleLayout, 
			void* swarm, 
			Particle_Index newParticle_I,
			void* particle);

#endif /* __StgDomain_Swarm_ManualParticleLayout_h__ */

