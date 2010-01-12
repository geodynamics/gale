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
 **	Abstract class faciliting how particles of a cell are arranged. 
 **	Swarms should only reference this... the base class, not the implementations.
 **
 ** Assumptions:
 **
 ** Comments:
 **
 ** $Id: ParticleLayout.h 3629 2006-06-13 08:36:51Z PatrickSunter $
 **
 **~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgDomain_Swarm_ParticleLayout_h__
#define __StgDomain_Swarm_ParticleLayout_h__

/** @see ParticleLayout_SetInitialCounts */
typedef void (ParticleLayout_SetInitialCountsFunction)    ( void* particleLayout, void* swarm );
/** @see ParticleLayout_InitialiseParticles */ 
typedef void (ParticleLayout_InitialiseParticlesFunction) ( void* particleLayout, void* swarm );
	
extern const Type ParticleLayout_Type;
	
/* ParticleLayout information */
#define __ParticleLayout \
    __Stg_Component \
    AbstractContext*										context; \
    /* virtual information */ \
    ParticleLayout_SetInitialCountsFunction*		_setInitialCounts; \
    ParticleLayout_InitialiseParticlesFunction*	_initialiseParticles; \
    /* member data */ \
    Stream*													debug; /**< For sending log messages to */ \
    CoordSystem											coordSystem; /**< Set once by subclass init */ \
    Bool														weightsInitialisedAtStartup; /**< Set once by subclass init */

struct ParticleLayout { __ParticleLayout };
	
/* Creation implementation / Virtual constructor */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define PARTICLELAYOUT_DEFARGS \
                STG_COMPONENT_DEFARGS, \
                ParticleLayout_SetInitialCountsFunction*               _setInitialCounts, \
                ParticleLayout_InitialiseParticlesFunction*         _initialiseParticles, \
                CoordSystem                                                  coordSystem, \
                Bool                                         weightsInitialisedAtStartup

	#define PARTICLELAYOUT_PASSARGS \
                STG_COMPONENT_PASSARGS, \
	        _setInitialCounts,           \
	        _initialiseParticles,        \
	        coordSystem,                 \
	        weightsInitialisedAtStartup

ParticleLayout* _ParticleLayout_New(  PARTICLELAYOUT_DEFARGS  );

/* Initialise implementation */
void _ParticleLayout_Init( 
      void*                     particleLayout,
      AbstractContext*          context,
      CoordSystem               coordSystem,
      Bool                      weightsInitialisedAtStartup );
	
/* Stg_Class_Delete a ParticleLayout construct */
void _ParticleLayout_Delete( void* particleLayout );
	
/* Print a ParticleLayout construct */
void _ParticleLayout_Print( void* particleLayout, Stream* stream );

void _ParticleLayout_Destroy( void* particleLayout, void* data );
	
/* Copy */
#define ParticleLayout_Copy( self )                                     \
    (ParticleLayout*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
#define ParticleLayout_DeepCopy( self )                                 \
    (ParticleLayout*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
	
void* _ParticleLayout_Copy( void* particleLayout, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
void _ParticleLayout_AssignFromXML( void* particleLayout, Stg_ComponentFactory *cf, void* data );
	
/** Setup the particle tables in a Swarm Class, both local and shadow values.
 *  Must setup cellParticleCountTbl, cellParticleSizeTbl, cellParticleTbl */
void ParticleLayout_SetInitialCounts( void* particleLayout, void* swarm );
	
/** Initialise startup positions of all particles in given swarm. */
void ParticleLayout_InitialiseParticles( void* particleLayout, void* swarm );
	
#endif /* __StgDomain_Swarm_ParticleLayout_h__ */

