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
**	Handles the communication of particles between processors when their coordinates are updated.
**
** Assumptions:
**
** Comments:
**	If we ever decide we need more than one strategy for doing this, this will probably become an
**	interface and we'll have separate instantiation classes. Can't see that need in the short
**	term though so we'll just use the one class for now.
**
** $Id: ParticleMovementHandler.h 4001 2007-02-09 01:26:29Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	
#ifndef __Domain_Swarm_ParticleMovementHandler_h__
#define __Domain_Swarm_ParticleMovementHandler_h__

	/** Textual name of this class */
	extern const Type ParticleMovementHandler_Type;

	#define __ParticleMovementHandler \
		__ParticleCommHandler \
		/* Virtual info */ \
		/* Member info */ \
		Index                           globalParticlesArrivingMyDomainCount; \
		Index                           globalParticlesOutsideDomainTotal; \
		Bool                            useGlobalFallbackCommStrategy; \
		Bool                            defensive;


	struct ParticleMovementHandler { __ParticleMovementHandler };	

	/* --- virtual functions --- */

	/** Constructor interface */
	void* ParticleMovementHandler_DefaultNew( Name name );
	
	ParticleMovementHandler* ParticleMovementHandler_New(
			Name name,
			Bool useGlobalFallbackCommStrategy
			);
	
	/** Private Constructor */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define PARTICLEMOVEMENTHANDLER_DEFARGS \
                PARTICLECOMMHANDLER_DEFARGS, \
                Bool  useGlobalFallbackCommStrategy

	#define PARTICLEMOVEMENTHANDLER_PASSARGS \
                PARTICLECOMMHANDLER_PASSARGS, \
	        useGlobalFallbackCommStrategy

	ParticleMovementHandler* _ParticleMovementHandler_New(  PARTICLEMOVEMENTHANDLER_DEFARGS  );
	
	/** Variable initialiser */
	void _ParticleMovementHandler_Init(
		ParticleMovementHandler*     self,
		Bool                     useGlobalFallbackCommStrategy
		);

	/** Stg_Class_Print() implementation */
	void _ParticleMovementHandler_Print( void* pCommsHandler, Stream* stream );
	
	void _ParticleMovementHandler_AssignFromXML( void* pCommsHandler, Stg_ComponentFactory* cf, void* data );
	
	void _ParticleMovementHandler_Build( void* pCommsHandler, void *data );
	
	void _ParticleMovementHandler_Initialise( void* pCommsHandler, void *data );
	
	void _ParticleMovementHandler_Execute( void* pCommsHandler, void *data );

	void _ParticleMovementHandler_Destroy( void* pCommsHandler, void *data );
	
	/** Copy */
	#define ParticleMovementHandler_Copy( self ) \
		(ParticleMovementHandler*)ParticleMovementHandler_CopyFunc( self, NULL, False, NULL, NULL )
	#define ParticleMovementHandler_DeepCopy( self ) \
		(ParticleMovementHandler*)ParticleMovementHandler_CopyFunc( self, NULL, True, NULL, NULL )
	
	void* _ParticleMovementHandler_CopyFunc( void* ParticleMovementHandler, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	/** Stg_Class_Delete() implementation */
	void _ParticleMovementHandler_Delete(void* pCommsHandler );
	
	/* --- Public virtual function interfaces --- */
	
	/** Handle particle movement between processors */
	void ParticleMovementHandler_HandleParticleMovementBetweenProcs( ParticleCommHandler* pCommsHandler );

	/* --- virtual function implementations --- */

	/* +++ Global fallback method related +++ */
	void ParticleMovementHandler_DoGlobalFallbackCommunication( ParticleMovementHandler* self );

	void ParticleMovementHandler_FindParticlesThatHaveMovedOutsideMyDomain( ParticleMovementHandler* self );

	void ParticleMovementHandler_GetCountOfParticlesOutsideDomainPerProcessor(
		ParticleMovementHandler*	self,
		Particle_Index**	globalParticlesOutsideDomainCountsPtr,
		Particle_Index*		maxGlobalParticlesOutsideDomainCountPtr,
		Particle_Index*		globalParticlesOutsideDomainTotalPtr );
		
	void ParticleMovementHandler_ShareAndUpdateParticlesThatHaveMovedOutsideDomains(
		ParticleMovementHandler* self,
		Particle_Index*      globalParticlesArrivingMyDomainCountPtr,
		Particle_Index*      globalParticlesOutsideDomainTotalPtr );

	void ParticleMovementHandler_EnsureParticleCountLeavingDomainsEqualsCountEnteringGlobally( ParticleMovementHandler* self );

	void ParticleMovementHandler_ZeroGlobalCommStrategyCounters( ParticleMovementHandler* self );

	/* --- private functions --- */
	/* +++ Managment of particle array insertions/deletions +++ */
	Particle_Index ParticleMovementHandler_FindFreeSlotAndPrepareForInsertion( ParticleCommHandler* self );

	void ParticleMovementHandler_FillRemainingHolesInLocalParticlesArray( ParticleCommHandler*	self );

	Particle_Index* ParticleMovementHandler_MergeListsOfUnfilledParticleSlots( ParticleCommHandler* self );

	void ParticleMovementHandler_PrintParticleSlotsYetToFill( ParticleCommHandler* self );
	void ParticleMovementHandler_FinishReceiveAndUpdateShadowParticlesEnteringMyDomain( ParticleCommHandler* self );

#endif

