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
**	Handles the shadowing of particles. So procs can use particles information from neighbouring procs.
**
** Assumptions:
**
** Comments:
**	
**
** $Id: ParticleShadowSync.h 4001 2007-02-09 01:26:29Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	
#ifndef __StgDomain_Swarm_ParticleShadowSync_h__
#define __StgDomain_Swarm_ParticleShadowSync_h__

	/** Textual name of this class */
	extern const Type ParticleShadowSync_Type;

	#define __ParticleShadowSync \
		__ParticleCommHandler 
		/* Virtual info */ 
		/* Member info */ 


	struct ParticleShadowSync { __ParticleShadowSync };	

	/* --- virtual functions --- */

	/** Constructor interface */
	void* ParticleShadowSync_DefaultNew( Name name );
	
	ParticleShadowSync* ParticleShadowSync_New(
			Name name,
			void* swarm
			);
	
	/** Private Constructor */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define PARTICLESHADOWSYNC_DEFARGS \
                PARTICLECOMMHANDLER_DEFARGS

	#define PARTICLESHADOWSYNC_PASSARGS \
                PARTICLECOMMHANDLER_PASSARGS

	ParticleShadowSync* _ParticleShadowSync_New(  PARTICLESHADOWSYNC_DEFARGS  );
	
	/** Variable initialiser */
	void _ParticleShadowSync_Init(
		ParticleShadowSync*     self );

	/** Stg_Class_Print() implementation */
	void _ParticleShadowSync_Print( void* pCommsHandler, Stream* stream );
	
	void _ParticleShadowSync_AssignFromXML( void* pCommsHandler, Stg_ComponentFactory* cf, void* data );
	
	void _ParticleShadowSync_Build( void* pCommsHandler, void *data );
	
	void _ParticleShadowSync_Initialise( void* pCommsHandler, void *data );
	
	void _ParticleShadowSync_Execute( void* pCommsHandler, void *data );

	void _ParticleShadowSync_Destroy( void* pCommsHandler, void *data );
	
	/** Copy */
	#define ParticleShadowSync_Copy( self ) \
		(ParticleShadowSync*)ParticleShadowSync_CopyFunc( self, NULL, False, NULL, NULL )
	#define ParticleShadowSync_DeepCopy( self ) \
		(ParticleShadowSync*)ParticleShadowSync_CopyFunc( self, NULL, True, NULL, NULL )
	
	void* _ParticleShadowSync_CopyFunc( void* ParticleShadowSync, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	/** Stg_Class_Delete() implementation */
	void _ParticleShadowSync_Delete(void* pCommsHandler );
	
	/* --- Public virtual function interfaces --- */
	
	/** Handle particle movement between processors */
	void ParticleShadowSync_HandleParticleMovementBetweenProcs( ParticleCommHandler* pCommsHandler );

	/* --- virtual function implementations --- */

	/* +++ Global fallback method related +++ */
	void ParticleShadowSync_DoGlobalFallbackCommunication( ParticleShadowSync* self );

	void ParticleShadowSync_FindParticlesThatHaveMovedOutsideMyDomain( ParticleShadowSync* self );

	void ParticleShadowSync_GetCountOfParticlesOutsideDomainPerProcessor(
		ParticleShadowSync*	self,
		Particle_Index**	globalParticlesOutsideDomainCountsPtr,
		Particle_Index*		maxGlobalParticlesOutsideDomainCountPtr,
		Particle_Index*		globalParticlesOutsideDomainTotalPtr );
		
	void ParticleShadowSync_EnsureParticleCountLeavingDomainsEqualsCountEnteringGlobally( ParticleShadowSync* self );

	void _ParticleShadowSync_FinishReceiveOfIncomingParticleCounts( ParticleCommHandler* self );
	void _ParticleShadowSync_BeginReceiveOfIncomingParticles( ParticleCommHandler* self );
	void _ParticleShadowSync_SendParticleTotalsInShadowCellsToNbrs( ParticleCommHandler* self );
	void _ParticleShadowSync_SendShadowParticles( ParticleCommHandler* self );
	void _ParticleShadowSync_FinishReceiveOfIncomingParticles( ParticleCommHandler* pCommHandler );

	/* --- private functions --- */

#endif

