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
** $Id: ParticleCommHandler.h 4119 2007-05-22 07:35:46Z RaquibulHassan $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	
#ifndef __Discretisation_Swarm_ParticleCommHandler_h__
#define __Discretisation_Swarm_ParticleCommHandler_h__

	
	/** Textual name of this class */
	extern const Type ParticleCommHandler_Type;

	typedef void (ParticleCommHandler_CommFunction)	( ParticleCommHandler * );
	typedef void (ParticleCommHandler_AllocateOutgoingCountArrays) ( ParticleCommHandler * );
	typedef void (ParticleCommHandler_AllocateOutgoingParticleArrays) ( ParticleCommHandler * );
	typedef void (ParticleCommHandler_FreeOutgoingArrays) ( ParticleCommHandler * );
	
	typedef void (ParticleCommHandler_AllocateIncomingCountArrays) ( ParticleCommHandler * );
	typedef void (ParticleCommHandler_AllocateIncomingParticleArrays) ( ParticleCommHandler * );
	typedef void (ParticleCommHandler_FreeIncomingArrays) ( ParticleCommHandler * );
	
	typedef void (ParticleCommHandler_BeginReceiveOfIncomingParticleCounts) ( ParticleCommHandler * );
	typedef void (ParticleCommHandler_FinishReceiveOfIncomingParticleCounts) ( ParticleCommHandler * );
	
	typedef void (ParticleCommHandler_BeginReceiveOfIncomingParticles) ( ParticleCommHandler * );
	typedef void (ParticleCommHandler_FinishReceiveOfIncomingParticlesAndUpdateIndices) ( ParticleCommHandler * );
	typedef void (ParticleCommHandler_SendOutgoingParticleCounts) ( ParticleCommHandler * );
	typedef void (ParticleCommHandler_BeginSendingParticles) ( ParticleCommHandler * );
	typedef void (ParticleCommHandler_ConfirmOutgoingSendsCompleted) ( ParticleCommHandler * );


	#define __ParticleCommHandler \
		__Stg_Component \
		/* Virtual info */ \
		ParticleCommHandler_AllocateOutgoingCountArrays						*allocateOutgoingCountArrays; \
		ParticleCommHandler_AllocateOutgoingParticleArrays						*allocateOutgoingParticleArrays; \
		ParticleCommHandler_FreeOutgoingArrays									*freeOutgoingArrays; \
		ParticleCommHandler_AllocateIncomingCountArrays						*allocateIncomingCountArrays; \
		ParticleCommHandler_AllocateIncomingParticleArrays						*allocateIncomingParticleArrays; \
		ParticleCommHandler_FreeIncomingArrays									*freeIncomingArrays; \
		ParticleCommHandler_BeginReceiveOfIncomingParticleCounts				*beginReceiveOfIncomingParticleCounts; \
		ParticleCommHandler_FinishReceiveOfIncomingParticleCounts				*finishReceiveOfIncomingParticleCounts; \
		ParticleCommHandler_BeginReceiveOfIncomingParticles					*beginReceiveOfIncomingParticles; \
		ParticleCommHandler_FinishReceiveOfIncomingParticlesAndUpdateIndices	*finishReceiveOfIncomingParticlesAndUpdateIndices; \
		ParticleCommHandler_SendOutgoingParticleCounts							*sendOutgoingParticleCounts; \
		ParticleCommHandler_BeginSendingParticles								*beginSendingParticles; \
		ParticleCommHandler_ConfirmOutgoingSendsCompleted						*confirmOutgoingSendsCompleted; \
		ParticleCommHandler_CommFunction										*_commFunction;\
		/* Member info */ \
		Stream*				debug; \
		Swarm*				swarm; \
		Particle_Index* 		shadowParticlesLeavingMeIndices; \
		Index				shadowParticlesLeavingMeTotalCount; \
		Index				shadowParticlesLeavingMeUnfilledCount; \
		Index				currShadowParticleLeavingMeIndex; \
		Index				currParticleLeavingMyDomainIndex; \
		\
		Particle_Index* 		particlesOutsideDomainIndices; \
		Index				particlesOutsideDomainTotalCount; \
		Index				particlesOutsideDomainUnfilledCount; \
		\
		/** cnts of [nbr][st_cell] outgoing particles */ \
		Particle_Index**		shadowParticlesLeavingMeCountsPerCell; \
		/** cnts of [nbr] total outgoing particles via my shadow cells to nbr procs */ \
		Particle_Index*			shadowParticlesLeavingMeTotalCounts; \
		/** transfer array [nbr] of particles to send */ \
		Particle**			shadowParticlesLeavingMe; \
		MPI_Request**			shadowParticlesLeavingMeHandles; \
		/** cnts of [nbr][st_cell] incoming particles */ \
		Particle_Index**		particlesArrivingFromNbrShadowCellCounts; \
		MPI_Request**			particlesArrivingFromNbrShadowCellCountsHandles; \
		/** cnts of [nbr] total incoming particles via shadow cells of nbr procs */\
		Particle_Index*			particlesArrivingFromNbrShadowCellsTotalCounts; \
		/** transfer array [nbr] of particles to recv */ \
		Particle**			particlesArrivingFromNbrShadowCells; \
		MPI_Request**			particlesArrivingFromNbrShadowCellsHandles;

	struct ParticleCommHandler { __ParticleCommHandler };	

	/* --- virtual functions --- */

	/** Constructor interface */
	
	/** Private Constructor */
	ParticleCommHandler* _ParticleCommHandler_New( 
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
		Name                                                            name,
		Bool                                                            initFlag,
		ParticleCommHandler_AllocateOutgoingCountArrays						*_allocateOutgoingCountArrays,
		ParticleCommHandler_AllocateOutgoingParticleArrays						*_allocateOutgoingParticleArrays,
		ParticleCommHandler_FreeOutgoingArrays									*_freeOutgoingArrays,
		ParticleCommHandler_AllocateIncomingCountArrays						*_allocateIncomingCountArrays,
		ParticleCommHandler_AllocateIncomingParticleArrays						*_allocateIncomingParticleArrays,
		ParticleCommHandler_FreeIncomingArrays									*_freeIncomingArrays,
		ParticleCommHandler_BeginReceiveOfIncomingParticleCounts				*_beginReceiveOfIncomingParticleCounts,
		ParticleCommHandler_FinishReceiveOfIncomingParticleCounts				*_finishReceiveOfIncomingParticleCounts,
		ParticleCommHandler_BeginReceiveOfIncomingParticles					*_beginReceiveOfIncomingParticles,
		ParticleCommHandler_FinishReceiveOfIncomingParticlesAndUpdateIndices	*_finishReceiveOfIncomingParticlesAndUpdateIndices,
		ParticleCommHandler_SendOutgoingParticleCounts							*_sendOutgoingParticleCounts,
		ParticleCommHandler_BeginSendingParticles								*_beginSendingParticle,
		ParticleCommHandler_ConfirmOutgoingSendsCompleted						*_confirmOutgoingSendsCompleted,
		ParticleCommHandler_CommFunction										*_commFunction );
	
	/** Variable initialiser */
	void _ParticleCommHandler_Init(
		ParticleCommHandler*     self
		);

	/** Stg_Class_Print() implementation */
	void _ParticleCommHandler_Print( void* pCommsHandler, Stream* stream );
	
	void _ParticleCommHandler_Construct( void* pCommsHandler, Stg_ComponentFactory* cf, void* data );
	
	void _ParticleCommHandler_Build( void* pCommsHandler, void *data );
	
	void _ParticleCommHandler_Initialise( void* pCommsHandler, void *data );
	
	void _ParticleCommHandler_Execute( void* pCommsHandler, void *data );

	void _ParticleCommHandler_Destroy( void* pCommsHandler, void *data );
	
	/** Copy */
	#define ParticleCommHandler_Copy( self ) \
		(ParticleCommHandler*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define ParticleCommHandler_DeepCopy( self ) \
		(ParticleCommHandler*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
	
	void* _ParticleCommHandler_Copy( void* particleCommHandler, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	/** Stg_Class_Delete() implementation */
	void _ParticleCommHandler_Delete(void* pCommsHandler );
	
	/* --- Public virtual function interfaces --- */
	
	/* --- virtual function implementations --- */

	/* --- private functions --- */

	/* +++ Comms via shadow cell related +++ */

	void _ParticleCommHandler_AllocateOutgoingCountArrays( ParticleCommHandler* self );
	void _ParticleCommHandler_AllocateOutgoingParticlesArrays( ParticleCommHandler* self );
	void _ParticleCommHandler_FreeOutgoingArrays( ParticleCommHandler* self );

	void _ParticleCommHandler_AllocateIncomingCountArrays( ParticleCommHandler* self );
	void _ParticleCommHandler_AllocateIncomingParticlesArrays( ParticleCommHandler* self );
	void _ParticleCommHandler_FreeIncomingArrays( ParticleCommHandler* self );

	void _ParticleCommHandler_BeginReceiveOfIncomingParticleCounts( ParticleCommHandler* self );
	void _ParticleCommHandler_FinishReceiveOfIncomingParticleCounts( ParticleCommHandler* self );

	void _ParticleCommHandler_BeginReceiveOfIncomingParticles( ParticleCommHandler* self );
	void _ParticleCommHandler_FinishReceiveAndUpdateShadowParticlesEnteringMyDomain( ParticleCommHandler* self );

	void _ParticleCommHandler_SendParticleTotalsInShadowCellsToNbrs( ParticleCommHandler* self );

	void _ParticleCommHandler_BeginSendingParticlesInShadowCellsToNbrs( ParticleCommHandler* self );
	void _ParticleCommHandler_ConfirmOutgoingSendsCompleted( ParticleCommHandler* self );

	void _ParticleCommHandler_ZeroShadowCommStrategyCounters( ParticleCommHandler* self );


	/* +++ Statistics Printing +++ */
	void _ParticleCommHandler_PrintCommunicationVolumeStats( ParticleCommHandler* self, double startTime, Stream* info );

#endif
