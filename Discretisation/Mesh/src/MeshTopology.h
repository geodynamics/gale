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
**
** Assumptions:
**
** Invariants:
**
** Comments:
**
** $Id: MeshTopology.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Discretisaton_Mesh_MeshTopology_h__
#define __Discretisaton_Mesh_MeshTopology_h__

	/** Textual name of this class */
	extern const Type MeshTopology_Type;

	/** Virtual function types */

	/** MeshTopology class contents */
	typedef enum {
		MT_VERTEX = 0, 
		MT_EDGE = 1, 
		MT_FACE = 2, 
		MT_VOLUME = 3
	} MeshTopology_Dim;

	#define __MeshTopology				\
		/* General info */			\
		__Stg_Component				\
							\
		/* Virtual info */			\
							\
		/* MeshTopology info */			\
		unsigned		nDims;		\
		unsigned		nTDims;		\
		Decomp_Sync**		domains;	\
		unsigned		nBndVerts;	\
		unsigned*		bndVerts;	\
		unsigned		shadowDepth;	\
		unsigned***		nIncEls;	\
		unsigned****		incEls;

	struct MeshTopology { __MeshTopology };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define MESHTOPOLOGY_DEFARGS	\
		STG_COMPONENT_DEFARGS

	#define MESHTOPOLOGY_PASSARGS	\
		STG_COMPONENT_PASSARGS

	MeshTopology* MeshTopology_New( Name name );
	MeshTopology* _MeshTopology_New( MESHTOPOLOGY_DEFARGS );
	void _MeshTopology_Init( MeshTopology* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _MeshTopology_Delete( void* topo );
	void _MeshTopology_Print( void* topo, Stream* stream );
	void _MeshTopology_Construct( void* topo, Stg_ComponentFactory* cf, void* data );
	void _MeshTopology_Build( void* topo, void* data );
	void _MeshTopology_Initialise( void* topo, void* data );
	void _MeshTopology_Execute( void* topo, void* data );
	void _MeshTopology_Destroy( void* topo, void* data );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void MeshTopology_SetDimSize( void* topo, unsigned nDims );
	void MeshTopology_SetElements( void* topo, MeshTopology_Dim dim, unsigned nElements, unsigned* elements );
	void MeshTopology_SetIncidence( void* topo, MeshTopology_Dim fromDim, MeshTopology_Dim toDim, 
					unsigned* nIncEls, unsigned** incEls );
	void MeshTopology_SetBoundaryVertices( void* topo, unsigned nVerts, unsigned* verts );
	void MeshTopology_SetShadowDepth( void* topo, unsigned depth );
	void MeshTopology_SetSync( void* topo, MeshTopology_Dim dim, Decomp_Sync* sync );

	unsigned MeshTopology_GetGlobalSize( void* mesh, MeshTopology_Dim dim );
	unsigned MeshTopology_GetLocalSize( void* meshTopology, MeshTopology_Dim dim );
	unsigned MeshTopology_GetRemoteSize( void* meshTopology, MeshTopology_Dim dim );
	unsigned MeshTopology_GetDomainSize( void* meshTopology, MeshTopology_Dim dim );
	unsigned MeshTopology_GetSharedSize( void* meshTopology, MeshTopology_Dim dim );
	void MeshTopology_GetLocalElements( void* meshTopology, MeshTopology_Dim dim, unsigned* nEls, unsigned** els );
	void MeshTopology_GetRemoteElements( void* meshTopology, MeshTopology_Dim dim, unsigned* nEls, unsigned** els );
	unsigned MeshTopology_GetOwner( void* meshTopology, MeshTopology_Dim dim, unsigned remote );
	void MeshTopology_GetSharers( void* meshTopology, MeshTopology_Dim dim, unsigned shared, 
				      unsigned* nSharers, unsigned** sharers );
	CommTopology* MeshTopology_GetCommTopology( void* meshTopology, MeshTopology_Dim dim );
	Decomp_Sync* MeshTopology_GetSync( void* meshTopology, MeshTopology_Dim dim );

	void MeshTopology_Complete( void* topo );
	void MeshTopology_Invert( void* topo, MeshTopology_Dim fromDim, MeshTopology_Dim toDim );
	void MeshTopology_Cascade( void* topo, MeshTopology_Dim fromDim, MeshTopology_Dim toDim );
	void MeshTopology_Neighbourhood( void* topo, MeshTopology_Dim dim );

	Bool MeshTopology_GlobalToDomain( void* meshTopology, MeshTopology_Dim dim, unsigned global, unsigned* domain );
	unsigned MeshTopology_DomainToGlobal( void* meshTopology, MeshTopology_Dim dim, unsigned domain );
	Bool MeshTopology_DomainToShared( void* meshTopology, MeshTopology_Dim dim, unsigned domain, unsigned* shared );
	unsigned MeshTopology_SharedToDomain( void* meshTopology, MeshTopology_Dim dim, unsigned shared );

	Bool MeshTopology_HasIncidence( void* meshTopology, MeshTopology_Dim fromDim, MeshTopology_Dim toDim );
	unsigned MeshTopology_GetIncidenceSize( void* meshTopology, MeshTopology_Dim fromDim, unsigned fromInd, 
						MeshTopology_Dim toDim );
	void MeshTopology_GetIncidence( void* meshTopology, MeshTopology_Dim fromDim, unsigned fromInd, MeshTopology_Dim toDim, 
					unsigned* nInc, unsigned** inc );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	void MeshTopology_ClearElements( MeshTopology* self, MeshTopology_Dim dim );
	void MeshTopology_ClearIncidence( MeshTopology* self, MeshTopology_Dim dim );
	void MeshTopology_CommUnion( MeshTopology* self );
	void MeshTopology_BuildShadows( MeshTopology* self );
	void MeshTopology_BuildShadowedElements( MeshTopology* self, MeshTopology_Dim dim, 
						 unsigned** nShadowedEls, unsigned*** shadowedEls, 
						 unsigned** nExternalEls, unsigned*** externalEls );
	void MeshTopology_BuildTopShadowedElements( MeshTopology* self, unsigned** nShadowedEls, unsigned*** shadowedEls, 
						    unsigned** nExternalEls, unsigned*** externalEls );
	void MeshTopology_ExpandShadows( MeshTopology* self, IndexSet* shadows, IndexSet* boundary );
	void MeshTopology_BuildMidShadowedElements( MeshTopology* self, MeshTopology_Dim dim, 
						    unsigned** nShadowedEls, unsigned*** shadowedEls, 
						    unsigned** nExternalEls, unsigned*** externalEls );
	void MeshTopology_BuildShadowedIncidence( MeshTopology* self, MeshTopology_Dim dim, 
						  unsigned* nShadowedEls, unsigned** shadowedEls, 
						  unsigned* nShadowedInc, Stg_Byte** shadowedInc );
	void MeshTopology_ExchangeShadowedElements( MeshTopology* self, MeshTopology_Dim dim, 
						    unsigned* nShadowedEls, unsigned** shadowedEls );
	void MeshTopology_ExchangeExternalElements( MeshTopology* self, MeshTopology_Dim dim, 
						    unsigned* nExternalEls, unsigned** externalEls, 
						    unsigned** nShadowedEls, unsigned*** shadowedEls );
	void MeshTopology_ExchangeShadowedIncidence( MeshTopology* self, MeshTopology_Dim dim, 
						     unsigned* nShadowedInc, Stg_Byte** shadowedInc );
	void MeshTopology_ExchangeElements( MeshTopology* self, MeshTopology_Dim dim, 
					    unsigned* nSendElements, unsigned** sendElements, 
					    unsigned* nRecvElements, unsigned** recvElements );

	void MeshTopology_PickleIncidence( MeshTopology* self, MeshTopology_Dim fromDim, MeshTopology_Dim toDim, 
					   unsigned nEls, unsigned* els, 
					   unsigned* nBytes, Stg_Byte** bytes );
	void MeshTopology_UnpickleIncidence( MeshTopology* self, MeshTopology_Dim fromDim, MeshTopology_Dim toDim, 
					     unsigned nBytes, Stg_Byte* bytes );

	void MeshTopology_Destruct( MeshTopology* self );

#ifndef NDEBUG
	Bool MeshTopology_ValidateElements( MeshTopology* self, MeshTopology_Dim dim, unsigned nEls, unsigned* els );
#endif

#endif /* __Discretisaton_Mesh_MeshTopology_h__ */
