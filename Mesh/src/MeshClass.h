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
** $Id: MeshClass.h 4184 2007-09-25 07:54:17Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Domain_Mesh_Mesh_h__
#define __Domain_Mesh_Mesh_h__

	/** Textual name of this class */
	extern const Type Mesh_Type;

	/** Virtual function types */

	/** Class contents */
	#define __Mesh						\
		/* General info */				\
		__Stg_Component					\
								\
		/* Virtual info */				\
								\
		/* Mesh info */					\
		MeshTopology*			topo;		\
		double**			verts;		\
								\
		List*				vars;		\
								\
		double				minSep;		\
		double*				minAxialSep;	\
		double*				minLocalCrd;	\
		double*				maxLocalCrd;	\
		double*				minDomainCrd;	\
		double*				maxDomainCrd;	\
		double*				minGlobalCrd;	\
		double*				maxGlobalCrd;	\
								\
		Mesh_Algorithms*		algorithms;	\
		unsigned			nElTypes;	\
		Mesh_ElementType**		elTypes;	\
		unsigned*			elTypeMap;	\
								\
		UIntMap*			topoDataSizes;	\
		ExtensionManager**		topoDataInfos;	\
		void**				topoDatas;	\
		ExtensionManager*		info;		\
								\
		MeshGenerator*			generator;	\
		ExtensionManager_Register*	emReg;

	struct Mesh { __Mesh };

	/* Checkpoint file version enum */
	typedef enum MeshCheckpointFileVersion {
                MeshCHECKPOINT_V1 = 1,      /** Original checkpointing format   */
                MeshCHECKPOINT_V2           /** * No longer store nodes within checkpoint files
                                            * now store attributes including
                                                 checkpoint version
                                                 number of dimensions
                                                 cartesion mesh size
                                            * now store mesh connectivity for xmdf generator
                                            * mesh vertex hdf group renamed to 'vertices' (from 'data') */
	} MeshCheckpointFileVersion;

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define MESH_DEFARGS		\
		STG_COMPONENT_DEFARGS

	#define MESH_PASSARGS		\
		STG_COMPONENT_PASSARGS

	Mesh* Mesh_New( Name name );
	Mesh* _Mesh_New( MESH_DEFARGS );
	void _Mesh_Init( Mesh* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _Mesh_Delete( void* mesh );
	void _Mesh_Print( void* mesh, Stream* stream );
	void _Mesh_Construct( void* mesh, Stg_ComponentFactory* cf, void* data );
	void _Mesh_Build( void* mesh, void* data );
	void _Mesh_Initialise( void* mesh, void* data );
	void _Mesh_Execute( void* mesh, void* data );
	void _Mesh_Destroy( void* mesh, void* data );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void Mesh_SetGenerator( void* mesh, void* generator );
	void Mesh_SetAlgorithms( void* mesh, void* algorithms );

	void Mesh_SetExtensionManagerRegister( void* mesh, void* extMgrReg );
	void Mesh_SetTopologyDataSize( void* mesh, MeshTopology_Dim dim, unsigned size );
	void* Mesh_GetTopologyData( void* mesh, MeshTopology_Dim dim );

	unsigned Mesh_GetDimSize( void* mesh );
	unsigned Mesh_GetGlobalSize( void* mesh, MeshTopology_Dim dim );
	unsigned Mesh_GetLocalSize( void* mesh, MeshTopology_Dim dim );
	unsigned Mesh_GetRemoteSize( void* mesh, MeshTopology_Dim dim );
	unsigned Mesh_GetDomainSize( void* mesh, MeshTopology_Dim dim );
	unsigned Mesh_GetSharedSize( void* mesh, MeshTopology_Dim dim );
	MeshTopology* Mesh_GetTopology( void* mesh );
	Sync* Mesh_GetSync( void* mesh, MeshTopology_Dim dim );

	Bool Mesh_GlobalToDomain( void* mesh, MeshTopology_Dim dim, unsigned global, unsigned* domain );
	unsigned Mesh_DomainToGlobal( void* mesh, MeshTopology_Dim dim, unsigned domain );
	Bool Mesh_LocalToShared( void* meshTopology, MeshTopology_Dim dim, unsigned domain, unsigned* shared );
	unsigned Mesh_SharedToLocal( void* meshTopology, MeshTopology_Dim dim, unsigned shared );

	unsigned Mesh_GetOwner( void* mesh, MeshTopology_Dim dim, unsigned remote );
	void Mesh_GetSharers( void* mesh, MeshTopology_Dim dim, unsigned shared, 
			      unsigned* nSharers, unsigned** sharers );

	Bool Mesh_HasIncidence( void* mesh, MeshTopology_Dim fromDim, MeshTopology_Dim toDim );
	unsigned Mesh_GetIncidenceSize( void* mesh, MeshTopology_Dim fromDim, unsigned fromInd, 
					MeshTopology_Dim toDim );
	void Mesh_GetIncidence( void* mesh, MeshTopology_Dim fromDim, unsigned fromInd, MeshTopology_Dim toDim, 
				IArray* inc );

	unsigned Mesh_NearestVertex( void* mesh, double* point );
	Bool Mesh_Search( void* mesh, double* point, 
			  MeshTopology_Dim* dim, unsigned* ind );

	Bool Mesh_SearchElements( void* mesh, double* point, unsigned* elInd );
	/* Mesh_SearchElements (
	 * mesh -- is a mesh
	 * point -- is a global coordinate
	 * elInd -- will be filled in by a local elementID
	 * )
	 * returns:
	 * False if the point is not in the DOMAIN space of the proc 
	 * True if the point is in the DOMAIN space
	 */

	Bool Mesh_ElementHasPoint( void* mesh, unsigned element, double* point, 
				   MeshTopology_Dim* dim, unsigned* ind );
	Mesh_ElementType* Mesh_GetElementType( void* mesh, unsigned element );

	Comm* Mesh_GetCommTopology( void* mesh, MeshTopology_Dim dim );
	double* Mesh_GetVertex( void* mesh, unsigned domain );

	Bool Mesh_HasExtension( void* mesh, const char* name );
	#define Mesh_GetExtension( mesh, type, name ) \
		(type)_Mesh_GetExtension( mesh, name )
	void* _Mesh_GetExtension( void* mesh, const char* name );

	void Mesh_GetMinimumSeparation( void* mesh, double* minSep, double* axial );
	void Mesh_GetLocalCoordRange( void* mesh, double* min, double* max );
	void Mesh_GetDomainCoordRange( void* mesh, double* min, double* max );
	void Mesh_GetGlobalCoordRange( void* mesh, double* min, double* max );

	void Mesh_DeformationUpdate( void* mesh );
	void Mesh_Sync( void* mesh );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	void Mesh_Destruct( Mesh* self );

#endif /* __Domain_Mesh_Mesh_h__ */
