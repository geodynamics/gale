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
** $Id: Mesh_CentroidAlgorithms.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Domain_Mesh_Mesh_CentroidAlgorithms_h__
#define __Domain_Mesh_Mesh_CentroidAlgorithms_h__

	/** Textual name of this class */
	extern const Type Mesh_CentroidAlgorithms_Type;

	/** Virtual function types */

	/** Class contents */
	#define __Mesh_CentroidAlgorithms		\
		/* General info */			\
		__Mesh_Algorithms			\
							\
		/* Virtual info */			\
							\
		/* Mesh_CentroidAlgorithms info */	\
		Mesh*			elMesh;

	struct Mesh_CentroidAlgorithms { __Mesh_CentroidAlgorithms };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define MESH_HEXALGORITHMS_DEFARGS \
		MESH_ALGORITHMS_DEFARGS

	#define MESH_HEXALGORITHMS_PASSARGS \
		MESH_ALGORITHMS_PASSARGS

	Mesh_CentroidAlgorithms* Mesh_CentroidAlgorithms_New( Name name, AbstractContext* context );
	Mesh_CentroidAlgorithms* _Mesh_CentroidAlgorithms_New( MESH_HEXALGORITHMS_DEFARGS );
	void _Mesh_CentroidAlgorithms_Init( Mesh_CentroidAlgorithms* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _Mesh_CentroidAlgorithms_Delete( void* centroidAlgorithms );
	void _Mesh_CentroidAlgorithms_Print( void* centroidAlgorithms, Stream* stream );
	void _Mesh_CentroidAlgorithms_AssignFromXML( void* centroidAlgorithms, Stg_ComponentFactory* cf, void* data );
	void _Mesh_CentroidAlgorithms_Build( void* centroidAlgorithms, void* data );
	void _Mesh_CentroidAlgorithms_Initialise( void* centroidAlgorithms, void* data );
	void _Mesh_CentroidAlgorithms_Execute( void* centroidAlgorithms, void* data );
	void _Mesh_CentroidAlgorithms_Destroy( void* centroidAlgorithms, void* data );

	void Mesh_CentroidAlgorithms_Update( void* centroidAlgorithms );
	unsigned Mesh_CentroidAlgorithms_NearestVertex( void* centroidAlgorithms, double* point );
	Bool Mesh_CentroidAlgorithms_Search( void* centroidAlgorithms, double* point, 
					     MeshTopology_Dim* dim, unsigned* ind );
	Bool Mesh_CentroidAlgorithms_SearchElements( void* centroidAlgorithms, double* point, 
						     unsigned* elInd );
	void Mesh_CentroidAlgorithms_GetLocalCoordRange( void* algorithms, double* min, double* max );
	void Mesh_CentroidAlgorithms_GetDomainCoordRange( void* algorithms, double* min, double* max );
	void Mesh_CentroidAlgorithms_GetGlobalCoordRange( void* algorithms, double* min, double* max );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void Mesh_CentroidAlgorithms_SetElementMesh( void* centroidAlgorithms, void* mesh );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

#endif /* __Domain_Mesh_Mesh_CentroidAlgorithms_h__ */
