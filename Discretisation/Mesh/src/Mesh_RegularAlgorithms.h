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
** $Id: Mesh_RegularAlgorithms.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Discretisaton_Mesh_Mesh_RegularAlgorithms_h__
#define __Discretisaton_Mesh_Mesh_RegularAlgorithms_h__

	/** Textual name of this class */
	extern const Type Mesh_RegularAlgorithms_Type;

	/** Virtual function types */

	/** Class contents */
	#define __Mesh_RegularAlgorithms		\
		/* General info */			\
		__Mesh_Algorithms			\
							\
		/* Virtual info */			\
							\
		/* Mesh_RegularAlgorithms info */	\
		double*		sep;			\
		double*		minCrd;			\
		double*		maxCrd;

	struct Mesh_RegularAlgorithms { __Mesh_RegularAlgorithms };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define MESH_REGULARALGORITHMS_DEFARGS \
		MESH_ALGORITHMS_DEFARGS

	#define MESH_REGULARALGORITHMS_PASSARGS \
		MESH_ALGORITHMS_PASSARGS

	Mesh_RegularAlgorithms* Mesh_RegularAlgorithms_New( Name name );
	Mesh_RegularAlgorithms* _Mesh_RegularAlgorithms_New( MESH_REGULARALGORITHMS_DEFARGS );
	void _Mesh_RegularAlgorithms_Init( Mesh_RegularAlgorithms* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _Mesh_RegularAlgorithms_Delete( void* algorithms );
	void _Mesh_RegularAlgorithms_Print( void* algorithms, Stream* stream );
	void _Mesh_RegularAlgorithms_Construct( void* algorithms, Stg_ComponentFactory* cf, void* data );
	void _Mesh_RegularAlgorithms_Build( void* algorithms, void* data );
	void _Mesh_RegularAlgorithms_Initialise( void* algorithms, void* data );
	void _Mesh_RegularAlgorithms_Execute( void* algorithms, void* data );
	void _Mesh_RegularAlgorithms_Destroy( void* algorithms, void* data );

	void Mesh_RegularAlgorithms_SetMesh( void* algorithms, void* mesh );
	void Mesh_RegularAlgorithms_Update( void* algorithms );
	Bool Mesh_RegularAlgorithms_Search( void* algorithms, void* mesh, double* point, 
					    MeshTopology_Dim* dim, unsigned* ind );
	Bool Mesh_RegularAlgorithms_SearchElements( void* algorithms, double* point, unsigned* elInd );
	double Mesh_RegularAlgorithms_GetMinimumSeparation( void* algorithms, void* mesh, double* perDim );
	void Mesh_RegularAlgorithms_GetLocalCoordRange( void* algorithms, void* mesh, double* min, double* max );
	void Mesh_RegularAlgorithms_GetDomainCoordRange( void* algorithms, void* mesh, double* min, double* max );
	void Mesh_RegularAlgorithms_GetGlobalCoordRange( void* algorithms, void* mesh, double* min, double* max );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	void Mesh_RegularAlgorithms_Destruct( Mesh_RegularAlgorithms* self );

#endif /* __Discretisaton_Mesh_Mesh_RegularAlgorithms_h__ */
