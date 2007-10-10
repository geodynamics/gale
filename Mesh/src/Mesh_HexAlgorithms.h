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
** $Id: Mesh_HexAlgorithms.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Domain_Mesh_Mesh_HexAlgorithms_h__
#define __Domain_Mesh_Mesh_HexAlgorithms_h__

	/** Textual name of this class */
	extern const Type Mesh_HexAlgorithms_Type;

	/** Virtual function types */

	/** Class contents */
	#define __Mesh_HexAlgorithms			\
		/* General info */			\
		__Mesh_Algorithms			\
							\
		/* Virtual info */			\
							\
		/* Mesh_HexAlgorithms info */

	struct Mesh_HexAlgorithms { __Mesh_HexAlgorithms };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define MESH_HEXALGORITHMS_DEFARGS \
		MESH_ALGORITHMS_DEFARGS

	#define MESH_HEXALGORITHMS_PASSARGS \
		MESH_ALGORITHMS_PASSARGS

	Mesh_HexAlgorithms* Mesh_HexAlgorithms_New( Name name );
	Mesh_HexAlgorithms* _Mesh_HexAlgorithms_New( MESH_HEXALGORITHMS_DEFARGS );
	void _Mesh_HexAlgorithms_Init( Mesh_HexAlgorithms* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _Mesh_HexAlgorithms_Delete( void* hexAlgorithms );
	void _Mesh_HexAlgorithms_Print( void* hexAlgorithms, Stream* stream );
	void _Mesh_HexAlgorithms_Construct( void* hexAlgorithms, Stg_ComponentFactory* cf, void* data );
	void _Mesh_HexAlgorithms_Build( void* hexAlgorithms, void* data );
	void _Mesh_HexAlgorithms_Initialise( void* hexAlgorithms, void* data );
	void _Mesh_HexAlgorithms_Execute( void* hexAlgorithms, void* data );
	void _Mesh_HexAlgorithms_Destroy( void* hexAlgorithms, void* data );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

#endif /* __Domain_Mesh_Mesh_HexAlgorithms_h__ */
