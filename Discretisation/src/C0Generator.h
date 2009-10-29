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
** $Id: C0Generator.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Discretisaton_Mesh_C0Generator_h__
#define __Discretisaton_Mesh_C0Generator_h__

	/** Textual name of this class */
	extern const Type C0Generator_Type;

	/** Virtual function types */

	/** C0Generator class contents */
	#define __C0Generator		\
		/* General info */		\
		__MeshGenerator			\
						\
		/* Virtual info */		\
						\
		/* C0Generator info */		\
		Mesh*		elMesh;

	struct C0Generator { __C0Generator };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define C0GENERATOR_DEFARGS	\
		MESHGENERATOR_DEFARGS

	#define C0GENERATOR_PASSARGS	\
		MESHGENERATOR_PASSARGS

	C0Generator* C0Generator_New( Name name );
	C0Generator* _C0Generator_New( C0GENERATOR_DEFARGS );
	void _C0Generator_Init( C0Generator* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _C0Generator_Delete( void* generator );
	void _C0Generator_Print( void* generator, Stream* stream );
	void _C0Generator_AssignFromXML( void* generator, Stg_ComponentFactory* cf, void* data );
	void _C0Generator_Build( void* generator, void* data );
	void _C0Generator_Initialise( void* generator, void* data );
	void _C0Generator_Execute( void* generator, void* data );
	void _C0Generator_Destroy( void* generator, void* data );

	void C0Generator_Generate( void* generator, void* _mesh );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void C0Generator_SetElementMesh( void* generator, void* mesh );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	void C0Generator_BuildTopology( C0Generator* self, FeMesh* mesh );
	void C0Generator_BuildGeometry( C0Generator* self, FeMesh* mesh );
	void C0Generator_BuildElementTypes( C0Generator* self, FeMesh* mesh );

#endif /* __Discretisaton_Mesh_C0Generator_h__ */
