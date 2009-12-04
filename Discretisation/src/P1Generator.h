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
** $Id: P1Generator.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_Discretisaton_P1Generator_h__
#define __StgFEM_Discretisaton_P1Generator_h__

	/** Textual name of this class */
	extern const Type P1Generator_Type;

	/** Virtual function types */

	/** P1Generator class contents */
	#define __P1Generator			\
		/* General info */		\
		__MeshGenerator			\
						\
		/* Virtual info */		\
						\
		/* P1Generator info */		\
		FeMesh*		elMesh;

	struct P1Generator { __P1Generator };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define P1GENERATOR_DEFARGS	\
		MESHGENERATOR_DEFARGS

	#define P1GENERATOR_PASSARGS	\
		MESHGENERATOR_PASSARGS

	P1Generator* P1Generator_New( Name name );
	P1Generator* _P1Generator_New( P1GENERATOR_DEFARGS );
	void _P1Generator_Init( P1Generator* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _P1Generator_Delete( void* generator );
	void _P1Generator_Print( void* generator, Stream* stream );
	void _P1Generator_AssignFromXML( void* generator, Stg_ComponentFactory* cf, void* data );
	void _P1Generator_Build( void* generator, void* data );
	void _P1Generator_Initialise( void* generator, void* data );
	void _P1Generator_Execute( void* generator, void* data );
	void _P1Generator_Destroy( void* generator, void* data );

	void P1Generator_Generate( void* generator, void* _mesh );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void P1Generator_SetElementMesh( void* generator, void* mesh );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	void P1Generator_BuildTopology( P1Generator* self, FeMesh* mesh );
	void P1Generator_BuildGeometry( P1Generator* self, FeMesh* mesh );
	void P1Generator_BuildElementTypes( P1Generator* self, FeMesh* mesh );

#endif /* __StgFEM_Discretisaton_P1Generator_h__ */
