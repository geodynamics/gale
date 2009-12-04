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
** $Id: Inner2DGenerator.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_Discretisaton_Inner2DGenerator_h__
#define __StgFEM_Discretisaton_Inner2DGenerator_h__

	/** Textual name of this class */
	extern const Type Inner2DGenerator_Type;

	/** Virtual function types */

	/** Inner2DGenerator class contents */
	#define __Inner2DGenerator		\
		/* General info */		\
		__MeshGenerator			\
						\
		/* Virtual info */		\
						\
		/* Inner2DGenerator info */		\
		Mesh*		elMesh;

	struct Inner2DGenerator { __Inner2DGenerator };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/



	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define INNER2DGENERATOR_DEFARGS \
                MESHGENERATOR_DEFARGS

	#define INNER2DGENERATOR_PASSARGS \
                MESHGENERATOR_PASSARGS

	Inner2DGenerator* Inner2DGenerator_New( Name name, AbstractContext* context );
	Inner2DGenerator* _Inner2DGenerator_New(  INNER2DGENERATOR_DEFARGS  );
	void _Inner2DGenerator_Init( Inner2DGenerator* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _Inner2DGenerator_Delete( void* generator );
	void _Inner2DGenerator_Print( void* generator, Stream* stream );
	void _Inner2DGenerator_AssignFromXML( void* generator, Stg_ComponentFactory* cf, void* data );
	void _Inner2DGenerator_Build( void* generator, void* data );
	void _Inner2DGenerator_Initialise( void* generator, void* data );
	void _Inner2DGenerator_Execute( void* generator, void* data );
	void _Inner2DGenerator_Destroy( void* generator, void* data );

	void Inner2DGenerator_Generate( void* generator, void* _mesh );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void Inner2DGenerator_SetElementMesh( void* generator, void* mesh );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	void Inner2DGenerator_BuildTopology( Inner2DGenerator* self, FeMesh* mesh );
	void Inner2DGenerator_BuildGeometry( Inner2DGenerator* self, FeMesh* mesh );
	void Inner2DGenerator_BuildElementTypes( Inner2DGenerator* self, FeMesh* mesh );

#endif /* __StgFEM_Discretisaton_Inner2DGenerator_h__ */

