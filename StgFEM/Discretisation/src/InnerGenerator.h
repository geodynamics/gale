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
** $Id: InnerGenerator.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_Discretisaton_InnerGenerator_h__
#define __StgFEM_Discretisaton_InnerGenerator_h__

	/** Textual name of this class */
	extern const Type InnerGenerator_Type;

	/** Virtual function types */

	/** InnerGenerator class contents */
	#define __InnerGenerator		\
		/* General info */		\
		__MeshGenerator			\
						\
		/* Virtual info */		\
						\
		/* InnerGenerator info */		\
		Mesh*		elMesh;

	struct InnerGenerator { __InnerGenerator };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/



	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define INNERGENERATOR_DEFARGS \
                MESHGENERATOR_DEFARGS

	#define INNERGENERATOR_PASSARGS \
                MESHGENERATOR_PASSARGS

	InnerGenerator* InnerGenerator_New( Name name, AbstractContext* context );
	InnerGenerator* _InnerGenerator_New(  INNERGENERATOR_DEFARGS  );
	void _InnerGenerator_Init( InnerGenerator* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _InnerGenerator_Delete( void* generator );
	void _InnerGenerator_Print( void* generator, Stream* stream );
	void _InnerGenerator_AssignFromXML( void* generator, Stg_ComponentFactory* cf, void* data );
	void _InnerGenerator_Build( void* generator, void* data );
	void _InnerGenerator_Initialise( void* generator, void* data );
	void _InnerGenerator_Execute( void* generator, void* data );
	void _InnerGenerator_Destroy( void* generator, void* data );

	void InnerGenerator_Generate( void* generator, void* _mesh );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void InnerGenerator_SetElementMesh( void* generator, void* mesh );
        void InnerGenerator_SetCoordinates( InnerGenerator* self, FeMesh* mesh );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	void InnerGenerator_BuildTopology( InnerGenerator* self, FeMesh* mesh );
	void InnerGenerator_BuildGeometry( InnerGenerator* self, FeMesh* mesh );
	void InnerGenerator_BuildElementTypes( InnerGenerator* self, FeMesh* mesh );

#endif /* __StgFEM_Discretisaton_InnerGenerator_h__ */

