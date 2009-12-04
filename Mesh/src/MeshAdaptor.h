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
** $Id: MeshAdaptor.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgDomain_Mesh_MeshAdaptor_h__
#define __StgDomain_Mesh_MeshAdaptor_h__

	/** Textual name of this class */
	extern const Type MeshAdaptor_Type;

	/** Virtual function types */

	/** MeshAdaptor class contents */
	#define __MeshAdaptor				\
		/* General info */			\
		__MeshGenerator				\
							\
		/* Virtual info */			\
							\
		/* MeshAdaptor info */			\
		MeshGenerator*		generator;	\
		Mesh*			srcMesh;

	struct MeshAdaptor { __MeshAdaptor };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/



	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define MESHADAPTOR_DEFARGS \
                MESHGENERATOR_DEFARGS

	#define MESHADAPTOR_PASSARGS \
                MESHGENERATOR_PASSARGS

	MeshAdaptor* _MeshAdaptor_New(  MESHADAPTOR_DEFARGS  );
	void _MeshAdaptor_Init( MeshAdaptor* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _MeshAdaptor_Delete( void* adaptor );
	void _MeshAdaptor_Print( void* adaptor, Stream* stream );
	void _MeshAdaptor_AssignFromXML( void* adaptor, Stg_ComponentFactory* cf, void* data );
	void _MeshAdaptor_Build( void* adaptor, void* data );
	void _MeshAdaptor_Initialise( void* adaptor, void* data );
	void _MeshAdaptor_Execute( void* adaptor, void* data );
	void _MeshAdaptor_Destroy( void* adaptor, void* data );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void MeshAdaptor_SetGenerator( void* adaptor, void* generator );
	void MeshAdaptor_SetSourceMesh( void* adaptor, void* mesh );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

#endif /* __StgDomain_Mesh_MeshAdaptor_h__ */

