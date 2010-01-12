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
** $Id: CompressionAdaptor.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgDomain_Mesh_CompressionAdaptor_h__
#define __StgDomain_Mesh_CompressionAdaptor_h__

	/** Textual name of this class */
	extern const Type CompressionAdaptor_Type;

	/** Virtual function types */
	#define __CompressionAdaptor			\
		/* General info */			\
		__MeshAdaptor				\
							\
		/* Virtual info */			\
							\
		/* CompressionAdaptor info */		\
		double			compressionfactor;

	struct CompressionAdaptor { __CompressionAdaptor };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/



	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define COMPRESSIONADAPTOR_DEFARGS \
                MESHADAPTOR_DEFARGS

	#define COMPRESSIONADAPTOR_PASSARGS \
                MESHADAPTOR_PASSARGS

	CompressionAdaptor* CompressionAdaptor_New( Name name, AbstractContext* context );
	CompressionAdaptor* _CompressionAdaptor_New(  COMPRESSIONADAPTOR_DEFARGS  );
	void _CompressionAdaptor_Init( CompressionAdaptor* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _CompressionAdaptor_Delete( void* adaptor );
	void _CompressionAdaptor_Print( void* adaptor, Stream* stream );
	void _CompressionAdaptor_AssignFromXML( void* adaptor, Stg_ComponentFactory* cf, void* data );
	void _CompressionAdaptor_Build( void* adaptor, void* data );
	void _CompressionAdaptor_Initialise( void* adaptor, void* data );
	void _CompressionAdaptor_Execute( void* adaptor, void* data );
	void _CompressionAdaptor_Destroy( void* adaptor, void* data );

	void CompressionAdaptor_Generate( void* adaptor, void* _mesh, void* data );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

#endif /* __StgDomain_Mesh_CompressionAdaptor_h__ */

