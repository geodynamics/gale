/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd,
** 110 Victoria Street, Melbourne, 3053, Australia.
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
**  Role:
**
** Assumptions:
**
** Invariants:
**
** Comments:
**
** $Id: Remesher.h 2225 1970-01-02 13:48:23Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StGermain_Domain_Utils_Remesher_h__
#define __StGermain_Domain_Utils_Remesher_h__

	/* Textual name of this class. */
	extern const Type Remesher_Type;

	/* Virtual function types. */
	typedef void (Remesher_RemeshFunc)( void* _self );

	/* Class contents. */
	#define __Remesher \
		__Stg_Component \
		AbstractContext*		context; \
		Remesher_RemeshFunc*	remeshFunc; \
		Mesh*						mesh;

	struct Remesher { __Remesher };



	/* Constructors */

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define REMESHER_DEFARGS \
                STG_COMPONENT_DEFARGS, \
                Remesher_RemeshFunc*  remeshFunc

	#define REMESHER_PASSARGS \
                STG_COMPONENT_PASSARGS, \
	        remeshFunc

	Remesher* _Remesher_New(  REMESHER_DEFARGS  );

	void _Remesher_Init( void* remeshser, AbstractContext* context, Mesh* mesh );

	/* Virtual functions */

	void _Remesher_Delete( void* remesher );

	void _Remesher_Print( void* remesher, Stream* stream );

	Remesher* _Remesher_DefaultNew( Name name );

	void _Remesher_AssignFromXML( void* remesher, Stg_ComponentFactory* cf, void* data );

	void _Remesher_Build( void* remesher, void* data );

	void _Remesher_Initialise( void* remesher, void* data );

	void _Remesher_Execute( void* remesher, void* data );

	void _Remesher_Destroy( void* remesher, void* data );

	/* Public functions */

	#define Remesher_Remesh( self ) \
		(self)->remeshFunc( self )

#endif

