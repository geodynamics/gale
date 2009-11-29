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
** $Id: RegularTrilinear.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Discretisaton_Mesh_RegularTrilinear_h__
#define __Discretisaton_Mesh_RegularTrilinear_h__

	/** Textual name of this class */
	extern const Type RegularTrilinear_Type;

	/** Virtual function types */

	/** RegularTrilinear class contents */
	#define __RegularTrilinear		\
		/* General info */		\
		__TrilinearElementType		\
						\
		/* Virtual info */		\
						\
		/* RegularTrilinear info */

	struct RegularTrilinear { __RegularTrilinear };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/



	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define REGULARTRILINEAR_DEFARGS \
                TRILINEARELEMENTTYPE_DEFARGS

	#define REGULARTRILINEAR_PASSARGS \
                TRILINEARELEMENTTYPE_PASSARGS

	RegularTrilinear* RegularTrilinear_New( Name name );
	RegularTrilinear* _RegularTrilinear_New(  REGULARTRILINEAR_DEFARGS  );
	void _RegularTrilinear_Init( RegularTrilinear* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _RegularTrilinear_Delete( void* elementType );
	void _RegularTrilinear_Print( void* elementType, Stream* stream );
	void _RegularTrilinear_AssignFromXML( void* elementType, Stg_ComponentFactory* cf, void* data );
	void _RegularTrilinear_Build( void* elementType, void* data );
	void _RegularTrilinear_Initialise( void* elementType, void* data );
	void _RegularTrilinear_Execute( void* elementType, void* data );
	void _RegularTrilinear_Destroy( void* elementType, void* data );

	void RegularTrilinear_ConvertGlobalCoordToElLocal( void* elementType, void* mesh, unsigned element, 
							   const double* globalCoord, double* localCoord );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

#endif /* __Discretisaton_Mesh_RegularTrilinear_h__ */

