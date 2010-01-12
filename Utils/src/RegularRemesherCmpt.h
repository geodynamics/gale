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
 ** $Id: RegularRemesherCmpt.h 2225 1970-01-02 13:48:23Z LukeHodkinson $
 **
 **~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgDomain_Utils_RegularRemesherCmpt_h__
#define __StgDomain_Utils_RegularRemesherCmpt_h__

#ifdef HAVE_PETSC

	/* Textual name of this class. */
	extern const Type RegularRemesherCmpt_Type;

	/* Virtual function types. */

	/* Class contents. */
	#define __RegularRemesherCmpt \
		/* General info */ \
		__Remesher \
		\
		/* Virtual info */ \
		\
		/* RegularRemesherCmpt info ... */ \
		RegularRemesher*	regRemesh;

	struct RegularRemesherCmpt { __RegularRemesherCmpt };



	/*-----------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	/* Create a RegularRemesherCmpt */
	RegularRemesherCmpt* RegularRemesherCmpt_New( Name name, AbstractContext* context, Mesh* mesh, RegularRemesher* regRemesh );

	/* Creation implementation */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define REGULARREMESHERCMPT_DEFARGS \
                REMESHER_DEFARGS

	#define REGULARREMESHERCMPT_PASSARGS \
                REMESHER_PASSARGS

	RegularRemesherCmpt* _RegularRemesherCmpt_New(  REGULARREMESHERCMPT_DEFARGS  );

	/* Initialisation implementation functions */
	void _RegularRemesherCmpt_Init( void* remesher, RegularRemesher* regRemesh );


	/*-----------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _RegularRemesherCmpt_Delete( void* remesher );

	void _RegularRemesherCmpt_Print( void* remesher, Stream* stream );

	RegularRemesherCmpt* _RegularRemesherCmpt_DefaultNew( Name name );

	void _RegularRemesherCmpt_AssignFromXML( void* remesher, Stg_ComponentFactory* cf, void* data );

	void _RegularRemesherCmpt_Build( void* remesher, void* data );

	void _RegularRemesherCmpt_Initialise( void* remesher, void* data );

	void _RegularRemesherCmpt_Execute( void* remesher, void* data );

	void _RegularRemesherCmpt_Destroy( void* remesher, void* data );


/*-----------------------------------------------------------------------------------------------------------------------------
** Public functions
*/


/*-----------------------------------------------------------------------------------------------------------------------------
** Private Member functions
*/

#endif

#endif

