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
** $Id: P1.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_Discretisaton_P1_h__
#define __StgFEM_Discretisaton_P1_h__

	/** Textual name of this class */
	extern const Type P1_Type;

	/** Virtual function types */

	/** P1 class contents */
	#define __P1			\
		/* General info */	\
		__ElementType		\
					\
		/* Virtual info */	\
					\
		/* P1 info */

	struct P1 { __P1 };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/



	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define P1_DEFARGS \
                ELEMENTTYPE_DEFARGS

	#define P1_PASSARGS \
                ELEMENTTYPE_PASSARGS

	P1* P1_New( Name name );
	P1* _P1_New(  P1_DEFARGS  );
	void _P1_Init( P1* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _P1_Delete( void* elementType );
	void _P1_Print( void* elementType, Stream* stream );
	void _P1_AssignFromXML( void* elementType, Stg_ComponentFactory* cf, void* data );
	void _P1_Build( void* elementType, void* data );
	void _P1_Initialise( void* elementType, void* data );
	void _P1_Execute( void* elementType, void* data );
	void _P1_Destroy( void* elementType, void* data );

	void P1_EvalBasis( void* elementType, const double* localCoord, double* derivs );
	void P1_EvalLocalDerivs( void* elementType, const double* localCoord, double** derivs );
	int _P1_SurfaceNormal( void* elementType, unsigned element_I, unsigned dim, double* xi, double* normal );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

#endif /* __StgFEM_Discretisaton_P1_h__ */

