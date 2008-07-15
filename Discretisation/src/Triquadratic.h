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
** $Id: Biquadratic.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Discretisaton_Mesh_Triquadratic_h__
#define __Discretisaton_Mesh_Triquadratic_h__

	/** Textual name of this class */
	extern const Type Triquadratic_Type;

	/** Virtual function types */

	/** Triquadratic class contents */
	#define __Triquadratic		\
		/* General info */	\
		__ElementType		\
					\
		/* Virtual info */	\
					\
		/* Triquadratic info */

	struct Triquadratic { __Triquadratic };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define TRIQUADRATIC_DEFARGS \
		ELEMENTTYPE_DEFARGS

	#define TRIQUADRATIC_PASSARGS \
		ELEMENTTYPE_PASSARGS

	Triquadratic* Triquadratic_New( Name name );
	Triquadratic* _Triquadratic_New( TRIQUADRATIC_DEFARGS );
	void _Triquadratic_Init( Triquadratic* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _Triquadratic_Delete( void* elementType );
	void _Triquadratic_Print( void* elementType, Stream* stream );
	void _Triquadratic_Construct( void* elementType, Stg_ComponentFactory* cf, void* data );
	void _Triquadratic_Build( void* elementType, void* data );
	void _Triquadratic_Initialise( void* elementType, void* data );
	void _Triquadratic_Execute( void* elementType, void* data );
	void _Triquadratic_Destroy( void* elementType, void* data );

	void Triquadratic_EvalBasis( void* elementType, const double* localCoord, double* derivs );
	void Triquadratic_EvalLocalDerivs( void* elementType, const double* localCoord, double** derivs );
	double Triquadratic_JacobianDeterminantSurface( void* elementType, void * mesh, const double localCoord[],
		       					unsigned* nodes, unsigned norm );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

#endif /* __Discretisaton_Mesh_Triquadratic_h__ */
