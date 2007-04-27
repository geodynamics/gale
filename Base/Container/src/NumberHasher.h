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
** $Id: NumberHasher.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Base_Container_NumberHasher_h__
#define __Base_Contianer_NumberHasher_h__

	/** Textual name of this class */
	extern const Type NumberHasher_Type;

	/** Virtual function types */

	/** Class contents */
	#define __NumberHasher			\
		/* General info */		\
		__Hasher			\
						\
		/* Virtual info */		\
						\
		/* NumberHasher info */

	struct NumberHasher { __NumberHasher };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define NUMBERHASHER_DEFARGS	\
		HASHER_DEFARGS

	#define NUMBERHASHER_PASSARGS	\
		HASHER_PASSARGS

	NumberHasher* NumberHasher_New();
	NumberHasher* _NumberHasher_New( NUMBERHASHER_DEFARGS );
	void _NumberHasher_Init( NumberHasher* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _NumberHasher_Delete( void* hasher );
	void _NumberHasher_Print( void* hasher, Stream* stream );

	unsigned _NumberHasher_Hash( void* hasher, void* key );
	void _NumberHasher_CalcTableSize( void* hasher );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

#endif /* __Base_Container_NumberHasher_h__ */
