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
** $Id: Hasher.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Base_Container_Hasher_h__
#define __Base_Contianer_Hasher_h__

	/** Textual name of this class */
	extern const Type Hasher_Type;

	/** Virtual function types */
	typedef unsigned (Hasher_HashFunc)( void* hasher, void* key );

	/** Class contents */
	#define __Hasher				\
		/* General info */			\
		__Stg_Class				\
							\
		/* Virtual info */			\
		Hasher_HashFunc*	hashFunc;	\
							\
		/* Hasher info */			\
		unsigned		keySize;

	struct Hasher { __Hasher };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define HASHER_DEFARGS				\
		STG_CLASS_DEFARGS,			\
		Hasher_HashFunc*	hashFunc

	#define HASHER_PASSARGS			\
		STG_CLASS_PASSARGS, hashFunc

	Hasher* Hasher_New();
	Hasher* _Hasher_New( HASHER_DEFARGS );
	void _Hasher_Init( Hasher* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _Hasher_Delete( void* hasher );
	void _Hasher_Print( void* hasher, Stream* stream );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void Hasher_SetKeySize( void* hasher, unsigned keySize );

	#define Hasher_Hash( hasher, key )				\
		(assert( (hasher) && ((Hasher*)hasher)->hashFunc ),	\
		 ((Hasher*)hasher)->hashFunc( hasher, key ))

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	void Hasher_Destruct( Hasher* self );

#endif /* __Base_Container_Hasher_h__ */
