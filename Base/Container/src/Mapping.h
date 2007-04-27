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
** $Id: Mapping.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Base_Container_Mapping_h__
#define __Base_Contianer_Mapping_h__

	/** Textual name of this class */
	extern const Type		Mapping_Type;
	extern const unsigned long	Mapping_VacantIndex;
	extern const unsigned long	Mapping_CollidedIndex;

	/** Virtual function types */

	/** Class contents */
	typedef struct {
		void*	key;
		void*	value;
	} Mapping_Tuple;

	#define __Mapping				\
		/* General info */			\
		__Stg_Class				\
							\
		/* Virtual info */			\
							\
		/* Mapping info */			\
		unsigned		maxItems;	\
		unsigned		delta;		\
		unsigned		nItems;		\
		unsigned		keySize;	\
		unsigned		valSize;	\
		Stg_Byte*		keys;		\
		Stg_Byte*		vals;		\
							\
		Hasher*			hasher;		\
		unsigned long*		itemInds;

	struct Mapping { __Mapping };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define MAPPING_DEFARGS		\
		STG_CLASS_DEFARGS

	#define MAPPING_PASSARGS	\
		STG_CLASS_PASSARGS

	Mapping* Mapping_New();
	Mapping* _Mapping_New( MAPPING_DEFARGS );
	void _Mapping_Init( Mapping* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _Mapping_Delete( void* mapping );
	void _Mapping_Print( void* mapping, Stream* stream );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void Mapping_SetTupleSizes( void* mapping, unsigned keySize, unsigned valueSize );
	void Mapping_SetMaxSize( void* mapping, unsigned maxSize );
	void Mapping_SetDelta( void* mapping, unsigned delta );

	void Mapping_Insert( void* mapping, void* key, void* value );
	void Mapping_Remove( void* mapping, void* key );
	void Mapping_Clear( void* mapping );

	unsigned Mapping_GetSize( void* mapping );
	Bool Mapping_Map( void* mapping, void* key, void** value );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	void Mapping_InitHasher( Mapping* self );
	void Mapping_InitItems( Mapping* self );
	void Mapping_Collision( Mapping* self, void* key, void* val, unsigned hashInd );
	void Mapping_DestructData( Mapping* self );
	void Mapping_Destruct( Mapping* self );

#endif /* __Base_Container_Mapping_h__ */
