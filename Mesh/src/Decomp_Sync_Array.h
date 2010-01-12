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
** $Id: Decomp_Sync_Array.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgDomain_Mesh_Decomp_Sync_Array_h__
#define __StgDomain_Mesh_Decomp_Sync_Array_h__

	/** Textual name of this class */
	extern const Type Decomp_Sync_Array_Type;

	/** Virtual function types */

	/** Mesh class contents */
	#define __Decomp_Sync_Array		\
		/* General info */		\
		__Stg_Class			\
						\
		/* Virtual info */		\
						\
		/* Decomp_Sync_Array info */	\
		Decomp_Sync*	sync;		\
						\
		void*		snkArray;	\
		unsigned	snkStride;	\
		unsigned*	snkDisps;	\
		unsigned*	snkSizes;	\
		unsigned*	snkOffs;	\
						\
		void*		srcArray;	\
		unsigned	srcStride;	\
		unsigned*	srcDisps;	\
		unsigned*	srcSizes;	\
		unsigned*	srcOffs;	\
						\
		size_t		itemSize;

	struct Decomp_Sync_Array { __Decomp_Sync_Array };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define DECOMP_SYNC_ARRAY_DEFARGS	\
		STG_CLASS_DEFARGS

	#define DECOMP_SYNC_ARRAY_PASSARGS	\
		STG_CLASS_PASSARGS

	Decomp_Sync_Array* Decomp_Sync_Array_New();
	Decomp_Sync_Array* _Decomp_Sync_Array_New( DECOMP_SYNC_ARRAY_DEFARGS );
	void _Decomp_Sync_Array_Init( Decomp_Sync_Array* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _Decomp_Sync_Array_Delete( void* array );
	void _Decomp_Sync_Array_Print( void* array, Stream* stream );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void Decomp_Sync_Array_SetSync( void* array, Decomp_Sync* sync );
	void Decomp_Sync_Array_SetMemory( void* array, 
					  void* localArray, void* remoteArray, 
					  size_t localStride, size_t remoteStride, 
					  size_t itemSize );
	void Decomp_Sync_Array_Sync( void* array );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	void Decomp_Sync_Array_BuildArray( Decomp_Sync_Array* self );
	void Decomp_Sync_Array_Destruct( Decomp_Sync_Array* self );

#endif /* __StgDomain_Mesh_Decomp_Sync_Array_h__ */
