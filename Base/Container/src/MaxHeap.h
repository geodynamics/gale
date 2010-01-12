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
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/** \file
** <b>Role:</b>
**	Binary Tree class for objects.
**
** <b>Assumptions:</b>
**	None
**
** <b>Comments:</b>
**	None
**
** $Id: MaxHeap.h 2087 2005-2-25 02:28:44Z RaquibulHassan $
**
**/

#ifndef __StGermain_Base_Container_MaxHeap_h__
#define __StGermain_Base_Container_MaxHeap_h__

	/** Textual name for List class. */
	extern const Type MaxHeap_Type;
	
	/** \def __List See __List */
	#define __MaxHeap \
		/* General info */ \
		__Heap
	struct MaxHeap { __MaxHeap };
	
	MaxHeap* MaxHeap_New( void **keys, int keyLength, int numArrayElements, 
						Heap_KeySwapFunction *keySwap,
						Heap_CompareFunction *keyCompare,
						Heap_ExtendArrayFunc *extendArray );

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define MAXHEAP_DEFARGS \
                HEAP_DEFARGS

	#define MAXHEAP_PASSARGS \
                HEAP_PASSARGS

	MaxHeap* _MaxHeap_New(  MAXHEAP_DEFARGS  );

	
	void _MaxHeap_Init( MaxHeap* self );
	
	void MaxHeap_Heapify( _Heap *heap, int index );

	void* MaxHeap_Extract( void *heap );

	/** Stg_Class_Delete interface. */
		/** Stg_Class delete function */
	void _MaxHeap_DeleteFunc ( void *heap );
	
	/** Print interface. */
		/** Stg_Class print function */
	void _MaxHeap_PrintFunc ( void *heap, Stream* stream );

	/** Public functions */
	
	void MaxHeap_InsertFunc( _Heap *maxHeap, void *key );
	
	#define MaxHeap_Insert( maxHeap, key ) \
		MaxHeap_InsertFunc( (_Heap*)maxHeap, key );

	/** Private Functions */
	
#endif /* __StGermain_Base_Container_MaxHeap_h__ */


