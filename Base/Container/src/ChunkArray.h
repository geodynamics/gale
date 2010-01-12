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
** $Id: ChunkArray.h 2087 2004-11-15 02:28:44Z RaquibulHassan $
**
**/


#ifndef __StGermain_Base_Container_ChunkArray_h__
#define __StGermain_Base_Container_ChunkArray_h__

#define CHUNK_ARRAY_DELTA 10
#define INVALID -1
#define SUCCESS 1
#define FAILURE 0
#define TWO_EXP16 65535

extern const Type ChunkArray_Type;

    struct Chunk{
        char *memory;
        int numFree;
        int chunkId;
        char **freeList;
    };

    typedef void ( ChunkArray_ResizeCallbackFunc ) ( void * );
    
	#define __ChunkArray \
		__Stg_Class \
        	SizeT                  elementSize; \
        	int                    numElementsPerChunk; \
       		int                    numChunks; \
        	struct Chunk           *chunks; \
        	int                    maxChunkEntries; \
        	int                    chunkToUse;
    
	struct ChunkArray {__ChunkArray};
    

ChunkArray*
ChunkArray_NewFunc
(
    int         elementSize,
    int         numElementsPerChunk
);

ChunkArray *_ChunkArray_New(
		SizeT					_sizeOfSelf,
		Type					type,
		Stg_Class_DeleteFunction*	_delete,
		Stg_Class_PrintFunction*	_print,
		Stg_Class_CopyFunction*		_copy
		);
	
void _ChunkArray_Init( ChunkArray* self );

#define ChunkArray_New(type, numElementsPerChunk) \
    ChunkArray_NewFunc(sizeof(type), numElementsPerChunk);

int
ChunkArray_GetChunkWithFreeSlots
(
    ChunkArray      *chunkArray
);

int
ChunkArray_GetFreeChunkSlot
(
    ChunkArray      *chunkArray
);

void
_ChunkArray_Delete
(
    void    *chunkArray
);

void _ChunkArray_Print( void *self, Stream *myStream );

int
ChunkArray_CreateChunk
(
    ChunkArray      *chunkArray,
    int         pos
);

void *
ChunkArray_NewObjectFunc
(
    SizeT           elementSize,
    ChunkArray      *chunkArray
);

unsigned int
ChunkArray_NewObjectIDFunc
(
    SizeT           elementSize,
    ChunkArray      *chunkArray
);
        
int
ChunkArray_DeleteObject
(
       ChunkArray       *chunkArray,
       void             *object
);

int
ChunkArray_DeleteObjectID
(
       ChunkArray       *chunkArray,
       unsigned int          objectId
);

void
ChunkArray_Shrink
(
    ChunkArray      *chunkArray
);

    /** Public functions */
#define ChunkArray_NewObject( type, chunkArray ) \
    (type*)ChunkArray_NewObjectFunc( sizeof(type), chunkArray )

#define ChunkArray_NewObjectID( type, chunkArray ) \
    ChunkArray_NewObjectIDFunc( sizeof(type), chunkArray )

/*#define ChunkArray_ObjectAt(chunkArray, objectId) \
        (char*)(chunkArray->chunks[objectId >> 16].freeList[objectId & TWO_EXP16])*/

char* ChunkArray_ObjectAt(ChunkArray *chunkArray, unsigned int objectId);

#endif
