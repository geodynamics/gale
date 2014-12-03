/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd,
** 110 Victoria Street, Melbourne, 3053, Australia.
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
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/** \file
 ** <b>Role:</b>
 **	Abstract class faciliting how class inheritance is done.
 **
 ** <b>Assumptions:</b>
 **	None
 **
 ** <b>Comments:</b>
 **	None
 **
 ** $Id: Class.h 3904 2006-12-14 00:52:06Z LukeHodkinson $
 **
 **/

#ifndef __StGermain_Base_Foundation_NewClass_h__
#define __StGermain_Base_Foundation_NewClass_h__

/* Forward declaring Stream and PtrMap. */
struct Stream;
struct PtrMap;

/* Some macros to do memory recording. */
#define Class_Alloc( self, itmType )			\
   (assert( self ),					\
    ((NewClass*)self)->curAllocd += sizeof(itmType),	\
    MemAlloc( itmType, ((NewClass*)self)->type ))

#define Class_Array( self, itmType, nItms )			\
   (assert( self ),						\
    ((NewClass*)self)->curAllocd += sizeof(itmType) * (nItms),	\
    MemArray( itmType, nItms, ((NewClass*)self)->type ))

#define Class_Array2D( self, itmType, nItmsX, nItmsY )			\
   (assert( self ),							\
    ((NewClass*)self)->curAllocd += sizeof(itmType) * 			\
    (nItmsX) * (nItmsY) + sizeof(itmType*) * (nItmsX),			\
    MemArray2D( itmType, nItmsX, nItmsY, ((NewClass*)self)->type ))

#define Class_Rearray( self, ptr, itmType, nItms )			\
   (assert( self ),							\
    (ptr) ? ((NewClass*)self)->curAllocd -= Memory_SizeGet( ptr ) : 0,	\
    ((NewClass*)self)->curAllocd += sizeof(itmType) * (nItms),		\
    MemRearray( ptr, itmType, nItms, ((NewClass*)self)->type ))

#define Class_Free( self, ptr )						\
   (assert( self ),							\
    (ptr) ? ((NewClass*)self)->curAllocd -= Memory_SizeGet( ptr ) : 0,	\
    MemFree( ptr ))

/* Macros to handle type checked type casting. */
#define Class_IsSuper( ptr, typeName )				\
   (assert( ptr ), assert( ((NewClass*)ptr)->isSuperFunc ),	\
    ((NewClass*)ptr)->isSuperFunc( typeName##_Type ))

#define Class_Cast( ptr, typeName )		\
   (assert( Class_IsSuper( ptr, typeName ) ),	\
    (typeName*)ptr )

#define Class_ConstCast( ptr, typeName )	\
   (assert( Class_IsSuper( ptr, typeName ) ),	\
    (const typeName*)ptr )

#include "StGermain/Base/Foundation/ClassClear.h"
#define CLASSDIR StGermain/Base/Foundation
#define CLASSNAME NewClass
#include "StGermain/Base/Foundation/ClassHdr.h"

#endif /* __StGermain_Base_Foundation_NewClass_h__ */
