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
** $Id: RangeSet.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Base_Container_RangeSet_h__
#define __Base_Container_RangeSet_h__

	/** Textual name of this class */
	extern const Type RangeSet_Type;

	/** Virtual function types */

	/** Mesh class contents */
	#define __RangeSet				\
		/* General info */			\
		__Stg_Class				\
							\
		/* Virtual info */			\
							\
		/* RangeSet info */			\
		unsigned		nInds;		\
		BTree*			btree;

	struct RangeSet { __RangeSet };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define RANGESET_DEFARGS \
		STG_CLASS_DEFARGS

	#define RANGESET_PASSARGS \
		STG_CLASS_PASSARGS

	RangeSet* RangeSet_New();
	RangeSet* _RangeSet_New( RANGESET_DEFARGS );
	void _RangeSet_Init( RangeSet* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _RangeSet_Delete( void* rangeSet );
	void _RangeSet_Print( void* rangeSet, Stream* stream );

	#define RangeSet_Copy( self )					\
		(RangeSet*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define RangeSet_DeepCopy( self )				\
		(RangeSet*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _RangeSet_Copy( void* rangeSet, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void RangeSet_SetIndices( void* rangeSet, unsigned nInds, unsigned* inds );
	void RangeSet_AddIndices( void* rangeSet, unsigned nInds, unsigned* inds );
	void RangeSet_SetRange( void* rangeSet, unsigned begin, unsigned end, unsigned step );
	void RangeSet_Clear( void* rangeSet );

	void RangeSet_GetIndices( void* rangeSet, unsigned* nInds, unsigned** inds );
	unsigned RangeSet_GetSize( void* rangeSet );
	unsigned RangeSet_GetNumRanges( void* rangeSet );
	RangeSet_Range* RangeSet_GetRange( void* rangeSet, unsigned index );
	Bool RangeSet_HasIndex( void* rangeSet, unsigned index );

	void RangeSet_Union( void* rangeSet, RangeSet* rSet );
	void RangeSet_Intersection( void* rangeSet, RangeSet* rSet );
	void RangeSet_Subtraction( void* rangeSet, RangeSet* rSet );

	void RangeSet_Pickle( void* rangeSet, unsigned* nBytes, Stg_Byte** bytes );
	void RangeSet_Unpickle( void* rangeSet, unsigned nBytes, Stg_Byte* bytes );

	void RangeSet_Range_Intersection( RangeSet_Range* left, RangeSet_Range* right, RangeSet_Range* result );
	Bool RangeSet_Range_HasIndex( RangeSet_Range* self, unsigned index );
	unsigned RangeSet_Range_GetNumIndices( RangeSet_Range* self );
	void RangeSet_Range_GetIndices( RangeSet_Range* self, unsigned* nInds, unsigned** inds );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	void RangeSet_GetIndicesParse( void* data, void* _parseStruct );
	void RangeSet_GetRangeParse( void* data, void* _parseStruct );
	void RangeSet_IntersectionParse( void* data, void* _parse );
	void RangeSet_SubtractionParse( void* data, void* _parse );
	void RangeSet_PickleParse( void* data, void* _parse );
#if 0
	void RangeSet_RangeIntersectionParse( BTreeNode* node, RangeSet_ParseStruct* parse );
#endif
	int RangeSet_SortCmp( const void* itema, const void* itemb );
	int RangeSet_DataCompare( void* left, void* right );
	void RangeSet_DataCopy( void** dstData, void* data, SizeT size );
	void RangeSet_DataDelete( void* data );

	void RangeSet_Destruct( RangeSet* self );

#endif /* __Base_Container_RangeSet_h__ */
