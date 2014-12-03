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
** $Id: Set.c 2192 2004-10-15 02:45:38Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include "Base/Foundation/Foundation.h"
#include "Base/IO/IO.h"

#include "types.h"
#include "BTreeNode.h"
#include "BTree.h"
#include "Set.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>


/* Textual name of this class */
const Type Set_Type = "Set";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

Set* Set_New_all(
		Dictionary*					dictionary, 
		SizeT						elementSize, 
		BTree_compareFunction*				compareFunc, 
		BTree_dataCopyFunction*				dataCopyFunc, 
		BTree_dataDeleteFunction*			dataDeleteFunc )
{
	/* Variables set in this function */
	SizeT                            _sizeOfSelf = sizeof(Set);
	Type                                    type = Set_Type;
	Stg_Class_DeleteFunction*            _delete = _Set_Delete;
	Stg_Class_PrintFunction*              _print = _Set_Print;
	Stg_Class_CopyFunction*                _copy = NULL;
	Set_UnionFunc*                    _unionFunc = _Set_Union;
	Set_IntersectionFunc*      _intersectionFunc = _Set_Intersection;
	Set_SubtractionFunc*        _subtractionFunc = _Set_Subtraction;

	return _Set_New(  SET_PASSARGS  );
}


void Set_Init(
		Set*						self,
		Dictionary*					dictionary, 
		SizeT						elementSize, 
		BTree_compareFunction*				compareFunc, 
		BTree_dataCopyFunction*				dataCopyFunc, 
		BTree_dataDeleteFunction*			dataDeleteFunc )
{
	/* General info */
	self->type = Set_Type;
	self->_sizeOfSelf = sizeof(Set);
	self->_deleteSelf = False;
	self->dictionary = dictionary;
	
	/* Virtual info */
	self->_delete = _Set_Delete;
	self->_print = _Set_Print;
	self->_copy = NULL;
	self->_unionFunc = _Set_Union;
	self->_intersectionFunc = _Set_Intersection;
	self->_subtractionFunc = _Set_Subtraction;
	_Stg_Class_Init( (Stg_Class*)self );
	
	/* Set info */
	_Set_Init( self, elementSize, compareFunc, dataCopyFunc, dataDeleteFunc );
}


Set* _Set_New(  SET_DEFARGS  )
{
	Set*	self;
	
	/* allocate memory */
	assert( _sizeOfSelf >= sizeof(Set) );
	self = (Set*)_Stg_Class_New(  STG_CLASS_PASSARGS  );
	
	/* general info */
	self->dictionary = dictionary;
	
	/* virtual info */
	self->_unionFunc = _unionFunc;
	self->_intersectionFunc = _intersectionFunc;
	self->_subtractionFunc = _subtractionFunc;
	
	/* Set info */
	_Set_Init( self, elementSize, compareFunc, dataCopyFunc, dataDeleteFunc );
	
	return self;
}


void _Set_Init(
		Set*						self, 
		SizeT						elementSize, 
		BTree_compareFunction*				compareFunc, 
		BTree_dataCopyFunction*				dataCopyFunc, 
		BTree_dataDeleteFunction*			dataDeleteFunc )
{
	/* TODO: convert to journal */
	assert( compareFunc && dataDeleteFunc );

	self->_elSize = elementSize;
	self->_btree = BTree_New( compareFunc, dataCopyFunc, dataDeleteFunc, NULL, BTREE_ALLOW_DUPLICATES );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Set_Delete( void* set ) {
	Set*	self = (Set*)set;
	
	/* delete the class itself */
	Stg_Class_Delete( self->_btree );
	
	/* delete parent */
	_Stg_Class_Delete( self );
}


void _Set_Print( void* set, Stream* stream ) {
	Set*		self = (Set*)set;
	Stream*		myStream = Journal_Register( InfoStream_Type, "SetStream" );

	/* print parent */
	_Stg_Class_Print( self, stream );
	
	/* general info */
	Journal_Printf( myStream, "Set (ptr): (%p)\n", self );
	Journal_Printf( myStream, "\tdictionary (ptr): %p\n", self->dictionary );
	
	/* virtual info */
	
	/* Set info */
}


void* _Set_Union( void* set, void* operand ) {
	Set*	self = (Set*)set;
	Set*	results;

	results = Set_New_all( self->dictionary, 
			       self->_elSize, 
			       self->_btree->compareFunction, 
			       self->_btree->dataCopyFunction, 
			       self->_btree->dataDeleteFunction );

	BTree_ParseTree( self->_btree, _Set_BTreeUnion, results );
	BTree_ParseTree( ((Set*)operand)->_btree, _Set_BTreeUnion, results );

	return results;
}


void* _Set_Intersection( void* set, void* operand ) {
	Set*	self = (Set*)set;
	Set*	results;
	void*	pack[2];

	results = Set_New_all( self->dictionary, 
			       self->_elSize, 
			       self->_btree->compareFunction, 
			       self->_btree->dataCopyFunction, 
			       self->_btree->dataDeleteFunction );

	pack[0] = operand;
	pack[1] = results;
	BTree_ParseTree( self->_btree, _Set_BTreeIntersection, pack );

	return results;
}


void* _Set_Subtraction( void* set, void* operand ) {
	Set*	self = (Set*)set;
	Set*	results;
	void*	pack[2];

	results = Set_New_all( self->dictionary, 
			       self->_elSize, 
			       self->_btree->compareFunction, 
			       self->_btree->dataCopyFunction, 
			       self->_btree->dataDeleteFunction );

	pack[0] = operand;
	pack[1] = results;
	BTree_ParseTree( self->_btree, _Set_BTreeSubtraction, pack );

	return results;
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

Bool Set_Insert( void* set, void* data ) {
	Set*      self = (Set*)set;

	if( BTree_FindNode( self->_btree, data ) != NULL ) {
		return False;
	}

	BTree_InsertNode( self->_btree, data, self->_elSize );

	return True;
}


void Set_Traverse( void* set, BTree_parseFunction* func, void* args ) {
	Set*      self = (Set*)set;

	BTree_ParseTree( self->_btree, func, args );
}


Bool Set_Exists( void* set, void* data ) {
	Set*      self = (Set*)set;

	return (BTree_FindNode( self->_btree, data ) != NULL) ? True : False;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void _Set_BTreeUnion( void* data, void* set ) {
	Set_Insert( (Set*)set, data );
}


void _Set_BTreeIntersection( void* data, void* pack ) {
	if( BTree_FindNode( ((Set**)pack)[0]->_btree, data ) != NULL ) {
		Set_Insert( ((Set**)pack)[1], data );
	}
}


void _Set_BTreeSubtraction( void* data, void* pack ) {
	if( BTree_FindNode( ((Set**)pack)[0]->_btree, data ) == NULL ) {
		Set_Insert( ((Set**)pack)[1], data );
	}
}


