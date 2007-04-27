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
** $Id: MeshVariable.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include "Base/Base.h"
#include "Discretisation/Geometry/Geometry.h"

#include "Mesh.h"


/* Textual name of this class */
const Type MeshVariable_Type = "MeshVariable";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

MeshVariable* MeshVariable_New( Name name ) {
	return _MeshVariable_New( sizeof(MeshVariable), 
				  MeshVariable_Type, 
				  _MeshVariable_Delete, 
				  _MeshVariable_Print, 
				  NULL, 
				  (void* (*)(Name))_MeshVariable_New, 
				  _MeshVariable_Construct, 
				  _MeshVariable_Build, 
				  _MeshVariable_Initialise, 
				  _MeshVariable_Execute, 
				  _MeshVariable_Destroy, 
				  name, 
				  False, 
				  0, 
				  NULL, 
				  NULL, 
				  NULL, 
				  NULL, 
				  NULL, 
				  NULL, 
				  NULL, 
				  NULL );
}

MeshVariable* _MeshVariable_New( MESHVARIABLE_DEFARGS ) {
	MeshVariable* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(MeshVariable) );
	self = (MeshVariable*)_Variable_New( VARIABLE_PASSARGS );

	/* Virtual info */

	/* MeshVariable info */
	_MeshVariable_Init( self );

	return self;
}

void _MeshVariable_Init( MeshVariable* self ) {
	self->mesh = NULL;
	self->topoDim = MT_VERTEX;
	self->meshArraySize = 0;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _MeshVariable_Delete( void* meshVariable ) {
	MeshVariable*	self = (MeshVariable*)meshVariable;

	MeshVariable_Destruct( self );

	/* Delete the parent. */
	_Variable_Delete( self );
}

void _MeshVariable_Print( void* meshVariable, Stream* stream ) {
	MeshVariable*	self = (MeshVariable*)meshVariable;
	
	/* Set the Journal for printing informations */
	Stream* meshVariableStream;
	meshVariableStream = Journal_Register( InfoStream_Type, "MeshVariableStream" );

	/* Print parent */
	Journal_Printf( stream, "MeshVariable (ptr): (%p)\n", self );
	_Variable_Print( self, stream );
}

void _MeshVariable_Construct( void* meshVariable, Stg_ComponentFactory* cf, void* data ) {
	MeshVariable*		self = (MeshVariable*)meshVariable;
	SizeT			    dataOffsets[]     = { 0 };
	Variable_DataType	dataTypes[]       = { 0 };		/* Init value later */
	Index			    dataTypeCounts[]  = { 1 };
	Dictionary *        componentDict     = NULL;
	Dictionary *        thisComponentDict = NULL;
	Name                dataTypeName      = NULL;
	Name                rankName          = NULL;
	void *              variableRegister  = NULL;
	void *              pointerRegister   = NULL;
	Name*               names             = NULL;
	Stream*             error             = Journal_Register( Error_Type, self->type );
	Mesh*			mesh;
	
	assert( self );

	componentDict = cf->componentDict;
	assert( componentDict );
	thisComponentDict = Dictionary_GetDictionary( componentDict, self->name );
	assert( thisComponentDict );
	
	/* Grab Registers */
	variableRegister = Stg_ObjectList_Get( cf->registerRegister, "Variable_Register" );
	assert( variableRegister );
	pointerRegister = Stg_ObjectList_Get( cf->registerRegister, "Pointer_Register" );
	assert( pointerRegister );

	/* Construct the mesh. */
	mesh = Stg_ComponentFactory_ConstructByKey( cf, self->name, "mesh", Mesh, True, data );
	MeshVariable_SetMesh( self, mesh );

	/* Get the topological element we're intereseted in. */
	self->topoDim = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "topologicalDim", 0 );
			
	/* Get Type of Variable */
	dataTypeName = Dictionary_GetString( thisComponentDict, "DataType" );
	if ( !strcasecmp( dataTypeName, "Double" ) )
		dataTypes[0] = Variable_DataType_Double;
	else if ( !strcasecmp( dataTypeName, "Float" ) )
		dataTypes[0] = Variable_DataType_Float;
	else if ( !strcasecmp( dataTypeName, "Int" ) )
		dataTypes[0] = Variable_DataType_Int;
	else if ( !strcasecmp( dataTypeName, "Char" ) )
		dataTypes[0] = Variable_DataType_Char;
	else if ( !strcasecmp( dataTypeName, "Short" ) )
		dataTypes[0] = Variable_DataType_Short;
	else 
		Journal_Firewall( False, error, "Variable '%s' cannot understand data type '%s'\n", self->name, dataTypeName );

	/* Get Rank of Variable - i.e. Scalar or Vector */
	rankName = Dictionary_GetString( thisComponentDict, "Rank" );
	if( !strcasecmp( rankName, "Scalar" ) ){
		dataTypeCounts[0] = 1;
	}
	else if ( !strcasecmp( rankName, "Vector" ) ){
		Dictionary_Entry_Value* list;
		Index                   nameCount = 0;

		/* Get Names from list */
		if (( list = Dictionary_Get( thisComponentDict, "names" ) )) {
			Index entry_I;

			nameCount = Dictionary_Entry_Value_GetCount( list );
			names = Memory_Alloc_Array( Name, nameCount, "Variable Names" );

			for ( entry_I = 0 ; entry_I < nameCount ; entry_I++ )
				names[ entry_I ] = Dictionary_Entry_Value_AsString( Dictionary_Entry_Value_GetElement(list, entry_I ) );
		}
		dataTypeCounts[0] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "VectorComponentCount", nameCount );

		Journal_Firewall( nameCount >= dataTypeCounts[0], error,
				"Variable '%s' has too few names in list for %d vector components.\n", self->name, dataTypeCounts[0] );
	}
	else
		Journal_Firewall( False, error, "Variable '%s' cannot understand rank '%s'\n", self->name, rankName );

	_Variable_Init( (Variable*)self, 
			1, 
			dataOffsets, 
			dataTypes, 
			dataTypeCounts, 
			names, 
			0, 
			&self->meshArraySize, 
			(void**)&self->arrayPtr,
			True, 
			variableRegister );

	/* Clean Up */
	if (names)
		Memory_Free(names);

#if 0
	MeshVariable*		self = (MeshVariable*)meshVariable;
	Variable_DataType	dataType;
	unsigned		dataRank;
	Dictionary*		dict;
	char*			dataTypeName;
	char*			rankName;
	unsigned		nDataNames;
	char**			names;
	Stream*			error;
	Mesh*			mesh;

	assert( self );
	assert( cf );

	/* Register streams. */
	error = Journal_Register( Error_Type, self->type );

	/* Shortcuts. */
	dict = Dictionary_GetDictionary( cf->componentDict, self->name );

	/* Construct the mesh. */
	mesh = Stg_ComponentFactory_ConstructByKey( cf, self->name, "mesh", Mesh, True );
	MeshVariable_SetMesh( self, mesh );

	/* Get the topological element we're intereseted in. */
	self->topoDim = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "topologicalDim", 0 );

	/* Get Type of Variable */
	dataTypeName = Stg_ComponentFactory_GetString( cf, self->name, "DataType", "" );
	if( !strcasecmp( dataTypeName, "Double" ) )
		dataType = Variable_DataType_Double;
	else if( !strcasecmp( dataTypeName, "Float" ) )
		dataType = Variable_DataType_Float;
	else if( !strcasecmp( dataTypeName, "Int" ) )
		dataType = Variable_DataType_Int;
	else if( !strcasecmp( dataTypeName, "Char" ) )
		dataType = Variable_DataType_Char;
	else if( !strcasecmp( dataTypeName, "Short" ) )
		dataType = Variable_DataType_Short;
	else {
		Journal_Firewall( False, error, "Variable '%s' cannot understand data type '%s'\n", 
				  self->name, dataTypeName );
	}

	/* Get Rank of Variable - i.e. Scalar or Vector */
	rankName = Stg_ComponentFactory_GetString( cf, self->name, "Rank", "" );
	if( !strcasecmp( rankName, "Scalar" ) ){
		dataRank = 1;
	}
	else if( !strcasecmp( rankName, "Vector" ) ) {
		Dictionary_Entry_Value* list;

		/* Get Names from list */
		if( (list = Dictionary_Get( dict, "names" ) )) {
			unsigned	n_i;

			nDataNames = Dictionary_Entry_Value_GetCount( list );
			dataNames = Memory_Alloc_Array_Unnamed( char*, nDataNames );

			for ( n_i = 0 ; n_i < nDataNames; n_i++ ) {
				Dictionary_Entry_Value*	tmp;

				tmp = Dictionary_Entry_Value_GetElement( list, n_i );
				dataNames[n_i] = Dictionary_Entry_Value_AsString( tmp );
			}
		}

		dataRank = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "VectorComponentCount", 
								  nNames );
		Journal_Firewall( nNames >= dataRank, error, 
				  "Variable '%s' has too few names in list for %d vector components.\n", 
				  self->name, dataRank );
	}
	else {
		Journal_Firewall( False, error, "Variable '%s' cannot understand rank '%s'\n", 
				  self->name, rankName );
	}

	/* Set the data type. */
	Variable_SetDataType( self, dataType, dataRank, dataNames );

	/* Free name array. */
	FreeArray( dataNames );
#endif
}

void _MeshVariable_Build( void* meshVariable, void* data ) {
	MeshVariable*	self = (MeshVariable*)meshVariable;

	assert( self );

	Build( self->mesh, data, False );

	self->meshArraySize = Mesh_GetDomainSize( self->mesh, self->topoDim );
	_Variable_Build( self, data );
}

void _MeshVariable_Initialise( void* meshVariable, void* data ) {
}

void _MeshVariable_Execute( void* meshVariable, void* data ) {
}

void _MeshVariable_Destroy( void* meshVariable, void* data ) {
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void MeshVariable_SetMesh( void* meshVariable, void* _mesh ) {
	MeshVariable*	self = (MeshVariable*)meshVariable;
	Mesh*		mesh = (Mesh*)_mesh;

	assert( self );

	MeshVariable_Destruct( self );

	self->mesh = mesh;
	if( mesh )
		List_Append( mesh->vars, self );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void MeshVariable_Destruct( MeshVariable* self ) {
}
