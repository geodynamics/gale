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
** $Id: OperatorFunction.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdarg.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include <StgFEM/SLE/SLE.h>

#include "types.h"
#include "OperatorFunction.h"
#include "StiffnessMatrix.h"

#include <stdio.h>
#include <assert.h>


static void ConvertTypedStrValue( Type type, const char* inValue, OperatorFunction_Datum* outValue );


const Type OperatorFunction_Type = "OperatorFunction";


/*--------------------------------------------------------------------------------------------------------------------------
** Constructor
*/

OperatorFunction* _OperatorFunction_New( 
		SizeT						_sizeOfSelf,
		Type						type,
	        Stg_Class_DeleteFunction*			_delete,
		Stg_Class_PrintFunction*			_print,
		Stg_Class_CopyFunction*				_copy,
		Stg_Component_DefaultConstructorFunction*	_defaultConstructor,
		Stg_Component_ConstructFunction*		_construct,
		Stg_Component_BuildFunction*			_build,
		Stg_Component_InitialiseFunction*		_initialise,
		Stg_Component_ExecuteFunction*			_execute,
		Stg_Component_DestroyFunction*			_destroy,
		OperatorFunction_ApplyFunc*			_applyMatrix,
		OperatorFunction_ApplyFunc*			_applyNeumann,
		OperatorFunction_ApplyRHSFunc*			_applyRHS,
		Index						dataCount,
		Name*						dataName,
		Bool*						dataIsComponent,
		Type*						dataType,
		Bool*						dataIsRequired,
		void**						values,
		Name						name )
{
	OperatorFunction* 		self;

	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(OperatorFunction) );
	self = (OperatorFunction*)_Stg_Component_New( 
				_sizeOfSelf, 
				type,
				_delete,
				_print,
				_copy,
				_defaultConstructor,
				_construct,
				_build,
				_initialise,
				_execute,
				_destroy,
				name,
		       		NON_GLOBAL );

	_OperatorFunction_Init(
		self, 
		_applyMatrix,
	        _applyNeumann,
		_applyRHS,	
		dataCount,
		dataName,
		dataIsComponent,
		dataType,
		dataIsRequired,
		values );

		return self;
}

void _OperatorFunction_Init(
		void*				operatorFunction, 
		OperatorFunction_ApplyFunc*	applyMatrix, 
		OperatorFunction_ApplyFunc*	applyNeumann, 
		OperatorFunction_ApplyRHSFunc*	applyRHS, 
		Index				dataCount,
		Name*				dataName,
		Bool*				dataIsComponent,
		Type*				dataType,
		Bool*				dataIsRequired,
		void**				values )
{
	OperatorFunction* 	self 	= (OperatorFunction*)operatorFunction;
	
	self->applyMatrix = applyMatrix;
	self->applyNeumann = applyNeumann;
	self->applyRHS = applyRHS;
	self->_dataCount = dataCount;
	if( self->_dataCount ) {
		Index i;

		/* Copy names */
		self->_dataName = Memory_Alloc_Array( Name, self->_dataCount, "OperatorFunction->_dataName" );
		for( i = 0; i < self->_dataCount; i++ ) {
			self->_dataName[i] = StG_Strdup( dataName[i] );
		}
		
		/* Copy isComponent */
		self->_dataIsComponent = Memory_Alloc_Array( Bool, self->_dataCount, "OperatorFunction->_dataIsComponent" );
		memcpy( self->_dataIsComponent, dataIsComponent, self->_dataCount * sizeof( Bool ) );

		/* Copy dataType */
		self->_dataType = Memory_Alloc_Array( Type, self->_dataCount, "OperatorFunction->_dataType" );
		for( i = 0; i < self->_dataCount; i++ ) {
			self->_dataType[i] = StG_Strdup( dataType[i] );
		}

		/* Copy isRequired */
		self->_dataIsRequired = Memory_Alloc_Array( Bool, self->_dataCount, "OperatorFunction->_dataIsRequired" );
		memcpy( self->_dataIsRequired, dataIsRequired, self->_dataCount * sizeof( Bool ) );

		/* There is no way of telling how to destroy the "values" (void**), so we wont copy the items... 
		   we just assume ownership when it comes to deleting. */
	}
	else {
		self->_dataName = NULL;
		self->_dataIsComponent = NULL;
		self->_dataType = NULL;
		self->_dataIsRequired = NULL;
	}
	self->_values = values; /*See note above */
}

OperatorFunction* OperatorFunction_New(
		OperatorFunction_ApplyFunc*	applyMatrix, 
		OperatorFunction_ApplyFunc*	applyNeumann, 
		OperatorFunction_ApplyRHSFunc*	applyRHS, 
		Name				name,
		Index				dataCount,
		Name*				dataName,
		Bool*				dataIsComponent,
		Type*				dataType,
		Bool*				dataIsRequired,
		void**				values )
{
	return _OperatorFunction_New( 
		sizeof( OperatorFunction ), 
		OperatorFunction_Type, 
		_OperatorFunction_Delete,
		_OperatorFunction_Print, 
		NULL,
	        _OperatorFunction_DefaultNew,
		_OperatorFunction_Construct,
		_OperatorFunction_Build,
		_OperatorFunction_Initialise,
		_OperatorFunction_Execute,
		_OperatorFunction_Destroy,	
		applyMatrix,
		applyNeumann,
		applyRHS,
		dataCount,
		dataName,
		dataIsComponent,
		dataType,
		dataIsRequired,
		values,
		name );
}

void OperatorFunction_Init(
		OperatorFunction*		self, 
		OperatorFunction_ApplyFunc*	applyMatrix, 
		OperatorFunction_ApplyFunc*	applyNeumann, 
		OperatorFunction_ApplyRHSFunc*	applyRHS, 
		Index				dataCount,
		Name*				dataName,
		Bool*				dataIsComponent,
		Type*				dataType,
		Bool*				dataIsRequired,
		void**				values )
{
	/* General info */
	self->type = OperatorFunction_Type;
	self->_sizeOfSelf = sizeof( OperatorFunction );
	self->_deleteSelf = False;
	
	/* Virtual info */
	self->_delete = _OperatorFunction_Delete;
	self->_print = _OperatorFunction_Print;
	self->_copy = NULL;
	self->_defaultConstructor = _OperatorFunction_DefaultNew;
	self->_construct = _OperatorFunction_Construct;
	self->_build = _OperatorFunction_Build;
	self->_initialise = _OperatorFunction_Initialise;
	self->_execute = _OperatorFunction_Execute;
	self->_destroy = _OperatorFunction_Destroy;

	_Stg_Component_Init( (Stg_Component*)self );
	
	/* Stg_Class info */
	_OperatorFunction_Init(
		self, 
		applyMatrix, 
		applyNeumann, 
		applyRHS, 
		dataCount,
		dataName,
		dataIsComponent,
		dataType,
		dataIsRequired,
		values );
}








/*
OperatorFunction* OperatorFunction_New(
		OperatorFunction_ApplyFunc*	apply, 
		Name				name,
		Index				dataCount,
		Name*				dataName,
		Bool*				dataIsComponent,
		Type*				dataType,
		Bool*				dataIsRequired,
		void**				values )
{
	return _OperatorFunction_New( 
		sizeof( OperatorFunction ), 
		OperatorFunction_Type, 
		_OperatorFunction_Delete,
		_OperatorFunction_Print, 
		NULL, 
		apply, 
		name,
		dataCount,
		dataName,
		dataIsComponent,
		dataType,
		dataIsRequired,
		values );
}*/

/*
void OperatorFunction_Init(
		OperatorFunction*		self, 
		OperatorFunction_ApplyFunc*	apply, 
		Name				name,
		Index				dataCount,
		Name*				dataName,
		Bool*				dataIsComponent,
		Type*				dataType,
		Bool*				dataIsRequired,
		void**				values )
{*/
	/* General info */
/*
	self->type = OperatorFunction_Type;
	self->_sizeOfSelf = sizeof( OperatorFunction );
	self->_deleteSelf = False;
*/	
	/* Virtual info */
/*
	self->_delete = _OperatorFunction_Delete;
	self->_print = _OperatorFunction_Print;
	self->_copy = NULL;
	
	_Stg_Class_Init( (Stg_Class*)self );
*/	
	/* Stg_Class info */
/*
	_OperatorFunction_Init(
		self, 
		apply, 
		name,
		dataCount,
		dataName,
		dataIsComponent,
		dataType,
		dataIsRequired,
		values );
}*/

/*
OperatorFunction* _OperatorFunction_New( 
		SizeT				_sizeOfSelf, 
		Type				type,
		Stg_Class_DeleteFunction*	_delete,
		Stg_Class_PrintFunction*	_print, 
		Stg_Class_CopyFunction*		_copy, 
		OperatorFunction_ApplyFunc*	apply, 
		Name				name,
		Index				dataCount,
		Name*				dataName,
		Bool*				dataIsComponent,
		Type*				dataType,
		Bool*				dataIsRequired,
		void**				values )
{
	OperatorFunction* self;
*/	
	/* Allocate memory */
/*
	assert( _sizeOfSelf >= sizeof( OperatorFunction ) );
	self = (OperatorFunction*)_Stg_Class_New( _sizeOfSelf, type, _delete, _print, _copy );
*/	
	/* General info */
	
	/* Virtual functions */
	
	/* Stg_Class info */
/*
	_OperatorFunction_Init(
		self, 
		apply, 
		name,
		dataCount,
		dataName,
		dataIsComponent,
		dataType,
		dataIsRequired,
		values );
	
	return self;
}
*/

/*
void _OperatorFunction_Init(
		void*				operatorFunction, 
		OperatorFunction_ApplyFunc*	apply, 
		Name				name,
		Index				dataCount,
		Name*				dataName,
		Bool*				dataIsComponent,
		Type*				dataType,
		Bool*				dataIsRequired,
		void**				values )
{
	OperatorFunction* self = (OperatorFunction*)operatorFunction;
	
	self->apply = apply;
	self->name = name;
	self->_dataCount = dataCount;
	if( self->_dataCount ) {
		Index i;*/

		/* Copy names */
/*
		self->_dataName = Memory_Alloc_Array( Name, self->_dataCount, "OperatorFunction->_dataName" );
		for( i = 0; i < self->_dataCount; i++ ) {
			self->_dataName[i] = StG_Strdup( dataName[i] );
		}
*/		
		/* Copy isComponent */
/*
		self->_dataIsComponent = Memory_Alloc_Array( Bool, self->_dataCount, "OperatorFunction->_dataIsComponent" );
		memcpy( self->_dataIsComponent, dataIsComponent, self->_dataCount * sizeof( Bool ) );
*/
		/* Copy dataType */
/*
		self->_dataType = Memory_Alloc_Array( Type, self->_dataCount, "OperatorFunction->_dataType" );
		for( i = 0; i < self->_dataCount; i++ ) {
			self->_dataType[i] = StG_Strdup( dataType[i] );
		}
*/
		/* Copy isRequired */
/*
		self->_dataIsRequired = Memory_Alloc_Array( Bool, self->_dataCount, "OperatorFunction->_dataIsRequired" );
		memcpy( self->_dataIsRequired, dataIsRequired, self->_dataCount * sizeof( Bool ) );
*/
		/* There is no way of telling how to destroy the "values" (void**), so we wont copy the items... 
		   we just assume ownership when it comes to deleting. */
/*
	}
	else {
		self->_dataName = NULL;
		self->_dataIsComponent = NULL;
		self->_dataType = NULL;
		self->_dataIsRequired = NULL;
	}
	self->_values = values;*/ /*See note above */
/*
}*/


OperatorFunction* OperatorFunction_NewFromMeta( OperatorFunction_ApplyFunc* 	applyMatrix, 
						OperatorFunction_ApplyFunc* 	applyNeumann,
						OperatorFunction_ApplyRHSFunc* 	applyRHS,
						Name name ) {
	OperatorFunction* 	self;
	Dictionary* 		meta;
	Dictionary_Entry_Value* paramList;
	Dictionary_Entry_Value* dependList;

	meta = Stg_ComponentRegister_GetMetadata( stgComponentRegister, name, 0/* not sure about this?? */ );

	if( meta ) {
		Index paramCount;
		Index usesCount;
		Index i;
		Index count;
		Name* names;
		Bool* isComponents;
		Type* types;
		Bool* isRequired;
		void** values;

		paramList = Dictionary_Get( meta, "Params" ); /* was Params before */
		paramCount = Dictionary_Entry_Value_GetCount( paramList ); /* Stg_Meta_GetParameterCount() ?? */
		dependList = Dictionary_Get( meta, "Dependencies" ); /* was Dependencies before */
		usesCount = Dictionary_Entry_Value_GetCount( dependList );
		count = paramCount + usesCount;
		names = Memory_Alloc_Array( Name, count, "OperatorFunction_NewFromMeta: names" );
		isComponents = Memory_Alloc_Array( Bool, count, "OperatorFunction_NewFromMeta: isComponents" );
		types = Memory_Alloc_Array( Type, count, "OperatorFunction_NewFromMeta: types" );
		isRequired = Memory_Alloc_Array( Bool, count, "OperatorFunction_NewFromMeta: isRequired" );

		values = Memory_Alloc_Array( void*, count, "OperatorFunction_NewFromMeta: values" );
		memset( values, 0, count * sizeof(count) );
		
		/* Obtain the parameters to the operator from the meta file information */
		for( i = 0; i < paramCount; i++ ) {
			/*Stg_ComponentMeta_Value* param;*/
			char* defaultValue;

			/* Note: we're just ptr copying the strings here because we'll delete "meta" after we've created the operator function 
				(i.e. the copying of values occurs there.) */
			names[i] = Stg_Meta_GetParameterName( meta, i );
			isComponents[i] = False;
			types[i] = Stg_Meta_GetParameterType( meta, i );
			isRequired[i] = False; /* not sure about this?? */
			defaultValue = Stg_Meta_GetParameterDefault( meta, i );

			/* Evaluate and store the default value for the parameter if provided */
			if( defaultValue ) {
				OperatorFunction_Datum v;
				ConvertTypedStrValue( types[i], defaultValue, &v );

				if( strcasecmp( types[i], "String" ) == 0 ) {
					values[i] = v._ptr;
				}
				else if( strcasecmp( types[i], "Int" ) == 0 ) {
					values[i] = Memory_Alloc( int, "OperatorFunction_NewFromMeta: values[]" );
					*((int*)values[i]) = v._int;
					/* Note: this is the same technique used by Dictionary_Entry_Value 22/06/2008 */
				}
				else if( strcasecmp( types[i], "Unsigned" ) == 0 ) {
					values[i] = Memory_Alloc( unsigned, "OperatorFunction_NewFromMeta: values[]" );
					*((int*)values[i]) = v._unsigned;
					/* Note: this is the same technique used by Dictionary_Entry_Value 22/06/2008 */
				}
				else if( strcasecmp( types[i], "Double" ) == 0 ) {
					values[i] = Memory_Alloc( double, "OperatorFunction_NewFromMeta: values[]" );
					*((double*)(values[i])) = v._double;
					/* Note: this is the same technique used by Dictionary_Entry_Value 22/06/2008 */
				}
				else {
					Journal_Firewall( 
						(0), 
						Journal_Register( Error_Type, OperatorFunction_Type ),
						"In func %s: For operator \"%s\", the \"%s\" parameter's type \"%s\" not yet implemented!",
						__func__,
						name,						
						names[i],
						types[i] );
				}
			}
		}

		/* Obtain the "uses" of the operator from the meta file information */
		for( i = 0; i < usesCount; i++ ) { 
			/*Stg_ComponentMeta_Value* uses;*/

			/* Note: we're just ptr copying the strings here because we'll delete "meta" after we've created the operator function 
				(i.e. the copying of values occurs there. */
			names[paramCount + i] = Stg_Meta_GetAssociationName( meta, i );
			isComponents[paramCount + i] = True;
			types[paramCount + i] = (Type)Stg_Meta_GetAssociationType( meta, i );

			Dictionary* associations = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( meta, "associations" ) );
			Dictionary* association  = Dictionary_Entry_Value_AsDictionary( Dictionary_GetByIndex( associations, i ) );
			isRequired[paramCount + i] = Dictionary_Entry_Value_AsBool( Dictionary_Get( association, "essential" ) );

			/* allocate memory for the dependencies */
			/* need to call the actual function, not the macro, as the type is being passed in as a string... dave - 02.02.09 */
			values[paramCount + i] = _Memory_InternalMalloc( sizeof(types[paramCount + i]) );
			
		}

		self = OperatorFunction_New( applyMatrix, applyNeumann, applyRHS, name, count, names, isComponents, types, isRequired, values );

		Memory_Free( isRequired );
		Memory_Free( types );
		Memory_Free( isComponents );
		Memory_Free( names );
	}
	else {
		self = OperatorFunction_New( applyMatrix, applyNeumann, applyRHS, name, 0, NULL, NULL, NULL, NULL, NULL );
	}

	return self;
}

/*--------------------------------------------------------------------------------------------------------------------------
** General virtual functions
*/

void _OperatorFunction_Delete( void* operatorFunction ) {
	OperatorFunction* 	self 	= (OperatorFunction*)operatorFunction;
	
	/* As a Operator Function... */
	if( self->_values ) { /* Assume ownership and destroy all. Assume list members are can be fully destroyed by Memory_Free. */
		Index i;

		for( i = 0; i < self->_dataCount; i++ ) {
			if( self->_values[i]._ptr ) {
				Memory_Free( self->_values[i]._ptr );
				self->_values[i]._ptr = NULL;
			}
		}
		Memory_Free( self->_values );
		self->_values = NULL;
	}
	if( self->_dataIsRequired ) {
		Memory_Free( self->_dataIsRequired );
		self->_dataIsRequired = NULL;
	}
	if( self->_dataType ) {
		Index i;

		for( i = 0; i < self->_dataCount; i++ ) {
			if( self->_dataType[i] ) {
				Memory_Free( self->_dataType[i] );
				self->_dataType[i] = NULL;
			}
		}
		Memory_Free( self->_dataType );
		self->_dataType = NULL;
	}
	if( self->_dataIsComponent ) {
		Memory_Free( self->_dataIsComponent );
		self->_dataIsComponent = NULL;
	}
	if( self->_dataName ) {
		Index i;

		for( i = 0; i < self->_dataCount; i++ ) {
			if( self->_dataName[i] ) {
				Memory_Free( self->_dataName[i] );
				self->_dataName[i] = NULL;
			}
		}
		Memory_Free( self->_dataName );
		self->_dataName = NULL;
	}
	self->_dataCount = 0;

	/* Delete parent class */
	_Stg_Class_Delete( self );
}


void _OperatorFunction_Print( void* operatorFunction, Stream* stream ) {
	OperatorFunction*	self 	= (OperatorFunction*)operatorFunction;
	
	/* Set the Journal for printing informations */
	Stream* operatorFunctionStream = stream;
	
	/* General info */
	Journal_Printf( operatorFunctionStream, "OperatorFunction (ptr): %p\n", self );
	
	/* Virtual info */
	
	/* Stg_Class info */
	Journal_Printf( operatorFunctionStream, "\tapply Matrix (func ptr): %p\n", self->applyMatrix );
	Journal_Printf( operatorFunctionStream, "\tapply Neumann (func ptr): %p\n", self->applyNeumann );
	Journal_Printf( operatorFunctionStream, "\tapply RHS (func ptr): %p\n", self->applyRHS );
	Journal_Printf( operatorFunctionStream, "\tname (ptr): %p\n", self->name );
	if( self->name ) {
		Journal_Printf( operatorFunctionStream, "\t\tname: %s\n", self->name );
	}

	Journal_Printf( operatorFunctionStream, "\t_dataCount: %d\n", self->_dataCount );
	Journal_Printf( operatorFunctionStream, "\t_dataName (ptr): %p\n", self->_dataName );
	if( self->_dataName ) {
		Index i;

		for( i = 0; i < self->_dataCount; i++ ) {
			Journal_Printf( operatorFunctionStream, "\t\t_dataName[%d]: %s\n", i, self->_dataName[i] );
		}
	}
	Journal_Printf( operatorFunctionStream, "\t_dataIsComponent (ptr): %p\n", self->_dataIsComponent );
	if( self->_dataIsComponent ) {
		Index i;

		for( i = 0; i < self->_dataCount; i++ ) {
			Journal_Printf( operatorFunctionStream, "\t\t_dataIsComponent[%d]: %s\n", i, self->_dataIsComponent[i] ? "Yes" : "No" );
		}
	}
	Journal_Printf( operatorFunctionStream, "\t_dataType (ptr): %p\n", self->_dataType );
	if( self->_dataType ) {
		Index i;

		for( i = 0; i < self->_dataCount; i++ ) {
			Journal_Printf( operatorFunctionStream, "\t\t_dataType: %s\n", i, self->_dataType[i] );
		}
	}
	Journal_Printf( operatorFunctionStream, "\t_dataIsRequired (ptr): %p\n", self->_dataIsRequired );
	if( self->_dataIsRequired ) {
		Index i;

		for( i = 0; i < self->_dataCount; i++ ) {
			Journal_Printf( operatorFunctionStream, "\t\t_dataIsRequired[%d]: %s\n", i, self->_dataIsRequired[i] ? "Yes" : "No" );
		}
	}
	Journal_Printf( operatorFunctionStream, "\t_values (ptr): %p\n", self->_values );
	if( self->_values ) {
		Index i;

		for( i = 0; i < self->_dataCount; i++ ) {
			Journal_Printf( operatorFunctionStream, "\t\t_values[%d]: (ptr) %p\n", i, self->_values[i] );
		}
	}

	/* Print parent class */
	_Stg_Class_Print( self, operatorFunctionStream );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Macros
*/


/*--------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void* _OperatorFunction_DefaultNew( Name name ) {
	return _OperatorFunction_New( 
		sizeof( OperatorFunction ), 
		OperatorFunction_Type, 
		_OperatorFunction_Delete,
		_OperatorFunction_Print, 
		NULL,
	        _OperatorFunction_DefaultNew,
		_OperatorFunction_Construct,
		_OperatorFunction_Build,
		_OperatorFunction_Initialise,
		_OperatorFunction_Execute,
		_OperatorFunction_Destroy,	
		NULL, /* _assembleMatrix function */
		NULL, /* _assembleNeumann function */
		NULL, /* _assembleRHS function */
		0,    /* dataCount */
		NULL, /* dataName */
		NULL, /* dataIsComponent */
		NULL, /* dataType */
		NULL, /* dataIsRequired */
		NULL, /* values */
		name );
}

void _OperatorFunction_Construct( void* operatorFunction, Stg_ComponentFactory* cf, void* data ) {
	OperatorFunction* 	self 		= (OperatorFunction*) operatorFunction;
	Index			i;
	Stream*			errStream 	= Journal_Register( ErrorStream_Type, OperatorFunction_Type );

	Stg_Component_Construct( self, cf, data, False );

	for( i = 0; i < self->_dataCount; i++ ) {
		if( self->_dataIsComponent[i] )
			/* need to call the actual function, not the macro, as the type is being passed in as a string, not a type... dave - 02.02.09 */
			self->_values[i]._ptr = _Stg_ComponentFactory_ConstructByKey( cf, self->name, self->_dataName[i], 
								(Type)self->_dataType[i], self->_dataIsRequired[i], data );
		else if( !strcasecmp( self->_dataType[i], "Double" ) )
			//*((double*)self->_values[i]) = Stg_ComponentFactory_GetDouble( cf, self->name, self->_dataName[i], *((double*)(self->_values[i])) );
			self->_values[i]._double = Stg_ComponentFactory_GetDouble( cf, self->name, self->_dataName[i], self->_values[i]._double );
		else if( !strcasecmp( self->_dataType[i], "String" ) )
			self->_values[i]._ptr = Stg_ComponentFactory_GetString( cf, self->name, self->_dataName[i], (char*)self->_values[i]._ptr );
		else if( !strcasecmp( self->_dataType[i], "Int" ) )
			self->_values[i]._int = Stg_ComponentFactory_GetInt( cf, self->name, self->_dataName[i], self->_values[i]._int );
		else if( !strcasecmp( self->_dataType[i], "UnsignedInt" ) )
			self->_values[i]._unsigned = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, self->_dataName[i], self->_values[i]._unsigned );
		else {
			Journal_Printf( errStream, "the type: %s not currently handled by the OperatorFunction\n", self->_dataType[i] );
			assert( 0 );
		}
	}

	/* assign the StiffnessMatrix (same as for the StiffnessMatrixTerms...) */
	self->sm = Stg_ComponentFactory_ConstructByKey( cf, self->name, "StiffnessMatrix", StiffnessMatrix, True, data );
	StiffnessMatrix_AddOperatorFunction( self->sm, self );

	self->integrationSwarm 	= Stg_ComponentFactory_ConstructByKey( cf, self->name, "IntegrationSwarm", Swarm, True, data );
	self->borderSwarm	= Stg_ComponentFactory_ConstructByKey( cf, self->name, "BorderSwarm", Swarm, False, data );
}

void _OperatorFunction_Build( void* operatorFunction, void* data ) {

}

void _OperatorFunction_Initialise( void* operatorFunction, void* data ) {

}

void _OperatorFunction_Execute( void* operatorFunction, void* data ) {

}

void _OperatorFunction_Destroy( void* operatorFunction, void* data ) {

}

/*--------------------------------------------------------------------------------------------------------------------------
** Build functions
*/


/*--------------------------------------------------------------------------------------------------------------------------
** Functions
*/

/* OperatorFunction_Apply is implemented as a macro */


OperatorFunction_Datum* OperatorFunction_CreateDataFromDictionary( 
		void*				operatorFunction,
		Stg_ComponentFactory*		cf, 
		Dictionary*			dictionary ) 
{
	OperatorFunction*	self = (OperatorFunction*)operatorFunction;
	OperatorFunction_Datum*	data;
	Index 			i;

	data = (OperatorFunction_Datum*)Memory_Alloc_Array( 
		OperatorFunction_Datum, 
		self->_dataCount,
		"OperatorFunction_CreateDataFromDictionary:data" );

	/* For each data item defined for this operator, populate the "data list" with actual values provided by the input file /
	    dictionary. */ 
	for( i = 0; i < self->_dataCount; i++ ) {
		/* Look the data item up in the dictionary, and ensure its present if required. */
		Dictionary_Entry_Value* dev = Dictionary_Get( dictionary, self->_dataName[i] );
		Journal_Firewall( 
			(!dev && self->_dataIsRequired), 
			Journal_Register( Error_Type, self->type ),
			"In func %s: For operator \"%s\", the \"%s\" parameter/dependancy/uses is required/essential, but not provided in input/model",
			__func__,
			self->name,
			self->_dataName[i] );

		/* If the data item is a dependancy/uses, look it up in the live component register. */
		if( self->_dataIsComponent[i] ) { 
			data[i]._ptr = LiveComponentRegister_Get( cf->LCRegister, Dictionary_Entry_Value_AsString( dev ) );
			Journal_Firewall( 
				(!data[i]._ptr && self->_dataIsRequired), 
				Journal_Register( Error_Type, self->type ),
				"In func %s: For operator \"%s\", the \"%s\" dependancy/uses is required/essential, but component not found to exist",
				__func__,
				self->name,
				self->_dataName[i] );
			/* TODO: make sure its the right type! */
		}
		/* Else its a parameter, and if its not provided in the input file, use the defined default. */
		else { 
			Journal_Firewall( 
				(!dev && (!self->_values[i]._int || !self->_values[i]._unsigned || fabs(self->_values[i]._double) > 1.0E-8)), 
				Journal_Register( Error_Type, self->type ),
				"In func %s: For operator \"%s\", the \"%s\" parameter was not provided, and there is no default",
				__func__,
				self->name,
				self->_dataName[i] );
			if( dev ) {
				ConvertTypedStrValue( self->_dataType[i], Dictionary_Entry_Value_AsString( dev ), &data[i] );
			}
			else if( self->_values[i]._ptr || self->_values[i]._int || self->_values[i]._unsigned || 
					fabs(self->_values[i]._double) > 1.0E-8 ) {
				memcpy( &data[i], &self->_values[i], sizeof(OperatorFunction_Datum) );
			}
		}
	}

	return data;
}




void ConvertTypedStrValue( Type type, const char* inValue, OperatorFunction_Datum* outValue ) {
	if( strcasecmp( type, "String" ) == 0 ) {
		outValue->_ptr = StG_Strdup( inValue );
	}
	else if( strcasecmp( type, "Int" ) == 0 ) {
		outValue->_int = strtoul( inValue, 0, 0 );
		/* Note: this is the same technique used by Dictionary_Entry_Value 22/06/2008 */
	}
	else if( strcasecmp( type, "Double" ) == 0 ) {
		outValue->_double = strtod( inValue, 0 );
		/* Note: this is the same technique used by Dictionary_Entry_Value 22/06/2008 */
	}
	else if( strcasecmp( type, "Unsigned" ) == 0 ) {
		outValue->_unsigned = strtod( inValue, 0 );
		/* Note: this is the same technique used by Dictionary_Entry_Value 22/06/2008 */
	}
	else {
		Journal_Firewall( 
			(0), 
			Journal_Register( Error_Type, OperatorFunction_Type ),
			"In func %s: The parameter type \"%s\" not yet implemented!",
			__func__,
			type );
	}
}


