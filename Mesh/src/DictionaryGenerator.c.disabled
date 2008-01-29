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
** $Id: DictionaryGenerator.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>
#include <StGermain/StGermain.h>

#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>

#include "Mesh.h"


/* Textual name of this class */
const Type DictionaryGenerator_Type = "DictionaryGenerator";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

DictionaryGenerator* DictionaryGenerator_New( Name name ) {
	return _DictionaryGenerator_New( sizeof(DictionaryGenerator), 
					DictionaryGenerator_Type, 
					_DictionaryGenerator_Delete, 
					_DictionaryGenerator_Print, 
					 NULL, 
					(void* (*)(Name))_DictionaryGenerator_New, 
					_DictionaryGenerator_Construct, 
					_DictionaryGenerator_Build, 
					_DictionaryGenerator_Initialise, 
					_DictionaryGenerator_Execute, 
					_DictionaryGenerator_Destroy, 
					name, 
					NON_GLOBAL, 
					DictionaryGenerator_Generate );
}

DictionaryGenerator* _DictionaryGenerator_New( DICTIONARYGENERATOR_DEFARGS ) {
	DictionaryGenerator* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(DictionaryGenerator) );
	self = (DictionaryGenerator*)_MeshGenerator_New( MESHGENERATOR_PASSARGS );

	/* Virtual info */

	/* DictionaryGenerator info */
	_DictionaryGenerator_Init( self );

	return self;
}

void _DictionaryGenerator_Init( DictionaryGenerator* self ) {
	self->dict = NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _DictionaryGenerator_Delete( void* generator ) {
	DictionaryGenerator*	self = (DictionaryGenerator*)generator;

	/* Delete the parent. */
	_MeshGenerator_Delete( self );
}

void _DictionaryGenerator_Print( void* generator, Stream* stream ) {
	DictionaryGenerator*	self = (DictionaryGenerator*)generator;
	
	/* Set the Journal for printing informations */
	Stream* generatorStream;
	generatorStream = Journal_Register( InfoStream_Type, "DictionaryGeneratorStream" );

	/* Print parent */
	Journal_Printf( stream, "DictionaryGenerator (ptr): (%p)\n", self );
	_MeshGenerator_Print( self, stream );
}

void _DictionaryGenerator_Construct( void* generator, Stg_ComponentFactory* cf, void* data ) {
	DictionaryGenerator*	self = (DictionaryGenerator*)generator;
	Dictionary*		dict;

	assert( self );
	assert( cf );

	/* Call parent construct. */
	_MeshGenerator_Construct( self, cf, data );

	/* Set the dictionary to the component dictionary. */
	DictionaryGenerator_SetDictionary( self, cf->componentDict );
}

void _DictionaryGenerator_Build( void* generator, void* data ) {
	_MeshGenerator_Build( generator, data );
}

void _DictionaryGenerator_Initialise( void* generator, void* data ) {
}

void _DictionaryGenerator_Execute( void* generator, void* data ) {
}

void _DictionaryGenerator_Destroy( void* generator, void* data ) {
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void DictionaryGenerator_Generate( void* generator, void* _mesh ) {
	DictionaryGenerator*	self = (DictionaryGenerator*)generator;
	Mesh*			mesh = (Mesh*)_mesh;

	/* Sanity check. */
	assert( self );
	assert( self->dict );

	/* For now, this can only work in serial. */
	
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
