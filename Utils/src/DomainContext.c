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
** $Id: Context.c 2454 2004-12-21 04:50:42Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
                                                                                                                                    
#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>
                                                                                                                                    
#include "types.h"
#include "FieldVariable.h"
#include "FieldVariable_Register.h"
#include "DomainContext.h"
                                                                                                                                    
#include <stdlib.h>
#include <string.h>
#include <assert.h>


/* Textual name of this class */
const Type DomainContext_Type = "DomainContext";

DomainContext* DomainContext_New( 
	Name			name,
	double		start,
	double		stop,
	MPI_Comm		communicator,
	Dictionary*	dictionary )
{
	DomainContext* self = _DomainContext_DefaultNew( name );

	self->isConstructed = True;
	_AbstractContext_Init( (AbstractContext*) self );
	_DomainContext_Init( self );

	return self;
}

DomainContext* _DomainContext_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(DomainContext);
	Type                                                      type = DomainContext_Type;
	Stg_Class_DeleteFunction*                              _delete = _DomainContext_Delete;
	Stg_Class_PrintFunction*                                _print = _DomainContext_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = NULL;
	Stg_Component_ConstructFunction*                    _construct = _DomainContext_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _AbstractContext_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _AbstractContext_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _AbstractContext_Execute;
	Stg_Component_DestroyFunction*                        _destroy = (Stg_Component_DestroyFunction*)_DomainContext_Destroy;
	AllocationType                              nameAllocationType = NON_GLOBAL;
	AbstractContext_SetDt*                                  _setDt = _DomainContext_SetDt;
	double                                               startTime = 0;
	double                                                stopTime = 0;
	MPI_Comm                                          communicator = MPI_COMM_WORLD;
	Dictionary*                                         dictionary = NULL;

	return _DomainContext_New(  DOMAINCONTEXT_PASSARGS  );
}

DomainContext* _DomainContext_New(  DOMAINCONTEXT_DEFARGS  ) {
	DomainContext* self;
	
	/* Allocate memory */
	self = (DomainContext*)_AbstractContext_New(  ABSTRACTCONTEXT_PASSARGS  );
	
	/* General info */

	/* Virtual info */
	
	return self;
}

void _DomainContext_Init( DomainContext* self ) {

	self->fieldVariable_Register = FieldVariable_Register_New();
}


/* Virtual Functions -------------------------------------------------------------------------------------------------------------*/

void _DomainContext_AssignFromXML( void* context, Stg_ComponentFactory* cf, void* data ) {
	DomainContext* self = (DomainContext*)context;

	_AbstractContext_AssignFromXML( context, cf, data );

	/* Check if we have been provided a constant to multiply our calculated dt values by. */
	self->dtFactor = Dictionary_GetDouble_WithDefault( self->dictionary, (Dictionary_Entry_Key)"timestepFactor", 1.0  );
	self->dim = Dictionary_GetUnsignedInt_WithDefault( self->dictionary, "dim", 2 );

	_DomainContext_Init( self );
}

void _DomainContext_Delete( void* context ) {
	DomainContext* self = (DomainContext*)context;
	
	Journal_DPrintf( self->debug, "In: %s()\n", __func__ );

	Journal_DPrintfL( self->debug, 2, "Deleting the FieldVariable register (and hence all FieldVariables).\n" );

	/* Stg_Class_Delete parent */
	_AbstractContext_Delete( self );
}

void _DomainContext_Destroy( void* context ) {
	DomainContext* self = (DomainContext*)context;
	
	Stg_Class_Delete( self->fieldVariable_Register ); 
	_AbstractContext_Destroy( self, 0 );
}

void _DomainContext_Print( void* context, Stream* stream ) {
	DomainContext* self = (DomainContext*)context;
	
	/* General info */
	Journal_Printf( (void*) stream, "DomainContext (ptr): %p\n", self );
	
	/* Print parent */
	_AbstractContext_Print( self, stream );

	Journal_Printf( (void*) stream, "\tfieldVariables (ptr): %p\n", self->fieldVariable_Register );
	Stg_Class_Print( self->fieldVariable_Register, stream );
}


void _DomainContext_SetDt( void* context, double dt ) {
}


