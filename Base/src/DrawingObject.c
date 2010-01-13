/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005, Monash Cluster Computing 
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
** 		* Redistributions of source code must retain the above copyright notice, 
** 			this list of conditions and the following disclaimer.
** 		* Redistributions in binary form must reproduce the above copyright 
**			notice, this list of conditions and the following disclaimer in the 
**			documentation and/or other materials provided with the distribution.
** 		* Neither the name of the Monash University nor the names of its contributors 
**			may be used to endorse or promote products derived from this software 
**			without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
** THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
** PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS 
** BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
** CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
** SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
** HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
** LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
** OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
**
** Contact:
*%		Cecile Duboz - Cecile.Duboz@sci.monash.edu.au
*%
** Contributors:
*+		Cecile Duboz
*+		Robert Turnbull
*+		Alan Lo
*+		Louis Moresi
*+		David Stegman
*+		David May
*+		Stevan Quenette
*+		Patrick Sunter
*+		Greg Watson
*+
** $Id: DrawingObject.c 740 2007-10-11 08:05:31Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "types.h"
#include "ColourMap.h"
#include "ViewportInfo.h"
#include "Window.h"
#include "DrawingObject.h"
#include "Init.h"

#include <assert.h>
#include <string.h>

const Type lucDrawingObject_Type = "lucDrawingObject";

lucDrawingObject* _lucDrawingObject_New(  LUCDRAWINGOBJECT_DEFARGS  )
{
	lucDrawingObject*    self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucDrawingObject) );
	/* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
	/* This means that any values of these parameters that are passed into this function are not passed onto the parent function
	   and so should be set to ZERO in any children of this class. */
	nameAllocationType = NON_GLOBAL;

	self = (lucDrawingObject*) _Stg_Component_New(  STG_COMPONENT_PASSARGS  );

	self->_setup   = _setup;
	self->_draw    = _draw;
	self->_cleanUp = _cleanUp;
	
	return self;
}

void _lucDrawingObject_Init( lucDrawingObject* self ) {
	self->isConstructed = True;

	self->infoStream  = Journal_MyStream( Info_Type,  self );
	self->errorStream = Journal_MyStream( Error_Type, self );
	self->debugStream = Journal_MyStream( Debug_Type, self );
}

void lucDrawingObject_InitAll( void* drawingObject ) {
	lucDrawingObject* self        = drawingObject;

	_lucDrawingObject_Init( self );
}

void _lucDrawingObject_Delete( void* drawingObject ) {
	lucDrawingObject* self        = drawingObject;
	
	_Stg_Component_Delete( self );
}

void _lucDrawingObject_Print( void* drawingObject, Stream* stream ) {
	lucDrawingObject* self        = drawingObject;
	
	Journal_Printf( stream, "lucDrawingObject: %s\n", self->name );

	Stream_Indent( stream );

	/* Print Parent */
	_Stg_Component_Print( self, stream );

	Journal_PrintPointer( stream, self->_setup );
	Journal_PrintPointer( stream, self->_draw );
	
	Stream_UnIndent( stream );
}

void* _lucDrawingObject_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	lucDrawingObject* self        = drawingObject;
	lucDrawingObject* newDrawingObject;

	newDrawingObject = _Stg_Component_Copy( self, dest, deep, nameExt, ptrMap );

	return (void*) newDrawingObject;
}

void _lucDrawingObject_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data ) {
	lucDrawingObject*        self            = (lucDrawingObject*) drawingObject ;

	self->context = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Context", AbstractContext, False, data );
	if( !self->context  ) 
		self->context = Stg_ComponentFactory_ConstructByName( cf, (Name)"context", AbstractContext, True, data );
	
	_lucDrawingObject_Init( self  );
}

void _lucDrawingObject_Build( void* camera, void* data ) { }
void _lucDrawingObject_Initialise( void* camera, void* data ) { }
void _lucDrawingObject_Execute( void* camera, void* data ) { }
void _lucDrawingObject_Destroy( void* camera, void* data ) { }

void lucDrawingObject_Setup( void* drawingObject, void* context ) {
	lucDrawingObject*   self       = (lucDrawingObject*) drawingObject ;
   double time;

	lucDebug_PrintFunctionBegin( self, 2 );

	time = MPI_Wtime();

	if ( self->needsToSetup ) 
		self->_setup( self, context );
	else {
		Journal_DPrintfL( lucDebug, 2, "Already setup\n" );
	}

	self->needsToSetup   = False;
	self->needsToCleanUp = True;

	time = MPI_Wtime() - time;
	Journal_DPrintfL( lucDebug, 2, "Setup took %f seconds, ", time );

	lucDebug_PrintFunctionEnd( self, 2 );
}

void lucDrawingObject_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* context ) {
	lucDrawingObject*   self       = (lucDrawingObject*) drawingObject ;
   double time;
   
	lucDebug_PrintFunctionBegin( self, 2 );
	Journal_DPrintfL( lucDebug, 2, "Drawing (%s),  ", self->name );

	lucDrawingObject_Setup( self, context );
	

	time = MPI_Wtime();
	
	self->_draw( self, window, viewportInfo, context );
	
	time = MPI_Wtime() - time;
	Journal_DPrintfL( lucDebug, 2, "Drawing took %f seconds\n", time );

	lucDebug_PrintFunctionEnd( self, 2 );
}


void lucDrawingObject_CleanUp( void* drawingObject, void* context ) {
	lucDrawingObject*   self       = (lucDrawingObject*) drawingObject ;

	lucDebug_PrintFunctionBegin( self, 2 );

	if ( self->needsToCleanUp ) {
		self->_cleanUp( self, context );

		self->needsToCleanUp = False;
	}
	else {
		Journal_DPrintfL( lucDebug, 2, "%s '%s' has already cleaned up.\n", self->type, self->name );
	}

	lucDebug_PrintFunctionEnd( self, 2 );
}

void lucDrawingObjectMask_Construct( lucDrawingObjectMask* self, Name drawingObjectName, Stg_ComponentFactory* cf, void* mask ) {
	Name                   maskTypeName;

	self->value     = Stg_ComponentFactory_GetDouble( cf, drawingObjectName, (Dictionary_Entry_Key)"maskValue", 0.0  );
	self->tolerance = Stg_ComponentFactory_GetDouble( cf, drawingObjectName, (Dictionary_Entry_Key)"maskTolerance", 0.001  );

	maskTypeName = Stg_ComponentFactory_GetString( cf, drawingObjectName, (Dictionary_Entry_Key)"maskType", "GreaterThan"  );
	if ( strcasecmp( maskTypeName, "GreaterThan" ) == 0 ) 
		self->type = GreaterThan;
	else if ( strcasecmp( maskTypeName, "LesserThan" ) == 0 || strcasecmp( maskTypeName, "LessThan" ) == 0 )
		self->type = LessThan;
	else if ( strcasecmp( maskTypeName, "EqualTo" ) == 0 )
		self->type = EqualTo;
	else {
		Journal_Printf( lucError, "In func %s: Cannot understand 'maskType' '%s'.\n", __func__, maskTypeName );
		abort();
	}
}

Bool lucDrawingObjectMask_Test( lucDrawingObjectMask* self, double value ) {
	double maskValue = self->value;

	switch (self->type) {
		case GreaterThan:
			if (value > maskValue) 
				return True;
			return False;
		case LessThan:
			if (value < maskValue) 
				return True;
			return False;
		case EqualTo:
			if (fabs( maskValue - value ) < self->tolerance )
				return True;
			return False;
	}
	return True;
}


