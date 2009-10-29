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
** $Id: Viewport.c 791 2008-09-01 02:09:06Z JulianGiordani $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "types.h"
#include "ColourMap.h"
#include "Window.h"
#include "DrawingObject_Register.h"
#include "Light_Register.h"
#include "Light.h"


#include "Viewport.h"

#include "DrawingObject.h"
#include "Camera.h"
#include "Init.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>


const Type lucViewport_Type = "lucViewport";

MPI_Datatype lucViewport_MPI_Datatype;

lucViewport* lucViewport_New(
		Name                                               name,
		lucCamera*                                         camera,
		lucDrawingObject**                                 drawingObjectList,
		DrawingObject_Index                                drawingObjectCount,
		lucLight**          				   lightList,
	  	Light_Index                                        lightCount,
		Bool                                               drawTitle,
		Bool                                               drawTime,
		Bool                                               compositeEachObject,
		double                                             nearClipPlane,
		double                                             farClipPlane )
{
	lucViewport* self = _lucViewport_DefaultNew( name );

	lucViewport_InitAll( self, camera, drawingObjectList, drawingObjectCount, lightList, lightCount, drawTitle, drawTime, compositeEachObject, nearClipPlane, farClipPlane );

	return self;
}

lucViewport* _lucViewport_New(
		SizeT                                              sizeOfSelf,
		Type                                               type,
		Stg_Class_DeleteFunction*                          _delete,
		Stg_Class_PrintFunction*                           _print,
		Stg_Class_CopyFunction*                            _copy, 
		Stg_Component_DefaultConstructorFunction*          _defaultConstructor,
		Stg_Component_ConstructFunction*                   _construct,
		Stg_Component_BuildFunction*                       _build,
		Stg_Component_InitialiseFunction*                  _initialise,
		Stg_Component_ExecuteFunction*                     _execute,
		Stg_Component_DestroyFunction*                     _destroy,		
		Name                                               name )
{
	lucViewport*    self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucViewport) );
	self = (lucViewport*) _Stg_Component_New( 
			sizeOfSelf,
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

	return self;
}

void _lucViewport_Init( 
		lucViewport*                                       self, 
		lucCamera*                                         camera, 
		lucDrawingObject**                                 drawingObjectList, 
		DrawingObject_Index                                drawingObjectCount,
		lucLight**          				   lightList,
	        Light_Index                                        lightCount,
		Bool                                               drawTitle,
		Bool                                               drawTime,
		Bool                                               compositeEachObject,
		double                                             nearClipPlane,
		double                                             farClipPlane )
{
	DrawingObject_Index object_I;
	Light_Index light_I;
	GLfloat lightPosition[4];
	GLfloat lmodel_ambient[4] = {0.2, 0.2, 0.2, 1.0};
	GLfloat spotCutOff = 180.0;
	GLfloat spotDirection[3] = {0.0, 0.0, -1.0};
	
	lightPosition[0]= LUC_LIGHT_DEFAULT_POS_X;
	lightPosition[1]= LUC_LIGHT_DEFAULT_POS_Y;
	lightPosition[2]= LUC_LIGHT_DEFAULT_POS_Z;
	lightPosition[3]= LUC_LIGHT_DEFAULT_POS_W;
	
	self->camera                   = camera;
	self->drawTitle                = drawTitle;
	self->drawTime                 = drawTime;
	self->nearClipPlane            = nearClipPlane;
	self->farClipPlane             = farClipPlane;
	self->compositeEachObject      = compositeEachObject;

	self->drawingObject_Register = lucDrawingObject_Register_New();

	for ( object_I = 0 ; object_I < drawingObjectCount ; object_I++ )
		lucDrawingObject_Register_Add( self->drawingObject_Register, drawingObjectList[ object_I ] );
		
	/* Setup light register stuff*/	
	self->light_Register = lucLight_Register_New();

	for ( light_I = 0 ; light_I < lightCount ; light_I++ )
		lucLight_Register_Add( self->light_Register, lightList[ light_I ] );

       	if(lightCount == 0){
       		self->defaultLight = lucLight_New( "defaultLight", 0, GL_LIGHT_MODEL_TWO_SIDE,  GL_AMBIENT_AND_DIFFUSE, lightPosition, lmodel_ambient, spotCutOff, spotDirection);
		lucLight_Register_Add( self->light_Register, self->defaultLight );
	}


}

void lucViewport_InitAll( 
		void*                                              viewport,
		lucCamera*                                         camera, 
		lucDrawingObject**                                 drawingObjectList, 
		DrawingObject_Index                                drawingObjectCount,
	        lucLight**          			           lightList,
	 	Light_Index                                        lightCount,
		Bool                                               drawTitle,
		Bool                                               drawTime,
		Bool                                               compositeEachObject,
		double                                             nearClipPlane,
		double                                             farClipPlane )
{
	lucViewport* self        = viewport;

	_lucViewport_Init( self, camera, drawingObjectList, drawingObjectCount, lightList, lightCount, drawTitle, drawTime, compositeEachObject, nearClipPlane, farClipPlane );
}

	
void _lucViewport_Delete( void* viewport ) {
	lucViewport* self        = viewport;
	
	_Stg_Component_Delete( self );
}

void _lucViewport_Print( void* viewport, Stream* stream ) {
	lucViewport*          self        = viewport;
	
	Journal_Printf( stream, "lucViewport: %s\n", self->name );

	Stream_Indent( stream );

	/* Print Parent */
	_Stg_Component_Print( self, stream );

	lucDrawingObject_Register_PrintAllObjects( self->drawingObject_Register, stream );

	Stg_Class_Print( self->camera, stream );

	Journal_PrintValue( stream, self->nearClipPlane );
	Journal_PrintValue( stream, self->farClipPlane );
	
	Journal_PrintBool( stream, self->drawTitle );
	Journal_PrintBool( stream, self->drawTime );
	Journal_PrintBool( stream, self->compositeEachObject );
	
	Stream_UnIndent( stream );
}

void* _lucViewport_Copy( void* viewport, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	lucViewport* self        = viewport;
	lucViewport* newViewport;

	newViewport = _Stg_Component_Copy( self, dest, deep, nameExt, ptrMap );
	if ( deep ) {
		newViewport->camera     = (lucCamera*)  Stg_Class_Copy( self->camera,     dest, deep, nameExt, ptrMap );
	}
	else {
		newViewport->camera        = self->camera;
	}

	newViewport->nearClipPlane       = self->nearClipPlane;
	newViewport->farClipPlane        = self->farClipPlane;
	newViewport->drawTitle           = self->drawTitle;
	newViewport->drawTime            = self->drawTime;
	newViewport->compositeEachObject = self->compositeEachObject;

	return (void*) newViewport;
}

void* _lucViewport_DefaultNew( Name name ) {
	return _lucViewport_New( 
			sizeof( lucViewport ),
			lucViewport_Type,
			_lucViewport_Delete,
			_lucViewport_Print,
			_lucViewport_Copy,
			_lucViewport_DefaultNew,
			_lucViewport_AssignFromXML,
			_lucViewport_Build,
			_lucViewport_Initialise,
			_lucViewport_Execute,
			_lucViewport_Destroy,
			name );
}

void _lucViewport_AssignFromXML( void* viewport, Stg_ComponentFactory* cf, void* data ) {
	lucViewport*        self               = (lucViewport*) viewport;
	DrawingObject_Index drawingObjectCount;
	lucDrawingObject**  drawingObjectList;
	lucLight**          lightList;
	Light_Index         lightCount;
	lucCamera*          camera;

	/* TODO Construct Parent */

	self->context = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Context", AbstractContext, False, data );
	if( !self->context ) 
		self->context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data );

	camera =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "Camera", lucCamera, True, data ) ;

	drawingObjectList = Stg_ComponentFactory_ConstructByList( 
		cf, 
		self->name, 
		"DrawingObject", 
		Stg_ComponentFactory_Unlimited, 
		lucDrawingObject, 
		True, 
		&drawingObjectCount,
		data );
	
	/* Grab a list of lights for this viewport */
	lightList = Stg_ComponentFactory_ConstructByList( 
		cf, 
		self->name, 
		"Light", 
		Stg_ComponentFactory_Unlimited, 
		lucLight, 
		False, 
		&lightCount,
		data );

	_lucViewport_Init(
			self,
			camera,
			drawingObjectList,
			drawingObjectCount,
			lightList,
			lightCount,
			Stg_ComponentFactory_GetBool( cf, self->name, "drawTitle", True ),
			Stg_ComponentFactory_GetBool( cf, self->name, "drawTime", False ),
			Stg_ComponentFactory_GetBool( cf, self->name, "compositeEachObject", False ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "nearClipPlane", 0.1 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "farClipPlane", 40.0 ) );

	Memory_Free( drawingObjectList );
        if(lightList)
		Memory_Free( lightList );
}

void _lucViewport_Build( void* camera, void* data ) { }
void _lucViewport_Initialise( void* camera, void* data ) { }
void _lucViewport_Execute( void* camera, void* data ) { }
void _lucViewport_Destroy( void* camera, void* data ) { }

void lucViewport_Draw( void* viewport, lucWindow* window, lucViewportInfo* viewportInfo, void* context ) {
	lucViewport*          self = (lucViewport*) viewport ;

	lucDebug_PrintFunctionBegin( self, 2 );

	/*Enables the lights */
	lucLight_Register_EnableAll( self->light_Register );

	lucDrawingObject_Register_DrawAll( self->drawingObject_Register, window, viewportInfo, context, self->compositeEachObject );

	lucDebug_PrintFunctionEnd( self, 2 );
}

void lucViewport_CleanUp( void* viewport, void* context ) {
	lucViewport*          self = (lucViewport*) viewport ;

	lucDebug_PrintFunctionBegin( self, 2 );

	lucDrawingObject_Register_CleanUpAll( self->drawingObject_Register, context );

	lucDebug_PrintFunctionEnd( self, 2 );
}

void lucViewport_SetNeedsToSetupFlag( void* viewport, Bool flag ) {
	lucViewport*       self = (lucViewport*) viewport ;
	
	lucDebug_PrintFunctionBegin( self, 2 );

	lucDrawingObject_Register_SetNeedsToSetupFlag( self->drawingObject_Register, flag );
	
	lucDebug_PrintFunctionEnd( self, 2 );
}

void lucViewport_Reset( void* viewport ) {
	lucViewport*       self = (lucViewport*) viewport ;

	if (self == NULL) 
		return;
	
	lucCamera_Reset( self->camera );
}

void lucViewport_Broadcast( void* viewport, int rootRank, MPI_Comm comm ) {
	lucViewport*       self      = (lucViewport*) viewport ;

	lucDebug_PrintFunctionBegin( self, 2 );

	lucCamera_Broadcast( self->camera, rootRank, comm );

	MPI_Bcast( self, 1, lucViewport_MPI_Datatype, rootRank, comm );
	
	lucDebug_PrintFunctionEnd( self, 2 );
}


#define lucViewport_TypesCount 5
void lucViewport_Create_MPI_Datatype() {
	MPI_Datatype        typeList[lucViewport_TypesCount]     = { MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT };
	int                 blocklen[lucViewport_TypesCount]     = {1, 1, 1, 1, 1};
	MPI_Aint            displacement[lucViewport_TypesCount];
	lucViewport         viewport;

	displacement[0] = GetOffsetOfMember( viewport, nearClipPlane );
	displacement[1] = GetOffsetOfMember( viewport, farClipPlane );
	displacement[2] = GetOffsetOfMember( viewport, drawTitle );
	displacement[3] = GetOffsetOfMember( viewport, drawTime );
	displacement[4] = GetOffsetOfMember( viewport, compositeEachObject );
	
	MPI_Type_struct( lucViewport_TypesCount, blocklen, displacement, typeList, &lucViewport_MPI_Datatype );
	MPI_Type_commit( & lucViewport_MPI_Datatype );
}
