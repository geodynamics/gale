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

lucViewport* _lucViewport_New(  LUCVIEWPORT_DEFARGS  )
{
	lucViewport*    self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucViewport) );
	/* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
	/* This means that any values of these parameters that are passed into this function are not passed onto the parent function
	   and so should be set to ZERO in any children of this class. */
	nameAllocationType = NON_GLOBAL;

	self = (lucViewport*) _Stg_Component_New(  STG_COMPONENT_PASSARGS  );

	return self;
}

void _lucViewport_Init( 
		lucViewport*                                       self, 
		lucCamera*                                         camera, 
		lucDrawingObject**                                 drawingObjectList, 
		DrawingObject_Index                                drawingObjectCount,
      lucLight**          				                     lightList,
      Light_Index                                        lightCount,
		Bool                                               drawTitle,
		Bool                                               drawTime,
		Bool                                               compositeEachObject,
		double                                             nearClipPlane,
		double                                             farClipPlane,
		double                                             scaleX,
		double                                             scaleY,
		double                                             scaleZ)
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
	self->scaleX                   = scaleX;
	self->scaleY                   = scaleY;
	self->scaleZ                   = scaleZ;
	self->compositeEachObject      = compositeEachObject;

	self->drawingObject_Register = lucDrawingObject_Register_New();

	for ( object_I = 0 ; object_I < drawingObjectCount ; object_I++ )
		lucDrawingObject_Register_Add( self->drawingObject_Register, drawingObjectList[ object_I ] );
		
	/* Setup light register stuff*/	
	self->light_Register = lucLight_Register_New();

	for ( light_I = 0 ; light_I < lightCount ; light_I++ )
		lucLight_Register_Add( self->light_Register, lightList[ light_I ] );

   if(lightCount == 0) {
      self->defaultLight = lucLight_New( "defaultLight", 0, GL_LIGHT_MODEL_TWO_SIDE,  GL_AMBIENT_AND_DIFFUSE, 
                                          lightPosition, lmodel_ambient, spotCutOff, spotDirection);
      lucLight_Register_Add( self->light_Register, self->defaultLight );
   }
}

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
		double                                             farClipPlane,
		double                                             scaleX,
		double                                             scaleY,
		double                                             scaleZ)
{
	lucViewport* self = _lucViewport_DefaultNew( name );

	_lucViewport_Init( self, camera, drawingObjectList, drawingObjectCount, lightList, lightCount, drawTitle, drawTime, compositeEachObject, nearClipPlane, farClipPlane, scaleX, scaleY, scaleZ);

	return self;
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
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof( lucViewport );
	Type                                                      type = lucViewport_Type;
	Stg_Class_DeleteFunction*                              _delete = _lucViewport_Delete;
	Stg_Class_PrintFunction*                                _print = _lucViewport_Print;
	Stg_Class_CopyFunction*                                  _copy = _lucViewport_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _lucViewport_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _lucViewport_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _lucViewport_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _lucViewport_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _lucViewport_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _lucViewport_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return _lucViewport_New(  LUCVIEWPORT_PASSARGS  );
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
			Stg_ComponentFactory_GetDouble( cf, self->name, "nearClipPlane", camera->focalLength / 10.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "farClipPlane", camera->focalLength * 10.0 ), 
			Stg_ComponentFactory_GetDouble( cf, self->name, "scaleX", 1.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "scaleY", 1.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "scaleZ", 1.0 ) );

	Memory_Free( drawingObjectList );
        if(lightList)
		Memory_Free( lightList );
}

void _lucViewport_Build( void* viewport, void* data ) { }
void _lucViewport_Initialise( void* viewport, void* data ) {}
void _lucViewport_Execute( void* viewport, void* data ) { }
void _lucViewport_Destroy( void* viewport, void* data ) { }

void lucViewport_Draw( void* viewport, lucWindow* window, lucViewportInfo* viewportInfo, void* context ) {
	lucViewport*          self = (lucViewport*) viewport ;

	lucDebug_PrintFunctionBegin( self, 2 );

	/* Enables the lights */
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


