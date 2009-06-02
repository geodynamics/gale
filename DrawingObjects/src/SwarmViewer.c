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
** $Id: SwarmViewer.c 791 2008-09-01 02:09:06Z JulianGiordani $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#ifdef GLUCIFER_USE_PICELLERATOR
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#endif

#include <glucifer/Base/Base.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>

#include "types.h"

#include "OpenGLDrawingObject.h"
#include "SwarmViewerBase.h"
#include "SwarmViewer.h"

#include <assert.h>
#ifdef HAVE_OPENGL_FRAMEWORK
	#include <OpenGL/gl.h>
	#include <OpenGL/glu.h>
#else
	#include <gl.h>
	#include <glu.h>
#endif
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucSwarmViewer_Type = "lucSwarmViewer";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucSwarmViewer* _lucSwarmViewer_New( 
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
	lucDrawingObject_SetupFunction*                    _setup,
	lucDrawingObject_DrawFunction*                     _draw,
	lucDrawingObject_CleanUpFunction*                  _cleanUp,
	lucOpenGLDrawingObject_BuildDisplayListFunction*   _buildDisplayList,
	lucSwarmViewerBase_PlotParticleFunction*           _plotParticle,
	lucSwarmViewerBase_SetParticleColourFunction*      _setParticleColour,
	Name                                               name ) 
{
	lucSwarmViewer*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucSwarmViewer) );
	self = (lucSwarmViewer*) _lucSwarmViewerBase_New( 
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
		_setup,
		_draw,
		_cleanUp,
		_buildDisplayList,
		_plotParticle,
		_setParticleColour,
		name );

	return self;
}

void _lucSwarmViewer_Init( 
	lucSwarmViewer*                                              self,
	float                                                        pointSize, 
	Bool													     pointSmoothing )
{
	self->pointSize           = pointSize;
	self->pointSmoothing      = pointSmoothing;
}

void _lucSwarmViewer_Delete( void* drawingObject ) {
	lucSwarmViewer*  self = (lucSwarmViewer*)drawingObject;

	_lucOpenGLDrawingObject_Delete( self );
}

void _lucSwarmViewer_Print( void* drawingObject, Stream* stream ) {
	lucSwarmViewer*  self = (lucSwarmViewer*)drawingObject;

	_lucOpenGLDrawingObject_Print( self, stream );
}

void* _lucSwarmViewer_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucSwarmViewer*  self = (lucSwarmViewer*)drawingObject;
	lucSwarmViewer* newDrawingObject;

	newDrawingObject = _lucOpenGLDrawingObject_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newDrawingObject;
}


void* _lucSwarmViewer_DefaultNew( Name name ) {
	return (void*) _lucSwarmViewer_New(
		sizeof(lucSwarmViewer),
		lucSwarmViewer_Type,
		_lucSwarmViewer_Delete,
		_lucSwarmViewer_Print,
		NULL,
		_lucSwarmViewer_DefaultNew,
		_lucSwarmViewer_Construct,
		_lucSwarmViewer_Build,
		_lucSwarmViewer_Initialise,
		_lucSwarmViewer_Execute,
		_lucSwarmViewer_Destroy,
		_lucSwarmViewer_Setup,
		_lucSwarmViewer_Draw,
		_lucSwarmViewer_CleanUp,
		_lucSwarmViewer_BuildDisplayList,
		_lucSwarmViewer_PlotParticle,
		_lucSwarmViewerBase_SetParticleColourDefault,
		name );
}

void _lucSwarmViewer_Construct( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucSwarmViewer*         self = (lucSwarmViewer*)drawingObject;

	/* Construct Parent */
	_lucSwarmViewerBase_Construct( self, cf, data );
	
	_lucSwarmViewer_Init( 
		self, 
		(float) Stg_ComponentFactory_GetDouble( cf, self->name, "pointSize",      1.0 ),
		(Bool ) Stg_ComponentFactory_GetBool(   cf, self->name, "pointSmoothing", 1 )
	  );
}

void _lucSwarmViewer_Build( void* drawingObject, void* data ) {
	_lucSwarmViewerBase_Build( drawingObject, data );
}

void _lucSwarmViewer_Initialise( void* drawingObject, void* data ) {
	_lucSwarmViewerBase_Initialise( drawingObject, data );
}


void _lucSwarmViewer_Execute( void* drawingObject, void* data ) {}
void _lucSwarmViewer_Destroy( void* drawingObject, void* data ) {}

void _lucSwarmViewer_Setup( void* drawingObject, void* _context ) {
	lucSwarmViewer*          self                = (lucSwarmViewer*)drawingObject;
	
	_lucSwarmViewerBase_Setup( self, _context );
}

void _lucSwarmViewer_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	lucSwarmViewer*          self          = (lucSwarmViewer*)drawingObject;
	lucCamera*               camera        = viewportInfo->viewport->camera;
	XYZ                      normal;

	/* Draw the particles with a normal facing the camera 
	 * these lines have to be here because the camera can move after building the display list */
	StGermain_VectorSubtraction( normal, camera->coord, camera->focalPoint, 3 );
	glNormal3dv(normal);

	_lucSwarmViewerBase_Draw( self, window, viewportInfo, _context ); 
}


void _lucSwarmViewer_CleanUp( void* drawingObject, void* context ) {
	lucSwarmViewer*          self          = (lucSwarmViewer*)drawingObject;
	
	_lucSwarmViewerBase_CleanUp( self, context );
}

void _lucSwarmViewer_BuildDisplayList( void* drawingObject, void* _context ) {
	lucSwarmViewer*          self                = (lucSwarmViewer*)drawingObject;

	
	/* lighting of small particle objects like particle dots seems to do more harm than good - 
	   the gl lighting system doesn't seem to deal with lighting such tiny objects well. Lighting
	   of particles doesn't really help the user anyway, so I've disabled it for now.
	   PatrickSunter - 8 Jun 2006 */
			
	glDisable(GL_LIGHTING);
		
	if(self->pointSmoothing) { 
        /* Round, smooth points */
		glEnable(GL_POINT_SMOOTH);
		/* Point smoothing will not work correctly with depth testing enabled*/
	    glDepthFunc(GL_ALWAYS);
    }
	else 
		glDisable(GL_POINT_SMOOTH);
		
	glPointSize( self->pointSize );
	
	glBegin( GL_POINTS );
	_lucSwarmViewerBase_BuildDisplayList( self, _context );
	glEnd( );

	/* Put back lighting / smoothing settings to low-impact options */
	glEnable(GL_LIGHTING);	
		
	if(self->pointSmoothing) {
		glDisable(GL_POINT_SMOOTH);
	    glDepthFunc(GL_LESS);
	}
}


void _lucSwarmViewer_PlotParticle( void* drawingObject, void* _context, Particle_Index lParticle_I ) {
	lucSwarmViewer*          self                = (lucSwarmViewer*)drawingObject;
	DomainContext*   context             = (DomainContext*) _context;
	GlobalParticle*          particle            = (GlobalParticle*)Swarm_ParticleAt( self->swarm, lParticle_I );
	double*                  coord               = particle->coord;
	float                    offset              = 0.001; 

	if (context->dim == 2)
		glVertex3f( (float)coord[0], (float)coord[1], offset);
	else   
		glVertex3dv( coord );
}

