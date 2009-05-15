#ifdef HAVE_VTK


#include <vtkCamera.h>

#ifdef __cplusplus
extern "C" {
#endif



#include <mpi.h>
#include <StGermain/StGermain.h>
#include <glucifer/Base/Base.h>
#include "types.h"

#include "RenderingEngineVTK.h"

#include <assert.h>
#include <string.h>

#ifndef MASTER
	#define MASTER 0
#endif

/* Textual name of this class */
const Type lucRenderingEngineVTK_Type = "lucRenderingEngineVTK";

/* Creation implementation / Virtual constructor */
lucRenderingEngineVTK* _lucRenderingEngineVTK_New( 
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
		lucRenderingEngine_RenderFunction*                 _render,
		lucRenderingEngine_ClearFunction*             	   _clear,
		lucRenderingEngine_GetPixelDataFunction*           _getPixelData,
		lucRenderingEngine_CompositeViewportFunction*      _compositeViewport,
		Name                                               name ) 
{
	lucRenderingEngineVTK*					self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(lucRenderingEngineVTK) );
	self = (lucRenderingEngineVTK*) _lucRenderingEngine_New( 
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
			_render,
			_clear,
			_getPixelData,
			_compositeViewport,
			name );
	
	return self;
}

void _lucRenderingEngineVTK_Init( 
		lucRenderingEngineVTK*                                      self  ) 
{
}

void _lucRenderingEngineVTK_Delete( void* renderingEngine ) {
	lucRenderingEngineVTK*  self = (lucRenderingEngineVTK*)renderingEngine;

	_lucRenderingEngine_Delete( self );
}

void _lucRenderingEngineVTK_Print( void* renderingEngine, Stream* stream ) {
	lucRenderingEngineVTK*  self = (lucRenderingEngineVTK*)renderingEngine;

	_lucRenderingEngine_Print( self, stream );
}

void* _lucRenderingEngineVTK_Copy( void* renderingEngine, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucRenderingEngineVTK*  self = (lucRenderingEngineVTK*)renderingEngine;
	lucRenderingEngineVTK* newRenderingEngine;

	newRenderingEngine = (lucRenderingEngineVTK*) _lucRenderingEngine_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newRenderingEngine;
}


void* _lucRenderingEngineVTK_DefaultNew( Name name ) {
	return (void*) _lucRenderingEngineVTK_New(
		sizeof(lucRenderingEngineVTK),
		lucRenderingEngineVTK_Type,
		_lucRenderingEngineVTK_Delete,
		_lucRenderingEngineVTK_Print,
		NULL,
		_lucRenderingEngineVTK_DefaultNew,
		_lucRenderingEngineVTK_Construct,
		_lucRenderingEngineVTK_Build,
		_lucRenderingEngineVTK_Initialise,
		_lucRenderingEngineVTK_Execute,
		_lucRenderingEngineVTK_Destroy,
		_lucRenderingEngineVTK_Render,
		_lucRenderingEngineVTK_Clear,
		_lucRenderingEngineVTK_GetPixelData,
		_lucRenderingEngineVTK_CompositeViewport_Stencil,
		name );
}

void _lucRenderingEngineVTK_Construct( void* renderingEngine, Stg_ComponentFactory* cf ){
	lucRenderingEngineVTK*  self = (lucRenderingEngineVTK*)renderingEngine;

	/* Construct Parent */
	_lucRenderingEngine_Construct( self, cf );
	
	_lucRenderingEngineVTK_Init( self );
}

void _lucRenderingEngineVTK_Build( void* renderingEngine, void* data ) {}
void _lucRenderingEngineVTK_Initialise( void* renderingEngine, void* data ) {}
void _lucRenderingEngineVTK_Execute( void* renderingEngine, void* data ) {}
void _lucRenderingEngineVTK_Destroy( void* renderingEngine, void* data ) {}

void lucRenderingEngineVTK_Camera_CXX( void* _context, lucViewport* viewport, vtkCamera* camera ) {
	lucCamera* myCamera = viewport->camera;

	/* Change viewer location */
	camera->SetPosition( myCamera->coord[ I_AXIS ], myCamera->coord[ J_AXIS ], myCamera->coord[ K_AXIS ] );
	camera->SetFocalPoint( myCamera->focalPoint[ I_AXIS ], myCamera->focalPoint[ J_AXIS ], myCamera->focalPoint[ K_AXIS ] );
	camera->SetViewUp( myCamera->upDirection[ I_AXIS ], myCamera->upDirection[ J_AXIS ], myCamera->upDirection[ K_AXIS ] );
}


void _lucRenderingEngineVTK_Render( void* renderingEngine, lucWindow* window, AbstractContext* context ) {
 /*	lucRenderingEngineVTK* self              = (lucRenderingEngineVTK*) renderingEngine;
	lucViewport*          viewport;
	Viewport_Index        viewport_I;
	Viewport_Index        viewportCount     = window->viewportCount;
	lucViewportInfo*      viewportInfo;

	Journal_DPrintfL( lucDebug, 2, "In func: %s for %s '%s'\n", __func__, self->type, self->name );
	Stream_Indent( lucDebug );

///TO RE INSERT

 AbstractContext* context = (AbstractContext*) _context;
	lucViewport *viewport;
	lucWindowingContextExtension* windowContextExt = (lucWindowingContextExtension*)ExtensionManager_Get( context->extensionMgr, context, lucWindowingContextExtensionHandle );
	 lucWindowingVTKContextExtension* VTKwindowContextExt = (lucWindowingVTKContextExtension*)ExtensionManager_Get( context->extensionMgr, context, lucWindowingVTKContextExtensionHandle);
	int nViewports = windowContextExt->nViewports;
	lucColour* bgColour = &windowContextExt->backgroundColour;
	vtkRenderer* renderer;
	vtkCamera* camera;
	int n;

	for ( n = 0 ; n < nViewports ; n++ ) {
		viewport = lucViewport_At( windowContextExt , n );
	
		// Create Renderer Window 
		//viewportExtension = ExtensionManager_Get( context->extensionMgr, context, lucWindowingContextExtensionHandle );
		renderer = vtkRenderer::New();
	
		//Add Renderer to Rendering Window 
		VTKwindowContextExt->renderer= renderer;
		VTKwindowContextExt->window->AddRenderer( renderer );

		//Set viewport 
		renderer->SetViewport( viewport->startx, viewport->starty, viewport->endx, viewport->endy );
		renderer->SetBackground( bgColour->red, bgColour->green, bgColour->blue );

		// Create camera 
		camera = vtkCamera::New();
		renderer->SetActiveCamera( camera );
		lucRenderingEngineVTK_Camera_CXX( context, viewport, camera );

		// Entry point for Visualisation plugins - All the VTK Actor stuff 
		lucWindowing_DrawViewport( context, viewport );
		
		
		 //DRAW TEST
	// Draw something scene 
	//TEST
	/vtkSphereSource *sphere= vtkSphereSource::New();
	sphere->SetRadius(0.5);
  sphere->SetThetaResolution(18);
  sphere->SetPhiResolution(18);
	
	// map to graphics library
  vtkPolyDataMapper *map = vtkPolyDataMapper::New();
  map->SetInput(sphere->GetOutput());

 
 	vtkActor *aSphere = vtkActor::New();
	aSphere->SetMapper(map);

  // an interactor (been created in lucWIndowingVTK/src/lucWindowingVTK.cx)
 // vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
//  iren->SetRenderWindow(VTKwindowContextExt->window);


	
	renderer->AddActor(aSphere);
	 renderer->SetBackground(0,0,1);
	 
///////	  // begin mouse interaction
  //iren->Start();
	
	//VTKwindowContextExt->interactor->Start();
	
	
	}*/
  
	
	// Render the scene 
	//VTKwindowContextExt->window->Render();
	//VTKwindowContextExt->window->MakeCurrent();
	//
///////////////////////////////////////////



/*	VTKPolygonMode(VTK_FRONT_AND_BACK, VTK_FILL);	
	VTKEnable(VTK_COLOR_MATERIAL);
	VTKColorMaterial(VTK_FRONT_AND_BACK, VTK_AMBIENT_AND_DIFFUSE);
	VTKLightModeli(VTK_LIGHT_MODEL_TWO_SIDE, VTK_TRUE);		
	VTKDrawBuffer(VTK_BACK_LEFT);

	/
	VTKEnable (VTK_BLEND);
	VTKBlendFunc (VTK_SRC_ALPHA, VTK_ONE_MINUS_SRC_ALPHA );
	VTKColorMask(VTK_TRUE,VTK_TRUE,VTK_TRUE,VTK_TRUE);
	
	VTKEnable(VTK_DEPTH_TEST);

	/
	VTKEnable(VTK_LIGHTING);
	VTKEnable(VTK_LIGHT0);*/
	
	///////////////////////////////////////////////////////

/*	lucWindow_Broadcast( window, 0, MPI_COMM_WORLD );
	lucWindow_CheckCameraFlag( window );
		
	for ( viewport_I = 0 ; viewport_I < viewportCount ; viewport_I++ ) {
		viewportInfo = &window->viewportInfoList[ viewport_I ];
		viewport = viewportInfo->viewport;

		Journal_DPrintfL( lucDebug, 2, "In loop for viewport '%s'.\n", viewport->name );
		Stream_Indent( lucDebug );
		
		 Set viewport 
		VTKViewport( viewportInfo->startx, viewportInfo->starty, viewportInfo->width, viewportInfo->height);
		VTKScissor( viewportInfo->startx, viewportInfo->starty, viewportInfo->width, viewportInfo->height);
		if ( ! viewportInfo->needsToDraw ) {
			Journal_DPrintfL( lucDebug, 2, "Viewport '%s' doesn't need to be redrawn.\n", viewport->name );
			Stream_UnIndent( lucDebug );
			continue;
		}

		//lucRenderingEngineVTK_Clear( self, window, False );
		self->_clear(self, window);

		if (context->rank == MASTER)
			lucRenderingEngineVTK_DrawTitle( self, window, viewportInfo );
			
			
		switch ( viewport->camera->stereoType ) {
			case lucMono:
				lucViewportInfo_SetOpenVTKCamera( viewportInfo );
				lucViewport_Draw( viewport, window, viewportInfo, context );
				break;
			case lucStereoToeIn: case lucStereoAsymmetric:
				VTKDrawBuffer(VTK_BACK_RIGHT);
				viewport->camera->buffer = lucRight;

				lucViewportInfo_SetOpenVTKCamera( viewportInfo );
				lucViewport_Draw( viewport, window, viewportInfo, context );

				VTKDrawBuffer(VTK_BACK_LEFT);
				viewport->camera->buffer = lucLeft;
				
				lucViewportInfo_SetOpenVTKCamera( viewportInfo );
				lucViewport_Draw( viewport, window, viewportInfo, context );
		}
		
		viewportInfo->needsToDraw = False;

		Stream_UnIndent( lucDebug );
		Journal_DPrintfL( lucDebug, 2, "Finised loop.\n" );
	}
	
	Stream_UnIndent( lucDebug );
	Journal_DPrintfL( lucDebug, 2, "Leaving func %s\n", __func__ );*/
}

void _lucRenderingEngineVTK_GetPixelData( void* renderingEngine, lucWindow* window, lucPixel* buffer ) {
/*	VTKsizei width  = window->width;
	VTKsizei height = window->height;
	
	VTKPixelStorei(VTK_PACK_ALIGNMENT,1);
	
	if ( lucWindow_HasStereoCamera( window ) ) {
		if ( window->currStereoBuffer == lucRight )
			VTKReadBuffer( VTK_FRONT_RIGHT );
		else
			VTKReadBuffer( VTK_FRONT_LEFT );
	}

	
	VTKReadPixels(0, 0, width, height, VTK_RGB, VTK_UNSIGNED_BYTE, buffer); */
}

void lucRenderingEngineVTK_DrawTitle( void* renderingEngine, lucWindow* window, lucViewportInfo* viewportInfo ) {
	lucViewport* viewport = viewportInfo->viewport;

	/* Print Title 
	if (viewport->drawTitle) {
		int stringWidth = lucStringWidth( viewport->name );

			VTKPushMatrix();
		VTKMatrixMode(VTK_PROJECTION);
		VTKLoadIdentity();
		VTKuOrtho2D((VTKfloat) 0.0, (VTKfloat) viewportInfo->width, (VTKfloat) 0.0, (VTKfloat) viewportInfo->height );
		VTKMatrixMode(VTK_MODELVIEW);
		VTKLoadIdentity();	

		lucColour_SetComplimentaryOpenVTKColour( &window->backgroundColour );

		VTKRasterPos2i( viewportInfo->width/2 - stringWidth/2, viewportInfo->height - 13 );
		lucPrintString( viewport->name );
		
		VTKPopMatrix();
	}*/
}

void _lucRenderingEngineVTK_Clear( void* renderingEngineVTK, lucWindow* window, Bool clearAll ) {
	}

Index lucRenderingEngineVTK_MapBufferIdToRank( void* renderingEngineVTK, Index bufferId, Index mergeCount ) {
}

void _lucRenderingEngineVTK_CompositeViewport_Stencil( 
		void*                                              renderingEngine, 
		lucViewportInfo*                                   viewportInfo, 
		AbstractContext*                                   context, 
		Bool                                               broadcast ) 
{

}

void lucRenderingEngineVTK_CombineToMaster( 
		void*                                              renderingEngine,
		lucViewportInfo*                                   viewportInfo,
		AbstractContext*                                   context,
		lucPixel*                                          imageBuffer, 
		float*                                             depthBuffer )
{
}

void _lucRenderingEngineVTK_CompositeViewport_Manual( 
		void*                                              renderingEngine, 
		lucViewportInfo*                                   viewportInfo, 
		AbstractContext*                                   context, 
		Bool                                               broadcast ) 
{
}

#ifdef __cplusplus
}
#endif

#endif
