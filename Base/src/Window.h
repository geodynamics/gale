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
** $Id: Window.h 754 2008-01-11 05:41:53Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __lucWindow_h__
#define __lucWindow_h__

	typedef void (lucWindow_DisplayFunction) ( void* object );
	typedef int (lucWindow_EventsWaitingFunction) ( void* object );
	typedef Bool (lucWindow_EventProcessorFunction) ( void* object );
	typedef void (lucWindow_ResizeFunction) ( void* object );

	extern const Type lucWindow_Type;

	extern MPI_Datatype lucWindow_MPI_Datatype;
	extern MPI_Datatype lucViewportInfo_MPI_Datatype;

	#define __lucWindow														\
		__Stg_Component                                                     \
		/* Virtual Functions */ 											\
		lucWindow_DisplayFunction*			_displayWindow;					\
		lucWindow_EventsWaitingFunction*	_eventsWaiting;					\
		lucWindow_EventProcessorFunction*  	_eventProcessor;				\
		lucWindow_ResizeFunction*  	        _resizeWindow;				    \
		\
		/* Other Info */ \
		lucRenderingEngine*                	renderingEngine;				\
		lucViewportInfo*              		viewportInfoList;				\
		Viewport_Index                		viewportCount;					\
        Viewport_Index                      verticalCount;                  \
		lucOutputFormat_Register*     		outputFormat_Register;			\
		lucWindowInteraction_Register*		windowInteraction_Register;		\
		lucWindowInteraction*				defaultWindowInteraction;		\
		Pixel_Index							width;							\
		Pixel_Index							height;							\
        Bool                                resized;                        \
		Bool								interactive;					\
		Bool								continuous;						\
		lucColour							backgroundColour;				\
		lucStereoBuffer						currStereoBuffer;				\
		Bool								quitEventLoop;					\
		Bool								toggleApplicationQuit;			\
		AbstractContext*					context;						\
		Bool								isMaster;						\
		Bool								isTimedOut;         			\
		double								maxIdleTime;        			\
		double								idleTime;						\
		Pixel_Index							startx;							\
		Pixel_Index							starty;							\
		char*								title;            	  			\
			
	struct lucWindow {__lucWindow};

	lucWindow* _lucWindow_New(
		SizeT                                    	sizeOfSelf,
		Type                                     	type,
		Stg_Class_DeleteFunction*                	_delete,
		Stg_Class_PrintFunction*                 	_print,
		Stg_Class_CopyFunction*                  	_copy, 
		Stg_Component_DefaultConstructorFunction*	_defaultConstructor,
		Stg_Component_ConstructFunction*            _construct,
		Stg_Component_BuildFunction*				_build,
		Stg_Component_InitialiseFunction*			_initialise,
		Stg_Component_ExecuteFunction*				_execute,
		Stg_Component_DestroyFunction*				_destroy,
		lucWindow_DisplayFunction*					_displayWindow,	
		lucWindow_EventsWaitingFunction*			_eventsWaiting,	
		lucWindow_EventProcessorFunction*			_eventProcessor,	
		lucWindow_ResizeFunction*					_resizeWindow,	
		Name										name );

	void _lucWindow_Delete( void* window ) ;
	void _lucWindow_Print( void* window, Stream* stream ) ;
	void* _lucWindow_Copy( void* window, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) ;

	/* Stg_Component Virtual Functions */
	void* _lucWindow_DefaultNew( Name name ) ;
	void _lucWindow_Construct( void* window, Stg_ComponentFactory* cf, void* data ) ;
	void _lucWindow_Build( void* window, void* data ) ;
	void _lucWindow_Initialise( void* window, void* data ) ;
	void _lucWindow_Execute( void* window, void* data ) ;
	void _lucWindow_Destroy( void* window, void* data ) ;

	/* Window Virtual Functions */
	void lucWindow_Display( void* window );		
	int lucWindow_EventsWaiting( void* window);	
	Bool lucWindow_EventProcessor( void* window );	
    void lucWindow_Resize( void* window);

	/* +++ Public Functions +++ */
	void lucWindow_Dump( void* window, AbstractContext* context ) ;
	void lucWindow_CleanUp( void* window, void* context ) ;

	/** Functions for interactivity */
	void lucWindow_MouseMotion( void* window, lucMouseButton button, Pixel_Index xpos, Pixel_Index ypos) ;
	void lucWindow_MouseClick( void* window, lucMouseButton button, lucMouseState state, Pixel_Index xpos, Pixel_Index ypos) ;
	void lucWindow_KeyboardEvent( void* window, char key, Pixel_Index xpos, Pixel_Index ypos) ;

	Bool lucWindow_HasStereoCamera( lucWindow* self ) ;
	void lucWindow_CheckCameraFlag( void* window ) ;
	void lucWindow_CheckLightFlag( void* window ) ;


	void lucWindow_SetViewportNeedsToSetupFlag( void* window, Bool flag ) ;
	void lucWindow_SetViewportNeedsToDrawFlag( void* window, Bool flag ) ;

	lucViewportInfo* lucWindow_GetViewportInfoByPixel( void* window, Pixel_Index xPixel, Pixel_Index yPixel ) ;
	lucViewport* lucWindow_GetViewportByPixel( void* window, Pixel_Index xPixel, Pixel_Index yPixel ) ;
	void lucWindow_Broadcast( void* window, int rootRank, MPI_Comm comm ) ;

	void lucViewportInfo_Broadcast( lucViewportInfo* self, int rootRank, MPI_Comm comm ) ;
	Bool lucViewportInfo_IsPixelInside( lucViewportInfo* viewportInfo, Pixel_Index xPixel, Pixel_Index yPixel ) ;

	#define lucViewportInfo_AspectRatio( viewportInfo ) \
		( (double) (viewportInfo)->width / (double) (viewportInfo)->height )

	lucViewportInfo* lucWindow_ConstructViewportInfoList( 
		lucWindow* self, 
		Stg_ComponentFactory* cf, 
		Pixel_Index width, 
		Pixel_Index height, 
		Viewport_Index* viewportCount,
		void* data ) ;

	void lucWindow_InteractionHelpMessage( void* window, Stream* stream ) ;

	void lucWindow_Create_MPI_Datatype() ;
	void lucViewportInfo_Create_MPI_Datatype() ;

	void lucWindow_ChangeInteractiveMode( void* window );
	void lucWindow_ChangeContinuousMode( void* window );

	void lucWindow_QuitEventLoop( void* window ) ;
	void lucWindow_ToggleApplicationQuit( void* window ) ;
	Bool lucWindow_SetSize( void* window, Pixel_Index newWidth, Pixel_Index newHeight ) ;
	void lucWindow_IdleReset(void *window);
	void lucWindow_IdleCheck(void *window);

#endif
