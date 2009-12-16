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
** $Id: Viewport.h 628 2006-10-12 08:23:07Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __lucViewport_h__
#define __lucViewport_h__

	extern const Type lucViewport_Type;
	extern MPI_Datatype lucViewport_MPI_Datatype;

	#define __lucViewport                                         \
		__Stg_Component                                           \
		AbstractContext*				   context; 		     \
      lucCamera*                                         camera;                   \
      lucDrawingObject_Register*                         drawingObject_Register;   \
      lucLight_Register*                         	      light_Register;           \
      lucLight*                                          defaultLight;             \
		double                                             nearClipPlane;            \
		double                                             farClipPlane;             \
		double                                             scaleX;              \
		double                                             scaleY;              \
		double                                             scaleZ;              \
		Bool                                               drawTitle;                \
		Bool                                               drawTime;                 \
		Bool                                               compositeEachObject;

	struct lucViewport {__lucViewport};

	lucViewport* lucViewport_New(
		Name                                               name,
		lucCamera*                                         camera,
		lucDrawingObject**                                 drawingObjectList,
		DrawingObject_Index                                drawingObjectCount,
		lucLight**                                         lightList,
      Light_Index                                        lightCount,
      Bool                                               drawTitle,
		Bool                                               drawTime,
		Bool                                               compositeEachObject,
		double                                             nearClipPlane,
		double                                             farClipPlane,
		double                                             scaleX,
		double                                             scaleY,
		double                                             scaleZ );

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define LUCVIEWPORT_DEFARGS \
                STG_COMPONENT_DEFARGS

	#define LUCVIEWPORT_PASSARGS \
                STG_COMPONENT_PASSARGS

	lucViewport* _lucViewport_New(  LUCVIEWPORT_DEFARGS  );

	void _lucViewport_Delete( void* viewport ) ;
	void _lucViewport_Print( void* viewport, Stream* stream ) ;
	void* _lucViewport_Copy( void* viewport, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) ;

	void* _lucViewport_DefaultNew( Name name ) ;
   void _lucViewport_AssignFromXML( void* viewport, Stg_ComponentFactory* cf, void* data ) ;
	void _lucViewport_Build( void* viewport, void* data );
	void _lucViewport_Initialise( void* viewport, void* data );
	void _lucViewport_Execute( void* viewport, void* data );
	void _lucViewport_Destroy( void* viewport, void* data );

	/* +++ Public Functions +++ */
	void lucViewport_SetNeedsToSetupFlag( void* viewport, Bool flag ) ;
	void lucViewport_Draw( void* viewport, lucWindow* window, lucViewportInfo* viewportInfo, void* context ) ;
	void lucViewport_CleanUp( void* viewport, void* context ) ;
	void lucViewport_Reset( void* viewport ) ;
	void lucViewport_Broadcast( void* viewport, int rootRank, MPI_Comm comm ) ;
	void lucViewport_Create_MPI_Datatype() ;

#endif

