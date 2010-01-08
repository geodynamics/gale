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
** $Id: Camera.h 628 2006-10-12 08:23:07Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/



#ifndef __lucCamera_h__
#define __lucCamera_h__

	extern const Type lucCamera_Type;
	extern MPI_Datatype lucCamera_MPI_Datatype;

	#define __lucCamera \
		__Stg_Component \
		AbstractContext*	context; \
		lucCamera*			originalCamera; \
		FieldVariable*		centreFieldVariable; \
		Coord					coord; \
		Coord					focalPoint; \
		Coord					rotationCentre; \
		XYZ					upDirection; \
		double				focalLength; \
		double				aperture; \
		double				eyeSeparation; \
		lucStereoBuffer	buffer; \
		lucStereoType		stereoType; \
		Bool					needsToDraw; 

	struct lucCamera {__lucCamera};

	/** Constructors */
	lucCamera* lucCamera_New( 
		Name					name,
		Coord					coord, 
		Coord					focalPoint,
		Coord					rotationCentre,
		XYZ					upDirection,
		double				focalLength,
		double				aperture,
		double				eyeSeparation,
		lucStereoType		stereoType,
		FieldVariable*		centreFieldVariable );

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define LUCCAMERA_DEFARGS \
                STG_COMPONENT_DEFARGS

	#define LUCCAMERA_PASSARGS \
                STG_COMPONENT_PASSARGS

	lucCamera* _lucCamera_New(  LUCCAMERA_DEFARGS  );

	/** Virtual Functions */
	void _lucCamera_Init(                                                                                                  
		void*          camera,                                                                                             
		Coord          coord,                                                                                              
		Coord          focalPoint,                                                                                         
		Coord          rotationCentre,                                                                                     
		XYZ            upDirection,                                                                                        
		double         focalLength,                                                                                        
		double         aperture,                                                                                           
		double         eyeSeparation,                                                                                      
		lucStereoType  stereoType,                                                                                         
   FieldVariable* centreFieldVariable );
      
	void _lucCamera_Delete( void* camera );

	void _lucCamera_Print( void* camera, Stream* stream );

	void _lucCamera_Copy( void* camera, void* dest );

	void* _lucCamera_DefaultNew( Name name );

	void _lucCamera_AssignFromXML( void* camera, Stg_ComponentFactory* cf, void* data );

	void _lucCamera_Build( void* camera, void* data );

	void _lucCamera_Initialise( void* camera, void* data );

	void _lucCamera_Execute( void* camera, void* data );

	void _lucCamera_Destroy( void* camera, void* data );

	/** Public Functions */
	void lucCamera_Zoom( void* camera, double zoomFactor );

	void lucCamera_RotateAroundUpDirection( void* camera, double deltaTheta );

	void lucCamera_RotateTowardsUpDirection( void* camera, double deltaTheta );

	void lucCamera_GetFocusDirection( void* camera, XYZ focusDirection );

	void lucCamera_GetLeftDirection( void* camera, XYZ leftDirection );

	void lucCamera_Reset( void* camera );

	void lucCamera_SetOriginal( void* camera );

	void lucCamera_Broadcast( void* camera, int rootRank, MPI_Comm comm );

	void lucCamera_SwapStereoBuffer( void* camera ) ;
	void lucCamera_CurrentEyePosition( void* camera, Coord currEyePos ) ;
	void lucCamera_CentreFromFieldVariable( void* camera ) ;
	void lucCamera_SetNeedsToDraw( void * camera );
	void lucCamera_ChangeFocalPoint( void* camera, Pixel_Index startx, Pixel_Index starty, Pixel_Index posx, Pixel_Index posy );
	void lucCamera_Pickle( void* camera, Stream* stream ) ;
	
	void lucCamera_Create_MPI_Datatype() ;
#endif

