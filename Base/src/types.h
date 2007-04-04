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
** $Id: types.h 540 2006-04-27 05:26:37Z CecileDuboz $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __lucBase_types_h__
#define __lucBase_types_h__

	typedef struct lucColourMap                lucColourMap;
	typedef struct lucColour                   lucColour;
	
	typedef struct lucCamera                   lucCamera;

	typedef struct lucWindow                   lucWindow;
	
	typedef struct lucWindowInteraction_Register        lucWindowInteraction_Register;	
	typedef struct lucWindowInteraction        lucWindowInteraction;	
	
	typedef struct lucRenderingEngine          lucRenderingEngine;

	typedef struct lucViewport                 lucViewport;
	typedef struct lucViewportInfo             lucViewportInfo;

	typedef struct lucOutputFormat_Register    lucOutputFormat_Register;
	typedef struct lucOutputFormat             lucOutputFormat;	
	
	typedef struct lucInputFormat_Register     lucInputFormat_Register;
	typedef struct lucInputFormat             lucInputFormat;

	typedef struct lucDrawingObject_Register   lucDrawingObject_Register;
	typedef struct lucDrawingObject            lucDrawingObject;
	
	typedef struct lucLight_Register           lucLight_Register;
  typedef struct lucLight                    lucLight;
	
	/* types for clarity */
	typedef Index Colour_Index;
	typedef Index Camera_Index;
	typedef Index Viewport_Index;
	typedef Index DrawingObject_Index;
	typedef Index Pixel_Index;
	typedef Index Window_Index;
	typedef Index OutputFormat_Index;
	typedef Index InputFormat_Index;
	typedef Index WindowInteraction_Index;
	typedef Index Light_Index;
	
	typedef unsigned char lucPixel[3];
	typedef unsigned char lucAlphaPixel[4];

	typedef enum {
		lucMono,
		lucStereoToeIn,
		lucStereoAsymmetric
	} lucStereoType;
		
	typedef enum {
		lucLeft,
		lucRight
	} lucStereoBuffer;

	typedef enum {
		lucButtonPress = 4,
		lucButtonRelease = 5
	} lucMouseState;

	typedef enum {
		lucLeftButton   = 1,
		lucMiddleButton = 2,
		lucRightButton  = 3,
		lucWheelUp      = 4,
		lucWheelDown    = 5
	} lucMouseButton;



#endif
