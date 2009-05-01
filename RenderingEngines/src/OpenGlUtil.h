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
** $Id: OpenGlUtil.h 768 2008-04-21 03:20:07Z JohnMansour $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __lucRenderingEngines_OpenGlUtil_h__
#define __lucRenderingEngines_OpenGlUtil_h__

#ifdef HAVE_GL2PS
/* call to gl2pText is required for text output using vector formats, 
 * as no text is stored in the GL feedback buffer */  
	#define lucPrintString( A )       gl2psText( A, "Times-Roman", 16); glCallLists( strlen( A ), GL_UNSIGNED_BYTE, A )
#else
	#define lucPrintString( A )       glListBase(2000); glCallLists( strlen( A ), GL_UNSIGNED_BYTE, A ); glListBase(0)
#endif

/* Big hack !! - Assumes that font width is 10 pixels - this is hardly ever true */
#define lucStringWidth( A )     (strlen( A ) * 10) 

/* In OpenGL you cannot set the raster position to be outside the viewport - 
 * but this is a hack which OpenGL allows (and advertises) - That you can get the raster position to a legal value
 * and then move the raster position by calling glBitmap with a NULL bitmap */
#define lucMoveRaster( deltaX, deltaY ) \
	glBitmap( 0,0,0.0,0.0, (float)(deltaX), (float)(deltaY), NULL )

void lucColour_SetOpenGLColour( lucColour* colour ) ;
void lucColour_SetComplimentaryOpenGLColour( lucColour* colour ) ;
void lucColourMap_SetOpenGLColourFromValue( lucColourMap* cmap, double value ) ;
void lucColourMap_SetOpenGLColourFromValue_ExplicitOpacity( lucColourMap* cmap, double value, float opacity ) ;
void lucColourMap_SetOpenGLColourFromRGB( double red, double green, double blue);
void lucColourMap_SetOpenGLColourFromRGB_ExplicitOpacity( double red, double green, double blue, float opacity );
void lucViewportInfo_GetCoordFromPixel( lucViewportInfo* viewportInfo, Pixel_Index xPos, Pixel_Index yPos, Coord coord ) ;
void lucViewportInfo_SetOpenGLCamera( lucViewportInfo* viewportInfo ) ;

void luc_OpenGlSquare( Dimension_Index dim, double* pos, double* normal, double* plane, double width) ;
void luc_OpenGlFatSquare( Dimension_Index dim, double* pos, double* normal, double* plane, double width, double thickness) ;
void luc_OpenGlCircle( Dimension_Index dim, double* pos, double* normal, double radius, Index nSides) ;
void luc_DrawVector( Dimension_Index dim , double* pos, double* vector, double scale, double headSize ) ;
void luc_DrawRod( Dimension_Index dim , double* pos, double* vector, double scale ) ;

void luc_OpenGlDebug_PrintAll( Stream* stream, int debugLevel);
void luc_OpenGlDebug( Stream* stream, int mode, int debugLevel);

#endif
