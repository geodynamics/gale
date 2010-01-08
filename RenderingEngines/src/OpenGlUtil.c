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
** $Id: OpenGlUtil.c 791 2008-09-01 02:09:06Z JulianGiordani $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <glucifer/Base/Base.h>

#include "types.h"
#include "OpenGlUtil.h"
#include <gl.h>
#include <glu.h>
#include <string.h>

#include  "font.h"

unsigned int fontbase = -1, fontcharset = FONT_DEFAULT, texture;

void lucViewport2d(Bool enabled, lucViewportInfo* viewportInfo)
{
	if (enabled)
	{
		/* Set up 2D Viewer the size of the viewport */
		glDisable( GL_DEPTH_TEST );
		glPushMatrix();
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho((GLfloat) 0.0, (GLfloat) viewportInfo->width, (GLfloat) viewportInfo->height, (GLfloat) 0.0, -1.0f,1.0f);
		
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		
		/* Disable lighting because we don't want a 3D effect */
		glDisable(GL_LIGHTING);
		/* Disable line smoothing in 2d mode */
		glDisable(GL_LINE_SMOOTH);
	}
	else
	{
		/* Put back settings */
		glEnable(GL_LIGHTING);
		glEnable( GL_DEPTH_TEST );
		glEnable(GL_LINE_SMOOTH);
		glPopMatrix();
		
		/*Set back the viewport to what it should be to render any other object */
		/* If this is not done, than any object displayed after the colour bar will not appear,*/
		/* because the projection matrix and lookAt point have been altered */
		lucViewportInfo_SetOpenGLCamera( viewportInfo );
	}
}

void lucPrintString(const char* str)
{
	if (fontbase < 0)								/* Load font if not yet done */
		lucSetupRasterFont();
	
	if (fontcharset > FONT_LARGE || fontcharset < FONT_FIXED)		/* Character set valid? */
		fontcharset = FONT_FIXED;	

   	glDisable( GL_CULL_FACE );
	   glDisable( GL_LIGHTING );

    glActiveTexture(GL_TEXTURE0);
	glEnable(GL_TEXTURE_2D);				 					/* Enable Texture Mapping */
	glBindTexture(GL_TEXTURE_2D, texture);
	glListBase(fontbase - 32 + (96 * fontcharset));		/* Choose the font and charset */
	glCallLists(strlen(str),GL_UNSIGNED_BYTE, str);		/* Display */
	glDisable(GL_TEXTURE_2D);									/* Disable Texture Mapping */

	#ifdef HAVE_GL2PS
	/* call to gl2pText is required for text output using vector formats, 
	 * as no text is stored in the GL feedback buffer */  
		gl2psText( A, "Times-Roman", 16);
	#endif
}

void lucPrintf(int x, int y, const char *fmt, ...)
{
    char    text[512];
    va_list ap;                 /* Pointer to arguments list */
    if (fmt == NULL) return;    /* No format string */
    va_start(ap, fmt);          /* Parse format string for variables */
    vsprintf(text, fmt, ap);    /* Convert symbols */
    va_end(ap);

    lucPrint(x, y, text);       /* Print result string */
}

void lucPrint(int x, int y, const char *str)
{
	glPushMatrix();
	glLoadIdentity();
	glTranslated(x, y, 0);
	lucPrintString(str);
	glPopMatrix();
}

void lucPrint3d(double x, double y, double z, const char *str)
{
   /* Calculate projected screen coords in viewport */
   GLdouble modelMatrix[16];
   GLdouble projMatrix[16];
   GLint    viewportArray[4];
   double  xPos, yPos, depth;

   glGetDoublev( GL_MODELVIEW_MATRIX, modelMatrix );
   glGetDoublev( GL_PROJECTION_MATRIX, projMatrix );
   glGetIntegerv( GL_VIEWPORT, viewportArray );

   gluProject(x, y, z, 
      modelMatrix, projMatrix, viewportArray,
      &xPos, &yPos, &depth	);

   /* Switch to ortho view with 1 unit = 1 pixel and print using calculated screen coords */
   glMatrixMode(GL_PROJECTION);
   glPushMatrix();
   glLoadIdentity();

   glOrtho((GLfloat) 0.0, viewportArray[2], viewportArray[3], 0.0, -1.0f,1.0f);

   glMatrixMode(GL_MODELVIEW);
   glPushMatrix();
   glLoadIdentity();   

   glDisable( GL_DEPTH_TEST );

   /* Print at calculated position, compensating for viewport offset */
   int xs, ys;
   xs = (int)(xPos - viewportArray[0]);
   ys = (int)(viewportArray[3] - yPos - viewportArray[1]);
   if (depth <= 1.0 && xs > 0 && ys > 0 && xs < viewportArray[2] && ys < viewportArray[3]) 
      lucPrint(xs, ys, str);  /* Print if in view */

   /* Restore state */
   glPopMatrix();
   glMatrixMode(GL_PROJECTION);
   glPopMatrix();
   glMatrixMode(GL_MODELVIEW);

   glEnable( GL_DEPTH_TEST );
}

void lucSetFontCharset(int charset)
{
	fontcharset = charset;
}

/* String width calc */ 
int lucStringWidth(const char *string)
{
	/* Sum character widths in string */
	int i, len = 0, slen = strlen(string);
	for (i = 0; i < slen; i++)
		len += font_charwidths[string[i]-32 + (96 * fontcharset)];

	return len + slen;	/* Additional pixel of spacing for each character */
}

void lucSetupRasterFont() {
	/* Load font bitmaps and Convert To Textures */
	int i, j;
 	unsigned char pixel_data[IMAGE_HEIGHT * IMAGE_WIDTH * IMAGE_BYTES_PER_PIXEL];
   	unsigned char fontdata[IMAGE_HEIGHT][IMAGE_WIDTH];	/* font texture data */
	
	/* Get font pixels from source data - only need alpha channel */
	IMAGE_RUN_LENGTH_DECODE(pixel_data, IMAGE_RLE_PIXEL_DATA, IMAGE_WIDTH * IMAGE_HEIGHT, IMAGE_BYTES_PER_PIXEL);
	for (i = 0; i < IMAGE_HEIGHT; i++)
		for (j = 0; j < IMAGE_WIDTH; j++)
				fontdata[ i ][ j ] = pixel_data[ IMAGE_BYTES_PER_PIXEL * (IMAGE_WIDTH * i + j) + 3 ];
	
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);

	/* create and bind texture */
    glActiveTexture(GL_TEXTURE0);
	glGenTextures(1, &texture);	
	glBindTexture(GL_TEXTURE_2D, texture);
	glEnable(GL_COLOR_MATERIAL);
	/* use linear filtering */
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	/* generate the texture from bitmap alpha data */
	glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA, IMAGE_WIDTH, IMAGE_HEIGHT, 0, GL_ALPHA, GL_UNSIGNED_BYTE, fontdata);
	fontbase = glGenLists(GLYPHS);

	/* Build font display lists */
	lucBuildFont(16, 16, 0, 384);		/* 16x16 glyphs, 16 columns - glyphs 0-383 = first 4 fonts */
	lucBuildFont(18, 14, 384, 512);	/* 18x18 glyphs, 14 columns - glyphs 384-511 = last font */
}

void lucBuildFont(int glyphsize, int columns, int startidx, int stopidx)
{
	/* Build font display lists */
	int i;
	static float yoffset;
	float divX = IMAGE_WIDTH / (float)glyphsize;
	float divY = IMAGE_HEIGHT / (float)glyphsize;
	float glyphX = 1 / divX;	/* Width & height of a glyph in texture coords */
	float glyphY = 1 / divY;
	GLfloat cx, cy;         /* the character coordinates in our texture */
    if (startidx == 0) yoffset = 0;
	glBindTexture(GL_TEXTURE_2D, texture);
	for (i = 0; i < (stopidx - startidx); i++)
	{
		cx = (float) (i % columns) / divX;
		cy = yoffset + (float) (i / columns) / divY;
		glNewList(fontbase + startidx + i, GL_COMPILE);
			glBegin(GL_QUADS);
				glTexCoord2f(cx, cy + glyphY);
				glVertex2i(0, glyphsize);
				glTexCoord2f(cx + glyphX, cy + glyphY);
				glVertex2i(glyphsize, glyphsize);
				glTexCoord2f(cx + glyphX, cy);
				glVertex2i(glyphsize, 0);
				glTexCoord2f(cx, cy);
				glVertex2i(0, 0);
			glEnd();
			/* Shift right width of character + 1 */
			glTranslated(font_charwidths[startidx + i]+1, 0, 0);
		glEndList();
	}
	/* Save vertical offset to resume from */
	yoffset = cy + glyphY;
}

void lucDeleteFont()
{
	/* Delete fonts */
    if (fontbase >= 0) glDeleteLists(fontbase, GLYPHS);
    fontbase = -1;
    if (texture) glDeleteTextures(1, &texture);
    texture = 0;
}

void lucColour_SetOpenGLColour( lucColour* colour ) {
	glColor4f(
			colour->red,
			colour->green,
			colour->blue,
			colour->opacity );
}
void lucColour_SetComplimentaryOpenGLColour( lucColour* colour ) {
	glColor4f(
			1.0 - colour->red,
			1.0 - colour->green,
			1.0 - colour->blue,
			colour->opacity );
}

void lucColourMap_SetOpenGLColourFromValue( lucColourMap* cmap, double value ) {
	lucColour colour;

	lucColourMap_GetColourFromValue( cmap, value, &colour);
	lucColour_SetOpenGLColour( &colour );
}

void lucColourMap_SetOpenGLColourFromScaledValue( lucColourMap* cmap, float scaledValue ) {
	lucColour colour;
 
	lucColourMap_GetColourFromScaledValue( cmap, scaledValue, &colour);
	lucColour_SetOpenGLColour( &colour );
}

void lucColourMap_SetOpenGLColourFromValue_ExplicitOpacity( lucColourMap* cmap, double value, float opacity ) {
	lucColour colour;

	lucColourMap_GetColourFromValue( cmap, value, &colour);
	colour.opacity = opacity;
	lucColour_SetOpenGLColour( &colour );
}

void lucColourMap_SetOpenGLColourFromRGB( double red, double green, double blue) {
	glColor4f(
			red,
			green,
			blue,
			1.0 );
}
void lucColourMap_SetOpenGLColourFromRGB_ExplicitOpacity( double red, double green, double blue, float opacity ) {
		glColor4f(
			red,
			green,
			blue,
			opacity );
}
void lucViewportInfo_GetCoordFromPixel( lucViewportInfo* viewportInfo, Pixel_Index xPos, Pixel_Index yPos, Coord coord ) {
	GLdouble modelMatrix[16];
	GLdouble projMatrix[16];
	GLint    viewportArray[4];
	float    depth;

	glViewport( viewportInfo->startx, viewportInfo->starty, viewportInfo->width, viewportInfo->height);
	lucViewportInfo_SetOpenGLCamera( viewportInfo );

	glGetDoublev( GL_MODELVIEW_MATRIX, modelMatrix );
	glGetDoublev( GL_PROJECTION_MATRIX, projMatrix );
	glGetIntegerv( GL_VIEWPORT, viewportArray );

	glReadBuffer(GL_FRONT);
	glReadPixels( xPos, yPos, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth );

	gluUnProject(
		(double) xPos,
		(double) yPos,
		(double) depth,
		modelMatrix,
		projMatrix,
		viewportArray,
		&coord[I_AXIS],
		&coord[J_AXIS],
		&coord[K_AXIS]
	);
}

void lucViewportInfo_SetOpenGLCamera( lucViewportInfo* viewportInfo ) {
	lucViewport* viewport  = viewportInfo->viewport;
	lucCamera*   camera      = viewport->camera;
	Coord        currEyePos;
	
	/*Declarations for the lucStereoAsymetric part*/
	double ratio, radians, wd2, ndfl;
	double left, right, top, bottom;
	#define DTOR 0.0174532925
	GLenum pname = GL_STEREO;
	GLboolean stereo_enabled;

	switch ( camera->stereoType ) {
		case lucMono: case lucStereoToeIn:
			/* set up projection transform */
			glMatrixMode(GL_PROJECTION);  
			glLoadIdentity();

			gluPerspective( camera->aperture, lucViewportInfo_AspectRatio(viewportInfo), 
					viewport->nearClipPlane, viewport->farClipPlane );
			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();

			lucCamera_CurrentEyePosition( camera, currEyePos );

			/* Change viewer location */
			gluLookAt( currEyePos[ I_AXIS ],          currEyePos[ J_AXIS ],          currEyePos[ K_AXIS ], 
					   camera->focalPoint[ I_AXIS ],  camera->focalPoint[ J_AXIS ],  camera->focalPoint[ K_AXIS ], 
					   camera->upDirection[ I_AXIS ], camera->upDirection[ J_AXIS ], camera->upDirection[ K_AXIS ] ); 
			break;
		
		
		case lucStereoAsymmetric:
		
			/* test if the strereo can be enabled */
			
			glGetBooleanv(pname, &stereo_enabled);
			if( stereo_enabled == GL_FALSE)
				Journal_DPrintfL( lucDebug, 2, " Stereo not supported\n");
			else 
				Journal_DPrintfL( lucDebug, 2, " Stereo is supported\n");
				
				
			radians = DTOR * camera->aperture/2.0;
			wd2 = viewport->nearClipPlane * tan(radians);
			
			ndfl = (viewport->nearClipPlane ) / camera->focalLength;
			
			/* View vector for each camera is parrallel and glFrustum() is used to describe the perspective projection*/
			lucCamera_CurrentEyePosition( camera, currEyePos );
			
			if(camera->buffer == lucRight){
				glMatrixMode (GL_PROJECTION);
				glLoadIdentity();
			 	left = -ratio * wd2 - 0.5* camera->eyeSeparation * ndfl;
				right = ratio * wd2 - 0.5* camera->eyeSeparation * ndfl;
				top = wd2;
				bottom = -wd2;
				glFrustum(left, right, bottom, top, viewport->nearClipPlane, viewport->farClipPlane);
		
				glMatrixMode(GL_MODELVIEW);
				/*glDrawBuffer(GL_BACK_RIGHT);*/
				glDrawBuffer(GL_BACK);
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
				glLoadIdentity();
			
				gluLookAt(currEyePos[ I_AXIS ], currEyePos[ J_AXIS ], currEyePos[ K_AXIS ], 
								currEyePos[ I_AXIS ] + camera->focalPoint[ I_AXIS ], 
								currEyePos[ J_AXIS ] + camera->focalPoint[ J_AXIS ],  
								currEyePos[ K_AXIS ] + camera->focalPoint[ K_AXIS ],
								camera->upDirection[ I_AXIS ], camera->upDirection[ J_AXIS ], camera->upDirection[ K_AXIS ] ); 
								
			}
			else if (camera->buffer == lucLeft){
			 	
				glMatrixMode (GL_PROJECTION);
				glLoadIdentity();
				left = -ratio * wd2 + 0.5* camera->eyeSeparation * ndfl;
				right = ratio * wd2 + 0.5* camera->eyeSeparation * ndfl;
				top = wd2;
				bottom = -wd2;
				glFrustum(left, right, bottom, top, viewport->nearClipPlane, viewport->farClipPlane);
		
				glMatrixMode(GL_MODELVIEW);
				glDrawBuffer(GL_BACK);
				/*glDrawBuffer(GL_BACK_RIGHT);*/
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
				glLoadIdentity();
			
				gluLookAt(currEyePos[ I_AXIS ], currEyePos[ J_AXIS ], currEyePos[ K_AXIS ], 
								currEyePos[ I_AXIS ] + camera->focalPoint[ I_AXIS ], 
								currEyePos[ J_AXIS ] + camera->focalPoint[ J_AXIS ],  
								currEyePos[ K_AXIS ] + camera->focalPoint[ K_AXIS ],
								camera->upDirection[ I_AXIS ], camera->upDirection[ J_AXIS ], camera->upDirection[ K_AXIS ] ); 
								
			
			}			
			else {}					
		/*	abort();*/
			break;
	}

   /* Apply scaling factors */
   if (viewport->scaleX != 1.0 || viewport->scaleY != 1.0 || viewport->scaleZ != 1.0) {
      glScalef(viewport->scaleX, viewport->scaleY, viewport->scaleZ);
      /* Enable automatic rescaling of normal vectors when scaling is turned on */
      //glEnable(GL_RESCALE_NORMAL);
      glEnable(GL_NORMALIZE);
   } 
}

void luc_OpenGlSquare( Dimension_Index dim, double* pos, double* normal, double* orientation, double width) {
	double plane2[3], corner[3];
	double imposedOrientation[3];
	double *plane;
	
	if(width == 0.0)
		return;
	
	if(orientation==NULL) {  /* Note the danger - if normal points precisely along z axis this fails !*/
		imposedOrientation[0] = -normal[1];
		imposedOrientation[1] =  normal[0];
		imposedOrientation[2] = 0.0; 
		plane = imposedOrientation;	
	}
	else
		plane = orientation;
	
	if (dim == 2) {
		luc_DrawVector( dim, pos, plane, width, 0.0 );
		return;
	}
	
	if(normal[1] < 0.0) { /* Planes can face (abitrarily) up or down - choose the upward facing normal for consistent shading */
		normal[0] *= -1.0;
		normal[1] *= -1.0;
		normal[2] *= -1.0;
	}
	
	StGermain_VectorCrossProduct(plane2, normal, plane);

	corner[0] = pos[0] - plane[0]*width/2.0 - plane2[0]*width/2.0;
	corner[1] = pos[1] - plane[1]*width/2.0 - plane2[1]*width/2.0;
	corner[2] = pos[2] - plane[2]*width/2.0 - plane2[2]*width/2.0;

	glBegin(GL_QUADS);
		glNormal3dv(normal);
		glVertex3d (corner[0], corner[1], corner[2]);
		glVertex3d (corner[0] + plane2[0]*width, corner[1] + plane2[1]*width, corner[2] + plane2[2]*width);
		glVertex3d (corner[0] + plane[0]*width + plane2[0]*width, corner[1] + plane[1]*width + plane2[1]*width, 
			corner[2] + plane[2]*width + plane2[2]*width);
		glVertex3d (corner[0] + plane[0]*width, corner[1] + plane[1]*width, corner[2] + plane[2]*width);
		
		glVertex3d (corner[0], corner[1], corner[2]);
		glVertex3d (corner[0] + plane2[0]*width, corner[1] + plane2[1]*width, corner[2] + plane2[2]*width);
		glVertex3d (corner[0] + plane[0]*width + plane2[0]*width, corner[1] + plane[1]*width + plane2[1]*width, 
			corner[2] + plane[2]*width + plane2[2]*width);
		glVertex3d (corner[0] + plane[0]*width, corner[1] + plane[1]*width, corner[2] + plane[2]*width);
		
		
	glEnd();
}

void luc_OpenGlFatSquare( Dimension_Index dim, double* pos, double* normal, double* orientation, double width, double thickness) {
	double plane2[3], corner[3], corner2[3];
	double imposedOrientation[3];
	double vertex[8][3];
	double *plane;
	
	if(width == 0.0 || thickness == 0.0)
		return;
	
	if(orientation==NULL) {  /* Note the danger - if normal points precisely along z axis this fails !*/
		imposedOrientation[0] = -normal[1];
		imposedOrientation[1] =  normal[0];
		imposedOrientation[2] = 0.0; 
		plane = imposedOrientation;	
	}
	else
		plane = orientation;
	
	if (dim == 2) {
		luc_DrawVector( dim, pos, plane, width, 0.0 );
		return;
	}
	
	if(normal[1] < 0.0) { /* Planes can face (abitrarily) up or down - choose the upward facing normal for consistent shading */
		normal[0] *= -1.0;
		normal[1] *= -1.0;
		normal[2] *= -1.0;
	}
		
	StGermain_VectorCrossProduct(plane2, normal, plane);

	corner[0] = pos[0] - plane[0]*width/2.0 - plane2[0]*width/2.0 - 0.5 * normal[0] * thickness;
	corner[1] = pos[1] - plane[1]*width/2.0 - plane2[1]*width/2.0 - 0.5 * normal[1] * thickness;
	corner[2] = pos[2] - plane[2]*width/2.0 - plane2[2]*width/2.0 - 0.5 * normal[2] * thickness;

	corner2[0] = corner[0] + normal[0] * thickness;
	corner2[1] = corner[1] + normal[1] * thickness;
	corner2[2] = corner[2] + normal[2] * thickness;
	
	vertex[0][0] = corner[0];
	vertex[0][1] = corner[1];
	vertex[0][2] = corner[2];

	vertex[1][0] = corner[0] + plane2[0]*width;
	vertex[1][1] = corner[1] + plane2[1]*width;
	vertex[1][2] = corner[2] + plane2[2]*width;

	vertex[2][0] = corner[0] + plane[0]*width;
	vertex[2][1] = corner[1] + plane[1]*width;
	vertex[2][2] = corner[2] + plane[2]*width;
	
	vertex[3][0] = corner[0] + plane[0]*width + plane2[0]*width;
	vertex[3][1] = corner[1] + plane[1]*width + plane2[1]*width;
	vertex[3][2] = corner[2] + plane[2]*width + plane2[2]*width;
	
	vertex[4][0] = corner2[0];
	vertex[4][1] = corner2[1];
	vertex[4][2] = corner2[2];

	vertex[5][0] = corner2[0] + plane2[0]*width;
	vertex[5][1] = corner2[1] + plane2[1]*width;
	vertex[5][2] = corner2[2] + plane2[2]*width;

	vertex[6][0] = corner2[0] + plane[0]*width;
	vertex[6][1] = corner2[1] + plane[1]*width;
	vertex[6][2] = corner2[2] + plane[2]*width;
	
	vertex[7][0] = corner2[0] + plane[0]*width + plane2[0]*width;
	vertex[7][1] = corner2[1] + plane[1]*width + plane2[1]*width;
	vertex[7][2] = corner2[2] + plane[2]*width + plane2[2]*width;
	
	
	glBegin(GL_QUADS);
		glNormal3d(-normal[0],-normal[1],-normal[2]);

		glVertex3dv(vertex[0]);
		glVertex3dv(vertex[1]);
		glVertex3dv(vertex[3]);
		glVertex3dv(vertex[2]);
		                    
		glNormal3d(normal[0],normal[1],normal[2]);
		
		glVertex3dv(vertex[4]);
		glVertex3dv(vertex[5]);
		glVertex3dv(vertex[7]);
		glVertex3dv(vertex[6]);
		
		// Keep constant shading for the edges so they don't glint in the sun and distract the viewer 
		
		glVertex3dv(vertex[0]);
		glVertex3dv(vertex[1]);
		glVertex3dv(vertex[5]);
		glVertex3dv(vertex[4]);
				
		glVertex3dv(vertex[2]);
		glVertex3dv(vertex[3]);
		glVertex3dv(vertex[7]);
		glVertex3dv(vertex[6]);
		
		glVertex3dv(vertex[0]);
		glVertex3dv(vertex[2]);
		glVertex3dv(vertex[6]);
		glVertex3dv(vertex[4]);
		
		glVertex3dv(vertex[1]);
		glVertex3dv(vertex[3]);
		glVertex3dv(vertex[7]);
		glVertex3dv(vertex[5]);

		
	glEnd();
	
}

void luc_OpenGlCircle( Dimension_Index dim, double* pos, double* normal, double radius, Index nSides) {
	double vector_polar[3], vector_rect[3];
	double normal_polar[3];
	int theta;

	if (dim == 2) {
		double plane[] = {normal[1], normal[0]};
		luc_DrawVector( dim, pos, plane, radius, 0.0 );
		return;
	}

	/* Convert Normal into spherical coordinates */
	StGermain_RectangularToSpherical( normal_polar, normal, 3);
	

	/* Draw circle */
	glBegin(GL_POLYGON) ; 
		glNormal3dv(normal);
		for (theta = 0 ; theta < 360 ; theta = theta + 360/nSides) {
			/* Find vector from centre of circle to edge of circle */
			vector_polar[0] = radius;
			vector_polar[1] = StGermain_DegreeToRadian( theta );
			vector_polar[2] = atan(			-1.0/
						(tan(normal_polar[2])*cos( normal_polar[1] - vector_polar[1] ) ) );
			/* Correct for arctan domain */
			if ( vector_polar[2] < 0.0) vector_polar[2] = vector_polar[2] + M_PI;
			
			/*Convert back to cartesian coordinates */
			StGermain_SphericalToRectangular( vector_rect, vector_polar, 3 );
			
			/* Find position vectors of edge of circle */
			vector_rect[0] = pos[0] + vector_rect[0];
			vector_rect[1] = pos[1] + vector_rect[1];
			vector_rect[2] = pos[2] + vector_rect[2];
			
			/* Plot them */
			glVertex3dv( vector_rect );
		}
	glEnd();

}

void luc_DrawVector3d( double* pos, double* vector, double scale, double headSize, int segment_count )
{
    /* Length of the drawn vector = vector magnitude * scaling factor */
    float length = scale * StGermain_VectorMagnitude( vector, 3 );
    static int segments;    /* Saves segment count */
    static float *x_coords, *y_coords;  /* Saves arrays of x,y points on circle for set segment count */

    /* Recalc required? Only done first time called and when segment count changes */
    if (segments != segment_count)
    {
        /* Calculate unit circle points when divided into specified segments
         * and store in static variable to re-use every time a vector with the
         * same segment count is drawn */
        segments = segment_count;
        if (x_coords != NULL) Memory_Free(x_coords);
        if (y_coords != NULL) Memory_Free(y_coords);

    	x_coords = Memory_Alloc_Array( float, (segment_count + 1), "Unit Circle X Coords" );
    	y_coords = Memory_Alloc_Array( float, (segment_count + 1), "Unit Circle Y Coords" );

        GLfloat angle; 
        //fprintf(stderr, "Vector Arrow -- Point %f,%f,%f\n", C[0], C[1], C[2]);
        // Loop around in a circle and specify even points along the circle
        // as the vertices for the triangle fan cone, cone base and arrow shaft
        float angle_inc = 2*M_PI / (float)segment_count;
        int idx;
        for(idx = 0; idx <= segments; idx++)
        {
            angle = angle_inc * (float)idx;
            // Calculate x and y position of the next vertex and cylinder normals (unit circle coords)
            x_coords[idx] = sin(angle);
            y_coords[idx] = cos(angle);
        }
    }
    
    /* Render a 3d arrow, cone with base for head, cylinder for shaft */
    glDisable(GL_CULL_FACE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    /* Radius of head */
    double radius = 0.8 * headSize * length;

    glPushMatrix();
    /* Vector is centered on pos[x,y,z]
     * Translate to the point of arrow -> position + vector/2 */
    glTranslated(pos[0] + scale * 0.5 * vector[0], 
                 pos[1] + scale * 0.5 * vector[1], 
                 pos[2] + scale * 0.5 * vector[2]);

    //Rotate to orient the cone
    //...Want to align our z-axis to point along arrow vector: 
    // axis of rotation = (z x vec) cosine of angle = (z . vec)
    double z[3] = {0.0, 0.0, 1.0}, axis[3];
    double rangle;
    StGermain_VectorCrossProduct( axis, z, vector );
    rangle = acos(StGermain_VectorDotProduct( z, vector, 3 ));
    rangle = StGermain_RadianToDegree( rangle );
    glRotated(rangle, axis[0], axis[1], axis[2]);
    //fprintf(stderr, "Axis %f,%f,%f angle %f\n", axis[0], axis[1], axis[2], rangle);

    //Translate back from point by size of arrowhead
    //thus our working coordinate system origin is at the base of head with z-axis aligned with head
    glTranslated(0.0, 0.0, -radius*2);

    /* Render the arrowhead cone and base with two triangle fans */
    /* Don't bother drawing head for tiny vectors */
    double normal[3];
    if ( radius >= 1.0e-10 ) 
    {
        /* Pinnacle vertex is at point of arrow */
        double pinnacle[3] = {0, 0, radius * 2};

        /* First pair of vertices on circle define a triangle when combined with pinnacle */
        double vertex0[3] = {radius * x_coords[0], radius * y_coords[0], 0.0};
        double vertex1[3] = {radius * x_coords[1], radius * y_coords[1], 0.0};
        
        /* Set normal for first triangle */
        StGermain_NormalToPlane( normal, pinnacle, vertex0, vertex1);

        int i;
        for (i=0; i<2; i++)
        {
            glBegin(GL_TRIANGLE_FAN);
            //Pinnacle vertex of cone / centre of base circle
            glNormal3dv(normal);
            glVertex3dv(pinnacle);

            //Subsequent vertices describe outer edges of cone base
            int v;
            for (v=0; v <= segments; v++)
            {
                /* Calc next vertex from unit circle normal */
                vertex1[0] = radius * x_coords[v];
                vertex1[1] = radius * y_coords[v];

                if (i==0 && v > 0)
                {
                    StGermain_NormalToPlane( normal, pinnacle, vertex0, vertex1);
                    glNormal3dv(normal);
                }
                
                /* Draw vertex */
                glVertex2dv(vertex1);

                /* Save previous vertex */
		        memcpy( vertex0, vertex1, sizeof(double) * 3 );
            }
            glEnd();
            
            //Flatten cone for circle base -> set common point to share z-coord
            pinnacle[2] = 0;
            //Normal for back of cone
            normal[0] = 0; normal[1] = 0; normal[2] = -1;
        }
    }

    /* Render a cylinder quad strip for shaft */
    glBegin(GL_QUAD_STRIP);
    int v;
    double shaft_vertex[3][3];
    for (v=0; v <= segments; v++)
    {
        /* Top of shaft */
        shaft_vertex[0][0] = 0.2 * radius * x_coords[v];
        shaft_vertex[0][1] = 0.2 * radius * y_coords[v];
        shaft_vertex[0][2] = 0.0;
        /* Base of shaft */
        shaft_vertex[1][0] = shaft_vertex[0][0];
        shaft_vertex[1][1] = shaft_vertex[0][1];
        shaft_vertex[1][2] = -length + radius*2; /* Shaft length to base = vector length - head size */

        normal[0] = -x_coords[v]; normal[1] = -y_coords[v]; normal[2] = 0;
        glNormal3dv(normal);
        glVertex3dv(shaft_vertex[0]);
        glVertex3dv(shaft_vertex[1]);
        //fprintf(stderr, " idx %d vertex %f,%f,%f --> ", v, shaft_vertex[0][0], shaft_vertex[0][1], shaft_vertex[0][2]);
        //fprintf(stderr, "%f,%f,%f\n", shaft_vertex[1][0], shaft_vertex[1][1], shaft_vertex[1][2]);
    }
    glEnd();
    
    glPopMatrix();
}

void luc_DrawVector( Dimension_Index dim , double* pos, double* vector, double scale, double headSize ) {
	if ( dim == 2 ) {
   	double pos3d[3], vector3d[3];
      pos3d[0] = pos[0]; pos3d[1] = pos[1]; pos3d[2] = 0.0;
      vector3d[0] = vector[0]; vector3d[1] = vector[1]; vector3d[2] = 0.0;
      luc_DrawVector3d(pos3d, vector3d, scale, headSize, 16.0);
	}
	else 
      luc_DrawVector3d(pos, vector, scale, headSize, 16.0);
}

//Leaving this here for now as is a useful debugging function for plotting visable surface normals
#define DEG2RAD (M_PI/180.0)
#define RAD2DEG (180.0/M_PI)
#define crossProduct(a,b,c) \
   (a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
   (a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
   (a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];

#define vecmag(v) \
   sqrt((v)[0] * (v)[0] + (v)[1] * (v)[1] + (v)[2] * (v)[2])
#define dotProduct(v,q) \
   ((v)[0] * (q)[0] + \
   (v)[1] * (q)[1] + \
   (v)[2] * (q)[2])

#define vectorSubtract(a, b, c) \
   (a)[0] = (b)[0] - (c)[0]; \
   (a)[1] = (b)[1] - (c)[1]; \
   (a)[2] = (b)[2] - (c)[2]; \
/* Debugging function, draws a small white/red line representing normal */
void luc_DrawNormalVector( double* pos, double* vector, double scale )
{
	glDisable(GL_LIGHTING);
    /* Length of the drawn vector = vector magnitude * scaling factor */
    double length = scale * sqrt(dotProduct(vector,vector));

    glPushMatrix();
    /* Vector is centered on pos[x,y,z]
     * Translate to the point of arrow -> position + vector/2 */
    glTranslated(pos[0] + scale * 0.5 * vector[0], 
                 pos[1] + scale * 0.5 * vector[1], 
                 pos[2] + scale * 0.5 * vector[2]);

    //Rotate to orient the cone
    //...Want to align our z-axis to point along arrow vector: 
    // axis of rotation = (z x vec) cosine of angle = (z . vec)
    double z[3] = {0.0, 0.0, 1.0}, axis[3];
    double rangle;
    crossProduct( axis, z, vector );
    rangle = acos(dotProduct( z, vector ));
    rangle = RAD2DEG * rangle;
    glRotated(rangle, axis[0], axis[1], axis[2]);
    //fprintf(stderr, "Axis %f,%f,%f angle %f radius %f\n", axis[0], axis[1], axis[2], rangle, radius);

    //Translate back from point by half length
    //thus our working coordinate system origin is halfway down vec with z-axis aligned with head
    glTranslated(0.0, 0.0, -length*0.5);

   /* Render vector as two lines, blue base, red at tip */
   glColor3f(1.0, 0.0, 0.0);
   glBegin(GL_LINES);
      glVertex3d(0, 0, 0);
      glVertex3d(0, 0, length*0.5);
   glEnd();

   glColor3f(0.0, 0.0, 1.0);
   glBegin(GL_LINES);
      glVertex3d(0, 0, 0);
      glVertex3d(0, 0, -length*0.5);
   glEnd();

   glPopMatrix();
	glEnable(GL_LIGHTING);
}

#define OFFSET2D 0.01

void luc_DrawRod( Dimension_Index dim , double* pos, double* vector, double scale ) {
	double magnitude;
	double normal[3], i[] = {1.0, 0.0, 0.0}, k[] = {0.0,0.0,1.0};
	double angle;

	magnitude = StGermain_VectorMagnitude( vector, dim );
	scale *= magnitude;
	angle = StGermain_RadianToDegree( acos(vector[0] / magnitude) );
	glPushMatrix();
	
	if (dim == 2) {
		if (vector[1] < 0) angle = -angle;
		glTranslated(pos[0], pos[1], 0.0);
		glRotated( angle, k[0], k[1], k[2] ); 
		glScaled( scale, scale, scale );
		glDisable(GL_LIGHTING);
		glBegin(GL_LINES);
			glVertex3d(-0.5,0.0, OFFSET2D);
			glVertex3d(0.5, 0.0, OFFSET2D);
		glEnd();
		glEnable(GL_LIGHTING);
	}
	else {
		
		StGermain_VectorCrossProduct( normal, i, vector );
		
		glTranslated( pos[0], pos[1], pos[2] );
		glRotated( angle, normal[0], normal[1], normal[2] ); 
		glScaled( scale, scale, scale );

		glDisable(GL_LIGHTING);
		glBegin(GL_LINES);
			glVertex3f(-0.5, 0.0, 0.0);
			glVertex3f(0.5, 0.0, 0.0);
		glEnd();
		glEnable(GL_LIGHTING);
	}
	glPopMatrix();	
}

void luc_OpenGlDebug_PrintAll( Stream* stream, int debugLevel ) {
    GLint r, g, b, a, depth, sten, aux, dbl;
    glGetIntegerv(GL_RED_BITS, &r);
    glGetIntegerv(GL_GREEN_BITS, &g);
    glGetIntegerv(GL_BLUE_BITS, &b);
    glGetIntegerv(GL_ALPHA_BITS, &a);
    glGetIntegerv(GL_DEPTH_BITS, &depth);
    glGetIntegerv(GL_STENCIL_BITS, &sten);
    glGetIntegerv(GL_AUX_BUFFERS, &aux);
    glGetIntegerv(GL_DOUBLEBUFFER, &dbl);
    Journal_PrintfL(stream, 2, "framebuffer: rgba %d %d %d %d, depth %d, stencil %d, aux bfr %d, double bfr %d\n",
        r, g, b, a, depth, sten, aux, dbl);

    GLint a_test, a_func; float a_ref;
    glGetIntegerv(GL_ALPHA_TEST, &a_test);
    glGetIntegerv(GL_ALPHA_TEST_FUNC, &a_func);
    glGetFloatv(GL_ALPHA_TEST_REF, &a_ref);
    Journal_PrintfL(stream, 2, "alpha testing: %d, func %04x ref %f\n", a_test, a_func, a_ref);

    GLint blend, b_src, b_equ, b_dst;
    GLfloat b_c[4];
    glGetIntegerv(GL_BLEND, &blend);
    glGetIntegerv(GL_BLEND_SRC, &b_src);
    glGetIntegerv(GL_BLEND_EQUATION, &b_equ);
    glGetIntegerv(GL_BLEND_DST, &b_dst);
    glGetFloatv(GL_BLEND_COLOR, b_c);
    Journal_PrintfL(stream, 2, "blending: %d, src %04x equ %04x dst %04x color %.2f %.2f %.2f %.2f\n",
        blend, b_src, b_equ, b_dst, b_c[0], b_c[1], b_c[2], b_c[3]);

    GLint logic, logmode, mask[4], cull, cullmode, front, draw;
    glGetIntegerv(GL_COLOR_LOGIC_OP, &logic);
    glGetIntegerv(GL_LOGIC_OP_MODE, &logmode);
    glGetIntegerv(GL_COLOR_WRITEMASK, mask);
    glGetIntegerv(GL_CULL_FACE, &cull);
    glGetIntegerv(GL_CULL_FACE_MODE, &cullmode);
    glGetIntegerv(GL_FRONT_FACE, &front);
    glGetIntegerv(GL_DRAW_BUFFER, &draw);
    Journal_PrintfL(stream, 2, "draw env: logic %d logmode %04x, mask %d %d %d %d, cull %d cullmode %04x frontface %04x, drawbuffer %04x\n",
        logic, logmode, mask[0], mask[1], mask[2], mask[3], cull, cullmode, front, draw);
        
    Journal_PrintfL(stream, 2, "GL err was %04x\n", glGetError());
    Journal_PrintfL(stream, 2, "***\n");


	luc_OpenGlDebug( stream, GL_LIST_INDEX, debugLevel);
	luc_OpenGlDebug( stream, GL_LIST_MODE, debugLevel);
	luc_OpenGlDebug( stream, GL_MATRIX_MODE, debugLevel);
	luc_OpenGlDebug( stream, GL_MODELVIEW_MATRIX, debugLevel);
	luc_OpenGlDebug( stream, GL_PROJECTION_MATRIX, debugLevel);
	luc_OpenGlDebug( stream, GL_CURRENT_RASTER_POSITION, debugLevel);
	luc_OpenGlDebug( stream, GL_CURRENT_RASTER_TEXTURE_COORDS, debugLevel);
	luc_OpenGlDebug( stream, GL_CURRENT_RASTER_POSITION_VALID, debugLevel);
	luc_OpenGlDebug( stream, GL_TEXTURE_ENV_MODE, debugLevel);
	luc_OpenGlDebug( stream, GL_TEXTURE_ENV_COLOR, debugLevel);
}

void luc_OpenGlDebug( Stream* stream, int mode, int debugLevel){
	GLfloat    modelviewMatrix[16];
	GLfloat    projectionMatrix[16];
	GLint      listIndex;
	GLfloat    listMode;
	GLfloat    matrixMode;
	GLfloat    currentRasterPosition[4];
	GLfloat    currentRasterTextureCoords[4];
	GLfloat    currentRasterPositionValid;
	GLfloat    textureEnvColor[4];
	
	Journal_PrintfL( stream, 2, "In func %s OpenglMode is is %d\n", __func__, mode);
	Journal_PrintfL(stream, 2, "GL err was %04x\n", glGetError());
	
	if(mode == GL_LIST_INDEX){
		glGetIntegerv(GL_LIST_INDEX, &listIndex);
		Journal_PrintfL( stream, debugLevel, "In func '%s' list index is %d\n", __func__, listIndex );	
	} 
	
	if(mode == GL_LIST_MODE){
		glGetFloatv(GL_LIST_MODE, &listMode);
		Journal_PrintfL( stream, debugLevel, "In func '%s' list index is %d\n", __func__, listMode);
	}  
	
	if(mode == GL_MATRIX_MODE){
		glGetFloatv(GL_MATRIX_MODE, &matrixMode);
		Journal_PrintfL( stream, debugLevel, "In func '%s' list index is %d\n", __func__, matrixMode);
	}
		
    if(mode == GL_MODELVIEW_MATRIX){
		glGetFloatv(GL_MODELVIEW_MATRIX, modelviewMatrix); 
		Journal_PrintfL( stream, debugLevel, "In func '%s' Model View matrix is \n",__func__);
		
		Journal_PrintfL( stream, debugLevel, " %f  %f  %f  %f\n", modelviewMatrix[0], modelviewMatrix[1], modelviewMatrix[2],modelviewMatrix[3]);
		Journal_PrintfL( stream, debugLevel, " %f  %f  %f  %f\n", modelviewMatrix[4], modelviewMatrix[5], modelviewMatrix[6],modelviewMatrix[7]);
		
		Journal_PrintfL( stream, debugLevel, " %f  %f  %f  %f\n", modelviewMatrix[8], modelviewMatrix[9], modelviewMatrix[10],modelviewMatrix[11]);
		Journal_PrintfL( stream, debugLevel, " %f  %f  %f  %f\n", modelviewMatrix[12], modelviewMatrix[13], modelviewMatrix[14],modelviewMatrix[15]);
	}
	
	if(mode == GL_PROJECTION_MATRIX){
		glGetFloatv(GL_PROJECTION_MATRIX, projectionMatrix); 
		Journal_PrintfL( stream, debugLevel, "In func '%s' projection matrix is \n",__func__);
		
		Journal_PrintfL( stream, debugLevel, " %f  %f  %f  %f\n", projectionMatrix[0], projectionMatrix[1], projectionMatrix[2],projectionMatrix[3]);
		Journal_PrintfL( stream, debugLevel, " %f  %f  %f  %f\n", projectionMatrix[4], projectionMatrix[5], projectionMatrix[6],projectionMatrix[7]); 
		Journal_PrintfL( stream, debugLevel, " %f  %f  %f  %f\n", projectionMatrix[8], projectionMatrix[9], projectionMatrix[10],projectionMatrix[11]);
		Journal_PrintfL( stream, debugLevel, " %f  %f  %f  %f\n", projectionMatrix[10], projectionMatrix[13], projectionMatrix[14],projectionMatrix[15]);
	}		 
	
	if(mode == GL_CURRENT_RASTER_POSITION){
		glGetFloatv(GL_CURRENT_RASTER_POSITION, currentRasterPosition); 
		Journal_PrintfL( stream, debugLevel, "In func %s current raster position is \n %f %f %f %f  \n",__func__, currentRasterPosition[0], currentRasterPosition[1], currentRasterPosition[2], currentRasterPosition[3]);
	}
	
	if(mode == GL_CURRENT_RASTER_TEXTURE_COORDS){
		glGetFloatv(GL_CURRENT_RASTER_TEXTURE_COORDS, currentRasterTextureCoords);
		Journal_PrintfL( stream, debugLevel, "In func %s current raster texture coords is \n, %f %f %f %f  \n",__func__, currentRasterTextureCoords[0], currentRasterTextureCoords[1],currentRasterTextureCoords[2], currentRasterTextureCoords[3]);
	}
	
	if(mode == GL_CURRENT_RASTER_POSITION_VALID){
		glGetFloatv(GL_CURRENT_RASTER_POSITION_VALID, &currentRasterPositionValid);
		Journal_PrintfL( stream, debugLevel, "In func %s current raster position valid  is \n %f  \n",__func__, currentRasterPositionValid);
	}
	
	if(mode == GL_TEXTURE_ENV_COLOR){
		glGetFloatv(GL_TEXTURE_ENV_COLOR, textureEnvColor);
		Journal_PrintfL( stream, debugLevel, "In func %s Texture Env Color is \n %f %f %f %f  \n",__func__, textureEnvColor[0], textureEnvColor[1], textureEnvColor[2], textureEnvColor[3]);
	} 

}


