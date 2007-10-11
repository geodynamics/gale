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
** $Id: DummyOpenGL.c 740 2007-10-11 08:05:31Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "gl.h"
#include "glu.h"
#include "stdio.h"

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

const Type DummyOpenGL_Type = "DummyOpenGL";

void glAccum (GLenum op, GLfloat value){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) op, (double) value );
}

void glAlphaFunc (GLenum func, GLclampf ref){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) func, (double) ref );
}

GLboolean glAreTexturesResident (GLsizei n, const GLuint *textures, GLboolean *residences){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, ptr, ptr );\n", __func__, (double) n );
	return 0;
}

void glArrayElement (GLint i){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d );\n", __func__, (int) i );
}

void glBegin (GLenum mode){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d );\n", __func__, (int) mode );
	Stream_Indent( stream );
}

void glBindTexture (GLenum target, GLuint texture){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) target, (double) texture );
}

void glBitmap (GLsizei width, GLsizei height, GLfloat xorig, GLfloat yorig, GLfloat xmove, GLfloat ymove, const GLubyte *bitmap){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, ptr );\n", __func__, (double) width, (double) height, (double) xorig, (double) yorig, (double) xmove, (double) ymove );
}

void glBlendFunc (GLenum sfactor, GLenum dfactor){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d, %d );\n", __func__, (int) sfactor, (int) dfactor );
}

void glCallList (GLuint list){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) list );
}

void glCallLists (GLsizei n, GLenum type, const GLvoid *lists){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d, %d, ptr );\n", __func__, (int) n, (int) type );
}

void glClear (GLbitfield mask){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %u );\n", __func__, (unsigned) mask );
}

void glClearAccum (GLfloat red, GLfloat green, GLfloat blue, GLfloat alpha){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) red, (double) green, (double) blue, (double) alpha );
}

void glClearColor (GLclampf red, GLclampf green, GLclampf blue, GLclampf alpha){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) red, (double) green, (double) blue, (double) alpha );
}

void glClearDepth (GLclampd depth){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) depth );
}

void glClearIndex (GLfloat c){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) c );
}

void glClearStencil (GLint s){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d );\n", __func__, (int) s );
}

void glClipPlane (GLenum plane, const GLdouble *equation){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d, ptr );\n", __func__, (int) plane );
}

void glColor3b (GLbyte red, GLbyte green, GLbyte blue){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) red, (double) green, (double) blue );
}

void glColor3bv (const GLbyte *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glColor3d (GLdouble red, GLdouble green, GLdouble blue){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) red, (double) green, (double) blue );
}

void glColor3dv (const GLdouble *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glColor3f (GLfloat red, GLfloat green, GLfloat blue){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) red, (double) green, (double) blue );
}

void glColor3fv (const GLfloat *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glColor3i (GLint red, GLint green, GLint blue){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) red, (double) green, (double) blue );
}

void glColor3iv (const GLint *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glColor3s (GLshort red, GLshort green, GLshort blue){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) red, (double) green, (double) blue );
}

void glColor3sv (const GLshort *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glColor3ub (GLubyte red, GLubyte green, GLubyte blue){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) red, (double) green, (double) blue );
}

void glColor3ubv (const GLubyte *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glColor3ui (GLuint red, GLuint green, GLuint blue){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) red, (double) green, (double) blue );
}

void glColor3uiv (const GLuint *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glColor3us (GLushort red, GLushort green, GLushort blue){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) red, (double) green, (double) blue );
}

void glColor3usv (const GLushort *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glColor4b (GLbyte red, GLbyte green, GLbyte blue, GLbyte alpha){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) red, (double) green, (double) blue, (double) alpha );
}

void glColor4bv (const GLbyte *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glColor4d (GLdouble red, GLdouble green, GLdouble blue, GLdouble alpha){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) red, (double) green, (double) blue, (double) alpha );
}

void glColor4dv (const GLdouble *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glColor4f (GLfloat red, GLfloat green, GLfloat blue, GLfloat alpha){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) red, (double) green, (double) blue, (double) alpha );
}

void glColor4fv (const GLfloat *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glColor4i (GLint red, GLint green, GLint blue, GLint alpha){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) red, (double) green, (double) blue, (double) alpha );
}

void glColor4iv (const GLint *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glColor4s (GLshort red, GLshort green, GLshort blue, GLshort alpha){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) red, (double) green, (double) blue, (double) alpha );
}

void glColor4sv (const GLshort *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glColor4ub (GLubyte red, GLubyte green, GLubyte blue, GLubyte alpha){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) red, (double) green, (double) blue, (double) alpha );
}

void glColor4ubv (const GLubyte *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glColor4ui (GLuint red, GLuint green, GLuint blue, GLuint alpha){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) red, (double) green, (double) blue, (double) alpha );
}

void glColor4uiv (const GLuint *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glColor4us (GLushort red, GLushort green, GLushort blue, GLushort alpha){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) red, (double) green, (double) blue, (double) alpha );
}

void glColor4usv (const GLushort *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glColorMask (GLboolean red, GLboolean green, GLboolean blue, GLboolean alpha){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) red, (double) green, (double) blue, (double) alpha );
}

void glColorMaterial (GLenum face, GLenum mode){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d, %d );\n", __func__, (int) face, (int) mode );
}

void glColorPointer (GLint size, GLenum type, GLsizei stride, const GLvoid *pointer){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, ptr );\n", __func__, (double) size, (double) type, (double) stride );
}

void glCopyPixels (GLint x, GLint y, GLsizei width, GLsizei height, GLenum type){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) x, (double) y, (double) width, (double) height, (double) type );
}

void glCopyTexImage1D (GLenum target, GLint level, GLenum internalFormat, GLint x, GLint y, GLsizei width, GLint border){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) target, (double) level, (double) internalFormat, (double) x, (double) y, (double) width, (double) border );
}

void glCopyTexImage2D (GLenum target, GLint level, GLenum internalFormat, GLint x, GLint y, GLsizei width, GLsizei height, GLint border){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) target, (double) level, (double) internalFormat, (double) x, (double) y, (double) width, (double) height, (double) border );
}

void glCopyTexSubImage1D (GLenum target, GLint level, GLint xoffset, GLint x, GLint y, GLsizei width){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) target, (double) level, (double) xoffset, (double) x, (double) y, (double) width );
}

void glCopyTexSubImage2D (GLenum target, GLint level, GLint xoffset, GLint yoffset, GLint x, GLint y, GLsizei width, GLsizei height){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) target, (double) level, (double) xoffset, (double) yoffset, (double) x, (double) y, (double) width, (double) height );
}

void glCullFace (GLenum mode){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) mode );
}

void glDeleteLists (GLuint list, GLsizei range){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) list, (double) range );
}

void glDeleteTextures (GLsizei n, const GLuint *textures){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, ptr );\n", __func__, (double) n );
}

void glDepthFunc (GLenum func){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) func );
}

void glDepthMask (GLboolean flag){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) flag );
}

void glDepthRange (GLclampd zNear, GLclampd zFar){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) zNear, (double) zFar );
}

void glDisable (GLenum cap){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d );\n", __func__, (int) cap );
}

void glDisableClientState (GLenum array){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) array );
}

void glDrawArrays (GLenum mode, GLint first, GLsizei count){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) mode, (double) first, (double) count );
}

void glDrawBuffer (GLenum mode){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d );\n", __func__, (int) mode );
}

void glDrawElements (GLenum mode, GLsizei count, GLenum type, const GLvoid *indices){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, ptr );\n", __func__, (double) mode, (double) count, (double) type );
}

void glDrawPixels (GLsizei width, GLsizei height, GLenum format, GLenum type, const GLvoid *pixels){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d, %d, %d, %d, ptr );\n", __func__, (int) width, (int) height, (int) format, (int) type );
}

void glEdgeFlag (GLboolean flag){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) flag );
}

void glEdgeFlagPointer (GLsizei stride, const GLvoid *pointer){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, ptr );\n", __func__, (double) stride );
}

void glEdgeFlagv (const GLboolean *flag){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glEnable (GLenum cap){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d );\n", __func__, (int) cap );
}

void glEnableClientState (GLenum array){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d );\n", __func__, (int) array );
}

void glEnd (){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Stream_UnIndent( stream );
	Journal_Printf( stream, "%s( );\n", __func__ );
}

void glEndList (){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Stream_UnIndent( stream );
	Journal_Printf( stream, "%s( );\n", __func__ );
}

void glEvalCoord1d (GLdouble u){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) u );
}

void glEvalCoord1dv (const GLdouble *u){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glEvalCoord1f (GLfloat u){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) u );
}

void glEvalCoord1fv (const GLfloat *u){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glEvalCoord2d (GLdouble u, GLdouble v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) u, (double) v );
}

void glEvalCoord2dv (const GLdouble *u){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glEvalCoord2f (GLfloat u, GLfloat v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) u, (double) v );
}

void glEvalCoord2fv (const GLfloat *u){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glEvalMesh1 (GLenum mode, GLint i1, GLint i2){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) mode, (double) i1, (double) i2 );
}

void glEvalMesh2 (GLenum mode, GLint i1, GLint i2, GLint j1, GLint j2){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) mode, (double) i1, (double) i2, (double) j1, (double) j2 );
}

void glEvalPoint1 (GLint i){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) i );
}

void glEvalPoint2 (GLint i, GLint j){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) i, (double) j );
}

void glFeedbackBuffer (GLsizei size, GLenum type, GLfloat *buffer){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) size, (double) type );
}

void glFinish (){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( );\n", __func__ );
}

void glFlush (){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( );\n", __func__ );
}

void glFogf (GLenum pname, GLfloat param){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) pname, (double) param );
}

void glFogfv (GLenum pname, const GLfloat *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, ptr );\n", __func__, (double) pname );
}

void glFogi (GLenum pname, GLint param){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) pname, (double) param );
}

void glFogiv (GLenum pname, const GLint *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, ptr );\n", __func__, (double) pname );
}

void glFrontFace (GLenum mode){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d );\n", __func__, (int) mode );
}

void glFrustum (GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble zNear, GLdouble zFar){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) left, (double) right, (double) bottom, (double) top, (double) zNear, (double) zFar );
}

GLuint glGenLists (GLsizei range){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) range );
	return 0;
}

void glGenTextures (GLsizei n, GLuint *textures){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, ptr );\n", __func__, (double) n );
}

void glGetBooleanv (GLenum pname, GLboolean *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, ptr );\n", __func__, (double) pname );
}

void glGetClipPlane (GLenum plane, GLdouble *equation){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, ptr );\n", __func__, (double) plane );
}

void glGetDoublev (GLenum pname, GLdouble *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, ptr );\n", __func__, (double) pname );
}

GLenum glGetError (){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( );\n", __func__ );
	return 0;
}

void glGetFloatv (GLenum pname, GLfloat *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, ptr );\n", __func__, (double) pname );
}

void glGetIntegerv (GLenum pname, GLint *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, ptr );\n", __func__, (double) pname );
}

void glGetLightfv (GLenum light, GLenum pname, GLfloat *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) light, (double) pname );
}

void glGetLightiv (GLenum light, GLenum pname, GLint *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) light, (double) pname );
}

void glGetMapdv (GLenum target, GLenum query, GLdouble *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) target, (double) query );
}

void glGetMapfv (GLenum target, GLenum query, GLfloat *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) target, (double) query );
}

void glGetMapiv (GLenum target, GLenum query, GLint *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) target, (double) query );
}

void glGetMaterialfv (GLenum face, GLenum pname, GLfloat *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) face, (double) pname );
}

void glGetMaterialiv (GLenum face, GLenum pname, GLint *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) face, (double) pname );
}

void glGetPixelMapfv (GLenum map, GLfloat *values){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, ptr );\n", __func__, (double) map );
}

void glGetPixelMapuiv (GLenum map, GLuint *values){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, ptr );\n", __func__, (double) map );
}

void glGetPixelMapusv (GLenum map, GLushort *values){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, ptr );\n", __func__, (double) map );
}

void glGetPointerv (GLenum pname, GLvoid* *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, ptr );\n", __func__, (double) pname );
}

void glGetPolygonStipple (GLubyte *mask){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

const GLubyte * glGetString (GLenum name){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) name );
	return 0;
}

void glGetTexEnvfv (GLenum target, GLenum pname, GLfloat *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) target, (double) pname );
}

void glGetTexEnviv (GLenum target, GLenum pname, GLint *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) target, (double) pname );
}

void glGetTexGendv (GLenum coord, GLenum pname, GLdouble *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) coord, (double) pname );
}

void glGetTexGenfv (GLenum coord, GLenum pname, GLfloat *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) coord, (double) pname );
}

void glGetTexGeniv (GLenum coord, GLenum pname, GLint *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) coord, (double) pname );
}

void glGetTexImage (GLenum target, GLint level, GLenum format, GLenum type, GLvoid *pixels){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, ptr );\n", __func__, (double) target, (double) level, (double) format, (double) type );
}

void glGetTexLevelParameterfv (GLenum target, GLint level, GLenum pname, GLfloat *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, ptr );\n", __func__, (double) target, (double) level, (double) pname );
}

void glGetTexLevelParameteriv (GLenum target, GLint level, GLenum pname, GLint *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, ptr );\n", __func__, (double) target, (double) level, (double) pname );
}

void glGetTexParameterfv (GLenum target, GLenum pname, GLfloat *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) target, (double) pname );
}

void glGetTexParameteriv (GLenum target, GLenum pname, GLint *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) target, (double) pname );
}

void glHint (GLenum target, GLenum mode){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) target, (double) mode );
}

void glIndexMask (GLuint mask){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) mask );
}

void glIndexPointer (GLenum type, GLsizei stride, const GLvoid *pointer){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) type, (double) stride );
}

void glIndexd (GLdouble c){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) c );
}

void glIndexdv (const GLdouble *c){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glIndexf (GLfloat c){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) c );
}

void glIndexfv (const GLfloat *c){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glIndexi (GLint c){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) c );
}

void glIndexiv (const GLint *c){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glIndexs (GLshort c){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) c );
}

void glIndexsv (const GLshort *c){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glIndexub (GLubyte c){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) c );
}

void glIndexubv (const GLubyte *c){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glInitNames (){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( );\n", __func__ );
}

void glInterleavedArrays (GLenum format, GLsizei stride, const GLvoid *pointer){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) format, (double) stride );
}

GLboolean glIsEnabled (GLenum cap){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) cap );
	return 0;
}

GLboolean glIsList (GLuint list){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) list );
	return 0;
}

GLboolean glIsTexture (GLuint texture){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) texture );
	return 0;
}

void glLightModelf (GLenum pname, GLfloat param){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) pname, (double) param );
}

void glLightModelfv (GLenum pname, const GLfloat *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, ptr );\n", __func__, (double) pname );
}

void glLightModeli (GLenum pname, GLint param){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d, %d );\n", __func__, (int) pname, (int) param );
}

void glLightModeliv (GLenum pname, const GLint *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, ptr );\n", __func__, (double) pname );
}

void glLightf (GLenum light, GLenum pname, GLfloat param){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) light, (double) pname, (double) param );
}

void glLightfv (GLenum light, GLenum pname, const GLfloat *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) light, (double) pname );
}

void glLighti (GLenum light, GLenum pname, GLint param){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) light, (double) pname, (double) param );
}

void glLightiv (GLenum light, GLenum pname, const GLint *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) light, (double) pname );
}

void glLineStipple (GLint factor, GLushort pattern){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) factor, (double) pattern );
}

void glLineWidth (GLfloat width){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) width );
}

void glListBase (GLuint base){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) base );
}

void glLoadIdentity (){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( );\n", __func__ );
}

void glLoadMatrixd (const GLdouble *m){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glLoadMatrixf (const GLfloat *m){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glLoadName (GLuint name){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) name );
}

void glLogicOp (GLenum opcode){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) opcode );
}

void glMap1d (GLenum target, GLdouble u1, GLdouble u2, GLint stride, GLint order, const GLdouble *points){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g, ptr );\n", __func__, (double) target, (double) u1, (double) u2, (double) stride, (double) order );
}

void glMap1f (GLenum target, GLfloat u1, GLfloat u2, GLint stride, GLint order, const GLfloat *points){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g, ptr );\n", __func__, (double) target, (double) u1, (double) u2, (double) stride, (double) order );
}

void glMap2d (GLenum target, GLdouble u1, GLdouble u2, GLint ustride, GLint uorder, GLdouble v1, GLdouble v2, GLint vstride, GLint vorder, const GLdouble *points){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, ptr );\n", __func__, (double) target, (double) u1, (double) u2, (double) ustride, (double) uorder, (double) v1, (double) v2, (double) vstride, (double) vorder );
}

void glMap2f (GLenum target, GLfloat u1, GLfloat u2, GLint ustride, GLint uorder, GLfloat v1, GLfloat v2, GLint vstride, GLint vorder, const GLfloat *points){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, ptr );\n", __func__, (double) target, (double) u1, (double) u2, (double) ustride, (double) uorder, (double) v1, (double) v2, (double) vstride, (double) vorder );
}

void glMapGrid1d (GLint un, GLdouble u1, GLdouble u2){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) un, (double) u1, (double) u2 );
}

void glMapGrid1f (GLint un, GLfloat u1, GLfloat u2){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) un, (double) u1, (double) u2 );
}

void glMapGrid2d (GLint un, GLdouble u1, GLdouble u2, GLint vn, GLdouble v1, GLdouble v2){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) un, (double) u1, (double) u2, (double) vn, (double) v1, (double) v2 );
}

void glMapGrid2f (GLint un, GLfloat u1, GLfloat u2, GLint vn, GLfloat v1, GLfloat v2){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) un, (double) u1, (double) u2, (double) vn, (double) v1, (double) v2 );
}

void glMaterialf (GLenum face, GLenum pname, GLfloat param){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) face, (double) pname, (double) param );
}

void glMaterialfv (GLenum face, GLenum pname, const GLfloat *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) face, (double) pname );
}

void glMateriali (GLenum face, GLenum pname, GLint param){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) face, (double) pname, (double) param );
}

void glMaterialiv (GLenum face, GLenum pname, const GLint *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) face, (double) pname );
}

void glMatrixMode (GLenum mode){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d );\n", __func__, (int) mode );
}

void glMultMatrixd (const GLdouble *m){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glMultMatrixf (const GLfloat *m){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glNewList (GLuint list, GLenum mode){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %u, %d );\n", __func__, (unsigned) list, (int) mode );
	Stream_Indent( stream );
}

void glNormal3b (GLbyte nx, GLbyte ny, GLbyte nz){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) nx, (double) ny, (double) nz );
}

void glNormal3bv (const GLbyte *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( (pointer argument - %.6g, %.6g, %.6g) );\n", __func__, (double) v[0], (double) v[1], (double) v[2] );
}

void glNormal3d (GLdouble nx, GLdouble ny, GLdouble nz){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) nx, (double) ny, (double) nz );
}

void glNormal3dv (const GLdouble *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( (pointer argument - %.6g, %.6g, %.6g) );\n", __func__, (double) v[0], (double) v[1], (double) v[2] );
}

void glNormal3f (GLfloat nx, GLfloat ny, GLfloat nz){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) nx, (double) ny, (double) nz );
}

void glNormal3fv (const GLfloat *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( (pointer argument - %.6g, %.6g, %.6g) );\n", __func__, (double) v[0], (double) v[1], (double) v[2] );
}

void glNormal3i (GLint nx, GLint ny, GLint nz){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) nx, (double) ny, (double) nz );
}

void glNormal3iv (const GLint *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( (pointer argument - %.6g, %.6g, %.6g) );\n", __func__, (double) v[0], (double) v[1], (double) v[2] );
}

void glNormal3s (GLshort nx, GLshort ny, GLshort nz){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) nx, (double) ny, (double) nz );
}

void glNormal3sv (const GLshort *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( (pointer argument - %.6g, %.6g, %.6g) );\n", __func__, (double) v[0], (double) v[1], (double) v[2] );
}

void glNormalPointer (GLenum type, GLsizei stride, const GLvoid *pointer){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) type, (double) stride );
}

void glOrtho (GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble zNear, GLdouble zFar){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) left, (double) right, (double) bottom, (double) top, (double) zNear, (double) zFar );
}

void glPassThrough (GLfloat token){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) token );
}

void glPixelMapfv (GLenum map, GLsizei mapsize, const GLfloat *values){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) map, (double) mapsize );
}

void glPixelMapuiv (GLenum map, GLsizei mapsize, const GLuint *values){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) map, (double) mapsize );
}

void glPixelMapusv (GLenum map, GLsizei mapsize, const GLushort *values){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) map, (double) mapsize );
}

void glPixelStoref (GLenum pname, GLfloat param){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) pname, (double) param );
}

void glPixelStorei (GLenum pname, GLint param){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d, %d );\n", __func__, (int) pname, (int) param );
}

void glPixelTransferf (GLenum pname, GLfloat param){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) pname, (double) param );
}

void glPixelTransferi (GLenum pname, GLint param){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) pname, (double) param );
}

void glPixelZoom (GLfloat xfactor, GLfloat yfactor){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) xfactor, (double) yfactor );
}

void glPointSize (GLfloat size){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) size );
}

void glPolygonMode (GLenum face, GLenum mode){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d, %d );\n", __func__, (int) face, (int) mode );
}

void glPolygonOffset (GLfloat factor, GLfloat units){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) factor, (double) units );
}

void glPolygonStipple (const GLubyte *mask){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glPopAttrib (){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( );\n", __func__ );
}

void glPopClientAttrib (){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( );\n", __func__ );
}

void glPopMatrix (){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( );\n", __func__ );
}

void glPopName (){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( );\n", __func__ );
}

void glPrioritizeTextures (GLsizei n, const GLuint *textures, const GLclampf *priorities){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, ptr, ptr );\n", __func__, (double) n );
}

void glPushAttrib (GLbitfield mask){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %u );\n", __func__, (unsigned int) mask );
}

void glPushClientAttrib (GLbitfield mask){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %u );\n", __func__, (unsigned int) mask );
}

void glPushMatrix (){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( );\n", __func__ );
}

void glPushName (GLuint name){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) name );
}

void glRasterPos2d (GLdouble x, GLdouble y){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) x, (double) y );
}

void glRasterPos2dv (const GLdouble *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glRasterPos2f (GLfloat x, GLfloat y){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) x, (double) y );
}

void glRasterPos2fv (const GLfloat *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glRasterPos2i (GLint x, GLint y){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d, %d );\n", __func__, (int) x, (int) y );
}

void glRasterPos2iv (const GLint *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glRasterPos2s (GLshort x, GLshort y){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) x, (double) y );
}

void glRasterPos2sv (const GLshort *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glRasterPos3d (GLdouble x, GLdouble y, GLdouble z){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) x, (double) y, (double) z );
}

void glRasterPos3dv (const GLdouble *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glRasterPos3f (GLfloat x, GLfloat y, GLfloat z){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) x, (double) y, (double) z );
}

void glRasterPos3fv (const GLfloat *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glRasterPos3i (GLint x, GLint y, GLint z){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) x, (double) y, (double) z );
}

void glRasterPos3iv (const GLint *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glRasterPos3s (GLshort x, GLshort y, GLshort z){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) x, (double) y, (double) z );
}

void glRasterPos3sv (const GLshort *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glRasterPos4d (GLdouble x, GLdouble y, GLdouble z, GLdouble w){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) x, (double) y, (double) z, (double) w );
}

void glRasterPos4dv (const GLdouble *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glRasterPos4f (GLfloat x, GLfloat y, GLfloat z, GLfloat w){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) x, (double) y, (double) z, (double) w );
}

void glRasterPos4fv (const GLfloat *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glRasterPos4i (GLint x, GLint y, GLint z, GLint w){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) x, (double) y, (double) z, (double) w );
}

void glRasterPos4iv (const GLint *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glRasterPos4s (GLshort x, GLshort y, GLshort z, GLshort w){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) x, (double) y, (double) z, (double) w );
}

void glRasterPos4sv (const GLshort *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glReadBuffer (GLenum mode){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) mode );
}

void glReadPixels (GLint x, GLint y, GLsizei width, GLsizei height, GLenum format, GLenum type, GLvoid *pixels){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d, %d, %d, %d, %d, %d, ptr );\n", __func__, (int) x, (int) y, (int) width, (int) height, (int) format, (int) type );
}

void glRectd (GLdouble x1, GLdouble y1, GLdouble x2, GLdouble y2){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) x1, (double) y1, (double) x2, (double) y2 );
}

void glRectdv (const GLdouble *v1, const GLdouble *v2){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, ptr );\n", __func__ );
}

void glRectf (GLfloat x1, GLfloat y1, GLfloat x2, GLfloat y2){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) x1, (double) y1, (double) x2, (double) y2 );
}

void glRectfv (const GLfloat *v1, const GLfloat *v2){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, ptr );\n", __func__ );
}

void glRecti (GLint x1, GLint y1, GLint x2, GLint y2){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d, %d, %d, %d );\n", __func__, (int) x1, (int) y1, (int) x2, (int) y2 );
}

void glRectiv (const GLint *v1, const GLint *v2){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, ptr );\n", __func__ );
}

void glRects (GLshort x1, GLshort y1, GLshort x2, GLshort y2){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) x1, (double) y1, (double) x2, (double) y2 );
}

void glRectsv (const GLshort *v1, const GLshort *v2){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, ptr );\n", __func__ );
}

GLint glRenderMode (GLenum mode){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) mode );
	return 0;
}

void glRotated (GLdouble angle, GLdouble x, GLdouble y, GLdouble z){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) angle, (double) x, (double) y, (double) z );
}

void glRotatef (GLfloat angle, GLfloat x, GLfloat y, GLfloat z){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) angle, (double) x, (double) y, (double) z );
}

void glScaled (GLdouble x, GLdouble y, GLdouble z){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) x, (double) y, (double) z );
}

void glScalef (GLfloat x, GLfloat y, GLfloat z){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) x, (double) y, (double) z );
}

void glScissor (GLint x, GLint y, GLsizei width, GLsizei height){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d, %d, %d, %d );\n", __func__, (int) x, (int) y, (int) width, (int) height );
}

void glSelectBuffer (GLsizei size, GLuint *buffer){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d, ptr );\n", __func__, (int) size );
}

void glShadeModel (GLenum mode){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) mode );
}

void glStencilFunc (GLenum func, GLint ref, GLuint mask){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) func, (double) ref, (double) mask );
}

void glStencilMask (GLuint mask){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) mask );
}

void glStencilOp (GLenum fail, GLenum zfail, GLenum zpass){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) fail, (double) zfail, (double) zpass );
}

void glTexCoord1d (GLdouble s){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) s );
}

void glTexCoord1dv (const GLdouble *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glTexCoord1f (GLfloat s){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) s );
}

void glTexCoord1fv (const GLfloat *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glTexCoord1i (GLint s){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) s );
}

void glTexCoord1iv (const GLint *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glTexCoord1s (GLshort s){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) s );
}

void glTexCoord1sv (const GLshort *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glTexCoord2d (GLdouble s, GLdouble t){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) s, (double) t );
}

void glTexCoord2dv (const GLdouble *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glTexCoord2f (GLfloat s, GLfloat t){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) s, (double) t );
}

void glTexCoord2fv (const GLfloat *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glTexCoord2i (GLint s, GLint t){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) s, (double) t );
}

void glTexCoord2iv (const GLint *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glTexCoord2s (GLshort s, GLshort t){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) s, (double) t );
}

void glTexCoord2sv (const GLshort *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glTexCoord3d (GLdouble s, GLdouble t, GLdouble r){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) s, (double) t, (double) r );
}

void glTexCoord3dv (const GLdouble *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glTexCoord3f (GLfloat s, GLfloat t, GLfloat r){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) s, (double) t, (double) r );
}

void glTexCoord3fv (const GLfloat *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glTexCoord3i (GLint s, GLint t, GLint r){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) s, (double) t, (double) r );
}

void glTexCoord3iv (const GLint *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glTexCoord3s (GLshort s, GLshort t, GLshort r){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) s, (double) t, (double) r );
}

void glTexCoord3sv (const GLshort *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glTexCoord4d (GLdouble s, GLdouble t, GLdouble r, GLdouble q){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) s, (double) t, (double) r, (double) q );
}

void glTexCoord4dv (const GLdouble *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glTexCoord4f (GLfloat s, GLfloat t, GLfloat r, GLfloat q){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) s, (double) t, (double) r, (double) q );
}

void glTexCoord4fv (const GLfloat *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glTexCoord4i (GLint s, GLint t, GLint r, GLint q){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) s, (double) t, (double) r, (double) q );
}

void glTexCoord4iv (const GLint *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glTexCoord4s (GLshort s, GLshort t, GLshort r, GLshort q){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) s, (double) t, (double) r, (double) q );
}

void glTexCoord4sv (const GLshort *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void glTexCoordPointer (GLint size, GLenum type, GLsizei stride, const GLvoid *pointer){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, ptr );\n", __func__, (double) size, (double) type, (double) stride );
}

void glTexEnvf (GLenum target, GLenum pname, GLfloat param){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) target, (double) pname, (double) param );
}

void glTexEnvfv (GLenum target, GLenum pname, const GLfloat *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) target, (double) pname );
}

void glTexEnvi (GLenum target, GLenum pname, GLint param){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) target, (double) pname, (double) param );
}

void glTexEnviv (GLenum target, GLenum pname, const GLint *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) target, (double) pname );
}

void glTexGend (GLenum coord, GLenum pname, GLdouble param){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) coord, (double) pname, (double) param );
}

void glTexGendv (GLenum coord, GLenum pname, const GLdouble *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) coord, (double) pname );
}

void glTexGenf (GLenum coord, GLenum pname, GLfloat param){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) coord, (double) pname, (double) param );
}

void glTexGenfv (GLenum coord, GLenum pname, const GLfloat *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) coord, (double) pname );
}

void glTexGeni (GLenum coord, GLenum pname, GLint param){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) coord, (double) pname, (double) param );
}

void glTexGeniv (GLenum coord, GLenum pname, const GLint *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, ptr );\n", __func__, (double) coord, (double) pname );
}

#ifdef __APPLE__
void glTexImage1D (GLenum target, GLint level, GLenum internalformat, GLsizei width, GLint border, GLenum format, GLenum type, const GLvoid *pixels){
#else
void glTexImage1D (GLenum target, GLint level, GLint internalformat, GLsizei width, GLint border, GLenum format, GLenum type, const GLvoid *pixels){
#endif
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d, %d, %d, %d, %d, %d, %d, ptr );\n", __func__, (int) target, (int) level, (int) internalformat, (int) width, (int) border, (int) format, (int) type );
}

#ifdef __APPLE__
void glTexImage2D (GLenum target, GLint level, GLenum internalformat, GLsizei width, GLsizei height, GLint border, GLenum format, GLenum type, const GLvoid *pixels){
#else
void glTexImage2D (GLenum target, GLint level, GLint internalformat, GLsizei width, GLsizei height, GLint border, GLenum format, GLenum type, const GLvoid *pixels){
#endif
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d, %d, %d, %d, %d, %d, %d, %d, ptr );\n", __func__, (int) target, (int) level, (int) internalformat, (int) width, (int) height, (int) border, (int) format, (int) type );
}

void glTexParameterf (GLenum target, GLenum pname, GLfloat param){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d, %d, %.6g );\n", __func__, (int) target, (int) pname, (double) param );
}

void glTexParameterfv (GLenum target, GLenum pname, const GLfloat *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d, %d, ptr );\n", __func__, (int) target, (int) pname );
}

void glTexParameteri (GLenum target, GLenum pname, GLint param){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d, %d, %d );\n", __func__, (int) target, (int) pname, (int) param );
}

void glTexParameteriv (GLenum target, GLenum pname, const GLint *params){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d, %d, ptr );\n", __func__, (int) target, (int) pname );
}

void glTexSubImage1D (GLenum target, GLint level, GLint xoffset, GLsizei width, GLenum format, GLenum type, const GLvoid *pixels){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, ptr );\n", __func__, (double) target, (double) level, (double) xoffset, (double) width, (double) format, (double) type );
}

void glTexSubImage2D (GLenum target, GLint level, GLint xoffset, GLint yoffset, GLsizei width, GLsizei height, GLenum format, GLenum type, const GLvoid *pixels){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, ptr );\n", __func__, (double) target, (double) level, (double) xoffset, (double) yoffset, (double) width, (double) height, (double) format, (double) type );
}

void glTranslated (GLdouble x, GLdouble y, GLdouble z){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) x, (double) y, (double) z );
}

void glTranslatef (GLfloat x, GLfloat y, GLfloat z){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) x, (double) y, (double) z );
}

void glVertex2d (GLdouble x, GLdouble y){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) x, (double) y );
}

void glVertex2dv (const GLdouble *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( (pointer argument - %.6g, %.6g) );\n", __func__, (double) v[0], (double) v[1] );
}

void glVertex2f (GLfloat x, GLfloat y){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) x, (double) y );
}

void glVertex2fv (const GLfloat *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( (pointer argument - %.6g, %.6g) );\n", __func__, (double) v[0], (double) v[1] );
}

void glVertex2i (GLint x, GLint y){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) x, (double) y );
}

void glVertex2iv (const GLint *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( (pointer argument - %.6g, %.6g) );\n", __func__, (double) v[0], (double) v[1] );
}

void glVertex2s (GLshort x, GLshort y){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g );\n", __func__, (double) x, (double) y );
}

void glVertex2sv (const GLshort *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( (pointer argument - %.6g, %.6g) );\n", __func__, (double) v[0], (double) v[1] );
}

void glVertex3d (GLdouble x, GLdouble y, GLdouble z){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) x, (double) y, (double) z );
}

void glVertex3dv (const GLdouble *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( (pointer argument - %.6g, %.6g, %.6g) );\n", __func__, (double) v[0], (double) v[1], (double) v[2] );
}

void glVertex3f (GLfloat x, GLfloat y, GLfloat z){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) x, (double) y, (double) z );
}

void glVertex3fv (const GLfloat *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( (pointer argument - %.6g, %.6g, %.6g) );\n", __func__, (double) v[0], (double) v[1], (double) v[2] );
}

void glVertex3i (GLint x, GLint y, GLint z){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) x, (double) y, (double) z );
}

void glVertex3iv (const GLint *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( (pointer argument - %.6g, %.6g, %.6g) );\n", __func__, (double) v[0], (double) v[1], (double) v[2] );
}

void glVertex3s (GLshort x, GLshort y, GLshort z){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g );\n", __func__, (double) x, (double) y, (double) z );
}

void glVertex3sv (const GLshort *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( (pointer argument - %.6g, %.6g, %.6g) );\n", __func__, (double) v[0], (double) v[1], (double) v[2] );
}

void glVertex4d (GLdouble x, GLdouble y, GLdouble z, GLdouble w){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) x, (double) y, (double) z, (double) w );
}

void glVertex4dv (const GLdouble *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( (pointer argument - %.6g, %.6g, %.6g, %.6g) );\n", __func__, (double) v[0], (double) v[1], (double) v[2], (double) v[3] );
}

void glVertex4f (GLfloat x, GLfloat y, GLfloat z, GLfloat w){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) x, (double) y, (double) z, (double) w );
}

void glVertex4fv (const GLfloat *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( (pointer argument - %.6g, %.6g, %.6g, %.6g) );\n", __func__, (double) v[0], (double) v[1], (double) v[2], (double) v[3] );
}

void glVertex4i (GLint x, GLint y, GLint z, GLint w){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) x, (double) y, (double) z, (double) w );
}

void glVertex4iv (const GLint *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( (pointer argument - %.6g, %.6g, %.6g, %.6g) );\n", __func__, (double) v[0], (double) v[1], (double) v[2], (double) v[3] );
}

void glVertex4s (GLshort x, GLshort y, GLshort z, GLshort w){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) x, (double) y, (double) z, (double) w );
}

void glVertex4sv (const GLshort *v){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( (pointer argument - %.6g, %.6g, %.6g, %.6g) );\n", __func__, (double) v[0], (double) v[1], (double) v[2], (double) v[3] );
}

void glVertexPointer (GLint size, GLenum type, GLsizei stride, const GLvoid *pointer){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, ptr );\n", __func__, (double) size, (double) type, (double) stride );
}

void glViewport (GLint x, GLint y, GLsizei width, GLsizei height){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %d, %d, %d, %d );\n", __func__, (int) x, (int) y, (int) width, (int) height );
}




void gluBeginCurve (GLUnurbs* nurb){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void gluBeginPolygon (GLUtesselator* tess){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void gluBeginSurface (GLUnurbs* nurb){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void gluBeginTrim (GLUnurbs* nurb){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

GLint gluBuild1DMipmapLevels (GLenum target, GLint internalFormat, GLsizei width, GLenum format, GLenum type, GLint level, GLint base, GLint max, const void *data){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, ptr );\n", __func__, (double) target, (double) internalFormat, (double) width, (double) format, (double) type, (double) level, (double) base, (double) max );
	return 0;
}

GLint gluBuild1DMipmaps (GLenum target, GLint internalFormat, GLsizei width, GLenum format, GLenum type, const void *data){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g, ptr );\n", __func__, (double) target, (double) internalFormat, (double) width, (double) format, (double) type );
	return 0;
}

GLint gluBuild2DMipmapLevels (GLenum target, GLint internalFormat, GLsizei width, GLsizei height, GLenum format, GLenum type, GLint level, GLint base, GLint max, const void *data){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, ptr );\n", __func__, (double) target, (double) internalFormat, (double) width, (double) height, (double) format, (double) type, (double) level, (double) base, (double) max );
	return 0;
}

GLint gluBuild2DMipmaps (GLenum target, GLint internalFormat, GLsizei width, GLsizei height, GLenum format, GLenum type, const void *data){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, ptr );\n", __func__, (double) target, (double) internalFormat, (double) width, (double) height, (double) format, (double) type );
	return 0;
}

GLint gluBuild3DMipmapLevels (GLenum target, GLint internalFormat, GLsizei width, GLsizei height, GLsizei depth, GLenum format, GLenum type, GLint level, GLint base, GLint max, const void *data){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, ptr );\n", __func__, (double) target, (double) internalFormat, (double) width, (double) height, (double) depth, (double) format, (double) type, (double) level, (double) base, (double) max );
	return 0;
}

GLint gluBuild3DMipmaps (GLenum target, GLint internalFormat, GLsizei width, GLsizei height, GLsizei depth, GLenum format, GLenum type, const void *data){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, ptr );\n", __func__, (double) target, (double) internalFormat, (double) width, (double) height, (double) depth, (double) format, (double) type );
	return 0;
}

GLboolean gluCheckExtension (const GLubyte *extName, const GLubyte *extString){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, ptr );\n", __func__ );
	return 0;
}

void gluCylinder (GLUquadric* quad, GLdouble base, GLdouble top, GLdouble height, GLint slices, GLint stacks){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, %.6g, %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) base, (double) top, (double) height, (double) slices, (double) stacks );
}

void gluDeleteNurbsRenderer (GLUnurbs* nurb){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void gluDeleteQuadric (GLUquadric* quad){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void gluDeleteTess (GLUtesselator* tess){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void gluDisk (GLUquadric* quad, GLdouble inner, GLdouble outer, GLint slices, GLint loops){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) inner, (double) outer, (double) slices, (double) loops );
}

void gluEndCurve (GLUnurbs* nurb){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void gluEndPolygon (GLUtesselator* tess){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void gluEndSurface (GLUnurbs* nurb){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void gluEndTrim (GLUnurbs* nurb){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

const GLubyte * gluErrorString (GLenum error){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) error );
	return 0;
}

void gluGetNurbsProperty (GLUnurbs* nurb, GLenum property, GLfloat* data){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, %.6g, ptr );\n", __func__, (double) property );
}

const GLubyte * gluGetString (GLenum name){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g );\n", __func__, (double) name );
	return 0;
}

void gluGetTessProperty (GLUtesselator* tess, GLenum which, GLdouble* data){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, %.6g, ptr );\n", __func__, (double) which );
}

void gluLoadSamplingMatrices (GLUnurbs* nurb, const GLfloat *model, const GLfloat *perspective, const GLint *view){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, ptr, ptr, ptr );\n", __func__ );
}

void gluLookAt (GLdouble eyeX, GLdouble eyeY, GLdouble eyeZ, GLdouble centerX, GLdouble centerY, GLdouble centerZ, GLdouble upX, GLdouble upY, GLdouble upZ){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) eyeX, (double) eyeY, (double) eyeZ, (double) centerX, (double) centerY, (double) centerZ, (double) upX, (double) upY, (double) upZ );
}

GLUnurbs* gluNewNurbsRenderer (){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( );\n", __func__ );
	return 0;
}

GLUquadric* gluNewQuadric (){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( );\n", __func__ );
	return 0;
}

GLUtesselator* gluNewTess (){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( );\n", __func__ );
	return 0;
}

void gluNextContour (GLUtesselator* tess, GLenum type){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, %.6g );\n", __func__, (double) type );
}

void gluNurbsCallbackData (GLUnurbs* nurb, GLvoid* userData){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, ptr );\n", __func__ );
}

void gluNurbsCallbackDataEXT (GLUnurbs* nurb, GLvoid* userData){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, ptr );\n", __func__ );
}

void gluNurbsCurve (GLUnurbs* nurb, GLint knotCount, GLfloat *knots, GLint stride, GLfloat *control, GLint order, GLenum type){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, %.6g, ptr, %.6g, ptr, %.6g, %.6g );\n", __func__, (double) knotCount, (double) stride, (double) order, (double) type );
}

void gluNurbsProperty (GLUnurbs* nurb, GLenum property, GLfloat value){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, %.6g, %.6g );\n", __func__, (double) property, (double) value );
}

void gluNurbsSurface (GLUnurbs* nurb, GLint sKnotCount, GLfloat* sKnots, GLint tKnotCount, GLfloat* tKnots, GLint sStride, GLint tStride, GLfloat* control, GLint sOrder, GLint tOrder, GLenum type){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, %.6g, ptr, %.6g, ptr, %.6g, %.6g, ptr, %.6g, %.6g, %.6g );\n", __func__, (double) sKnotCount, (double) tKnotCount, (double) sStride, (double) tStride, (double) sOrder, (double) tOrder, (double) type );
}

void gluOrtho2D (GLdouble left, GLdouble right, GLdouble bottom, GLdouble top){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) left, (double) right, (double) bottom, (double) top );
}

void gluPartialDisk (GLUquadric* quad, GLdouble inner, GLdouble outer, GLint slices, GLint loops, GLdouble start, GLdouble sweep){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) inner, (double) outer, (double) slices, (double) loops, (double) start, (double) sweep );
}

void gluPerspective (GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g );\n", __func__, (double) fovy, (double) aspect, (double) zNear, (double) zFar );
}

void gluPickMatrix (GLdouble x, GLdouble y, GLdouble delX, GLdouble delY, GLint *viewport){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, ptr );\n", __func__, (double) x, (double) y, (double) delX, (double) delY );
}

GLint gluProject (GLdouble objX, GLdouble objY, GLdouble objZ, const GLdouble *model, const GLdouble *proj, const GLint *view, GLdouble* winX, GLdouble* winY, GLdouble* winZ){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, ptr, ptr, ptr, ptr, ptr, ptr );\n", __func__, (double) objX, (double) objY, (double) objZ );
	return 0;
}

void gluPwlCurve (GLUnurbs* nurb, GLint count, GLfloat* data, GLint stride, GLenum type){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, %.6g, ptr, %.6g, %.6g );\n", __func__, (double) count, (double) stride, (double) type );
}

void gluQuadricDrawStyle (GLUquadric* quad, GLenum draw){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, %.6g );\n", __func__, (double) draw );
}

void gluQuadricNormals (GLUquadric* quad, GLenum normal){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, %.6g );\n", __func__, (double) normal );
}

void gluQuadricOrientation (GLUquadric* quad, GLenum orientation){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, %.6g );\n", __func__, (double) orientation );
}

void gluQuadricTexture (GLUquadric* quad, GLboolean texture){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, %.6g );\n", __func__, (double) texture );
}

GLint gluScaleImage (GLenum format, GLsizei wIn, GLsizei hIn, GLenum typeIn, const void *dataIn, GLsizei wOut, GLsizei hOut, GLenum typeOut, GLvoid* dataOut){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, ptr, %.6g, %.6g, %.6g, ptr );\n", __func__, (double) format, (double) wIn, (double) hIn, (double) typeIn, (double) wOut, (double) hOut, (double) typeOut );
	return 0;
}

void gluSphere (GLUquadric* quad, GLdouble radius, GLint slices, GLint stacks){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, %.6g, %.6g, %.6g );\n", __func__, (double) radius, (double) slices, (double) stacks );
}

void gluTessBeginContour (GLUtesselator* tess){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void gluTessBeginPolygon (GLUtesselator* tess, GLvoid* data){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, ptr );\n", __func__ );
}

void gluTessEndContour (GLUtesselator* tess){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void gluTessEndPolygon (GLUtesselator* tess){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr );\n", __func__ );
}

void gluTessNormal (GLUtesselator* tess, GLdouble valueX, GLdouble valueY, GLdouble valueZ){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, %.6g, %.6g, %.6g );\n", __func__, (double) valueX, (double) valueY, (double) valueZ );
}

void gluTessProperty (GLUtesselator* tess, GLenum which, GLdouble data){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, %.6g, %.6g );\n", __func__, (double) which, (double) data );
}

void gluTessVertex (GLUtesselator* tess, GLdouble *location, GLvoid* data){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr, ptr, ptr );\n", __func__ );
}

GLint gluUnProject (GLdouble winX, GLdouble winY, GLdouble winZ, const GLdouble *model, const GLdouble *proj, const GLint *view, GLdouble* objX, GLdouble* objY, GLdouble* objZ){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, ptr, ptr, ptr, ptr, ptr, ptr );\n", __func__, (double) winX, (double) winY, (double) winZ );
	return 0;
}

GLint gluUnProject4 (GLdouble winX, GLdouble winY, GLdouble winZ, GLdouble clipW, const GLdouble *model, const GLdouble *proj, const GLint *view, GLdouble near, GLdouble far, GLdouble* objX, GLdouble* objY, GLdouble* objZ, GLdouble* objW){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( %.6g, %.6g, %.6g, %.6g, ptr, ptr, ptr, %.6g, %.6g, ptr, ptr, ptr, ptr );\n", __func__, (double) winX, (double) winY, (double) winZ, (double) clipW, (double) near, (double) far );
	return 0;
}


#if 0
/* These lines are commented out because it wont compile on the mac otherwise */
void gluTessCallback (GLUtesselator* tess, GLenum which, _GLUfuncptr CallBackFunc){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr %.6g );\n", __func__, (double) which );
}

void gluNurbsCallback (GLUnurbs* nurb, GLenum which, _GLUfuncptr CallBackFunc){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr %.6g );\n", __func__, (double) which );
}

void gluQuadricCallback (GLUquadric* quad, GLenum which, _GLUfuncptr CallBackFunc){
	Stream* stream = Journal_Register( Info_Type, DummyOpenGL_Type ); 
	Journal_Printf( stream, "%s( ptr %.6g );\n", __func__, (double) which );
}
#endif
