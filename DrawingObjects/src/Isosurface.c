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
** $Id: Isosurface.c 791 2008-09-01 02:09:06Z JulianGiordani $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include <glucifer/Base/Base.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>

#include "types.h"
#include "OpenGLDrawingObject.h"
#include "Isosurface.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>
#include <string.h>

#ifndef MASTER
	#define MASTER 0
#endif

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucIsosurface_Type = "lucIsosurface";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucIsosurface* _lucIsosurface_New(  LUCISOSURFACE_DEFARGS  ) 
{
	lucIsosurface*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucIsosurface) );
	self = (lucIsosurface*) _lucOpenGLDrawingObject_New(  LUCOPENGLDRAWINGOBJECT_PASSARGS  );
	
	return self;
}

void _lucIsosurface_Init( 
		lucIsosurface*                                     self,
		FieldVariable*                                     isosurfaceField,
		double                                             isovalue,
		IJK                                                resolution,
		Bool                                               drawWalls,
		Bool                                               wireframe,
		Bool                                               cullFrontFace,
		Bool                                               cullBackFace,
		Name                                               colourName,
		lucColourMap*                                      colourMap,
		FieldVariable*                                     colourField,
		FieldVariable*                                     maskField,
		lucDrawingObjectMask*                              mask )
{
	self->isosurfaceField = isosurfaceField;
	self->isovalue        = isovalue;
	lucColour_FromString( &self->colour, colourName );
	memcpy( self->resolution, resolution, sizeof(IJK) );
	self->drawWalls       = drawWalls;
	self->wireframe       = wireframe;
	self->cullFrontFace   = cullFrontFace;
	self->cullBackFace    = cullBackFace;
	self->colourMap       = colourMap;
	self->colourField     = colourField;
	self->maskField       = maskField;
	memcpy( &self->mask, mask, sizeof(lucDrawingObjectMask) );

	self->trianglesAlloced = 100;
	self->triangleList = Memory_Alloc_Array( Surface_Triangle, self->trianglesAlloced, "triangleList" );
}

void _lucIsosurface_Delete( void* drawingObject ) {
	lucIsosurface*  self = (lucIsosurface*)drawingObject;

	Memory_Free( self->triangleList );

	_lucOpenGLDrawingObject_Delete( self );
}

void _lucIsosurface_Print( void* drawingObject, Stream* stream ) {
	lucIsosurface*  self = (lucIsosurface*)drawingObject;

	_lucOpenGLDrawingObject_Print( self, stream );
}

void* _lucIsosurface_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucIsosurface*  self = (lucIsosurface*)drawingObject;
	lucIsosurface* newDrawingObject;

	newDrawingObject = _lucOpenGLDrawingObject_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newDrawingObject;
}


void* _lucIsosurface_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(lucIsosurface);
	Type                                                             type = lucIsosurface_Type;
	Stg_Class_DeleteFunction*                                     _delete = _lucIsosurface_Delete;
	Stg_Class_PrintFunction*                                       _print = _lucIsosurface_Print;
	Stg_Class_CopyFunction*                                         _copy = NULL;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _lucIsosurface_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _lucIsosurface_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _lucIsosurface_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _lucIsosurface_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _lucIsosurface_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _lucIsosurface_Destroy;
	lucDrawingObject_SetupFunction*                                _setup = _lucIsosurface_Setup;
	lucDrawingObject_DrawFunction*                                  _draw = _lucIsosurface_Draw;
	lucDrawingObject_CleanUpFunction*                            _cleanUp = _lucIsosurface_CleanUp;
	lucOpenGLDrawingObject_BuildDisplayListFunction*    _buildDisplayList = _lucIsosurface_BuildDisplayList;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	return (void*) _lucIsosurface_New(  LUCISOSURFACE_PASSARGS  );
}

void _lucIsosurface_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucIsosurface*         self               = (lucIsosurface*)drawingObject;
	FieldVariable*         isosurfaceField;
	FieldVariable*         colourField;
	FieldVariable*         maskField;
	lucColourMap*          colourMap;
	Index                  defaultResolution;
	IJK                    resolution;
	double                 isovalue;
	lucDrawingObjectMask   mask;

	/* Construct Parent */
	_lucOpenGLDrawingObject_AssignFromXML( self, cf, data );

	isosurfaceField = Stg_ComponentFactory_ConstructByKey( cf, self->name, "IsosurfaceField", FieldVariable, True,  data );
	colourMap       = Stg_ComponentFactory_ConstructByKey( cf, self->name, "ColourMap",       lucColourMap,  False, data );
	colourField     = Stg_ComponentFactory_ConstructByKey( cf, self->name, "ColourField",     FieldVariable, False, data );
	maskField       = Stg_ComponentFactory_ConstructByKey( cf, self->name, "MaskField",       FieldVariable, False, data );

	defaultResolution = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolution", 64 );
	resolution[ I_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionX", defaultResolution );
	resolution[ J_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionY", defaultResolution );
	resolution[ K_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionZ", defaultResolution );

	lucDrawingObjectMask_Construct( &mask, self->name, cf, data );

	isovalue = Stg_ComponentFactory_GetDouble( cf, self->name, "isovalue", 0.0 );
	
	_lucIsosurface_Init( 
			self, 
			isosurfaceField,
			isovalue,
			resolution,
			Stg_ComponentFactory_GetBool( cf, self->name, "drawWalls", False ),
			Stg_ComponentFactory_GetBool( cf, self->name, "wireframe", False ),
			Stg_ComponentFactory_GetBool( cf, self->name, "cullFrontFace", False ),
			Stg_ComponentFactory_GetBool( cf, self->name, "cullBackFace", False ),
			Stg_ComponentFactory_GetString( cf, self->name, "colour", "White" ),
			colourMap,
			colourField,
			maskField,
			&mask );
}

void _lucIsosurface_Build( void* drawingObject, void* data ) {}
void _lucIsosurface_Initialise( void* drawingObject, void* data ) {}
void _lucIsosurface_Execute( void* drawingObject, void* data ) {}
void _lucIsosurface_Destroy( void* drawingObject, void* data ) {}

void _lucIsosurface_Setup( void* drawingObject, void* _context ) {
	lucIsosurface*           self            = (lucIsosurface*)drawingObject;
	DomainContext*   context         = (DomainContext*) _context;
	FieldVariable*           isosurfaceField = self->isosurfaceField;
	int                      i, j, k;
	int                      nx, ny, nz;
	double                   dx, dy, dz;
	Vertex***                vertex;
	Coord                    pos;
	Coord                    min;
	Coord                    max;
	Dimension_Index          dim             = context->dim;

	lucOpenGLDrawingObject_SyncShadowValues( self, self->isosurfaceField );

	/* Initialise Variables */
	self->triangleCount = 0;
	nx = self->resolution[ I_AXIS ];
	ny = self->resolution[ J_AXIS ];
	nz = self->resolution[ K_AXIS ];
	
	FieldVariable_GetMinAndMaxLocalCoords( isosurfaceField, min, max );
	dx = (max[0] - min[0])/((double) nx - 1);
	dy = (max[1] - min[1])/((double) ny - 1);
	dz = (max[2] - min[2])/((double) nz - 1);
	if (dim == 2) {
		dz = 0.0;
		nz = 1;
	}
	
	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	vertex = Memory_Alloc_3DArray( Vertex , nx, ny, nz , "Vertex array" );

	/* Sample Field in in regular grid */
	for ( i = 0 ; i < nx ; i++ ) {
		for ( j = 0 ; j < ny ; j++ ) {
			for ( k = 0 ; k < nz ; k++ ) {
				pos[ I_AXIS ] = min[ I_AXIS ] + dx * (double) i;
				pos[ J_AXIS ] = min[ J_AXIS ] + dy * (double) j;
				pos[ K_AXIS ] = min[ K_AXIS ] + dz * (double) k;

				memcpy( vertex[i][j][k].pos, pos, 3 * sizeof(double) );
				
				if ( i == 0 )
					pos[ I_AXIS ] = min[ I_AXIS ] + 0.001 * dx; 
				if ( j == 0 )
					pos[ J_AXIS ] = min[ J_AXIS ] + 0.001 * dy; 
				if ( k == 0 )
					pos[ K_AXIS ] = min[ K_AXIS ] + 0.001 * dz; 
				if ( i == nx - 1 )
					pos[ I_AXIS ] = max[ I_AXIS ] - 0.001 * dx;
				if ( j == ny - 1 )
					pos[ J_AXIS ] = max[ J_AXIS ] - 0.001 * dy;
				if ( k == nz - 1 )
					pos[ K_AXIS ] = max[ K_AXIS ] - 0.001 * dz;
				FieldVariable_InterpolateValueAt( isosurfaceField, pos, &vertex[i][j][k].value );
			}
		}
	}
	if (dim == 3) {
		/* Find normals */
		lucIsosurface_Normals( self, vertex );
		
		/* Find Surface with Marching Cubes */
		lucIsosurface_MarchingCubes( self, vertex );
	}
	if (dim == 2 || self->drawWalls)
		lucIsosurface_DrawWalls( self, vertex );

	/* Free memory */
	Memory_Free( vertex );

	/* Call parents setup function */
	_lucOpenGLDrawingObject_Setup( self, _context );
}



	#define lucIsosurface_GetColourForPos( self, pos, min, max, fudgeFactor ) \
		do { \
			double          __value; \
			Coord           __interpolationPos; \
			if ( self->colourField && self->colourMap ) { \
				memcpy( __interpolationPos, pos, sizeof(Coord) ); \
				if ( __interpolationPos[ I_AXIS ] >= max[ I_AXIS ] ) \
					__interpolationPos[ I_AXIS ] = max[ I_AXIS ] - fudgeFactor[ I_AXIS ]; \
				if ( __interpolationPos[ J_AXIS ] >= max[ J_AXIS ] ) \
					__interpolationPos[ J_AXIS ] = max[ J_AXIS ] - fudgeFactor[ J_AXIS ]; \
				if ( __interpolationPos[ K_AXIS ] >= max[ K_AXIS ] ) \
					__interpolationPos[ K_AXIS ] = max[ K_AXIS ] - fudgeFactor[ K_AXIS ]; \
				FieldVariable_InterpolateValueAt( self->colourField, __interpolationPos, &__value ); \
				lucColourMap_SetOpenGLColourFromValue( self->colourMap, __value ); \
			} \
		} while (0)
	
void _lucIsosurface_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	lucIsosurface*           self          = (lucIsosurface*)drawingObject;
	FieldVariable*           colourField   = self->colourField;
	lucColourMap*            colourMap     = self->colourMap;
	Coord                    min;
	Coord                    max;
	XYZ                      fudgeFactor   = { 0.0, 0.0, 0.0 };
	Index                    dof_I;

	lucOpenGLDrawingObject_SyncShadowValues( self, self->isosurfaceField );
	lucOpenGLDrawingObject_SyncShadowValues( self, colourField );
	lucOpenGLDrawingObject_SyncShadowValues( self, self->maskField );

	/* Calibrate Colour Map using Colour Variable */
	if ( colourMap && colourField ) {
		lucColourMap_CalibrateFromFieldVariable( colourMap, colourField );
		FieldVariable_GetMinAndMaxGlobalCoords( colourField, min, max );
		for ( dof_I = 0 ; dof_I < 3 ; dof_I++ )
			fudgeFactor[ dof_I ] = 0.001 * (max[ dof_I ] - min[ dof_I ]);
		FieldVariable_GetMinAndMaxLocalCoords( colourField, min, max );
	}
	else if ( colourMap ) 
		lucColourMap_CalibrateFromFieldVariable( colourMap, self->isosurfaceField );	
	glEnable(GL_LIGHTING);
	_lucOpenGLDrawingObject_Draw( self, window, viewportInfo, _context );
}


void _lucIsosurface_CleanUp( void* drawingObject, void* _context ) {
	lucIsosurface*            self            = (lucIsosurface*) drawingObject;
	
	_lucOpenGLDrawingObject_CleanUp( self, _context );
}

void _lucIsosurface_BuildDisplayList( void* drawingObject, void* _context ) {
	lucIsosurface*           self          = (lucIsosurface*)drawingObject;
	FieldVariable*           colourField   = self->colourField;
	FieldVariable*           maskField     = self->maskField;
	Surface_Triangle*        triangleList  = self->triangleList;
	Surface_Triangle*        currentTriangle;
	Index                    triangle_I;
	Index                    triangleCount = self->triangleCount;
	lucColourMap*            colourMap     = self->colourMap;
	Coord                    min;
	Coord                    max;
	XYZ                      fudgeFactor   = { 0.0, 0.0, 0.0 };
	Index                    dof_I;

	lucOpenGLDrawingObject_SyncShadowValues( self, self->isosurfaceField );
	lucOpenGLDrawingObject_SyncShadowValues( self, colourField );
	lucOpenGLDrawingObject_SyncShadowValues( self, maskField );

	/* Give option to draw surface as wireframe */
	if (self->wireframe) 
		glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
	else 
		glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );	

	/* Set up Face Culling */
	if (self->cullBackFace) {
		if (self->cullFrontFace) 
			glCullFace( GL_FRONT_AND_BACK );
		else 
			glCullFace( GL_BACK );
		glEnable( GL_CULL_FACE );
	}
	else if( self->cullFrontFace ) {
		glCullFace( GL_FRONT );
		glEnable( GL_CULL_FACE );
	}

	/* Calibrate Colour Map using Colour Variable */
	if ( colourMap && colourField ) {
		lucColourMap_CalibrateFromFieldVariable( colourMap, colourField );
		FieldVariable_GetMinAndMaxGlobalCoords( colourField, min, max );
		for ( dof_I = 0 ; dof_I < 3 ; dof_I++ )
			fudgeFactor[ dof_I ] = 0.001 * (max[ dof_I ] - min[ dof_I ]);
		FieldVariable_GetMinAndMaxLocalCoords( colourField, min, max );
	}
	else if ( colourMap ) 
		lucColourMap_CalibrateFromFieldVariable( colourMap, self->isosurfaceField );


	glFrontFace(GL_CW);
	glBegin(GL_TRIANGLES);

	if ( ! self->colourMap || ! colourField ) {
		lucColour_SetOpenGLColour( &self->colour );
	}
	else if ( !colourField ) {
		lucColourMap_SetOpenGLColourFromValue( colourMap, self->isovalue );
	}
	
	for( triangle_I = 0 ; triangle_I < triangleCount ; triangle_I++) {
		currentTriangle = (Surface_Triangle*) &triangleList[triangle_I];

		/* Test for masked triangle */
		if ( maskField != NULL ) {
			if (lucIsosurface_TestMask( self, currentTriangle->pos1 ) == False) continue;
			if (lucIsosurface_TestMask( self, currentTriangle->pos2 ) == False) continue;
			if (lucIsosurface_TestMask( self, currentTriangle->pos3 ) == False) continue;
		}

		/* Plot First Vertex */
		lucIsosurface_GetColourForPos( self, currentTriangle->pos1, min, max, fudgeFactor );
		glNormal3dv(currentTriangle->normal1);
		glVertex3dv(currentTriangle->pos1);
		
		/* Plot Second Vertex */
		lucIsosurface_GetColourForPos( self, currentTriangle->pos3, min, max, fudgeFactor );
		glNormal3dv(currentTriangle->normal3);
		glVertex3dv(currentTriangle->pos3);	
		
		/* Plot Third Vertex */
		lucIsosurface_GetColourForPos( self, currentTriangle->pos2, min, max, fudgeFactor );
		glNormal3dv(currentTriangle->normal2);
		glVertex3dv(currentTriangle->pos2);
	}
	glEnd();
}

Bool lucIsosurface_TestMask( lucIsosurface* self, Coord pos ) {
	double value;

	FieldVariable_InterpolateValueAt( self->maskField, pos, &value );

	return lucDrawingObjectMask_Test( &self->mask, value );
}

/*
   Given a grid cell and an isolevel, calculate the triangular
   facets required to represent the isosurface through the cell.
   Return the number of triangular facets, the array "triangleList"
   will be loaded up with the vertices at most 5 triangular facets.
	0 will be returned if the grid cell is either totally above
   of totally below the isolevel.
*/

void VertexInterp(double isolevel, Vertex* point, Vertex* vertex1, Vertex* vertex2 ) ;
double getQuadraticDerive(double y1, double y2, double y3, double x1, double x2, double x3) ;

/* This algorithm for constructing an isosurface is taken from:
Lorensen, William and Harvey E. Cline. Marching Cubes: A High Resolution 3D Surface Construction Algorithm. Computer Graphics (SIGGRAPH 87 Proceedings) 21(4) July 1987, p. 163-170) http://www.cs.duke.edu/education/courses/fall01/cps124/resources/p163-lorensen.pdf 
The lookup table is taken from http://astronomy.swin.edu.au/~pbourke/modelling/polygonise/ 
*/
void lucIsosurface_MarchingCubes( lucIsosurface* self, Vertex*** vertex ) {
	double isolevel = self->isovalue;
	int i,j,k,n;
	int triangleCount = 0;
	int nx = self->resolution[ I_AXIS ];
	int ny = self->resolution[ J_AXIS ];
	int nz = self->resolution[ K_AXIS ];
	int cubeindex;
	Vertex points[12];

int edgeTable[256]={
0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };
int triTable[256][16] =
{{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
{3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
{3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
{3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
{9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
{9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
{2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
{8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
{9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
{4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
{3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
{1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
{4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
{4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
{5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
{2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
{9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
{0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
{2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
{10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
{5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
{5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
{9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
{0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
{1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
{10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
{8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
{2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
{7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
{2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
{11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
{5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
{11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
{11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
{1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
{9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
{5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
{2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
{5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
{6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
{3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
{6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
{5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
{1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
{10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
{6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
{8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
{7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
{3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
{5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
{0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
{9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
{8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
{0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
{10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
{10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
{1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
{0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
{10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
{3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
{9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
{3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
{6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
{10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
{10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
{7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
{7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
{1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
{8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
{0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
{7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
{10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
{2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
{7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
{2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
{10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
{0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
{7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
{6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
{8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
{6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
{4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
{8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
{0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
{1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
{8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
{10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
{4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
{10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
{5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
{9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
{7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
{7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
{3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
{9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
{4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
{6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
{0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
{6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
{11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
{5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
{9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
{1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
{1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
{0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
{5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
{11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
{9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
{2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
{9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
{1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
{9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
{9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
{0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
{2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
{0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
{9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
{5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
{3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
{8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
{0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
{9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
{1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
{4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
{11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
{2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
{3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
{1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
{4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
{3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
{0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
{9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
{1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};


	for ( i = 0 ; i < nx - 1  ; i++ ) {
		for ( j = 0 ; j < ny - 1 ; j++ ) {
			for ( k = 0 ; k < nz - 1 ; k++ ) {
				/* Determine the index into the edge table which tells us which vertices are inside of the surface */
				cubeindex = 0;
				if (vertex[i][j][k].value       < isolevel) cubeindex |= 1;	
				if (vertex[i+1][j][k].value     < isolevel) cubeindex |= 2;
				if (vertex[i+1][j][k+1].value   < isolevel) cubeindex |= 4;
				if (vertex[i][j][k+1].value     < isolevel) cubeindex |= 8;
				if (vertex[i][j+1][k].value     < isolevel) cubeindex |= 16;
				if (vertex[i+1][j+1][k].value   < isolevel) cubeindex |= 32;
				if (vertex[i+1][j+1][k+1].value < isolevel) cubeindex |= 64;
				if (vertex[i][j+1][k+1].value   < isolevel) cubeindex |= 128;
				
				/* Cube is entirely in/out of the surface */
				if (edgeTable[cubeindex] == 0) continue;
					
				/* Find the vertices where the surface intersects the cube */
				if (edgeTable[cubeindex] & 1)
					VertexInterp(isolevel, &points[0], &vertex[i][j][k] , &vertex[i+1][j][k]);
				if (edgeTable[cubeindex] & 2)
					VertexInterp(isolevel, &points[1], &vertex[i+1][j][k] , &vertex[i+1][j][k+1] );
				if (edgeTable[cubeindex] & 4)
					VertexInterp(isolevel, &points[2], &vertex[i+1][j][k+1] , &vertex[i][j][k+1] );
				if (edgeTable[cubeindex] & 8)
					VertexInterp(isolevel, &points[3], &vertex[i][j][k+1] , &vertex[i][j][k] );
				if (edgeTable[cubeindex] & 16)
					VertexInterp(isolevel, &points[4], &vertex[i][j+1][k] , &vertex[i+1][j+1][k] );
				if (edgeTable[cubeindex] & 32)
					VertexInterp(isolevel, &points[5], &vertex[i+1][j+1][k] , &vertex[i+1][j+1][k+1] );
				if (edgeTable[cubeindex] & 64)
					VertexInterp(isolevel, &points[6], &vertex[i+1][j+1][k+1] , &vertex[i][j+1][k+1] );
				if (edgeTable[cubeindex] & 128)
					VertexInterp(isolevel, &points[7], &vertex[i][j+1][k+1] , &vertex[i][j+1][k] );
				if (edgeTable[cubeindex] & 256)
					VertexInterp(isolevel, &points[8], &vertex[i][j][k] , &vertex[i][j+1][k] );
				if (edgeTable[cubeindex] & 512)
					VertexInterp(isolevel, &points[9], &vertex[i+1][j][k] , &vertex[i+1][j+1][k] );
				if (edgeTable[cubeindex] & 1024)
					VertexInterp(isolevel, &points[10], &vertex[i+1][j][k+1] , &vertex[i+1][j+1][k+1] );
				if (edgeTable[cubeindex] & 2048)
					VertexInterp(isolevel, &points[11], &vertex[i][j][k+1] , &vertex[i][j+1][k+1] );

				/* Create the triangle */
				for ( n = 0 ; triTable[cubeindex][n] != -1 ; n += 3 ) {
					if (triangleCount > self->trianglesAlloced - 1 ) {
						self->trianglesAlloced = self->trianglesAlloced + 100;
						self->triangleList = Memory_Realloc_Array( self->triangleList, Surface_Triangle, self->trianglesAlloced );
					}
					/* Get positions */
					memcpy( self->triangleList[triangleCount].pos1 , points[ triTable[cubeindex][n  ] ].pos , 3 * sizeof(double) );
					memcpy( self->triangleList[triangleCount].pos2 , points[ triTable[cubeindex][n+1] ].pos , 3 * sizeof(double) );
					memcpy( self->triangleList[triangleCount].pos3 , points[ triTable[cubeindex][n+2] ].pos , 3 * sizeof(double) );

					/* Get Normals */
					memcpy( self->triangleList[triangleCount].normal1 , points[ triTable[cubeindex][n  ] ].normal , 3 * sizeof(double) );
					memcpy( self->triangleList[triangleCount].normal2 , points[ triTable[cubeindex][n+1] ].normal , 3 * sizeof(double) );
					memcpy( self->triangleList[triangleCount].normal3 , points[ triTable[cubeindex][n+2] ].normal , 3 * sizeof(double) );

					triangleCount++;
				}
			}
		}
	}
	self->triangleCount = triangleCount;
}

/*
   Linearly interpolate the position where an isosurface cuts
   an edge between two vertices, each with their own scalar value
*/
void VertexInterp(double isolevel, Vertex* point, Vertex* vertex1, Vertex* vertex2 ) {
	double mu;

	if (fabs(isolevel - vertex1->value) < 0.00001) {
		memcpy( point, vertex1, sizeof(Vertex) );
		return;
	}
	if (fabs(isolevel - vertex2->value) < 0.00001) {
		memcpy( point, vertex2, sizeof(Vertex) );
		return;
	}
	if (fabs(vertex1->value - vertex2->value) < 0.00001) {
		memcpy( point, vertex1, sizeof(Vertex) );
		return;
	}
		
	mu = (isolevel - vertex1->value) / (vertex2->value - vertex1->value);
	
	point->pos[0] = vertex1->pos[0] + mu * (vertex2->pos[0] - vertex1->pos[0]);
	point->pos[1] = vertex1->pos[1] + mu * (vertex2->pos[1] - vertex1->pos[1]);
	point->pos[2] = vertex1->pos[2] + mu * (vertex2->pos[2] - vertex1->pos[2]);

	point->normal[0] = vertex1->normal[0] + mu * (vertex2->normal[0] - vertex1->normal[0]);
	point->normal[1] = vertex1->normal[1] + mu * (vertex2->normal[1] - vertex1->normal[1]);
	point->normal[2] = vertex1->normal[2] + mu * (vertex2->normal[2] - vertex1->normal[2]);
}

void lucIsosurface_Normals( lucIsosurface* self, Vertex*** vertex ) {
	int i, j, k;
	int nx = self->resolution[ I_AXIS ];
	int ny = self->resolution[ J_AXIS ];
	int nz = self->resolution[ K_AXIS ];

	for ( i = 0 ; i < nx ; i++ ) {
		for ( j = 0 ; j < ny ; j++ ) {
			for ( k = 0 ; k < nz ; k++ ) {
				/* Set up component in x direction */
				if (i == 0) {
					vertex[i][j][k].normal[0] = (vertex[i+1][j][k].value - vertex[i][j][k].value) / 
									(vertex[i+1][j][k].pos[0] - vertex[i][j][k].pos[0] ) ;
				}
				else if (i == nx - 1)  {
					vertex[i][j][k].normal[0] = (vertex[i][j][k].value - vertex[i-1][j][k].value) / 
									(vertex[i][j][k].pos[0] - vertex[i-1][j][k].pos[0] ) ;
				}
				else {	
					vertex[i][j][k].normal[0] = getQuadraticDerive(	
					                              vertex[i-1][j][k].value, vertex[i][j][k].value, vertex[i+1][j][k].value,
					                              vertex[i-1][j][k].pos[0], vertex[i][j][k].pos[0], vertex[i+1][j][k].pos[0]);
				}
				
				/* Set up component in y direction */
				if (j == 0) {
					vertex[i][j][k].normal[1] = (vertex[i][j+1][k].value - vertex[i][j][k].value) / 
									(vertex[i][j+1][k].pos[1] - vertex[i][j][k].pos[1] ) ;
				}
				else if (j == ny - 1)  {
					vertex[i][j][k].normal[1] = (vertex[i][j][k].value - vertex[i][j-1][k].value) / 
									(vertex[i][j][k].pos[1] - vertex[i][j-1][k].pos[1] ) ;
				}
				else {	
					vertex[i][j][k].normal[1] = getQuadraticDerive(	
					                              vertex[i][j-1][k].value, vertex[i][j][k].value, vertex[i][j+1][k].value,
					                              vertex[i][j-1][k].pos[1], vertex[i][j][k].pos[1], vertex[i][j+1][k].pos[1]);
				}			
				
				/* Set up component in z direction */
				if (k == 0) {
					vertex[i][j][k].normal[2] = (vertex[i][j][k+1].value - vertex[i][j][k].value) / 
									(vertex[i][j][k+1].pos[2] - vertex[i][j][k].pos[2] ) ;
				}
				else if (k == nz - 1)  {
					vertex[i][j][k].normal[2] = (vertex[i][j][k].value - vertex[i][j][k-1].value) / 
									(vertex[i][j][k].pos[2] - vertex[i][j][k-1].pos[2] ) ;
				}
				else {	
					vertex[i][j][k].normal[2] = getQuadraticDerive(	
					                              vertex[i][j][k-1].value, vertex[i][j][k].value, vertex[i][j][k+1].value,
					                              vertex[i][j][k-1].pos[2], vertex[i][j][k].pos[2], vertex[i][j][k+1].pos[2]);
				}														
				StGermain_VectorNormalise(vertex[i][j][k].normal, 3);	
			}
		}
	}
}

double getQuadraticDerive(double y1, double y2, double y3, double x1, double x2, double x3) {
	double dy12, dy32;
	double dx12sq, dx32sq, dx12, dx32;
	double a, b;
	double result;

	if (fabs(y2 - y1) < 0.0001 && fabs(y3 - y1) < 0.0001) {
		return 0.0;
	}
	//  Get the problem variables
	dy12 = y1-y2;
	dy32 = y3-y2;
	dx12 = x1-x2;
	dx32 = x3-x2;
	dx12sq = (x1*x1)-(x2*x2);
	dx32sq = (x3*x3)-(x2*x2);
	
	//  Get a
	a = dx12 * dy32 / dx32;
	a -= dy12;
	a /= ((dx32sq/dx32) + dx12sq);
	
	//  Get b
	b = dy32;
	b -= a * dx32sq;
	b /= dx32;
	
	//  Now get the derivative
	result = 2*a*x2 + b;
	
	return result;
}


#define gLucifer_CW	0
#define gLucifer_CCW	1

#define LEFT_BOTTOM     0
#define RIGHT_BOTTOM    1
#define LEFT_TOP        2
#define RIGHT_TOP       3
#define LEFT            4
#define RIGHT           5
#define TOP	            6
#define BOTTOM          7

void lucIsosurface_DrawWalls( lucIsosurface* self, Vertex ***array ) {
	int nx = self->resolution[ I_AXIS ];
	int ny = self->resolution[ J_AXIS ];
	int nz = self->resolution[ K_AXIS ];
	int i, j, k;
	Vertex ** points;
	Vertex * midVerticies;
	char order;

	/* Allocate Memory */
	points = Memory_Alloc_Array( Vertex* , 8, "array for marching squares");
	midVerticies = Memory_Alloc_Array( Vertex , 4, "array for marching squares");
	points[LEFT] = &midVerticies[0];
	points[RIGHT] = &midVerticies[1];
	points[TOP] = &midVerticies[2];
	points[BOTTOM] = &midVerticies[3];
	
	for ( i = 0 ; i < nx - 1 ; i++ ) {
		for ( j = 0 ; j < ny - 1 ; j++ ) {
			k = 0;
			order = gLucifer_CCW;
			lucIsosurface_SetupPointsZ( points, array, i, j, k );
			lucIsosurface_WallElement( self, points, order, K_AXIS );

			k = nz - 1;
			order = gLucifer_CW;
			lucIsosurface_SetupPointsZ( points, array, i, j, k );
			lucIsosurface_WallElement( self, points, order, K_AXIS );
		}
	}
	for ( k = 0 ; k < nz - 1 ; k++ ) {
		for ( j = 0 ; j < ny - 1 ; j++ ) {
			i = 0;
			order = gLucifer_CW;
			lucIsosurface_SetupPointsX( points, array, i, j, k );
			lucIsosurface_WallElement( self, points, order, I_AXIS );

			i = nx - 1;
			order = gLucifer_CCW;
			lucIsosurface_SetupPointsX( points, array, i, j, k );
			lucIsosurface_WallElement( self, points, order, I_AXIS );

		}
	}
	for ( i = 0 ; i < nx - 1 ; i++ ) {
		for ( k = 0 ; k < nz - 1 ; k++ ) {
			j = 0;
			order = gLucifer_CW;
			lucIsosurface_SetupPointsY( points, array, i, j, k );
			lucIsosurface_WallElement( self, points, order, J_AXIS );

			j = ny - 1;
			order = gLucifer_CCW;
			lucIsosurface_SetupPointsY( points, array, i, j, k );
			lucIsosurface_WallElement( self, points, order, J_AXIS );
		}
	}
	Memory_Free( points );
	Memory_Free( midVerticies );
}

void lucIsosurface_SetupPointsX( Vertex** points, Vertex*** array, Index i, Index j, Index k ){
	points[LEFT_BOTTOM]  = &array[i][ j ][ k ];
	points[RIGHT_BOTTOM] = &array[i][j+1][ k ];
	points[LEFT_TOP]     = &array[i][ j ][k+1];
	points[RIGHT_TOP]    = &array[i][j+1][k+1];
}

void lucIsosurface_SetupPointsY( Vertex** points, Vertex*** array, Index i, Index j, Index k ){
	points[LEFT_BOTTOM]  = &array[ i ][j][ k ];
	points[RIGHT_BOTTOM] = &array[i+1][j][ k ];
	points[LEFT_TOP]     = &array[ i ][j][k+1];
	points[RIGHT_TOP]    = &array[i+1][j][k+1];
}

void lucIsosurface_SetupPointsZ( Vertex** points, Vertex*** array, Index i, Index j, Index k ){
	points[LEFT_BOTTOM]  = &array[ i ][ j ][k];
	points[RIGHT_BOTTOM] = &array[i+1][ j ][k];
	points[LEFT_TOP]     = &array[ i ][j+1][k];
	points[RIGHT_TOP]    = &array[i+1][j+1][k];
}

void lucIsosurface_WallElement( lucIsosurface* self, Vertex** points, char order, Dimension_Index axis ) {
	double value  = self->isovalue;
	char   cubeType = 0;

	/* find cube type */
	if (points[LEFT_BOTTOM]->value  > value) cubeType += 1;
	if (points[RIGHT_BOTTOM]->value > value) cubeType += 2;
	if (points[LEFT_TOP]->value     > value) cubeType += 4;
	if (points[RIGHT_TOP]->value    > value) cubeType += 8;

	/* Create Points */
	lucIsosurface_CreateIntermediatePoints( self, points, axis );
	
	lucIsosurface_MarchingRectangles( self, points, cubeType, order );
}	

void lucIsosurface_AddWallTriangle( lucIsosurface* self, int a , int b, int c, Vertex** points, char order) {
	int n = self->triangleCount;
	if ( n > self->trianglesAlloced - 1 ) {
		self->trianglesAlloced = self->trianglesAlloced + 100;
		self->triangleList = Memory_Realloc_Array( self->triangleList, Surface_Triangle, self->trianglesAlloced );
	}
	
	if (order == gLucifer_CCW) {
		memcpy( self->triangleList[n].pos1, points[a]->pos, 3*sizeof(double) );
		memcpy( self->triangleList[n].pos2, points[c]->pos, 3*sizeof(double) );
		memcpy( self->triangleList[n].pos3, points[b]->pos, 3*sizeof(double) );
	}
	else {
		memcpy( self->triangleList[n].pos1, points[a]->pos, 3*sizeof(double) );
		memcpy( self->triangleList[n].pos2, points[b]->pos, 3*sizeof(double) );
		memcpy( self->triangleList[n].pos3, points[c]->pos, 3*sizeof(double) );
	}

	/* Calculate Normal */
	StGermain_NormalToPlane( self->triangleList[n].normal1 , 
		self->triangleList[n].pos1, self->triangleList[n].pos2, self->triangleList[n].pos3 );

	memcpy( self->triangleList[n].normal2, self->triangleList[n].normal1 , 3 * sizeof(double) );
	memcpy( self->triangleList[n].normal3, self->triangleList[n].normal1 , 3 * sizeof(double) );
	self->triangleCount++;
	
}

void lucIsosurface_CreateIntermediatePoints( lucIsosurface* self, Vertex **points, Dimension_Index axis ) {
	double value = self->isovalue;
	Dimension_Index A = ( axis == I_AXIS ? J_AXIS : I_AXIS );
	Dimension_Index B = ( axis == K_AXIS ? J_AXIS : K_AXIS );
	double dA = points[RIGHT_TOP]->pos[A] - points[LEFT_BOTTOM]->pos[A];
	double dB = points[RIGHT_TOP]->pos[B] - points[LEFT_BOTTOM]->pos[B];

	memcpy( points[ LEFT ]->pos, points[ LEFT_BOTTOM ]->pos, 3*sizeof(double) );
	points[ LEFT ]->pos[B] += dB * (value - points[LEFT_BOTTOM]->value) 
				/ (points[LEFT_TOP]->value - points[LEFT_BOTTOM]->value ) ;
	
	memcpy( points[ RIGHT ]->pos, points[ RIGHT_BOTTOM ]->pos, 3*sizeof(double) );
	points[ RIGHT ]->pos[B] += dB * (value - points[RIGHT_BOTTOM]->value) 
				/ (points[RIGHT_TOP]->value - points[RIGHT_BOTTOM]->value ) ;
	
	memcpy( points[ BOTTOM ]->pos, points[ LEFT_BOTTOM ]->pos, 3*sizeof(double) );
	points[ BOTTOM ]->pos[A] += dA * (value - points[LEFT_BOTTOM]->value) 
				/ (points[RIGHT_BOTTOM]->value - points[LEFT_BOTTOM]->value ) ;
	
	memcpy( points[ TOP ]->pos, points[ LEFT_TOP ]->pos, 3*sizeof(double) );
	points[ TOP ]->pos[A] += dA * (value - points[LEFT_TOP]->value) 
				/ (points[RIGHT_TOP]->value - points[LEFT_TOP]->value ) ;
}		

void lucIsosurface_MarchingRectangles( lucIsosurface* self, Vertex** points, char cubeType, char order) {
	switch (cubeType) {
		case 0:
			/*  @@  */
			/*  @@  */
			break;
		case 1:		
			/*  @@  */
			/*  #@  */
			lucIsosurface_AddWallTriangle( self, LEFT_BOTTOM , LEFT, BOTTOM, points, order);
			break;
		case 2:
			/*  @@  */
			/*  @#  */	
			lucIsosurface_AddWallTriangle(self, BOTTOM , RIGHT, RIGHT_BOTTOM, points, order);
			break;
		case 3:
			/*  @@  */
			/*  ##  */	
			lucIsosurface_AddWallTriangle( self, LEFT_BOTTOM , LEFT, RIGHT, points, order);
			lucIsosurface_AddWallTriangle(self, LEFT_BOTTOM , RIGHT, RIGHT_BOTTOM, points, order);
			break;
		case 4:
			/*  #@  */
			/*  @@  */
			lucIsosurface_AddWallTriangle(self, LEFT_TOP , TOP, LEFT, points, order);
			break;
		case 5:
			/*  #@  */
			/*  #@  */
			lucIsosurface_AddWallTriangle(self, LEFT_TOP , TOP, LEFT_BOTTOM, points, order);
			lucIsosurface_AddWallTriangle(self, TOP , BOTTOM, LEFT_BOTTOM, points, order);
			break;
		case 6:
			/*  #@  */
			/*  @#  */
			lucIsosurface_AddWallTriangle(self, LEFT_TOP , TOP, LEFT, points, order);
			lucIsosurface_AddWallTriangle(self, BOTTOM , RIGHT, RIGHT_BOTTOM, points, order);
			break;
		case 7:
			/*  #@  */
			/*  ##  */
			lucIsosurface_AddWallTriangle(self, RIGHT , RIGHT_BOTTOM, TOP, points, order);
			lucIsosurface_AddWallTriangle(self, TOP , RIGHT_BOTTOM, LEFT_TOP, points, order);
			lucIsosurface_AddWallTriangle(self, LEFT_TOP , RIGHT_BOTTOM, LEFT_BOTTOM, points, order);
			break;
		case 8:
			/*  @#  */
			/*  @@  */
			lucIsosurface_AddWallTriangle(self, TOP , RIGHT_TOP, RIGHT, points, order);
			break;
		case 9:
			/*  @#  */
			/*  #@  */
			lucIsosurface_AddWallTriangle(self, TOP , RIGHT_TOP, RIGHT, points, order);
			lucIsosurface_AddWallTriangle(self, LEFT_BOTTOM , LEFT, BOTTOM, points, order);
			break;
		case 10:
			/*  @#  */
			/*  @#  */
			lucIsosurface_AddWallTriangle(self, TOP , RIGHT_TOP, RIGHT_BOTTOM, points, order);
			lucIsosurface_AddWallTriangle(self, BOTTOM , TOP, RIGHT_BOTTOM, points, order);
			
			break;
		case 11:
			/*  @#  */
			/*  ##  */
			lucIsosurface_AddWallTriangle(self, TOP , LEFT_BOTTOM, LEFT, points, order);
			lucIsosurface_AddWallTriangle(self, RIGHT_TOP , LEFT_BOTTOM, TOP, points, order);
			lucIsosurface_AddWallTriangle(self, RIGHT_BOTTOM , LEFT_BOTTOM, RIGHT_TOP, points, order);
			break;
		case 12:
			/*  ##  */
			/*  @@  */
			lucIsosurface_AddWallTriangle(self, LEFT , LEFT_TOP, RIGHT, points, order);
			lucIsosurface_AddWallTriangle(self, RIGHT , LEFT_TOP, RIGHT_TOP, points, order);
			break;
		case 13:
			/*  ##  */
			/*  #@  */
			lucIsosurface_AddWallTriangle(self, BOTTOM , RIGHT_TOP, RIGHT, points, order);
			lucIsosurface_AddWallTriangle(self, LEFT_BOTTOM , RIGHT_TOP, BOTTOM, points, order);
			lucIsosurface_AddWallTriangle(self, LEFT_TOP , RIGHT_TOP, LEFT_BOTTOM, points, order);
			break;
		case 14:
			/*  ##  */
			/*  @#  */
			lucIsosurface_AddWallTriangle(self, LEFT , LEFT_TOP, BOTTOM, points, order);
			lucIsosurface_AddWallTriangle(self, BOTTOM , LEFT_TOP, RIGHT_BOTTOM, points, order);
			lucIsosurface_AddWallTriangle(self, RIGHT_BOTTOM , LEFT_TOP, RIGHT_TOP, points, order);
			break;
		case 15:
			/*  ##  */
			/*  ##  */
			lucIsosurface_AddWallTriangle(self, LEFT_TOP , RIGHT_TOP, RIGHT_BOTTOM, points, order);
			lucIsosurface_AddWallTriangle(self, RIGHT_BOTTOM , LEFT_BOTTOM, LEFT_TOP, points, order);
			break;
		default:
			Journal_Printf( self->errorStream, "In func %s: Cannot understand cube type %d\n", __func__, cubeType );
			abort();
	}
}


