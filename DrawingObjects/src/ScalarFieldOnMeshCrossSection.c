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
** $Id: ScalarFieldCrossSection.c 568 2006-06-02 06:21:50Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#ifdef GLUCIFER_USE_PICELLERATOR
	#include <StgFEM/StgFEM.h>
	#include <PICellerator/PICellerator.h>
#endif

#include <glucifer/Base/Base.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>

#include "types.h"
#include "CrossSection.h"
#include "ScalarFieldOnMeshCrossSection.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>
#include <string.h>
#include <ctype.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucScalarFieldOnMeshCrossSection_Type = "lucScalarFieldOnMeshCrossSection";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucScalarFieldOnMeshCrossSection* _lucScalarFieldOnMeshCrossSection_New(  LUCSCALARFIELDONMESHCROSSSECTION_DEFARGS  ) 
{
	lucScalarFieldOnMeshCrossSection*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucScalarFieldOnMeshCrossSection) );
	self = (lucScalarFieldOnMeshCrossSection*) _lucCrossSection_New(  LUCCROSSSECTION_PASSARGS  );
	
	return self;
}

void _lucScalarFieldOnMeshCrossSection_Init( 
		lucScalarFieldOnMeshCrossSection*                            self,
		lucColourMap*                                                colourMap,
		XYZ                                                          minCropValues,
		XYZ                                                          maxCropValues,
      Bool                                                         wireFrame  ) 
{
	self->colourMap = colourMap;
	self->wireFrame = wireFrame;
	memcpy( self->minCropValues, minCropValues, sizeof(XYZ) );
	memcpy( self->maxCropValues, maxCropValues, sizeof(XYZ) );
}

void _lucScalarFieldOnMeshCrossSection_Delete( void* drawingObject ) {
	lucScalarFieldOnMeshCrossSection*  self = (lucScalarFieldOnMeshCrossSection*)drawingObject;

	_lucCrossSection_Delete( self );
}

void _lucScalarFieldOnMeshCrossSection_Print( void* drawingObject, Stream* stream ) {
	lucScalarFieldOnMeshCrossSection*  self = (lucScalarFieldOnMeshCrossSection*)drawingObject;

	_lucCrossSection_Print( self, stream );
}

void* _lucScalarFieldOnMeshCrossSection_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(lucScalarFieldOnMeshCrossSection);
	Type                                                             type = lucScalarFieldOnMeshCrossSection_Type;
	Stg_Class_DeleteFunction*                                     _delete = _lucScalarFieldOnMeshCrossSection_Delete;
	Stg_Class_PrintFunction*                                       _print = _lucScalarFieldOnMeshCrossSection_Print;
	Stg_Class_CopyFunction*                                         _copy = NULL;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _lucScalarFieldOnMeshCrossSection_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _lucScalarFieldOnMeshCrossSection_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _lucScalarFieldOnMeshCrossSection_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _lucScalarFieldOnMeshCrossSection_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _lucScalarFieldOnMeshCrossSection_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _lucScalarFieldOnMeshCrossSection_Destroy;
	lucDrawingObject_SetupFunction*                                _setup = _lucScalarFieldOnMeshCrossSection_Setup;
	lucDrawingObject_DrawFunction*                                  _draw = _lucCrossSection_Draw;
	lucDrawingObject_CleanUpFunction*                            _cleanUp = _lucOpenGLDrawingObject_CleanUp;
	lucOpenGLDrawingObject_BuildDisplayListFunction*    _buildDisplayList = _lucScalarFieldOnMeshCrossSection_BuildDisplayList;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _lucScalarFieldOnMeshCrossSection_New(  LUCSCALARFIELDONMESHCROSSSECTION_PASSARGS  );
}

void _lucScalarFieldOnMeshCrossSection_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucScalarFieldOnMeshCrossSection*     self = (lucScalarFieldOnMeshCrossSection*)drawingObject;
	lucColourMap*    colourMap;
	XYZ              minCropValues;
	XYZ              maxCropValues;
   Bool             wireFrame;

	/* Construct Parent */
	_lucCrossSection_AssignFromXML( self, cf, data );

	colourMap = Stg_ComponentFactory_ConstructByKey( cf, self->name, "ColourMap", lucColourMap, True, data );

	/* Get Values with which to crop the cross section */
	minCropValues[ I_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, "minCropX", -HUGE_VAL );
	minCropValues[ J_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, "minCropY", -HUGE_VAL );
	minCropValues[ K_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, "minCropZ", -HUGE_VAL );
	maxCropValues[ I_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, "maxCropX", +HUGE_VAL );
	maxCropValues[ J_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, "maxCropY", +HUGE_VAL );
	maxCropValues[ K_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, "maxCropZ", +HUGE_VAL );
	wireFrame = Stg_ComponentFactory_GetBool( cf, self->name, "wireFrame", False);
	
	_lucScalarFieldOnMeshCrossSection_Init( 
			self, 
			colourMap,
			minCropValues,
			maxCropValues,
         wireFrame );
}

void _lucScalarFieldOnMeshCrossSection_Build( void* drawingObject, void* data ) {
	lucScalarFieldOnMeshCrossSection*     self    = (lucScalarFieldOnMeshCrossSection*)drawingObject;
	FeVariable*                     feVariable;
	Mesh*                           mesh;
	Stream*                         errorStream = Journal_Register( Error_Type, self->type );
	
   /* Build field variable in parent */
   _lucCrossSection_Build(self, data);

	Journal_Firewall( self->fieldVariable->fieldComponentCount == 1, errorStream,
		"Error - in %s(): provided FieldVariable \"%s\" has %u components - but %s Component "
		"can only visualise FieldVariables with 1 component. Did you mean to visualise the "
		"magnitude of the given field?\n", __func__, self->fieldVariable->name,
		self->fieldVariable->fieldComponentCount, self->type );

   feVariable = (FeVariable*)self->fieldVariable;
	mesh = (Mesh*)feVariable->feMesh;

	/* Store the Vertex Grid */
	self->vertexGridHandle = ExtensionManager_GetHandle( mesh->info, "vertexGrid" );
	if ( self->vertexGridHandle == (ExtensionInfo_Index)-1 )

	Journal_Firewall( self->vertexGridHandle != (ExtensionInfo_Index)-1, errorStream,
		"Error - in %s(): provided FieldVariable \"%s\" doesn't have a Vertex Grid.\n"
		"Try visualising with lucScalarField instead.\n", __func__, self->fieldVariable->name );
		
}

void _lucScalarFieldOnMeshCrossSection_Initialise( void* drawingObject, void* data ) {}
void _lucScalarFieldOnMeshCrossSection_Execute( void* drawingObject, void* data ) {}
void _lucScalarFieldOnMeshCrossSection_Destroy( void* drawingObject, void* data ) {}

void _lucScalarFieldOnMeshCrossSection_Setup( void* drawingObject, void* _context ) {
	lucScalarFieldOnMeshCrossSection*       self            = (lucScalarFieldOnMeshCrossSection*)drawingObject;

	lucColourMap_CalibrateFromFieldVariable( self->colourMap, self->fieldVariable );
	_lucOpenGLDrawingObject_Setup( self, _context );
}
	
void _lucScalarFieldOnMeshCrossSection_BuildDisplayList( void* drawingObject, void* _context ) {
	lucScalarFieldOnMeshCrossSection*       self            = (lucScalarFieldOnMeshCrossSection*)drawingObject;
	lucScalarFieldOnMeshCrossSection_DrawCrossSection( self );
}

void lucScalarFieldOnMeshCrossSection_PlotColouredNode( void* drawingObject, MeshVertex* vert) {
	lucScalarFieldOnMeshCrossSection*       self            = (lucScalarFieldOnMeshCrossSection*)drawingObject;

	/* Get colour for vertex */
	lucColourMap_SetOpenGLColourFromValue( self->colourMap, vert->value);

   if (!self->wireFrame) glNormal3dv(vert->normal);

	/* Plot vertex */
	if ( self->fieldVariable->dim == 2 )
		glVertex2dv( vert->pos );
	else 
		glVertex3dv( vert->pos );
}

void lucScalarFieldOnMeshCrossSection_DrawCrossSection( void* drawingObject ) {
   lucScalarFieldOnMeshCrossSection*       self            = (lucScalarFieldOnMeshCrossSection*)drawingObject;
   FeVariable*          fieldVariable = (FeVariable*) self->fieldVariable;
   Mesh*                mesh          = (Mesh*) fieldVariable->feMesh;
   Node_LocalIndex      crossSection_I = self->value;
   Axis                 aAxis, bAxis;
   Grid*                vertGrid;
   IJK                  node_ijk;
   Node_GlobalIndex     node_gI;
   Node_DomainIndex     node_dI;
   Node_DomainIndex     nDomainNodes;
   int i,j;

	/* Get Axis Directions */
	aAxis = ( self->axis == I_AXIS ? J_AXIS : I_AXIS );
	bAxis = ( self->axis == K_AXIS ? J_AXIS : K_AXIS );
	
	vertGrid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, self->vertexGridHandle );
	
	Journal_DPrintfL( self->debugStream, 2, 
			"%s called on field %s, with axis of cross section as %d, crossSection_I as %d\n",
			__func__, fieldVariable->name, self->axis, crossSection_I );
	
	nDomainNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );

   /* Visible from both sides by default and wound clockwise */
   glDisable(GL_CULL_FACE);
   glFrontFace(GL_CW);
   if (self->wireFrame) {
      /* Edges only */
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	   glDisable(GL_LINE_SMOOTH);
	   glDisable(GL_LIGHTING);
   } else {
      /* Filled */
	   glEnable(GL_LIGHTING);
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
   }

	/* Get mesh cross section vertices */
   MeshVertex* vertices[ vertGrid->sizes[ aAxis ] ] [ vertGrid->sizes[ bAxis ] ];
	node_ijk[ self->axis ] = crossSection_I;
	for ( i = 0 ; i < vertGrid->sizes[ aAxis ]; i++ ) {
      node_ijk[ aAxis ] = i;
	   for ( j = 0 ; j < vertGrid->sizes[ bAxis ]; j++ ) {
         node_ijk[ bAxis ] = j;
			node_gI = Grid_Project( vertGrid, node_ijk );
			/* Get Local Node Index */
			if( !Mesh_GlobalToDomain( mesh, MT_VERTEX, node_gI, &node_dI ) || node_dI >= nDomainNodes ) {
            vertices[i][j] = NULL; /* Flag not on this proc */
            continue;
         }

         vertices[i][j] = Memory_Alloc(MeshVertex, "Mesh Vertex");
         FeVariable_GetValueAtNode( fieldVariable, node_dI, &vertices[i][j]->value  );
		   memcpy( &vertices[i][j]->pos, fieldVariable->feMesh->verts[node_dI], 3*sizeof(double) );
         vertices[i][j]->normal[0] = vertices[i][j]->normal[1] = vertices[i][j]->normal[2] = 0; /* Zero normal */
		}
	}

   /* Calc normals for irregular meshes per vertex by averaging four surrounding triangle normals */
   if (!self->wireFrame) { 
      for ( i = 0 ; i < vertGrid->sizes[ aAxis ]; i++ ) {
         for ( j = 0 ; j < vertGrid->sizes[ bAxis ]; j++ ) {
            /* Get sum of normal vectors */
            if (vertices[i][j] ==  NULL) continue;
            int n, num = 0;
            double tempnormal[3];

            if (i > 0) {
               if (j > 0 && vertices[i][j-1] != NULL && vertices[i-1][j] != NULL) {
                  /* Look back in both axis  \|  */
                  StGermain_NormalToPlane( tempnormal, vertices[i][j]->pos, vertices[i][j-1]->pos, vertices[i-1][j]->pos);
                  for (n=0; n<3; n++) vertices[i][j]->normal[n] += tempnormal[n];
                  num++;
               }

               if (j < vertGrid->sizes[ bAxis ] - 1 && vertices[i-1][j] != NULL && vertices[i][j+1] != NULL) {
                  /* Look back in aAxis, forward in bAxis  /|  */
                  StGermain_NormalToPlane( tempnormal, vertices[i][j]->pos, vertices[i-1][j]->pos, vertices[i][j+1]->pos);
                  for (n=0; n<3; n++) vertices[i][j]->normal[n] += tempnormal[n];
                  num++;
               }
            }

            if (i <  vertGrid->sizes[ aAxis ] - 1) {
               if (j > 0 && vertices[i+1][j] != NULL && vertices[i][j-1] != NULL) {
                  /* Look forward in aAxis, back in bAxis  |/  */
                  StGermain_NormalToPlane( tempnormal, vertices[i][j]->pos, vertices[i+1][j]->pos, vertices[i][j-1]->pos);
                  for (n=0; n<3; n++) vertices[i][j]->normal[n] += tempnormal[n];
                  num++;
               }

               if (j < vertGrid->sizes[ bAxis ] - 1 && vertices[i][j+1] != NULL && vertices[i+1][j] != NULL) {
                  /* Look forward both axis  |\  */
                  StGermain_NormalToPlane( tempnormal, vertices[i][j]->pos, vertices[i][j+1]->pos, vertices[i+1][j]->pos);
                  for (n=0; n<3; n++) vertices[i][j]->normal[n] += tempnormal[n];
                  num++;
               }
            }

            //StGermain_VectorNormalise(vertices[i][j]->normal, 3);
            for (n=0; n<3; n++) vertices[i][j]->normal[n] /= num;
         }
      }
   }

   /* Plot quad strip vertices */
	for ( i = 0 ; i < vertGrid->sizes[ aAxis ]; i++ ) {
		glBegin(GL_QUAD_STRIP);
	   for ( j = 0 ; j < vertGrid->sizes[ bAxis ]; j++ ) {
         if (vertices[i][j] !=  NULL) {
            if (i+1 < vertGrid->sizes[ aAxis ] && vertices[i+1][j] != NULL) {
               lucScalarFieldOnMeshCrossSection_PlotColouredNode( self, vertices[i][j]);
               lucScalarFieldOnMeshCrossSection_PlotColouredNode( self, vertices[i+1][j]);
			      /* TODO Cropping - implement cropping for all in generic CrossSection object? */
            }
            Memory_Free(vertices[i][j]);
         }
      }
		glEnd();
   }

   glEnable(GL_LIGHTING);
   glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
   glEnable(GL_CULL_FACE);
}

