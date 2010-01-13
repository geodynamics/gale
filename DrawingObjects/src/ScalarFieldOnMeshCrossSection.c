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
      Bool                                                         wireFrame, 
		Bool                                                         cullFace )
{
	self->colourMap = colourMap;
	self->cullFace = cullFace;
	self->wireFrame = wireFrame;
	memcpy( self->minCropValues, minCropValues, sizeof(XYZ) );
	memcpy( self->maxCropValues, maxCropValues, sizeof(XYZ) );
   self->flipNormals = False;
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
   Bool             wireFrame, cullFace;

	/* Construct Parent */
	_lucCrossSection_AssignFromXML( self, cf, data );

	colourMap = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"ColourMap", lucColourMap, True, data  );

	/* Get Values with which to crop the cross section */
	minCropValues[ I_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"minCropX", -HUGE_VAL  );
	minCropValues[ J_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"minCropY", -HUGE_VAL  );
	minCropValues[ K_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"minCropZ", -HUGE_VAL  );
	maxCropValues[ I_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"maxCropX", +HUGE_VAL  );
	maxCropValues[ J_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"maxCropY", +HUGE_VAL  );
	maxCropValues[ K_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"maxCropZ", +HUGE_VAL  );
	wireFrame = Stg_ComponentFactory_GetBool( cf, self->name, (Dictionary_Entry_Key)"wireFrame", False );
   cullFace = Stg_ComponentFactory_GetBool( cf, self->name, (Dictionary_Entry_Key)"cullFace", True  );
	
	_lucScalarFieldOnMeshCrossSection_Init( 
			self, 
			colourMap,
			minCropValues,
			maxCropValues,
         wireFrame,
         cullFace );
}

void _lucScalarFieldOnMeshCrossSection_Build( void* drawingObject, void* data ) {
	lucScalarFieldOnMeshCrossSection*     self    = (lucScalarFieldOnMeshCrossSection*)drawingObject;
	FeVariable*                     feVariable;
	Mesh*                           mesh;
	Stream*                         errorStream = Journal_Register( Error_Type, (Name)self->type  );
	
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
	self->vertexGridHandle = ExtensionManager_GetHandle( mesh->info, (Name)"vertexGrid" );
	if ( self->vertexGridHandle == (ExtensionInfo_Index)-1 )

	Journal_Firewall( self->vertexGridHandle != (ExtensionInfo_Index )-1, errorStream,
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
	lucScalarFieldOnMeshCrossSection_DrawCrossSection( self, GL_CCW);
}

void lucScalarFieldOnMeshCrossSection_PlotColouredNode( void* drawingObject, MeshVertex* vert) {
	lucScalarFieldOnMeshCrossSection*       self            = (lucScalarFieldOnMeshCrossSection*)drawingObject;

	/* Get colour for vertex */
	lucColourMap_SetOpenGLColourFromValue( self->colourMap, vert->value);

   if (!self->wireFrame && self->fieldVariable->dim == 3)
      glNormal3dv(vert->normal);

	/* Plot vertex */
	if ( self->fieldVariable->dim == 2 )
		glVertex2dv( vert->pos );
	else 
		glVertex3dv( vert->pos );
}

Bool lucMeshVertex_SumNormal( double normal[3], MeshVertex* v1, MeshVertex* v2, MeshVertex* v3, Bool flip) {
   /* Utility function for calculating per vertex normals, pass vertices of four triangles
    * surrounding the vertex you want to calculate the normal for and they will be summed */
   double tempnormal[3];
   if (v1 == NULL || v2 == NULL || v3 == NULL) return False;
   StGermain_NormalToPlane( tempnormal, v1->pos, v2->pos, v3->pos);
   if (flip == False) {
      normal[0] += tempnormal[0];
      normal[1] += tempnormal[1];
      normal[2] += tempnormal[2];
   } else {
      normal[0] -= tempnormal[0];
      normal[1] -= tempnormal[1];
      normal[2] -= tempnormal[2];
   }
   return True;
}

void lucScalarFieldOnMeshCrossSection_DrawCrossSection( void* drawingObject, int direction ) {
   lucScalarFieldOnMeshCrossSection*       self            = (lucScalarFieldOnMeshCrossSection*)drawingObject;
   FeVariable*          fieldVariable = (FeVariable*) self->fieldVariable;
   Mesh*                mesh          = (Mesh*) fieldVariable->feMesh;
   Node_LocalIndex      crossSection_I = self->value;
   Grid*                vertGrid;
   IJK                  node_ijk;
   Node_GlobalIndex     node_gI;
   Node_DomainIndex     node_dI;
   Node_DomainIndex     nDomainNodes;
   int i,j;

	vertGrid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, self->vertexGridHandle );
	
	Journal_DPrintfL( self->debugStream, 2, 
			"%s called on field %s, with axis of cross section as %d, crossSection_I as %d\n",
			__func__, fieldVariable->name, self->axis, crossSection_I );
	
	nDomainNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );

   /* set polygon face winding */
	glFrontFace(direction); 
   /* Visible from both sides? */
	if ( self->cullFace )
      glEnable(GL_CULL_FACE);
	else
      glDisable(GL_CULL_FACE);

   if (self->wireFrame) {
      /* Edges only */
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	   glDisable(GL_LINE_SMOOTH);
	   glDisable(GL_LIGHTING);
   } else {
      /* Filled */
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	   glEnable(GL_LIGHTING);
   }

	/* Get mesh cross section vertices */
   MeshVertex* vertices[ vertGrid->sizes[ self->axis1 ] ] [ vertGrid->sizes[ self->axis2 ] ];
	node_ijk[ self->axis ] = crossSection_I;
	for ( i = 0 ; i < vertGrid->sizes[ self->axis1 ]; i++ ) {
      node_ijk[ self->axis1 ] = i;
	   for ( j = 0 ; j < vertGrid->sizes[ self->axis2 ]; j++ ) {
         node_ijk[ self->axis2 ] = j;
			node_gI = Grid_Project( vertGrid, node_ijk );
			/* Get Local Node Index */
			if( !Mesh_GlobalToDomain( mesh, MT_VERTEX, node_gI, &node_dI ) || node_dI >= nDomainNodes ) {
            vertices[i][j] = NULL; /* Flag not on this proc */
            continue;
         }

         vertices[i][j] = Memory_Alloc(MeshVertex, "Mesh Vertex");
         FeVariable_GetValueAtNode( fieldVariable, node_dI, &vertices[i][j]->value  );
		   memcpy( &vertices[i][j]->pos, fieldVariable->feMesh->verts[node_dI], 3 * sizeof(double) );
         vertices[i][j]->normal[0] = vertices[i][j]->normal[1] = vertices[i][j]->normal[2] = 0; /* Zero normal */
		}
	}

   /* Calc normals for irregular meshes per vertex by averaging four surrounding triangle normals */
   if (!self->wireFrame && self->fieldVariable->dim == 3) { 
      for ( i = 0 ; i < vertGrid->sizes[ self->axis1 ]; i++ ) {
         for ( j = 0 ; j < vertGrid->sizes[ self->axis2 ]; j++ ) {
            /* Get sum of normal vectors */
            if (vertices[i][j] ==  NULL) continue;

            if (i > 0) {
               if (j > 0) 
                  /* Look back in both axis  \|  */
                  lucMeshVertex_SumNormal(vertices[i][j]->normal, vertices[i][j], vertices[i-1][j], vertices[i][j-1], self->flipNormals);

               if (j < vertGrid->sizes[ self->axis2 ] - 1) 
                  /* Look back in self->axis1, forward in self->axis2  /|  */
                  lucMeshVertex_SumNormal(vertices[i][j]->normal, vertices[i][j], vertices[i][j+1], vertices[i-1][j], self->flipNormals);
            }

            if (i <  vertGrid->sizes[ self->axis1 ] - 1) {
               if (j > 0) 
                  /* Look forward in self->axis1, back in self->axis2  |/  */
                  lucMeshVertex_SumNormal(vertices[i][j]->normal, vertices[i][j], vertices[i][j-1], vertices[i+1][j], self->flipNormals);

               if (j < vertGrid->sizes[ self->axis2 ] - 1) 
                  /* Look forward both axis  |\  */
                  lucMeshVertex_SumNormal(vertices[i][j]->normal, vertices[i][j], vertices[i+1][j], vertices[i][j+1], self->flipNormals);
            }

            StGermain_VectorNormalise(vertices[i][j]->normal, 3);
         }
      }
   } else {
      /* Default plane normal */
      double normal[3]; 
      normal[self->axis1] = 0.0;
      normal[self->axis2] = 0.0;
      normal[self->axis] = 1.0;
      if (self->flipNormals == True) normal[self->axis] = -1.0;
      glNormal3dv(normal);
   }

   /* Plot quad strip vertices */
	for ( i = 0 ; i < vertGrid->sizes[ self->axis1 ]; i++ ) {
		glBegin(GL_QUAD_STRIP);
	   for ( j = 0 ; j < vertGrid->sizes[ self->axis2 ]; j++ ) {
         if (vertices[i][j] !=  NULL) {
            if (i+1 < vertGrid->sizes[ self->axis1 ] && vertices[i+1][j] != NULL) {
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
   glFrontFace(GL_CCW);

   /* Plot normals - for debugging /
	for ( i = 0 ; i < vertGrid->sizes[ self->axis1 ]; i++ ) {
	   for ( j = 0 ; j < vertGrid->sizes[ self->axis2 ]; j++ ) {
         if (vertices[i][j] !=  NULL) {
            luc_DrawNormalVector( vertices[i][j]->pos, vertices[i][j]->normal, 1000.0);
            Memory_Free(vertices[i][j]);
         }
      }
   }*/
}

