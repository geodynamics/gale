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
** $Id: Arrhenius.c 78 2005-11-29 11:58:21Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __lucMeshViewer_h__
#define __lucMeshViewer_h__

	typedef void (vertexFuncType)( double* );
	typedef void (lucMeshViewer_RenderEdgesFunc)( lucMeshViewer* self, 
						      vertexFuncType* vertexFunc );

	/* Textual name of this class - This is a global pointer which is used for 
	   times when you need to refer to class and not a particular instance of a class */
	extern const Type lucMeshViewer_Type;

	/* Class contents - this is defined as a macro so that sub-classes of 
	   this class can use this macro at the start of the definition of their struct */
	#define __lucMeshViewer							\
		/* Macro defining parent goes here - This means you can */	\
		/* cast this class as its parent */				\
		__lucOpenGLDrawingObject					\
		/* Virtual functions go here */					\
	        Mesh*					mesh;			\
		/* Other info */						\
		/* Colour Stuff */						\
	        lucColour				localColour;		\
		lucColour				shadowColour;		\
		lucColour				vacantColour;		\
	        /* Other Stuff */						\
		/* Stg_Class info */						\
		unsigned				nEdges;			\
		unsigned**				edges;			\
		lucMeshViewer_RenderEdgesFunc*		renderEdges;		\
										\
		Bool                    		nodeNumbers;		\
		Bool                    		elementNumbers;		\
		Bool                    		displayNodes;				\
		float                         lineWidth;        \
		Bool				skipXedges;			\
		Bool				skipYedges;			\
		Bool				skipZedges;			\

	struct lucMeshViewer { __lucMeshViewer };

	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define LUCMESHVIEWER_DEFARGS \
                LUCOPENGLDRAWINGOBJECT_DEFARGS

	#define LUCMESHVIEWER_PASSARGS \
                LUCOPENGLDRAWINGOBJECT_PASSARGS

	lucMeshViewer* _lucMeshViewer_New(  LUCMESHVIEWER_DEFARGS  );

	void _lucMeshViewer_Delete( void* drawingObject ) ;
	void _lucMeshViewer_Print( void* drawingObject, Stream* stream ) ;
	void* _lucMeshViewer_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) ;

	/* 'Stg_Component' implementations */
	void* _lucMeshViewer_DefaultNew( Name name ) ;
	void _lucMeshViewer_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data );
	void _lucMeshViewer_Build( void* drawingObject, void* data ) ;
	void _lucMeshViewer_Initialise( void* drawingObject, void* data ) ;
	void _lucMeshViewer_Execute( void* drawingObject, void* data );
	void _lucMeshViewer_Destroy( void* drawingObject, void* data ) ;
	
	void _lucMeshViewer_Setup( void* drawingObject, void* _context ) ;
	void _lucMeshViewer_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) ;
	void _lucMeshViewer_CleanUp( void* drawingObject, void* context ) ;

	void _lucMeshViewer_BuildDisplayList( void* drawingObject, void* _context ) ;

	void lucMeshViewer_UpdateVariables( void* drawingObject ) ;
	
	void lucMeshViewer_RenderGlobal( void* drawingObject );
	void lucMeshViewer_RenderLocal( void* drawingObject );
	void lucMeshViewer_RenderShadow( void* drawingObject, Partition_Index rank );
	void lucMeshViewer_RenderVacant( void* drawingObject, Partition_Index rank );
	void lucMeshViewer_Render( void* drawingObject );

        void lucMeshViewer_PrintNodeNumber( void* drawingObject, Coord coord, int* nodeNumber );
	void lucMeshViewer_PrintElementNumber( void* drawingObject, Coord coord, int* elementNumber );

	void lucMeshViewer_BuildEdges( lucMeshViewer* self );
	void lucMeshViewer_RenderEdges_WithInc( lucMeshViewer* self, vertexFuncType* vertexFunc );
	void lucMeshViewer_RenderEdges( lucMeshViewer* self, vertexFuncType* vertexFunc );
   Bool EdgeSkip(lucMeshViewer* self, double* v1, double* v2);

        void lucMeshViewer_PrintAllNodesNumber( void* drawingObject );

#endif

