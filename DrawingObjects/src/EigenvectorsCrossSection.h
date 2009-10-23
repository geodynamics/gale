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
** $Id: EigenvectorsCrossSection.h 628 2006-10-12 08:23:07Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <glucifer/Base/CrossSection.h>

#ifndef __lucEigenvectorsCrossSection_h__
#define __lucEigenvectorsCrossSection_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type lucEigenvectorsCrossSection_Type;
		
	/** Class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __lucEigenvectorsCrossSection \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__lucOpenGLDrawingObject \
		/* Virtual functions go here */ \
		/* Other info */\
		FieldVariable*                                     tensorField;              \
		lucColour                                          colour[3];                \
		/* Colour used to display negative EigenValues */                            \
		lucColour                                          colourForNegative[3];     \
		IJK                                                resolution;               \
		double                                             arrowHeadSize;            \
		double                                             lengthScale;              \
		float                                              lineWidth;                \
      lucCrossSection*                                   crossSection;           \
		/* Specifies if the eigenvalue is used to draw the vector - default true */  \
		Bool 						   useEigenValue;            \
		/* Value used to draw the vector if the eigenvalue is not used */            \
		double 				                   notEigenValue;            \
		/* Specifies if the EigenVector and/or EigenValues are to be drawn */        \
		/* Default is True for EigenVector, False for EigenValues */                 \
                Bool                                               plotEigenVector;          \
		Bool                                               plotEigenValue;           \
		/* Used to scale the EigenValue if needed */                                 \
		double 					           scaleEigenValue;


	struct lucEigenvectorsCrossSection { __lucEigenvectorsCrossSection };
	
	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	lucEigenvectorsCrossSection* _lucEigenvectorsCrossSection_New( 
		SizeT                                              sizeOfSelf,
		Type                                               type,
		Stg_Class_DeleteFunction*                          _delete,
		Stg_Class_PrintFunction*                           _print,
		Stg_Class_CopyFunction*                            _copy, 
		Stg_Component_DefaultConstructorFunction*          _defaultConstructor,
		Stg_Component_ConstructFunction*                   _construct,
		Stg_Component_BuildFunction*                       _build,
		Stg_Component_InitialiseFunction*                  _initialise,
		Stg_Component_ExecuteFunction*                     _execute,
		Stg_Component_DestroyFunction*                     _destroy,
		lucDrawingObject_SetupFunction*                    _setup,
		lucDrawingObject_DrawFunction*                     _draw,
		lucDrawingObject_CleanUpFunction*                  _cleanUp,
		lucOpenGLDrawingObject_BuildDisplayListFunction*   _buildDisplayList,
		Name                                               name );

	void _lucEigenvectorsCrossSection_Delete( void* drawingObject ) ;
	void _lucEigenvectorsCrossSection_Print( void* drawingObject, Stream* stream ) ;
	void* _lucEigenvectorsCrossSection_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) ;

	/* 'Stg_Component' implementations */
	void* _lucEigenvectorsCrossSection_DefaultNew( Name name ) ;
	void _lucEigenvectorsCrossSection_Construct( void* drawingObject, Stg_ComponentFactory* cf, void* data );
	void _lucEigenvectorsCrossSection_Build( void* drawingObject, void* data ) ;
	void _lucEigenvectorsCrossSection_Initialise( void* drawingObject, void* data ) ;
	void _lucEigenvectorsCrossSection_Execute( void* drawingObject, void* data );
	void _lucEigenvectorsCrossSection_Destroy( void* drawingObject, void* data ) ;
	
	void _lucEigenvectorsCrossSection_Setup( void* drawingObject, void* _context ) ;
	void _lucEigenvectorsCrossSection_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) ;
	void _lucEigenvectorsCrossSection_CleanUp( void* drawingObject, void* _context ) ;

	void _lucEigenvectorsCrossSection_BuildDisplayList( void* drawingObject, void* _context ) ;
   void _lucEigenvectorsCrossSection_DrawCrossSection( void* drawingObject, Dimension_Index dim, lucCrossSection* crossSection );

#endif
