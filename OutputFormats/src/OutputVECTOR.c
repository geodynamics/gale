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
** 
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#ifdef HAVE_GL2PS

#include <gl2ps.h>
#include <mpi.h>

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>

#include "types.h"
#include "OutputVECTOR.h"

#include <assert.h>
#include <string.h>

#ifndef MASTER
	#define MASTER 0
#endif

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucOutputVECTOR_Type = "lucOutputVECTOR";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucOutputVECTOR* _lucOutputVECTOR_New(  LUCOUTPUTVECTOR_DEFARGS  ) 
{
	lucOutputVECTOR*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucOutputVECTOR) );
	self = (lucOutputVECTOR*) _lucOutputFormat_New(  LUCOUTPUTFORMAT_PASSARGS  );
	
	return self;
}

void _lucOutputVECTOR_Init( lucOutputVECTOR* self, Stg_ComponentFactory* cf ){
	Name             formatName;
	Index            buffersize;
	
	formatName = Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"Format", "ps"  );

	if ( strcasecmp( formatName, "ps" ) == 0 ) {
		self->format            = "ps";
		self->gl2psFormatIndex = GL2PS_PS;
	}
	else if ( strcasecmp( formatName, "eps" ) == 0 ) {
		self->format            = "eps";
		self->gl2psFormatIndex = GL2PS_EPS;
	}
	else if ( strcasecmp( formatName, "svg" ) == 0 ) {
		self->format            = "svg";
		self->gl2psFormatIndex = GL2PS_SVG;
	}
	else if ( strcasecmp( formatName, "pdf" ) == 0 ) {
		self->format            = "pdf";
		self->gl2psFormatIndex = GL2PS_PDF;
	}
	else
	Journal_Firewall( 
			False,
			Journal_MyStream( Error_Type, self ),
			"\n Error:  Vector image output format '%s' appears to be incorrect or is unsupported. \n \n \
			 Supported formats are: \n \
			    svg  -  scalable vector graphics\n \
			    ps   -  postscript \n \
			    eps  -  encapsulated postscript \n \
			    pdf  -  portable document format \n \n", \
			formatName);	

	buffersize = Stg_ComponentFactory_GetInt( cf, self->name, (Dictionary_Entry_Key)"Buffersize", 4096*4096 );
	self->buffersize = buffersize;
}

void _lucOutputVECTOR_Delete( void* outputFormat ) {
	lucOutputVECTOR*  self = (lucOutputVECTOR*)outputFormat;

	_lucOutputFormat_Delete( self  );
}

void _lucOutputVECTOR_Print( void* outputFormat, Stream* stream ) {
	lucOutputVECTOR*  self = (lucOutputVECTOR*)outputFormat;

	_lucOutputFormat_Print( self, stream );
}

void* _lucOutputVECTOR_Copy( void* outputFormat, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucOutputVECTOR*  self = (lucOutputVECTOR*)outputFormat;
	lucOutputVECTOR* newOutputFormat;

	newOutputFormat = _lucOutputFormat_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newOutputFormat;
}


void* _lucOutputVECTOR_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(lucOutputVECTOR);
	Type                                                      type = lucOutputVECTOR_Type;
	Stg_Class_DeleteFunction*                              _delete = _lucOutputVECTOR_Delete;
	Stg_Class_PrintFunction*                                _print = _lucOutputVECTOR_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _lucOutputVECTOR_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _lucOutputVECTOR_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _lucOutputVECTOR_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _lucOutputVECTOR_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _lucOutputVECTOR_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _lucOutputVECTOR_Destroy;
	lucOutputFormat_OutputFunction*                        _output = _lucOutputVECTOR_Output;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _lucOutputVECTOR_New(  LUCOUTPUTVECTOR_PASSARGS  );
}

void _lucOutputVECTOR_AssignFromXML( void* outputFormat, Stg_ComponentFactory* cf, void* data ){
	lucOutputVECTOR*  self = (lucOutputVECTOR*)outputFormat;
	AbstractContext* context = Stg_ComponentFactory_ConstructByName( cf, (Name)"context", AbstractContext, True, data ) ;
	
	if(context->rank == MASTER )
		Journal_Firewall( context->nproc == 1, Journal_MyStream( Error_Type, self ), "\n \n     Vector outputting is not supported in parallel.\n     Please choose an alternate output format.\n\n");

	_lucOutputVECTOR_Init( self, cf );
	
	/* Construct Parent */
	_lucOutputFormat_AssignFromXML( outputFormat, cf, data);
	
}

void _lucOutputVECTOR_Build( void* outputFormat, void* data ) {}
void _lucOutputVECTOR_Initialise( void* outputFormat, void* data ) {}
void _lucOutputVECTOR_Execute( void* outputFormat, void* data ) {}
void _lucOutputVECTOR_Destroy( void* outputFormat, void* data ) {}

void _lucOutputVECTOR_Output( void* outputFormat, lucWindow* window, AbstractContext* context, void* pixelData ) {
	lucOutputVECTOR* self       = (lucOutputVECTOR*) outputFormat;
	Pixel_Index   width        = window->width;
	Pixel_Index   height       = window->height;
	/*FILE*         file         = fopen("test", "wb"); */
	FILE*         file         = lucOutputFormat_OpenFile( self, window, context, "wb");
	GLint         viewport[4];
	GLint         state;
	Name          filename;

	viewport[0] = 0;
	viewport[1] = 0;
	viewport[2] = width;
	viewport[3] = height;
	
	/* get filename for file headers */	
	filename = lucOutputFormat_GetImageFilename( self, window, context );
	/* remove directory from output filename */
	filename = (filename + strlen(context->outputPath) + 1);

	/* call to gl2ps which sets up glRenderMode(GL_FEEDBACK) and parses feedback to required format */
	gl2psBeginPage(filename, "gLucifer", viewport, self->gl2psFormatIndex, GL2PS_SIMPLE_SORT,GL2PS_NONE,  
		 GL_RGBA, 0, NULL, 0, 0, 0, self->buffersize, file, NULL);
	
	/* cleanup pre-existing scene (which was possibly used for current timestep raster images) */
	lucWindow_CleanUp( window, context );
	/* setup scene to be rendered again now that feedback mode is enabled */
	lucWindow_SetViewportNeedsToSetupFlag( window, True );
	lucWindow_SetViewportNeedsToDrawFlag( window, True );
	lucWindow_Draw( window, context );

	/* return to glRenderMode(GL_RENDER), and complete writing output file */
	state = gl2psEndPage();
	if(state == 5)
		Journal_Printf( Journal_MyStream( Error_Type, self ), "\nError. Insufficient GL feedback buffer size. \
								       \nConsider increasing the OutputVECTOR buffersize. \
								      \nVector image will not be created correctly.\n\n" );

	/* Clean Up */
	fclose(file);	

}

#endif /* HAVE_GL2PS */


