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
** $Id: EncoderLibfame.c 740 2007-10-11 08:05:31Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef HAVE_LIBFAME

#include <fame.h>

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>

#include "types.h"
#include "EncoderLibfame.h"

#include <assert.h>
#include <string.h>


#ifndef MASTER
	#define MASTER 0
#endif

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucEncoderLibfame_Type = "lucEncoderLibfame";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucEncoderLibfame* _lucEncoderLibfame_New(  LUCENCODERLIBFAME_DEFARGS  ) 
{
	lucEncoderLibfame*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucEncoderLibfame) );
	self = (lucEncoderLibfame*) _lucOutputFormat_New(  LUCOUTPUTFORMAT_PASSARGS  );
	
	return self;
}

void _lucEncoderLibfame_Init( 
		lucEncoderLibfame*                                                self,
		lucWindow*                                                        window,
		char*                                                             outputPath,
		Bool                                                             loadFromCheckPoint,
		unsigned int                                                      keyframe,
		unsigned int                                                      quality,
		unsigned int                                                      framesPerSecond,
		char*                                                             profile,
		Bool                                                              includeFrame0)
{

	Pixel_Index              width      = window->width;
	Pixel_Index              height     = window->height;
	unsigned int numpixels, numbytes, quarterpixels, linebytes;
	fame_context_t*          famecontext;
	fame_parameters_t*       fameparameters;
	fame_frame_statistics_t* famestatistics;
	fame_object_t*           fameprofile;
	fame_yuv_t*              fameyuv;
	Name                     filename;
	Index                    i;

   /* Output should only run on root node */
	if (self->context->rank != MASTER) return;

	/* Initialise the inlcudeFrame0 parameter */
	self->includeFrame0 = includeFrame0;

	/* Check to make sure width and height of window are multiples of 16 */
	Journal_Firewall( width % 16 == 0, lucError, "In func %s: Width %u is not multiple of 16.", __func__, width );
	Journal_Firewall( height % 16 == 0, lucError, "In func %s: Height %u is not multiple of 16.", __func__, height );

	lucOutputFormat_Register_Add( window->outputFormat_Register, self );

	/* Setup paramaters on context */
	numpixels = self->numpixels = width * height;
	numbytes = self->numbytes = self->numpixels * 3;
	linebytes = self->linebytes = width * 3;
	quarterpixels = self->quarterpixels = self->numpixels / 4;

	/* Allocate Memory */
	fameparameters = self->fameparameters = Memory_Alloc( fame_parameters_t, "Fame Parameters" );
	famestatistics = self->famestatistics = Memory_Alloc( fame_frame_statistics_t, "Fame Statistics" );
	fameyuv = self->fameyuv = Memory_Alloc( fame_yuv_t, "Fame YUV" );
	self->buffer = Memory_Alloc_Array( unsigned char , numbytes, "Buffer" );

	/* Set up YUV */
	fameyuv->w = width;
	fameyuv->h = height;
	fameyuv->p = 0;
	fameyuv->y = Memory_Alloc_Array( unsigned char, numpixels, "Y" );
	fameyuv->u = Memory_Alloc_Array( unsigned char, quarterpixels, "U" );
	fameyuv->v = Memory_Alloc_Array( unsigned char, quarterpixels, "V" );

	/* Setup keyframe */
	self->coding = Memory_Alloc_Array( char, keyframe + 1, "Coding" );
	sprintf( self->coding, "I" );
	for ( i = 1 ; i < keyframe ; i++ )
		strcat( self->coding, "P" );

	/* Set up Fame Parameters */
	memset(fameparameters, 0, sizeof(fame_parameters_t));
	fameparameters->width = width;
	fameparameters->height = height;
	fameparameters->coding = self->coding;
	fameparameters->quality = quality;
	fameparameters->bitrate = 0;
	fameparameters->slices_per_frame = 1;
	fameparameters->frames_per_sequence = 0xffffffff;
	fameparameters->frame_rate_num = framesPerSecond;
	fameparameters->frame_rate_den = 1;
	fameparameters->shape_quality = 100;
	fameparameters->search_range = 0;
	fameparameters->verbose = 0;
	fameparameters->profile = self->profile = StG_Strdup( profile );
	fameparameters->total_frames = 0;

	/* Open Output File */
	Stg_asprintf( &filename, "%s/%s.mpeg", outputPath, window->name );
	
	if(!loadFromCheckPoint ) 
	    self->stream = fopen( filename, "w" );
	else  {
		self->stream = fopen( filename, "a" );
	}
	Memory_Free( filename );

	/* Fame Initialisation */
	famecontext = self->famecontext = fame_open();
	fameprofile = self->fameprofile = fame_get_object(famecontext, fameparameters->profile );
	fame_register(famecontext, "profile", fameprofile);
	fame_init(famecontext, fameparameters, self->buffer, numbytes);
}

void _lucEncoderLibfame_Delete( void* outputFormat ) {
	lucEncoderLibfame*  self         = (lucEncoderLibfame*)outputFormat;
	unsigned int        framebytes;

   /* Output should only run on root node */
	if (self->context->rank == MASTER) {
      /* Finish writing mpeg and close file*/
      framebytes = fame_close(self->famecontext);
      fwrite(self->buffer, framebytes, 1, self->stream);
      fflush(self->stream);
      fclose(self->stream);

      /* Free Memory */
      Memory_Free( self->fameyuv->y );
      Memory_Free( self->fameyuv->u );
      Memory_Free( self->fameyuv->v );

      Memory_Free( self->fameyuv );
      Memory_Free( self->famestatistics );
      Memory_Free( self->fameparameters );
      Memory_Free( self->buffer );
      
      Memory_Free( self->coding );
      Memory_Free( self->profile );
   }

	_lucOutputFormat_Delete( self );
}

void _lucEncoderLibfame_Print( void* outputFormat, Stream* stream ) {
	lucEncoderLibfame*  self = (lucEncoderLibfame*)outputFormat;

	_lucOutputFormat_Print( self, stream );
}

void* _lucEncoderLibfame_Copy( void* outputFormat, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucEncoderLibfame*  self = (lucEncoderLibfame*)outputFormat;
	lucEncoderLibfame* newOutputFormat;

	newOutputFormat = _lucOutputFormat_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newOutputFormat;
}


void* _lucEncoderLibfame_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(lucEncoderLibfame);
	Type                                                      type = lucEncoderLibfame_Type;
	Stg_Class_DeleteFunction*                              _delete = _lucEncoderLibfame_Delete;
	Stg_Class_PrintFunction*                                _print = _lucEncoderLibfame_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _lucEncoderLibfame_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _lucEncoderLibfame_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _lucEncoderLibfame_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _lucEncoderLibfame_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _lucEncoderLibfame_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _lucEncoderLibfame_Destroy;
	lucOutputFormat_OutputFunction*                        _output = _lucEncoderLibfame_Output;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _lucEncoderLibfame_New(  LUCENCODERLIBFAME_PASSARGS  );
}

void _lucEncoderLibfame_AssignFromXML( void* outputFormat, Stg_ComponentFactory* cf, void* data ){
	lucEncoderLibfame*  self = (lucEncoderLibfame*)outputFormat;
	lucWindow*          window;
	AbstractContext*    context;

	/* Construct Parent */
	_lucOutputFormat_AssignFromXML( outputFormat, cf, data);
	lucOutputFormat_InitAll( self, "mpeg" );

	window =  Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Window", lucWindow, True, data   ) ;
	context = Stg_ComponentFactory_ConstructByName( cf, (Name)"context", AbstractContext, True, data  ) ;

	_lucEncoderLibfame_Init( 
			self,
			window,
			context->outputPath,
			context->loadFromCheckPoint,
			Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"keyframe", 4  ),
			Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"quality", 93  ),
			Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"framesPerSecond", 30  ),
			Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"profile", "profile/mpeg1"  ),
			Stg_ComponentFactory_GetBool( cf, self->name, (Dictionary_Entry_Key)"includeFrame0", False)  );
}

void _lucEncoderLibfame_Build( void* outputFormat, void* data ) {}
void _lucEncoderLibfame_Initialise( void* outputFormat, void* data ) {}
void _lucEncoderLibfame_Execute( void* outputFormat, void* data ) {}
void _lucEncoderLibfame_Destroy( void* outputFormat, void* data ) {}

void _lucEncoderLibfame_Output( void* outputFormat, lucWindow* window, AbstractContext* context, lucPixel* pixelData ) {
	lucEncoderLibfame*            self            = (lucEncoderLibfame*) outputFormat;
	Pixel_Index                   width           = window->width;
	unsigned int                  quarterpixels   = self->quarterpixels;
	unsigned int                  numpixels       = self->numpixels;
	unsigned int                  numbytes        = self->numbytes;
	unsigned int                  linebytes       = self->linebytes;
	fame_context_t*               famecontext     = self->famecontext;
	fame_yuv_t*                   fameyuv         = self->fameyuv;
	fame_frame_statistics_t*      famestatistics  = self->famestatistics; 
	unsigned int                  i, j, k, tmp;
	int                           framebytes;
	unsigned char                 red, green, blue;
	unsigned char*                rgbframe        = (unsigned char*) pixelData;


    /* preparing yuv12 format */
    memset(fameyuv->y, numpixels, 0);
    memset(fameyuv->u, quarterpixels, 0);
    memset(fameyuv->v, quarterpixels, 0);

    i = 0;   //position in yuv.
    j = numpixels - (width*2);   //position in rgb.
    k = 0;   //position in scanline.
    while(i < quarterpixels) {
		tmp = j*3;
		
		red = (rgbframe[tmp] + rgbframe[tmp+3] + rgbframe[tmp+linebytes]
			+ rgbframe[tmp+linebytes+3]) / 4;
		green = (rgbframe[tmp+1] + rgbframe[tmp+4]
			+ rgbframe[tmp+linebytes+1] + rgbframe[tmp+linebytes+4]) / 4;
		blue = (rgbframe[tmp+2] + rgbframe[tmp+5]
			+ rgbframe[tmp+linebytes+2] + rgbframe[tmp+linebytes+5]) / 4;
			
		fameyuv->u[i] = ((-38*red-74*green+112*blue+128)>>8)+128;
		fameyuv->v[i] = ((112*red-94*green-18*blue+128)>>8)+128;
		
		j += 2;
		k += 2;
		
		if(k == width) {
			k = 0;
			j -= (width*3);
		}
	
		i++;
	}

    i = 0;
    j = numpixels - width;
    k = 0;
    while(i < numbytes) {
		red = rgbframe[i];
		green = rgbframe[i+1];
		blue = rgbframe[i+2];
		
		fameyuv->y[j] = ((66*red+129*green+25*blue+128)>>8)+16;
		
		j++;
		k++;
		if(k == width) {
			k = 0;
			j -= (width*2);
		}
		i += 3;
	}
	
	/* Encode and Write frame 0 only if includeFrame0 is True and not in restart mode */

	if( context->timeStep == 0 ){
		if( self->includeFrame0 && (!context->loadFromCheckPoint) ) {	
			/* Initialise memory */
			memset(self->buffer,0, self->numbytes);
			memset(famestatistics, 0, sizeof(fame_frame_statistics_t));
			
			/* Encode */
			fame_start_frame(famecontext, fameyuv, NULL);
			framebytes = fame_encode_slice(famecontext);
			fame_end_frame(famecontext, famestatistics);
			
			/* Write encoded data to file */
			fwrite(self->buffer, framebytes, 1, self->stream);
			fflush(self->stream);
		}	
	}
	else{
		/* Initialise memory */
		memset(self->buffer,0, self->numbytes);
		memset(famestatistics, 0, sizeof(fame_frame_statistics_t));
		
		/* Encode */
		fame_start_frame(famecontext, fameyuv, NULL);
		framebytes = fame_encode_slice(famecontext);
		fame_end_frame(famecontext, famestatistics);
			
		/* Write encoded data to file */
		fwrite(self->buffer, framebytes, 1, self->stream);
		fflush(self->stream);
	}
}

#endif /* HAVE_LIBFAME */



