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
** $Id: EncoderLibavcodec.c 740 2007-10-11 08:05:31Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef HAVE_LIBAVCODEC

#include <ffmpeg/avcodec.h>

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>

#include "types.h"
#include "EncoderLibavcodec.h"

#include <assert.h>
#include <string.h>

#define MAX_BUFFER_SIZE 100000

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucEncoderLibavcodec_Type = "lucEncoderLibavcodec";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucEncoderLibavcodec* _lucEncoderLibavcodec_New(  LUCENCODERLIBAVCODEC_DEFARGS  ) 
{
	lucEncoderLibavcodec*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucEncoderLibavcodec) );
	self = (lucEncoderLibavcodec*) _lucOutputFormat_New(  LUCOUTPUTFORMAT_PASSARGS  );
	
	return self;
}

void _lucEncoderLibavcodec_Init( 
		lucEncoderLibavcodec*                                             self,
		lucWindow*                                                        window,
		char*                                                             outputPath,
		Bool                                                              loadFromCheckPoint,
		Name                                                              codecName,
		unsigned int                                                      framesPerSecond,
		unsigned int                                                      bitRate,
		Bool                                                              includeFrame0)
{
	Pixel_Index              width      = window->width;
	Pixel_Index              height     = window->height;
	Pixel_Index              pixelCount = width * height;
	Name                     filename;
	AVCodecContext*          codecContext;
	AVCodec*                 codec;
	AVFrame*                 frame;
	Stream*                  errorStream = Journal_MyStream( Error_Type, self );

	/*set the bool to display the Frame 0 or not */
	self->includeFrame0 = includeFrame0;

	/* Check to make sure width and height of window are multiples of 16 */
	Journal_Firewall( width % 2 == 0, lucError, "In func %s: Width %u is not multiple of 2.", __func__, width );
	Journal_Firewall( height % 2 == 0, lucError, "In func %s: Height %u is not multiple of 2.", __func__, height );

	lucOutputFormat_Register_Add( window->outputFormat_Register, self );

	/* Create codec context */
	codecContext = self->codecContext = avcodec_alloc_context();

	/* Setup paramaters on context */
	codecContext->bit_rate        = bitRate;
	codecContext->width           = width;
	codecContext->height          = height;
	
	/* Set the frame rate of the movie
	 * NB: libavcodec changed the way to define the frame rate in April, 2005 */
#if LIBAVCODEC_BUILD > 4753
	codecContext->time_base.num   = 1;
	codecContext->time_base.den   = (int) framesPerSecond;
#else
	codecContext->frame_rate      = (int) framesPerSecond;
	codecContext->frame_rate_base = 1;
#endif
	codecContext->gop_size        = 10;

	/* Get Codec from name */
	codec = avcodec_find_encoder_by_name( codecName );
	if ( codec == NULL ) {
		Journal_Printf( errorStream, "Error in func %s for %s '%s' - Couldn't find codec '%s'. Available codecs are:\n",
				__func__, self->type, self->name, codecName );

		/* Go through linked list of all the codecs avaiable and print them out. */
		Stream_Indent( errorStream );
		for ( codec = first_avcodec ; codec != NULL ; codec = codec->next ) {
			if ( codec->encode != NULL ) 
				Journal_Printf( errorStream, "%s (id = %d)\n", codec->name, codec->id );
		}
		abort();
	}

	/* Open it */
	if (avcodec_open(codecContext, codec) < 0) {
		Journal_Printf( errorStream, "Error in func %s for %s '%s' - Cannot open codec '%s'.\n", 
				__func__, self->type, self->name, codecName );
		abort();
	}

	/* Open Output File */
	Stg_asprintf( &filename, "%s/%s.mpeg", outputPath, window->name );
	self->stream = Journal_MyStream( Dump_Type, self );
	
	if(!loadFromCheckPoint)
		Stream_RedirectFile( self->stream, filename );
	else
		Stream_AppendFile( self->stream, filename );

	Memory_Free( filename );

	/* Create 'frame' data structure */
	frame = self->frame = avcodec_alloc_frame();
	frame->data[0] = Memory_Alloc_Array( unsigned char, pixelCount, "Y" );
	frame->data[1] = Memory_Alloc_Array( unsigned char, pixelCount/4, "Cr" );
	frame->data[2] = Memory_Alloc_Array( unsigned char, pixelCount/4, "Cb" );

	frame->linesize[0] = width;
	frame->linesize[1] = width / 2; 
	frame->linesize[2] = width / 2;

	/* Create output buffer which is what we write to the file for each frame */
	self->outputBuffer = Memory_Alloc_Bytes(MAX_BUFFER_SIZE,char, "outputBufferForFrame");
}

void _lucEncoderLibavcodec_Delete( void* outputFormat ) {
	lucEncoderLibavcodec*  self              = (lucEncoderLibavcodec*) outputFormat;
	AVFrame*               frame             = (AVFrame*) self->frame;
	uint8_t                endCodeSequence[] = { 0x00, 0x00, 0x01, 0xb7 };
	int                    sizeToWrite;


	/* Write the delayed frames to the file */
	do {
		sizeToWrite = avcodec_encode_video(self->codecContext, self->outputBuffer, MAX_BUFFER_SIZE, NULL);
		Journal_Write( self->stream, self->outputBuffer, sizeToWrite, 1 );
	} while ( sizeToWrite != 0 );
	
	/* Write the end code sequence to get a real mpeg file */
	Journal_Write( self->stream, endCodeSequence, 4, 1 );

	/* Free Memory */
	Memory_Free( frame->data[0] );
	Memory_Free( frame->data[1] );
	Memory_Free( frame->data[2] );
	Memory_Free( frame );
	Memory_Free( self->codecContext );
	Memory_Free( self->outputBuffer );

	_lucOutputFormat_Delete( self );
}

void _lucEncoderLibavcodec_Print( void* outputFormat, Stream* stream ) {
	lucEncoderLibavcodec*  self = (lucEncoderLibavcodec*) outputFormat;

	_lucOutputFormat_Print( self, stream );
}

void* _lucEncoderLibavcodec_Copy( void* outputFormat, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucEncoderLibavcodec*  self = (lucEncoderLibavcodec*)outputFormat;
	lucEncoderLibavcodec* newOutputFormat;

	newOutputFormat = _lucOutputFormat_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newOutputFormat;
}


void* _lucEncoderLibavcodec_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(lucEncoderLibavcodec);
	Type                                                      type = lucEncoderLibavcodec_Type;
	Stg_Class_DeleteFunction*                              _delete = _lucEncoderLibavcodec_Delete;
	Stg_Class_PrintFunction*                                _print = _lucEncoderLibavcodec_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _lucEncoderLibavcodec_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _lucEncoderLibavcodec_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _lucEncoderLibavcodec_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _lucEncoderLibavcodec_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _lucEncoderLibavcodec_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _lucEncoderLibavcodec_Destroy;
	lucOutputFormat_OutputFunction*                        _output = _lucEncoderLibavcodec_Output;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	return (void*) _lucEncoderLibavcodec_New(  LUCENCODERLIBAVCODEC_PASSARGS  );
}

void _lucEncoderLibavcodec_AssignFromXML( void* outputFormat, Stg_ComponentFactory* cf, void* data ){
	lucEncoderLibavcodec*  self = (lucEncoderLibavcodec*)outputFormat;
	lucWindow*          window;
	AbstractContext*    context;

	/* Construct Parent */
	lucOutputFormat_InitAll( self, "mpeg" );

	window =  Stg_ComponentFactory_ConstructByKey(  cf,  self->name,  "Window", lucWindow,  True, data  ) ;
	context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data ) ;
        
	_lucEncoderLibavcodec_Init( 
			self,
			window,
			context->outputPath,
			context->loadFromCheckPoint,
			Stg_ComponentFactory_GetString( cf, self->name, "codec", "mpeg1video" ),
			Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "frameRate", 25 ),
			Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "bitRate", 400000),
			Stg_ComponentFactory_GetBool( cf, self->name, "includeFrame0", False)
			);
}

void _lucEncoderLibavcodec_Build( void* outputFormat, void* data ) {}
void _lucEncoderLibavcodec_Initialise( void* outputFormat, void* data ) {}
void _lucEncoderLibavcodec_Execute( void* outputFormat, void* data ) {}
void _lucEncoderLibavcodec_Destroy( void* outputFormat, void* data ) {}

void _lucEncoderLibavcodec_Output( void* outputFormat, lucWindow* window, AbstractContext* context, lucPixel* pixelData ) {
	lucEncoderLibavcodec*         self            = (lucEncoderLibavcodec*) outputFormat;
	AVFrame*                      frame           = (AVFrame*) self->frame;
	Pixel_Index                   width           = window->width;
	Pixel_Index                   height          = window->height;
	Pixel_Index                   pixel_I;
	Pixel_Index                   xPixel_I;
	Pixel_Index                   yPixel_I;
	unsigned char*                macroPixel0;
	unsigned char*                macroPixel1;
	unsigned char*                macroPixel2;
	unsigned char*                macroPixel3;
	float                         red, green, blue;
	int                           sizeToWrite;

	/* Setup the 'Y' part of the frame */
	for ( yPixel_I = 0 ; yPixel_I < height ; yPixel_I++ ) {
		for ( xPixel_I = 0 ; xPixel_I < width ; xPixel_I++ ) {
			pixel_I = (height - yPixel_I - 1) * width + xPixel_I;
			
			red   = (float) pixelData[pixel_I][0];
			green = (float) pixelData[pixel_I][1];
			blue  = (float) pixelData[pixel_I][2];

			frame->data[0][yPixel_I * frame->linesize[0] + xPixel_I] = 
				(unsigned char)((0.257 * red) + (0.504 * green) + (0.098 * blue) + 16);
		}
	}

	/* Setup the 'Cb (U)' and 'Cr (V)' part of the frame 
	 * loop over macro pixels - pixels twice the size as the normal ones */
	for ( yPixel_I = 0 ; yPixel_I < height/2 ; yPixel_I++ ) {
		for ( xPixel_I = 0 ; xPixel_I < width/2 ; xPixel_I++ ) {
			/* Find Four pixels in this macro pixel */
			macroPixel0 = pixelData[ (height - yPixel_I*2 - 1) * width + xPixel_I * 2 ];
			macroPixel1 = pixelData[ (height - yPixel_I*2 - 1) * width + xPixel_I * 2 + 1 ];
			macroPixel2 = pixelData[ (height - yPixel_I*2 - 2) * width + xPixel_I * 2 ];
			macroPixel3 = pixelData[ (height - yPixel_I*2 - 2) * width + xPixel_I * 2 + 1 ];

			/* Average red, green and blue for four pixels around point */
			red   = 0.25 * ((float) (macroPixel0[0] + macroPixel1[0] + macroPixel2[0] + macroPixel3[0] ));
			green = 0.25 * ((float) (macroPixel0[1] + macroPixel1[1] + macroPixel2[1] + macroPixel3[1] ));
			blue  = 0.25 * ((float) (macroPixel0[2] + macroPixel1[2] + macroPixel2[2] + macroPixel3[2] ));

			/* 'Cb (U)' Component */
			frame->data[1][yPixel_I * frame->linesize[1] + xPixel_I] = 
				(unsigned char) (-(0.148 * red) - (0.291 * green) + (0.439 * blue) + 128);

			/* 'Cr (V)' Component */
			frame->data[2][yPixel_I * frame->linesize[2] + xPixel_I] = 
				(unsigned char) ((0.439 * red) - (0.368 * green) - (0.071 * blue) + 128);
		}
	}
	
	/* Write data to file */
	/* Never write frame 0  if it's a restart. If it's a normal run, write it only if specified. */
	if( context->timeStep == 0 ){
		if( self->includeFrame0 && (!context->loadFromCheckPoint) ) {
			sizeToWrite = avcodec_encode_video(self->codecContext, self->outputBuffer, MAX_BUFFER_SIZE, frame);
			Journal_Write( self->stream, self->outputBuffer, sizeToWrite, 1 );
		}
	}

	else{
		sizeToWrite = avcodec_encode_video(self->codecContext, self->outputBuffer, MAX_BUFFER_SIZE, frame);
		Journal_Write( self->stream, self->outputBuffer, sizeToWrite, 1 );
	}
}

#endif /* HAVE_LIBAVCODEC */


