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
** $Id: OutputFormat.h 628 2006-10-12 08:23:07Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __lucOutputFormat_h__
#define __lucOutputFormat_h__

	extern const Type lucOutputFormat_Type;

	typedef void (lucOutputFormat_OutputFunction) ( void* outputFormat, lucWindow* window, AbstractContext* context, lucPixel* pixelData );

	#define __lucOutputFormat                                         \
		__Stg_Component                                           \
		AbstractContext*				   context;		     \
		/* Virtual Functions */ \
		lucOutputFormat_OutputFunction*                    _output;                  \
		/* Other Info */   \
		Name                                               extension;

	struct lucOutputFormat {__lucOutputFormat};

	lucOutputFormat* _lucOutputFormat_New(
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
		lucOutputFormat_OutputFunction*                    _output,
		Name                                               name );

	void lucOutputFormat_InitAll( 
		void*                                              outputFormat,
		Name                                               extension );

	void _lucOutputFormat_Delete( void* outputFormat ) ;
	void _lucOutputFormat_Print( void* outputFormat, Stream* stream ) ;
	void* _lucOutputFormat_Copy( void* outputFormat, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) ;

	void* _lucOutputFormat_DefaultNew( Name name ) ;
void _lucOutputFormat_AssignFromXML( void* outputFormat, Stg_ComponentFactory* cf, void* data ) ;
	void _lucOutputFormat_Build( void* outputFormat, void* data );
	void _lucOutputFormat_Initialise( void* outputFormat, void* data );
	void _lucOutputFormat_Execute( void* outputFormat, void* data );
	void _lucOutputFormat_Destroy( void* outputFormat, void* data );

	/* +++ Public Functions +++ */
	void lucOutputFormat_Output( void* outputFormat, lucWindow* window, AbstractContext* context, lucPixel* pixelData ) ;

	Name lucOutputFormat_GetImageFilename( void* outputFormat, lucWindow* window, void* _context ) ;
	FILE* lucOutputFormat_OpenFile( void* outputFormat, lucWindow* window, void* _context, const char *mode ) ;
	Stream* lucOutputFormat_OpenStream( void* outputFormat, lucWindow* window, void* context ) ;

#endif
