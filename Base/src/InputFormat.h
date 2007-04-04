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


#ifndef __lucInputFormat_h__
#define __lucInputFormat_h__

	extern const Type lucInputFormat_Type;

	typedef lucPixel* (lucInputFormat_InputFunction) (void* inputFormat,  Name imageName, Pixel_Index* width, Pixel_Index* height);

	#define __lucInputFormat                                         \
		__Stg_Component                                           \
		/* Virtual Functions */ \
		lucInputFormat_InputFunction*                    _input;   \
		/* Other Info */   \
		Name                                             imageName;\
		
		

	struct lucInputFormat {__lucInputFormat};

	lucInputFormat* _lucInputFormat_New(
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
		lucInputFormat_InputFunction*                      _input,
		Name                                               name );

	void lucInputFormat_InitAll( 
		void*                                              inputFormat,
		Name                                               extension );

	void _lucInputFormat_Delete( void* inputFormat ) ;
	void _lucInputFormat_Print( void* inputFormat, Stream* stream ) ;
	void* _lucInputFormat_Copy( void* inputFormat, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) ;

	void* _lucInputFormat_DefaultNew( Name name ) ;
void _lucInputFormat_Construct( void* inputFormat, Stg_ComponentFactory* cf, void* data ) ;
	void _lucInputFormat_Build( void* inputFormat, void* data );
	void _lucInputFormat_Initialise( void* inputFormat, void* data );
	void _lucInputFormat_Execute( void* inputFormat, void* data );
	void _lucInputFormat_Destroy( void* inputFormat, void* data );

	/* +++ Public Functions +++ */
	lucPixel* lucInputFormat_Input(  void* inputFormat, Name imageName, Pixel_Index* imageWidth, Pixel_Index* imageHeight ) ;

	FILE* lucInputFormat_OpenFile( void* inputFormat, lucWindow* window, void* _context, const char *mode ) ;

#endif
