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
** $Id: OutputJPEG.h 628 2006-10-12 08:23:07Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __lucOutputJPEG_h__
#define __lucOutputJPEG_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type lucOutputJPEG_Type;
		
	/** Class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __lucOutputJPEG \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__lucOutputFormat \
		/* Virtual functions go here */ \
		/* Other info */\
		int                              quality;                   \


	struct lucOutputJPEG { __lucOutputJPEG };
	
	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define LUCOUTPUTJPEG_DEFARGS \
                LUCOUTPUTFORMAT_DEFARGS

	#define LUCOUTPUTJPEG_PASSARGS \
                LUCOUTPUTFORMAT_PASSARGS

	lucOutputJPEG* _lucOutputJPEG_New(  LUCOUTPUTJPEG_DEFARGS  );

	void _lucOutputJPEG_Delete( void* outputFormat ) ;
	void _lucOutputJPEG_Print( void* outputFormat, Stream* stream ) ;
	void* _lucOutputJPEG_Copy( void* outputFormat, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) ;

	/* 'Stg_Component' implementations */
	void* _lucOutputJPEG_DefaultNew( Name name ) ;
	void _lucOutputJPEG_AssignFromXML( void* outputFormat, Stg_ComponentFactory* cf, void* data );
	void _lucOutputJPEG_Build( void* outputFormat, void* data ) ;
	void _lucOutputJPEG_Initialise( void* outputFormat, void* data ) ;
	void _lucOutputJPEG_Execute( void* outputFormat, void* data );
	void _lucOutputJPEG_Destroy( void* outputFormat, void* data ) ;
	
	void _lucOutputJPEG_Output( void* outputFormat, lucWindow* window, AbstractContext* context, void* pixelData ) ;

#endif

