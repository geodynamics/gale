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



#ifndef __lucInputFormat_Register_h__
#define __lucInputFormat_Register_h__

	extern const Type lucInputFormat_Register_Type;

	#define __lucInputFormat_Register                                  \
		__Stg_ComponentRegister                                 \

	struct lucInputFormat_Register {__lucInputFormat_Register};

	/** Constructors */
	lucInputFormat_Register* lucInputFormat_Register_New(  );

	lucInputFormat_Register* _lucInputFormat_Register_New(
		SizeT                                              sizeOfSelf,
		Type                                               type,
		Stg_Class_DeleteFunction*                          _delete,
		Stg_Class_PrintFunction*                           _print,
		Stg_Class_CopyFunction*                            _copy );

	void lucInputFormat_Register_InitAll(  ) ;
	
	/** Virtual Functions */
	void _lucInputFormat_Register_Delete( void* inputFormat_Register ) ;
	void _lucInputFormat_Register_Print( void* inputFormat_Register, Stream* stream ) ;
	void* _lucInputFormat_Register_Copy( void* inputFormat_Register, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) ;
	#define lucInputFormat_Register_Copy( self ) \
		(lucInputFormat_Register*) Stg_Class_Copy( self, NULL, False, NULL, NULL )

	#define lucInputFormat_Register_Add( self, extension, componentType, version, func ) \
		Stg_ComponentRegister_AddFunc( (Stg_ComponentRegister*)(self), (extension), (version), (func), NULL ); \
		Stg_ComponentRegister_Add( (Stg_ComponentRegister*)(self), componentType, (version), (func) );

	/* This goes through each input format and finds the which one works with the extension for the filename passed in */
	lucInputFormat* lucInputFormat_Register_CreateFromFileName( void* inputFormat_Register, Name filename );

#endif
