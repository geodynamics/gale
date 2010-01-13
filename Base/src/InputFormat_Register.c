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

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "types.h"
#include "InputFormat_Register.h"
#include "Init.h"

#include <string.h>
#include <assert.h>

const Type lucInputFormat_Register_Type = "lucInputFormat_Register";

lucInputFormat_Register* lucInputFormat_Register_New() {
	/* Variables set in this function */
	SizeT                      _sizeOfSelf = sizeof( lucInputFormat_Register );
	Type                              type = lucInputFormat_Register_Type;
	Stg_Class_DeleteFunction*      _delete = _lucInputFormat_Register_Delete;
	Stg_Class_PrintFunction*        _print = _lucInputFormat_Register_Print;
	Stg_Class_CopyFunction*          _copy = _lucInputFormat_Register_Copy;

	lucInputFormat_Register* self;

	self = _lucInputFormat_Register_New(  LUCINPUTFORMAT_REGISTER_PASSARGS  );

	lucInputFormat_Register_InitAll( self );

	return self;
}

lucInputFormat_Register* _lucInputFormat_Register_New(  LUCINPUTFORMAT_REGISTER_DEFARGS  )
{
	lucInputFormat_Register*    self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucInputFormat_Register) );
	self = (lucInputFormat_Register*) _Stg_ComponentRegister_New(  STG_COMPONENTREGISTER_PASSARGS  );
	
	return self;
}

void _lucInputFormat_Register_Init( void* inputFormat_Register ) {
}

void lucInputFormat_Register_InitAll( 
		void*                                              inputFormat_Register ) 
{
	lucInputFormat_Register* self        = inputFormat_Register;

	Stg_ComponentRegister_Init( ( Stg_ComponentRegister*)self );	
	_lucInputFormat_Register_Init( self );
}

void _lucInputFormat_Register_Delete( void* inputFormat_Register ) {
	lucInputFormat_Register* self        = inputFormat_Register;
	
	_Stg_ComponentRegister_Delete( self );
}

void _lucInputFormat_Register_Print( void* inputFormat_Register, Stream* stream ) {
	lucInputFormat_Register* self        = (lucInputFormat_Register*) inputFormat_Register;
	
	_Stg_ComponentRegister_Print( self, stream );
}

void* _lucInputFormat_Register_Copy( void* inputFormat_Register, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	lucInputFormat_Register* self        = inputFormat_Register;
	lucInputFormat_Register* newInputFormat_Register;

	newInputFormat_Register = (lucInputFormat_Register*) _Stg_Component_Copy( self, dest, deep, nameExt, ptrMap );

	return (void*) newInputFormat_Register;
}

/* This goes through each input format and finds the which one works with the extension for the filename passed in */
lucInputFormat* lucInputFormat_Register_CreateFromFileName( void* inputFormat_Register, Name imageName ) {
	lucInputFormat_Register*                  self        = (lucInputFormat_Register*) inputFormat_Register;
	Stg_Component_DefaultConstructorFunction* defaultNewFunctionPtr;
	Name                                      extension;
	Stream* errorStream = Journal_Register( Error_Type, (Name)lucInputFormat_Register_Type  );

	Journal_Firewall( imageName != NULL, errorStream, "In func '%s for %s Image file name %s cannot be found \n", 
			__func__, self->type, imageName );

	/* Find extension of image name */ 
	extension = strrchr( imageName, '.' );
	Journal_Firewall( extension != (char*) 1, errorStream, 
			"In func %s for %s - Cannot find extension for filename '%s'.\n", 
			__func__, self->type, imageName );

	defaultNewFunctionPtr = Stg_ComponentRegister_AssertGet( (Stg_ComponentRegister*)self, extension, "0" );

	return (lucInputFormat*) defaultNewFunctionPtr( imageName );
}
	


