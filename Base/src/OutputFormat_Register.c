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
** $Id: OutputFormat_Register.c 740 2007-10-11 08:05:31Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "types.h"
#include "OutputFormat_Register.h"
#include "OutputFormat.h"

const Type lucOutputFormat_Register_Type = "lucOutputFormat_Register";

lucOutputFormat_Register*	lucOutputFormat_Register_New( void ) {
	/* Variables set in this function */
	SizeT                      _sizeOfSelf = sizeof(lucOutputFormat_Register);
	Type                              type = lucOutputFormat_Register_Type;
	Stg_Class_DeleteFunction*      _delete = _NamedObject_Register_Delete;
	Stg_Class_PrintFunction*        _print = _NamedObject_Register_Print;
	Stg_Class_CopyFunction*          _copy = _NamedObject_Register_Copy;

	lucOutputFormat_Register* self;
	
	self = (lucOutputFormat_Register*) _NamedObject_Register_New(  NAMEDOBJECT_REGISTER_PASSARGS  );

	return self;
}

void lucOutputFormat_Register_OutputAll( void* outputFormat_Register, lucWindow* window, AbstractContext* context, lucPixel* pixelData, lucAlphaPixel* alphaPixelData ) {
	lucOutputFormat_Register* self          = (lucOutputFormat_Register*) outputFormat_Register;
	OutputFormat_Index        object_I;
	OutputFormat_Index        objectCount   = lucOutputFormat_Register_GetCount( self );
	lucOutputFormat*          object;

	for ( object_I = 0 ; object_I < objectCount ; object_I++ ) {
		object = lucOutputFormat_Register_GetByIndex( self, object_I );
      if (object->transparent)
   		lucOutputFormat_Output( object, window, context, alphaPixelData );
      else
   		lucOutputFormat_Output( object, window, context, pixelData );
	}
}




