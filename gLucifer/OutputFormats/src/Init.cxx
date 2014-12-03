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
** $Id: Init.c 768 2008-04-21 03:20:07Z JohnMansour $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <glucifer/Base/Base.h>

#include "OutputFormats.h"
#ifdef HAVE_LIBAVCODEC
	#include "ffmpeg/avcodec.h"
#endif

Bool lucOutputFormats_Init() {
	Stg_ComponentRegister* componentRegister = Stg_ComponentRegister_Get_ComponentRegister();

	Journal_Printf( Journal_Register( DebugStream_Type, (Name)"Context"  ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */

	Stg_ComponentRegister_Add( componentRegister, lucOutputPPM_Type, (Name)"0", _lucOutputPPM_DefaultNew  );
	RegisterParent( lucOutputPPM_Type,         lucOutputFormat_Type );

	#ifdef HAVE_GL2PS
	  Stg_ComponentRegister_Add( componentRegister, lucOutputVECTOR_Type, (Name)"0", _lucOutputVECTOR_DefaultNew  );
	  RegisterParent( lucOutputVECTOR_Type,         lucOutputFormat_Type );
	#endif

	#ifdef HAVE_LIBPNG
	   Stg_ComponentRegister_Add( componentRegister, lucOutputPNG_Type, (Name)"0", _lucOutputPNG_DefaultNew  );
	   RegisterParent( lucOutputPNG_Type,         lucOutputFormat_Type );
	#endif	
	
	#ifdef HAVE_LIBJPEG
    	Stg_ComponentRegister_Add( componentRegister, lucOutputJPEG_Type, (Name)"0", _lucOutputJPEG_DefaultNew  );
    	RegisterParent( lucOutputJPEG_Type,        lucOutputFormat_Type );
    #endif
    
    #ifdef HAVE_TIFF
    	Stg_ComponentRegister_Add( componentRegister, lucOutputTIFF_Type, (Name)"0", _lucOutputTIFF_DefaultNew  );
    	RegisterParent( lucOutputTIFF_Type,        lucOutputFormat_Type );
    #endif
    
    #ifdef HAVE_LIBFAME	
    	Stg_ComponentRegister_Add( componentRegister, lucEncoderLibfame_Type, (Name)"0", _lucEncoderLibfame_DefaultNew  );
    	RegisterParent( lucEncoderLibfame_Type,    lucOutputFormat_Type );
    #endif
    
    #ifdef HAVE_LIBAVCODEC	
    	Stg_ComponentRegister_Add( componentRegister, lucEncoderLibavcodec_Type, (Name)"0", _lucEncoderLibavcodec_DefaultNew  );
    	RegisterParent( lucEncoderLibavcodec_Type,    lucOutputFormat_Type );
		avcodec_init();
		avcodec_register_all();
    #endif

	return True;
}


