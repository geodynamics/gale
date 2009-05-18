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
** $Id: Init.c 740 2007-10-11 08:05:31Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/



#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <glucifer/Base/Base.h>

#include "Windowing.h"

const Type lucDefaultWindow_Type = "lucDefaultWindow";

Dictionary* lucDefaultWindow_MetaAsDictionary() {
	return Dictionary_New();
}

Dictionary* lucDefaultWindow_Type_MetaAsDictionary() {
	return lucDefaultWindow_MetaAsDictionary();
}

Bool lucWindowing_Init() {
	Stg_ComponentRegister* componentRegister = Stg_ComponentRegister_Get_ComponentRegister();

	Journal_Printf( Journal_Register( DebugStream_Type, "Context" ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */

    /* Order of priority for default output window: SDL, OSMesa, Carbon, X11, VTK */
    /* SDL will work with OSMesa to allow on/off-screen rendering so if present is first choice */
    /* If OSMesa is linked, Carbon and X11 will NOT work as they depend on the system OpenGL */
    /* library, which OSMesa replaces. */ 

	#ifdef HAVE_SDL
		Stg_ComponentRegister_Add( componentRegister, lucSDLWindow_Type,     "0", _lucSDLWindow_DefaultNew );
		RegisterParent( lucSDLWindow_Type, lucWindow_Type );
		if ( !Stg_ComponentRegister_Get( componentRegister, lucDefaultWindow_Type, "0" ) )
			Stg_ComponentRegister_Add( componentRegister, lucDefaultWindow_Type, "0", _lucSDLWindow_DefaultNew );		
	#endif		
	
	#ifdef HAVE_OSMESA  
		Stg_ComponentRegister_Add( componentRegister, lucOSMesaWindow_Type,     "0", _lucOSMesaWindow_DefaultNew );
		RegisterParent( lucOSMesaWindow_Type, lucWindow_Type );
		if ( !Stg_ComponentRegister_Get( componentRegister, lucDefaultWindow_Type, "0" ) )
			Stg_ComponentRegister_Add( componentRegister, lucDefaultWindow_Type, "0", _lucOSMesaWindow_DefaultNew );
	#endif
	
	#ifdef HAVE_CARBON
		Stg_ComponentRegister_Add( componentRegister, lucCarbonWindow_Type,     "0", _lucCarbonWindow_DefaultNew );
		RegisterParent( lucCarbonWindow_Type, lucWindow_Type );
		if ( !Stg_ComponentRegister_Get( componentRegister, lucDefaultWindow_Type, "0" ) )
			Stg_ComponentRegister_Add( componentRegister, lucDefaultWindow_Type, "0", _lucCarbonWindow_DefaultNew );

	#endif
		
	#ifdef HAVE_X11
		Stg_ComponentRegister_Add( componentRegister, lucX11Window_Type,     "0", _lucX11Window_DefaultNew );
		RegisterParent( lucX11Window_Type, lucWindow_Type );
		if ( !Stg_ComponentRegister_Get( componentRegister, lucDefaultWindow_Type, "0" ) )
			Stg_ComponentRegister_Add( componentRegister, lucDefaultWindow_Type, "0", _lucX11Window_DefaultNew );
	#endif	

	#ifdef HAVE_VTK
		Stg_ComponentRegister_Add( componentRegister, lucVTKWindow_Type,     "0", _lucVTKWindow_DefaultNew );
		RegisterParent( lucVTKWindow_Type, lucWindow_Type );
		if ( !Stg_ComponentRegister_Get( componentRegister, lucDefaultWindow_Type, "0" ) )
			Stg_ComponentRegister_Add( componentRegister, lucDefaultWindow_Type, "0", _lucVTKWindow_DefaultNew );		
	#endif	
	
	return True;
}
