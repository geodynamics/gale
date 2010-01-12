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
** $Id: Light.h 510 2006-02-17 04:33:32Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __lucLight_h__
#define __lucLight_h__

	extern const Type lucLight_Type;

	extern const double LUC_LIGHT_DEFAULT_POS_X;
	extern const double LUC_LIGHT_DEFAULT_POS_Y;
	extern const double LUC_LIGHT_DEFAULT_POS_Z;
	extern const double LUC_LIGHT_DEFAULT_POS_W;

	#define __lucLight                                \
		__Stg_Component                                \
		AbstractContext*		            context; \
		Light_Index 		               index;\
		int                              model; \
		int                              material;\
		float	                           position[4];\
		float                            lmodel_ambient[4];\
		float                            spotCutOff;\
		float                            spotDirection[3];\
		Bool                             needsToDraw;          \

		
	struct lucLight {__lucLight};

	/** Constructors */
	lucLight* lucLight_New( 
		Name                                name,
		Light_Index            				   index,
		int                                 model,
		int                                 material,
		float						               position[4],
		float                               lmodel_ambient[4],
		float                               spotCutOff,
		float                               spotDirection[3]
	 );

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define LUCLIGHT_DEFARGS \
                STG_COMPONENT_DEFARGS

	#define LUCLIGHT_PASSARGS \
                STG_COMPONENT_PASSARGS

	lucLight* _lucLight_New(  LUCLIGHT_DEFARGS  );

	/** Virtual Functions */
	void _lucLight_Delete( void* light ) ;
	void _lucLight_Print( void* light, Stream* stream ) ;
	void* _lucLight_Copy( void* light, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) ;
	#define lucLight_Copy( self ) \
		(lucLight*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	void* _lucLight_DefaultNew( Name name ) ;
	void _lucLight_AssignFromXML( void* light, Stg_ComponentFactory* cf, void* data ) ;
	void _lucLight_Build( void* light, void* data );
	void _lucLight_Initialise( void* light, void* data );
	void _lucLight_Execute( void* light, void* data );
	void _lucLight_Destroy( void* light, void* data );

	/** Public Functions */
	void lucLight_Position( void * light, int index, float posX, float posY, float posZ, float posW);
	void lucLight_Material( int material);



#endif

