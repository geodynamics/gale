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
** $Id: ColourMap.h 784 2008-08-18 13:54:31Z LukeHodkinson $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/



#ifndef __lucColour_ColourMap_h__
#define __lucColour_ColourMap_h__

	extern const Type lucColourMap_Type;

	struct lucColour {
		float red;
		float green;
		float blue;
		float opacity;
	};

	#define __lucColourMap                \
		__Stg_Component			\
		AbstractContext*				   context;      \
		Colour_Index                                       colourCount;  \
		lucColour*                                         colourList;   \
		double                                             minimum;      \
		double                                             maximum;      \
		Bool                                               logScale;     \
		Bool                                               dynamicRange; \
        Bool 											   discrete; \
		Bool											   centreOnFixedValue; \
		double											   centringValue;

	struct lucColourMap {__lucColourMap};

	/** Constructors */
	lucColourMap* lucColourMap_New( 
		Name                                               name,
		char*                                              _colourMapString, 
		double                                             minimum,
		double                                             maximum,
		Bool                                               logScale,
		Bool                                               dynamicRange,
		Bool											   centreOnFixedValue,
		double											   centringValue
		 );

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define LUCCOLOURMAP_DEFARGS \
                STG_COMPONENT_DEFARGS

	#define LUCCOLOURMAP_PASSARGS \
                STG_COMPONENT_PASSARGS

	lucColourMap* _lucColourMap_New(  LUCCOLOURMAP_DEFARGS  );

	void lucColourMap_InitAll( 
		void*                                              colourMap, 
		char*                                              _colourMapString, 
		double                                             minimum,
		double                                             maximum,
		Bool                                               logScale,
		Bool                                               dynamicRange,
		Bool											   centreOnFixedValue,
		double											   centringValue
		 );

	/** Virtual Functions */
	void _lucColourMap_Delete( void* colourMap ) ;
	void _lucColourMap_Print( void* colourMap, Stream* stream ) ;
	void* _lucColourMap_Copy( void* colourMap, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) ;

	void* _lucColourMap_DefaultNew( Name name ) ;
    void _lucColourMap_AssignFromXML( void* colourMap, Stg_ComponentFactory* cf, void* data ) ;
	void _lucColourMap_Build( void* colourMap, void* data ) ;
	void _lucColourMap_Initialise( void* colourMap, void* data ) ;
	void _lucColourMap_Execute( void* colourMap, void* data ) ;
	void _lucColourMap_Destroy( void* colourMap, void* data ) ;

	/** Public Functions */
	void lucColourMap_GetColourFromValue( void* colourMap, double value, lucColour* colour ) ;
	float lucColourMap_ScaleValue( void* colourMap, double value ) ;
	void lucColourMap_GetColourFromScaledValue( void* colourMap, float scaledValue, lucColour* colour ) ;

	#define lucColourMap_GetColourFromList( self, colour_I ) \
		(&(self)->colourList[ colour_I ])

	void lucColourMap_SetMinMax( void* colourMap, double min, double max ) ;
	void lucColourMap_CalibrateFromVariable( void* colourMap, void* _variable ) ;
	void lucColourMap_CalibrateFromFieldVariable( void* colourMap, void* _fieldVariable ) ;
	void lucColourMap_CalibrateFromSwarmVariable( void* colourMap, void* swarmVariable ) ;

	void lucColour_FromHSV( lucColour* self, float hue, float saturation, float value, float opacity ) ;
	void lucColour_FromString( lucColour* self, char* string ) ;

	void lucColour_SetTerminalColours( lucColour* textColour, lucColour* backgroundColour, Stream* stream ) ;

#endif

