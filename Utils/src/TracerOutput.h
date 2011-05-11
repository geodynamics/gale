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
*%		Louis Moresi - Louis.Moresi@sci.monash.edu.au
*%
** Contributors:
*+		Robert Turnbull
*+		Vincent Lemiale
*+		Louis Moresi
*+		David May
*+		David Stegman
*+		Mirko Velic
*+		Patrick Sunter
*+		Julian Giordani
*+
** $Id: TracerOutput.h 354 2006-10-12 08:19:27Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __Underworld_Utils_TracerOutput_h__
#define __Underworld_Utils_TracerOutput_h__

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
extern const Type TracerOutput_Type;

/* TracerOutput information */
#define __TracerOutput                                                  \
  /* Macro defining parent goes here - This means you can cast this class as its parent */ \
  __SwarmOutput                                                         \
  /* Virtual Info */                                                    \
  FeVariable* pressureField;      \
  FeVariable** fields;   \
  unsigned int num_fields; \

struct TracerOutput { __TracerOutput };
	
/*---------------------------------------------------------------------------------------------------------------------
** Constructors
*/
	
#ifndef ZERO
#define ZERO 0
#endif

#define TRACEROUTPUT_DEFARGS       \
  SWARMOUTPUT_DEFARGS

#define TRACEROUTPUT_PASSARGS      \
  SWARMOUTPUT_PASSARGS

TracerOutput* _TracerOutput_New(  TRACEROUTPUT_DEFARGS  );

/* Stg_Class_Delete TracerOutput implementation */
void _TracerOutput_Delete( void* swarmOutput );
void _TracerOutput_Print( void* swarmOutput, Stream* stream );
#define TracerOutput_Copy( self )                                       \
  (TracerOutput*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
#define TracerOutput_DeepCopy( self )                                   \
  (TracerOutput*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
void* _TracerOutput_Copy( const void* swarmOutput, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
void* _TracerOutput_DefaultNew( Name name ) ;
void _TracerOutput_AssignFromXML( void* shape, Stg_ComponentFactory* cf, void* data ) ;
void _TracerOutput_Build( void* swarmOutput, void* data ) ;
void _TracerOutput_Initialise( void* swarmOutput, void* data ) ;
void _TracerOutput_Execute( void* swarmOutput, void* data );
void _TracerOutput_Destroy( void* swarmOutput, void* data ) ;
	
void _TracerOutput_PrintHeader( void* swarmOutput, Stream* stream, Particle_Index lParticle_I, void* context );
void _TracerOutput_PrintData( void* swarmOutput, Stream* stream, Particle_Index lParticle_I, void* context );
		
/*---------------------------------------------------------------------------------------------------------------------
** Private functions
*/
	
/*---------------------------------------------------------------------------------------------------------------------
** Entry Point Hooks
*/
	
/*---------------------------------------------------------------------------------------------------------------------
** Public functions
*/
	
#endif 

