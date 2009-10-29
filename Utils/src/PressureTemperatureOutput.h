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
** $Id: PressureTemperatureOutput.h 354 2006-10-12 08:19:27Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __Underworld_Utils_PressureTemperatureOutput_h__
#define __Underworld_Utils_PressureTemperatureOutput_h__

	/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type PressureTemperatureOutput_Type;

	/* PressureTemperatureOutput information */
	#define __PressureTemperatureOutput \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__SwarmOutput \
		/* Virtual Info */\
		FeVariable*                                       pressureField;                \
		FeVariable*                                       temperatureField;             \

	struct PressureTemperatureOutput { __PressureTemperatureOutput };
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/
	PressureTemperatureOutput* _PressureTemperatureOutput_New(
		SizeT                                              _sizeOfSelf, 
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
		SwarmOutput_PrintHeaderFunction*                   _printHeader,		
		SwarmOutput_PrintDataFunction*                     _printData,		
		Name                                               name );

	/* Stg_Class_Delete PressureTemperatureOutput implementation */
	void _PressureTemperatureOutput_Delete( void* swarmOutput );
	void _PressureTemperatureOutput_Print( void* swarmOutput, Stream* stream );
	#define PressureTemperatureOutput_Copy( self ) \
		(PressureTemperatureOutput*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define PressureTemperatureOutput_DeepCopy( self ) \
		(PressureTemperatureOutput*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _PressureTemperatureOutput_Copy( void* swarmOutput, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	void* _PressureTemperatureOutput_DefaultNew( Name name ) ;
void _PressureTemperatureOutput_AssignFromXML( void* shape, Stg_ComponentFactory* cf, void* data ) ;
	void _PressureTemperatureOutput_Build( void* swarmOutput, void* data ) ;
	void _PressureTemperatureOutput_Initialise( void* swarmOutput, void* data ) ;
	void _PressureTemperatureOutput_Execute( void* swarmOutput, void* data );
	void _PressureTemperatureOutput_Destroy( void* swarmOutput, void* data ) ;
	
	void _PressureTemperatureOutput_PrintHeader( void* swarmOutput, Stream* stream, Particle_Index lParticle_I, void* context );
	void _PressureTemperatureOutput_PrintData( void* swarmOutput, Stream* stream, Particle_Index lParticle_I, void* context );
		
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
