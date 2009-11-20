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
** $Id: Underworld_SwarmOutput.h 354 2006-10-12 08:19:27Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __Underworld_Utils_Underworld_SwarmOutput_h__
#define __Underworld_Utils_Underworld_SwarmOutput_h__

	/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type Underworld_SwarmOutput_Type;

	typedef void (Underworld_SwarmOutput_GetFeVariableValues)(Underworld_SwarmOutput* uwSwarmOutput, FeVariable* feVariable, MaterialPointsSwarm* swarm, FILE* outputFile);
	typedef void (Underworld_SwarmOutput_PrintStandardFormat)( MaterialPoint* particle, double* result, unsigned fieldComponentCount, FILE* outputFile );

	/* Underworld_SwarmOutput information */
	#define __Underworld_SwarmOutput \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__Stg_Component \
		UnderworldContext*        context; \
		/* Virtual Info */\
		Underworld_SwarmOutput_GetFeVariableValues* _getFeValuesFunc; \
		Underworld_SwarmOutput_PrintStandardFormat* _printFunc; \
		MaterialPointsSwarm*      materialSwarm; \
		FeVariable**              feVariableList; \
		Index                     sizeList; \
		Stream*                   infoStream; \
		Stream*                   errorStream; \

	struct Underworld_SwarmOutput { __Underworld_SwarmOutput };
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/
	Underworld_SwarmOutput* _Underworld_SwarmOutput_New(
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
		Name                                               name );

	/* Stg_Class_Delete Underworld_SwarmOutput implementation */
	void _Underworld_SwarmOutput_Delete( void* swarmOutput );
	void _Underworld_SwarmOutput_Print( void* swarmOutput, Stream* stream );
	#define Underworld_SwarmOutput_Copy( self ) \
		(Underworld_SwarmOutput*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define Underworld_SwarmOutput_DeepCopy( self ) \
		(Underworld_SwarmOutput*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _Underworld_SwarmOutput_Copy( void* swarmOutput, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	void* _Underworld_SwarmOutput_DefaultNew( Name name ) ;
   void _Underworld_SwarmOutput_AssignFromXML( void* shape, Stg_ComponentFactory* cf, void* data ) ;
   void _Underworld_SwarmOutput_Init( Underworld_SwarmOutput* self, PICelleratorContext* context, MaterialPointsSwarm* materialSwarm, unsigned int listCount, FeVariable** feVariableList );
   
	void _Underworld_SwarmOutput_Build( void* swarmOutput, void* data ) ;
	void _Underworld_SwarmOutput_Initialise( void* swarmOutput, void* data ) ;
	void _Underworld_SwarmOutput_Execute( void* swarmOutput, void* data );
	void _Underworld_SwarmOutput_Destroy( void* swarmOutput, void* data ) ;
	
	void _Underworld_SwarmOutput_PrintHeader( void* swarmOutput, Stream* stream, Particle_Index lParticle_I, void* context );
	void _Underworld_SwarmOutput_PrintData( void* swarmOutput, Stream* stream, Particle_Index lParticle_I, void* context );
		
	/*---------------------------------------------------------------------------------------------------------------------
	** Private functions
	*/
	
	void _Underworld_SwarmOutput_GetFeVariableValues(Underworld_SwarmOutput* uwSwarmOutput, FeVariable* feVariable, MaterialPointsSwarm* swarm, FILE* outputFile );
	void _Underworld_SwarmOutput_PrintStandardFormat( MaterialPoint* particle, double* result, unsigned fieldComponentCount, FILE* outputFile  );
	/*---------------------------------------------------------------------------------------------------------------------
	** Entry Point Hooks
	*/
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/
	
#endif 
