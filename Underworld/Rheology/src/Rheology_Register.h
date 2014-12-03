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
** $Id: Rheology_Register.h 169 2006-04-21 07:08:28Z AlanLo $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __Underworld_Rheology_Rheology_Register_h__
#define __Underworld_Rheology_Rheology_Register_h__
	
	extern const Type Rheology_Register_Type;
	
	#define __Rheology_Register \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__NamedObject_Register \
		\
		/* Virtual functions go here */ \
		\
		/* Class info */ \

	struct Rheology_Register { __Rheology_Register };
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructor
	*/
	
	Rheology_Register*	Rheology_Register_New( void );
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** General virtual functions
	*/
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/
	#define Rheology_Register_Add NamedObject_Register_Add

	#define Rheology_Register_GetIndex NamedObject_Register_GetIndex 

	#define Rheology_Register_GetByName( self, rheologyName ) \
		( (Rheology*) NamedObject_Register_GetByName( self, rheologyName ) ) 

	#define Rheology_Register_GetByIndex( self, rheology_I ) \
		( (Rheology*) NamedObject_Register_GetByIndex( self, rheology_I ) )

	#define Rheology_Register_GetCount( self ) \
		(self)->objects->count
	
	#define Rheology_Register_DeleteAllObjects( self ) \
		NamedObjectList_DeleteAllObjects( (self)->objects )
	#define Rheology_Register_PrintAllObjects( self, stream ) \
		NamedObjectList_PrintAllObjects( (self)->objects, stream )

	void Rheology_Register_RunRheologies( 	
		Rheology_Register*                                 self,
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               swarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint, 
		Coord                                              xi );
#endif 
