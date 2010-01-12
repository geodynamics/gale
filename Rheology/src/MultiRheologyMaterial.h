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
** $Id: MultiRheologyMaterial.h 446 2007-03-04 09:55:46Z PatrickSunter $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __Underworld_Rheology_MultiRheologyMaterial_h__
#define __Underworld_Rheology_MultiRheologyMaterial_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type MultiRheologyMaterial_Type;
		
	/** Class contents */
	#define __MultiRheologyMaterial \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__RheologyMaterial \
		/* Virtual functions go here */ \
		/* Material Parameters */\
		Rheology_Register**	rheology_RegisterList; \
		Index						rheology_RegisterCount;

	struct MultiRheologyMaterial { __MultiRheologyMaterial };
	
	/** Public "New" C++-Style constructor */
	MultiRheologyMaterial* MultiRheologyMaterial_New( 
		Name						name,
		PICelleratorContext*	context,
		Stg_Shape*				shape,
		Dictionary*				materialDictionary,
		Materials_Register*	materialRegister,
		Rheology**				rheologyList,
		Rheology_Index			rheologyCount,
		Compressible*			compressible,
		Rheology***				rheologyListList,
		Rheology_Index*		rheologyCountList, 
		Index						rheologyListCount );

	void* _MultiRheologyMaterial_DefaultNew( Name name ) ;

	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define MULTIRHEOLOGYMATERIAL_DEFARGS \
                RHEOLOGYMATERIAL_DEFARGS

	#define MULTIRHEOLOGYMATERIAL_PASSARGS \
                RHEOLOGYMATERIAL_PASSARGS

	MultiRheologyMaterial* _MultiRheologyMaterial_New(  MULTIRHEOLOGYMATERIAL_DEFARGS  );

	void _MultiRheologyMaterial_AssignFromXML( void* material, Stg_ComponentFactory* cf, void* data );

	void _MultiRheologyMaterial_Init( 
		MultiRheologyMaterial*  self, 
		Rheology***             rheologyListList,
		Rheology_Index*         rheologyCountList, 
		Index                   rheologyListCount );

	/* 'Stg_Component' implementations */
	void _MultiRheologyMaterial_Delete( void* material );
   void _MultiRheologyMaterial_Destroy( void* material, void* data );
   void _MultiRheologyMaterial_Init( 
		MultiRheologyMaterial*  self, 
		Rheology***             rheologyListList,
		Rheology_Index*         rheologyCountList, 
		Index                   rheologyListCount );

	void _MultiRheologyMaterial_Destroy( void* material, void* data );

	void _MultiRheologyMaterial_RunRheologies( 	
		void*						material,
		ConstitutiveMatrix*	constitutiveMatrix,
		MaterialPointsSwarm*	swarm,
		Element_LocalIndex	lElement_I,
		MaterialPoint*			materialPoint,
		Coord						xi );

#endif

