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
** $Id: RheologyMaterial.h 446 2007-03-04 09:55:46Z PatrickSunter $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __Underworld_Rheology_RheologyMaterial_h__
#define __Underworld_Rheology_RheologyMaterial_h__

	/* typedefs for virtual functions: */
	typedef void  (RheologyMaterial_RunRheologiesFunction)            (
		void*                                              rheologyMaterial,
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               swarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi );

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type RheologyMaterial_Type;
	
	/** Rheology class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __RheologyMaterial \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__Material \
		/* Virtual functions go here */ \
		RheologyMaterial_RunRheologiesFunction*            _runRheologies;    \
		/* Other info */ \
		Rheology_Register*                                 rheology_Register; \
		Compressible*                                      compressible;

	struct RheologyMaterial { __RheologyMaterial };

	RheologyMaterial* RheologyMaterial_New( 
		Name                             name,
		Stg_Shape*                       shape,
		Dictionary*                      materialDictionary,
		Materials_Register*              materialRegister,
		Rheology**                       rheologyList,
		Rheology_Index                   rheologyCount,
		Compressible*                    compressible );

	void* _RheologyMaterial_DefaultNew( Name name ) ;

	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	RheologyMaterial* _RheologyMaterial_New( 
		SizeT                                              sizeOfSelf,
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
		RheologyMaterial_RunRheologiesFunction*            _runRheologies,
		Name                                               name,
		Stg_Shape*                                         shape,
		Dictionary*                                        materialDictionary,
		Materials_Register*                                materialRegister,
		Rheology**                                         rheologyList,
		Rheology_Index                                     rheologyCount,
		Compressible*                                      compressible );

	void _RheologyMaterial_Construct( void* rheologyMaterial, Stg_ComponentFactory* cf, void* data );

	void _RheologyMaterial_Init(
		void*                                              rheologyMaterial,
		Rheology**                                         rheologyList,
		Rheology_Index                                     rheologyCount,
		Compressible*                                      compressible );

	/* 'Stg_Class' implementations */
	void _RheologyMaterial_Delete( void* rheologyMaterial );
	void _RheologyMaterial_Print( void* rheologyMaterial, Stream* stream );
	#define RheologyMaterial_Copy( self ) \
		(RheologyMaterial*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define RheologyMaterial_DeepCopy( self ) \
		(RheologyMaterial*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _RheologyMaterial_Copy( void* rheologyMaterial, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	/* 'Stg_Component' implementations */
	void _RheologyMaterial_Build( void* rheologyMaterial, void* data );
	void _RheologyMaterial_Initialise( void* rheologyMaterial, void* data );
	void _RheologyMaterial_Execute( void* rheologyMaterial, void* data );
	void _RheologyMaterial_Destroy( void* rheologyMaterial, void* data );
	
	void RheologyMaterial_RunRheologies( 	
		void*                                              rheologyMaterial,
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               swarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi );

	void _RheologyMaterial_RunRheologies( 	
		void*                                              rheologyMaterial,
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               swarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi );

	Bool RheologyMaterial_IsNonLinear( void* rheologyMaterial ) ;

	Rheology* RheologyMaterial_GetRheologyByType( void* rheologyMaterial, Type type ) ;
	Bool RheologyMaterial_HasRheology( void* rheologyMaterial, void* rheology ) ;

#endif
