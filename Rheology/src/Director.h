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
*+      Kathleen Humble
*+
** $Id: Director.h 354 2006-10-12 08:19:27Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __Underworld_Rheology_Director_h__
#define __Underworld_Rheology_Director_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type Director_Type;
    typedef enum InitialDirectionType { INIT_DIR_GLOBAL=0, INIT_DIR_RANDOM=1, INIT_DIR_PER_MAT=2 } InitialDirectionType;
	
	typedef struct {
		XYZ            director;
		Particle_Bool  dontUpdateParticle;
	}  Director_ParticleExt; 

	/** Rheology class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __Director \
		/* Parent info */ \
 		__TimeIntegratee \
		/* Virtual functions go here */ \
		/* General Info */\
		SwarmVariable*                                      directorSwarmVariable;                 \
		SwarmVariable*                                      dontUpdateParticle;                 \		
		ExtensionInfo_Index                                 particleExtHandle;                     \
		/* Param passed in */ \
		FeVariable*                                         velGradField;                          \
		MaterialPointsSwarm*                                materialPointsSwarm;                  \
		Materials_Register*                                 materialsRegister;                  \
		InitialDirectionType                                initialDirectionType;               \
		XYZ                                                 globalInitialDirection;                      \
		unsigned int                                        randomInitialDirectionSeed;             \
		Bool                                                dontUpdate; 
				
	struct Director { __Director };

   /** Public Constructor */
   Director* Director_New( 
		Name                   name,
		DomainContext*         context,
		TimeIntegrator*        timeIntegrator, 
		Variable*              variable,
		Index                  dataCount, 
		Stg_Component**        data,
		Bool                   allowFallbackToFirstOrder,
		FeVariable*            velGradField,
		MaterialPointsSwarm*   materialPointsSwarm,
		InitialDirectionType   initialDirectionType,
		double                 globalInitialDirectionX,
		double                 globalInitialDirectionY,
		double                 globalInitialDirectionZ,
		int                    randomInitialDirectionSeed,
		Bool                   dontUpdate );

	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	Director* _Director_New( 
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
		TimeIntegratee_CalculateTimeDerivFunction*         _calculateTimeDeriv,
		TimeIntegratee_IntermediateFunction*               _intermediate,
		Name                                               name ) ;
	
	/* 'Stg_Component' implementations */
	void* _Director_DefaultNew( Name name );
	void _Director_AssignFromXML( void* director, Stg_ComponentFactory* cf, void* data );
	void _Director_Build( void* director, void* data );
	void _Director_Initialise( void* director, void* data );
   void _Director_Init(
         Director*                                          self,
         FeVariable*                                        velGradField,
         MaterialPointsSwarm*                               materialPointsSwarm,
         InitialDirectionType                               initialDirectionType,
         double                                             globalInitialDirectionX,
         double                                             globalInitialDirectionY,
         double                                             globalInitialDirectionZ,
         int                                                randomInitialDirectionSeed,
         Bool                                               dontUpdate );
	void _Director_Destroy( void* _self, void* data ) ;
	
	Bool _Director_TimeDerivative( void* _director, Index lParticle_I, double* timeDeriv );
	void _Director_Intermediate( void* director, Index lParticle_I );

	void Director_UpdateVariables( void* timeIntegrator, Director* self );

	void Director_GetNormal( void* director, void* particle, XYZ normal );
	
	void Director_SetNormal( void* director, void* particle, XYZ normal );

	double* Director_GetNormalPtr( void* director, void* particle);

	void Director_SetDontUpdateParticleFlag( void* director, void* particle, Particle_Bool dontUpdateFlag );
	
#endif
