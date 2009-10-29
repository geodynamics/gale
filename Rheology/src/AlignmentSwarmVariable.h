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
** $Id: AlignmentSwarmVariable.h 354 2006-10-12 08:19:27Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __Underworld_Rheology_AlignmentSwarmVariable_h__
#define __Underworld_Rheology_AlignmentSwarmVariable_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type AlignmentSwarmVariable_Type;

	/** Class contents */
	#define __AlignmentSwarmVariable \
		/* Parent info */ \
 		__SwarmVariable \
		/* Virtual functions go here */ \
		/* General Info */\
		Director*                                           director;                              \
		FeVariable*                                         velocityField;
				
	struct AlignmentSwarmVariable { __AlignmentSwarmVariable };

	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	AlignmentSwarmVariable* _AlignmentSwarmVariable_New( 
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
		SwarmVariable_ValueAtFunction*                     _valueAt,
		SwarmVariable_GetGlobalValueFunction*              _getMinGlobalMagnitude,
		SwarmVariable_GetGlobalValueFunction*              _getMaxGlobalMagnitude,
		Name                                               name ) ;
	
	/* 'Stg_Component' implementations */
	void* _AlignmentSwarmVariable_DefaultNew( Name name ) ;
	void _AlignmentSwarmVariable_AssignFromXML( void* alignment, Stg_ComponentFactory* cf, void* data );
	void _AlignmentSwarmVariable_Initialise( void* alignment, void* data ) ;

	/* 'SwarmVariable Virtual Implementations */
	void _AlignmentSwarmVariable_ValueAt( void* alignment, Particle_Index lParticle_I, double* value ) ;
	double _AlignmentSwarmVariable_GetMinGlobalMagnitude( void* alignment ) ;
	double _AlignmentSwarmVariable_GetMaxGlobalMagnitude( void* alignment ) ;

#endif
