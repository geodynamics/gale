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
** $Id: testStiffnessMatrixAssembly.c 430 2007-02-07 00:10:36Z PatrickSunter $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include <Underworld/Rheology/Rheology.h>

#include <assert.h>

const Type Underworld_testStiffnessMatrixAssembly_Type = "Underworld_testStiffnessMatrixAssembly";

typedef struct {
	__Codelet
} Underworld_testStiffnessMatrixAssembly;

void Underworld_testStiffnessMatrixAssemblyFunction( FiniteElementContext* context ) {
	Stream*                     stream     = Journal_Register( Info_Type, "testConstitutiveMatrix" );
	StiffnessMatrix*            stiffnessMatrix;
	SystemLinearEquations*      sle;

	Stream_SetPrintingRank( stream, 0 );
	
	/* Assemble */
	stiffnessMatrix = (StiffnessMatrix*) LiveComponentRegister_Get( context->CF->LCRegister, "k_matrix" );
	sle = (SystemLinearEquations*)       LiveComponentRegister_Get( context->CF->LCRegister, "stokesEqn" );
	assert( stiffnessMatrix );
	StiffnessMatrix_Assemble( stiffnessMatrix, False, sle, context );
}

void _Underworld_testStiffnessMatrixAssembly_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	FiniteElementContext* context = (FiniteElementContext*)Stg_ComponentFactory_ConstructByName( 
		cf, 
		"context", 
		FiniteElementContext, 
		True,
		data );
	ContextEP_ReplaceAll( context, AbstractContext_EP_Solve, Underworld_testStiffnessMatrixAssemblyFunction );
}

void* _Underworld_testStiffnessMatrixAssembly_DefaultNew( Name name ) {
	return Codelet_New(
		Underworld_testStiffnessMatrixAssembly_Type,
		_Underworld_testStiffnessMatrixAssembly_DefaultNew,
		_Underworld_testStiffnessMatrixAssembly_Construct,
		_Codelet_Build,
		_Codelet_Initialise,
		_Codelet_Execute,
		_Codelet_Destroy,
		name );
}

Index Underworld_testStiffnessMatrixAssembly_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Underworld_testStiffnessMatrixAssembly_Type, "0",
	_Underworld_testStiffnessMatrixAssembly_DefaultNew );
}
