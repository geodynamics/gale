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
** $Id: Vrms.c 182 2006-05-01 12:32:01Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <StgFEM/FrequentOutput/FrequentOutput.h>

const Type Underworld_MaxVelocity_Type = "Underworld_MaxVelocity";
void Underworld_MaxVelocity_PrintHeaderToFile( void* context );
void Underworld_MaxVelocity_Output( void* _context );

void _Underworld_MaxVelocity_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	UnderworldContext*  context;

	context = Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data );

	Underworld_MaxVelocity_PrintHeaderToFile( context );
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, Underworld_MaxVelocity_Output );
}

void* _Underworld_MaxVelocity_DefaultNew( Name name ) {
	return Codelet_New(
		Underworld_MaxVelocity_Type,
		_Underworld_MaxVelocity_DefaultNew,
		_Underworld_MaxVelocity_Construct,
		_Codelet_Build,
		_Codelet_Initialise,
		_Codelet_Execute,
		_Codelet_Destroy,
		name );
}

Index Underworld_MaxVelocity_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Underworld_MaxVelocity_Type, "0", _Underworld_MaxVelocity_DefaultNew );
}

void Underworld_MaxVelocity_Output( void* _context ) {
	UnderworldContext* context       = (UnderworldContext*) _context;
	FeVariable*        velocityFe = context->velocityField;
	double             maxVel;

	maxVel = _FeVariable_GetMaxGlobalFieldMagnitude( velocityFe );
	StgFEM_FrequentOutput_PrintValue( context, maxVel );
}

void Underworld_MaxVelocity_PrintHeaderToFile( void* context ) {
	StgFEM_FrequentOutput_PrintString( context, "MaxVelocity" );
}

