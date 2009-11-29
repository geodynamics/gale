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
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <StgFEM/FrequentOutput/FrequentOutput.h>

/* Every plugin or component should have a Type (a class name), which specifies the name of the data structure it is */ 
const Type Underworld_MaxVelocity_Type = "Underworld_MaxVelocity";

void Underworld_MaxVelocity_PrintHeaderToFile( void* context );
void Underworld_MaxVelocity_Output( void* _context );
void _Underworld_MaxVelocity_AssignFromXML( void* plugin, Stg_ComponentFactory* cf, void* data );
void* _Underworld_MaxVelocity_DefaultNew( Name name );


Index Underworld_MaxVelocity_Register( PluginsManager* pluginsManager ) {
/*
 * 	Purpose:
 * 		Registers the Plugin, by linking the function _Underworld_MaxVelocity_DefaultNew
 * 		to the place holder (string type) Underworld_MaxVelocity_Type.
 * 		When the place holder is read from the input file the plugin manager needs to
 * 		identify the place holder string to the creation of a plugin data structure. 
 *
 * 	Inputs:
 * 		Inputs are all automatic and this code never needs to be called explicitly
 *
 * 	Interactions:
 * 		It's VERY IMPORTANT that the name of this function (prefix of the function name)
 * 		is spelt exactly the same as the Type place holder string.
 * 		In this case Underworld_MaxVelocity_Register, where the
 * 		prefix = "Underworld_MaxVelocity", is exactly the same as 
 * 		Underworld_MaxVelocity_Type = "Underworld_MaxVelocity".
 * 		
 * 		The function PluginsManager_Submit, needs to be given an Instantiation function
 * 		as it's final argument. This argument should define the functionality of the plugin
 */

	return PluginsManager_Submit( pluginsManager, Underworld_MaxVelocity_Type, "0", _Underworld_MaxVelocity_DefaultNew );
}

void* _Underworld_MaxVelocity_DefaultNew( Name name ) {
/*
 * 	Purpose:
 * 		Registers the 'key' StGermain functions that will be used by this plugin.
 * 		The 'key' functions are called during the different StGermain 'phases'.
 * 		The important 'phases' are: Construct, Build, Initialise, Exectute, Destroy.
 *
 * 	Inputs:
 * 		Inputs are all automatic.
 *
 * 	Interactions:
 * 		The Instantiation funtion calls a New function that contains either 
 * 		functions defined within this plugin OR a parent plugin. Depending on the level
 * 		of complexity of the plugin you may need to defined custom _Build or _Initialise
 * 		functions; in this case the plugin is straight forward and no overriding is
 * 		needed  
 */
	return Codelet_New(
		Underworld_MaxVelocity_Type,
		_Underworld_MaxVelocity_DefaultNew,
		_Underworld_MaxVelocity_AssignFromXML,
		_Codelet_Build,
		_Codelet_Initialise,
		_Codelet_Execute,
		_Codelet_Destroy,
		name );
}

void _Underworld_MaxVelocity_AssignFromXML( void* plugin, Stg_ComponentFactory* cf, void* data ) {
/*
 * 	Purpose:
 * 		This function is called on the 'Construct phase' as defined by the
 * 		Codelet_New(.....) above. The construct phase is called at the beginning of the 
 * 		code. 
 * 		In this function all data structures (objects or components) that this 
 * 		plugin needs should be collected. To be more specific, the objects that
 * 		are defined in the inputfile (and thus are static) should be collected here.
 *
 * 		Also the user functionily (as opposed to non StGermain functionality) 
 * 		should be defined in here. 
 * 		So users defined functions should be attached to EntryPoints. These user
 * 		defined functions are called Hooks. So StGermain executes entry points, which
 * 		kick off user defined Hooks.
 *
 * 	Input:
 * 		Inputs are all automatic. 
 * 		1st args is the plugin pointer itself. (Not used here)
 * 		2nd is the component factory. Will be used to collect components from
 * 			the inputfile.
 * 		3rd Ask Steve.
 *
 * 	Interactions:
 * 		The context is generally gathered first from the ComponentFactory in
 * 		every plugin's construction phase. This allows the plugin we're constructing
 * 		to have access to a wide range of data. (Note this is not the behaviour of components;
 * 		they only access the data structures they need).
 * 		Once the context has been gathered we append a hook to the EP called 
 * 		AbstractContext_EP_FrequentOutput. A list of StGermain entry points can
 * 		be found in StGermain/Base/Context/src/AbstractContext.c
 * 		
 *
 */
	UnderworldContext*  context;

	/* Gather context from the ComponentFactory here */
	context = Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data );

	/* Simple user defined function which belong to no EntryPoint */
	Underworld_MaxVelocity_PrintHeaderToFile( context );

	/* Append user defined function, Underworld_MaxVelocity_Output, onto the EntryPoint
	 * AbstractContext_EP_FrequentOutput */
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, Underworld_MaxVelocity_Output );
}



void Underworld_MaxVelocity_Output( void* _context ) {
/*
 * 	Purpose:
 * 		User defined function. In this case is just prints out the max velocity to 
 * 		the FrequentOutput.dat file.
 *
 * 	Inputs:
 * 		The context, which can be used to control any other object in the code.
 * 		As this function is on the AbstractContext_EP_FrequentOutput it
 * 		receives the input which is defined for that entry point. Other entry points
 * 		have a different set of input args so be sure your user defined plugin is receiving
 * 		the right data structure as the void*. 
 * 		(Note some entry point specify multiple input args).
 *
 * 	Interactions:
 * 		Using the velocityField, which is pointered to from the context, the max velocity
 * 		component value is gathered and piped to the frequent output file. 
 */

	UnderworldContext* context       = (UnderworldContext*) _context;
	FeVariable*        velocityFe    = (FeVariable*) LiveComponentRegister_Get( context->CF->LCRegister, "VelocityField" );
	double             maxVel;

	/* Find the max field component */
	maxVel = _FeVariable_GetMaxGlobalFieldMagnitude( velocityFe );
	/* Print to the FrequentOutput stream */
	StgFEM_FrequentOutput_PrintValue( context, maxVel );
}

void Underworld_MaxVelocity_PrintHeaderToFile( void* context ) {
	StgFEM_FrequentOutput_PrintString( context, "MaxVelocity" );
}



