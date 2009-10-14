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
** $Id: IncompressibleExtensionBC.c 728 2008-05-12 02:29:30Z LouisMoresi $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <assert.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include "IncompressibleExtensionBC.h"

/* Each component type needs a unique identifier (as a string for us to read and as an integer for realtime comparisions) */
const Type Underworld_IncompressibleExtensionBC_Type = "Underworld_IncompressibleExtensionBC_Type";

	/*

	                       ^ V_e
	                      /
	         ____________/__________________________
	        /                                      /|
	       /                                      / |
	      /                                      /  |
	     /                  ^ V_c               /   |
	    /                   |                  /    |
	   /          ^ V_f                       /     |
	  /          /                           /      |
	 /__________/___________________________/       |
	|     ^                                |        |
	|     |                                |        |
	|     h_1                              |  ->V_a |
	|     |                                |        |
	|     v                                |        |
	|--------------------------------------|       /
	|         ^                            |      /
	|         |                            |     /
	|-> V_b   |                            |    /
	|         |h_2                         |   / 
	|         |                            |  /
	|         |      ^V_d                  | /
	|________________|_____________________|/

	*/

double GetLeftWallVelocity( UnderworldContext* context ) {
	return Dictionary_GetDouble_WithDefault( context->dictionary,  "leftWallVelocity", 0.0 );
}	

double GetLeftWallShearVelocity( UnderworldContext* context ) {
		return Dictionary_GetDouble_WithDefault( context->dictionary,  "leftWallShearVelocity", 0.0 );
	}
	
double GetRightWallVelocity( UnderworldContext* context ) {
		return Dictionary_GetDouble_WithDefault( context->dictionary,  "rightWallVelocity", 0.0 );
	}	
	
double GetRightWallShearVelocity( UnderworldContext* context ) {
			return Dictionary_GetDouble_WithDefault( context->dictionary,  "rightWallShearVelocity", 0.0 );
}

double GetBackWallVelocity( UnderworldContext* context ) {
	if ( context->dim == 2 )
		return 0.0;
	
	return Dictionary_GetDouble_WithDefault( context->dictionary,  "backWallVelocity", 0.0 );
}
double GetFrontWallVelocity( UnderworldContext* context ) {
	if ( context->dim == 2 )
		return 0.0;
	
	return Dictionary_GetDouble_WithDefault( context->dictionary,  "frontWallVelocity", 0.0 );
}
double GetReferenceHeight( UnderworldContext* context ) {
	return Dictionary_GetDouble_WithDefault( context->dictionary,  "constantHeight", 0.0 );
}

double GetTopWallVelocity( UnderworldContext* context ) {
	FeVariable*         velocityField = (FeVariable*) LiveComponentRegister_Get( context->CF->LCRegister, "velocityField" );
	double              y   = GetReferenceHeight( context );
	double              V_a = GetRightWallVelocity( context );
	double              V_b = GetLeftWallVelocity( context );
	double              V_e = GetBackWallVelocity( context );
	double              V_f = GetFrontWallVelocity( context );
	double              V_c;
	double              h_1;
	XYZ                 min, max;
	double              width;
	double              depth;
	
	/* Calculate Width and Height */
	FieldVariable_GetMinAndMaxGlobalCoords( velocityField, min, max );
	width  = (max[ I_AXIS ] - min[ I_AXIS ]);
	depth  = context->dim == 3 ? (max[ K_AXIS ] - min[ K_AXIS ]) : 1.0; /* if only 3D depth cancels in division */

	h_1 = max[ J_AXIS ] - y;

	/* Calculate velocity at the top and at the bottom of the nodes */
	V_c = - h_1 * ( (V_a - V_b) * depth + (V_f - V_e) * width )/(width*depth);
	
	return V_c;
}

double GetBottomWallVelocity( UnderworldContext* context ) {
	FeVariable*         velocityField = (FeVariable*) LiveComponentRegister_Get( context->CF->LCRegister, "velocityField" );
	double              y   = GetReferenceHeight( context );
	double              V_a = GetRightWallVelocity( context );
	double              V_b = GetLeftWallVelocity( context );
	double              V_e = GetBackWallVelocity( context );
	double              V_f = GetFrontWallVelocity( context );
	double              V_d;
	double              h_2;
	XYZ                 min, max;
	double              width;
	double              depth;
	
	/* Calculate Width and Height */
	FieldVariable_GetMinAndMaxGlobalCoords( velocityField, min, max );
	width  = (max[ I_AXIS ] - min[ I_AXIS ]);
	depth  = context->dim == 3 ? (max[ K_AXIS ] - min[ K_AXIS ]) : 1.0; /* if only 3D depth cancels in division */

	h_2 = y - min[ J_AXIS ];
	
	/* Calculate velocity at the top and at the bottom of the nodes */
	V_d =   h_2 * ( (V_a - V_b) * depth + (V_f - V_e) * width )/(width*depth);

	return V_d;
}

void IncompressibleExtensionBC_RightCondition( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	UnderworldContext* context = (UnderworldContext*) _context;
	double*            result  = (double*) _result;

	*result = GetRightWallVelocity( context );
}

void IncompressibleExtensionBC_RightShearCondition( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	UnderworldContext* context = (UnderworldContext*) _context;
	double*            result  = (double*) _result;

	*result = GetRightWallShearVelocity( context );
}


void IncompressibleExtensionBC_LeftCondition( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	UnderworldContext* context = (UnderworldContext*) _context;
	double*            result  = (double*) _result;

	*result = GetLeftWallVelocity( context );
}

void IncompressibleExtensionBC_LeftShearCondition( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	UnderworldContext* context = (UnderworldContext*) _context;
	double*            result  = (double*) _result;

	*result = GetLeftWallShearVelocity( context );
}


void IncompressibleExtensionBC_BackCondition( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	UnderworldContext* context = (UnderworldContext*) _context;
	double*            result  = (double*) _result;

	*result = GetBackWallVelocity( context );
}
void IncompressibleExtensionBC_FrontCondition( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	UnderworldContext* context = (UnderworldContext*) _context;
	double*            result  = (double*) _result;

	*result = GetFrontWallVelocity( context );
}

void IncompressibleExtensionBC_TopCondition( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	UnderworldContext* context = (UnderworldContext*) _context;
	double*            result  = (double*) _result;

	*result = GetTopWallVelocity( context );
}

void IncompressibleExtensionBC_BottomCondition( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	UnderworldContext* context = (UnderworldContext*) _context;
	double*            result  = (double*) _result;

	*result = GetBottomWallVelocity( context );
}

void _Underworld_IncompressibleExtensionBC_Construct( void* self, Stg_ComponentFactory* cf, void* data ) {
	UnderworldContext*  context  = Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data );
	ConditionFunction*  condFunc;

	condFunc = ConditionFunction_New( IncompressibleExtensionBC_TopCondition, "IncompressibleExtensionBC_TopCondition" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( IncompressibleExtensionBC_BottomCondition, "IncompressibleExtensionBC_BottomCondition" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( IncompressibleExtensionBC_LeftCondition, "IncompressibleExtensionBC_LeftCondition" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( IncompressibleExtensionBC_LeftShearCondition, "IncompressibleExtensionBC_LeftShearCondition" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( IncompressibleExtensionBC_RightCondition, "IncompressibleExtensionBC_RightCondition" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( IncompressibleExtensionBC_RightShearCondition, "IncompressibleExtensionBC_RightShearCondition" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( IncompressibleExtensionBC_FrontCondition, "IncompressibleExtensionBC_FrontCondition" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( IncompressibleExtensionBC_BackCondition, "IncompressibleExtensionBC_BackCondition" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

}


/* This function will provide StGermain the abilty to instantiate (create) this codelet on demand. */
void* _Underworld_IncompressibleExtensionBC_DefaultNew( Name name ) {
	return Codelet_New(
			Underworld_IncompressibleExtensionBC_Type,
			_Underworld_IncompressibleExtensionBC_DefaultNew,
			_Underworld_IncompressibleExtensionBC_Construct, /* SQ NOTE: Used to be a construct extensions. */
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}
	
/* This function is automatically run by StGermain when this plugin is loaded. The name must be "<plugin-name>_Register". */
Index Underworld_IncompressibleExtensionBC_Register( PluginsManager* pluginsManager ) {
	/* A plugin is only properly registered once it returns the handle provided when submitting a codelet to StGermain. */
	return PluginsManager_Submit( pluginsManager, Underworld_IncompressibleExtensionBC_Type, "0", _Underworld_IncompressibleExtensionBC_DefaultNew );
}
