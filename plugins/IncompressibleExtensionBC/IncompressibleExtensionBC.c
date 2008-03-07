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
** $Id: IncompressibleExtensionBC.c 675 2008-03-07 06:32:38Z JulianGiordani $
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
	                 ^ V_c
	 ________________|_____________________ 
	|     ^                                | --
	|     |                                |   |
	|     h_1                              |    }
	|     |                                |    }
	|     v                                |   |
	|--------------------------------------| --
	|         ^                            |
	|         |                            |
	|-> V_b   |                            |->V_a
	|         |                            |
	|         |                            |
	|         |      ^V_d                  |
	|________________|_____________________|

	*/

double GetLeftWallVelocity( FeVariable* velocityField ) {
	FeMesh*             mesh = velocityField->feMesh;
	IJK                 globalIJK = {0,0,0};
	Node_GlobalIndex    global_I;
	XYZ                 velocity;
	Grid*		    vertGrid;
	
	vertGrid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

	/* Grab global indicies for these nodes */
	global_I = Grid_Project( vertGrid, globalIJK );
	
	/* Grab velocities on these walls */
	FeVariable_GetValueAtNodeGlobal( velocityField, global_I, velocity );

	/* Get components in horizontal direction */
	return velocity[ I_AXIS ]; 
}
double GetRightWallVelocity( FeVariable* velocityField ) {
	FeMesh*             mesh = velocityField->feMesh;
	IJK                 globalIJK = {0,0,0};
	Node_GlobalIndex    global_I;
	XYZ                 velocity;
	Grid*		    vertGrid;
	
	vertGrid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

	globalIJK[ I_AXIS ] = vertGrid->sizes[ I_AXIS ] - 1;

	/* Grab global indicies for these nodes */
	global_I = Grid_Project( vertGrid, globalIJK );
	
	/* Grab velocities on these walls */
	FeVariable_GetValueAtNodeGlobal( velocityField, global_I, velocity );
	
	/* Get components in horizontal direction */
	return velocity[ I_AXIS ]; 
}

double GetTopWallVelocity( FeVariable* velocityField, double y ) {
	double              V_b = GetLeftWallVelocity( velocityField );
	double              V_a = GetRightWallVelocity( velocityField );
	double              V_c;
	double              h_1;
	XYZ                 min, max;
	double              width;
	
	/* Calculate Width and Height */
	FieldVariable_GetMinAndMaxGlobalCoords( velocityField, min, max );
	width  = (max[ I_AXIS ] - min[ I_AXIS ]);
	h_1 = max[ J_AXIS ] - y;

	/* Calculate velocity at the top and at the bottom of the nodes */
	V_c = - h_1/width * ( V_a - V_b );

	return V_c;
}

double GetBottomWallVelocity( FeVariable* velocityField, double y ) {
	double              V_b = GetLeftWallVelocity( velocityField );
	double              V_a = GetRightWallVelocity( velocityField );
	double              V_d;
	double              h_2;
	XYZ                 min, max;
	double              width;
	
	/* Calculate Width and Height */
	FieldVariable_GetMinAndMaxGlobalCoords( velocityField, min, max );
	width  = (max[ I_AXIS ] - min[ I_AXIS ]);
	h_2 = y - min[ J_AXIS ];
	
	/* Calculate velocity at the top and at the bottom of the nodes */
	V_d =   h_2/width * ( V_a - V_b );
	
	return V_d;
}

void GetVelocity( FeVariable* velocityField, double y, Coord coord, double* velocity ) {
	double              V_a = GetRightWallVelocity( velocityField );
	double              V_b = GetLeftWallVelocity( velocityField );
	double              V_c = GetTopWallVelocity( velocityField, y );
	double              V_d = GetBottomWallVelocity( velocityField, y );
	XYZ                 min, max;
	double              width;
	double              height;
	
	printf(" %g %g \n", coord[0], coord[1] );

	FieldVariable_GetMinAndMaxGlobalCoords( velocityField, min, max );
	width  = (max[ I_AXIS ] - min[ I_AXIS ]);
	height  = (max[ J_AXIS ] - min[ J_AXIS ]);

	velocity[ I_AXIS ] = ( coord[ I_AXIS ] - min[ I_AXIS ] ) / width  * ( V_a - V_b ) + V_b;
	velocity[ J_AXIS ] = ( coord[ J_AXIS ] - min[ J_AXIS ] ) / height * ( V_c - V_d ) + V_d;
	if ( velocityField->dim == 3 )
		velocity[ K_AXIS ] = 0.0;
}

void IncompressibleExtensionBC_TopCondition( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	UnderworldContext* context = (UnderworldContext*) _context;
	double*            result  = (double*) _result;
	double             y       = Dictionary_GetDouble_WithDefault( context->dictionary, "constantHeight", 0.0 );

	*result = GetTopWallVelocity( context->velocityField, y );
}

void IncompressibleExtensionBC_BottomCondition( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	UnderworldContext* context = (UnderworldContext*) _context;
	double*            result  = (double*) _result;
	double             y       = Dictionary_GetDouble_WithDefault( context->dictionary, "constantHeight", 0.0 );

	*result = GetBottomWallVelocity( context->velocityField, y );
}

void _Underworld_IncompressibleExtensionBC_Construct( void* self, Stg_ComponentFactory* cf, void* data ) {
	UnderworldContext*  context  = Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data );
	ConditionFunction*  condFunc;

	condFunc = ConditionFunction_New( IncompressibleExtensionBC_TopCondition, "IncompressibleExtensionBC_TopCondition" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( IncompressibleExtensionBC_BottomCondition, "IncompressibleExtensionBC_BottomCondition" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );

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
