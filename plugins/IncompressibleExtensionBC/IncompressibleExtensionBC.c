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
** $Id: IncompressibleExtensionBC.c 466 2007-04-27 06:24:33Z LukeHodkinson $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <assert.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

/* Each component type needs a unique identifier (as a string for us to read and as an integer for realtime comparisions) */
const Type Underworld_IncompressibleExtensionBC_Type = "Underworld_IncompressibleExtensionBC_Type";


void CalculateVelocities( UnderworldContext* context, double* V_c, double* V_d )  {
	FeVariable*         velocityField  = context->velocityField;
	FeMesh* mesh           = velocityField->feMesh;
	double              V_a, V_b;
	double              h_1, h_2;
	double              x;
	IJK                 leftWall_global_IJK;
	IJK                 rightWall_global_IJK;
	Node_GlobalIndex    leftWall_globalIndex;
	Node_GlobalIndex    rightWall_globalIndex;
	Node_LocalIndex     leftWall_localIndex;
	Node_LocalIndex     rightWall_localIndex;
	XYZ                 velocity;
	XYZ                 min, max;
	double              width;
	Grid*		    vertGrid;

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

	x = Dictionary_GetDouble_WithDefault( context->dictionary, "constantHeight", 0.0 );

	vertGrid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
						  ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

	/* Grab V_a and V_b (i.e. velocities at left and right walls) */
	leftWall_global_IJK[ I_AXIS ] = 0;
	leftWall_global_IJK[ J_AXIS ] = 0;
	leftWall_global_IJK[ K_AXIS ] = 0;

	rightWall_global_IJK[ I_AXIS ] = vertGrid->sizes[ I_AXIS ] - 1;
	rightWall_global_IJK[ J_AXIS ] = 0;
	rightWall_global_IJK[ K_AXIS ] = 0;
	
	/* Grab global indicies for these nodes */
	leftWall_globalIndex = Grid_Project( vertGrid, leftWall_global_IJK );
	rightWall_globalIndex = Grid_Project( vertGrid, rightWall_global_IJK );
	
	/* Grab local indicies for these nodes */
	insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, leftWall_globalIndex, &leftWall_localIndex ) );
	insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, rightWall_globalIndex, &rightWall_localIndex ) );

	/* Check if the left wall is on processor */
	if ( leftWall_localIndex < FeMesh_GetNodeLocalSize( mesh ) ) {
		FeVariable_GetValueAtNode( velocityField, leftWall_localIndex, velocity );
		V_b = velocity[ I_AXIS ];
	}
	/* Check if the right wall is on processor */
	if ( rightWall_localIndex < FeMesh_GetNodeLocalSize( mesh ) ) {
		FeVariable_GetValueAtNode( velocityField, rightWall_localIndex, velocity );
		V_a = velocity[ I_AXIS ];
	}

	/* Send these values across processors  - TODO */
	
	/* Calculate Width and Height */
	FieldVariable_GetMinAndMaxGlobalCoords( velocityField, min, max );
	width  = (max[ I_AXIS ] - min[ I_AXIS ]);
	h_1 = max[ J_AXIS ] - x;
	h_2 = x - min[ J_AXIS ];

	/* Calculate velocity at the top and at the bottom of the nodes */
	*V_c = - h_1/width * ( V_a - V_b );
	*V_d =   h_2/width * ( V_a - V_b );

}

void IncompressibleExtensionBC_TopCondition( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	UnderworldContext* context = (UnderworldContext*) _context;
	double* result = (double*) _result;
	double tmp;

	CalculateVelocities( context, result, &tmp );
}

void IncompressibleExtensionBC_BottomCondition( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	UnderworldContext* context = (UnderworldContext*) _context;
	double* result = (double*) _result;
	double tmp;

	CalculateVelocities( context, &tmp, (double*) result );
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
