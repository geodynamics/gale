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

#include <stdlib.h>
#include <string.h>
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
	return Dictionary_GetDouble_WithDefault( context->dictionary, (Dictionary_Entry_Key)"leftWallVelocity", 0.0 );
}	

double GetLeftWallShearVelocity( UnderworldContext* context  ) {
		return Dictionary_GetDouble_WithDefault( context->dictionary, (Dictionary_Entry_Key)"leftWallShearVelocity", 0.0 );
	}
	
double GetRightWallVelocity( UnderworldContext* context  ) {
		return Dictionary_GetDouble_WithDefault( context->dictionary, (Dictionary_Entry_Key)"rightWallVelocity", 0.0 );
	}	
	
double GetRightWallShearVelocity( UnderworldContext* context  ) {
			return Dictionary_GetDouble_WithDefault( context->dictionary, (Dictionary_Entry_Key)"rightWallShearVelocity", 0.0 );
}

double GetBackWallVelocity( UnderworldContext* context ) {
	if ( context->dim == 2  )
		return 0.0;
	
	return Dictionary_GetDouble_WithDefault( context->dictionary, (Dictionary_Entry_Key)"backWallVelocity", 0.0 );
}
double GetFrontWallVelocity( UnderworldContext* context ) {
	if ( context->dim == 2  )
		return 0.0;
	
	return Dictionary_GetDouble_WithDefault( context->dictionary, (Dictionary_Entry_Key)"frontWallVelocity", 0.0 );
}
double GetReferenceHeight( UnderworldContext* context  ) {
	return Dictionary_GetDouble_WithDefault( context->dictionary, (Dictionary_Entry_Key)"constantHeight", 0.0 );
}

double GetTopWallVelocity( UnderworldContext* context ) {
	FeVariable*         velocityField = (FeVariable* ) LiveComponentRegister_Get( context->CF->LCRegister, (Name)"VelocityField" );
	double              y   = GetReferenceHeight( context );
	double              V_a = GetRightWallVelocity( context );
	double              V_b = GetLeftWallVelocity( context );
	double              V_e = GetBackWallVelocity( context );
	double              V_f = GetFrontWallVelocity( context  );
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
	FeVariable*         velocityField = (FeVariable*) LiveComponentRegister_Get( context->CF->LCRegister, (Name)"VelocityField" );
	double              y   = GetReferenceHeight( context );
	double              V_a = GetRightWallVelocity( context );
	double              V_b = GetLeftWallVelocity( context );
	double              V_e = GetBackWallVelocity( context );
	double              V_f = GetFrontWallVelocity( context  );
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

void Underworld_IncompressibleExtensionBC_Remesh( TimeIntegrator* timeIntegrator, IncExtBC* self ) {
    FeVariable* velocityField = (FeVariable*) LiveComponentRegister_Get( self->context->CF->LCRegister, (Name)"VelocityField" );
    FeVariable* pressureField = (FeVariable* ) LiveComponentRegister_Get( self->context->CF->LCRegister, (Name)"PressureField"  );
    FeMesh *mesh;
    Grid *nodeGrid;
    double dt;
    double top, bottom, left, right;
    double minCrd[3], maxCrd[3];
    double nodeWidth[3];
    int numNodes;
    int nodeInds[3];
    int ii;

    mesh = velocityField->feMesh;

    nodeGrid = *Mesh_GetExtension(mesh, Grid**, "vertexGrid");

    dt = AbstractContext_Dt(self->ctx);
    top = GetTopWallVelocity(self->ctx);
    bottom = GetBottomWallVelocity(self->ctx);
    left = GetLeftWallVelocity(self->ctx);
    right = GetRightWallVelocity(self->ctx);
    Mesh_GetGlobalCoordRange(mesh, minCrd, maxCrd);

    minCrd[0] += dt*left;
    minCrd[1] += dt*bottom;
    maxCrd[0] += dt*right;
    maxCrd[1] += dt*top;

    nodeWidth[0] = (maxCrd[0] - minCrd[0])/(double)(nodeGrid->sizes[0] - 1);
    nodeWidth[1] = (maxCrd[1] - minCrd[1])/(double)(nodeGrid->sizes[1] - 1);

    numNodes = FeMesh_GetNodeLocalSize(mesh);
    for(ii = 0; ii < numNodes; ii++) {
		Grid_Lift(nodeGrid, FeMesh_NodeDomainToGlobal(mesh, ii), nodeInds);
		Mesh_GetVertex(mesh, ii)[0] = minCrd[0] + ((double)nodeInds[0])*nodeWidth[0];
		Mesh_GetVertex(mesh, ii)[1] = minCrd[1] + ((double)nodeInds[1])*nodeWidth[1];
    }

    Mesh_Sync(mesh);
    Mesh_DeformationUpdate(mesh);

    if(velocityField->feMesh != pressureField->feMesh) {
	if(!strcmp(pressureField->feMesh->type, "CartesianGenerator")) {
	    FeMesh *pmesh;
	    Grid *pnodeGrid;
	    int lind;

	    pmesh = pressureField->feMesh;
	    pnodeGrid = *Mesh_GetExtension(pmesh, Grid**, "vertexGrid");

	    numNodes = FeMesh_GetNodeLocalSize(pmesh);
	    for(ii = 0; ii < numNodes; ii++) {
		Grid_Lift(pnodeGrid, FeMesh_NodeDomainToGlobal(pmesh, ii), nodeInds);
		nodeInds[0] *= 2;
		nodeInds[1] *= 2;
		insist(FeMesh_NodeGlobalToDomain(mesh, Grid_Project(nodeGrid, nodeInds), &lind), != 0);
		memcpy(Mesh_GetVertex(pmesh, ii), Mesh_GetVertex(mesh, lind), 2*sizeof(double));
	    }
	}
    }
}

void _Underworld_IncompressibleExtensionBC_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data ) {
    IncExtBC* self = (IncExtBC*)_self;
	UnderworldContext*  context  = Stg_ComponentFactory_ConstructByName( cf, (Name)"context", UnderworldContext, True, data );
	ConditionFunction*  condFunc;
   TimeIntegrator* timeIntegrator = (TimeIntegrator* )  LiveComponentRegister_Get( context->CF->LCRegister, (Name)"timeIntegrator"  );

        self->context   = context;

	condFunc = ConditionFunction_New( IncompressibleExtensionBC_TopCondition, (Name)"IncompressibleExtensionBC_TopCondition"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( IncompressibleExtensionBC_BottomCondition, (Name)"IncompressibleExtensionBC_BottomCondition"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( IncompressibleExtensionBC_LeftCondition, (Name)"IncompressibleExtensionBC_LeftCondition"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( IncompressibleExtensionBC_LeftShearCondition, (Name)"IncompressibleExtensionBC_LeftShearCondition"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( IncompressibleExtensionBC_RightCondition, (Name)"IncompressibleExtensionBC_RightCondition"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( IncompressibleExtensionBC_RightShearCondition, (Name)"IncompressibleExtensionBC_RightShearCondition"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( IncompressibleExtensionBC_FrontCondition, (Name)"IncompressibleExtensionBC_FrontCondition"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( IncompressibleExtensionBC_BackCondition, (Name)"IncompressibleExtensionBC_BackCondition"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

        if( Stg_ComponentFactory_PluginGetBool( cf, self, (Dictionary_Entry_Key)"Remesh", False )  ) {
	    TimeIntegrator_PrependFinishEP( 
		timeIntegrator, "Underworld_IncompressibleExtensionBC_Remesh", Underworld_IncompressibleExtensionBC_Remesh, 
		CURR_MODULE_NAME, self );
	}
}


/* This function will provide StGermain the abilty to instantiate (create) this codelet on demand. */
void* _Underworld_IncompressibleExtensionBC_DefaultNew( Name name ) {
   return _Codelet_New(
			sizeof( IncExtBC ),
			Underworld_IncompressibleExtensionBC_Type,
			_Codelet_Delete,
			_Codelet_Print, 
			_Codelet_Copy,
                        name,
                        NON_GLOBAL,
			_Underworld_IncompressibleExtensionBC_DefaultNew,
			_Underworld_IncompressibleExtensionBC_AssignFromXML,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy );
}
	
/* This function is automatically run by StGermain when this plugin is loaded. The name must be "<plugin-name>_Register". */
Index Underworld_IncompressibleExtensionBC_Register( PluginsManager* pluginsManager ) {
	/* A plugin is only properly registered once it returns the handle provided when submitting a codelet to StGermain. */
	return PluginsManager_Submit( pluginsManager, Underworld_IncompressibleExtensionBC_Type, (Name)"0", _Underworld_IncompressibleExtensionBC_DefaultNew  );
}


