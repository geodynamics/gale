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
** $Id: MovingMeshEnergyCorrection.c 466 2007-04-27 06:24:33Z LukeHodkinson $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>


/* Each component type needs a unique identifier (as a string for us to read and as an integer for realtime comparisions) */
const Type Underworld_MovingMeshEnergyCorrection_Type = "Underworld_MovingMeshEnergyCorrection_Type";

ExtensionInfo_Index Underworld_MovingMeshEnergyCorrection_ContextExtHandle;


typedef struct {
	/* Save the "old" execute function pointer (which we replace) */ 
	Stg_Component_ExecuteFunction* energySolverExecute;
} Underworld_MovingMeshEnergyCorrection_ContextExt;



void FeVariable_GetSideWallValues( void* feVariable, Dimension_Index planeAxis, double** valueMinSideWall, double** valueMaxSideWall ) 
{
	FeVariable*                self          = (FeVariable*)  feVariable;
	Mesh*                      mesh          = (Mesh*)        self->feMesh;
	Grid*			   vertGrid;
	/* Assume IJK Topology */

	double*                    valueMinSideWallLocal;
	double*                    valueMaxSideWallLocal;
	double*                    valueMinSideWallGlobal;
	double*                    valueMaxSideWallGlobal;
	Dimension_Index            aAxis              = ( planeAxis == I_AXIS ? J_AXIS : I_AXIS );
	Dimension_Index            bAxis              = ( planeAxis == K_AXIS ? J_AXIS : K_AXIS );
	IJK                        minNodeIJK;
	IJK                        maxNodeIJK;
	Index                      array_I;
	Node_Index                 aNode_I;
	Node_Index                 bNode_I;
	Node_GlobalIndex           minNodeGlobal_I;
	Node_GlobalIndex           maxNodeGlobal_I;
	Node_LocalIndex            minNodeLocal_I;
	Node_LocalIndex            maxNodeLocal_I;
	Node_Index                 sideNodeCount;
	MPI_Comm                   comm;
	Dof_Index                  dof            = self->fieldComponentCount;
	Dof_Index                  dof_I;

	comm = CommTopology_GetComm( Mesh_GetCommTopology( mesh, MT_VERTEX ) );
	vertGrid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
						  ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );
	
	sideNodeCount = vertGrid->sizes[ aAxis ] * vertGrid->sizes[ bAxis ];

	valueMinSideWallLocal = Memory_Alloc_Array( double, sideNodeCount*dof, "valueMinSideWallLocal" );
	valueMaxSideWallLocal = Memory_Alloc_Array( double, sideNodeCount*dof, "valueMaxSideWallLocal" );
	valueMinSideWallGlobal = Memory_Alloc_Array( double, sideNodeCount*dof, "valueMinSideWallGlobal" );
	valueMaxSideWallGlobal = Memory_Alloc_Array( double, sideNodeCount*dof, "valueMaxSideWallGlobal" );

	/* Loop over all nodes on the minimum side wall and the maximum side wall */
	minNodeIJK[ planeAxis ] = 0;
	maxNodeIJK[ planeAxis ] = vertGrid->sizes[ planeAxis ] - 1 ;
	array_I = 0;
	for ( aNode_I = 0 ; aNode_I < vertGrid->sizes[ aAxis ] ; aNode_I++ ) {
		minNodeIJK[ aAxis ] = maxNodeIJK[ aAxis ] = aNode_I;
		for ( bNode_I = 0 ; bNode_I < vertGrid->sizes[ bAxis ] ; bNode_I++ ) {
			minNodeIJK[ bAxis ] = maxNodeIJK[ bAxis ] = bNode_I;

			/* Get Global Node Indicies for the node at these walls */
			minNodeGlobal_I = Grid_Project( vertGrid, minNodeIJK );
			maxNodeGlobal_I = Grid_Project( vertGrid, maxNodeIJK );

			/* Check to see whether this min node is on this processor */
			if ( Mesh_GlobalToDomain( mesh, MT_VERTEX, minNodeGlobal_I, &minNodeGlobal_I ) && 
			     minNodeLocal_I < Mesh_GetLocalSize( mesh, MT_VERTEX ) )
			{
				FeVariable_GetValueAtNode( feVariable, minNodeLocal_I, &valueMinSideWallLocal[ array_I*dof ] );
			}
			else {
				for ( dof_I = 0 ; dof_I < dof ; dof_I++ ) 
					/* Initialise all values to some really large value so that 
					 * MPI_Allreduce can find value from correct processor */
					valueMinSideWallLocal[ array_I*dof + dof_I ] = HUGE_VAL;
			}
			
			/* Check to see whether this max node is on this processor */
			if ( Mesh_GlobalToDomain( mesh, MT_VERTEX, maxNodeGlobal_I, &maxNodeGlobal_I ) && 
			     maxNodeLocal_I < Mesh_GetLocalSize( mesh, MT_VERTEX ) )
			{
				FeVariable_GetValueAtNode( feVariable, maxNodeLocal_I, &valueMaxSideWallLocal[ array_I*dof ] );
			}
			else {
				for ( dof_I = 0 ; dof_I < dof ; dof_I++ ) 
					/* Initialise all values to some really large value so that 
					 * MPI_Allreduce can find value from correct processor */
					valueMaxSideWallLocal[ array_I*dof + dof_I ] = HUGE_VAL;
			}
			
			/* Increment counter for side wall node index */
			array_I++;
		}
	}

	MPI_Allreduce( valueMinSideWallLocal, valueMinSideWallGlobal, sideNodeCount*dof, MPI_DOUBLE, MPI_MIN, comm );
	MPI_Allreduce( valueMaxSideWallLocal, valueMaxSideWallGlobal, sideNodeCount*dof, MPI_DOUBLE, MPI_MIN, comm );

	Memory_Free( valueMinSideWallLocal );
	Memory_Free( valueMaxSideWallLocal );
	
	/* Set return values */
	*valueMinSideWall = valueMinSideWallGlobal;
	*valueMaxSideWall = valueMaxSideWallGlobal;
}

typedef enum {
	ADD,
	MINUS
} CorrectionFlag;

void MovingMeshEnergyCorrection_AddCorrection( FeVariable* velocityField, double* xLeftVelocities, double* xRightVelocities, CorrectionFlag flag ) {
	Mesh*                      mesh               = (Mesh*)        velocityField->feMesh;
	Grid*			   vertGrid;
	Node_LocalIndex            nodeGlobalCount;
	Node_LocalIndex            nodeLocalCount;
	Dof_Index                  dof                = velocityField->fieldComponentCount;
	double                     min[3];
	double                     max[3];
	Node_LocalIndex            localNode_I;
	Node_LocalIndex            globalNode_I;
	double*                    coord;
	IJK                        globalNodeIJK      = { 0, 0, 0 };
	Index                      array_I;
	Dof_Index                  dof_I;
	XYZ                        movingMeshVelocity = { 0.0, 0.0, 0.0 };
	Variable*                  currVariable;
	double*                    nodeDataPtr;

	nodeGlobalCount = Mesh_GetGlobalSize( mesh, MT_VERTEX );
	nodeLocalCount = Mesh_GetLocalSize( mesh, MT_VERTEX );
	Mesh_GetGlobalCoordRange( mesh, min, max );
	vertGrid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
						  ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

	for ( globalNode_I = 0 ; globalNode_I < nodeGlobalCount ; globalNode_I++ ) {
		/* Check if node is on processor */
		if ( !Mesh_GlobalToDomain( mesh, MT_VERTEX, globalNode_I, &localNode_I ) || 
		     localNode_I >= nodeLocalCount )
		{
			continue;
		}

		coord = Mesh_GetVertex( mesh, localNode_I );
		Grid_Lift( vertGrid, globalNode_I, globalNodeIJK );

		array_I = globalNodeIJK[ J_AXIS ] * dof; /* TODO GET TO WORK IN 3D */
		movingMeshVelocity[ I_AXIS ] = 
			(coord[ I_AXIS ] - min[ I_AXIS ])/(max[ I_AXIS ] - min[ I_AXIS ]) * 
			( xRightVelocities[ array_I ] - xLeftVelocities[ array_I ] ) + xLeftVelocities[ array_I ];

		for ( dof_I = 0 ; dof_I < dof ; dof_I++ ) {
			currVariable = DofLayout_GetVariable( velocityField->dofLayout, localNode_I, dof_I );
			nodeDataPtr = Variable_GetPtrDouble( currVariable, localNode_I );

			if ( flag == ADD )
				*nodeDataPtr += movingMeshVelocity[ dof_I ];
			else 
				*nodeDataPtr -= movingMeshVelocity[ dof_I ];
		}
	}
}

void MovingMeshEnergyCorrection( void* sle, void* data ) {
	UnderworldContext*                                      context                 = (UnderworldContext*) data;
	Underworld_MovingMeshEnergyCorrection_ContextExt*       contextExt;
	FeVariable*                                             velocityField           = context->velocityField;
	double*                                                 xLeftVelocities;
	double*                                                 xRightVelocities;
	
	
	contextExt = ExtensionManager_Get( 
		context->extensionMgr, 
		context, 
		Underworld_MovingMeshEnergyCorrection_ContextExtHandle );
	
	/* Get Side Wall velocities */
	FeVariable_GetSideWallValues( velocityField, I_AXIS, &xLeftVelocities, &xRightVelocities );

	/* Correct velocity accordingly */
	MovingMeshEnergyCorrection_AddCorrection( velocityField, xLeftVelocities, xRightVelocities, MINUS );

	/* Solve Energy equation */
	contextExt->energySolverExecute( sle, context );

	/* Reverse correction */
	MovingMeshEnergyCorrection_AddCorrection( velocityField, xLeftVelocities, xRightVelocities, ADD );

	Memory_Free( xLeftVelocities );
	Memory_Free( xRightVelocities );
}

void _Underworld_MovingMeshEnergyCorrection_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	UnderworldContext*                                      context = 
	Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data ); 
	Underworld_MovingMeshEnergyCorrection_ContextExt*       contextExt;
	
	Journal_DFirewall( 
		(Bool)context, 
		Journal_Register( Error_Type, Underworld_MovingMeshEnergyCorrection_Type ), 
		"No context found\n" );
	Journal_DFirewall( 
		(Bool)context->energySLE, 
		Journal_Register( Error_Type, Underworld_MovingMeshEnergyCorrection_Type ), 
		"The required energy SLE component has not been created or placed on the context.\n");	
	
	/* Add the extension to the context */
	Underworld_MovingMeshEnergyCorrection_ContextExtHandle = ExtensionManager_Add(
		context->extensionMgr, 
		Underworld_MovingMeshEnergyCorrection_Type, 
		sizeof( Underworld_MovingMeshEnergyCorrection_ContextExt ) );
	contextExt = ExtensionManager_Get( 
		context->extensionMgr, 
		context, 
		Underworld_MovingMeshEnergyCorrection_ContextExtHandle );
	
	/* Replace the energy SLE's execute with this one. Save the old value for use later. */
	contextExt->energySolverExecute = context->energySLE->_execute;
	context->energySLE->_execute = MovingMeshEnergyCorrection;
}

/* This function will provide StGermain the abilty to instantiate (create) this codelet on demand. */
void* _Underworld_MovingMeshEnergyCorrection_DefaultNew( Name name ) {
	return Codelet_New(
			Underworld_MovingMeshEnergyCorrection_Type,
			_Underworld_MovingMeshEnergyCorrection_DefaultNew,
			_Underworld_MovingMeshEnergyCorrection_Construct, /* SQ NOTE: Used to be a construct extensions. */
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}

/* This function is automatically run by StGermain when this plugin is loaded. The name must be "<plugin-name>_Register". */
Index Underworld_MovingMeshEnergyCorrection_Register( PluginsManager* pluginsManager ) {
	/* A plugin is only properly registered once it returns the handle provided when submitting a codelet to StGermain. */
	return PluginsManager_Submit( 
		pluginsManager, 
		Underworld_MovingMeshEnergyCorrection_Type, 
		"0", 
		_Underworld_MovingMeshEnergyCorrection_DefaultNew );
}
