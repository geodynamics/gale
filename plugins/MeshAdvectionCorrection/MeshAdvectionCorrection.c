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
** $Id: MeshAdvectionCorrection.c 610 2007-10-11 08:09:29Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include <string.h>

/* Each component type needs a unique identifier (as a string for us to read and as an integer for realtime comparisions) */
const Type Underworld_MeshAdvectionCorrection_Type = "Underworld_MeshAdvectionCorrection_Type";

ExtensionInfo_Index Underworld_MeshAdvectionCorrection_ContextExtHandle;

typedef struct {
	/* Save the "old" execute function pointer (which we replace) */ 
	Stg_Component_ExecuteFunction* energySolverExecute;
	void (*_findArtificalVelocity)(void* feVariable, double* artVelocity, void* data);
	Axis extAxis;
} Underworld_MeshAdvectionCorrection_ContextExt;


void MeshAdVectionCorrection_SingleAxisExtension( void* feVariable, double* artVelocity, void* data ) {
	/* This function assumes the two plane walls in a single axis only are moving. */
	FeVariable*                self          = (FeVariable*)  feVariable;
	Mesh*                      mesh          = (Mesh*)        self->feMesh;
	int                        planeAxis     = *(int*)data; /* recast data */
	Grid*                      vertGrid;

	double*                    valueMinSideWallLocal;
	double*                    valueMaxSideWallLocal;
	double*                    valueMinSideWallGlobal;
	double*                    valueMaxSideWallGlobal;
	double                     minCoord[3], maxCoord[3];
	double*                    coord;
	int                        globalNodeIJK[3];
 	double                     artificialVelocity[3] = {0,0,0};
	Dimension_Index            aAxis, bAxis;
	IJK                        minNodeIJK, maxNodeIJK;
	Index                      array_I;
	Node_Index                 aNode_I, bNode_I;
	Node_GlobalIndex           minNodeGlobal_I, maxNodeGlobal_I;
	Node_LocalIndex            minNodeLocal_I, maxNodeLocal_I;
	Node_Index                 sideNodeCount, lNodeCount, lNode_I, gNode_I;
	//MPI_Comm                   comm;
	Dof_Index                  dof           = self->fieldComponentCount;
	Dof_Index                  dof_I;

	//comm = Mesh_GetCommTopology( mesh, MT_VERTEX ) ;
	vertGrid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
						  ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );
	
	lNodeCount = FeMesh_GetNodeLocalSize( mesh );

	aAxis =  ( planeAxis == I_AXIS ? J_AXIS : I_AXIS );
	sideNodeCount = vertGrid->sizes[ aAxis ];

	if( self->dim == 3 ) {
		bAxis = ( planeAxis == K_AXIS ? I_AXIS : K_AXIS );
		sideNodeCount = vertGrid->sizes[ aAxis ] * vertGrid->sizes[ bAxis ];
	}

	valueMinSideWallLocal = Memory_Alloc_Array( double, sideNodeCount*dof, "valueMinSideWallLocal" );
	valueMaxSideWallLocal = Memory_Alloc_Array( double, sideNodeCount*dof, "valueMaxSideWallLocal" );
	valueMinSideWallGlobal = Memory_Alloc_Array( double, sideNodeCount*dof, "valueMinSideWallGlobal" );
	valueMaxSideWallGlobal = Memory_Alloc_Array( double, sideNodeCount*dof, "valueMaxSideWallGlobal" );

	/* Loop over all nodes on the minimum side wall and the maximum side wall */
	minNodeIJK[ planeAxis ] = 0;
	maxNodeIJK[ planeAxis ] = vertGrid->sizes[ planeAxis ] - 1 ;
	array_I = 0;

	if( self->dim == 3 ) {
		for ( aNode_I = 0 ; aNode_I < vertGrid->sizes[ aAxis ] ; aNode_I++ ) {
			minNodeIJK[ aAxis ] = maxNodeIJK[ aAxis ] = aNode_I;
			for ( bNode_I = 0 ; bNode_I < vertGrid->sizes[ bAxis ] ; bNode_I++ ) {
				minNodeIJK[ bAxis ] = maxNodeIJK[ bAxis ] = bNode_I;

				/* Get Global Node Indicies for the node at these walls */
				minNodeGlobal_I = Grid_Project( vertGrid, minNodeIJK );
				maxNodeGlobal_I = Grid_Project( vertGrid, maxNodeIJK );

				/* Check to see whether this min node is on this processor */
				if ( Mesh_GlobalToDomain( mesh, MT_VERTEX, minNodeGlobal_I, &minNodeLocal_I ) )
					FeVariable_GetValueAtNode( feVariable, minNodeLocal_I, &valueMinSideWallLocal[ array_I*dof ] );
				else {
					for ( dof_I = 0 ; dof_I < dof ; dof_I++ ) 
						/* Initialise all values to some really large value so that 
						 * MPI_Allreduce can find value from correct processor */
						valueMinSideWallLocal[ array_I*dof + dof_I ] = HUGE_VAL;
				}
				
				/* Check to see whether this max node is on this processor */
				if ( Mesh_GlobalToDomain( mesh, MT_VERTEX, maxNodeGlobal_I, &maxNodeLocal_I ) )
					FeVariable_GetValueAtNode( feVariable, maxNodeLocal_I, &valueMaxSideWallLocal[ array_I*dof ] );
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
	} else {
		for( aNode_I = 0 ; aNode_I < vertGrid->sizes[ aAxis ] ; aNode_I++ ) {
			minNodeIJK[ aAxis ] = maxNodeIJK[ aAxis ] = aNode_I;
			/* Get Global Node Indicies for the node at these walls */
			minNodeGlobal_I = Grid_Project( vertGrid, minNodeIJK );
			maxNodeGlobal_I = Grid_Project( vertGrid, maxNodeIJK );

			/* Check to see whether this min node is on this processor */
			if ( Mesh_GlobalToDomain( mesh, MT_VERTEX, minNodeGlobal_I, &minNodeLocal_I ) )
				FeVariable_GetValueAtNode( feVariable, minNodeLocal_I, &valueMinSideWallLocal[ array_I*dof ] );
			else {
				for ( dof_I = 0 ; dof_I < dof ; dof_I++ ) 
					/* Initialise all values to some really large value so that 
					 * MPI_Allreduce can find value from correct processor */
					valueMinSideWallLocal[ array_I*dof + dof_I ] = HUGE_VAL;
			}
			
			/* Check to see whether this max node is on this processor */
			if ( Mesh_GlobalToDomain( mesh, MT_VERTEX, maxNodeGlobal_I, &maxNodeLocal_I ) )
				FeVariable_GetValueAtNode( feVariable, maxNodeLocal_I, &valueMaxSideWallLocal[ array_I*dof ] );
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


	MPI_Allreduce( valueMinSideWallLocal, valueMinSideWallGlobal, sideNodeCount*dof, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
	MPI_Allreduce( valueMaxSideWallLocal, valueMaxSideWallGlobal, sideNodeCount*dof, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );

	Memory_Free( valueMinSideWallLocal );
	Memory_Free( valueMaxSideWallLocal );
	
	/** Set artificial velcoity field **/

	Mesh_GetGlobalCoordRange( mesh, minCoord, maxCoord );

	for( lNode_I = 0 ; lNode_I < lNodeCount ; lNode_I++ ) {
		/* get node coord */
		coord = Mesh_GetVertex( mesh, lNode_I );

		/* get IJK params for node */
		gNode_I = Mesh_DomainToGlobal( mesh, MT_VERTEX, lNode_I );
		Grid_Lift( vertGrid, gNode_I, globalNodeIJK );

		/* assumes 2D only, and extension in X */
		array_I = globalNodeIJK[aAxis] * dof + planeAxis;
		/* 3D */
		if( self->dim==3 ) 
			array_I = (globalNodeIJK[aAxis]*vertGrid->sizes[bAxis] + globalNodeIJK[bAxis]) * dof + planeAxis;

		artificialVelocity[ planeAxis ] = 
			(coord[ planeAxis ] - minCoord[ planeAxis ])/(maxCoord[ planeAxis ] - minCoord[ planeAxis ]) * 
			( valueMaxSideWallGlobal[ array_I ] - valueMinSideWallGlobal[ array_I ] ) + valueMinSideWallGlobal[ array_I ];

		memcpy( &artVelocity[lNode_I*dof] , artificialVelocity, dof*sizeof(double) );
	}
}

typedef enum {
	ADD,
	MINUS
} CorrectionFlag;

void MeshAdvectionCorrection_AddCorrection( FeVariable* feVariable, double* artVelocity, CorrectionFlag flag ) {
	/* This function offset the VelocityField by the artificialVelocityField */
	FeVariable*  self          = (FeVariable*)  feVariable;
	FeMesh*      mesh          = self->feMesh;
	int          lNodeCount    = FeMesh_GetNodeLocalSize( mesh );
	int          dof           = self->fieldComponentCount;
	int          dim           = self->dim;
	double       oldVec[3], currentVec[3], newVec[3];
	int lNode_I;

	for( lNode_I = 0 ; lNode_I < lNodeCount ; lNode_I++ ) {
		memcpy( oldVec, &artVelocity[lNode_I*dof], dof*sizeof(double) );

		FeVariable_GetValueAtNode( feVariable, lNode_I, currentVec );

		if ( flag == ADD ) 
			StGermain_VectorAddition(newVec, currentVec, oldVec, dim );
		else 
			StGermain_VectorSubtraction(newVec, currentVec, oldVec, dim );

		FeVariable_SetValueAtNode( feVariable, lNode_I, newVec );
	}
}

void MeshAdvectionCorrection( void* sle, void* data ) {
	UnderworldContext*                                      context                 = (UnderworldContext*) data;
	Underworld_MeshAdvectionCorrection_ContextExt*          plugin;
	FeVariable*                                             feVariable           = context->velocityField;
	double*                                                 artVelocity;
	int lNodeCount;
	
	lNodeCount = FeMesh_GetNodeLocalSize( feVariable->feMesh );

	artVelocity = Memory_Alloc_Array( double, 
			lNodeCount * feVariable->fieldComponentCount, 
			"artificial nodal velocities" );

	plugin = ExtensionManager_Get( 
		context->extensionMgr, 
		context, 
		Underworld_MeshAdvectionCorrection_ContextExtHandle );

	/* Fine articfical nodal velocities */
	/* IF plugin->extAxis is valid then */
	plugin->_findArtificalVelocity( feVariable, artVelocity, (void*)&plugin->extAxis );
	/* ELSE
	 * 	yet to be implemented */

	/* Correct velocity and re-sync shadow space */
	MeshAdvectionCorrection_AddCorrection( feVariable, artVelocity, MINUS );
	FeVariable_SyncShadowValues( feVariable );

	/* Solve Energy equation */
	plugin->energySolverExecute( sle, context );

	/* Reverse correction and re-sync */
	MeshAdvectionCorrection_AddCorrection( feVariable, artVelocity, ADD );
	FeVariable_SyncShadowValues( feVariable );

	Memory_Free( artVelocity );
}

void _Underworld_MeshAdvectionCorrection_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	UnderworldContext*                                      context = 
	Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data ); 
	Underworld_MeshAdvectionCorrection_ContextExt*       plugin;
	char* extAxisName = NULL;
	
	Journal_DFirewall( 
		(Bool)context, 
		Journal_Register( Error_Type, Underworld_MeshAdvectionCorrection_Type ), 
		"No context found\n" );
	Journal_DFirewall( 
		(Bool)context->energySLE, 
		Journal_Register( Error_Type, Underworld_MeshAdvectionCorrection_Type ), 
		"The required energy SLE component has not been created or placed on the context.\n");	
	
	/* Add the extension to the context */
	Underworld_MeshAdvectionCorrection_ContextExtHandle = ExtensionManager_Add(
		context->extensionMgr, 
		Underworld_MeshAdvectionCorrection_Type, 
		sizeof( Underworld_MeshAdvectionCorrection_ContextExt ) );
	plugin = ExtensionManager_Get( 
		context->extensionMgr, 
		context, 
		Underworld_MeshAdvectionCorrection_ContextExtHandle );

	extAxisName = Stg_ComponentFactory_GetRootDictString( cf, "MeshAdvectionCorrect_wallExtensionAxis", "x" );
	switch ( extAxisName[0] ) {
		case 'x': case 'X': case 'i': case 'I': case '0':
			plugin->extAxis = I_AXIS; break;
		case 'y': case 'Y': case 'j': case 'J': case '1':
			plugin->extAxis = J_AXIS; break;
		case 'z': case 'Z': case 'k': case 'K': case '2':
			plugin->extAxis = K_AXIS; break;
		default:
			Journal_Firewall( False, Journal_Register( Error_Type, Underworld_MeshAdvectionCorrection_Type ),
				"Error: MeshAdvectionCorrect doesn't know how to estimate wallVC\n"
				"For example add:\n"
				"<param name=\"MeshAdvectionCorrect_wallExtensionAxis\">x</param>\n");
	}
	
	plugin->_findArtificalVelocity = MeshAdVectionCorrection_SingleAxisExtension;
	/* Replace the energy SLE's execute with this one. Save the old value for use later. */
	plugin->energySolverExecute = context->energySLE->_execute;
	context->energySLE->_execute = MeshAdvectionCorrection;
}

/* This function will provide StGermain the abilty to instantiate (create) this codelet on demand. */
void* _Underworld_MeshAdvectionCorrection_DefaultNew( Name name ) {
	return Codelet_New(
			Underworld_MeshAdvectionCorrection_Type,
			_Underworld_MeshAdvectionCorrection_DefaultNew,
			_Underworld_MeshAdvectionCorrection_Construct, /* SQ NOTE: Used to be a construct extensions. */
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}

/* This function is automatically run by StGermain when this plugin is loaded. The name must be "<plugin-name>_Register". */
Index Underworld_MeshAdvectionCorrection_Register( PluginsManager* pluginsManager ) {
	/* A plugin is only properly registered once it returns the handle provided when submitting a codelet to StGermain. */
	return PluginsManager_Submit( 
		pluginsManager, 
		Underworld_MeshAdvectionCorrection_Type, 
		"0", 
		_Underworld_MeshAdvectionCorrection_DefaultNew );
}
