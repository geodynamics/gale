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
** $Id: MovingMesh.c 734 2008-05-15 23:15:29Z LukeHodkinson $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include "MovingMesh.h"
#include <assert.h>


/* Each component type needs a unique identifier (as a string for us to read and as an integer for realtime comparisions) */
const Type Underworld_MovingMesh_Type = "Underworld_MovingMesh";

void _Underworld_MovingMesh_Construct( void* meshExtender, Stg_ComponentFactory* cf, void* data ) {
	UnderworldContext*  context = Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data );
	MeshExtender*       self = (MeshExtender*) meshExtender;
	Dimension_Index     dim_I = 0;
	Dimension_Index     axisToRemeshOnTotal = 0;
	
	Journal_Firewall( 
		(Bool)context, 
		Journal_Register( Error_Type, Underworld_MovingMesh_Type ), 
		"No context found\n" );

	self->context = (AbstractContext*)context;	
	self->velocityField = context->velocityField;

	Journal_Firewall( (Bool)self->velocityField, 
			  Journal_Register( Error_Type, Underworld_MovingMesh_Type ), 
			  "The required velocity field component has not been created or placed on the context.\n");

	self->remeshAccordingToAxis[I_AXIS] = 		
		Dictionary_GetBool_WithDefault( context->dictionary, "remeshAccordingToIAxis", True );
	self->remeshAccordingToAxis[J_AXIS] = 		
		Dictionary_GetBool_WithDefault( context->dictionary, "remeshAccordingToJAxis", True );
	if ( self->velocityField->dim == 2 ) {
		self->remeshAccordingToAxis[K_AXIS] = False;
	}
	else {
		self->remeshAccordingToAxis[K_AXIS] = 		
			Dictionary_GetBool_WithDefault( context->dictionary, "remeshAccordingToKAxis", True );
	}

	axisToRemeshOnTotal = 0;
	for ( dim_I = 0; dim_I < self->velocityField->dim; dim_I++ ) {
		if ( self->remeshAccordingToAxis[dim_I] == True ) {
			axisToRemeshOnTotal++;
		}
	}

	Journal_Firewall( 
		axisToRemeshOnTotal > 0,
		Journal_Register( Error_Type, Underworld_MovingMesh_Type ), 
		"Error: in %s: you have disabled remeshing on all axis. Please set at least one axis to remesh on "
		"using the remeshAccordingToIAxis, remeshAccordingToJAxis, remeshAccordingToKAxis dictionary "
		"parameters, if you wish to use the remesher.\n",
		__func__ );

	/*  */ 
	/* TODO: this timeIntegrator interface seems weird. Would have thought you'd pass a ptr to self. */
	TimeIntegrator_PrependFinishEP( 
			context->timeIntegrator, "Underworld_MovingMesh_Remesh", Underworld_MovingMesh_Remesh, 
			CURR_MODULE_NAME, self );

	/* Add to Live Stg_Component Register so it gets automatically built, initialised and deleted */
	LiveComponentRegister_Add( context->CF->LCRegister, (Stg_Component*) self );
}

void Underworld_MovingMesh_Build( void* meshExtender, void* data ) {
	MeshExtender*       self = (MeshExtender*)meshExtender;
	Mesh*		    mesh;

	mesh = (Mesh*)self->velocityField->feMesh;

	Journal_Firewall( ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) != (unsigned)-1, 
			  Journal_Register( Error_Type, Underworld_MovingMesh_Type ), 
			  "Error: in %s: provided Velocity field's mesh doesn't have a regular decomposition.\n",
			  __func__ );
}


void Underworld_MovingMesh_Remesh( TimeIntegrator* timeIntegrator, MeshExtender* self ) {
	FeVariable*           velocityField = self->velocityField;
	Mesh*                 mesh          = (Mesh*) velocityField->feMesh;
	unsigned	      rank;
	Stream*               debug = Journal_Register( Debug_Type, Underworld_MovingMesh_Type );
	Stream*               info = Journal_Register( Info_Type, self->type );
	double                remeshTime, remeshTimeStart, remeshTimeEnd;
	Dimension_Index       dim_I = 0;
	Dimension_Index       axisToRemeshOnTotal = 0;
	double                dt = 0.0;
	MPI_Comm	      comm;
	#if DEBUG
	Index                 element_dI;
	#endif

	/*comm = CommTopology_GetComm( Mesh_GetCommTopology( mesh, MT_VERTEX ) );*/
	comm = Comm_GetMPIComm( Mesh_GetCommTopology( mesh, MT_VERTEX ) );
	MPI_Comm_rank( comm, (int*)&rank );

	Journal_Printf( info, "%d: starting remeshing of mesh \"%s\":\n", rank, mesh->name );
	axisToRemeshOnTotal = 0;
	Journal_Printf( info, "Axis set to remesh on are : " );
	for ( dim_I = 0; dim_I < velocityField->dim; dim_I++ ) {
		if ( self->remeshAccordingToAxis[dim_I] == True ) {
                  /* Journal_Printf( info, "%c, ", IJKTopology_DimNumToDimLetter[dim_I] );*/
			axisToRemeshOnTotal++;
		}
	}
	Journal_Printf( info, "\n" );
	Journal_DPrintf( debug, "In %s(): about to remesh:\n", __func__ );
	Stream_Indent( debug );

	Journal_Firewall( 
		axisToRemeshOnTotal > 0,
		Journal_Register( Error_Type, Underworld_MovingMesh_Type ), 
		"Error: in %s: you have disabled remeshing on all axis. Please set at least one axis to remesh on "
		"using the remeshAccordingToIAxis, remeshAccordingToJAxis, remeshAccordingToKAxis dictionary "
		"parameters, if you wish to use the remesher.\n",
		__func__ );

	dt = AbstractContext_Dt( self->context );
	/* TODO: the if statement is here since we seem to be using last timestep's dt to update. Is this correct?
	 PatrickSunter - 5 June 2006 */
	if ( self->context->timeStep > 1 ) {
		Journal_Firewall( 
			dt > 0.0,
			Journal_Register( Error_Type, Underworld_MovingMesh_Type ), 
			"Error- in %s: in timeStep %u, provided Context \"%s\"'s dt <= 0.\n",
			__func__, self->context->timeStep, self->context->name );
	}	
	
	remeshTimeStart = MPI_Wtime();

	Journal_Firewall( ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) != (unsigned)-1, 
			  Journal_Register( Error_Type, Underworld_MovingMesh_Type ), 
			  "Error: in %s: provided Velocity field's mesh doesn't have a regular decomposition.\n",
			  __func__ );

	/*
	Journal_Firewall( 
		True == Stg_Class_IsInstance( velocityField->feMesh->layout->elementLayout, ParallelPipedHexaEL_Type ),
		Journal_Register( Error_Type, Underworld_MovingMesh_Type ), 
		"Error: in %s: provided Velocity field's mesh doesn't have a %s elementLayout.\n",
		__func__, ParallelPipedHexaEL_Type );
	*/

	Underworld_MovingMesh_RemeshAccordingToSidewalls( self, velocityField );
	
	Journal_DPrintf( debug, "remeshing complete.\n" );

	/*
	#if DEBUG
	if ( Stream_IsEnable( debug ) && Stream_IsPrintableLevel( debug, 3 ) ) {
		Journal_DPrintfL( debug, 3, "Updated node coords:\n" );
		Stream_Indent( debug );
		for ( element_dI = 0; element_dI < Mesh_GetDomainSize( velocityField->feMesh, MT_VERTEX ); element_dI++ ) {
			Journal_DPrintfL( debug, 3, "In Element %3d: ", element_dI );
			Mesh_PrintNodeCoordsOfElement( velocityField->feMesh, element_dI, debug );
		}	
		Stream_UnIndent( debug );
	}	
	#endif
	*/

	remeshTimeEnd = MPI_Wtime();
	remeshTime = remeshTimeEnd - remeshTimeStart;
	Journal_Printf( info, "%d: finished remeshing: took %f sec.\n", rank, remeshTime );

	Stream_UnIndent( debug );

	return;
}


void Underworld_MovingMesh_RemeshAccordingToSidewalls( MeshExtender* self, FeVariable* velocityField ) {
	Mesh*                 mesh          = (Mesh*) velocityField->feMesh;
	Dimension_Index       remeshAxis_I;
	
	for ( remeshAxis_I = 0 ; remeshAxis_I < velocityField->dim ; remeshAxis_I++ ) {
		if ( self->remeshAccordingToAxis[remeshAxis_I] == True ) {
			Underworld_MovingMesh_RemeshAccordingToSidewall_SingleAxis( self, velocityField, remeshAxis_I );
		}	
	}

	/* Update mesh details. */
	Mesh_DeformationUpdate( mesh );
}


void Underworld_MovingMesh_RemeshAccordingToSidewall_SingleAxis(
		MeshExtender*    self,
		FeVariable*      velocityField,
		Dimension_Index  remeshAxis )
{
	Mesh*                      mesh          = (Mesh*) velocityField->feMesh;
	double			   maxCrd[3], minCrd[3];
	Grid*			   vertGrid;
	/* Assume IJK Topology */
	double*                    minGlobal;
	double*                    maxGlobal;
	Node_Index                 sideWallNodeCount;
	Dimension_Index            otherAxisA = (remeshAxis + 1) % 3;
	Dimension_Index            otherAxisB = (remeshAxis + 2) % 3;
	IJK                        nodeIJK;
	Node_Index                 sideWallNode_I;
	Node_Index                 aNode_I;
	Node_Index                 bNode_I;
	Node_Index                 remeshAxisNode_I;
	Node_GlobalIndex           nodeGlobal_I;
	Node_DomainIndex           nodeDomain_I;
	double                     newElementWidthInRemeshAxis;
	double                     minGlobalCoord =  HUGE_VAL;
	double                     maxGlobalCoord = -HUGE_VAL;
	/*Bool                       nodeG2DBuiltTemporarily = False;*/
	double                     newCoordInRemeshAxis = 0;
	double                     tolerance;

	Stream* debug = Journal_Register( Debug_Type, Underworld_MovingMesh_Type );

	/*Journal_DPrintf( debug, "In %s(): for remeshAxis %c\n", __func__, IJKTopology_DimNumToDimLetter[remeshAxis] );*/
	Stream_Indent( debug );	

	Mesh_GetGlobalCoordRange( mesh, minCrd, maxCrd );
	tolerance = (maxCrd[remeshAxis] - minCrd[remeshAxis]) * 1e-9;

	vertGrid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
						  ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

	sideWallNodeCount = vertGrid->sizes[ otherAxisA ] * vertGrid->sizes[ otherAxisB ];
	minGlobal = Memory_Alloc_Array( double, sideWallNodeCount, "min node coords" );
	maxGlobal = Memory_Alloc_Array( double, sideWallNodeCount, "max node coords" );

	/* As long as memory not too tight, create temporary tables to massively speed up
	   necessary global to local calculations. */
	/*
	if ( (mesh->buildTemporaryGlobalTables == True) && (mesh->nodeG2D == 0) ) {
		mesh->nodeG2D = MeshDecomp_BuildNodeGlobalToDomainMap( meshDecomp );
		nodeG2DBuiltTemporarily = True;
	}
	*/
	
	/* calc min and max wall new coords given velocity field */
	Underworld_MovingMesh_CalculateMinOrMaxCoordsOnSidewall(
		self, velocityField, remeshAxis, MIN_COORDS, minGlobal );
	Underworld_MovingMesh_CalculateMinOrMaxCoordsOnSidewall(
		self, velocityField, remeshAxis, MAX_COORDS, maxGlobal );

	/* calc overall min and max wall new coords */
	for ( sideWallNode_I = 0; sideWallNode_I < sideWallNodeCount; sideWallNode_I++ ) {
		if ( minGlobal[ sideWallNode_I ] < minGlobalCoord )
			minGlobalCoord = minGlobal[ sideWallNode_I ];
		if ( maxGlobal[ sideWallNode_I ] > maxGlobalCoord )
			maxGlobalCoord = maxGlobal[ sideWallNode_I ];
	}	

	if ( ( fabs( minGlobalCoord - minCrd[ remeshAxis ] ) < tolerance ) &&
		( fabs( maxGlobalCoord - maxCrd[ remeshAxis ] ) < tolerance ) )
	{
          /*Journal_DPrintfL( debug, 1, "Found no extension/compression in remeshAxis %c: thus not "
            "updating node coords.\n", IJKTopology_DimNumToDimLetter[remeshAxis] );*/
		Memory_Free( minGlobal );
		Memory_Free( maxGlobal );
		/*
		if ( nodeG2DBuiltTemporarily ) {
			Memory_Free( mesh->nodeG2D );
			mesh->nodeG2D = NULL;
		}
		*/
		Stream_UnIndent( debug );	
		return;
	}

	/* Change Geometry */
	/*
	Journal_DPrintfL( debug, 3, "Updating mesh's BlockGeometry \"%s\" with new global "
		"min and max in remesh axis %c: min=%.3f, max=%.3f\n",
		geometry->name, IJKTopology_DimNumToDimLetter[remeshAxis], minGlobalCoord, maxGlobalCoord );

	minCrd[ remeshAxis ] = minGlobalCoord;
	maxCrd[ remeshAxis ] = maxGlobalCoord;
	*/

	/*Journal_DPrintfL( debug, 2, "Looping over all nodes, updating coordinates in remeshAxis %c based on "
          "stretch/compression on the relevant side walls.\n", IJKTopology_DimNumToDimLetter[remeshAxis] );*/
	Stream_Indent( debug );	
	sideWallNode_I = 0;
	for ( aNode_I = 0 ; aNode_I < vertGrid->sizes[ otherAxisA ] ; aNode_I++ ) {
		nodeIJK[ otherAxisA ] = aNode_I;
		for ( bNode_I = 0 ; bNode_I < vertGrid->sizes[ otherAxisB ] ; bNode_I++ ) {
			nodeIJK[ otherAxisB ] = bNode_I;

			newElementWidthInRemeshAxis = (maxGlobal[ sideWallNode_I ] - minGlobal[ sideWallNode_I ])
				/ ( (double)vertGrid->sizes[ remeshAxis ] - 1.0 );

			/*Journal_DPrintfL( debug, 3, "For slice with aNode_I=%u in otherAxisA=%c and bNode_I=%u in "
				"otherAxisB=%c: calculated new element width in remesh axis= %.3f\n",
				aNode_I, IJKTopology_DimNumToDimLetter[otherAxisA],
				bNode_I, IJKTopology_DimNumToDimLetter[otherAxisB],
				newElementWidthInRemeshAxis );*/
			Stream_Indent( debug );	

			for ( remeshAxisNode_I = 0 ; remeshAxisNode_I < vertGrid->sizes[ remeshAxis ] ; remeshAxisNode_I++ ) {
				nodeIJK[ remeshAxis ] = remeshAxisNode_I;

				/* Find this local coord */
				nodeGlobal_I = Grid_Project( vertGrid, nodeIJK );

				/* Adjust position if this node is on my processor */
				if ( Mesh_GlobalToDomain( mesh, MT_VERTEX, nodeGlobal_I, &nodeDomain_I ) && 
				     nodeDomain_I < Mesh_GetDomainSize( mesh, MT_VERTEX ) )
				{
					newCoordInRemeshAxis = minGlobal[ sideWallNode_I ] +
						newElementWidthInRemeshAxis * (double) remeshAxisNode_I;
					Journal_DPrintfL( debug, 3, "updating node coord[%u] at IJK(%u,%u,%u) (domain "
						"node I %u) from %.3f to %.3f\n", remeshAxis, nodeIJK[0], nodeIJK[1], nodeIJK[2],
						nodeDomain_I, Mesh_GetVertex( mesh, nodeDomain_I )[ remeshAxis ],
						newCoordInRemeshAxis );
					Mesh_GetVertex( mesh, nodeDomain_I )[ remeshAxis ] = newCoordInRemeshAxis;
				}		
			}
			Stream_UnIndent( debug );	

			/* Increment counter for side wall node index */
			sideWallNode_I++;
		}
	}
	Stream_UnIndent( debug );	

	/*
	if ( nodeG2DBuiltTemporarily ) {
		Memory_Free( mesh->nodeG2D );
		mesh->nodeG2D = NULL;
	}
	*/

	/* Update the element layout's partition info based on geometry changes... */
	/*
	ParallelPipedHexaEL_UpdateGeometryPartitionInfo( mesh->layout->elementLayout, mesh->layout->decomp );
	*/

	/* Update mesh information. */
	Mesh_Sync( mesh );
	Mesh_DeformationUpdate( mesh );

	Memory_Free( minGlobal );
	Memory_Free( maxGlobal );
	Stream_UnIndent( debug );	
}


void Underworld_MovingMesh_CalculateMinOrMaxCoordsOnSidewall(
		MeshExtender*    self,
		FeVariable*      velocityField,
		Dimension_Index  remeshAxis,
		MinOrMaxFlag     minOrMaxFlag,
		double*          newWallCoordsInRemeshAxisGlobal )
{
	Mesh*                      mesh          = (Mesh*) velocityField->feMesh;
	Grid*			   vertGrid;
	double			   maxCrd[3], minCrd[3];
	IJK                        wallNodeIJK;
	Node_Index                 sideWallNode_I;
	Node_LocalIndex            wallNode_lI;
	Node_Index                 aNode_I;
	Node_Index                 bNode_I;
	Node_Index                 sideWallNodeCount;
	Dimension_Index            otherAxisA = (remeshAxis + 1) % 3;
	Dimension_Index            otherAxisB = (remeshAxis + 2) % 3;
	MPI_Op                     mpiOperation;
	double*                    newWallCoordsInRemeshAxis = NULL;
	Node_LocalIndex            lastCalculatedWallNode_lI = 0;
	Node_Index                 lastCalculatedWallNode_sideI = 0;
	double                     lastCalculatedNewWallCoordInRemeshAxis;
	double*                    currNodeCoord;
	double                     velVector[3];
	/* Somewhat arbitrary tolerance to check that extension of the mesh is not irregular in the current axis */
	double                     tolerance;
	double                     dt = AbstractContext_Dt( self->context );
	Bool                       atLeastOneLocalSideNodeProcessed = False;
	MPI_Comm		   comm;

	Mesh_GetGlobalCoordRange( mesh, minCrd, maxCrd );
	tolerance = (maxCrd[remeshAxis] - minCrd[remeshAxis]) * 1e-9;

	vertGrid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
						  ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

	sideWallNodeCount = vertGrid->sizes[ otherAxisA ] * vertGrid->sizes[ otherAxisB ];
	newWallCoordsInRemeshAxis = Memory_Alloc_Array( double, sideWallNodeCount, "newWallCoordsInRemeshAxis" );

	if ( minOrMaxFlag == MIN_COORDS ) {
		/* Loop over all nodes on the minimum side wall */
		wallNodeIJK[ remeshAxis ] = 0;
	}
	else {
		/* Loop over all nodes on the maximum side wall */
		wallNodeIJK[ remeshAxis ] = vertGrid->sizes[ remeshAxis ] - 1 ;
	}

	sideWallNode_I = 0;
	for ( aNode_I = 0 ; aNode_I < vertGrid->sizes[ otherAxisA ] ; aNode_I++ ) {
		wallNodeIJK[ otherAxisA ] = aNode_I;
		for ( bNode_I = 0 ; bNode_I < vertGrid->sizes[ otherAxisB ] ; bNode_I++ ) {
			wallNodeIJK[ otherAxisB ] = bNode_I;

			/* Convert to Local Index */
			wallNode_lI = RegularMeshUtils_Node_3DTo1D( mesh, wallNodeIJK );

			/* Check to see whether this wall node is on this processor */
			if ( Mesh_GlobalToDomain( mesh, MT_VERTEX, wallNode_lI, &wallNode_lI ) && 
			     wallNode_lI < Mesh_GetLocalSize( mesh, MT_VERTEX ) )
			{
				/* Record the current coord at this node */
				currNodeCoord = Mesh_GetVertex( mesh, wallNode_lI );
				newWallCoordsInRemeshAxis[sideWallNode_I] = currNodeCoord[remeshAxis];
				/* Now add the velocity BC / solved velocity on this wall to see where it will move to */
				FeVariable_GetValueAtNode( velocityField, wallNode_lI, velVector );
				newWallCoordsInRemeshAxis[sideWallNode_I] += velVector[remeshAxis] * dt;
				/* make sure that all the side wall stretches are the same if this is for a
				 * parallelPipedHexa mesh */
				/*
				if ( ( True == atLeastOneLocalSideNodeProcessed ) &&
					Stg_Class_IsInstance( mesh->layout->elementLayout, ParallelPipedHexaEL_Type ) )
				{
					double                     differenceBetweenCurrAndPrev = 0;
					Stream*    errorStream = Journal_Register( Error_Type, self->type );
					differenceBetweenCurrAndPrev = fabs( newWallCoordsInRemeshAxis[ sideWallNode_I ]
						- lastCalculatedNewWallCoordInRemeshAxis );

					Journal_Firewall( differenceBetweenCurrAndPrev < tolerance,
						errorStream,
						"Error - in %s(): While calculating new coords to remesh to in axis %c for mesh "
						"\"%s\" with element layout type %s, found side wall node %u (with IJK={%u,%u,%u}, "
						"local index %u and coord={%.2f,%.2f,%.2f}) has "
						"new coord in remeshAxis of %.2f (= coord %.2f + new vel %.2f in remesh "
						"axis * dt=%.2f), but "
						"already calculated side wall node %u (local index %u) has new wall coord in remesh "
						"axis of %.2f. "
						"This would cause irregular deformation, which is not supported for this "
						"mesh type. Check your Boundary Condition specifications.\n",
						__func__, IJKTopology_DimNumToDimLetter[remeshAxis], mesh->name,
						ParallelPipedHexaEL_Type,
						sideWallNode_I, wallNodeIJK[0], wallNodeIJK[1], wallNodeIJK[2],
						wallNode_lI,
						currNodeCoord[0], currNodeCoord[1], currNodeCoord[2],
						newWallCoordsInRemeshAxis[sideWallNode_I], currNodeCoord[remeshAxis],
						velVector[remeshAxis], dt,
						lastCalculatedWallNode_sideI, lastCalculatedWallNode_lI,
						lastCalculatedNewWallCoordInRemeshAxis );
				}
				*/
				/* Copy values for comparison */
				lastCalculatedWallNode_sideI = sideWallNode_I;
				lastCalculatedWallNode_lI = wallNode_lI;
				lastCalculatedNewWallCoordInRemeshAxis = newWallCoordsInRemeshAxis[sideWallNode_I];
				atLeastOneLocalSideNodeProcessed = True;
			}	
			else {
				if ( minOrMaxFlag == MIN_COORDS ) {
					newWallCoordsInRemeshAxis[sideWallNode_I] = HUGE_VAL;
				}	
				else {
					newWallCoordsInRemeshAxis[ sideWallNode_I ] = -HUGE_VAL;
				}
			}	
			
			/* Increment counter for side wall node index */
			sideWallNode_I++;
		}
	}

	if ( minOrMaxFlag == MIN_COORDS ) {
		mpiOperation = MPI_MIN;
	}	
	else {
		mpiOperation = MPI_MAX;
	}	

	/*comm = CommTopology_GetComm( Mesh_GetCommTopology( mesh, MT_VERTEX ) );*/
	comm = Comm_GetMPIComm( Mesh_GetCommTopology( mesh, MT_VERTEX ) );
	MPI_Allreduce( newWallCoordsInRemeshAxis, newWallCoordsInRemeshAxisGlobal, sideWallNodeCount,
		       MPI_DOUBLE, mpiOperation, comm );

	/* TODO : should really do another iteration over all the nodes, checking identical 	 */

	Memory_Free( newWallCoordsInRemeshAxis );	
}

		
/* This function will provide StGermain the abilty to instantiate (create) this codelet on demand. */
void* _Underworld_MovingMesh_DefaultNew( Name name ) {
	return _Codelet_New(
			sizeof( MeshExtender ),
			Underworld_MovingMesh_Type,
			_Codelet_Delete,
			_Codelet_Print, 
			_Codelet_Copy,
			_Underworld_MovingMesh_DefaultNew,
			_Underworld_MovingMesh_Construct, /* SQ NOTE: Used to be a construct extensions. */
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );		
}
	
/* This function is automatically run by StGermain when this plugin is loaded. The name must be "<plugin-name>_Register". */
Index Underworld_MovingMesh_Register( PluginsManager* pluginsManager ) {
	/* A plugin is only properly registered once it returns the handle provided when submitting a codelet to StGermain. */
	return PluginsManager_Submit( pluginsManager, Underworld_MovingMesh_Type, "0", _Underworld_MovingMesh_DefaultNew );
}
