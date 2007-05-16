/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
** Copyright (c) 2005-2006, Monash Cluster Computing, Building 28, Monash University Clayton Campus,
**	Victoria, 3800, Australia
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**
** Contributors:
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Patrick D. Sunter, Software Engineer, VPAC. (patrick@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	David Stegman, Postdoctoral Fellow, Monash University. (david.stegman@sci.monash.edu.au)
**	Wendy Sharples, PhD Student, Monash University (wendy.sharples@sci.monash.edu.au)
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Lesser General Public
**  License as published by the Free Software Foundation; either
**  version 2.1 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License along with this library; if not, write to the Free Software
**  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**
** $Id: TempNodeHeight.c 340 2006-06-27 00:54:55Z AlanLo $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <assert.h>

const Type ExperimentalUnderworld_CylinderNodeProfiling = "ExperimentalUnderworld_CylinderNodeProfiling";

void ExperimentalUnderworld_NodeTempProfile( PICelleratorContext* context ) {
#if 0
	static FeVariable*          temperatureFe;
	FiniteElement_Mesh*         mesh;
	double*                     nodeCoord;
	double                      lmaxTemp, nodeTemp;
	Index                       lNode_I, newNodeID;
	int                         startProfileNodeID; 
	static Bool                 beenHere              = False;
	Stream*                     stream                = Journal_Register( Info_Type, "TempNodeHeight" );
	char*                       filename;
	double*                     maxTempList = Memory_Alloc_Array( double, context->nproc, "Hold the max temperature of each array");
	unsigned                    rootProc;
	
	rootProc = 0;
	if (!beenHere) {
		Name                 tempFeName;
		tempFeName = Dictionary_GetString_WithDefault( context->dictionary, "TemperatureField", "TemperatureField" );
		
		/* Get TemperatureField FeVariable */
		temperatureFe = (FeVariable*) LiveComponentRegister_Get( context->CF->LCRegister, tempFeName );
		assert( temperatureFe );
	
		/* Set up stream */
		Stg_asprintf( &filename, "CylinderNodeProfiling.%dof%d.dat", context->rank, context->nproc );
		Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, filename );
		Memory_Free( filename );
		Stream_SetAutoFlush( stream, True );

		/* Print header to stream */
		Journal_Printf( stream, "Timestep\tCentralPosition\t\tNodeTemperature\n" );
		beenHere = True;
	}
	mesh = temperatureFe->feMesh;
	
	for( lNode_I = 0; lNode_I < temperatureFe->feMesh->nodeLocalCount ; lNode_I++ ) {
		nodeTemp = FeVariable_GetScalarAtNode( temperatureFe , lNode_I );
		if( lNode_I == 0 ) {
			lmaxTemp = nodeTemp;
			startProfileNodeID = lNode_I;
		}
		if ( nodeTemp > lmaxTemp ) {
			lmaxTemp = nodeTemp;
			startProfileNodeID = lNode_I;
		}
	}
	MPI_Gather( &lmaxTemp, 1, MPI_DOUBLE, maxTempList, 1, MPI_DOUBLE, rootProc, context->communicator );

	Journal_Printf( stream, "%.6g   ******************************\n", (double)context->timeStep );
	nodeCoord = Mesh_CoordAt( mesh , startProfileNodeID );
	nodeTemp  = FeVariable_GetScalarAtNode( temperatureFe, startProfileNodeID );
	Journal_Printf( stream, "\t\tNode Profile in Y direction, starting at (%.6g, %.6g, %.6g)\tat temperature = %g\n",
				nodeCoord[0], nodeCoord[1], nodeCoord[2], nodeTemp );

	newNodeID = mesh->nodeNeighbourTbl[ startProfileNodeID ][ 1 ]; /* 1 = +y direction, 2 = +z direction */
	if( newNodeID >= mesh->nodeLocalCount ) { /* is node NOT a local, thus stop profiling */
		Journal_Printf( stream, "\t\tProfiling has reached the boundary of the local processor\n");
	} else {  /* continue profiling */
		nodeTemp = FeVariable_GetScalarAtNode( temperatureFe, newNodeID );
		while( nodeTemp > 0.0 ) {
			nodeCoord = Mesh_CoordAt( mesh, newNodeID );
			Journal_Printf( stream, "\t\t (%.6g, %.6g, %.6g)\tat temperature = %g\n", nodeCoord[0], nodeCoord[1], nodeCoord[2], nodeTemp );
			newNodeID = mesh->nodeNeighbourTbl[ newNodeID ][ 1 ];
			if( newNodeID >= mesh->nodeLocalCount ) {
                          /* is node not a local i.e. not on processor, tus stop profiling	 */
				Journal_Printf( stream, "\t\tProfiling has reached the boundary of the local processor\n");
				break;
			}
			nodeTemp = FeVariable_GetScalarAtNode( temperatureFe, newNodeID );
		}
	}
	Memory_Free( maxTempList );
	MPI_Barrier( context->communicator );
#endif
}

void _ExperimentalUnderworld_CylinderNodeProfiling_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {

	AbstractContext* context;
	FeVariable*  temperatureField  = Stg_ComponentFactory_ConstructByName( cf, "TemperatureField", FeVariable, True, data );
	FeMesh* mesh     = temperatureField->feMesh;

	context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data );
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, ExperimentalUnderworld_NodeTempProfile );
}


void* _ExperimentalUnderworld_CylinderNodeProfiling_DefaultNew( Name name ) {
	return _Codelet_New(
			sizeof( Codelet ),
			ExperimentalUnderworld_CylinderNodeProfiling,
			_Codelet_Delete,
			_Codelet_Print,
			_Codelet_Copy,
			_ExperimentalUnderworld_CylinderNodeProfiling_DefaultNew,
			_ExperimentalUnderworld_CylinderNodeProfiling_Construct,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}


Index ExperimentalUnderworld_CylinderNodeProfiling_Register( PluginsManager* pluginsManager ) {
	Index result;

	result = PluginsManager_Submit( pluginsManager, ExperimentalUnderworld_CylinderNodeProfiling, "0",
		_ExperimentalUnderworld_CylinderNodeProfiling_DefaultNew );

	return result;
}
