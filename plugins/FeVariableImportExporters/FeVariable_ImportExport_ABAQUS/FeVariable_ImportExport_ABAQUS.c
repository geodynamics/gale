/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**	Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
** Contributors:
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
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
** $Id: FeVariable_ImportExport_ABAQUS.c 742 2007-02-13 02:44:01Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>

#include <string.h>

const Type   StgFEM_FeVariable_ImportExport_ABAQUS_Type = "StgFEM_FeVariable_ImportExport_ABAQUS";
const char*  ABAQUS_ImportExportType = "ABAQUS";


void FeVariable_ReadNodalValuesFromFile_ABAQUS( void* _feVariable, const char* prefixStr, unsigned int timeStep ) {
	FeVariable*        feVariable = (FeVariable*)_feVariable;
	char*              filename;
	Node_LocalIndex    lNode_I = 0;
	Node_GlobalIndex   gNode_I = 0;
	Node_GlobalIndex   gNodeCount_I = 0;
	Dof_Index          dof_I;
	Dof_Index          dofAtEachNodeCount;
	FILE*              inputFile;
	double             variableVal;
	char               lineString[1000];
	const unsigned int MAX_LINE_LENGTH = 1000;
	Processor_Index    proc_I=0;
	Dimension_Index    dim_I=0;
	BlockGeometry*     geometry = (BlockGeometry*)feVariable->feMesh->layout->elementLayout->geometry;
	char*              matchString;
	Index              currentFileLine = 0;
	Coord              localGeometryMin;
	Coord              localGeometryMax;
	Stream*            debugStream = Journal_Register( Debug_Type, StgFEM_FeVariable_ImportExport_ABAQUS_Type );
	Bool               nodeG2LBuiltTemporarily = False;
	
	Journal_DPrintf( debugStream, "In %s(): for FeVariable \"%s\"\n", __func__, feVariable->name );
	Stream_Indent( debugStream );
	
	/* Necessary for now since we need to update the geometry min and max - see comment below */
	geometry = Stg_CheckType( geometry, BlockGeometry );

	/*                                                prefix            feVariable->name        . 00000 .  dat \0 */
	filename = Memory_Alloc_Array_Unnamed( char, strlen(prefixStr) + strlen(feVariable->name) + 1 + 5 + 1 + 3 + 1 );
	sprintf( filename, "%s%s.%.5u.rpt", prefixStr, feVariable->name, timeStep );

	/* TODO May need/want to change to MPI file stuff */
	
	/* This loop used to stop 2 processors trying to open the file at the same time, which
	  * seems to cause problems */
	for ( proc_I = 0; proc_I < feVariable->feMesh->layout->decomp->nproc; proc_I++ ) {
		MPI_Barrier( feVariable->feMesh->layout->decomp->communicator );
		if ( proc_I == feVariable->feMesh->layout->decomp->rank ) {	
			inputFile = fopen( filename, "r" );
		}
	}

	if ( False == inputFile ) {
		Stream*    errorStr = Journal_Register( Error_Type, feVariable->type );
		Journal_Printf( errorStr, "Error- in %s(), for feVariable \"%s\": Couldn't find checkpoint file with "
			"prefix \"%s\", timestep %d - thus full filename \"%s\" - aborting.\n", __func__, feVariable->name,
			prefixStr, timeStep, filename );
		exit(EXIT_FAILURE);	
	}

	/* This is where we skip over the ABAQUS header stuff */
	while ( !feof(inputFile) ) {
		currentFileLine++;
		fgets( lineString, MAX_LINE_LENGTH, inputFile );
		matchString = strstr( lineString, "            Node" );
		if ( matchString != NULL ) {
			/* Grab the "Label" and the "----" lines */ 
			fgets( lineString, MAX_LINE_LENGTH, inputFile );
			fgets( lineString, MAX_LINE_LENGTH, inputFile );
			currentFileLine += 2;
			/* Ok, now we're ready to start reading the actual field values */
			break;
		}
	}	
	Journal_DPrintf( debugStream, "Skipped %u lines of ABAQUS header info...\n", currentFileLine );

	dofAtEachNodeCount = feVariable->fieldComponentCount;

	/* Need to re-set the geometry here, in case we're loading from a checkpoint that had compression/squashing BCs,
		and hence ended up with a smaller mesh than the original */
	for ( dim_I = 0; dim_I < 3; dim_I++ ) {
		localGeometryMin[dim_I] = HUGE_VAL;
		localGeometryMax[dim_I] = -HUGE_VAL;
	}

	/* Modification: build the temporary global tables for speed here. Important for parallel runs.
 	 *  -- PatrickSunter, 13 Feb 2007 */
	if ( (feVariable->feMesh->buildTemporaryGlobalTables == True) && (feVariable->feMesh->nodeG2L == 0) ) {
		feVariable->feMesh->nodeG2L = MeshDecomp_BuildNodeGlobalToLocalMap( feVariable->feMesh->layout->decomp );
		nodeG2LBuiltTemporarily = True;
	}

	Journal_DPrintf( debugStream, "Processing Nodal Data...\n", currentFileLine );
	Stream_Indent( debugStream );
	/* Note: in ABAQUS, the number of lines of nodal values from hereon is always == the number of global nodes */
	for ( gNodeCount_I = 0; gNodeCount_I < feVariable->feMesh->nodeGlobalCount; gNodeCount_I++ ) {
		fscanf( inputFile, "%u ", &gNode_I );
		/* Note: ABAQUS has same layout of global node indices as StgFEM, except it indexes starting from 1 - thus we
		 * need to subtract 1 here */
		gNode_I -= 1;
		
		lNode_I = Mesh_NodeMapGlobalToLocal( feVariable->feMesh, gNode_I );
		Journal_DPrintfL( debugStream, 3, "Found info for global node %u, local node %u:\n", gNode_I, lNode_I );
		Stream_Indent( debugStream );

		if ( lNode_I != Mesh_Node_Invalid( feVariable->feMesh ) ) {
			/* Note: until we have proper mesh geometry, topology etc checkpointing, we re-load the 
			node co-ords from the feVariable file - and also update the geometry */
			fscanf( inputFile, "%lg %lg ",
				&feVariable->feMesh->nodeCoord[lNode_I][0],
				&feVariable->feMesh->nodeCoord[lNode_I][1] );
			if ( feVariable->fieldComponentCount == 3 ) {	
				fscanf( inputFile, "%lg",
					&feVariable->feMesh->nodeCoord[lNode_I][2] );
			}	
			else {
				feVariable->feMesh->nodeCoord[lNode_I][2] = 0.0;
			}

			Journal_DPrintfL( debugStream, 3, "read coord (%.3f, %.3f, %.3f)\n", 
				feVariable->feMesh->nodeCoord[lNode_I][0],
				feVariable->feMesh->nodeCoord[lNode_I][1],
				feVariable->feMesh->nodeCoord[lNode_I][2] );

			for ( dof_I = 0; dof_I < dofAtEachNodeCount; dof_I++ ) {
				fscanf( inputFile, "%lg ", &variableVal );
				DofLayout_SetValueDouble( feVariable->dofLayout, lNode_I, dof_I, variableVal );
				Journal_DPrintfL( debugStream, 3, "read dof %u: %g\n", dof_I, variableVal ); 
				/* TODO: Hack for now - ABAQUS only uses initial coords, so assume this is displacement and add */
				feVariable->feMesh->nodeCoord[lNode_I][dof_I] += variableVal;
				Journal_DPrintfL( debugStream, 3, "TODO: using HACK assumption of disp.:- updating nodeCoord[%u] to %.3f\n",
					dof_I, feVariable->feMesh->nodeCoord[lNode_I][dof_I] );
			}

			for ( dim_I = 0; dim_I < 3; dim_I++ ) {
				if ( feVariable->feMesh->nodeCoord[lNode_I][dim_I] < localGeometryMin[dim_I] ) {
					localGeometryMin[dim_I] = feVariable->feMesh->nodeCoord[lNode_I][dim_I];
				}
				else if ( feVariable->feMesh->nodeCoord[lNode_I][dim_I] > localGeometryMax[dim_I] ) {
					localGeometryMax[dim_I] = feVariable->feMesh->nodeCoord[lNode_I][dim_I];
				}
			}
			
		}
		else {
			Journal_DPrintfL( debugStream, 3, "not on current proc -> skipping\n" );
			fgets( lineString, MAX_LINE_LENGTH, inputFile );
		}
		Stream_UnIndent( debugStream );
		currentFileLine++;
	}			
	fclose( inputFile );
	Stream_UnIndent( debugStream );

	if ( nodeG2LBuiltTemporarily ) {
		Memory_Free( feVariable->feMesh->nodeG2L );
		feVariable->feMesh->nodeG2L = NULL;
	}

	/* Since we could be loading in parallel, need to find global min and max of geometry */
	for ( dim_I = 0; dim_I < 3; dim_I++ ) {
		MPI_Allreduce( localGeometryMin, geometry->min, 3, MPI_DOUBLE, MPI_MIN, 
			feVariable->feMesh->layout->decomp->communicator );
		MPI_Allreduce( localGeometryMax, geometry->max, 3, MPI_DOUBLE, MPI_MAX, 
			feVariable->feMesh->layout->decomp->communicator );
	}

	Journal_DPrintf( debugStream, "Recalculated global field min coord as (%.3f, %.3f, %.3f)\n",
		geometry->min[0], geometry->min[1], geometry->min[2] );
	Journal_DPrintf( debugStream, "Recalculated global field max coord as (%.3f, %.3f, %.3f)\n",
		geometry->max[0], geometry->max[1], geometry->max[2] );

	if ( Stg_Class_IsInstance( feVariable->feMesh->layout->elementLayout, ParallelPipedHexaEL_Type ) ) {
		ParallelPipedHexaEL_UpdateGeometryPartitionInfo( feVariable->feMesh->layout->elementLayout, feVariable->feMesh->layout->decomp );
	}

	Stream_UnIndent( debugStream );
}			


void FeVariable_SaveNodalValuesToFile_ABAQUS( void* _feVariable, const char* prefixStr, unsigned int timeStep ) {
	FeVariable*      feVariable = (FeVariable*)_feVariable;
	Stream*          errorStream = Journal_Register( Error_Type, StgFEM_FeVariable_ImportExport_ABAQUS_Type );

	Journal_Firewall( 0, errorStream, "Error - in %s(), for FeVariable \"%s\": function not implemented yet.\n",
		__func__, feVariable->name );
}


void _StgFEM_FeVariable_ImportExport_ABAQUS_Construct( void* componment, Stg_ComponentFactory* cf, void* data ) {
	AbstractContext* context;
	FeVariable_ImportExportInfo*    abaqusImportExportInfo = NULL;

	context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data );
	
	abaqusImportExportInfo = Memory_Alloc( FeVariable_ImportExportInfo, "ABAQUS_ImportExportInfo" );
	abaqusImportExportInfo->readNodalValuesFromFile = FeVariable_ReadNodalValuesFromFile_ABAQUS;
	abaqusImportExportInfo->saveNodalValuesToFile = FeVariable_SaveNodalValuesToFile_ABAQUS;

	Stg_ObjectList_PointerAppend( 
		FeVariable_FileFormatImportExportList, 
		abaqusImportExportInfo,
		ABAQUS_ImportExportType,
		FeVariable_ImportExportInfo_Delete,
		FeVariable_ImportExportInfo_Print,
		FeVariable_ImportExportInfo_Copy );
}

void* _StgFEM_FeVariable_ImportExport_ABAQUS_DefaultNew( Name name ) {
	return _Codelet_New(
			sizeof( Codelet ),
			StgFEM_FeVariable_ImportExport_ABAQUS_Type,
			_Codelet_Delete,
			_Codelet_Print,
			_Codelet_Copy,
			_StgFEM_FeVariable_ImportExport_ABAQUS_DefaultNew,
			_StgFEM_FeVariable_ImportExport_ABAQUS_Construct,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}
   
Index StgFEM_FeVariable_ImportExport_ABAQUS_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, StgFEM_FeVariable_ImportExport_ABAQUS_Type, "0", _StgFEM_FeVariable_ImportExport_ABAQUS_DefaultNew );
}


