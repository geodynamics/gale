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
** $Id:  $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include "Context.h"
#include "MultiGrid.h"


const Type		StgFEM_MultiGrid_Type = "MultiGrid";
ExtensionInfo_Index	MultiGrid_ContextHandle;

const char*	MULTIGRID_PLUGIN_TAG = "MultiGrid";
const char*	SLES_TO_CONVERT_LIST_TAG = "SLEsToConvert";


Index _StgFEM_MultiGrid_Register( PluginsManager* pluginsMgr ) {
	return PluginsManager_Submit( pluginsMgr, 
				      StgFEM_MultiGrid_Type, 
				      "0", 
				      _StgFEM_MultiGrid_DefaultNew );
}


void* _StgFEM_MultiGrid_DefaultNew( Name name ) {
	return _Codelet_New( sizeof(Codelet), 
			     StgFEM_MultiGrid_Type, 
			     _Codelet_Delete, 
			     _Codelet_Print, 
			     _Codelet_Copy, 
			     _StgFEM_MultiGrid_DefaultNew, 
			     _StgFEM_MultiGrid_Construct, 
			     _StgFEM_MultiGrid_Build, 
			     _Codelet_Initialise, 
			     _Codelet_Execute, 
			     _StgFEM_MultiGrid_Destroy, 
			     name );
}


void _StgFEM_MultiGrid_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	FiniteElementContext*	feCtx;
	MultiGrid_Context*	mgCtx;

	assert( component );
	assert( cf );

	/* Retrieve context. */
	feCtx = (FiniteElementContext*)Stg_ComponentFactory_ConstructByName( cf, "context", FiniteElementContext, True, data ); 


	/*
	** Register all entry points, add extensions to the FE context and whatever else needs modifying.
	*/

	Journal_DPrintf( StgFEM_Debug, "In: %s( void* )\n", __func__ );

	/* Add extensions to the context and get them. */
	MultiGrid_ContextHandle = ExtensionManager_Add( feCtx->extensionMgr, 
							StgFEM_MultiGrid_Type, 
							sizeof(MultiGrid_Context) );

	mgCtx = ExtensionManager_Get( feCtx->extensionMgr, 
				      feCtx, 
				      MultiGrid_ContextHandle );


	/*
	** Clear the context.
	*/

	mgCtx->ctx.nMappings = 0;
	mgCtx->ctx.mappings = NULL;
	mgCtx->ctx.nInfos = 0;
	mgCtx->ctx.infos = NULL;
	mgCtx->nSLEs = 0;
	mgCtx->sles = NULL;
	mgCtx->nLevels = NULL;
	mgCtx->nUpIts = NULL;
	mgCtx->nDownIts = NULL;
	mgCtx->nCycles = NULL;
	mgCtx->nFinalIts = NULL;
	mgCtx->upSmoothers = NULL;
	mgCtx->downSmoothers = NULL;
	
	
	/*
	** Read in all required parameters to the context from the dictionary.
	*/
	
	{
		Dictionary_Entry_Value*	mgEntry=NULL;
		Dictionary_Entry_Value*	slesToConvertList=NULL;
		Index				nSLEsToConvert=0;
		Index				sle_i;
		Stream*				warningStr = Journal_Register( Error_Type, CURR_MODULE_NAME );
		
		/* Parse in the name of the SLE to multi-grid-ize. */
		mgEntry = Dictionary_Get( feCtx->dictionary, (char*)MULTIGRID_PLUGIN_TAG );
		if( mgEntry == NULL ) {
			Journal_Printf( warningStr, 
					"Warning - in %s: Plugin \"%s\" loaded, but no \"%s\" tag found "
					"in dictionary. Not converting any SLEs to multigrid.\n", 
					__func__, CURR_MODULE_NAME, MULTIGRID_PLUGIN_TAG );
			return;
		}
		
		slesToConvertList = Dictionary_Entry_Value_GetMember( mgEntry, (char*)SLES_TO_CONVERT_LIST_TAG );
		if( slesToConvertList == NULL ) {
			Journal_Printf( warningStr, 
					"Warning - in %s: Plugin \"%s\" loaded, \"%s\" tag found, "
					"but it doesn't contain a list \"%s\" of SLEs to convert to multigrid.\n",
					__func__, CURR_MODULE_NAME, MULTIGRID_PLUGIN_TAG, SLES_TO_CONVERT_LIST_TAG );
			return;
		}
		
		nSLEsToConvert = Dictionary_Entry_Value_GetCount( slesToConvertList );
		if( nSLEsToConvert == 0 ) {
			Journal_Printf( warningStr, 
					"Warning - in %s: Plugin \"%s\" loaded, \"%s\" tag found, " 
					"but list \"%s\" has 0 entries.\n", 
					__func__, CURR_MODULE_NAME, MULTIGRID_PLUGIN_TAG, SLES_TO_CONVERT_LIST_TAG );
			return;
		}
		
		for( sle_i = 0; sle_i < nSLEsToConvert; sle_i++ ) {
			Dictionary_Entry_Value*	currSLE_ConversionRequest = NULL;
			Dictionary*		currSLE_ConversionRequestDict = NULL;
			Name			targetedSLE_Name = NULL;
			void*			targetedSLE = NULL;
			unsigned int		mgLevel = 0;
			unsigned*		upIts;
			unsigned*		downIts;
			unsigned*		cycles;
			unsigned		finalIts;
			unsigned		level_i;
			unsigned		nAllDownIts, nAllUpIts;
			char**			upSmoothers;
			char**			downSmoothers;
			char*			allUpSmoother;
			char*			allDownSmoother;
			
			
			currSLE_ConversionRequest = Dictionary_Entry_Value_GetElement( slesToConvertList, sle_i );
			currSLE_ConversionRequestDict = Dictionary_Entry_Value_AsDictionary( currSLE_ConversionRequest );
			
			/* Find the SLE name. */
			targetedSLE_Name = Dictionary_GetString( currSLE_ConversionRequestDict, "SLE" );
			if( targetedSLE_Name == NULL ) {
				Journal_Printf( warningStr, 
						"Warning - in %s(): One of the multigrid conversion requests doesn't " 
						"have an \"SLE\" parameter describing which SLE should have multigrid applied to it. " 
						"Skipping this entry.\n", 
						__func__ );
				continue;
			}
			
			/* Search the FE context for that SLE. */
			targetedSLE = Stg_ObjectList_Get( feCtx->slEquations, targetedSLE_Name );
			if( targetedSLE == NULL ) {
				Journal_Printf( warningStr, 
						"Warning - in %s(): You requested appling MultiGrid to SLE \"%s\", but this " 
						"SLE wasn't found in the SLE register on the context.", 
						__func__, targetedSLE_Name );
				Journal_Printf( warningStr, 
						" Did you remember to load the %s plugin _after_ the one that constructs the App? " 
						"Skipping this conversion.\n", 
						CURR_MODULE_NAME );
				Journal_Printf( warningStr, 
						"(SLEs currently loaded in the FE context are: " );
				Stg_ObjectList_PrintAllEntryNames( feCtx->slEquations, warningStr );	
				Journal_Printf( warningStr, ")\n" );
				continue;
			}
			
			/* Find out the level requested (and other things). */
			mgLevel = Dictionary_GetUnsignedInt_WithDefault( currSLE_ConversionRequestDict, "level", 1 );  /* ambiguous documentation level or levels ??!! */
			if( Dictionary_GetUnsignedInt_WithDefault( currSLE_ConversionRequestDict, "levels", (unsigned)-1 ) != (unsigned)-1 ) {
				mgLevel = Dictionary_GetUnsignedInt_WithDefault( currSLE_ConversionRequestDict, "levels", 1 );
			}
			finalIts = Dictionary_GetUnsignedInt_WithDefault( currSLE_ConversionRequestDict, "finalIts", 1 );
			nAllDownIts = Dictionary_GetUnsignedInt_WithDefault( currSLE_ConversionRequestDict, "downIts", 1 );
			nAllUpIts = Dictionary_GetUnsignedInt_WithDefault( currSLE_ConversionRequestDict, "upIts", 1 );
			allUpSmoother = Dictionary_GetString_WithDefault( currSLE_ConversionRequestDict, "upSmoother", "sor" );
			allDownSmoother = Dictionary_GetString_WithDefault( currSLE_ConversionRequestDict, "downSmoother", "sor" );
			
			downIts = Memory_Alloc_Array( unsigned, mgLevel + 1, "MultiGrid" );
			upIts = Memory_Alloc_Array( unsigned, mgLevel + 1, "MultiGrid" );
			cycles = Memory_Alloc_Array( unsigned, mgLevel + 1, "MultiGrid" );
			upSmoothers = Memory_Alloc_Array( char*, mgLevel + 1, "MultiGrid" );
			downSmoothers = Memory_Alloc_Array( char*, mgLevel + 1, "MultiGrid" );
			
			for( level_i = 0; level_i <= mgLevel; level_i++ ) {
				char		levelName[10];
				Dictionary*	levelEntry;
				
				/* This is the name of the dictionary entry for this level. */
				sprintf( levelName, "level_%d", level_i );
				levelEntry = Dictionary_Entry_Value_AsDictionary( 
					Dictionary_Get( currSLE_ConversionRequestDict, levelName ) );
				if( levelEntry ) {
					fprintf(stderr,"*** MGplugin: Reading additional parameters for level %d\n",level_i);
					downIts[level_i] = Dictionary_GetUnsignedInt_WithDefault( levelEntry, "downIts", nAllDownIts );
					upIts[level_i] = Dictionary_GetUnsignedInt_WithDefault( levelEntry, "upIts", nAllUpIts );
					cycles[level_i] = Dictionary_GetUnsignedInt_WithDefault( levelEntry, "cycles", 1 );
					upSmoothers[level_i] = StG_Strdup( Dictionary_GetString_WithDefault( levelEntry, 
													 "upSmoother", 
													 allUpSmoother ) );
					downSmoothers[level_i] = StG_Strdup( Dictionary_GetString_WithDefault( levelEntry, 
													   "downSmoother", 
													   allDownSmoother ) );
				}
				else {
					if( level_i < mgLevel ) {
						downIts[level_i] = nAllDownIts;
						upIts[level_i] = nAllUpIts;
						cycles[level_i] = 1;
						upSmoothers[level_i] = StG_Strdup( allUpSmoother );
						downSmoothers[level_i] = StG_Strdup( allDownSmoother );
					}
					else {
						downIts[level_i] = 1;
						upIts[level_i] = 0;
						cycles[level_i] = 1;
						upSmoothers[level_i] = StG_Strdup( "default" );
						downSmoothers[level_i] = StG_Strdup( "default" );
					}
				}
			}
			
			/* This is a "level_mg_coarse" lookup which sets coarse grid parameters and 
				override all other settings. This would let us give coarse grid settings which would not
				be changed if the user messes about with the grid resolution */
			
			{
				char		levelName[10];
				Dictionary*	levelEntry;
			
				level_i = mgLevel;
				sprintf( levelName, "level_mg_coarse");
				levelEntry = Dictionary_Entry_Value_AsDictionary( 
					Dictionary_Get( currSLE_ConversionRequestDict, levelName ) );
				if( levelEntry ) {
					fprintf(stderr,"*** MGplugin: Reading additional parameters for level %d\n",level_i);
					downIts[level_i] = Dictionary_GetUnsignedInt_WithDefault( levelEntry, "downIts", downIts[level_i] );
					upIts[level_i] = Dictionary_GetUnsignedInt_WithDefault( levelEntry, "upIts", upIts[level_i] );
					cycles[level_i] = Dictionary_GetUnsignedInt_WithDefault( levelEntry, "cycles", cycles[level_i] );
					upSmoothers[level_i] = StG_Strdup( Dictionary_GetString_WithDefault( levelEntry, 
													 "upSmoother", 
													 upSmoothers[level_i] ) );
					downSmoothers[level_i] = StG_Strdup( Dictionary_GetString_WithDefault( levelEntry, 
													   "downSmoother", 
													   downSmoothers[level_i] ) );
				}
			}	
				
			
			
			

			/* Add this info to the multi-grid context. */
			mgCtx->sles = Memory_Realloc_Array( mgCtx->sles, SystemLinearEquations*, mgCtx->nSLEs + 1 );
			mgCtx->sles[mgCtx->nSLEs] = targetedSLE;
			mgCtx->nLevels = Memory_Realloc_Array( mgCtx->nLevels, unsigned, mgCtx->nSLEs + 1 );
			mgCtx->nLevels[mgCtx->nSLEs] = mgLevel;
			mgCtx->nUpIts = Memory_Realloc_Array( mgCtx->nUpIts, unsigned*, mgCtx->nSLEs + 1 );
			mgCtx->nUpIts[mgCtx->nSLEs] = upIts;
			mgCtx->nDownIts = Memory_Realloc_Array( mgCtx->nDownIts, unsigned*, mgCtx->nSLEs + 1 );
			mgCtx->nDownIts[mgCtx->nSLEs] = downIts;
			mgCtx->nCycles = Memory_Realloc_Array( mgCtx->nCycles, unsigned*, mgCtx->nSLEs + 1 );
			mgCtx->nCycles[mgCtx->nSLEs] = cycles;
			mgCtx->nFinalIts = Memory_Realloc_Array( mgCtx->nFinalIts, unsigned, mgCtx->nSLEs + 1 );
			mgCtx->nFinalIts[mgCtx->nSLEs] = finalIts;
			mgCtx->upSmoothers = Memory_Realloc_Array( mgCtx->upSmoothers, char**, mgCtx->nSLEs + 1);
			mgCtx->upSmoothers[mgCtx->nSLEs] = upSmoothers;
			mgCtx->downSmoothers = Memory_Realloc_Array( mgCtx->downSmoothers, char**, mgCtx->nSLEs + 1);
			mgCtx->downSmoothers[mgCtx->nSLEs] = downSmoothers;
			mgCtx->nSLEs++;

			/* For this SLE, reach in and change the desired FeEquationNumbers to not remove their BCs. */
			{
				unsigned		nSMs;
				StiffnessMatrix**	sms;
				unsigned		sm_i;

				((SystemLinearEquations*)targetedSLE)->bcRemoveQuery = True;
				SystemLinearEquations_MG_SelectStiffMats( targetedSLE, &nSMs, &sms );
				for( sm_i = 0; sm_i < nSMs; sm_i++ ) {
					StiffnessMatrix*	sm = sms[sm_i];

					/* Disable BCs for this SM. */
					if( sm->rowVariable ) {
						sm->rowVariable->eqNum->removeBCs = True;
					}
					if( sm->columnVariable ) {
						sm->columnVariable->eqNum->removeBCs = True;
					}
				}
			}
		}
	}
	
	/* Enable the plugin's context. */
	MultiGrid_SetContext( &mgCtx->ctx );
}


void _StgFEM_MultiGrid_Build( void* component, void* data ) {
	FiniteElementContext*	feCtx;
	MultiGrid_Context*	mgCtx;
	unsigned		sle_i;

	assert( component );
	assert( data );

	/* Get the contexts. */
	feCtx = (FiniteElementContext*)data;
	mgCtx = ExtensionManager_Get( feCtx->extensionMgr, 
				      feCtx, 
				      MultiGrid_ContextHandle );

	
	/*
	** Convert all the target SLEs to use MG.
	*/

	for( sle_i = 0; sle_i < mgCtx->nSLEs; sle_i++ ) {
		MultiGrid_ConvertSLE( mgCtx->sles[sle_i], 
				      mgCtx->nLevels[sle_i], 
				      mgCtx->nUpIts[sle_i], mgCtx->nDownIts[sle_i], mgCtx->nCycles[sle_i], 
				      mgCtx->nFinalIts[sle_i], 
				      mgCtx->upSmoothers[sle_i], mgCtx->downSmoothers[sle_i] );
	}
}


void _StgFEM_MultiGrid_Destroy( void* component, void* data ) {
	FiniteElementContext*	feCtx;
	MultiGrid_Context*	mgCtx;

	assert( component );
	assert( data );

	/* Get the contexts. */
	feCtx = (FiniteElementContext*)data;
	mgCtx = ExtensionManager_Get( feCtx->extensionMgr, 
				      feCtx, 
				      MultiGrid_ContextHandle );


	/*
	** Delete the plugin's context.
	*/

	/* TODO
	   for( sle_i = 0; sle_i < mgCtx->nSLEs; sle_i++ ) {
	   FreeArray2D( mgCtx->nLevels[sle_i], mgCtx->upSmoothers[sle_i] );
	   FreeArray2D( mgCtx->nLevels[sle_i], mgCtx->downSmoothers[sle_i] );
	   } */

	MultiGrid_DeleteContext( &mgCtx->ctx );
	FreeArray( mgCtx->sles );
	FreeArray( mgCtx->nLevels );
	FreeArray2D( mgCtx->nSLEs, mgCtx->nUpIts );
	FreeArray2D( mgCtx->nSLEs, mgCtx->nDownIts );
	FreeArray2D( mgCtx->nSLEs, mgCtx->nCycles );
	FreeArray( mgCtx->nFinalIts );
	FreeArray( mgCtx->upSmoothers );
	FreeArray( mgCtx->downSmoothers );
}
