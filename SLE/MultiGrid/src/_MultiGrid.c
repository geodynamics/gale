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
#include <stdarg.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>
#include <StGermain/StGermain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/LinearAlgebra/LinearAlgebra.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"

#include "units.h"
#include "types.h"
#include "shortcuts.h"
#include "MGCoarsener.h"
#include "MGCoarsener_RegCartesian.h"
#include "_MultiGrid.h"


/* A global context for MG.  Probably bad coding, but it'll work for now.  Anyway, the MG plugin needs to set this 
   when it starts up. */
MGContext*	mgCtx;


/*
** DEBUG INCLUDES
*/

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscmg.h>
#include <StGermain/compatibility/petsccompat.h>


/*** debug */
void dumpBinMat( void* mat, const char* name, ... ) {
	va_list		ap;
        PetscViewer     viewer;
	char*		modName;

	va_start( ap, name );
	modName = Memory_Alloc_Array( char, strlen( name ) + 1, "MultiGrid" );
	vsprintf( modName, name, ap );

        PetscViewerBinaryOpen( PETSC_COMM_WORLD, modName, PETSC_FILE_CREATE, &viewer );
        MatView( (Mat)mat, viewer );
        PetscViewerDestroy( viewer );

	FreeArray( modName );
	va_end( ap );
}


/*** debug */
void dumpBinVec( void* vec, const char* name, ... ) {
	va_list		ap;
        PetscViewer     viewer;
	char*		modName;

	va_start( ap, name );
	modName = Memory_Alloc_Array( char, strlen( name ) + 1, "MultiGrid" );
	vsprintf( modName, name, ap );

        PetscViewerBinaryOpen( PETSC_COMM_WORLD, modName, PETSC_FILE_CREATE, &viewer );
        VecView( (Vec)vec, viewer );
        PetscViewerDestroy( viewer );

	FreeArray( modName );
	va_end( ap );
}


/*** debug */
void dumpAscMat( void* mat, const char* name, ... ) {
	va_list		ap;
        PetscViewer     viewer;
	char*		modName;

	va_start( ap, name );
	modName = Memory_Alloc_Array( char, strlen( name ) + 1, "MultiGrid" );
	vsprintf( modName, name, ap );

        PetscViewerASCIIOpen( PETSC_COMM_WORLD, modName, &viewer );
	MatView( (Mat)mat, viewer );
        PetscViewerDestroy( viewer );

	FreeArray( modName );
	va_end( ap );
}


/*** debug */
void dumpAscVec( void* vec, const char* name, ... ) {
	va_list		ap;
        PetscViewer     viewer;
	char*		modName;

	va_start( ap, name );
	modName = Memory_Alloc_Array( char, strlen( name ) + 1, "MultiGrid" );
	vsprintf( modName, name, ap );

        PetscViewerASCIIOpen( PETSC_COMM_WORLD, modName, &viewer );
        VecView( (Vec)vec, viewer );
        PetscViewerDestroy( viewer );

	FreeArray( modName );
	va_end( ap );
}


/*
** DEBUG FUNCTION
*/

void dumpVec( Vector* vec, const char* fn ) {
	FILE*		fp;
	PetscInt	size;
	PetscScalar*	array;
	unsigned	el_i;

	VecGetSize( (Vec)vec, &size );
	VecGetArray( (Vec)vec, &array );

	fp = fopen( fn, "w" );
	fprintf( fp, "%d\n", size );
	for( el_i = 0; el_i < size; el_i++ ) {
		fprintf( fp, "%g\n", array[el_i] );
	}
	fclose( fp );

	VecRestoreArray( (Vec)vec, &array );
}


/*
** DEBUG FUNCTION
*/

void dumpMat( Mat mat, const char* fn, int rank ) {
	FILE*	fp;
	MatInfo	matInfo;
	unsigned*	colInds;
	double*	entries;
	unsigned	nCols;
	unsigned	nGRows;
	unsigned	row_i, col_i;
	
	fp = fopen( fn, "w" );
	MatGetInfo( mat, MAT_GLOBAL_MAX, &matInfo );
	nGRows = matInfo.rows_global;
	MatGetInfo( mat, MAT_LOCAL, &matInfo );
	nGRows -= matInfo.rows_local;
	fprintf( fp, "%d %d\n", (unsigned)matInfo.rows_local, (unsigned)matInfo.columns_local );
	for( row_i = 0; row_i < matInfo.rows_local; row_i++ ) {
		MatGetRow( mat, row_i + rank * nGRows, &nCols, (const int**)&colInds, (const double**)&entries );
		fprintf( fp, "%d", nCols );
		for( col_i = 0; col_i < nCols; col_i++ ) {
			fprintf( fp, " %d %g", colInds[col_i], entries[col_i] );
		}
		fprintf( fp, "\n" );
		MatRestoreRow( mat, row_i, &nCols, (const int**)&colInds, (const double**)&entries );
	}
	fclose( fp );
}


void dumpMatrices( unsigned handle ) {
	MGInfo*		info = mgCtx->infos + handle;
	char		fn[10];
	unsigned	level_i;
	int	       	rank;

	MPI_Comm_rank( info->stiffMat->rowVariable->feMesh->layout->decomp->communicator, &rank );

	strcpy( fn, "mat.?.?" );
	
	/* Dump original. */
	fn[4] = '0' + rank;
	fn[6] = '0';
	dumpMat( (Mat)info->stiffMat->matrix, fn, rank );

	/* Dump levels. */
	for( level_i = 0; level_i < info->nLevels; level_i++ ) {
		fn[6] = '1' + level_i;
		dumpMat( (Mat)info->smoothers[level_i], fn, rank );
	}
	
	/* Dump restriction operators. */
	strcpy( fn, "rOp.?.?" );
	fn[4] ='0' + rank;
	for( level_i = 0; level_i < info->nLevels; level_i++ ) {
		fn[6] = '0' + level_i;
		dumpMat( (Mat)info->rOps[level_i], fn, rank );
	}

	/* Dump projection operators. */
	strcpy( fn, "pOp.?.?" );
	fn[4] = '0' + rank;
	for( level_i = 0; level_i < info->nLevels; level_i++ ) {
		fn[6] = '0' + level_i;
		dumpMat( (Mat)info->pOps[level_i], fn, rank );
	}
}


/*
** DEBUG FUNCTION
*/

void dumpState( MatrixSolver* solver ) {
	PC		pc;
	unsigned	nLevels;
	unsigned	level_i;
	
	KSPGetPC( (KSP)solver, &pc );
	PCMGGetLevels( pc, &nLevels );
	
	printf( "*** MG: State dump:\n" );
	printf( "*** MG: \tNumber of levels: %d\n", nLevels );
	
	if( nLevels ) {
		unsigned	nDownIts, nUpIts;
		KSP		down, up;
		
		/* Coarsest level. */
		MatrixSolver_MG_GetDownSmoother( solver, 0, (MatrixSolver**)&down );
		KSPGetTolerances( down, PETSC_NULL, PETSC_NULL, PETSC_NULL, &nDownIts );
		printf( "*** MG: \tState of level 0:\n" );
		printf( "*** MG: \t\tNumber of down iterations: %d\n", nDownIts );
		
		for( level_i = 1; level_i < nLevels; level_i++ ) {
			MatrixSolver_MG_GetDownSmoother( solver, level_i, (MatrixSolver**)&down );
			MatrixSolver_MG_GetUpSmoother( solver, level_i, (MatrixSolver**)&up );
			KSPGetTolerances( down, PETSC_NULL, PETSC_NULL, PETSC_NULL, &nDownIts );
			KSPGetTolerances( up, PETSC_NULL, PETSC_NULL, PETSC_NULL, &nUpIts );
			printf( "*** MG: \tState of level %d:\n", level_i );
			printf( "*** MG: \t\tNumber of down iterations: %d\n", nDownIts );
			printf( "*** MG: \t\tNumber of up iterations: %d\n", nUpIts );
		}
	}
}


void MultiGrid_SetContext( MGContext* ctx ) {
	mgCtx = ctx;
}


void MultiGrid_ConvertSLE( SystemLinearEquations* sle, 
			   unsigned nLevels, 
			   unsigned* nUpIts, unsigned* nDownIts, unsigned* nCycles, 
			   unsigned nFinalIts, 
			   char** upSmoothers, char** downSmoothers )
{
	unsigned			nSMs;
	StiffnessMatrix**	sms = NULL;
	
	
	assert( mgCtx );
	
	if( sle->mgEnabled ) {
		fprintf( stderr, "Warning: MG already enabled on this SLE.\n" );
		return;
	}
	
	
	/*
	** Set the required values on the SLE such that it will use MG.  We'll need to have it give us a list of its
	** stiffness matrices that are capable of having MG applied, then build the required number of levels on a grid
	** mapping.
	*/
	
	SystemLinearEquations_MG_SelectStiffMats( sle, &nSMs, &sms );
	
	/* If there are none to use then don't apply MG. */
	if( nSMs ) {
		unsigned	sm_i;
		
		for( sm_i = 0; sm_i < nSMs; sm_i++ ) {
			StiffnessMatrix*	sm = sms[sm_i];
			unsigned			handle;
			
			/* NOTE: I have no idea how to handle situations where the row and column variables aren't the same. */
			if( sm->rowVariable != sm->columnVariable ) {
				fprintf( stderr, "Warning: Bad stiffness matrix (row and column variables aren't the same \
							   and I'm uneducated).\n" );
				continue;
			}
			
			/* Create a new 'MGInfo' for this stiffness matrix. */
			{
				MGInfo*		info;
				MGGridMapping*	mapping;
				unsigned	map_i;
				
				mgCtx->infos = Memory_Realloc_Array( mgCtx->infos, MGInfo, mgCtx->nInfos + 1 );
				info = mgCtx->infos + mgCtx->nInfos;
				handle = mgCtx->nInfos;
				mgCtx->nInfos++;
				
				info->stiffMat = sm;
				info->nLevels = nLevels;
				info->nUpIts = Memory_Alloc_Array( unsigned, nLevels + 1, "MultiGrid" );
				memcpy( info->nUpIts, nUpIts, sizeof(unsigned) * (nLevels + 1) );
				info->nDownIts = Memory_Alloc_Array( unsigned, nLevels + 1, "MultiGrid" );
				memcpy( info->nDownIts, nDownIts, sizeof(unsigned) * (nLevels + 1) );
				info->nCycles = Memory_Alloc_Array( unsigned, nLevels + 1, "MultiGrid" );
				memcpy( info->nCycles, nCycles, sizeof(unsigned) * (nLevels + 1) );
				info->nFinalIts = nFinalIts;
				info->downSmoothers = downSmoothers;
				info->upSmoothers = upSmoothers;
				info->smoothers = NULL;
				info->rOps = NULL;
				info->pOps = NULL;
				info->iOps = NULL;
				info->nLocalRows = NULL;
				info->nLocalCols = NULL;
				info->rhsVecs = NULL;
				info->xVecs = NULL;
				info->rVecs = NULL;
				
				/* Attempt to locate an existing mapping, if non-existant create a new one. */
				mapping = NULL;
				for( map_i = 0; map_i < mgCtx->nMappings; map_i++ ) {
					if( mgCtx->mappings[map_i].feVar == sm->columnVariable ) {
						mapping = mgCtx->mappings + map_i;
						break;
					}
				}
				if( map_i == mgCtx->nMappings ) {
					mgCtx->mappings = Memory_Realloc_Array( mgCtx->mappings, MGGridMapping, mgCtx->nMappings + 1 );
					mapping = mgCtx->mappings + mgCtx->nMappings;
					mgCtx->nMappings++;
					
					mapping->feVar = sm->columnVariable;
					mapping->maxLevels = 0;
					mapping->maps = NULL;
				}
				
				/* Ensure the required number of levels have been built. */
				MultiGrid_BuildGridMappings( mapping, nLevels );
				
				/* Store the mapping. */
				info->mapping = mapping;
			}
			
			/* Create a new entry in the SLE's handle array. */
			sle->mgHandles = Memory_Realloc_Array( sle->mgHandles, unsigned, sle->nMGHandles + 1 );
			sle->mgHandles[sle->nMGHandles] = handle;
			sle->nMGHandles++;
		}
		
		/* Inform the SLE that multi-grid is to be used. */
		SystemLinearEquations_MG_Enable( sle );
	}
	else {
		fprintf( stderr, "Warning: There are no valid stiffness matrices on a SLE selected for multi-grid.\n" );
	}
	
	/* Free arrays. */
	FreeArray( sms );
}


void MultiGrid_InitMatrixSolver( unsigned handle, MatrixSolver* solver ) {
	MGInfo*		info;
	unsigned	nProcs;
	unsigned	level_i, level_j;

	void MultiGrid_SetSmootherType( MatrixSolver* smoother, const char* type, unsigned nProcs );
	
	
	assert( mgCtx );
	assert( handle < mgCtx->nInfos );
	
	fprintf(stderr,"*** MG: Initialize matrix solvers\n");
	
	
	info = mgCtx->infos + handle;
	nProcs = info->stiffMat->rowVariable->feMesh->layout->decomp->procsInUse;
	
	
	/*
	** Enable multi-grid on the matrix solver and setup all persistant values/options.
	*/

	MatrixSolver_SetKSP_Type( solver, "richardson" );
	MatrixSolver_SetPC_Type( solver, "mg" );
	MatrixSolver_MG_SetLevels( solver, info->nLevels + 1 );
	
	/* Set the iterations for each level. */
	for( level_i = 0, level_j = info->nLevels; level_i <= info->nLevels; level_i++, level_j-- ) {
		MatrixSolver*	down;
		MatrixSolver*	up;

		MatrixSolver_MG_GetDownSmoother( solver, level_j, &down );
		MatrixSolver_SetMaxIts( down, info->nDownIts[level_i] );
		MultiGrid_SetSmootherType( down, info->downSmoothers[level_i], nProcs );

		/* Setup the up smoother. */
		if( level_j > 0 ) {
			MatrixSolver_MG_GetUpSmoother( solver, level_j, &up );
			MatrixSolver_SetMaxIts( up, info->nDownIts[level_i] );
			MultiGrid_SetSmootherType( up, info->upSmoothers[level_i], nProcs );
		}

		MatrixSolver_MG_SetCycles( solver, level_j, info->nCycles[level_i] );
	}
}


void MultiGrid_UpdateMatrixSolver( unsigned handle, MatrixSolver* solver, SLE_Solver* sleSolver ) {
	unsigned	level_i, level_j;
	unsigned	nProcs;
	MGInfo*		info;
	
	assert( mgCtx );
	assert( handle < mgCtx->nInfos );
	
	info = mgCtx->infos + handle;
	nProcs = info->stiffMat->rowVariable->feMesh->layout->decomp->procsInUse;

	fprintf(stderr,"*** MG: Update matrix solvers\n");

	/*
	** Set the required options in the matrix solver for this solve step.
	*/
	
	for( level_i = 0, level_j = info->nLevels; level_i < info->nLevels + 1; level_i++, level_j-- ) {
		MatrixSolver*	down;
		MatrixSolver*	up;

		/* Set this level's smoother.  Remember to skip the coarse grid; it only has a down smoother. */
		MatrixSolver_MG_GetDownSmoother( solver, level_j, &down );
		if( level_j > 0 ) {
			MatrixSolver_MG_GetUpSmoother( solver, level_j, &up );
		}
		
		if( level_i > 0 ) {
			/* For all but the finest. */
			fprintf(stderr,"*** MG: Setting down smoother for level %d to %s \n",level_i,info->downSmoothers[level_i]);
			SLE_Solver_MG_SetupSmoother( sleSolver, down, info->smoothers[level_i - 1], level_i, info->nLevels + 1 );  /* !! i-1 */
			if( level_j > 0 ) {
				fprintf(stderr,"*** MG: Setting up smoother for level %d to %s \n",level_i,info->upSmoothers[level_i]);
				SLE_Solver_MG_SetupSmoother( sleSolver, up, info->smoothers[level_i - 1], level_i, info->nLevels + 1 );
			}
		}
		else {
			/* For the finest, send in the original stiffness matrix. */
			SLE_Solver_MG_SetupSmoother( sleSolver, down, info->stiffMat->matrix, level_i, info->nLevels + 1 );
			if( level_j > 0 ) {
				SLE_Solver_MG_SetupSmoother( sleSolver, up, info->stiffMat->matrix, level_i, info->nLevels + 1 );
			}
		}

		/* Setup the sub blocks if needed. */
		if( nProcs > 1 && level_j > 0 ) {
			MatrixSolver**	blocks;
			unsigned	nBlocks;
			unsigned	block_i;

			if( !strcmp( info->downSmoothers[level_i], "sor" ) ) {
				MatrixSolver_SetupKSP( down );
				MatrixSolver_PC_GetBJacobiSubBlocks( down, &blocks, &nBlocks );
				for( block_i = 0; block_i < nBlocks; block_i++ ) {
					MatrixSolver_SetMaxIts( blocks[block_i], info->nDownIts[level_i] );
					MatrixSolver_SetKSP_Type( blocks[block_i], "richardson" );
					MatrixSolver_SetPC_Type( blocks[block_i], "sor" );
				}
			}

			if( !strcmp( info->upSmoothers[level_i], "sor" ) ) {
				MatrixSolver_SetupKSP( up );
				MatrixSolver_PC_GetBJacobiSubBlocks( up, &blocks, &nBlocks );
				for( block_i = 0; block_i < nBlocks; block_i++ ) {
					MatrixSolver_SetMaxIts( blocks[block_i], info->nDownIts[level_i] );
					MatrixSolver_SetKSP_Type( blocks[block_i], "richardson" );
					MatrixSolver_SetPC_Type( blocks[block_i], "sor" );
				}
			}
		}
		
		/* Set this level's operators aswell as the residual function, except for the coarse level. */
		if( level_i < info->nLevels ) {
			MatrixSolver_MG_SetRestrictionOp( solver, level_j, info->rOps[level_i] );
			MatrixSolver_MG_SetInterpolationOp( solver, level_j, info->rOps[level_i] ); /*  ROPS <-> POPS */
		}

		/* Set this level's work vectors. */
		MatrixSolver_MG_SetWorkSpace( solver, level_j, 
					      info->rhsVecs[level_i], info->xVecs[level_i], info->rVecs[level_i] );

		/* Set this level's matrix, used for extracting residuals. */
		MatrixSolver_MG_SetLevelMatrix( solver, level_j, 
						(level_i > 0) ? info->smoothers[level_i - 1] : info->stiffMat->matrix );
	}
}


void MultiGrid_SetSmootherType( MatrixSolver* smoother, const char* type, unsigned nProcs ) {
	assert( smoother );
	assert( type );

	if( !strcmp( type, "sor" ) ) {
		MatrixSolver_SetKSP_Type( smoother, "richardson" );
		if( nProcs == 1 ) {
			MatrixSolver_SetPC_Type( smoother, "sor" );
		}
		else {
			MatrixSolver_SetPC_Type( smoother, "bjacobi" );
		}
	}
	else if( !strcmp( type, "jacobi" ) ) {
		MatrixSolver_SetKSP_Type( smoother, "richardson" );
		MatrixSolver_SetPC_Type( smoother, "jacobi" );
	}
	else if( !strcmp( type, "bjacobi" ) ) {
		MatrixSolver_SetKSP_Type( smoother, "richardson" );
		MatrixSolver_SetPC_Type( smoother, "bjacobi" );
	}
	else if( !strcmp( type, "gmres" ) ) {
		MatrixSolver_SetKSP_Type( smoother, "gmres" );
		MatrixSolver_SetPC_Type( smoother, "ilu" );  /* Which overrides the "default" LU on the coarse grid of one processor */
		if( nProcs > 1 ) {
			MatrixSolver_SetPC_Type( smoother, "bjacobi" );
		}
	}
	else if( !strcmp( type, "default" ) ) {
		/* Leave the smoother as whatever PETSc wants. */
	}
	else {
		fprintf( stderr, "*** MG: Unsupported smoother type: %s\n", type );
		exit( 1 );
	}
}


void MultiGrid_BuildGridMappings( MGGridMapping* mapping, unsigned nLevels ) {
	double			sTime;
	FiniteElement_Mesh*	mesh;
	MGMapping**		dstMaps;
	MGCoarsener*		coarsener;
	
	void MultiGrid_BuildLMTable( FeVariable* feVar, MGMapping* map );
	
	
	/*
	** Validate input and status.
	*/
	
	/* If the mapping already has enough levels built then return. */
	if( mapping->maxLevels >= nLevels ) {
		return;
	}
	
	assert( mapping->feVar );
	
	
	/*
	** Begin timing.
	*/
	
	fprintf( stderr, "*** MG: Begin building grid mappings...\n" );
	sTime = MPI_Wtime();
	
	
	/*
	** Using either the mesh hint or the mesh itself, determine the most appropriate method to coarsen the
	** existing mesh.
	**
	** NOTE: For now we are assuming the mesh is regular Cartesian.
	*/
	
	mesh = mapping->feVar->feMesh;
	dstMaps = &mapping->maps;
	coarsener = (MGCoarsener*)MGCoarsener_RegCartesian_New();
	MGCoarsener_SetMesh( coarsener, mesh );
	
	
	/*
	** Loop over the mesh map lists, generating each level's mapping as we go.
	*/
	
	{
		unsigned	level_i;
		
		/* Allocate the mesh maps and LM tables. */
		*dstMaps = Memory_Realloc_Array( *dstMaps, MGMapping, nLevels );
		
		/* Generate a mapping for each level of this mesh. */
		for( level_i = mapping->maxLevels; level_i < nLevels; level_i++ ) {
			MGMapping*	dstMap = (*dstMaps) + level_i;
			
			/* Use the coarsener to generate the top-level mapping. */
			MGCoarsener_Coarsen( coarsener, level_i, dstMap );
			
			/* Build this level's LM table based on the top-level. */
			MultiGrid_BuildLMTable( mapping->feVar, dstMap );
		}
	}
	
	/* Set the maximum levels generated. */
	mapping->maxLevels = nLevels;
	
	
	/*
	** End timing.
	*/
	
	fprintf( stderr, "*** MG: Ended, elapsed time is %g seconds.\n", MPI_Wtime() - sTime );
}


void MultiGrid_BuildLMTable( FeVariable* feVar, MGMapping* map ) {
	const int		tag_offset = 2222;
	unsigned		nNodes;
	unsigned*		nLDofs;
	unsigned**	      	lmTbl;
	unsigned		nLUCDofs;
	unsigned		gOffs = 0;
	unsigned		allDofs = 0;
	unsigned*		tuples;
	DofLayout*		dofLayout = feVar->dofLayout;
	FeEquationNumber*	eqNum = feVar->eqNum;
	MPI_Comm		comm = feVar->feMesh->layout->decomp->communicator;
	unsigned    		nProcsInUse = feVar->feMesh->layout->decomp->procsInUse;
	unsigned		rank = feVar->feMesh->layout->decomp->rank;
	unsigned		nOwnedEqNums;
	unsigned*		ownedEqNums;
	unsigned*		reqEqNums;
	unsigned		nReqEqNums;
	unsigned		node_i;

	int unsigned_cmp( const void* opA, const void* opB );
	int cmpEqNum( const void* opA, const void* opB );
	
	
	/*
	** Validate parameters and state.
	*/
	
	assert( feVar && map );


	/*
	** If this proc is not the first (proc 0) then receive from the previous proc a list of equation numbers
	** shared by some procs that are to be ignored from now on.  Also receive the global offset for this proc.
	*/

	if( rank > 0 ) {
		MPI_Status	status;

		/* Receive offset. */
		MPI_Recv( &gOffs, 1, MPI_UNSIGNED, rank - 1, tag_offset, comm, &status );
	}
	
	
	/*
	** Build a new LM table for a coarser level based on the finest level.  This involves a quick sort, will need to
	** be careful else could take forever.
	*/

	/* Build the local mapping. */
	nNodes = map->nNodes;
	nLDofs = Memory_Alloc_Array( unsigned, nNodes, "MultiGrid" );
	lmTbl = Memory_Alloc_Array( unsigned*, nNodes, "MultiGrid" );
	nLUCDofs = 0;
	ownedEqNums = NULL;
	nOwnedEqNums = 0;
	reqEqNums = NULL;
	nReqEqNums = 0;
	for( node_i = 0; node_i < nNodes; node_i++ ) {
		unsigned	tNodeInd = map->nodesTop[node_i];
		unsigned	dof_i;
		
		/* How many dofs on this node? */
		nLDofs[node_i] = dofLayout->dofCounts[tNodeInd];
		
		/* Add this node's DOF entries to the map. */
		lmTbl[node_i] = Memory_Alloc_Array( unsigned, nLDofs[node_i], "MultiGrid" );
		for( dof_i = 0; dof_i < nLDofs[node_i]; dof_i++ ) {
			/* We need to count the range this proc owns. */
			lmTbl[node_i][dof_i] = eqNum->destinationArray[tNodeInd][dof_i];
			if( lmTbl[node_i][dof_i] >= eqNum->firstOwnedEqNum && 
			    lmTbl[node_i][dof_i] <= eqNum->lastOwnedEqNum )
			{
				allDofs++;
				nLUCDofs++;

				/* Store this as a local. */
				ownedEqNums = Memory_Realloc_Array( ownedEqNums, unsigned, nOwnedEqNums + 1 );
				ownedEqNums[nOwnedEqNums] = lmTbl[node_i][dof_i];
				nOwnedEqNums++;
			}
			else if( lmTbl[node_i][dof_i] != -1 ) {
				/* We'll need to collect the equation number for this guy from another proc. */
				reqEqNums = Memory_Realloc_Array( reqEqNums, unsigned, nReqEqNums + 1 );
				reqEqNums[nReqEqNums] = lmTbl[node_i][dof_i];
				nReqEqNums++;

				allDofs++;
				lmTbl[node_i][dof_i] = -1;
			}
		}
	}

	/* Flatten the array. */
	tuples = Memory_Alloc_Array( unsigned, nLUCDofs * 3, "MultiGrid" );
	nLUCDofs = 0;
	for( node_i = 0; node_i < nNodes; node_i++ ) {
		unsigned	dof_i;
		
		/* Add this node's DOF entries to the map. */
		for( dof_i = 0; dof_i < nLDofs[node_i]; dof_i++ ) {
			if( lmTbl[node_i][dof_i] != -1 ) {
				tuples[nLUCDofs * 3 + 0] = lmTbl[node_i][dof_i];
				tuples[nLUCDofs * 3 + 1] = node_i;
				tuples[nLUCDofs * 3 + 2] = dof_i;
				nLUCDofs++;
			}
		}
	}

	/* Sort the array based on equation number. */
	qsort( tuples, nLUCDofs, sizeof(unsigned) * 3, cmpEqNum );

	/* Remove any holes in the equation numbers and reassemble. */
	for( node_i = 0; node_i < nLUCDofs; node_i++ ) {
		unsigned	ucDofInd = node_i * 3;

		lmTbl[tuples[ucDofInd + 1]][tuples[ucDofInd + 2]] = gOffs + node_i;

		/* We'll use the tuples array later. */
		tuples[ucDofInd] = gOffs + node_i;
	}


	/*
	** Send the next proc's global offset.
	*/

	if( rank < nProcsInUse - 1 ) {
		/* Send offset. */
		gOffs += nLUCDofs;
		MPI_Send( &gOffs, 1, MPI_UNSIGNED, rank + 1, tag_offset, comm );
	}


	/*
	** Now, we need to communicate the equation numbers stored on other procs that we need here.
	*/

	{
		Sync*		sync;

		/* First of all, sort the arrays we have, as 'tuples' is sorted by global eqNum. */
		qsort( ownedEqNums, nOwnedEqNums, sizeof(unsigned), unsigned_cmp );
		qsort( reqEqNums, nReqEqNums, sizeof(unsigned), unsigned_cmp );

		/* We'll need to expand the tuples array to accomodate the incoming values. */
		tuples = Memory_Realloc_Array( tuples, unsigned, (nLUCDofs + nReqEqNums) * 3 );

		/* Build the sync class, initialise. */
		sync = Sync_New( "Multigrid_LMSync" );
		Sync_Negotiate( sync, 
				eqNum->globalSumUnconstrainedDofs, 
				ownedEqNums, nOwnedEqNums, 
				NULL, 0, 
				reqEqNums, nReqEqNums, 
				comm );
		Sync_SetDomainArray( sync, 
				     sizeof(unsigned) * 3, 
				     sizeof(unsigned) * 3, tuples );

		/* Send and recieve. */
		Sync_SendRecv( sync );

#if 0
		/* Build the sync class. */
		sync = Sync_New( "Multigrid_LMSync" );
		Sync_Negotiate( sync, eqNum->globalSumUnconstrainedDofs, nOwnedEqNums, ownedEqNums, nReqEqNums, reqEqNums, comm );

		{
			unsigned	proc_i;
			unsigned	eqNum_i;
			unsigned	snkInd = 0;

			/* Build the sink and source offsets arrays. */
			snkOffs = Memory_Alloc_Array( unsigned, sync->netSink, "MultiGrid" );
			for( proc_i = 0; proc_i < sync->nProcs; proc_i++ ) {
				for( item_i = 0; item_i < sync->nSink[proc_i]; item_i++ ) {
					for( eqNum_i = 0; eqNum_i < nOwnedEqNums; eqNum_i++ ) {
						if( ownedEqNums[eqNum_i] == sync->sink[proc_i][item_i] ) {
							snkOffs[snkInd++] = eqNum_i;
						}
					}
				}
			}
		}

		srcOffs = Memory_Alloc_Array( unsigned, nReqEqNums, "MultiGrid" );
		for( item_i = 0; item_i < nReqEqNums; item_i++ ) {
			srcOffs[item_i] = item_i;
		}

		/* Initialise the send/recv process. */
		Sync_SendRecvInitialise( sync, sizeof(unsigned), snkOffs, sizeof(unsigned) * 3, srcOffs, sizeof(unsigned) );

		/* Build the destination array and send/recv. */
		ghostEqNums = Memory_Alloc_Array( unsigned, nReqEqNums, "MultiGrid" );
		Sync_SendRecv( sync, tuples, ghostEqNums );

		/* Free some stuff. */
		FreeArray( reqEqNums );
		FreeArray( ownedEqNums );
#endif

		/* Update the LM table. */
		for( node_i = 0; node_i < nNodes; node_i++ ) {
			unsigned	tNodeInd = map->nodesTop[node_i];
			unsigned	dof_i;

			for( dof_i = 0; dof_i < nLDofs[node_i]; dof_i++ ) {
				if( !(eqNum->destinationArray[tNodeInd][dof_i] >= eqNum->firstOwnedEqNum && 
				      eqNum->destinationArray[tNodeInd][dof_i] <= eqNum->lastOwnedEqNum) &&
				    eqNum->destinationArray[tNodeInd][dof_i] != -1 )
				{
					unsigned	gEqNum = eqNum->destinationArray[tNodeInd][dof_i];
					unsigned	dInd;

					dInd = Sync_MapGlobal( sync, gEqNum );
					assert( dInd >= nLUCDofs );
					assert( dInd < (nLUCDofs + nReqEqNums) );

					lmTbl[node_i][dof_i] = tuples[dInd * 3];

#if 0
					/* Hunt down the global eq num. */
					for( rEqNum_i = 0; rEqNum_i < nReqEqNums; rEqNum_i++ ) {
						if( reqEqNums[rEqNum_i] == gEqNum ) {
							lmTbl[node_i][dof_i] = ghostEqNums[rEqNum_i];
							break;
						}
					}
					assert( rEqNum_i < nReqEqNums );
#endif
				}
			}
		}
		
		/* Free more stuff. */
		FreeArray( tuples );
		Stg_Class_Delete( sync );
	}

	
	/* Store results. */
	map->nDofs = nLDofs;
	map->nUCDofs = allDofs;
	map->nOwnedEqNums = nLUCDofs;
	map->lmTable = lmTbl;
}


int unsigned_cmp( const void* opA, const void* opB ) {
	if( *(unsigned*)opA > *(unsigned*)opB ) {
		return 1;
	}
	else if( *(unsigned*)opA < *(unsigned*)opB ) {
		return -1;
	}
	else {
		return 0;
	}
}


int cmpEqNum( const void* opA, const void* opB ) {
	if( *((unsigned*)opA) > *((unsigned*)opB) ) {
		return 1;
	}
	else if( *((unsigned*)opA) < *((unsigned*)opB) ) {
		return -1;
	}
	else {
		return 0;
	}
}


void MultiGrid_BuildGridOps( unsigned handle ) {
	double	sTime;
	unsigned	level_i;
	MGInfo*	info;

	void MultiGrid_BuildProjectionOp( FeVariable* feVar, 
					  MGMapping* cMap, MGMapping* fMap, 
					  Matrix** mat, 
					  unsigned* nLocalRows, unsigned* nLocalCols );
	void MultiGrid_BuildRestrictionOp( FeVariable* feVar, 
					   MGMapping* cMap, MGMapping* fMap, 
					   Matrix** mat, 
					   unsigned* nLocalRows, unsigned* nLocalCols, 
					   Bool avgBCs );
#if 0					
	void MultiGrid_BuildInjectionOp( FeVariable* feVar, 
					 MGMapping* cMap, MGMapping* fMap, 
					 Matrix** mat, 
					 unsigned* nLocalRows, unsigned* nLocalCols );
#endif

	void MultiGrid_NormaliseOp( Matrix* op );
	
	
	assert( mgCtx );
	assert( handle < mgCtx->nInfos );
	
	info = mgCtx->infos + handle;
	
	
	/*
	** Start timing.
	*/
	
	fprintf( stderr, "*** MG: Begin building inter-grid operators...\n" );
	sTime = MPI_Wtime();
	
	
	/*
	** Validate status and input.
	*/
	
	assert( info->stiffMat );
	assert( info->mapping );

	/*
	** Build/update the restriction and interpolation operators.
	*/
	
	/* Allocate space for matrices. */
	if( !info->rOps ) {
		fprintf( stderr, "*** MG: Allocating restriction operator space\n");
		info->rOps = Memory_Alloc_Array( Matrix*, info->nLevels, "MultiGrid" );
		memset( info->rOps, 0, sizeof(Matrix*) * info->nLevels );
	}
	
#if 0	
	if( !info->pOps ) {
		fprintf( stderr, "*** MG: Allocating projection operator space\n");
		info->pOps = Memory_Alloc_Array( Matrix*, info->nLevels, "MultiGrid" );
		memset( info->pOps, 0, sizeof(Matrix*) * info->nLevels );
	}

	if( !info->iOps ) {
		info->iOps = Memory_Alloc_Array( Matrix*, info->nLevels, "MultiGrid" );
		memset( info->iOps, 0, sizeof(Matrix*) * info->nLevels );
	}
#endif
	
	/* Allocate space for the row/col counts. */
	if( !info->nLocalRows ) {
		info->nLocalRows = Memory_Alloc_Array( unsigned, info->nLevels, "MultiGrid" );
		memset( info->nLocalRows, 0, sizeof(unsigned) * info->nLevels );
	}
	if( !info->nLocalCols ) {
		info->nLocalCols = Memory_Alloc_Array( unsigned, info->nLevels, "MultiGrid" );
		memset( info->nLocalCols, 0, sizeof(unsigned) * info->nLevels );
	}
	
	/* Build an operator for each level. */
	for( level_i = 0; level_i < info->nLevels; level_i++ ) {
		/* If we already have operators in there, clear them out. */
		if( info->rOps[level_i] ) {
			Matrix_Destroy( info->rOps[level_i] );
		}
		
#if 0				
		if( info->pOps[level_i] && info->pOps[level_i] != info->rOps[level_i] ) {
			Matrix_Destroy( info->pOps[level_i] );
		}
		if( info->iOps[level_i] ) {
			Matrix_Destroy( info->iOps[level_i] );
		}
#endif
		
		fprintf( stderr, "*** MG: Building inter-grid operators for level %d\n", level_i );
		
		MultiGrid_BuildRestrictionOp( info->mapping->feVar, 
					      &info->mapping->maps[level_i], (level_i > 0) ? &info->mapping->maps[level_i - 1] : NULL, 
					      &info->rOps[level_i], 
					      &info->nLocalRows[level_i], &info->nLocalCols[level_i], 
					      False );
					
					
					
#if 0
		MultiGrid_BuildProjectionOp( info->mapping->feVar, 
					     &info->mapping->maps[level_i], (level_i > 0) ? &info->mapping->maps[level_i - 1] : NULL, 
					     &info->pOps[level_i], 
					     &info->nLocalRows[level_i], &info->nLocalCols[level_i] );
			
		MultiGrid_BuildInjectionOp( info->mapping->feVar, 
					    &info->mapping->maps[level_i], (level_i > 0) ? &info->mapping->maps[level_i - 1] : NULL, 
					    &info->iOps[level_i], 
					    &info->nLocalRows[level_i], &info->nLocalCols[level_i] );
#endif

		fprintf( stderr, "*** MG: Building inter-grid operators for level %d - completed\n", level_i );
	}
	
	
	/*
	** End timing.
	*/
	
	fprintf( stderr, "*** MG: Ended, elapsed time is %g seconds.\n", MPI_Wtime() - sTime );
}


void MultiGrid_BuildSmoothers( unsigned handle ) {
	double	sTime;
	unsigned	level_i;
	MGInfo*	info;
	
	
	assert( mgCtx );
	assert( handle < mgCtx->nInfos );
	
	info = mgCtx->infos + handle;


	/*
	** Start timing.
	*/
	
	fprintf( stderr, "*** MG: Begin building smoother operators...\n" );
	sTime = MPI_Wtime();
	
	
	/*
	** Validate status and input.
	*/
	
	assert( info );
	assert( info->stiffMat );
	assert( info->rOps );
	/*  assert( info->pOps ); */
	
	
	/*
	** Use the restriction operators to build a set of smoothers from the base stiffness matrix.
	*/
	
	/* Allocate space if needed. */
	if( !info->smoothers ) {
		info->smoothers = Memory_Alloc_Array( Matrix*, info->nLevels, "MultiGrid" );
		memset( info->smoothers, 0, sizeof(Matrix*) * info->nLevels );
	}
	
	for( level_i = 0; level_i < info->nLevels; level_i++ ) {
		/* Zero the smoother before filling. */
		if( info->smoothers[level_i] ) {
			Matrix_Destroy( info->smoothers[level_i] );
			info->smoothers[level_i] = NULL;
		}
		assert( info->rOps[level_i] );
		if( level_i == 0 ) {
			Matrix_PtAP( info->stiffMat->matrix, info->rOps[level_i], &info->smoothers[level_i], 0.0 );
		}
		else {
			Matrix_PtAP( info->smoothers[level_i - 1], info->rOps[level_i], &info->smoothers[level_i], 0.0 );
		}
	}
	
	
	/*
	** End timing.
	*/
	
	fprintf( stderr, "*** MG: Ended, elapsed time is %g seconds.\n", MPI_Wtime() - sTime );
}


void MultiGrid_BuildProjectionOp( FeVariable* feVar, 
				  MGMapping* cMap, MGMapping* fMap, 
				  Matrix** mat, 
				  unsigned* nLocalRows, unsigned* nLocalCols )
{
	FiniteElement_Mesh*	mesh = feVar->feMesh;
	MPI_Comm		comm = mesh->layout->decomp->communicator;
	unsigned		nRows;
	unsigned		nCols;
	unsigned		nzs;
	unsigned		maxNZs;
	unsigned*		rowNZs;
	unsigned**		indices;
	double**		entries;
	unsigned*		rowIndices;
	
	
	/*
	** Validate.
	*/
	
	assert( feVar );
	assert( cMap );
	assert( mat );
	assert( nLocalRows && nLocalCols );
	
	
	/*
	** Build the restriction operator in local memory first.  We need to do this so we can calculate the real number of non-zeros
	** there will be.
	*/
	
	{
		ElementLayout*	elLayout = mesh->layout->elementLayout;
		unsigned		nCNodes = cMap->nNodes;
		unsigned*		nRowDofs;
		unsigned*		nColDofs;
		unsigned**	rowEqNums;
		unsigned**	colEqNums;
		unsigned		curRow;
		unsigned		cNode_i;
		
		fprintf( stderr, "*** MG: Build projectiion operators\n");
		
		
		/* Setup LM table info. */
		if( fMap ) {
			nCols = fMap->nUCDofs;
		}
		else {
			unsigned		node_i;
			unsigned		dof_i;
			FeEquationNumber*	eqNum = feVar->eqNum;
			DofLayout*		dofLayout = feVar->dofLayout;

			nCols = 0;
			for( node_i = 0; node_i < mesh->nodeLocalCount; node_i++ ) {
				for( dof_i = 0; dof_i < dofLayout->dofCounts[node_i]; dof_i++ ) {
					if( eqNum->destinationArray[node_i][dof_i] != -1 ) {
						nCols++;
					}
				}
			}
		}
		nColDofs = fMap ? fMap->nDofs : feVar->dofLayout->dofCounts;
		colEqNums = fMap ? fMap->lmTable : (unsigned**)feVar->eqNum->destinationArray;
		nRows = cMap->nUCDofs;
		nRowDofs = cMap->nDofs;
		rowEqNums = cMap->lmTable;
		

		
		/* Allocate initial space for the matrix entries. */
		nzs = 0;
		maxNZs = 0;
		curRow = 0;
		rowNZs = Memory_Alloc_Array( unsigned, nRows, "MultiGrid" );
		memset( rowNZs, 0, sizeof(unsigned) * nRows );
		indices = Memory_Alloc_Array( unsigned*, nRows, "MultiGrid" );
		memset( indices, 0, sizeof(unsigned*) * nRows );
		entries = Memory_Alloc_Array( double*, nRows, "MultiGrid" );
		memset( entries, 0, sizeof(double*) * nRows );
		rowIndices = Memory_Alloc_Array( unsigned, nRows, "MultiGrid" );
		memset( indices, 0, sizeof(unsigned) * nRows );
		
		
		/* Loop over the nodes in the coarse mesh. */
		for( cNode_i = 0; cNode_i < nCNodes; cNode_i++ ) {
			unsigned	rDof_i;
			
			/* Loop over the DOFs on this coarse node (row DOFs). */
			for( rDof_i = 0; rDof_i < nRowDofs[cNode_i]; rDof_i++ ) {
				unsigned	rDofInd = rowEqNums[cNode_i][rDof_i];
				unsigned	cEl_i;
				
				if( rDofInd == -1 ) {
					continue;
				}
				
				/* Loop over the coarse elements incident on this coarse node. */
				for( cEl_i = 0; cEl_i < cMap->nIncElsLocal[cNode_i]; cEl_i++ ) {
					unsigned		cElInd = cMap->incElsLocal[cNode_i][cEl_i];
					unsigned		elNodeInd = cMap->elNodeInds[cNode_i][cEl_i];
					ElementType*	elType;
					double*		sfVals;
					Coord**		tNodeCoordPtrs;
					Coord		lCoord;
					unsigned		cNode_j;
					unsigned		pNode_i;
					
					/* Get this element's type and allocate space for the nodal shape function values. It is assumed that the 
					   finest mesh has the same element types for all elements covered by this coarse element. */
					elType = FiniteElement_Mesh_ElementTypeAt( mesh, cMap->incElTop[cElInd] );
					sfVals = Memory_Alloc_Array( double, elType->nodeCount, "MultiGrid" );
					
					/* For the time being we require there to be the same number of elType->nodeCount as incident nodes. */
					assert( elType->nodeCount == cMap->nIncNodesLocal[cElInd] );
					
					/* Collect a set of top-level node coordinate pointers. */
					tNodeCoordPtrs = Memory_Alloc_Array( Coord*, cMap->nIncNodesLocal[cElInd], "MultiGrid" );
					for( cNode_j = 0; cNode_j < cMap->nIncNodesLocal[cElInd]; cNode_j++ ) {
						unsigned	tNodeInd = cMap->nodesTop[cMap->incNodesLocal[cElInd][cNode_j]];
						
						tNodeCoordPtrs[cNode_j] = &mesh->nodeCoord[tNodeInd];
					}
					
					/* Loop over the previous-level nodes incident on this coarse element. */
					for( pNode_i = 0; pNode_i < cMap->nIncNodesPrev[cElInd]; pNode_i++ ) {
						unsigned	pNodeInd = cMap->incNodesPrev[cElInd][pNode_i];
						unsigned	tNodeInd = fMap ? fMap->nodesTop[pNodeInd] : pNodeInd;
						unsigned	cDof_i;
						
						/* Calculate the shape function values of this top-node with respect to the coarse node on this
						   element. */
						ElementType_ConvertGlobalCoordToElLocal( elType, 
														 elLayout, 
														 (const Coord**)tNodeCoordPtrs, 
														 mesh->nodeCoord[tNodeInd], 
														 lCoord );
						ElementType_EvaluateShapeFunctionsAt( elType, lCoord, sfVals );
						
						if( sfVals[elNodeInd] != 0.0 ) {
							/* Loop over the DOFs on the previous node (column DOFs). */
							for( cDof_i = 0; cDof_i < nColDofs[pNodeInd]; cDof_i++ ) {
								unsigned	cDofInd = colEqNums[pNodeInd][cDof_i];
								
								if( cDofInd == -1 || cDof_i != rDof_i ) {
									continue;
								}
								
								/* Expand the row arrays and store value. */
								indices[curRow] = Memory_Realloc_Array( indices[curRow], unsigned, rowNZs[curRow] + 1 );
								indices[curRow][rowNZs[curRow]] = cDofInd;
								entries[curRow] = Memory_Realloc_Array( entries[curRow], double, rowNZs[curRow] + 1 );
								entries[curRow][rowNZs[curRow]] = sfVals[elNodeInd];
								rowNZs[curRow]++;
							}
						}
					}
					
					/* Free some arrays. */
					FreeArray( sfVals );
					FreeArray( tNodeCoordPtrs );
				}
				
				/* Update the net number of non-zeros and max non-zeros. */
				nzs += rowNZs[curRow];
				maxNZs = (rowNZs[curRow] > maxNZs) ? rowNZs[curRow] : maxNZs;

				/* Update the current row. */
				rowIndices[curRow] = rDofInd;
				curRow++;
			}
		}
	}
	
	
	/*
	** Transfer the operator entries to a matrix.
	*/
	
	{
		unsigned	row_i;
		Matrix*		rOp;
		unsigned	nMatRows = fMap ? fMap->nOwnedEqNums : feVar->eqNum->localEqNumsOwnedCount;
		unsigned	nMatCols = cMap->nOwnedEqNums;

		/* Stash some info. */
		*nLocalRows = nMatRows;
		*nLocalCols = nMatCols;
		
		/* Allocate the matrix.  As it turns out, the non-zeros this function refers to is on a row by row basis.  Luckily, 
		   Galerkin style operators should have about the same number of non-zeros per row. */
		
		fprintf( stderr, "*** MG: Build projection matrix %d x %d (maxNZ = %d)\n",nMatRows,nMatCols,maxNZs);
		
		rOp = Matrix_New( comm, nMatRows, nMatCols, maxNZs );
		
		/* Transfer, row by row. */
		for( row_i = 0; row_i < nRows; row_i++ ) {
			if( rowNZs[row_i] ) {
				Matrix_Insert( rOp, rowNZs[row_i], indices[row_i], 1, rowIndices + row_i, entries[row_i] );
			}
		}
		
		/* Assemble the matrix. */
		Matrix_AssemblyBegin( rOp );
		Matrix_AssemblyEnd( rOp );

		/* Save. */
		*mat = rOp;
	}
	
	/* Free arrays. */
	FreeArray( rowNZs );
	FreeArray2D( nRows, indices );
	FreeArray2D( nRows, entries );
	FreeArray( rowIndices );
	
	fprintf( stderr, "*** MG: Build projectiion operators - done\n");
	
	
}


void MultiGrid_BuildRestrictionOp( FeVariable* feVar, 
				   MGMapping* cMap, MGMapping* fMap, 
				   Matrix** mat, 
				   unsigned* nLocalRows, unsigned* nLocalCols, 
				   Bool avgBCs )
{
	FiniteElement_Mesh*	mesh = feVar->feMesh;
	MPI_Comm		comm = mesh->layout->decomp->communicator;
	unsigned		nRows;
	unsigned		nCols;
	unsigned		nzs;
	unsigned		maxNZs;
	unsigned*		rowNZs;
	unsigned**		indices;
	double**		entries;
	unsigned*		rowIndices;
	
	
	/*
	** Validate.
	*/
	
	assert( feVar );
	assert( cMap );
	assert( mat );
	assert( nLocalRows && nLocalCols );
	
	
	/*
	** Build the restriction operator in local memory first.  We need to do this so we can calculate the real number of non-zeros
	** there will be.
	*/
	
	{
		ElementLayout*	elLayout = mesh->layout->elementLayout;
		unsigned	nCNodes = cMap->nNodes;
		unsigned*	nRowDofs;
		unsigned*	nColDofs;
		unsigned**	rowEqNums;
		unsigned**	colEqNums;
		unsigned	curRow;
		IndexSet*	bcSet;
		unsigned	cNode_i;
		
		double vanishinglySmallContribution = 1.0e-8;
		
		fprintf( stderr, "*** MG: Build restriction operators\n");
		
		
		
		/* Setup LM table info. */
		if( fMap ) {
			nCols = fMap->nUCDofs;
		}
		else {
			unsigned		node_i;
			unsigned		dof_i;
			FeEquationNumber*	eqNum = feVar->eqNum;
			DofLayout*		dofLayout = feVar->dofLayout;

			nCols = 0;
			for( node_i = 0; node_i < mesh->nodeLocalCount; node_i++ ) {
				for( dof_i = 0; dof_i < dofLayout->dofCounts[node_i]; dof_i++ ) {
					if( eqNum->destinationArray[node_i][dof_i] != -1 ) {
						nCols++;
					}
				}
			}
		}
		nColDofs = fMap ? fMap->nDofs : feVar->dofLayout->dofCounts;
		colEqNums = fMap ? fMap->lmTable : (unsigned**)feVar->eqNum->destinationArray;
		nRows = cMap->nUCDofs;
		nRowDofs = cMap->nDofs;
		rowEqNums = cMap->lmTable;
		
		/* Allocate initial space for the matrix entries. */
		nzs = 0;
		maxNZs = 0;
		curRow = 0;
		rowNZs = Memory_Alloc_Array( unsigned, nRows, "MultiGrid" );
		memset( rowNZs, 0, sizeof(unsigned) * nRows );
		indices = Memory_Alloc_Array( unsigned*, nRows, "MultiGrid" );
		memset( indices, 0, sizeof(unsigned*) * nRows );
		entries = Memory_Alloc_Array( double*, nRows, "MultiGrid" );
		memset( entries, 0, sizeof(double*) * nRows );
		rowIndices = Memory_Alloc_Array( unsigned, nRows, "MultiGrid" );
		memset( indices, 0, sizeof(unsigned) * nRows );

		/* Hack out the bc index set. */
		bcSet = feVar->eqNum->bcEqNums;

		/* Loop over the nodes in the coarse mesh. */
		for( cNode_i = 0; cNode_i < nCNodes; cNode_i++ ) {
			unsigned	rDof_i;
			
			/* Loop over the DOFs on this coarse node (row DOFs). */
			for( rDof_i = 0; rDof_i < nRowDofs[cNode_i]; rDof_i++ ) {
				unsigned	rDofInd = rowEqNums[cNode_i][rDof_i];
				unsigned	cEl_i;
				
				if( rDofInd == -1 ) {
					continue;
				}

				/* Loop over the coarse elements incident on this coarse node. */
				for( cEl_i = 0; cEl_i < cMap->nIncElsLocal[cNode_i]; cEl_i++ ) {
					unsigned	cElInd = cMap->incElsLocal[cNode_i][cEl_i];
					unsigned	elNodeInd = cMap->elNodeInds[cNode_i][cEl_i];
					ElementType*	elType;
					double*		sfVals;
					Coord**		tNodeCoordPtrs;
					Coord		lCoord;
					unsigned	cNode_j;
					unsigned	pNode_i;
					
					/* Get this element's type and allocate space for the nodal shape function values. It is assumed that the 
					   finest mesh has the same element types for all elements covered by this coarse element. */
					elType = FiniteElement_Mesh_ElementTypeAt( mesh, cMap->incElTop[cElInd] );
					sfVals = Memory_Alloc_Array( double, elType->nodeCount, "MultiGrid" );
					
					/* For the time being we require there to be the same number of elType->nodeCount as incident nodes. */
					assert( elType->nodeCount == cMap->nIncNodesLocal[cElInd] );
					
					/* Collect a set of top-level node coordinate pointers. */
					tNodeCoordPtrs = Memory_Alloc_Array( Coord*, cMap->nIncNodesLocal[cElInd], "MultiGrid" );
					for( cNode_j = 0; cNode_j < cMap->nIncNodesLocal[cElInd]; cNode_j++ ) {
						unsigned	tNodeInd = cMap->nodesTop[cMap->incNodesLocal[cElInd][cNode_j]];
						
						tNodeCoordPtrs[cNode_j] = &mesh->nodeCoord[tNodeInd];
					}
					
					/* Loop over the previous-level nodes incident on this coarse element. */
					for( pNode_i = 0; pNode_i < cMap->nIncNodesPrev[cElInd]; pNode_i++ ) {
						unsigned	pNodeInd = cMap->incNodesPrev[cElInd][pNode_i];
						unsigned	tNodeInd = fMap ? fMap->nodesTop[pNodeInd] : pNodeInd;
						unsigned	cDof_i;

						/* If we aren't looking at a direct mapping and 'avgBCs' is false, 
						   and this row is a bc, skip it. */
						if( !avgBCs && bcSet ) {
							unsigned	tDofInd;

							tDofInd = feVar->eqNum->destinationArray[cMap->nodesTop[cNode_i]][rDof_i];
							if( IndexSet_IsMember( bcSet, tDofInd ) && 
							    tNodeInd != cMap->nodesTop[cNode_i] )
							{
								continue;
							}
						}
						
						/* Calculate the shape function values of this top-node with respect to the coarse node on this
						   element. */
						ElementType_ConvertGlobalCoordToElLocal( elType, 
											 elLayout, 
											 (const Coord**)tNodeCoordPtrs, 
											 mesh->nodeCoord[tNodeInd], 
											 lCoord );
						ElementType_EvaluateShapeFunctionsAt( elType, lCoord, sfVals );
						
						if( sfVals[elNodeInd] > vanishinglySmallContribution ) {
							/* Loop over the DOFs on the previous node (column DOFs). */
							for( cDof_i = 0; cDof_i < nColDofs[pNodeInd]; cDof_i++ ) {
								unsigned	cDofInd = colEqNums[pNodeInd][cDof_i];
								
								if( cDofInd == -1 || cDof_i != rDof_i ) {
									continue;
								}

								/* Have we already covered this index? */
								{
									unsigned	col_i;

									for( col_i = 0; col_i < rowNZs[curRow]; col_i++ ) {
										if( indices[curRow][col_i] == cDofInd ) {
											break;
										}
									}
									if( col_i < rowNZs[curRow] ) {
										continue;
									}
								}
								
								/* Expand the row arrays and store value. */
								indices[curRow] = Memory_Realloc_Array( indices[curRow], unsigned, rowNZs[curRow] + 1 );
								indices[curRow][rowNZs[curRow]] = cDofInd;
								entries[curRow] = Memory_Realloc_Array( entries[curRow], double, rowNZs[curRow] + 1 );
								entries[curRow][rowNZs[curRow]] = sfVals[elNodeInd];
								rowNZs[curRow]++;
							}
						}
					}
					
					/* Free some arrays. */
					FreeArray( sfVals );
					FreeArray( tNodeCoordPtrs );
				}
				
				/* Update the net number of non-zeros and max non-zeros. */
				nzs += rowNZs[curRow];
				maxNZs = (rowNZs[curRow] > maxNZs) ? rowNZs[curRow] : maxNZs;

				/* Update the current row. */
				rowIndices[curRow] = rDofInd;
				curRow++;
			}
		}
	}
	
	
	/*
	** Transfer the operator entries to a matrix.
	*/
	
	{
		unsigned	row_i;
		Matrix*		rOp;
		unsigned	nMatRows = fMap ? fMap->nOwnedEqNums : feVar->eqNum->localEqNumsOwnedCount;
		unsigned	nMatCols = cMap->nOwnedEqNums;

		/* Stash some info. */
		*nLocalRows = nMatCols;
		*nLocalCols = nMatRows;
		
		/* Allocate the matrix.  As it turns out, the non-zeros this function refers to is on a row by row basis.  Luckily, 
		   Galerkin style operators should have about the same number of non-zeros per row. */
		
		
		fprintf( stderr, "*** MG: Build restriction matrix %d x %d (maxNZ = %d)\n",nMatRows,nMatCols,maxNZs);
		rOp = Matrix_New( comm, nMatRows, nMatCols, maxNZs );
		
		/* Transfer, row by row. */
		for( row_i = 0; row_i < nRows; row_i++ ) {
			if( rowNZs[row_i] ) {
				Matrix_Insert( rOp, rowNZs[row_i], indices[row_i], 1, rowIndices + row_i, entries[row_i] );
			}
		}
		
		/* Assemble the matrix. */
		Matrix_AssemblyBegin( rOp );
		Matrix_AssemblyEnd( rOp );

		/* Save. */
		*mat = rOp;
	}
	
	/* Free arrays. */
	FreeArray( rowNZs );
	FreeArray2D( nRows, indices );
	FreeArray2D( nRows, entries );
	FreeArray( rowIndices );
	
	fprintf( stderr, "*** MG: Build restriction operators - Done\n");
	
	
}


void MultiGrid_BuildInjectionOp( FeVariable* feVar, 
				 MGMapping* cMap, MGMapping* fMap, 
				 Matrix** mat, 
				 unsigned* nLocalRows, unsigned* nLocalCols )
{
	FiniteElement_Mesh*	mesh = feVar->feMesh;
	MPI_Comm		comm = mesh->layout->decomp->communicator;
	unsigned		nRows;
	unsigned		nCols;
	unsigned		nzs;
	unsigned		maxNZs;
	unsigned*		rowNZs;
	unsigned**		indices;
	double**		entries;
	unsigned*		rowIndices;
	
	
	/*
	** Validate.
	*/
	
	assert( feVar );
	assert( cMap );
	assert( mat );
	assert( nLocalRows && nLocalCols );
	
	
	/*
	** Build the restriction operator in local memory first.  We need to do this so we can calculate the real number of non-zeros
	** there will be.
	*/
	
	{
		ElementLayout*	elLayout = mesh->layout->elementLayout;
		unsigned	nCNodes = cMap->nNodes;
		unsigned*	nRowDofs;
		unsigned*	nColDofs;
		unsigned**	rowEqNums;
		unsigned**	colEqNums;
		unsigned	curRow;
		unsigned	cNode_i;
		
		/* Setup LM table info. */
		nCols = fMap ? fMap->nUCDofs : feVar->eqNum->localEqNumsOwnedCount;
		nColDofs = fMap ? fMap->nDofs : feVar->dofLayout->dofCounts;
		colEqNums = fMap ? fMap->lmTable : (unsigned**)feVar->eqNum->destinationArray;
		nRows = cMap->nUCDofs;
		nRowDofs = cMap->nDofs;
		rowEqNums = cMap->lmTable;
		
		/* Stash some info. */
		*nLocalRows = nRows;
		*nLocalCols = nCols;
		
		/* Allocate initial space for the matrix entries. */
		nzs = 0;
		maxNZs = 0;
		curRow = 0;
		rowNZs = Memory_Alloc_Array( unsigned, nRows, "MultiGrid" );
		memset( rowNZs, 0, sizeof(unsigned) * nRows );
		indices = Memory_Alloc_Array( unsigned*, nRows, "MultiGrid" );
		memset( indices, 0, sizeof(unsigned*) * nRows );
		entries = Memory_Alloc_Array( double*, nRows, "MultiGrid" );
		memset( entries, 0, sizeof(double*) * nRows );
		rowIndices = Memory_Alloc_Array( unsigned, nRows, "MultiGrid" );
		memset( indices, 0, sizeof(unsigned) * nRows );
		
		/* Loop over the nodes in the coarse mesh. */
		for( cNode_i = 0; cNode_i < nCNodes; cNode_i++ ) {
			unsigned	rDof_i;
			
			/* Loop over the DOFs on this coarse node (row DOFs). */
			for( rDof_i = 0; rDof_i < nRowDofs[cNode_i]; rDof_i++ ) {
				unsigned	rDofInd = rowEqNums[cNode_i][rDof_i];
				unsigned	cEl_i;
				
				if( rDofInd == -1 ) {
					continue;
				}

				/* Loop over the coarse elements incident on this coarse node. */
				for( cEl_i = 0; cEl_i < cMap->nIncElsLocal[cNode_i]; cEl_i++ ) {
					unsigned		cElInd = cMap->incElsLocal[cNode_i][cEl_i];
					unsigned		elNodeInd = cMap->elNodeInds[cNode_i][cEl_i];
					ElementType*	elType;
					double*		sfVals;
					Coord**		tNodeCoordPtrs;
					Coord		lCoord;
					unsigned		cNode_j;
					unsigned		pNode_i;
					
					/* Get this element's type and allocate space for the nodal shape function values. It is assumed that the 
					   finest mesh has the same element types for all elements covered by this coarse element. */
					elType = FiniteElement_Mesh_ElementTypeAt( mesh, cMap->incElTop[cElInd] );
					sfVals = Memory_Alloc_Array( double, elType->nodeCount, "MultiGrid" );
					
					/* For the time being we require there to be the same number of elType->nodeCount as incident nodes. */
					assert( elType->nodeCount == cMap->nIncNodesLocal[cElInd] );
					
					/* Collect a set of top-level node coordinate pointers. */
					tNodeCoordPtrs = Memory_Alloc_Array( Coord*, cMap->nIncNodesLocal[cElInd], "MultiGrid" );
					for( cNode_j = 0; cNode_j < cMap->nIncNodesLocal[cElInd]; cNode_j++ ) {
						unsigned	tNodeInd = cMap->nodesTop[cMap->incNodesLocal[cElInd][cNode_j]];
						
						tNodeCoordPtrs[cNode_j] = &mesh->nodeCoord[tNodeInd];
					}
					
					/* Loop over the previous-level nodes incident on this coarse element. */
					for( pNode_i = 0; pNode_i < cMap->nIncNodesPrev[cElInd]; pNode_i++ ) {
						unsigned	pNodeInd = cMap->incNodesPrev[cElInd][pNode_i];
						unsigned	tNodeInd = fMap ? fMap->nodesTop[pNodeInd] : pNodeInd;
						unsigned	cDof_i;

						/* If we aren't looking at a direct mapping, skip it. */
						if( tNodeInd != cMap->nodesTop[cNode_i] ) {
							continue;
						}
						
						/* Calculate the shape function values of this top-node with respect to the coarse node on this
						   element. */
						ElementType_ConvertGlobalCoordToElLocal( elType, 
											 elLayout, 
											 (const Coord**)tNodeCoordPtrs, 
											 mesh->nodeCoord[tNodeInd], 
											 lCoord );
						ElementType_EvaluateShapeFunctionsAt( elType, lCoord, sfVals );
						
						if( sfVals[elNodeInd] != 0.0 ) {
							/* Loop over the DOFs on the previous node (column DOFs). */
							for( cDof_i = 0; cDof_i < nColDofs[pNodeInd]; cDof_i++ ) {
								unsigned	cDofInd = colEqNums[pNodeInd][cDof_i];
								
								if( cDofInd == -1 || cDof_i != rDof_i ) {
									continue;
								}

								/* Have we already covered this index? */
								{
									unsigned	col_i;

									for( col_i = 0; col_i < rowNZs[curRow]; col_i++ ) {
										if( indices[curRow][col_i] == cDofInd ) {
											break;
										}
									}
									if( col_i < rowNZs[curRow] ) {
										continue;
									}
								}
								
								/* Expand the row arrays and store value. */
								indices[curRow] = Memory_Realloc_Array( indices[curRow], unsigned, rowNZs[curRow] + 1 );
								indices[curRow][rowNZs[curRow]] = cDofInd;
								entries[curRow] = Memory_Realloc_Array( entries[curRow], double, rowNZs[curRow] + 1 );
								entries[curRow][rowNZs[curRow]] = sfVals[elNodeInd];
								rowNZs[curRow]++;
								
							}
						}
					}
					
					/* Free some arrays. */
					FreeArray( sfVals );
					FreeArray( tNodeCoordPtrs );
				}
				
				/* Update the net number of non-zeros and max non-zeros. */
				nzs += rowNZs[curRow];
				maxNZs = (rowNZs[curRow] > maxNZs) ? rowNZs[curRow] : maxNZs;

				/* Update the current row. */
				rowIndices[curRow] = rDofInd;
				curRow++;
			}
		}
	}
	
	
	/*
	** Transfer the operator entries to a matrix.
	*/
	
	{
		unsigned	row_i;
		Matrix*		rOp;
		
		/* Allocate the matrix.  As it turns out, the non-zeros this function refers to is on a row by row basis.  Luckily, 
		   Galerkin style operators should have about the same number of non-zeros per row. */
		rOp = Matrix_New( comm, nRows, nCols, maxNZs );
		
		/* Transfer, row by row. */
		for( row_i = 0; row_i < nRows; row_i++ ) {
			if( rowNZs[row_i] ) {
				Matrix_Insert( rOp, 1, rowIndices + row_i, rowNZs[row_i], indices[row_i], entries[row_i] );
			}
		}
		
		/* Assemble the matrix. */
		Matrix_AssemblyBegin( rOp );
		Matrix_AssemblyEnd( rOp );

		/* Save. */
		*mat = rOp;
	}
	
	/* Free arrays. */
	FreeArray( rowNZs );
	FreeArray2D( nRows, indices );
	FreeArray2D( nRows, entries );
	FreeArray( rowIndices );
}


void MultiGrid_NormaliseOp( Matrix* op ) {
	unsigned	nRows, nCols;
	unsigned	row_i;

	assert( op );


	/*
	** We may need to set up a restriction inter-grid operator such that the scale is preserved.  Here we will
	** sum each row of the operator and divide each entry in the row by the total.
	*/

	Matrix_GetLocalSize( op, &nRows, &nCols );
	for( row_i = 0; row_i < nRows; row_i++ ) {
		unsigned     	nEntries;
		unsigned*	inds;
		unsigned*	indsCopy;
		double*		entries;
		double*		newEntries;

		/* Get the row details. */
		Matrix_GetRow( op, row_i, &nEntries, &inds, &entries );

		if( nEntries > 1 ) {
			double		sum = 0.0;
			unsigned	entry_i;

			/* We need new space to store modifications because the memory used for the entries
			   may be inside PETSc. We'll also need the indices.*/
			newEntries = Memory_Alloc_Array( double, nEntries, "MultiGrid" );
			indsCopy = Memory_Alloc_Array( unsigned, nEntries, "MultiGrid" );
			memcpy( indsCopy, inds, nEntries * sizeof(unsigned) );

			for( entry_i = 0; entry_i < nEntries; entry_i++ ) {
				sum += entries[entry_i];
			}
			sum = 1.0 / sum;
			for( entry_i = 0; entry_i < nEntries; entry_i++ ) {
				newEntries[entry_i] = entries[entry_i] * sum;
			}
		}
		else {
			newEntries = NULL;
			indsCopy = NULL;
		}

		/* Can't keep them... */
		Matrix_RestoreRow( op, row_i, nEntries, &inds, &entries );

		/* If there were mods, update the matrix. */
		if( nEntries > 1 ) {
			Matrix_Insert( op, 1, &row_i, nEntries, indsCopy, newEntries );
			FreeArray( indsCopy );
			FreeArray( newEntries );

			/* Re-assemble the matrix. */
			Matrix_AssemblyBegin( op );
			Matrix_AssemblyEnd( op );
		}
	}


}


void MultiGrid_BuildWorkVectors( unsigned handle ) {
	Vector*	tmplVec;
	unsigned	level_i;
	MGInfo*	info;
	
	
	assert( mgCtx );
	assert( handle < mgCtx->nInfos );
	
	info = mgCtx->infos + handle;
	
	
	/*
	** Validate status.
	*/
	
	assert( info->nLocalRows && info->nLocalCols );
	assert( info->stiffMat->rhs->vector );
	
	
	/*
	** Construct appropriately sized vectors for calculating residuals and solutions for each level (except the finest, 
	** which only needs a residual work vector).
	*/
	
	/* Allocate space for work vectors. */
	if( !info->rhsVecs ) {
		info->rhsVecs = Memory_Alloc_Array( Vector*, info->nLevels + 1, "MultiGrid" );
		memset( info->rhsVecs, 0, sizeof(Vector*) * (info->nLevels + 1) );
	}
	if( !info->rVecs ) {
		info->rVecs = Memory_Alloc_Array( Vector*, info->nLevels + 1, "MultiGrid" );
		memset( info->rVecs, 0, sizeof(Vector*) * (info->nLevels + 1) );
	}
	if( !info->xVecs ) {
		info->xVecs = Memory_Alloc_Array( Vector*, info->nLevels + 1, "MultiGrid" );
		memset( info->xVecs, 0, sizeof(Vector*) * (info->nLevels + 1) );
	}
	
	/* Set the finest level work vectors. */
	tmplVec = info->stiffMat->rhs->vector;
	info->rhsVecs[0] = NULL;
	info->xVecs[0] = NULL;
	if( !info->rVecs[0] ) {
		Vector_Duplicate( tmplVec, &info->rVecs[0] );
	}
	
	if( info->nLevels >= 1 ) {
		MPI_Comm	comm;
		
		/* Unfortunately, we still need the communicator to do a 'DupNewSize'. */
		comm = info->stiffMat->comm;
		
		/* Duplicate from template vector and resize accordingly. */
		if( info->rhsVecs[1] && Vector_LocalSize( info->rhsVecs[1] ) != info->nLocalRows[0] ) {
			Vector_Destroy( info->rhsVecs[1] );
			info->rhsVecs[1] = NULL;
		}
		if( !info->rhsVecs[1] ) {
			Vector_DupNewSize( comm, tmplVec, info->nLocalRows[0], &info->rhsVecs[1] );
		}
		
		if( info->xVecs[1] && Vector_LocalSize( info->xVecs[1] ) != info->nLocalRows[0] ) {
			Vector_Destroy( info->xVecs[1] );
			info->xVecs[1] = NULL;
		}
		if( !info->xVecs[1] ) {
			Vector_DupNewSize( comm, tmplVec, info->nLocalRows[0], &info->xVecs[1] );
		}
		
		/* Except for the coarse vector. */
		if( info->nLevels > 1 ) {
			if( info->rVecs[1] && Vector_LocalSize( info->rVecs[1] ) != info->nLocalRows[0] ) {
				Vector_Destroy( info->rVecs[1] );
				info->rVecs[1] = NULL;
			}
			if( !info->rVecs[1] ) {
				Vector_DupNewSize( comm, tmplVec, info->nLocalRows[0], &info->rVecs[1] );
			}
		}
		else if( info->rVecs[1] ) {
			Vector_Destroy( info->rVecs[1] );
			info->rVecs[1] = NULL;
		}
		
		/* Set the workspace vectors for all levels except the finest. */
		for( level_i = 2; level_i < info->nLevels + 1; level_i++ ) {
			/* Duplicate from previous vectors. */
			if( info->rhsVecs[level_i] && Vector_LocalSize( info->rhsVecs[level_i] ) != info->nLocalRows[level_i - 1] ) {
				Vector_Destroy( info->rhsVecs[level_i] );
				info->rhsVecs[level_i] = NULL;
			}
			if( !info->rhsVecs[level_i] ) {
				Vector_DupNewSize( comm, info->rhsVecs[level_i - 1], info->nLocalRows[level_i - 1], &info->rhsVecs[level_i] );
			}
			
			if( info->xVecs[level_i] && Vector_LocalSize( info->xVecs[level_i] ) != info->nLocalRows[level_i - 1] ) {
				Vector_Destroy( info->xVecs[level_i] );
				info->xVecs[level_i] = NULL;
			}
			if( !info->xVecs[level_i] ) {
				Vector_DupNewSize( comm, info->xVecs[level_i - 1], info->nLocalRows[level_i - 1], &info->xVecs[level_i] );
			}
			
			/* And we don't need a coarse residual. */
			if( level_i < info->nLevels ) {
				if( info->rVecs[level_i] && Vector_LocalSize( info->rVecs[level_i] ) != info->nLocalRows[level_i - 1] ) {
					Vector_Destroy( info->rVecs[level_i] );
					info->rVecs[level_i] = NULL;
				}
				if( !info->rVecs[level_i] ) {
					Vector_DupNewSize( comm, info->rVecs[level_i - 1], info->nLocalRows[level_i - 1], &info->rVecs[level_i] );
				}
			}
			else if( info->rVecs[level_i] ) {
				Vector_Destroy( info->rVecs[level_i] );
				info->rVecs[level_i] = NULL;
			}
		}
	}
}


void MultiGrid_DeleteContext( MGContext* ctx ) {
	
	/*
	** Free all the data structures contained in the context.
	*/
	
	/* Free the mappings. */
	{
		unsigned	map_i;
		
		for( map_i = 0; map_i < ctx->nMappings; map_i++ ) {
			MGGridMapping*	mapping = ctx->mappings + map_i;
			unsigned		level_i;
			

			
			for( level_i = 0; level_i < mapping->maxLevels; level_i++ ) {
				MGMapping*	map = mapping->maps + level_i;
				
				FreeArray( map->nodesTop );
				FreeArray( map->nIncNodesLocal );
				FreeArray2D( map->nEls, map->incNodesLocal );
				FreeArray( map->nIncNodesPrev );
				FreeArray2D( map->nEls, map->incNodesPrev );
				FreeArray( map->incElTop );
				FreeArray( map->nIncElsLocal );
				FreeArray2D( map->nNodes, map->incElsLocal );
				FreeArray2D( map->nNodes, map->elNodeInds );
				FreeArray( map->nDofs );
				FreeArray2D( map->nNodes, map->lmTable );
			}
			
			FreeArray( mapping->maps );
		}
		
		FreeArray( ctx->mappings );
	}
	
	/* Free the infos. */
	{
		unsigned	info_i;
		
		for( info_i = 0; info_i < ctx->nInfos; info_i++ ) {
			MGInfo*	info = ctx->infos + info_i;
			
			FreeArray( info->nDownIts );
			FreeArray( info->nUpIts );
			FreeArray( info->nCycles );
			
			if( info->smoothers ) {
				unsigned	smoother_i;
				
				for( smoother_i = 0; smoother_i < info->nLevels; smoother_i++ ) {
					if( info->smoothers[smoother_i] ) {
						Matrix_Destroy( info->smoothers[smoother_i] );
					}
				}
				FreeArray( info->smoothers );
			}
			
			if( info->rOps ) {
				unsigned	op_i;
				
				for( op_i = 0; op_i < info->nLevels; op_i++ ) {
					if( info->rOps[op_i] && !(info->pOps && (info->rOps[op_i] == info->pOps[op_i])) ) {
						Matrix_Destroy( info->rOps[op_i] );
					}
				}
				FreeArray( info->rOps );
			}
			
			if( info->pOps ) {
				unsigned	op_i;
				
				for( op_i = 0; op_i < info->nLevels; op_i++ ) {
					if( info->pOps[op_i] ) {
						Matrix_Destroy( info->pOps[op_i] );
					}
				}
				FreeArray( info->pOps );
			}

			if( info->iOps ) {
				unsigned	op_i;
				
				for( op_i = 0; op_i < info->nLevels; op_i++ ) {
					if( info->iOps[op_i] ) {
						Matrix_Destroy( info->iOps[op_i] );
					}
				}
				FreeArray( info->iOps );
			}
			
			FreeArray( info->nLocalRows );
			FreeArray( info->nLocalCols );
			
			if( info->rhsVecs ) {
				unsigned	vec_i;
				
				for( vec_i = 0; vec_i < info->nLevels; vec_i++ ) {
					if( info->rhsVecs[vec_i] ) {
						Vector_Destroy( info->rhsVecs[vec_i] );
					}
				}
				FreeArray( info->rhsVecs );
			}
			
			if( info->xVecs ) {
				unsigned	vec_i;
				
				for( vec_i = 0; vec_i < info->nLevels; vec_i++ ) {
					if( info->xVecs[vec_i] ) {
						Vector_Destroy( info->xVecs[vec_i] );
					}
				}
				FreeArray( info->xVecs );
			}
			
			if( info->rVecs ) {
				unsigned	vec_i;
				
				for( vec_i = 0; vec_i < info->nLevels; vec_i++ ) {
					if( info->rVecs[vec_i] ) {
						Vector_Destroy( info->rVecs[vec_i] );
					}
				}
				FreeArray( info->rVecs );
			}
		}
	}
}
