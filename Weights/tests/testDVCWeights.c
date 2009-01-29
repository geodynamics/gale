/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	Robert B. Turnbull, Monash Cluster Computing. (Robert.Turnbull@sci.monash.edu.au)
**	Kathleen M. Humble, Computational Scientist, VPAC. (khumble@vpac.org)
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
** $Id: testTensorMath.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PopulationControl/PopulationControl.h>

#include "PICellerator/Weights/Weights.h"

#include <stdio.h>
#include <stdlib.h>

int main( int argc, char* argv[] ) {
	MPI_Comm CommWorld;
	int rank;
	int numProcessors;
	int procToWatch;
	Stream*  stream;
	
	/* Initialise MPI, get world info */
	MPI_Init( &argc, &argv );
	MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
	MPI_Comm_size( CommWorld, &numProcessors );
	MPI_Comm_rank( CommWorld, &rank );
	
	StGermain_Init( &argc, &argv );
	StgDomain_Init( &argc, &argv );
	StgFEM_Init( &argc, &argv );
    PICellerator_Voronoi_Init( &argc, &argv );
    PICellerator_PopulationControl_Init( &argc, &argv );
	PICellerator_Weights_Init( &argc, &argv );

	/* TODO Not sure how to use the DVCWeights Init as it doesn't match the other ones */
	MPI_Barrier( CommWorld ); /* Ensures copyright info always come first in output */
	
	stream = Journal_Register( InfoStream_Type, "testDVCWeights" );

	stJournal->firewallProducesAssert = False;
	
	
	if( argc >= 2 ) {
		procToWatch = atoi( argv[1] );
	}
	else {
		procToWatch = 0;
	}

	if( rank == procToWatch ) {
	
		double dx,dy,dz,da;
		static struct cell *cells;/* the 3D connected grid */ 
		struct particle *pList;/* 3D particle List */
		struct chain *bchain;/*3D boundary chain */
			
		static struct cell2d *cells2D;/* the 2D connected grid */
		struct particle2d *pList2D;/*2D particle List */
		struct chain *bchain2D;/*2D boundary chain */
			
		int nump,numx,numy,numz;
		int px, py, pz;
		double BBXMIN = -1.0; /* the ranges of the local coordinates of a FEM cell.*/
		double BBXMAX = 1.0;
		double BBYMIN = -1.0;
		double BBYMAX = 1.0;
		double BBZMIN = -1.0;
		double BBZMAX = 1.0;
		int i,j,k,l;
		double particleDistance;		
		double particle0[3], particle1[3];
		
		/* Do the testing of the 3D functions*/

		Journal_Printf( stream, "-----------------------------\n");	
		Journal_Printf( stream, "3D functions:  \n");
		Journal_Printf( stream, "-----------------------------\n\n");	


		Journal_Printf( stream, "****************************\n");
		Journal_Printf( stream, "Test _DVCWeights_ConstructGrid function \n");
		Journal_Printf( stream, "****************************\n");
		
		/*Define the resolution */
		
		numx = 2;
		numy = 2;
		numz = 2;

		Journal_Printf( stream, "size of element:\n\t x = (%f, %f)\n\t y = (%f, %f) \n\t z = (%f, %f)\n",
			BBXMIN, BBXMAX, BBYMIN, BBYMAX, BBZMIN, BBZMAX);
		Journal_Printf( stream, "Resolution: \n\t (x, y, z) = (%d, %d, %d)\n", numx, numy, numz);
		
     		_DVCWeights_ConstructGrid(&cells,numz,numy,numx,BBXMIN,BBYMIN,BBZMIN,BBXMAX,BBYMAX,BBZMAX);		
		
		/* Print out the grid somehow */
		for (i = 0; i < (numx * numy * numz); i++ ) {
			Journal_Printf(stream, "cell[%d]:\tParticle Index: %d\n", 
					i, cells[i].p);
			Journal_Printf(stream, " \t\tValues: (N: %d, S: %d, E: %d, W: %d, U: %d, D: %d) \n", 
					cells[i].N, cells[i].S, cells[i].E, cells[i].W,cells[i].U, cells[i].D );
			Journal_Printf(stream, " \t\tCoords: (%f, %f,%f) \t Done = %d\n", 
					cells[i].x, cells[i].y, cells[i].z, cells[i].done);			
		}

		Journal_Printf( stream, "****************************\n");
		Journal_Printf( stream, "Test _DVCWeights_InitialiseStructs function \n");
		Journal_Printf( stream, "****************************\n");
		
		/*Define the number of particles */
		px = 2;
		py = 2; 
		pz = 2;
		nump = px * py * pz;
		
		_DVCWeights_InitialiseStructs( &bchain, &pList, nump);
		for (i = 0; i < nump; i++) {
			Journal_Printf( stream, "bchain[%d]: ", i);
			Journal_Printf( stream, "No of new_claimed_cells = %d, ", 
				bchain[i].new_claimed_cells_malloced);
			Journal_Printf( stream, "No of new_bound_cells = %d\n",
				bchain[i].new_bound_cells_malloced);
		}
		
		Journal_Printf( stream, "****************************\n");
		Journal_Printf( stream, "Test _DVCWeights_ResetGrid function \n");
		Journal_Printf( stream, "****************************\n");
		Journal_Printf( stream, "numz * numy * numx = %d\n", numz*numy*numx);
		
		_DVCWeights_ResetGrid(&cells,numz*numy*numx);

		for ( i = 0; i < numz*numy*numx; i++) {
			Journal_Printf( stream, "cells[%d].p = %d \t cells[%d].done = %d\n",
				i, cells[i].p, i, cells[i].done);
			
		}
		
		Journal_Printf( stream, "****************************\n");
		Journal_Printf( stream, "Test _DVCWeights_CreateVoronoi function \n");
		Journal_Printf( stream, "****************************\n");

		dx = (BBXMAX - BBXMIN)/numx;
		dy = (BBYMAX - BBYMIN)/numy;
		dz = (BBZMAX - BBZMIN)/numz;
		da = dx*dy*dz;
		/*Initialise particle coords */
		l = 0;
		for(i = 0; i < px ;i++){
	    	for ( j = 0; j < py ; j++) {
				for ( k = 0; k < pz; k++ ) { 
	      			pList[l].x = (1 + i) / (px + 1.0);
	      			pList[l].y = (1 + j) / ( py + 1.0);
	      			pList[l].z = (1 + k) / (pz + 1.0);
					l++;
				}
		  	}
	  	}
		for ( i = 0; i < nump; i++) {	    
			Journal_Printf( stream, "pList[%d]:", i);
			Journal_Printf( stream, "\t\t coords: (x, y, z) = (%f, %f, %f)\n",
				pList[i].x, pList[i].y, pList[i].z);

		}
		Journal_Printf( stream, "\n(dx, dy, dz) = (%f, %f, %f)	da = %f\n\n",
			dx, dy, dz, da);		
		_DVCWeights_CreateVoronoi( &bchain, &pList, &cells, dx, dy, dz, nump, numx, numy, numz, BBXMIN, BBXMAX, BBYMIN, BBYMAX, BBZMIN, BBZMAX);
		
		/* print out bchain changes */
		for (i = 0; i < nump; i++) {
			Journal_Printf( stream, "bchain[%d]: \t index = %d \n",
					i, bchain[i].index);
			Journal_Printf( stream, "\t\t sizeofboundary = %d \n\t\t numclaimed = %d \n",
					bchain[i].sizeofboundary, bchain[i].numclaimed);
			Journal_Printf( stream, "\t\t totalclaimed = %d\n", bchain[i].totalclaimed);
			Journal_Printf( stream, "\t\t new_bound_cells_malloced = %d \n",
					bchain[i].new_bound_cells_malloced);
			Journal_Printf( stream, "\t\t new_claimed_cells_malloced = %d \n",
					bchain[i].new_claimed_cells_malloced);
			Journal_Printf( stream, "\t\t done = %d\n", bchain[i].done);
		}
		/* Print out particle values */
		for (i = 0; i < nump; i++) {
		
			Journal_Printf( stream, "pList[%d]:\n", i);
			Journal_Printf( stream, "\t\t coords: (x, y, z) = (%f, %f, %f)\n",
				pList[i].x, pList[i].y, pList[i].z);
			Journal_Printf( stream, "\t\t centroids: (cx, cy, cz) = (%f, %f %f)\n",
				pList[i].cx, pList[i].cy, pList[i].cz);
			Journal_Printf( stream, "\t\t weight = %f\n", pList[i].w);
		
		}
		Journal_Printf( stream, "****************************\n");
		Journal_Printf( stream, "Test _DVCWeights_GetCentroids function \n");
		Journal_Printf( stream, "****************************\n");

		_DVCWeights_GetCentroids( cells, pList,numz,numy,numx,nump,da);
		for (i = 0; i < nump; i++) {
		
			Journal_Printf( stream, "pList[%d]:\n", i);
			Journal_Printf( stream, "\t\t coords: (x, y, z) = (%f, %f, %f)\n",
				pList[i].x, pList[i].y, pList[i].z);
			Journal_Printf( stream, "\t\t centroids: (cx, cy, cz) = (%f, %f %f)\n",
				pList[i].cx, pList[i].cy, pList[i].cz);
			Journal_Printf( stream, "\t\t weight = %f\n", pList[i].w);
		
		}
		Journal_Printf( stream, "****************************\n");
		Journal_Printf( stream, "Test _DVCWeights_DistanceSquared function \n");
		Journal_Printf( stream, "****************************\n");
		particle0[0] = 0.5;	particle0[1] = 0.5;	particle0[2] = 0.5;
		particle1[0] = 0.25; particle1[1] = 0.25; 	particle1[2] = 0; 		
		Journal_Printf( stream, "particle0:\n");
		Journal_Printf( stream, "\t\t coords: (x, y, z) = (%f, %f, %f)\n",
				particle0[0], particle0[1], particle0[2]);
		Journal_Printf( stream, "particle1:\n");
		Journal_Printf( stream, "\t\t coords: (x, y, z) = (%f, %f, %f)\n",
				particle1[0], particle1[1], particle1[2]);
		
		particleDistance = _DVCWeights_DistanceSquared(
			particle0[0], particle0[1], particle0[2],
			particle1[0], particle1[1], particle1[2]);
		Journal_Printf( stream, "calculated distance^2 between particles = %f \n", particleDistance);
		
		free(bchain);
		free(pList);
		
		Journal_Printf( stream, "\n-----------------------------\n");	
		Journal_Printf( stream, "2D functions:  \n");
		Journal_Printf( stream, "-----------------------------\n\n");

		/*********************************/
		/* Test construct grid 2D */

		Journal_Printf( stream, "****************************\n");
		Journal_Printf( stream, "Test _DVCWeights_ConstructGrid2D function \n");
		Journal_Printf( stream, "****************************\n");
		
		/*Define the resolution */
		
		numx = 2;
		numy = 2;

		Journal_Printf( stream, "size of element:\n\t x = (%f, %f)\n\t y = (%f, %f) \n",
			BBXMIN, BBXMAX, BBYMIN, BBYMAX);
		Journal_Printf( stream, "Resolution: \n\t (x, y) = (%d, %d)\n", numx, numy);
		
     		_DVCWeights_ConstructGrid2D(&cells2D,numy,numx, BBXMIN,BBYMIN,BBXMAX,BBYMAX);		
		
		/* Print out the grid somehow */
		for (i = 0; i < (numx * numy ); i++ ) {
			Journal_Printf(stream, "cells2d[%d]:\tParticle Index: %d \n", 
					i, cells2D[i].p);
			Journal_Printf(stream, " \t\tValues: (N: %d, S: %d, E: %d, W: %d) \n", 
					cells2D[i].N, cells2D[i].S, cells2D[i].E, cells2D[i].W );
			Journal_Printf(stream, " \t\tCoords: (%f, %f) \t Done = %d\n", 
					cells2D[i].x, cells2D[i].y, cells[i].done);			
		}

		Journal_Printf( stream, "****************************\n");
		Journal_Printf( stream, "Test _DVCWeights_InitialiseStructs2D function \n");
		Journal_Printf( stream, "****************************\n");
		
		/*Define size of swarm-to-be */
		px = 2;
		py = 2; 
		nump = px * py ;
		
		_DVCWeights_InitialiseStructs2D( &bchain2D, &pList2D, nump);
		for (i = 0; i < nump; i++) {
			Journal_Printf( stream, "bchain2D[%d]: ", i);
			Journal_Printf( stream, "No of new_claimed_cells = %d, ", 
				bchain2D[i].new_claimed_cells_malloced);
			Journal_Printf( stream, "No of new_bound_cells = %d\n",
				bchain2D[i].new_bound_cells_malloced);
		}
		
		Journal_Printf( stream, "****************************\n");
		Journal_Printf( stream, "Test _DVCWeights_ResetGrid2D function \n");
		Journal_Printf( stream, "****************************\n");
		Journal_Printf( stream, "numz * numy = %d\n", numz*numy);
		
		_DVCWeights_ResetGrid2D(&cells2D,numx*numy);

		for ( i = 0; i < numx*numy; i++) {
			Journal_Printf( stream, "cells2D[%d].p = %d \t cells2D[%d].done = %d\n",
				i, cells2D[i].p, i, cells2D[i].done);
			
		}
		
		Journal_Printf( stream, "****************************\n");
		Journal_Printf( stream, "Test _DVCWeights_CreateVoronoi2D function \n");
		Journal_Printf( stream, "****************************\n");

		dx = (BBXMAX - BBXMIN)/numx;
		dy = (BBYMAX - BBYMIN)/numy;
		da = dx*dy;
		/*Initialise particle coords */
		l = 0;
		for(i = 0; i < px ;i++){
	    	for ( j = 0; j < py ; j++) {
	      		pList2D[l].x = (1 + i) / (px + 1.0);
	      		pList2D[l].y = (1 + j) / ( py + 1.0);
				l++;
		  	}
	  	}
		for ( i = 0; i < nump; i++) {	    
			Journal_Printf( stream, "pList2D[%d]:", i);
			Journal_Printf( stream, "\t\t coords: (x, y) = (%f, %f)\n",
				pList2D[i].x, pList2D[i].y);

		}
		Journal_Printf( stream, "\n(dx, dy) = (%f, %f)	da = %f\n\n",
			dx, dy, da);		
		_DVCWeights_CreateVoronoi2D( &bchain2D, &pList2D, &cells2D, dx, dy, nump, numx, numy, BBXMIN, BBXMAX, BBYMIN, BBYMAX);
		
		/* print out bchain changes */
		for (i = 0; i < nump; i++) {
			Journal_Printf( stream, "bchain2D[%d]: \t  index = %d \n",
					i, bchain2D[i].index);
			Journal_Printf( stream, "\t\t sizeofboundary = %d \n\t\t numclaimed = %d \n",
					bchain2D[i].sizeofboundary, bchain2D[i].numclaimed);
			Journal_Printf( stream, "\t\t totalclaimed = %d\n", bchain2D[i].totalclaimed);
			Journal_Printf( stream, "\t\t new_bound_cells_malloced = %d \n",
					bchain2D[i].new_bound_cells_malloced);
			Journal_Printf( stream, "\t\t new_claimed_cells_malloced = %d \n",
					bchain2D[i].new_claimed_cells_malloced);
			Journal_Printf( stream, "\t\t done = %d\n", bchain2D[i].done);
		}
		/* Print out particle values */
		for (i = 0; i < nump; i++) {
		
			Journal_Printf( stream, "pList2D[%d]:\n", i);
			Journal_Printf( stream, "\t\t coords: (x, y) = (%f, %f)\n",
				pList2D[i].x, pList2D[i].y);
			Journal_Printf( stream, "\t\t centroids: (cx, cy) = (%f, %f)\n",
				pList2D[i].cx, pList2D[i].cy);
			Journal_Printf( stream, "\t\t weight = %f\n", pList2D[i].w);
		
		}
		Journal_Printf( stream, "****************************\n");
		Journal_Printf( stream, "Test _DVCWeights_GetCentroids2D function \n");
		Journal_Printf( stream, "****************************\n");

		_DVCWeights_GetCentroids2D( cells2D, pList2D,numy,numx,nump,da);
		for (i = 0; i < nump; i++) {
		
			Journal_Printf( stream, "pList2D[%d]:\n", i);
			Journal_Printf( stream, "\t\t coords: (x, y) = (%f, %f)\n",
				pList2D[i].x, pList2D[i].y);
			Journal_Printf( stream, "\t\t centroids: (cx, cy) = (%f, %f)\n",
				pList2D[i].cx, pList2D[i].cy);
			Journal_Printf( stream, "\t\t weight = %f\n", pList2D[i].w);
		
		}
		Journal_Printf( stream, "****************************\n");
		Journal_Printf( stream, "Test _DVCWeights_DistanceSquared2D function \n");
		Journal_Printf( stream, "****************************\n");
		
		Journal_Printf( stream, "particle0:\n");
		Journal_Printf( stream, "\t\t coords: (x, y) = (%f, %f)\n",
				particle0[0], particle0[1]);
		Journal_Printf( stream, "particle1:\n");
		Journal_Printf( stream, "\t\t coords: (x, y) = (%f, %f)\n",
				particle1[0], particle1[1]);
		
		particleDistance = _DVCWeights_DistanceSquared2D(
			particle0[0], particle0[1],
			particle1[0], particle1[1] );
		Journal_Printf( stream, "calculated distance^2 between particles = %f \n", particleDistance);
		
		
		
		free(bchain2D);
		free(pList2D);

	}
		
	PICellerator_Weights_Finalise();	
	PICellerator_PopulationControl_Finalise();
	PICellerator_Voronoi_Finalise();
	StgFEM_Finalise();
	StgDomain_Finalise();
	StGermain_Finalise();

	/* Close off MPI */
	MPI_Finalize();
	
	return 0;
}
