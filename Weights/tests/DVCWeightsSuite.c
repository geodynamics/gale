/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**   Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**   Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**   Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**   Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**   Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**   Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
** $Id: testList.c 2136 2004-09-30 02:47:13Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include "PICellerator/PopulationControl/PopulationControl.h"
#include "PICellerator/Weights/Weights.h"


/* the ranges of the local coordinates of a FEM cell.*/
const double BBXMIN = -1.0; 
const double BBXMAX = 1.0;
const double BBYMIN = -1.0;
const double BBYMAX = 1.0;
const double BBZMIN = -1.0;
const double BBZMAX = 1.0;

typedef struct {
   struct cell*         cells;   /* the 3D connected grid */ 
   struct particle*     pList;   /* 3D particle List */
   struct chain*        bchain;  /*3D boundary chain */

   struct cell2d*       cells2D; /* the 2D connected grid */
   struct particle2d*   pList2D; /*2D particle List */
   struct chain*        bchain2D;/*2D boundary chain */
   unsigned int         numx,numy,numz;
   unsigned int         nump;
   unsigned int         px, py, pz;
   double               dx,dy,dz;
   double               da;
} DVCWeightsSuiteData;


void DVCWeightsSuite_Setup( DVCWeightsSuiteData* data ) {
   data->cells = NULL;
   data->pList = NULL;
   data->bchain = NULL;
   data->cells2D = NULL;
   data->pList2D = NULL;
   data->bchain2D = NULL;

   /*Define the resolution */
   data->numx = 2;
   data->numy = 2;
   data->numz = 2;

   /*Define the number of particles */
   data->px = 2;
   data->py = 2; 
   data->pz = 2;
   data->nump = data->px * data->py * data->pz;

   data->dx = (BBXMAX - BBXMIN)/data->numx;
   data->dy = (BBYMAX - BBYMIN)/data->numy;
   data->dz = (BBZMAX - BBZMIN)/data->numz;
   data->da = data->dx*data->dy*data->dz;
}


void DVCWeightsSuite_Teardown( DVCWeightsSuiteData* data ) {
   if(data->bchain) { free(data->bchain); }
   if(data->pList) { free(data->pList); }
   if(data->cells) { free(data->cells); }

   if(data->bchain2D) { free(data->bchain2D); }
   if(data->pList2D) { free(data->pList2D); }
   if(data->cells2D) { free(data->cells2D); }
}


/* Do the testing of the 3D functions*/
void DVCWeightsSuite_TestConstructGrid( DVCWeightsSuiteData* data ) {
   unsigned int   ii=0;
   Stream*        stream = Journal_Register( Info_Type, "TestConstructGrid" );
   const char*    gridFilename = "DVCWeightsSuite_testConstructGrid.txt";
   char           expectedGridFilename[PCU_PATH_MAX];
   
   /* We'll use an expected file for this test, as visual checking of the correct parameters is probably
    *  most logical way to maintain it */
   Stream_RedirectFile( stream, gridFilename );

   _DVCWeights_ConstructGrid(&data->cells, data->numz, data->numy, data->numx, BBXMIN,BBYMIN,BBZMIN,BBXMAX,BBYMAX,BBZMAX);		
   
   /* Print out the grid */
   for (ii = 0; ii < (data->numx * data->numy * data->numz); ii++ ) {
      pcu_check_true( data->cells[ii].p == -1 );   /* Particle index: shouldn't be set up yet */
      pcu_check_true( data->cells[ii].done == 0 );

      Journal_Printf(stream, "cell[%d]\n", ii);
      Journal_Printf(stream, " \t\tValues: (N: %d, S: %d, E: %d, W: %d, U: %d, D: %d) \n", 
            data->cells[ii].N, data->cells[ii].S, data->cells[ii].E, data->cells[ii].W,
            data->cells[ii].U, data->cells[ii].D );
      Journal_Printf(stream, " \t\tCoords: (%f, %f,%f) \n", 
            data->cells[ii].x, data->cells[ii].y, data->cells[ii].z );
   }

   pcu_filename_expected( gridFilename, expectedGridFilename );
   pcu_check_fileEq( gridFilename, expectedGridFilename );
   remove( gridFilename );
}


void DVCWeightsSuite_TestInitialiseStructs( DVCWeightsSuiteData* data ) {
   unsigned int   ii;   

   
   _DVCWeights_InitialiseStructs( &data->bchain, &data->pList, data->nump);
   for (ii = 0; ii < data->nump; ii++) {
      pcu_check_true( data->bchain[ii].new_claimed_cells_malloced == DVC_INC );
      pcu_check_true( data->bchain[ii].new_bound_cells_malloced == DVC_INC );
   }
}
   

void DVCWeightsSuite_TestResetGrid( DVCWeightsSuiteData* data ) {
   unsigned int   i;
   
   _DVCWeights_ConstructGrid(&data->cells, data->numz, data->numy, data->numx, BBXMIN,BBYMIN,BBZMIN,BBXMAX,BBYMAX,BBZMAX);		
   _DVCWeights_ResetGrid(&data->cells, data->numz*data->numy*data->numx);

   for ( i = 0; i < data->numz*data->numy*data->numx; i++) {
      pcu_check_true( data->cells[i].p == -1 );
      pcu_check_true( data->cells[i].done == 0 );
   }
}


void DVCWeightsSuite_TestCreateVoronoi( DVCWeightsSuiteData* data ) {
   int            i,j,k,l;

   _DVCWeights_ConstructGrid(&data->cells, data->numz, data->numy, data->numx, BBXMIN,BBYMIN,BBZMIN,BBXMAX,BBYMAX,BBZMAX);		
   _DVCWeights_InitialiseStructs( &data->bchain, &data->pList, data->nump);

   /*Initialise particle coords */
   l = 0;
   for(i = 0; i < data->px ;i++){
      for ( j = 0; j < data->py ; j++) {
         for ( k = 0; k < data->pz; k++ ) { 
               data->pList[l].x = (1 + i) / (data->px + 1.0);
               data->pList[l].y = (1 + j) / ( data->py + 1.0);
               data->pList[l].z = (1 + k) / (data->pz + 1.0);
            l++;
         }
      }
   }
   for ( i = 0; i < data->nump; i++) {	    
      //Journal_Printf( stream, "data->pList[%d]:", i);
      //Journal_Printf( stream, "\t\t coords: (x, y, z) = (%f, %f, %f)\n",
      //   data->pList[i].x, data->pList[i].y, data->pList[i].z);

   }

   _DVCWeights_CreateVoronoi( &data->bchain, &data->pList, &data->cells, data->dx, data->dy, data->dz,
      data->nump, data->numx, data->numy, data->numz, BBXMIN, BBXMAX, BBYMIN, BBYMAX, BBZMIN, BBZMAX);
   
   /* data->bchain changes */
   for (i = 0; i < data->nump; i++) {
      pcu_check_true( data->bchain[i].index == 7 );
      pcu_check_true( data->bchain[i].sizeofboundary == 0 );
      pcu_check_true( data->bchain[i].numclaimed == 0 );
      pcu_check_true( data->bchain[i].totalclaimed == 1 );
      pcu_check_true( data->bchain[i].new_bound_cells_malloced == DVC_INC );
      pcu_check_true( data->bchain[i].new_claimed_cells_malloced == DVC_INC );
      pcu_check_true( data->bchain[i].done == 0 );
   }
   /* particle values */
   for (i = 0; i < data->nump; i++) {
      pcu_check_true( (data->pList[i].cx == 0) && (data->pList[i].cy == 0) && (data->pList[i].cz == 0));
      pcu_check_true( data->pList[i].w == 0 );
   }
}


void DVCWeightsSuite_TestGetCentroids( DVCWeightsSuiteData* data ) {
   unsigned int   i;
   Stream*        stream = Journal_Register( Info_Type, "TestGetCentroids" );
   const char*    centroidsFilename = "DVCWeightsSuite_testGetCentroids.txt";
   char           expectedCentroidsFilename[PCU_PATH_MAX];
   
   /* We'll use an expected file for this test, as visual checking of the correct parameters is probably
    *  most logical way to maintain it */
   Stream_RedirectFile( stream, centroidsFilename );

   _DVCWeights_ConstructGrid(&data->cells, data->numz, data->numy, data->numx, BBXMIN,BBYMIN,BBZMIN,BBXMAX,BBYMAX,BBZMAX);		
   _DVCWeights_InitialiseStructs( &data->bchain, &data->pList, data->nump);
   _DVCWeights_CreateVoronoi( &data->bchain, &data->pList, &data->cells, data->dx, data->dy, data->dz,
      data->nump, data->numx, data->numy, data->numz, BBXMIN, BBXMAX, BBYMIN, BBYMAX, BBZMIN, BBZMAX);

   _DVCWeights_GetCentroids( data->cells, data->pList,data->numz,data->numy,data->numx,data->nump,data->da);

   for (i = 0; i < data->nump; i++) {
      Journal_Printf( stream, "data->pList[%d]:\n", i);
      Journal_Printf( stream, "\t\t coords: (x, y, z) = (%f, %f, %f)\n",
         data->pList[i].x, data->pList[i].y, data->pList[i].z);
      Journal_Printf( stream, "\t\t centroids: (cx, cy, cz) = (%f, %f %f)\n",
         data->pList[i].cx, data->pList[i].cy, data->pList[i].cz);
      Journal_Printf( stream, "\t\t weight = %f\n", data->pList[i].w);
   }

   pcu_filename_expected( centroidsFilename, expectedCentroidsFilename );
   pcu_check_fileEq( centroidsFilename, expectedCentroidsFilename );
   remove( centroidsFilename );
}


void DVCWeightsSuite_TestDistanceSquared( DVCWeightsSuiteData* data ) {
   double particleDistance;		
   double particle0[3], particle1[3];

   particle0[0] = 0.5;	particle0[1] = 0.5;	particle0[2] = 0.5;
   particle1[0] = 0.25; particle1[1] = 0.25; 	particle1[2] = 0; 		

   particleDistance = _DVCWeights_DistanceSquared(
      particle0[0], particle0[1], particle0[2],
      particle1[0], particle1[1], particle1[2]);

   pcu_check_true( particleDistance == 0.375 );
}


#if 0

/* 2D Functions */

void DVCWeightsSuite_TestConstructGrid2D( DVCWeightsSuiteData* data ) {
   int data->numx,data->numy,data->numz;
   
   /*Define the resolution */
   
   data->numx = 2;
   data->numy = 2;

   Journal_Printf( stream, "size of element:\n\t x = (%f, %f)\n\t y = (%f, %f) \n",
      BBXMIN, BBXMAX, BBYMIN, BBYMAX);
   Journal_Printf( stream, "Resolution: \n\t (x, y) = (%d, %d)\n", data->numx, data->numy);
   
      _DVCWeights_ConstructGrid2D(&data->cells2D,data->numy,data->numx, BBXMIN,BBYMIN,BBXMAX,BBYMAX);		
   
   /* Print out the grid somehow */
   for (i = 0; i < (data->numx * data->numy ); i++ ) {
      Journal_Printf(stream, "data->cells2d[%d]:\tParticle Index: %d \n", 
            i, data->cells2D[i].p);
      Journal_Printf(stream, " \t\tValues: (N: %d, S: %d, E: %d, W: %d) \n", 
            data->cells2D[i].N, data->cells2D[i].S, data->cells2D[i].E, data->cells2D[i].W );
      Journal_Printf(stream, " \t\tCoords: (%f, %f) \t Done = %d\n", 
            data->cells2D[i].x, data->cells2D[i].y, data->cells[i].done);			
   }
}


void DVCWeightsSuite_TestInitialiseStructs2D( DVCWeightsSuiteData* data ) {
   int data->nump;
   int data->px, data->py;
   
   /*Define size of swarm-to-be */
   data->px = 2;
   data->py = 2; 
   data->nump = data->px * data->py ;
   
   _DVCWeights_InitialiseStructs2D( &data->bchain2D, &data->pList2D, data->nump);
   for (i = 0; i < data->nump; i++) {
      Journal_Printf( stream, "data->bchain2D[%d]: ", i);
      Journal_Printf( stream, "No of new_claimed_cells = %d, ", 
         data->bchain2D[i].new_claimed_cells_malloced);
      Journal_Printf( stream, "No of new_bound_cells = %d\n",
         data->bchain2D[i].new_bound_cells_malloced);
   }
}

   
void DVCWeightsSuite_TestResetGrid2D( DVCWeightsSuiteData* data ) {
   Journal_Printf( stream, "data->numz * data->numy = %d\n", data->numz*data->numy);
   
   _DVCWeights_ResetGrid2D(&data->cells2D,data->numx*data->numy);

   for ( i = 0; i < data->numx*data->numy; i++) {
      Journal_Printf( stream, "data->cells2D[%d].p = %d \t data->cells2D[%d].done = %d\n",
         i, data->cells2D[i].p, i, data->cells2D[i].done);
   }
}

   
void DVCWeightsSuite_TestCreateVoronoi2D( DVCWeightsSuiteData* data ) {
   double dx,dy,dz,da;
   int i,j,l;

   dx = (BBXMAX - BBXMIN)/data->numx;
   dy = (BBYMAX - BBYMIN)/data->numy;
   da = dx*dy;
   /*Initialise particle coords */
   l = 0;
   for(i = 0; i < data->px ;i++){
      for ( j = 0; j < data->py ; j++) {
            data->pList2D[l].x = (1 + i) / (data->px + 1.0);
            data->pList2D[l].y = (1 + j) / ( data->py + 1.0);
         l++;
      }
   }
   for ( i = 0; i < data->nump; i++) {	    
      Journal_Printf( stream, "data->pList2D[%d]:", i);
      Journal_Printf( stream, "\t\t coords: (x, y) = (%f, %f)\n",
         data->pList2D[i].x, data->pList2D[i].y);

   }
   Journal_Printf( stream, "\n(dx, dy) = (%f, %f)	da = %f\n\n",
      dx, dy, da);		
   _DVCWeights_CreateVoronoi2D( &data->bchain2D, &data->pList2D, &data->cells2D, dx, dy, data->nump, data->numx, data->numy, BBXMIN, BBXMAX, BBYMIN, BBYMAX);
   
   /* print out data->bchain changes */
   for (i = 0; i < data->nump; i++) {
      Journal_Printf( stream, "data->bchain2D[%d]: \t  index = %d \n",
            i, data->bchain2D[i].index);
      Journal_Printf( stream, "\t\t sizeofboundary = %d \n\t\t numclaimed = %d \n",
            data->bchain2D[i].sizeofboundary, data->bchain2D[i].numclaimed);
      Journal_Printf( stream, "\t\t totalclaimed = %d\n", data->bchain2D[i].totalclaimed);
      Journal_Printf( stream, "\t\t new_bound_cells_malloced = %d \n",
            data->bchain2D[i].new_bound_cells_malloced);
      Journal_Printf( stream, "\t\t new_claimed_cells_malloced = %d \n",
            data->bchain2D[i].new_claimed_cells_malloced);
      Journal_Printf( stream, "\t\t done = %d\n", data->bchain2D[i].done);
   }
   /* Print out particle values */
   for (i = 0; i < data->nump; i++) {
   
      Journal_Printf( stream, "data->pList2D[%d]:\n", i);
      Journal_Printf( stream, "\t\t coords: (x, y) = (%f, %f)\n",
         data->pList2D[i].x, data->pList2D[i].y);
      Journal_Printf( stream, "\t\t centroids: (cx, cy) = (%f, %f)\n",
         data->pList2D[i].cx, data->pList2D[i].cy);
      Journal_Printf( stream, "\t\t weight = %f\n", data->pList2D[i].w);
   }
}


void DVCWeightsSuite_TestGetCentroids2D( DVCWeightsSuiteData* data ) {

   _DVCWeights_GetCentroids2D( data->cells2D, data->pList2D,data->numy,data->numx,data->nump,da);
   for (i = 0; i < data->nump; i++) {
   
      Journal_Printf( stream, "data->pList2D[%d]:\n", i);
      Journal_Printf( stream, "\t\t coords: (x, y) = (%f, %f)\n",
         data->pList2D[i].x, data->pList2D[i].y);
      Journal_Printf( stream, "\t\t centroids: (cx, cy) = (%f, %f)\n",
         data->pList2D[i].cx, data->pList2D[i].cy);
      Journal_Printf( stream, "\t\t weight = %f\n", data->pList2D[i].w);
   
   }
}


void DVCWeightsSuite_TestDistanceSquared2D( DVCWeightsSuiteData* data ) {
   double particleDistance;		
   double particle0[3], particle1[3];
   
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
   
}
#endif


void DVCWeightsSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, DVCWeightsSuiteData );
   pcu_suite_setFixtures( suite, DVCWeightsSuite_Setup, DVCWeightsSuite_Teardown );
   pcu_suite_addTest( suite, DVCWeightsSuite_TestConstructGrid );
   pcu_suite_addTest( suite, DVCWeightsSuite_TestInitialiseStructs );
   pcu_suite_addTest( suite, DVCWeightsSuite_TestCreateVoronoi );
   pcu_suite_addTest( suite, DVCWeightsSuite_TestGetCentroids );
   pcu_suite_addTest( suite, DVCWeightsSuite_TestDistanceSquared );
#if 0
   pcu_suite_addTest( suite, DVCWeightsSuite_TestConstructGrid2D );
   pcu_suite_addTest( suite, DVCWeightsSuite_TestInitialiseStructs2D );
   pcu_suite_addTest( suite, DVCWeightsSuite_TestCreateVoronoi2D );
   pcu_suite_addTest( suite, DVCWeightsSuite_TestGetCentroids2D );
   pcu_suite_addTest( suite, DVCWeightsSuite_TestDistanceSquared2D );
#endif
}
