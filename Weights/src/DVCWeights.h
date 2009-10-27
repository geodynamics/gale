/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2006, Monash Cluster Computing 
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
** Author:
**              Mirko Velic - Mirko.Velic@sci.monash.edu.au
**
**  Assumptions:
**  	 I am assuming that the xi's (local coords) on the IntegrationPoint particles
**       are precalculated somewhere and get reset based on material PIC positions each time step.
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


/****************************************************************************************************************

  The algorithm here-in computes a discrete voronoi diagram per FEM cell given a set of local
  particle positions, in 3D and 2D. The volumes of the Voronoi regions are used as integration weights for
  the integration point swarm and the integration points are the centroids of the same volumes.

  For a description of this algorithm, see the article by Velic et.al.
     "A Fast Robust Algorithm for computing Discrete Voronoi Diagrams in N-dimensions"



*****************************************************************************************************************/

#ifndef __PICellerator_Weights_DVCWeightsClass_h__
#define __PICellerator_Weights_DVCWeightsClass_h__

/* Textual name of this class */
extern const Type DVCWeights_Type;

/* DVCWeights information */
#define __DVCWeights                            \
    /* General info */                          \
    __WeightsCalculator                         \
                                                \
    /* Virtual Info */                          \
    /* Parameters that are passed in */         \
    int resX;                                   \
    int resY;                                   \
    int resZ;

struct DVCWeights { __DVCWeights };
	
#define DVC_INC 150
	
struct cell{
    int p;/*particle index number*/
    int index;
    int N;
    int S;
    int E;
    int W;
    int U;
    int D;
    double x;
    double y;
    double z;
    int done;
};
struct cell2d{
    int p;/* particle index number */
    int index;
    int N;
    int S;
    int E;
    int W;
    double x;
    double y;
    int done;
};
struct chain{
    int p;/*particle index number*/
    int index;/*cell number in grid*/
    int sizeofboundary; /* number of cells on boundary so far */
    int numclaimed;
    int totalclaimed;
    int new_bound_cells_malloced;
    int new_claimed_cells_malloced; 
    int *new_bound_cells;
    int *new_claimed_cells;
    int done;
};
struct particle{
    double x;
    double y;
    double z;
    double cx;
    double cy;
    double cz;
    double w;
    int index;
};
struct particle2d{
    double x;
    double y;
    double cx;
    double cy;
    double w;
    int index;
};
/* Private function prototypes in 3D */
void   _DVCWeights_GetCentroids(	
    struct cell *cells,
    struct particle *pList,
    int n,
    int m, 
    int l,
    int nump,
    double vol);
void   _DVCWeights_ClaimCells(
    struct chain **bchain,
    struct cell **cells,
    struct particle **pList,
    int p_i);
void   _DVCWeights_UpdateBchain(
    struct chain **bchain,
    struct cell **cells,
    int p_i);
void   _DVCWeights_ResetGrid(struct cell **cells, int n );

double _DVCWeights_DistanceSquared(
    double x0, double y0, double z0, 
    double x1, double y1, double z1);
double _DVCWeights_DistanceTest(double x0, double y0, double z0, 
                                double x1, double y1, double z1,
                                double x2, double y2, double z2);	
void   _DVCWeights_ConstructGrid(
    struct cell **cell_list, 
    int n, int m, int l,
    double x0, double y0, double z0,
    double x1, double y1, double z1);
void   _DVCWeights_InitialiseStructs( 
    struct chain **bchain, 
    struct particle **pList, 
    int nump);
void   _DVCWeights_CreateVoronoi( 
    struct chain **bchain, 
    struct particle **pList, 
    struct cell **cells, 
    double dx, double dy, double dz,
    int nump,
    int numx, int numy, int numz, 
    double BBXMIN, double BBXMAX, 
    double BBYMIN, double BBYMAX,
    double BBZMIN, double BBZMAX);
/* Private function prototypes in 2D */
void   _DVCWeights_GetCentroids2D( 
    struct cell2d *cells, 
    struct particle2d *pList,
    int n, 
    int m, 
    int nump,
    double vol);
void   _DVCWeights_ClaimCells2D(
    struct chain **bchain,
    struct cell2d **cells,
    struct particle2d **pList,
    int p_i);
void   _DVCWeights_UpdateBchain2D(
    struct chain **bchain,
    struct cell2d **cells,
    int p_i);
void   _DVCWeights_ResetGrid2D(struct cell2d **cells, int n );
double _DVCWeights_DistanceSquared2D(double x0, double y0, double x1, double y1);
double _DVCWeights_DistanceTest2D(double x0, double y0,
                                  double x1, double y1,
                                  double x2, double y2);
void   _DVCWeights_ConstructGrid2D(
    struct cell2d **cell_list, 
    int m, int l,
    double x0, double y0,
    double x1, double y1);
void   _DVCWeights_InitialiseStructs2D( 
    struct chain **bchain, 
    struct particle2d **pList, 
    int nump);
void   _DVCWeights_CreateVoronoi2D( 
    struct chain **bchain, 
    struct particle2d **pList, 
    struct cell2d **cells, 
    double dx, double dy,
    int nump,
    int numx, int numy,
    double BBXMIN, double BBXMAX, 
    double BBYMIN, double BBYMAX);

void   _DVCWeights_Calculate2D( 
    void* dvcWeights, 
    void* _swarm, 
    Cell_LocalIndex lCell_I );
void   _DVCWeights_Calculate3D( 
    void* dvcWeights, 
    void* _swarm, 
    Cell_LocalIndex lCell_I );


/*---------------------------------------------------------------------------------------------------------------------
** Constructors
*/
DVCWeights* DVCWeights_New( Name name, Dimension_Index dim, int *res ) ;

DVCWeights* _DVCWeights_New(
    SizeT                                 _sizeOfSelf, 
    Type                                  type,
    Stg_Class_DeleteFunction*             _delete,
    Stg_Class_PrintFunction*              _print,
    Stg_Class_CopyFunction*               _copy, 
    Stg_Component_DefaultConstructorFunction* _defaultConstructor,
    Stg_Component_ConstructFunction*      _construct,
    Stg_Component_BuildFunction*          _build,
    Stg_Component_InitialiseFunction*     _initialise,
    Stg_Component_ExecuteFunction*        _execute,
    Stg_Component_DestroyFunction*        _destroy,		
    WeightsCalculator_CalculateFunction*  _calculate,
    Name                                  name,
    Bool                                  initFlag,
    int                                   dim,
    int*                                  res );

void _DVCWeights_Init( void* dvcWeights, int *res ) ;

/* Stg_Class_Delete DVCWeights implementation */
void _DVCWeights_Delete( void* dvcWeights );
void _DVCWeights_Print( void* dvcWeights, Stream* stream );
#define DVCWeights_Copy( self )                                         \
    (DVCWeights*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
#define DVCWeights_DeepCopy( self )                                     \
    (DVCWeights*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
void* _DVCWeights_Copy( void* dvcWeights, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
void* _DVCWeights_DefaultNew( Name name ) ;

void _DVCWeights_Construct( void* dvcWeights, Stg_ComponentFactory* cf, void* data ) ;

void _DVCWeights_Build( void* dvcWeights, void* data ) ;
void _DVCWeights_Initialise( void* dvcWeights, void* data ) ;
void _DVCWeights_Execute( void* dvcWeights, void* data );
void _DVCWeights_Destroy( void* dvcWeights, void* data ) ;
	
		
void _DVCWeights_Calculate( void* dvcWeights, void* _swarm, Cell_LocalIndex lCell_I ) ;

#endif
