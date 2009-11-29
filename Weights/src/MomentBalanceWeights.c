/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**      Melbourne, 3053, Australia.
** Copyright (c) 2005-2006, Monash Cluster Computing, Building 28, Monash University Clayton Campus,
**      Victoria, 3800, Australia
**
** Primary Contributing Organisations:
**      Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**      Australian Computational Earth Systems Simulator - http://www.access.edu.au
**      Monash Cluster Computing - http://www.mcc.monash.edu.au
**
** Contributors:
**      Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**      Patrick D. Sunter, Software Engineer, VPAC. (patrick@vpac.org)
**      Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**      Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**      David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**      Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**      Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**      Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**      Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**      Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**      David Stegman, Postdoctoral Fellow, Monash University. (david.stegman@sci.monash.edu.au)
**      Wendy Sharples, PhD Student, Monash University (wendy.sharples@sci.monash.edu.au)
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
** $Id: MomentBalanceWeights.c 518 2007-10-11 08:07:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include "types.h"
#include "WeightsCalculator.h"
#include "MomentBalanceWeights.h"
#include "ConstantWeights.h"

#include <assert.h>
#include <string.h>
#include <math.h>


/* Textual name of this class */
const Type MomentBalanceWeights_Type = "MomentBalanceWeights";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

MomentBalanceWeights* MomentBalanceWeights_New( Name name, Dimension_Index dim, WeightsCalculator* backupWeights ) {
    MomentBalanceWeights* self = _MomentBalanceWeights_DefaultNew( name );

    self->isConstructed = True;
    _WeightsCalculator_Init( self, dim );
    _MomentBalanceWeights_Init( self, backupWeights );
}

MomentBalanceWeights* _MomentBalanceWeights_New(  MOMENTBALANCEWEIGHTS_DEFARGS  ) {
    MomentBalanceWeights* self;

    /* Allocate memory */
    assert( _sizeOfSelf >= sizeof(MomentBalanceWeights) );
    self = (MomentBalanceWeights*)_WeightsCalculator_New(  WEIGHTSCALCULATOR_PASSARGS  ); 

    /* General info */

    /* Virtual Info */

    return self;
}

void _MomentBalanceWeights_Init( void* momentBalanceWeights, WeightsCalculator* backupWeights ) {
    MomentBalanceWeights* self = (MomentBalanceWeights*)momentBalanceWeights;

    if ( backupWeights ) {
        self->backupWeights = backupWeights;
        self->freeBackupWeights = False;
    }
    else {
        self->backupWeights = (WeightsCalculator*) ConstantWeights_New( "backupWeights", self->dim );
        self->freeBackupWeights = True;
    }
        
    Journal_Firewall( self->dim == 2, Journal_Register( Error_Type, self->type ), "%s only works in 2D.\n", self->type );
        
}

/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _MomentBalanceWeights_Delete( void* momentBalanceWeights ) {
    MomentBalanceWeights* self = (MomentBalanceWeights*)momentBalanceWeights;
        
    /* Delete parent */
    _WeightsCalculator_Delete( self );
}

void _MomentBalanceWeights_Print( void* momentBalanceWeights, Stream* stream ) {
    MomentBalanceWeights* self = (MomentBalanceWeights*)momentBalanceWeights;
        
    /* Print parent */
    _WeightsCalculator_Print( self, stream );
}

void* _MomentBalanceWeights_Copy( void* momentBalanceWeights, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
    MomentBalanceWeights*       self = (MomentBalanceWeights*)momentBalanceWeights;
    MomentBalanceWeights*       newMomentBalanceWeights;
        
    newMomentBalanceWeights = (MomentBalanceWeights*)_WeightsCalculator_Copy( self, dest, deep, nameExt, ptrMap );
        
    return (void*)newMomentBalanceWeights;
}

void* _MomentBalanceWeights_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(MomentBalanceWeights);
	Type                                                      type = MomentBalanceWeights_Type;
	Stg_Class_DeleteFunction*                              _delete = _MomentBalanceWeights_Delete;
	Stg_Class_PrintFunction*                                _print = _MomentBalanceWeights_Print;
	Stg_Class_CopyFunction*                                  _copy = _MomentBalanceWeights_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _MomentBalanceWeights_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _MomentBalanceWeights_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _MomentBalanceWeights_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _MomentBalanceWeights_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _MomentBalanceWeights_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _MomentBalanceWeights_Destroy;
	AllocationType                              nameAllocationType = NON_GLOBAL;
	WeightsCalculator_CalculateFunction*                _calculate = _MomentBalanceWeights_Calculate;

    return (void*) _MomentBalanceWeights_New(  MOMENTBALANCEWEIGHTS_PASSARGS  );
}

void _MomentBalanceWeights_AssignFromXML( void* momentBalanceWeights, Stg_ComponentFactory* cf, void* data ) {
    MomentBalanceWeights*            self          = (MomentBalanceWeights*) momentBalanceWeights;
    WeightsCalculator*           backupWeights;

    _WeightsCalculator_AssignFromXML( self, cf, data );

    backupWeights =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "BackupWeights", WeightsCalculator, False, data );

    _MomentBalanceWeights_Init( self, backupWeights );
}

void _MomentBalanceWeights_Build( void* momentBalanceWeights, void* data ) {
    MomentBalanceWeights*       self = (MomentBalanceWeights*)momentBalanceWeights;
    
    Stg_Component_Build( self->backupWeights, data, False );
    _WeightsCalculator_Build( self, data );
}
void _MomentBalanceWeights_Initialise( void* momentBalanceWeights, void* data ) {
    MomentBalanceWeights*       self = (MomentBalanceWeights*)momentBalanceWeights;

    Stg_Component_Initialise( self->backupWeights, data, False );
    _WeightsCalculator_Initialise( self, data );
}
void _MomentBalanceWeights_Execute( void* momentBalanceWeights, void* data ) {
    MomentBalanceWeights*       self = (MomentBalanceWeights*)momentBalanceWeights;
        
    _WeightsCalculator_Execute( self, data );
}
void _MomentBalanceWeights_Destroy( void* momentBalanceWeights, void* data ) {
    MomentBalanceWeights*       self = (MomentBalanceWeights*)momentBalanceWeights;

    if ( self->freeBackupWeights )
       Stg_Component_Destroy( self->backupWeights, data, False );
    _WeightsCalculator_Destroy( self, data );
}


/*-------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

typedef struct {
    double         x;
    double         y;
    double         z;
    double         weight;
    Particle_Index particleCount;
} RegionInfo;

typedef enum {
    REGION_A,
    REGION_B,
    REGION_C,
    REGION_D,
    REGION_E,
    REGION_F,
    REGION_G,
    REGION_H
} RegionIndex;


#define _MomentBalanceWeights_WhichRegion2D( xi )               \
    (( xi[ I_AXIS ] < 0.0 )                                     \
     ? (( xi[ J_AXIS ] < 0.0 ) ? REGION_C : REGION_B )          \
     : (( xi[ J_AXIS ] < 0.0 ) ? REGION_D : REGION_A ) )

void _MomentBalanceWeights_Calculate( void* momentBalanceWeights, void* _swarm, Cell_LocalIndex lCell_I ) {
    MomentBalanceWeights*        self            = (MomentBalanceWeights*)  momentBalanceWeights;
    Swarm*                       swarm           = (Swarm*) _swarm;
    Particle_InCellIndex         cParticleCount  = swarm->cellParticleCountTbl[lCell_I];
    Particle_InCellIndex         cParticle_I;
    RegionInfo                   region[4];
    RegionIndex                  region_I;
    double*                      xi;
    double                       n;
    double                       m;
    IntegrationPoint*            particle;
    double                       a,b,c;
                
    /* First Step -
     * Assume all weights are constant in quadrants - i.e.                  
     * -------------
     * |  B  |  A  |
     * |  C  |  D  |
     * -------------
     *
     * This gives us three constaint equations:
     * W_A N_A + W_B N_B + W_C N_C + W_D N_D = 4
     * W_A X_A + W_D X_D = W_B X_B + W_C X_C
     * W_A Y_A + W_B Y_B = W_C Y_C + W_D Y_D 
     *
     * where
     *
     * W_Q is the weight for each particle in region Q,
     * N_Q is the number of particles in region Q,
     * X_Q = \sum_{q=1}^{N_q} |\xi_{q} | and
     * Y_Q = \sum_{q=1}^{N_q} |\eta_{q}| */

    memset( region, 0, 4 * sizeof(RegionInfo) );

    /* Calculate region data from particle local coordinates */
    for ( cParticle_I = 0 ; cParticle_I < cParticleCount ; cParticle_I++ ) {
        particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, lCell_I, cParticle_I );
        xi       = particle->xi;

        region_I = _MomentBalanceWeights_WhichRegion2D( xi );
        region[ region_I ].x += fabs(xi[ I_AXIS ]);
        region[ region_I ].y += fabs(xi[ J_AXIS ]);

        region[ region_I ].particleCount++;
    }

    /* Make sure each region has particles - otherwise use back up */
    if ( region[ REGION_A ].particleCount == 0 ||
         region[ REGION_B ].particleCount == 0 ||
         region[ REGION_C ].particleCount == 0 ||
         region[ REGION_D ].particleCount == 0 ) 
    {
        WeightsCalculator_CalculateCell( self->backupWeights, swarm, lCell_I );
        return;
    }


    /** Let  n = W_B/W_C and
     *       m = W_D/W_C 
     *
     * Assume that the ratios of W_A/W_B and W_D/W_C are the same and
     * (which also means that W_A/W_D = W_B/W_C)
     *  
     * Therefore:
     *       nm = W_A/W_C 
     *
     * This gives us the constraint equations
     *
     * nm N_A + n N_B + N_C + m N_D = 4/W_C
     * nm X_A + m X_D = n X_B + X_C
     * nm Y_A + n Y_B = Y_C + m Y_D
     *
     * Therefore:
     * m = \frac{n X_B + X_C}{n X_A + X_D} = \frac{ Y_C - n Y_B }{ nY_A - Y_D }
     * and
     * (X_B Y_A + Y_B X_A)n^2 + (X_C Y_A + Y_B X_D - Y_D X_B - X_A Y_C)n - (X_C Y_D + Y_C X_D) = 0
     * which we can solve using the quadratic formula */

    a = ( region[ REGION_A ].y * region[ REGION_B ].x + region[ REGION_A ].x * region[ REGION_B ].y );
    b = ( region[ REGION_A ].y * region[ REGION_C ].x - region[ REGION_A ].x * region[ REGION_C ].y ) +
        ( region[ REGION_B ].y * region[ REGION_D ].x - region[ REGION_B ].x * region[ REGION_D ].y );
    c = - ( region[ REGION_D ].y * region[ REGION_C ].x + region[ REGION_D ].x * region[ REGION_C ].y );

    /* We take the solution which will always give a positive n */
    n = (-b + sqrt( b*b - 4.0 * a * c ))/(2.0 * a);
                
    m = ( n * region[ REGION_B ].x + region[ REGION_C ].x)/
        ( n * region[ REGION_A ].x + region[ REGION_D ].x);

    /* Calculate weights from intermediate values */
    region[ REGION_C ].weight = self->cellLocalVolume /
        ( (double) region[ REGION_A ].particleCount * n * m
          + (double) region[ REGION_B ].particleCount * n
          + (double) region[ REGION_C ].particleCount  
          + (double) region[ REGION_D ].particleCount * m );

    region[ REGION_A ].weight = n * m * region[ REGION_C ].weight;
    region[ REGION_B ].weight = n *     region[ REGION_C ].weight;
    region[ REGION_D ].weight =     m * region[ REGION_C ].weight;

    /* Assign weights to particles according to which region it is in */
    for ( cParticle_I = 0 ; cParticle_I < cParticleCount ; cParticle_I++ ) {
        particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, lCell_I, cParticle_I );
        xi       = particle->xi;

        region_I = _MomentBalanceWeights_WhichRegion2D( xi );
        particle->weight = region[ region_I ].weight;
    }
}


#define _MomentBalanceWeights_WhichRegion3D( xi )               \
    (( xi[ K_AXIS ] < 0.0 )                                     \
     ? ( (( xi[ I_AXIS ] < 0.0 )                                \
          ? (( xi[ J_AXIS ] < 0.0 ) ? REGION_C : REGION_B )     \
          : (( xi[ J_AXIS ] < 0.0 ) ? REGION_D : REGION_A ) ))  \
     : ( (( xi[ I_AXIS ] < 0.0 )                                \
          ? (( xi[ J_AXIS ] < 0.0 ) ? REGION_G : REGION_F )     \
          : (( xi[ J_AXIS ] < 0.0 ) ? REGION_H : REGION_E ) )))

void _MomentBalanceWeights_Calculate3D( void* momentBalanceWeights, void* _swarm, Cell_LocalIndex lCell_I ) {
    MomentBalanceWeights*        self            = (MomentBalanceWeights*)  momentBalanceWeights;
    Swarm*                       swarm           = (Swarm*) _swarm;
    Particle_InCellIndex         cParticleCount  = swarm->cellParticleCountTbl[lCell_I];
    Particle_InCellIndex         cParticle_I;
    RegionInfo                   region[8];
    RegionIndex                  region_I;
    double*                      xi;
    double                       n;
    double                       m;
/*      double                       q; */
    IntegrationPoint*            particle;
    double                       a,b,c;
                
    /* First Step -
     * Assume all weights are constant in octants - i.e.                  
     * -FRONT-FACE--
     * |  B  |  A  |
     * |  C  |  D  |
     * -------------
     *
     * --BACK-FACE--
     * |  F  |  E  |
     * |  G  |  H  |
     * -------------
     *
     * This gives us four constaint equations:
     * W_A N_A + W_B N_B + W_C N_C + W_D N_D + W_E N_E + W_F N_F + W_G N_G + W_H N_H= 8
     * W_A X_A + W_D X_D + W_E X_E + W_H X_H = W_B X_B + W_C X_C + W_F X_F + W_G X_G
     * W_A Y_A + W_B Y_B + W_F Y_F + W_E Y_E = W_C Y_C + W_D Y_D + W_G Y_G + W_H Y_H
     * W_A Z_A + W_B Z_B + W_C Z_C + W_D Z_D = W_E Z_E + W_F Z_F + W_G Z_G + W_H Z_H
     *
     * where
     *
     * W_Q is the weight for each particle in region Q,
     * N_Q is the number of particles in region Q,
     * X_Q = \sum_{q=1}^{N_q} |\xi_{q}  | and
     * Y_Q = \sum_{q=1}^{N_q} |\eta_{q} | 
     * Z_Q = \sum_{q=1}^{N_q} |\zeta_{q}| */

    memset( region, 0, 8 * sizeof(RegionInfo) );

    /* Calculate region data from particle local coordinates */
    for ( cParticle_I = 0 ; cParticle_I < cParticleCount ; cParticle_I++ ) {
        particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, lCell_I, cParticle_I );
        xi       = particle->xi;

        region_I = _MomentBalanceWeights_WhichRegion3D( xi );
        region[ region_I ].x += fabs(xi[ I_AXIS ]);
        region[ region_I ].y += fabs(xi[ J_AXIS ]);
        region[ region_I ].z += fabs(xi[ K_AXIS ]);

        region[ region_I ].particleCount++;
    }

    /* Make sure each region has particles - otherwise use back up */
    if ( region[ REGION_A ].particleCount == 0 ||
         region[ REGION_B ].particleCount == 0 ||
         region[ REGION_C ].particleCount == 0 ||
         region[ REGION_D ].particleCount == 0 ||
         region[ REGION_E ].particleCount == 0 ||
         region[ REGION_F ].particleCount == 0 ||
         region[ REGION_G ].particleCount == 0 ||
         region[ REGION_H ].particleCount == 0 ) 
    {
        WeightsCalculator_CalculateCell( self->backupWeights, swarm, lCell_I );
        return;
    }


    /** Let  n = W_B/W_C and
     *       m = W_D/W_C 
     *       q = W_G/W_C 
     *
     * To keep ratios even - as in 2D case assume:
     *  
     * Therefore:
     *       nm  = W_A/W_C 
     *       nmq = W_E/W_C 
     *       nq  = W_F/W_C 
     *       mq  = W_H/W_C 
     *
     * This gives us the constrain equations
     *
     * n X_B + X_C + nq X_F + q X_G      = nm X_A + m X_D + nmq X_E + mq X_H
     * n Y_B + nm Y_A + nq Y_F + nmq Y_E = Y_C + m Y_D + q Y_G + mq Y_H
     * nm Z_A + n Z_B + Z_C + m Z_D      = nmq Z_E + nq Z_F + q Z_G + mq Z_H
     *
     * Therefore, we obtain 3 equations for m:
     * m = \frac{ n X_B + X_C + nq X_F + q X_G }{ n X_A + X_D + nq X_E + q X_H }
     *   = \frac{ n Y_B + nq Y_F - Y_C - q Y_G }{ Y_D + q Y_H + - n Y_A - nq Y_E }
     *   TODO other one
     *
     * which gives two quadratics for n in terms of q:
     *
     *   n^2 ( q^2(X_E Y_F + X_F Y_E) + q(X_E Y_B + X_B Y_E + X_A Y_F + X_F Y_A) + (X_B Y_A + X_A Y_B) )
     * + n   ( q^2( X_G Y_E - X_E Y_G + X_H Y_F - X_F Y_H ) + 
     *                          + q ( X_H Y_B - X_B Y_H + X_D Y_F - X_F Y_D + X_C Y_E - Y_C X_E + X_G Y_A - X_A Y_G )
     *                          ( X_C Y_A - X_A Y_C ) + (X_D Y_B - X_B Y_D) )
     * - ( q^2( X_H Y_G + X_G Y_H ) + q( X_G Y_D + X_D Y_G + X_H Y_C + X_C Y_H ) + X_D Y_C + X_C Y_D ) = 0
     *
     * and
     *   TODO other one
     *
     * which we can solve using the quadratic formula:
     *
     * n = 
     * */

    a = ( region[ REGION_A ].y * region[ REGION_B ].x + region[ REGION_A ].x * region[ REGION_B ].y );
    b = ( region[ REGION_A ].y * region[ REGION_C ].x - region[ REGION_A ].x * region[ REGION_C ].y ) +
        ( region[ REGION_B ].y * region[ REGION_D ].x - region[ REGION_B ].x * region[ REGION_D ].y );
    c = - ( region[ REGION_D ].y * region[ REGION_C ].x + region[ REGION_D ].x * region[ REGION_C ].y );

    /* We take the solution which will always give a positive n */
    n = (-b + sqrt( b*b - 4.0 * a * c ))/(2.0 * a);
                
    m = ( n * region[ REGION_B ].x + region[ REGION_C ].x)/
        ( n * region[ REGION_A ].x + region[ REGION_D ].x);

    /* Calculate weights from intermediate values */
    region[ REGION_C ].weight = self->cellLocalVolume /
        ( (double) region[ REGION_A ].particleCount * n * m
          + (double) region[ REGION_B ].particleCount * n
          + (double) region[ REGION_C ].particleCount  
          + (double) region[ REGION_D ].particleCount * m );

    region[ REGION_A ].weight = n * m * region[ REGION_C ].weight;
    region[ REGION_B ].weight = n *     region[ REGION_C ].weight;
    region[ REGION_D ].weight =     m * region[ REGION_C ].weight;

    /* Assign weights to particles according to which region it is in */
    for ( cParticle_I = 0 ; cParticle_I < cParticleCount ; cParticle_I++ ) {
        particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, lCell_I, cParticle_I );
        xi       = particle->xi;

        region_I = _MomentBalanceWeights_WhichRegion2D( xi );
        particle->weight = region[ region_I ].weight;
    }
}

/*-------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/




