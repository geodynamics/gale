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
** $Id: IterativeWeights.c 518 2007-10-11 08:07:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include "types.h"
#include "WeightsCalculator.h"
#include "ConstantWeights.h"
#include "IterativeWeights.h"

#include <assert.h>
#include <string.h>
#include <math.h>

/** IterativeWeights is an implementation of Frederic Dufour's weights routine from Ellipsis which is described
 * in his PhD Thesis Section 2.6.1 pp. 63-66 */

/* Textual name of this class */
const Type IterativeWeights_Type = "IterativeWeights";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

IterativeWeights* _IterativeWeights_New(  ITERATIVEWEIGHTS_DEFARGS  ) {
    IterativeWeights* self;

    /* Allocate memory */
    assert( _sizeOfSelf >= sizeof(IterativeWeights) );
    self = (IterativeWeights*)_ConstantWeights_New(  CONSTANTWEIGHTS_PASSARGS  );
	
    /* General info */

    /* Virtual Info */
	
    return self;
}

void _IterativeWeights_Init( void* iterativeWeights, WeightsCalculator* initialWeights,
                             Iteration_Index maxIterations, double tolerance, double alpha )
{
    IterativeWeights* self = (IterativeWeights*)iterativeWeights;

    if ( initialWeights ) {
        self->initialWeights = initialWeights;
        self->freeInitialWeights = False;
    }
    else {
        self->initialWeights = (WeightsCalculator*) ConstantWeights_New( "initialWeights", self->dim );
        self->freeInitialWeights = True;
    }

    self->maxIterations = maxIterations;
    self->tolerance = tolerance;
    self->alpha = alpha;
	
    Journal_Firewall( self->dim == 2, Journal_Register( Error_Type, self->type ), "%s only works in 2D.\n", self->type );
}

/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _IterativeWeights_Delete( void* iterativeWeights ) {
    IterativeWeights* self = (IterativeWeights*)iterativeWeights;
	
    if (self->freeInitialWeights )
        Stg_Class_Delete( self->initialWeights );
    /* Delete parent */
    _ConstantWeights_Delete( self );
}


void _IterativeWeights_Print( void* iterativeWeights, Stream* stream ) {
    IterativeWeights* self = (IterativeWeights*)iterativeWeights;
	
    /* Print parent */
    _ConstantWeights_Print( self, stream );
}



void* _IterativeWeights_Copy( void* iterativeWeights, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
    IterativeWeights*	self = (IterativeWeights*)iterativeWeights;
    IterativeWeights*	newIterativeWeights;
	
    newIterativeWeights = (IterativeWeights*)_ConstantWeights_Copy( self, dest, deep, nameExt, ptrMap );
	
    return (void*)newIterativeWeights;
}

void* _IterativeWeights_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(IterativeWeights);
	Type                                                      type = IterativeWeights_Type;
	Stg_Class_DeleteFunction*                              _delete = _IterativeWeights_Delete;
	Stg_Class_PrintFunction*                                _print = _IterativeWeights_Print;
	Stg_Class_CopyFunction*                                  _copy = _IterativeWeights_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _IterativeWeights_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _IterativeWeights_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _IterativeWeights_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _IterativeWeights_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _IterativeWeights_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _IterativeWeights_Destroy;
	AllocationType                              nameAllocationType = NON_GLOBAL;
	WeightsCalculator_CalculateFunction*                _calculate = _IterativeWeights_Calculate;

    return (void*) _IterativeWeights_New(  ITERATIVEWEIGHTS_PASSARGS  );
}


void _IterativeWeights_AssignFromXML( void* iterativeWeights, Stg_ComponentFactory* cf, void* data ) {
    IterativeWeights*	     self          = (IterativeWeights*) iterativeWeights;
    WeightsCalculator*       initialWeights;
    Iteration_Index          maxIterations;
    double                   tolerance;
    double                   alpha;

    _ConstantWeights_AssignFromXML( self, cf, data );

    initialWeights =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "InitialWeights", WeightsCalculator, False, data );

    maxIterations = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "maxIterations", 10 );
    tolerance     = Stg_ComponentFactory_GetDouble( cf, self->name, "tolerance", 0.01 );
    alpha         = Stg_ComponentFactory_GetDouble( cf, self->name, "alpha", 0.8 ); /* 0.8 is default in Dufour p. 65 */
	
    _IterativeWeights_Init( self, initialWeights, maxIterations, tolerance, alpha );
}

void _IterativeWeights_Build( void* iterativeWeights, void* data ) {
    IterativeWeights*	self = (IterativeWeights*)iterativeWeights;

    Stg_Component_Build( self->initialWeights, data, False ); 
    _ConstantWeights_Build( self, data );
}
void _IterativeWeights_Initialise( void* iterativeWeights, void* data ) {
    IterativeWeights*	self = (IterativeWeights*)iterativeWeights;

    Stg_Component_Initialise( self->initialWeights, data, False ); 	
    _ConstantWeights_Initialise( self, data );
}
void _IterativeWeights_Execute( void* iterativeWeights, void* data ) {
    IterativeWeights*	self = (IterativeWeights*)iterativeWeights;
	
    _ConstantWeights_Execute( self, data );
}
void _IterativeWeights_Destroy( void* iterativeWeights, void* data ) {
    IterativeWeights*	self = (IterativeWeights*)iterativeWeights;

    if (self->freeInitialWeights )
       Stg_Component_Destroy( self->initialWeights, data, False ); 	
    _ConstantWeights_Destroy( self, data );
}

/*-------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
void _IterativeWeights_Calculate( void* iterativeWeights, void* _swarm, Cell_LocalIndex lCell_I ) {
    IterativeWeights*            self            = (IterativeWeights*)  iterativeWeights;
    Swarm*                       swarm           = (Swarm*) _swarm;
    Iteration_Index              iteration_I;
    double                       constaintLHS[3] = { 0.0,0.0,0.0 };
    double                       squareTerms[3]  = { 0.0,0.0,1.0 };
    double                       constaintError;
    Particle_InCellIndex         cParticle_I;
    Particle_InCellIndex         cParticleCount  = swarm->cellParticleCountTbl[lCell_I];
    IntegrationPoint*            particle;
    Dimension_Index              dim_I;
    Dimension_Index              dim             = self->dim;

    /* Initialise Weights Using Constant Weights routine - Dufour Equation 2.52 */
    WeightsCalculator_CalculateCell( self->initialWeights, swarm, lCell_I );
		
    for ( iteration_I = 0 ; iteration_I < self->maxIterations ; iteration_I++ ) {
        constaintLHS[ I_AXIS ] = WeightsCalculator_GetConstraintLHS( self, swarm, lCell_I, 1, 0, 0 );
        constaintLHS[ J_AXIS ] = WeightsCalculator_GetConstraintLHS( self, swarm, lCell_I, 0, 1, 0 );
        if ( self->dim == 3 )
            constaintLHS[ K_AXIS ] = WeightsCalculator_GetConstraintLHS( self, swarm, lCell_I, 0, 0, 1 );

        /* Test for convergence - See Dufour Equation 2.55 */
        constaintError = 
            fabs(constaintLHS[ I_AXIS ]) + 
            fabs(constaintLHS[ J_AXIS ]) + 
            fabs(constaintLHS[ K_AXIS ]);

        if ( constaintError < self->tolerance )
            /* Then this cell has converged - bail out of loop */
            break;

        /* Adjust weights according to Dufour Equation 2.54 */
        squareTerms[ I_AXIS ] = WeightsCalculator_GetLocalCoordSum( self, swarm, lCell_I, 2, 0, 0 );
        squareTerms[ J_AXIS ] = WeightsCalculator_GetLocalCoordSum( self, swarm, lCell_I, 0, 2, 0 );
        if ( self->dim == 3 )
            squareTerms[ K_AXIS ] = WeightsCalculator_GetLocalCoordSum( self, swarm, lCell_I, 0, 0, 2 );

        for ( cParticle_I = 0 ; cParticle_I < cParticleCount ; cParticle_I++ ) {
            particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, lCell_I, cParticle_I );
			
            for ( dim_I = 0 ; dim_I < dim ; dim_I++ )
                particle->weight -= self->alpha * constaintLHS[ dim_I ] / squareTerms[ dim_I ] * particle->xi[ dim_I ] ;
        }

        /* Scale Weights to ensure constant constaint */
        IterativeWeights_ScaleForConstantConstraint( self, swarm, lCell_I );
    }
}


/*-------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void IterativeWeights_ScaleForConstantConstraint( void* iterativeWeights, void* _swarm, Cell_LocalIndex lCell_I ) {
    IterativeWeights*            self            = (IterativeWeights*)  iterativeWeights;
    Swarm*                       swarm           = (Swarm*) _swarm;
    double                       weightsTotal    = 0.0;
    Particle_InCellIndex         cParticle_I;
    Particle_InCellIndex         cParticleCount  = swarm->cellParticleCountTbl[lCell_I];
    IntegrationPoint*            particle;

    weightsTotal = WeightsCalculator_SumCellWeights( self, swarm, lCell_I );

    for ( cParticle_I = 0 ; cParticle_I < cParticleCount ; cParticle_I++ ) {
        particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, lCell_I, cParticle_I );

        /* Scale weights so that sum of weights = cellLocalVolume */
        particle->weight *= self->cellLocalVolume/weightsTotal;
    }
}
	



