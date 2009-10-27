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
** $Id: ConstantWeights.c 518 2007-10-11 08:07:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include "types.h"
#include "WeightsCalculator.h"
#include "ConstantWeights.h"

#include <assert.h>
#include <string.h>
#include <math.h>


/* Textual name of this class */
const Type ConstantWeights_Type = "ConstantWeights";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

ConstantWeights* ConstantWeights_New( Name name, int dim ) {
    return _ConstantWeights_New(
        sizeof(ConstantWeights),
        ConstantWeights_Type,
        _ConstantWeights_Delete,
        _ConstantWeights_Print,
        _ConstantWeights_Copy,
        _ConstantWeights_DefaultNew,
        _ConstantWeights_Construct,
        _ConstantWeights_Build,
        _ConstantWeights_Initialise,
        _ConstantWeights_Execute,
        _ConstantWeights_Destroy,
        _ConstantWeights_Calculate,
        name,
        True,
        dim );
}

ConstantWeights* _ConstantWeights_New(
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
    int                                   dim )
{
    ConstantWeights* self;

    /* Allocate memory */
    assert( _sizeOfSelf >= sizeof(ConstantWeights) );
    self = (ConstantWeights*)_WeightsCalculator_New( 
        _sizeOfSelf,
        type,
        _delete,
        _print,
        _copy,
        _defaultConstructor,
        _construct,
        _build,
        _initialise,
        _execute,
        _destroy,		
        _calculate,
        name,
        initFlag,
        dim );

	
    /* General info */

    /* Virtual Info */

    if(initFlag)
        _ConstantWeights_Init( self );

    return self;
}

void _ConstantWeights_Init( void* constantWeights  ) {
    ConstantWeights* self = (ConstantWeights*)constantWeights;
}


/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _ConstantWeights_Delete( void* constantWeights ) {
    ConstantWeights* self = (ConstantWeights*)constantWeights;
	
    /* Delete parent */
    _WeightsCalculator_Delete( self );
}


void _ConstantWeights_Print( void* constantWeights, Stream* stream ) {
    ConstantWeights* self = (ConstantWeights*)constantWeights;
	
    /* Print parent */
    _WeightsCalculator_Print( self, stream );
}



void* _ConstantWeights_Copy( void* constantWeights, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
    ConstantWeights*	self = (ConstantWeights*)constantWeights;
    ConstantWeights*	newConstantWeights;
	
    newConstantWeights = (ConstantWeights*)_WeightsCalculator_Copy( self, dest, deep, nameExt, ptrMap );
	
    return (void*)newConstantWeights;
}

void* _ConstantWeights_DefaultNew( Name name ) {
    return (void*) _ConstantWeights_New(
        sizeof(ConstantWeights),
        ConstantWeights_Type,
        _ConstantWeights_Delete,
        _ConstantWeights_Print,
        _ConstantWeights_Copy,
        _ConstantWeights_DefaultNew,
        _ConstantWeights_Construct,
        _ConstantWeights_Build,
        _ConstantWeights_Initialise,
        _ConstantWeights_Execute,
        _ConstantWeights_Destroy,
        _ConstantWeights_Calculate,
        name,
        False,
        0 );
}


void _ConstantWeights_Construct( void* constantWeights, Stg_ComponentFactory* cf, void* data ) {
    ConstantWeights*	     self          = (ConstantWeights*) constantWeights;

    _WeightsCalculator_Construct( self, cf, data );
	
    _ConstantWeights_Init( self );
}

void _ConstantWeights_Build( void* constantWeights, void* data ) {
    ConstantWeights*	self = (ConstantWeights*)constantWeights;

    _WeightsCalculator_Build( self, data );
}
void _ConstantWeights_Initialise( void* constantWeights, void* data ) {
    ConstantWeights*	self = (ConstantWeights*)constantWeights;
	
    _WeightsCalculator_Initialise( self, data );
}
void _ConstantWeights_Execute( void* constantWeights, void* data ) {
    ConstantWeights*	self = (ConstantWeights*)constantWeights;
	
    _WeightsCalculator_Execute( self, data );
}
void _ConstantWeights_Destroy( void* constantWeights, void* data ) {
    ConstantWeights*	self = (ConstantWeights*)constantWeights;
	
    _WeightsCalculator_Destroy( self, data );
}

/*-------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
void _ConstantWeights_Calculate( void* constantWeights, void* _swarm, Cell_LocalIndex lCell_I ) {
    ConstantWeights*             self            = (ConstantWeights*)  constantWeights;
    Swarm*                       swarm           = (Swarm*) _swarm;
    double                       weight;
    Particle_InCellIndex         cParticleCount;
		
    cParticleCount = swarm->cellParticleCountTbl[lCell_I];
    weight = self->cellLocalVolume / (double) cParticleCount;
    WeightsCalculator_SetWeightsValueAllInCell( self, swarm, lCell_I, weight );
}

/*-------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/


