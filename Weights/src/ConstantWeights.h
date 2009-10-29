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
*/
/** \file
 **  Role:
 **
 ** Assumptions:
 **
 ** Comments:
 **
 ** $Id: ConstantWeights.h 374 2006-10-12 08:59:41Z SteveQuenette $
 **
 **~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __PICellerator_Weights_ConstantWeightsClass_h__
#define __PICellerator_Weights_ConstantWeightsClass_h__

/* Textual name of this class */
extern const Type ConstantWeights_Type;

/* ConstantWeights information */
#define __ConstantWeights                       \
    /* General info */                          \
    __WeightsCalculator                         \
                                                \
    /* Virtual Info */                          \
                                                \

struct ConstantWeights { __ConstantWeights };

        
/*---------------------------------------------------------------------------------------------------------------------
** Constructors
*/

#define CONSTANTWEIGHTS_DEFARGS \
    WEIGHTSCALCULATOR_DEFARGS

#define CONSTANTWEIGHTS_PASSARGS \
    WEIGHTSCALCULATOR_PASSARGS

ConstantWeights* ConstantWeights_New( Name name, int dim );

ConstantWeights* _ConstantWeights_New( CONSTANTWEIGHTS_DEFARGS );

void _ConstantWeights_Init( void* constantWeights );


/* Stg_Class_Delete ConstantWeights implementation */
void _ConstantWeights_Delete( void* constantWeights );
void _ConstantWeights_Print( void* constantWeights, Stream* stream );
#define ConstantWeights_Copy( self )                                    \
    (ConstantWeights*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
#define ConstantWeights_DeepCopy( self )                                \
    (ConstantWeights*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
void* _ConstantWeights_Copy( void* constantWeights, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
        
void* _ConstantWeights_DefaultNew( Name name ) ;
void _ConstantWeights_Construct( void* constantWeights, Stg_ComponentFactory* cf, void* data ) ;
void _ConstantWeights_Build( void* constantWeights, void* data ) ;
void _ConstantWeights_Initialise( void* constantWeights, void* data ) ;
void _ConstantWeights_Execute( void* constantWeights, void* data );
        
                
void _ConstantWeights_Calculate( void* constantWeights, void* _swarm, Cell_LocalIndex lCell_I ) ;
/*---------------------------------------------------------------------------------------------------------------------
** Private functions
*/
        
/*---------------------------------------------------------------------------------------------------------------------
** Public functions
*/
        
        
#endif 
