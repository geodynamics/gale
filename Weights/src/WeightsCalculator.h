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
 ** $Id: WeightsCalculator.h 374 2006-10-12 08:59:41Z SteveQuenette $
 **
 **~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __PICellerator_Weights_WeightsCalculator_h__
#define __PICellerator_Weights_WeightsCalculator_h__

typedef void (WeightsCalculator_CalculateFunction)( void* self, void* swarm, Cell_LocalIndex lCell_I );

/* Textual name of this class */
extern const Type WeightsCalculator_Type;

/* WeightsCalculator information */
#define __WeightsCalculator \
	/* General info */ \
	__Stg_Component \
	\
	/* Virtual Info */ \
	FiniteElementContext*						context; \
	WeightsCalculator_CalculateFunction*  _calculate; \
	/* Other Info */ \
	double                                cellLocalVolume; \
	Dimension_Index                       dim;

	struct WeightsCalculator { __WeightsCalculator };
        
	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define WEIGHTSCALCULATOR_DEFARGS \
	    STG_COMPONENT_DEFARGS, \
	        WeightsCalculator_CalculateFunction* _calculate, \
	        int dim

	#define WEIGHTSCALCULATOR_PASSARGS \
	    STG_COMPONENT_PASSARGS, \
	        _calculate, \
				 dim

	WeightsCalculator* _WeightsCalculator_New( WEIGHTSCALCULATOR_DEFARGS );

	void _WeightsCalculator_Init( void* self , int dim );

	/* Stg_Class_Delete WeightsCalculator implementation */
	void _WeightsCalculator_Delete( void* self );
	void _WeightsCalculator_Print( void* self, Stream* stream );
	#define WeightsCalculator_Copy( self )                                  \
	    (WeightsCalculator*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define WeightsCalculator_DeepCopy( self )                              \
	    (WeightsCalculator*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _WeightsCalculator_Copy( void* self, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
        
	void _WeightsCalculator_AssignFromXML( void* self, Stg_ComponentFactory* cf, void* data ) ;
	void _WeightsCalculator_Build( void* self, void* data ) ;
	void _WeightsCalculator_Initialise( void* self, void* data ) ;
	void _WeightsCalculator_Execute( void* self, void* data );
        
	/*---------------------------------------------------------------------------------------------------------------------
	** Public member functions
	*/
	void WeightsCalculator_CalculateCell( void* self, void* swarm, Cell_LocalIndex lCell_I ) ;
        
	void WeightsCalculator_CalculateAll( void* self, void* _swarm ) ;
	/*---------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	#define WeightsCalculator_ZeroWeights( self, swarm )                    \
	    WeightsCalculator_SetWeightsValueAll( (self), (swarm), 0.0 )
	#define WeightsCalculator_ZeroWeightsInCell( self, swarm, lCell_I )     \
	    WeightsCalculator_SetWeightsValueAllInCell( (self), (swarm), (lCell_I), 0.0 )
        
	void WeightsCalculator_SetWeightsValueAll( void* self, void* _swarm, double weight ) ;
	void WeightsCalculator_SetWeightsValueAllInCell( void* self, void* _swarm, Cell_LocalIndex lCell_I, double weight ) ;
	Constraint_Index WeightsCalculator_FindConstraintOrder( void* self, void* _swarm, Dimension_Index dim, Stream* stream ) ;
	double WeightsCalculator_TestConstraint( void* self, void* _swarm, Dimension_Index dim, Constraint_Index order ) ;
	double WeightsCalculator_TestConstraintOverCell( void* self, void* _swarm, Cell_LocalIndex lCell_I, Dimension_Index dim, Constraint_Index order ) ;
	double WeightsCalculator_GetConstraintLHS( void* self, void* _swarm, Cell_LocalIndex lCell_I, Index power_i, Index power_j, Index power_k ) ;
	double WeightsCalculator_GetLocalCoordSum( void* self, void* _swarm, Cell_LocalIndex lCell_I, Index power_i, Index power_j, Index power_k ) ;
	double WeightsCalculator_SumCellWeights( void* self, void* _swarm, Cell_LocalIndex lCell_I ) ;
	void WeightsCalculator_CheckEmptyCell( void* self, void* _swarm, Cell_LocalIndex lCell_I ) ;
        
#endif 
