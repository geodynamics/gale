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

	typedef void (WeightsCalculator_CalculateFunction)( void* weightsCalculator, void* swarm, Cell_LocalIndex lCell_I );

	/* Textual name of this class */
	extern const Type WeightsCalculator_Type;

	/* WeightsCalculator information */
	#define __WeightsCalculator \
		/* General info */ \
		__Stg_Component \
		/* Virtual Info */\
		FiniteElementContext*		      context;		 \
		WeightsCalculator_CalculateFunction*  _calculate;        \
		/* Other Info */\
		double                                cellLocalVolume;   \
		Dimension_Index                       dim;

	struct WeightsCalculator { __WeightsCalculator };
	
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/
	WeightsCalculator* _WeightsCalculator_New(
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
		Name                                  name );
	
	void _WeightsCalculator_Init( void* weightsCalculator , Dimension_Index dim) ;
	void WeightsCalculator_InitAll( void* weightsCalculator, Dimension_Index dim );

	/* Stg_Class_Delete WeightsCalculator implementation */
	void _WeightsCalculator_Delete( void* weightsCalculator );
	void _WeightsCalculator_Print( void* weightsCalculator, Stream* stream );
	#define WeightsCalculator_Copy( self ) \
		(WeightsCalculator*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define WeightsCalculator_DeepCopy( self ) \
		(WeightsCalculator*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _WeightsCalculator_Copy( void* weightsCalculator, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
void _WeightsCalculator_Construct( void* shape, Stg_ComponentFactory* cf, void* data ) ;
	void _WeightsCalculator_Build( void* weightsCalculator, void* data ) ;
	void _WeightsCalculator_Initialise( void* weightsCalculator, void* data ) ;
	void _WeightsCalculator_Execute( void* weightsCalculator, void* data );
	void _WeightsCalculator_Destroy( void* weightsCalculator, void* data ) ;
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Public member functions
	*/
	void WeightsCalculator_CalculateCell( void* weightsCalculator, void* swarm, Cell_LocalIndex lCell_I ) ;
	
	void WeightsCalculator_CalculateAll( void* weightsCalculator, void* _swarm ) ;
	/*---------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	#define WeightsCalculator_ZeroWeights( self, swarm ) \
		WeightsCalculator_SetWeightsValueAll( (self), (swarm), 0.0 )
	#define WeightsCalculator_ZeroWeightsInCell( self, swarm, lCell_I ) \
		WeightsCalculator_SetWeightsValueAllInCell( (self), (swarm), (lCell_I), 0.0 )
	
	void WeightsCalculator_SetWeightsValueAll( void* weightsCalculator, void* _swarm, double weight ) ;
	void WeightsCalculator_SetWeightsValueAllInCell( void* weightsCalculator, void* _swarm, Cell_LocalIndex lCell_I, double weight ) ;
	Constraint_Index WeightsCalculator_FindConstraintOrder( void* weightsCalculator, void* _swarm, Dimension_Index dim, Stream* stream ) ;
	double WeightsCalculator_TestConstraint( void* weightsCalculator, void* _swarm, Dimension_Index dim, Constraint_Index order ) ;
	double WeightsCalculator_TestConstraintOverCell( void* weightsCalculator, void* _swarm, Cell_LocalIndex lCell_I, Dimension_Index dim, Constraint_Index order ) ;
	double WeightsCalculator_GetConstraintLHS( void* weightsCalculator, void* _swarm, Cell_LocalIndex lCell_I, Index power_i, Index power_j, Index power_k ) ;
	double WeightsCalculator_GetLocalCoordSum( void* weightsCalculator, void* _swarm, Cell_LocalIndex lCell_I, Index power_i, Index power_j, Index power_k ) ;
	double WeightsCalculator_SumCellWeights( void* weightsCalculator, void* _swarm, Cell_LocalIndex lCell_I ) ;
	void WeightsCalculator_CheckEmptyCell( void* weightsCalculator, void* _swarm, Cell_LocalIndex lCell_I ) ;
	
#endif 
