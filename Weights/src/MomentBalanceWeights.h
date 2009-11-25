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
 ** $Id: MomentBalanceWeights.h 374 2006-10-12 08:59:41Z SteveQuenette $
 **
 **~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __PICellerator_Weights_MomentBalanceWeightsClass_h__
#define __PICellerator_Weights_MomentBalanceWeightsClass_h__

/* Textual name of this class */
extern const Type MomentBalanceWeights_Type;

	/* MomentBalanceWeights information */
	#define __MomentBalanceWeights \
		/* General info */ \
		__WeightsCalculator \
		\
		/* Virtual Info */ \
		WeightsCalculator*	backupWeights; \
		Bool						freeBackupWeights;

	struct MomentBalanceWeights { __MomentBalanceWeights };

	#define MOMENTBALANCEWEIGHTS_DEFARGS \
		WEIGHTSCALCULATOR_DEFARGS, \
			WeightsCalculator*	backupWeights

	#define MOMENTBALANCEWEIGHTS_PASSARGS \
		WEIGHTSCALCULATOR_PASSARG, \
			backupWeights

	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	MomentBalanceWeights* MomentBalanceWeights_New( Name name, Dimension_Index dim, WeightsCalculator* backupWeights );

	MomentBalanceWeights* _MomentBalanceWeights_New( MOMENTBALANCEWEIGHTS_DEFARGS );

	/* Stg_Class implementation */
	void _MomentBalanceWeights_Delete( void* momentBalanceWeights );

	void _MomentBalanceWeights_Print( void* momentBalanceWeights, Stream* stream );

	#define MomentBalanceWeights_Copy( self ) \
		(MomentBalanceWeights*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define MomentBalanceWeights_DeepCopy( self )                           \
		(MomentBalanceWeights*) Stg_Class_Copy( self, NULL, True, NULL, NULL )

	void _MomentBalanceWeights_Init( void* momentBalanceWeights, WeightsCalculator* backupWeights );

	void* _MomentBalanceWeights_Copy( void* momentBalanceWeights, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
        
	void* _MomentBalanceWeights_DefaultNew( Name name );

	void _MomentBalanceWeights_AssignFromXML( void* momentBalanceWeights, Stg_ComponentFactory* cf, void* data );

	void _MomentBalanceWeights_Build( void* momentBalanceWeights, void* data );

	void _MomentBalanceWeights_Initialise( void* momentBalanceWeights, void* data );

	void _MomentBalanceWeights_Destroy( void* momentBalanceWeights, void* data );

	void _MomentBalanceWeights_Execute( void* momentBalanceWeights, void* data );
                
	void _MomentBalanceWeights_Calculate( void* momentBalanceWeights, void* _swarm, Cell_LocalIndex lCell_I );

/*---------------------------------------------------------------------------------------------------------------------
** Private functions
*/
        
/*---------------------------------------------------------------------------------------------------------------------
** Public functions
*/
        
        
#endif 
