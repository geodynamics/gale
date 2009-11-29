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
** $Id: Everywhere.h 3851 2006-10-12 08:57:22Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StGermain_Domain_Shape_EverywhereClass_h__
#define __StGermain_Domain_Shape_EverywhereClass_h__

	/* Textual name of this class */
	extern const Type Everywhere_Type;

	/* Everywhere information */
	#define __Everywhere \
		/* General info */ \
		__Stg_Shape \
		/* Virtual Info */\
		\

	struct Everywhere { __Everywhere };
	
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/
	Everywhere* Everywhere_New(
		Name                                  name,
		Dimension_Index                       dim );

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define EVERYWHERE_DEFARGS \
                STG_SHAPE_DEFARGS

	#define EVERYWHERE_PASSARGS \
                STG_SHAPE_PASSARGS

	Everywhere* _Everywhere_New(  EVERYWHERE_DEFARGS  );
	
	void _Everywhere_Init( void* everywhere ) ;

	/* Stg_Class_Delete Everywhere implementation */
	void _Everywhere_Delete( void* everywhere );
	void _Everywhere_Print( void* everywhere, Stream* stream );
	#define Everywhere_Copy( self ) \
		(Everywhere*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define Everywhere_DeepCopy( self ) \
		(Everywhere*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _Everywhere_Copy( void* everywhere, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	void* _Everywhere_DefaultNew( Name name ) ;
	void _Everywhere_AssignFromXML( void* shape, Stg_ComponentFactory* cf, void* data ) ;
	void _Everywhere_Build( void* everywhere, void* data ) ;
	void _Everywhere_Initialise( void* everywhere, void* data ) ;
	void _Everywhere_Execute( void* everywhere, void* data );
	void _Everywhere_Destroy( void* everywhere, void* data ) ;
	
	Bool _Everywhere_IsCoordInside( void* everywhere, Coord coord ) ;
	double _Everywhere_CalculateVolume( void* everywhere );
	void _Everywhere_DistanceFromCenterAxis( void* shape, Coord coord, double* disVec );

	/*---------------------------------------------------------------------------------------------------------------------
	** Public member functions
	*/
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/
	
#endif 

