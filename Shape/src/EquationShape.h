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
** $Id: EquationShape.h 3851 2006-10-12 08:57:22Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgDomain_Shape_EquationShapeClass_h__
#define __StgDomain_Shape_EquationShapeClass_h__

	/* Textual name of this class */
	extern const Type EquationShape_Type;

	/* EquationShape information */
	#define __EquationShape \
		/* General info */ \
		__Stg_Shape \
		/* Virtual Info */\
		\
		char *equation;

	struct EquationShape { __EquationShape };
	
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/
	EquationShape* EquationShape_New(
		Name                                  name,
		Dimension_Index                       dim,
		XYZ                                   centre,
		double                                alpha,
		double                                beta,
		double                                gamma,
                char*                                 equation);

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define EQUATIONSHAPE_DEFARGS \
                STG_SHAPE_DEFARGS

	#define EQUATIONSHAPE_PASSARGS \
                STG_SHAPE_PASSARGS

	EquationShape* _EquationShape_New(  EQUATIONSHAPE_DEFARGS  );
	
	void _EquationShape_Init( void* sphere, char *equation ) ;

	/* Stg_Class_Delete EquationShape implementation */
	void _EquationShape_Delete( void* sphere );
	void _EquationShape_Print( void* sphere, Stream* stream );
	#define EquationShape_Copy( self ) \
		(EquationShape*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define EquationShape_DeepCopy( self ) \
		(EquationShape*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _EquationShape_Copy( const void* sphere, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	void* _EquationShape_DefaultNew( Name name ) ;
	void _EquationShape_AssignFromXML( void* shape, Stg_ComponentFactory* cf, void* data ) ;
	void _EquationShape_Build( void* sphere, void* data ) ;
	void _EquationShape_Initialise( void* sphere, void* data ) ;
	void _EquationShape_Execute( void* sphere, void* data );
	void _EquationShape_Destroy( void* sphere, void* data ) ;
	
	Bool _EquationShape_IsCoordInside( void* sphere, const Coord coord ) ;
	void _EquationShape_DistanceFromCenterAxis( void* sphere, const Coord coord, double* disVec );
	double _EquationShape_CalculateVolume( void* sphere );

	/*---------------------------------------------------------------------------------------------------------------------
	** Public member functions
	*/
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/
	
	
#endif 

