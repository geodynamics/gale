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
** $Id: Cylinder.h 3851 2006-10-12 08:57:22Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgDomain_Shape_CylinderClass_h__
#define __StgDomain_Shape_CylinderClass_h__

	/* Textual name of this class */
	extern const Type Cylinder_Type;

	/* Cylinder information */
	#define __Cylinder \
		/* General info */ \
		__Stg_Shape \
		/* Virtual Info */\
		\
		double                                radius;           \
		XYZ                                   start;            \
		XYZ                                   end;              \
		Axis                                  alongAxis;

	struct Cylinder { __Cylinder };
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/
	Cylinder* Cylinder_New(
		Name                                  name,
		Dimension_Index                       dim,
		XYZ                                   centre, 
		double                                alpha,
		double                                beta,
		double                                gamma,
		double                                radius, 
		XYZ                                   start, 
		XYZ                                   end, 
		Axis                                  perpendicularAxis );

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define CYLINDER_DEFARGS \
                STG_SHAPE_DEFARGS

	#define CYLINDER_PASSARGS \
                STG_SHAPE_PASSARGS

	Cylinder* _Cylinder_New(  CYLINDER_DEFARGS  );
	
	void _Cylinder_Init( Cylinder* self, double radius, XYZ start, XYZ end, Axis perpendicularAxis ) ;

	/* Stg_Class_Delete Cylinder implementation */
	void _Cylinder_Delete( void* cylinder );
	void _Cylinder_Print( void* cylinder, Stream* stream );
	#define Cylinder_Copy( self ) \
		(Cylinder*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define Cylinder_DeepCopy( self ) \
		(Cylinder*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _Cylinder_Copy( const void* cylinder, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	void* _Cylinder_DefaultNew( Name name ) ;
	void _Cylinder_AssignFromXML( void* shape, Stg_ComponentFactory* cf, void* data ) ;
	void _Cylinder_Build( void* cylinder, void* data ) ;
	void _Cylinder_Initialise( void* cylinder, void* data ) ;
	void _Cylinder_Execute( void* cylinder, void* data );
	void _Cylinder_Destroy( void* cylinder, void* data ) ;
	
	Bool _Cylinder_IsCoordInside( void* cylinder, Coord coord ) ;
	void _Cylinder_DistanceFromCenterAxis( void* cylinder, Coord coord, double* disVec );
	double _Cylinder_CalculateVolume( void* cylinder );

	/*---------------------------------------------------------------------------------------------------------------------
	** Public member functions
	*/
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/
	
	
#endif 

