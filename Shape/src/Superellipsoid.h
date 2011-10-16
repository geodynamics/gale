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
** $Id: Superellipsoid.h 3851 2006-10-12 08:57:22Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgDomain_Shape_SuperellipsoidClass_h__
#define __StgDomain_Shape_SuperellipsoidClass_h__

	/* Textual name of this class */
	extern const Type Superellipsoid_Type;

	/* Superellipsoid information */
	#define __Superellipsoid \
		/* General info */ \
		__Stg_Shape \
		/* Virtual Info */\
		\
		double                                epsilon1;            \
		double                                epsilon2;            \
		XYZ                                   radius;		   \

	struct Superellipsoid { __Superellipsoid };
	
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/
	Superellipsoid* Superellipsoid_New(
		Name                                  name,
		Dimension_Index                       dim,
		XYZ                                   centre, 
		double                                alpha,
		double                                beta,
		double                                gamma,
		double                                epsilon1,
		double                                epsilon2,   
		XYZ                                   radius );

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define SUPERELLIPSOID_DEFARGS \
                STG_SHAPE_DEFARGS

	#define SUPERELLIPSOID_PASSARGS \
                STG_SHAPE_PASSARGS

	Superellipsoid* _Superellipsoid_New(  SUPERELLIPSOID_DEFARGS  );
	
	void _Superellipsoid_Init( void* superellipsoid, double epsilon1, double epsilon2, XYZ radius ) ;

	/* Stg_Class_Delete Superellipsoid implementation */
	void _Superellipsoid_Delete( void* superellipsoid );
	void _Superellipsoid_Print( void* superellipsoid, Stream* stream );
	#define Superellipsoid_Copy( self ) \
		(Superellipsoid*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define Superellipsoid_DeepCopy( self ) \
		(Superellipsoid*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _Superellipsoid_Copy( const void* superellipsoid, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	void* _Superellipsoid_DefaultNew( Name name ) ;
	void _Superellipsoid_AssignFromXML( void* shape, Stg_ComponentFactory* cf, void* data ) ;
	void _Superellipsoid_Build( void* superellipsoid, void* data ) ;
	void _Superellipsoid_Initialise( void* superellipsoid, void* data ) ;
	void _Superellipsoid_Execute( void* superellipsoid, void* data );
	void _Superellipsoid_Destroy( void* superellipsoid, void* data ) ;
	
	Bool _Superellipsoid_IsCoordInside( void* superellipsoid, const Coord coord ) ;
	double _Superellipsoid_CalculateVolume( void* superellipsoid );
	void _Superellipsoid_DistanceFromCenterAxis( void* shape, const Coord coord, double* disVec );

	/*---------------------------------------------------------------------------------------------------------------------
	** Public member functions
	*/
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/
	
	
#endif 

