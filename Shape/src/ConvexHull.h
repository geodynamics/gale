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
** $Id: ConvexHull.h 3851 2006-10-12 08:57:22Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StGermain_Domain_Shape_ConvexHullClass_h__
#define __StGermain_Domain_Shape_ConvexHullClass_h__

	/* Textual name of this class */
	extern const Type ConvexHull_Type;

	/* ConvexHull information */
	#define __ConvexHull \
		/* General info */ \
		__Stg_Shape \
		/* Virtual Info */\
		\
		Coord_List              vertexList;    \
		Index                   vertexCount;   \
		XYZ*                    facesList;     \

	struct ConvexHull { __ConvexHull };
	
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/
	ConvexHull* ConvexHull_New(
		Name                                  name,
		Dimension_Index                       dim,
		XYZ                                   centre, 
		double                                alpha,
		double                                beta,
		double                                gamma,
		Coord_List                            vertexList,
		Index                                 vertexCount);
		
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define CONVEXHULL_DEFARGS \
                STG_SHAPE_DEFARGS

	#define CONVEXHULL_PASSARGS \
                STG_SHAPE_PASSARGS

	ConvexHull* _ConvexHull_New(  CONVEXHULL_DEFARGS  );
	
	void _ConvexHull_Init( void* convexHull, Coord_List vertexList, Index vertexCount);
		
	/* Stg_Class_Delete ConvexHull implementation */
	void _ConvexHull_Delete( void* convexHull );
	void _ConvexHull_Print( void* convexHull, Stream* stream );
	#define ConvexHull_Copy( self ) \
		(ConvexHull*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define ConvexHull_DeepCopy( self ) \
		(ConvexHull*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _ConvexHull_Copy( void* convexHull, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	void* _ConvexHull_DefaultNew( Name name ) ;
	void _ConvexHull_AssignFromXML( void* shape, Stg_ComponentFactory* cf, void* data ) ;
	void _ConvexHull_Build( void* convexHull, void* data ) ;
	void _ConvexHull_Initialise( void* convexHull, void* data ) ;
	void _ConvexHull_Execute( void* convexHull, void* data );
	void _ConvexHull_Destroy( void* convexHull, void* data ) ;
	
	Bool _ConvexHull_IsCoordInside( void* convexHull, Coord coord ) ;
	double _ConvexHull_CalculateVolume( void* convexHull );
	void _ConvecHull_DistanceFromCenterAxis( void* self, Coord coord, double* disVec );

	/*---------------------------------------------------------------------------------------------------------------------
	** Public member functions
	*/
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/
	
	
#endif 

