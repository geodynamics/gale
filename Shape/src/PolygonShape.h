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
** $Id: PolygonShape.h 4054 2007-03-28 06:46:32Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgDomain_Shape_PolygonShapeClass_h__
#define __StgDomain_Shape_PolygonShapeClass_h__

	/* Textual name of this class */
	extern const Type PolygonShape_Type;

	/* PolygonShape information */
	#define __PolygonShape \
		/* General info */ \
		__Stg_Shape \
		/* Virtual Info */\
		\
		Coord_List              vertexList;    \
		Index                   vertexCount;   \
		XYZ                     start;        \
		XYZ                     end;          \
		Axis                    perpendicularAxis; \

	struct PolygonShape { __PolygonShape };
	
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/
	PolygonShape* PolygonShape_New(
		Name                                  name,
		Dimension_Index                       dim,
		XYZ                                   centre, 
		double                                alpha,
		double                                beta,
		double                                gamma,
		Coord_List                            vertexList,
		Index                                 vertexCount,
		XYZ                                   start,
		XYZ                                   end,
		Axis                                  perpendicularAxis);
		
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define POLYGONSHAPE_DEFARGS \
                STG_SHAPE_DEFARGS

	#define POLYGONSHAPE_PASSARGS \
                STG_SHAPE_PASSARGS

	PolygonShape* _PolygonShape_New(  POLYGONSHAPE_DEFARGS  );
	
	void _PolygonShape_Init( void* polygon, Coord_List vertexList, Index vertexCount, XYZ start, XYZ end, Axis perpendicular ) ;
		
	/* Stg_Class_Delete PolygonShape implementation */
	void _PolygonShape_Delete( void* polygon );
	void _PolygonShape_Print( void* polygon, Stream* stream );
	#define PolygonShape_Copy( self ) \
		(PolygonShape*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define PolygonShape_DeepCopy( self ) \
		(PolygonShape*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _PolygonShape_Copy( const void* polygon, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	void* _PolygonShape_DefaultNew( Name name ) ;
	void _PolygonShape_AssignFromXML( void* shape, Stg_ComponentFactory* cf, void* data ) ;
	void _PolygonShape_Build( void* polygon, void* data ) ;
	void _PolygonShape_Initialise( void* polygon, void* data ) ;
	void _PolygonShape_Execute( void* polygon, void* data );
	void _PolygonShape_Destroy( void* polygon, void* data ) ;
	
	Bool _PolygonShape_IsCoordInside( void* polygon, const Coord coord ) ;
	double _PolygonShape_CalculateVolume( void* polygon );
	void _PolygonShape_DistanceFromCenterAxis( void* shape, const Coord coord, double* disVec );

	/*---------------------------------------------------------------------------------------------------------------------
	** Public member functions
	*/
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/
	
	
#endif 

