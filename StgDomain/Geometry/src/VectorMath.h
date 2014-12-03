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
**	Robert B. Turnbull, Monash Cluster Computing. (Robert.Turnbull@sci.monash.edu.au)
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
** Role:
**    Provides basic vector operations.
**
** Assumptions:
**    - Coord is an array of 3 doubles.
**
** Comments:
**    In any operation that involves two or more input vectors, those vectors 
**    may be one and the same; it may be assumed that such an occurence will be
**    handled.
**
** $Id: VectorMath.h 4081 2007-04-27 06:20:07Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgDomain_Geometry_VectorMath_h__
#define __StgDomain_Geometry_VectorMath_h__

	/*--------------------------------------------------------------------------------------------------------------------------
	** Macros
	*/

	#define Vec_Set2D( dst, src )			\
		((dst)[0] = (src)[0],			\
		 (dst)[1] = (src)[1], 0)

	#define Vec_Set3D( dst, src )			\
		(Vec_Set2D( dst, src ),			\
		 (dst)[2] = (src)[2], 0)

	#define Vec_SetScalar2D( dst, x, y )		\
		((dst)[0] = x,				\
		 (dst)[1] = y, 0)

	#define Vec_SetScalar3D( dst, x, y, z )		\
		(Vec_SetScalar2D( dst, x, y ),		\
		 (dst)[2] = z, 0)

	#define Vec_Add2D( dst, a, b )			\
		((dst)[0] = (a)[0] + (b)[0],		\
		 (dst)[1] = (a)[1] + (b)[1], 0)

	#define Vec_Add3D( dst, a, b )			\
		(Vec_Add2D( dst, a, b ),		\
		 (dst)[2] = (a)[2] + (b)[2], 0)

	#define Vec_Sub2D( dst, a, b )			\
		((dst)[0] = (a)[0] - (b)[0],		\
		 (dst)[1] = (a)[1] - (b)[1], 0)

	#define Vec_Sub3D( dst, a, b )			\
		(Vec_Sub2D( dst, a, b ),		\
		 (dst)[2] = (a)[2] - (b)[2], 0)

	#define Vec_Dot2D( a, b )			\
		((a)[0] * (b)[0] + (a)[1] * (b)[1])

	#define Vec_Dot3D( a, b )			\
		(Vec_Dot2D( a, b ) + (a)[2] * (b)[2])

	#define Vec_Scale2D( dst, a, s )		\
		((dst)[0] = (a)[0] * (s),		\
		 (dst)[1] = (a)[1] * (s), 0)

	#define Vec_Scale3D( dst, a, s )		\
		(Vec_Scale2D( dst, a, s ),		\
		 (dst)[2] = (a)[2] * (s), 0)

	#define Vec_MagSq2D( a )			\
		Vec_Dot2D( a, a )

	#define Vec_MagSq3D( a )			\
		Vec_Dot3D( a, a )

	#define Vec_Mag2D( a )				\
		sqrt( Vec_MagSq2D( a ) )

	#define Vec_Mag3D( a )				\
		sqrt( Vec_MagSq3D( a ) )

	#define Vec_Proj2D( dst, a, b )					\
		(Vec_Norm2D( dst, b ),					\
		 Vec_Scale2D( dst, dst, Vec_Dot2D( a, b ) ), 0)

	#define Vec_Proj3D( dst, a, b )					\
		(Vec_Norm3D( dst, b ),					\
		 Vec_Scale3D( dst, dst, Vec_Dot3D( a, b ) ), 0)


	/*--------------------------------------------------------------------------------------------------------------------------
	** Functions
	*/

	void Vec_Cross3D( double* dst, const double* a, const double* b );
	void Vec_Div2D( double* dst, const double* a, const double s );
	void Vec_Div3D( double* dst, const double* a, const double s );
	void Vec_Norm2D( double* dst, const double* a );
	void Vec_Norm3D( double* dst, const double* a );

	void StGermain_RotateVector(double* rotatedVector, const double* vector, const double* w, const double theta) ;
	void StGermain_RotateCoordinateAxis( double* rotatedVector, const double* vector, Index axis, const double theta ) ;
	void StGermain_VectorSubtraction(double* destination, const double* vector1, const double* vector2, Index dim) ;
	void StGermain_VectorAddition(double* destination, const double* vector1, const double* vector2, Index dim) ;
	double StGermain_VectorMagnitude(const double* vector, Index dim) ;
	double StGermain_VectorDotProduct(const double* vector1, const double* vector2, Index dim) ;
	void StGermain_VectorCrossProduct(double* destination, const double* vector1, const double* vector2) ;
	double StGermain_VectorCrossProductMagnitude( const double* vector1, const double* vector2, Dimension_Index dim ) ;
	double StGermain_ScalarTripleProduct( const double* vectorA, const double* vectorB, const double* vectorC ) ;

	void StGermain_VectorNormalise(double* vector, Index dim) ;
	double StGermain_AngleBetweenVectors( const double* vectorA, const double* vectorB, Index dim ) ;
	double StGermain_DistanceBetweenPoints( const double* pos1, const double* pos2, Index dim) ;
	void StGermain_NormalToPlane( double* normal, const double* pos0, const double* pos1, const double* pos2) ;

	void StGermain_TriangleCentroid( double* centroid, const double* pos0, const double* pos1, const double* pos2, Index dim) ;
	double StGermain_TriangleArea( const double* pos0, const double* pos1, const double* pos2, Index dim ) ;
	double StGermain_ConvexQuadrilateralArea( const double* vertexCoord1, const double* vertexCoord2, 
						  const double* vertexCoord3, const double* vertexCoord4, 
						  Dimension_Index dim ) ;
	double StGermain_ParallelepipedVolume( const double* coordLeftBottomFront, 
					       const double* coordRightBottomFront, 
					       const double* coordLeftTopFront, 
					       const double* coordLeftBottomBack );
	double StGermain_ParallelepipedVolumeFromCoordList( Coord_List list ) ;
	
	void StGermain_AverageCoord( double* coord, double** coordList, Index count, Dimension_Index dim ) ;
	void StGermain_PrintVector( Stream* stream, const double* vector, Index dim ) ;

	/** Print a named vector. Name comes from vector variable in file*/
	#define StGermain_PrintNamedVector(stream, vector, dim)		\
		do {							\
			Journal_Printf( stream, #vector " - " );	\
			StGermain_PrintVector( stream, vector, dim );	\
		} while(0)

#endif /* __StgDomain_Geometry_VectorMath_h__ */
