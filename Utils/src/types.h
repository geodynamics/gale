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
**	Basic framework types.
**
** Assumptions:
**	None as yet.
**
** Comments:
**	None as yet.
**
** $Id: types.h 4160 2007-07-30 06:17:06Z DavidLee $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgDomain_Utils_types_h__
#define __StgDomain_Utils_types_h__

typedef struct NewRemesher NewRemesher;
typedef struct RegularRemesher RegularRemesher;
typedef struct RegularRemesherCmpt RegularRemesherCmpt;
	
	/* Classes */
	typedef struct _AllElementsVC			AllElementsVC;
	typedef struct _AllNodesVC			AllNodesVC;
	typedef struct _WallVC				WallVC;
	typedef struct _ContactVC ContactVC;
	typedef struct _InnerWallVC			InnerWallVC;
	typedef struct _CornerVC			CornerVC;
	typedef struct MeshShapeVC				MeshShapeVC;
	typedef struct MeshBoundaryShape MeshBoundaryShape;
	typedef struct _DofLayout			DofLayout;
	typedef struct _EmbeddedSurface			EmbeddedSurface;
	typedef struct Operator      		Operator;
	typedef struct FieldVariable			FieldVariable;
	typedef struct FieldVariable_Register		FieldVariable_Register;
	typedef struct OperatorFieldVariable		OperatorFieldVariable;
	typedef struct DomainContext		DomainContext;
	typedef struct LinearRegression     	LinearRegression;
	typedef struct SemiRegDeform		SemiRegDeform;
	typedef struct SobolGenerator            SobolGenerator;
	typedef struct TimeIntegrand            TimeIntegrand;
	typedef struct TimeIntegrator            TimeIntegrator;
	typedef struct ShapeAdvector             ShapeAdvector;
	typedef struct StripRemesher		StripRemesher;
	typedef struct CellRemesher		CellRemesher;

	/* Degrees of Freedom */
	typedef unsigned int				Dof;
	typedef Index                       Bit_Index;

	/* WallVC types */
	typedef enum
	{
		WallVC_Wall_Back,
		WallVC_Wall_Left,
		WallVC_Wall_Bottom,
		WallVC_Wall_Right,
		WallVC_Wall_Top,
		WallVC_Wall_Front,
		WallVC_Wall_BottomLeft,
		WallVC_Wall_BottomRight,
		WallVC_Wall_Size
	} WallVC_Wall;

	/* CornerVC types */
	typedef enum
	{
		CornerVC_Corner_BottomLeftFront,
		CornerVC_Corner_BottomRightFront,
		CornerVC_Corner_TopLeftFront,
		CornerVC_Corner_TopRightFront,
		CornerVC_Corner_BottomLeftBack,
		CornerVC_Corner_BottomRightBack,
		CornerVC_Corner_TopLeftBack,
		CornerVC_Corner_TopRightBack,
		CornerVC_Corner_Size
	} CornerVC_Corner;

	/* InnerWallVC types */
	typedef enum
	{
		InnerWallVC_InnerWall_Back,
		InnerWallVC_InnerWall_Left,
		InnerWallVC_InnerWall_Bottom,
		InnerWallVC_InnerWall_Right,
		InnerWallVC_InnerWall_Top,
		InnerWallVC_InnerWall_Front,
		InnerWallVC_InnerWall_Size
	} InnerWallVC_InnerWall;


	typedef struct _WallVC_Entry			WallVC_Entry;
	typedef Index					WallVC_Entry_Index;
	typedef Index					MeshShapeVC_Entry_Index;

	typedef struct _CornerVC_Entry			CornerVC_Entry;
	typedef Index					CornerVC_Entry_Index;

	typedef struct _InnerWallVC_Entry		InnerWallVC_Entry;
	typedef Index					InnerWallVC_Entry_Index;
	
	typedef struct _AllElementsVC_Entry		AllElementsVC_Entry;
	typedef Index					AllElementsVC_Entry_Index;
	
	typedef struct _AllNodesVC_Entry		AllNodesVC_Entry;
	typedef Index					AllNodesVC_Entry_Index;
	
	typedef struct _AllVC_Entry			AllVC_Entry;
	typedef Index					AllVC_Entry_Index;
	
	typedef Index					FieldVariable_Index;
	
#endif /* __StgDomain_Utils_types_h__ */
