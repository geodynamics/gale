/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd,
** 110 Victoria Street, Melbourne, 3053, Australia.
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
**  Role:
**
** Assumptions:
**
** Invariants:
**
** Comments:
**
** $Id: MeshBoundaryShape.h 2225 1970-01-02 13:48:23Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgDomain_Utils_MeshBoundaryShape_h__
#define __StgDomain_Utils_MeshBoundaryShape_h__

/* Textual name of this class */
extern const Type MeshBoundaryShape_Type;

/* Class contents. */
#define __MeshBoundaryShape                     \
   __Stg_Shape                                  \
   Mesh* mesh;                                  \
   CartesianGenerator* gen;                     \
   int depth;                                   \
   Bool walls[6];

struct MeshBoundaryShape { __MeshBoundaryShape };


/*
** Constructors */



	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define MESHBOUNDARYSHAPE_DEFARGS \
                STG_SHAPE_DEFARGS

	#define MESHBOUNDARYSHAPE_PASSARGS \
                STG_SHAPE_PASSARGS

MeshBoundaryShape* MeshBoundaryShape_New( Name name );
MeshBoundaryShape* _MeshBoundaryShape_New(  MESHBOUNDARYSHAPE_DEFARGS  );
void _MeshBoundaryShape_Init( MeshBoundaryShape* _self );

void _MeshBoundaryShape_Delete( void* _self );
void _MeshBoundaryShape_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data );
void _MeshBoundaryShape_Build( void* _self, void* data );
void _MeshBoundaryShape_Initialise( void* _self, void* data );

Bool _MeshBoundaryShape_IsCoordInside( void* _self, const Coord coord ) ;
double _MeshBoundaryShape_CalculateVolume( void* _self );
void _MeshBoundaryShape_DistanceFromCenterAxis( void* _self, const Coord coord, double* disVec );

#endif

