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
** $Id: MeshBoundaryShape.c 2192 2004-10-15 02:45:38Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <string.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>

#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>

#include "types.h"
#include "MeshBoundaryShape.h"


/* Textual name of this class */
const Type MeshBoundaryShape_Type = "MeshBoundaryShape";


/*
** Constructors */

MeshBoundaryShape* MeshBoundaryShape_New( Name name ) {
	/* Variables set in this function */
	SizeT                                                  _sizeOfSelf = sizeof(MeshBoundaryShape);
	Type                                                          type = MeshBoundaryShape_Type;
	Stg_Class_DeleteFunction*                                  _delete = _MeshBoundaryShape_Delete;
	Stg_Class_PrintFunction*                                    _print = _Stg_Shape_Print;
	Stg_Class_CopyFunction*                                      _copy = _Stg_Shape_Copy;
	Stg_Component_DefaultConstructorFunction*      _defaultConstructor = (void*(*)(Name))MeshBoundaryShape_New;
	Stg_Component_ConstructFunction*                        _construct = _MeshBoundaryShape_AssignFromXML;
	Stg_Component_BuildFunction*                                _build = _MeshBoundaryShape_Build;
	Stg_Component_InitialiseFunction*                      _initialise = _MeshBoundaryShape_Initialise;
	Stg_Component_ExecuteFunction*                            _execute = _Stg_Shape_Execute;
	Stg_Component_DestroyFunction*                            _destroy = _Stg_Shape_Destroy;
	Stg_Shape_IsCoordInsideFunction*                    _isCoordInside = _MeshBoundaryShape_IsCoordInside;
	Stg_Shape_CalculateVolumeFunction*                _calculateVolume = _MeshBoundaryShape_CalculateVolume;
	Stg_Shape_DistanceFromCenterAxisFunction*  _distanceFromCenterAxis = _MeshBoundaryShape_DistanceFromCenterAxis;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

   return _MeshBoundaryShape_New(  MESHBOUNDARYSHAPE_PASSARGS  );
}

MeshBoundaryShape* _MeshBoundaryShape_New(  MESHBOUNDARYSHAPE_DEFARGS  ) {
   MeshBoundaryShape* self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof(MeshBoundaryShape) );
   self = (MeshBoundaryShape*)_Stg_Shape_New(  STG_SHAPE_PASSARGS  );

   _MeshBoundaryShape_Init( self );

   return self;
}

void _MeshBoundaryShape_Init( MeshBoundaryShape* self ) {
   self->mesh = NULL;
   self->depth = 0;
   memset( self->walls, 0, 6 * sizeof(Bool) );
}


/*
** Virtual functions */

void _MeshBoundaryShape_Delete( void* _self ) {
   MeshBoundaryShape* self = (MeshBoundaryShape*)_self;

   /* Delete parent */
   _Stg_Shape_Delete( self );
}

void _MeshBoundaryShape_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data ) {
   MeshBoundaryShape* self = (MeshBoundaryShape*)_self;
   Dictionary_Entry_Value* wallList;
   int ii;

   _Stg_Shape_AssignFromXML( self, cf, data );
   _MeshBoundaryShape_Init( self );

   /* Need a mesh with a cartesian generator. */
   self->mesh = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"mesh", Mesh, True, data  );

   /* Read in the walls to have friction applied. */
   wallList = _Stg_ComponentFactory_GetDictionaryValue( cf, self->name, "walls", NULL );
   if( wallList ) {
      int nWalls, curWall;
      char* name;

      /* List exists, use it. Extract the number of walls specified. */
      nWalls = Dictionary_Entry_Value_GetCount( wallList );

      /* Read in each wall's name. */
      curWall = 0;
      for( ii = 0; ii < nWalls; ii++ ) {

         /* Rip out the name from the list. */
         name = Dictionary_Entry_Value_AsString(
            Dictionary_Entry_Value_GetElement( wallList, ii ) );

         /* Store enabled wall. */
         if( !strcasecmp( name, "left" ) )
            self->walls[0] = True;
         else if( !strcasecmp( name, "right" ) )
            self->walls[1] = True;
         else if( !strcasecmp( name, "bottom" ) )
            self->walls[2] = True;
         else if( !strcasecmp( name, "top" ) )
            self->walls[3] = True;
         else if( !strcasecmp( name, "back" ) )
            self->walls[4] = True;
         else if( !strcasecmp( name, "front" ) )
            self->walls[5] = True;
         else
            abort();
      }
   }

   /* If no wall list, assume all walls. */
   else {
      for( ii = 0; ii < 6; ii++ )
         self->walls[ii] = True;
   }

}


/*
** Virtual functions */

void _MeshBoundaryShape_Build( void* _self, void* data ) {
   MeshBoundaryShape* self = (MeshBoundaryShape*)_self;

   _Stg_Shape_Build( self, data );
   Stg_Component_Build( self->mesh, data, False );
   Stg_Component_Build( self->gen, data, False );
   if( !self->mesh->generator || strcmp( self->mesh->generator->type, CartesianGenerator_Type ) )
      abort();
   self->gen = (CartesianGenerator*)self->mesh->generator;
}

void _MeshBoundaryShape_Destroy( void* _self, void* data ) {
   MeshBoundaryShape* self = (MeshBoundaryShape*)_self;

   Stg_Component_Destroy( self->mesh, data, False );
   Stg_Component_Destroy( self->gen, data, False );
   _Stg_Shape_Destroy( self, data );
}

void _MeshBoundaryShape_Initialise( void* _self, void* data ) {
   MeshBoundaryShape* self = (MeshBoundaryShape*)_self;

   _Stg_Shape_Initialise( self, data );
   Stg_Component_Initialise( self->mesh, data, False );
   Stg_Component_Initialise( self->gen, data, False );
}

Bool _MeshBoundaryShape_IsCoordInside( void* _self, Coord coord ) {
   MeshBoundaryShape* self = (MeshBoundaryShape*)_self;
   Coord newCoord;
#if 0
   Grid* grid;
   int inds[3], element, nDims, wallInd;
#endif
   int ii;

   /* Transform coordinate into canonical reference frame */
   Stg_Shape_TransformCoord( self, coord, newCoord );

   /* Easy, just check if the coord is in the boundary region. */
   for( ii = 0; ii < Mesh_GetDimSize( self->mesh ); ii++ ) {
      if( (self->walls[2 * ii] && coord[ii] < self->gen->crdMin[ii] + self->gen->contactGeom[ii]) ||
          (self->walls[2 * ii + 1] && coord[ii] > self->gen->crdMax[ii] - self->gen->contactGeom[ii]) )
      {
         return True;
      }
   }
   return False;

#if 0
   /* Get the element grid from the mesh. */
   grid = *(Grid**)Mesh_GetExtension( self->mesh, Grid**, "elementGrid" );
   assert( grid );

   /* Find which element the current coordinate is in. */
   if( !Mesh_SearchElements( self->mesh, newCoord, &element ) ) {

      /* Couldn't find the point inside the mesh. */
      return False;
   }

   /* Convert the element into n-dimensional indices. */
   Grid_Lift( grid, element, inds );

   /* Check if we're on any of the specified boundaries. */
   nDims = Mesh_GetDimSize( self->mesh );
   for( ii = 0; ii < nDims; ii++ ) {
      wallInd = 2 * ii;
      if( (self->walls[wallInd] && inds[ii] < self->depth) ||
          (self->walls[wallInd + 1] && inds[ii] > (grid->sizes[ii] - self->depth - 1)) )
      {
         if( coords[
         return True;
      }
   }

   return False;
#endif
}

double _MeshBoundaryShape_CalculateVolume( void* _self ) {
   abort();
   return 0.0;
}

void _MeshBoundaryShape_DistanceFromCenterAxis( void* _self, Coord coord, double* disVec ) {
   abort();
}



