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
** $Id: Mesh_HexType.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>

#include "types.h"
#include "shortcuts.h"
#include "Decomp.h"
#include "Sync.h"
#include "MeshTopology.h"
#include "IGraph.h"
#include "Mesh_ElementType.h"
#include "MeshClass.h"
#include "Mesh_HexType.h"


/* Textual name of this class */
const Type Mesh_HexType_Type = "Mesh_HexType";

const int max_vertices(27);

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

Mesh_HexType* Mesh_HexType_New() {
	/* Variables set in this function */
	SizeT                                                    _sizeOfSelf = sizeof(Mesh_HexType);
	Type                                                            type = Mesh_HexType_Type;
	Stg_Class_DeleteFunction*                                    _delete = _Mesh_HexType_Delete;
	Stg_Class_PrintFunction*                                      _print = _Mesh_HexType_Print;
	Stg_Class_CopyFunction*                                        _copy = NULL;
	Mesh_ElementType_UpdateFunc*                              updateFunc = Mesh_HexType_Update;
	Mesh_ElementType_ElementHasPointFunc*            elementHasPointFunc = Mesh_HexType_ElementHasPoint;
	Mesh_ElementType_GetMinimumSeparationFunc*  getMinimumSeparationFunc = Mesh_HexType_GetMinimumSeparation;
	Mesh_ElementType_GetCentroidFunc*                    getCentroidFunc = _Mesh_ElementType_GetCentroid;

	return _Mesh_HexType_New(  MESH_HEXTYPE_PASSARGS  );
}

Mesh_HexType* _Mesh_HexType_New(  MESH_HEXTYPE_DEFARGS  ) {
	Mesh_HexType* self;

	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(Mesh_HexType) );
	self = (Mesh_HexType*)_Mesh_ElementType_New(  MESH_ELEMENTTYPE_PASSARGS  );

	/* Virtual info */

	/* Mesh_HexType info */
	_Mesh_HexType_Init( self );

	return self;
}

void _Mesh_HexType_Init( Mesh_HexType* self ) {
	assert( self && Stg_CheckType( self, Mesh_HexType ) );

        self->num_simplexes[0]=2;
        self->num_simplexes[1]=10;
	self->vertMap = AllocArray( unsigned, max_vertices );
	self->inc = AllocArray( unsigned, max_vertices );
	Mesh_HexType_SetVertexMap( self, NULL );

	self->elementHasPoint = NULL;

	self->triInds = AllocArray2D( unsigned, 2, 3 );
	self->triInds[0][0] = 0; self->triInds[0][1] = 1; self->triInds[0][2] = 2;
	self->triInds[1][0] = 1; self->triInds[1][1] = 3; self->triInds[1][2] = 2;

	self->tetInds = AllocArray2D( unsigned, 10, 4 );
	self->tetInds[0][0] = 0; self->tetInds[0][1] = 1; self->tetInds[0][2] = 2; self->tetInds[0][3] = 4;
	self->tetInds[1][0] = 1; self->tetInds[1][1] = 2; self->tetInds[1][2] = 3; self->tetInds[1][3] = 7;
	self->tetInds[2][0] = 1; self->tetInds[2][1] = 4; self->tetInds[2][2] = 5; self->tetInds[2][3] = 7;
	self->tetInds[3][0] = 2; self->tetInds[3][1] = 4; self->tetInds[3][2] = 6; self->tetInds[3][3] = 7;
	self->tetInds[4][0] = 1; self->tetInds[4][1] = 2; self->tetInds[4][2] = 4; self->tetInds[4][3] = 7;
	self->tetInds[5][0] = 0; self->tetInds[5][1] = 1; self->tetInds[5][2] = 3; self->tetInds[5][3] = 5;
	self->tetInds[6][0] = 0; self->tetInds[6][1] = 4; self->tetInds[6][2] = 5; self->tetInds[6][3] = 6;
	self->tetInds[7][0] = 0; self->tetInds[7][1] = 2; self->tetInds[7][2] = 3; self->tetInds[7][3] = 6;
	self->tetInds[8][0] = 3; self->tetInds[8][1] = 5; self->tetInds[8][2] = 6; self->tetInds[8][3] = 7;
	self->tetInds[9][0] = 0; self->tetInds[9][1] = 3; self->tetInds[9][2] = 5; self->tetInds[9][3] = 6;

	self->incArray = IArray_New();
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Mesh_HexType_Delete( void* elementType ) {
	Mesh_HexType*	self = (Mesh_HexType*)elementType;

	FreeArray( self->vertMap );
	FreeArray( self->inc );
	FreeArray( self->triInds );
	FreeArray( self->tetInds );
	NewClass_Delete( self->incArray );

	/* Delete the parent. */
	_Mesh_ElementType_Delete( self );
}

void _Mesh_HexType_Print( void* elementType, Stream* stream ) {
	Mesh_HexType*	self = (Mesh_HexType*)elementType;
	Stream*			elementTypeStream;

	elementTypeStream = Journal_Register( InfoStream_Type, (Name)"Mesh_HexTypeStream"  );

	/* Print parent */
	Journal_Printf( stream, "Mesh_HexType (ptr): (%p)\n", self );
	_Mesh_ElementType_Print( self, stream );
}

void Mesh_HexType_Update( void* hexType ) {
	Mesh_HexType*	self = (Mesh_HexType*)hexType;
	unsigned	nDims;
	unsigned	d_i;

	assert( self && Stg_CheckType( self, Mesh_HexType ) );

	nDims = Mesh_GetDimSize( self->mesh );
	for( d_i = 0; d_i < nDims; d_i++ ) {
		if( Class_IsSuper( self->mesh->topo, IGraph ) && (!Mesh_GetGlobalSize( self->mesh, (MeshTopology_Dim)d_i ) || !Mesh_HasIncidence( self->mesh, (MeshTopology_Dim)nDims, (MeshTopology_Dim)d_i )) ) {
			break;
		}
	}

	if( Mesh_GetDimSize( self->mesh ) == 3 ) {
		if( d_i == nDims ) {
			self->elementHasPoint = (Mesh_ElementType_ElementHasPointFunc*)Mesh_HexType_ElementHasPoint3DWithIncidence;
		}
		else {
			self->elementHasPoint = (Mesh_ElementType_ElementHasPointFunc*)Mesh_HexType_ElementHasPoint3DGeneral;
		}
	}
	else if( Mesh_GetDimSize( self->mesh ) == 2 ) {
		if( d_i == nDims ) {
			self->elementHasPoint = (Mesh_ElementType_ElementHasPointFunc*)Mesh_HexType_ElementHasPoint2DWithIncidence;
		}
		else {
			self->elementHasPoint = (Mesh_ElementType_ElementHasPointFunc*)Mesh_HexType_ElementHasPoint2DGeneral;
		}
	}
	else {
		if( d_i == nDims ) {
			self->elementHasPoint = (Mesh_ElementType_ElementHasPointFunc*)Mesh_HexType_ElementHasPoint1DWithIncidence;
		}
		else {
			self->elementHasPoint = (Mesh_ElementType_ElementHasPointFunc*)Mesh_HexType_ElementHasPoint1DGeneral;
		}
	}
}

Bool Mesh_HexType_ElementHasPoint( void* hexType, unsigned elInd, double* point, MeshTopology_Dim* dim, unsigned* ind ) {
	Mesh_HexType*	self = (Mesh_HexType*)hexType;

	assert( self && Stg_CheckType( self, Mesh_HexType ) );
	assert( Mesh_GetDimSize( self->mesh ) <= 3 );
	assert( self->elementHasPoint );

	return self->elementHasPoint( self, elInd, point, dim, ind );
}

double Mesh_HexType_GetMinimumSeparation( void* hexType, unsigned elInd, double* perDim ) {
	Mesh_HexType*	self = (Mesh_HexType*)hexType;
	unsigned*		map = NULL;
	double			curSep = 0.0;
	double*			dimSep = NULL;
	unsigned		e_i, nInc = 0;
	int			*inc = NULL;

	assert( self );
	assert( elInd < Mesh_GetDomainSize( self->mesh, Mesh_GetDimSize( self->mesh ) ) );
	assert( Mesh_GetDimSize( self->mesh ) <= 3 );

	/*
	** We know we're a hexahedral element but we may not be regular.  This algorithm (originally from
	** FeVariable.c) doesn't calculate the exact separation but provides an answer that's pretty
	** close.
	*/

	dimSep = AllocArray( double, Mesh_GetDimSize( self->mesh ) );

	for( e_i = 0; e_i > (unsigned)Mesh_GetDimSize( self->mesh ); e_i++ )
		dimSep[e_i] = 0.0;

	Mesh_GetIncidence( self->mesh, Mesh_GetDimSize( self->mesh ), elInd, MT_VERTEX, self->incArray );
	nInc = IArray_GetSize( self->incArray );
	inc = IArray_GetPtr( self->incArray );
	map = self->vertMap;

	curSep = Mesh_GetVertex( self->mesh, inc[map[1]] )[0] - Mesh_GetVertex( self->mesh, inc[map[0]] )[0];
	dimSep[0] = curSep;	

	if( Mesh_GetDimSize( self->mesh ) >= 2 ) {
		curSep = Mesh_GetVertex( self->mesh, inc[map[3]] )[0] - Mesh_GetVertex( self->mesh, inc[map[2]] )[0];
		if( curSep < dimSep[0] )
			dimSep[0] = curSep;
	}

	if( Mesh_GetDimSize( self->mesh ) == 3 ) {
		curSep = Mesh_GetVertex( self->mesh, inc[map[5]] )[0] - Mesh_GetVertex( self->mesh, inc[map[4]] )[0];
		if( curSep < dimSep[0] )
			dimSep[0] = curSep;

		curSep = Mesh_GetVertex( self->mesh, inc[map[7]] )[0] - Mesh_GetVertex( self->mesh, inc[map[6]] )[0];
		if( curSep < dimSep[0] )
			dimSep[0] = curSep;
	}

	if( Mesh_GetDimSize( self->mesh ) >= 2 ) {
		dimSep[1] = Mesh_GetVertex( self->mesh, inc[map[2]] )[1] - Mesh_GetVertex( self->mesh, inc[map[0]] )[1];

		curSep = Mesh_GetVertex( self->mesh, inc[map[3]] )[1] - Mesh_GetVertex( self->mesh, inc[map[1]] )[1];
		if( curSep < dimSep[1] )
			dimSep[1] = curSep;
	}

	if( Mesh_GetDimSize( self->mesh ) == 3 ) {
		curSep = Mesh_GetVertex( self->mesh, inc[map[6]] )[1] - Mesh_GetVertex( self->mesh, inc[map[4]] )[1];
		if( curSep < dimSep[1] )
			dimSep[1] = curSep;

		curSep = Mesh_GetVertex( self->mesh, inc[map[7]] )[1] - Mesh_GetVertex( self->mesh, inc[map[5]] )[1];
		if( curSep < dimSep[1] )
			dimSep[1] = curSep;
	}

	if( Mesh_GetDimSize( self->mesh ) == 3 ) {
		curSep = Mesh_GetVertex( self->mesh, inc[map[4]] )[2] - Mesh_GetVertex( self->mesh, inc[map[0]] )[2];
		dimSep[2] = curSep;

		curSep = Mesh_GetVertex( self->mesh, inc[map[5]] )[2] - Mesh_GetVertex( self->mesh, inc[map[1]] )[2];
		if( curSep < dimSep[2] )
			dimSep[2] = curSep;

		curSep = Mesh_GetVertex( self->mesh, inc[map[6]] )[2] - Mesh_GetVertex( self->mesh, inc[map[2]] )[2];
		if ( curSep < dimSep[2] )
			dimSep[2] = curSep;

		curSep = Mesh_GetVertex( self->mesh, inc[map[7]] )[2] - Mesh_GetVertex( self->mesh, inc[map[3]] )[2];
		if ( curSep < dimSep[2] )
			dimSep[2] = curSep;
	}

	curSep = dimSep[0];
	if( Mesh_GetDimSize( self->mesh ) >= 2 ) {
		curSep = dimSep[1] < curSep ? dimSep[1] : curSep;
		if ( Mesh_GetDimSize( self->mesh ) == 3 )
			curSep = dimSep[2] < curSep ? dimSep[2] : curSep;
	}

	if( perDim )
		memcpy( perDim, dimSep, Mesh_GetDimSize( self->mesh ) * sizeof(double) );

	FreeArray( dimSep );

	return curSep;
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void Mesh_HexType_SetVertexMap( void* hexType, unsigned* map ) {
	Mesh_HexType*	self = (Mesh_HexType*)hexType;
	unsigned	v_i;

	assert( self && Stg_CheckType( self, Mesh_HexType ) );

	if( map ) {
		self->mapSize = (Mesh_GetDimSize( self->mesh ) == 3) ? 8 : 4;
		memcpy( self->vertMap, map, self->mapSize * sizeof(unsigned) );
	}
	else {
		self->mapSize = 0;
		for( v_i = 0; v_i < max_vertices; v_i++ )
			self->vertMap[v_i] = v_i;
	}
}

/* Modify triInds and tetInds for Q2 elements */

void Mesh_HexType_SetQ2Inds( void* hexType) {
  Mesh_HexType*	self = (Mesh_HexType*)hexType;

  unsigned index_map[]={0,1,3,4,9,10,12,13};
  unsigned start_index[]={0,1,3,4,9,10,12,13};

  /* Set vertmap so that MinimumSeparation will work */
  for(int i=0; i<max_vertices; i++)
    self->vertMap[i] = index_map[i];

  unsigned triInds[8][3], tetInds[80][4];

  for(int n=0;n<4;++n)
    for(int i=0;i<2;++i)
      for(int j=0;j<3;++j)
        {
          triInds[2*n+i][j]=index_map[self->triInds[i][j]]+start_index[n];
        }

  self->triInds = ReallocArray2D( self->triInds, unsigned, 8, 3 );
  for(int i=0;i<8;++i)
    for(int j=0;j<3;++j)
      self->triInds[i][j]=triInds[i][j];

  for(int n=0;n<8;++n)
    for(int i=0;i<10;++i)
      for(int j=0;j<4;++j)
        {
          tetInds[10*n+i][j]=index_map[self->tetInds[i][j]]+start_index[n];
        }

  self->tetInds = ReallocArray2D( self->tetInds, unsigned, 80, 4 );
  for(int i=0;i<80;++i)
    for(int j=0;j<4;++j)
      self->tetInds[i][j]=tetInds[i][j];

  self->num_simplexes[0]=8;
  self->num_simplexes[1]=80;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

Bool Mesh_HexType_ElementHasPoint3DGeneral( Mesh_HexType* self, unsigned elInd, double* point, 
					    MeshTopology_Dim* dim, unsigned* ind )
{
	Mesh*		mesh;
	unsigned	nInc;
	double		bc[4];
	unsigned	inside;
	const int	*inc;

	assert( self && Stg_CheckType( self, Mesh_HexType ) );
	assert( Mesh_GetDimSize( self->mesh ) == 3 );
	assert( elInd < Mesh_GetDomainSize( self->mesh, Mesh_GetDimSize( self->mesh ) ) );
	assert( point );
	assert( dim );
	assert( ind );

	/* Shortcuts. */
	mesh = self->mesh;

	/* Get element to vertex incidence. */
	Mesh_GetIncidence( mesh, Mesh_GetDimSize( mesh ), elInd, MT_VERTEX, self->incArray );
	nInc = IArray_GetSize( self->incArray );
	inc = IArray_GetPtr( self->incArray );

	/* Search for tetrahedra. */
	if( self->mapSize ) {
		unsigned	v_i;

		for( v_i = 0; v_i < self->mapSize; v_i++ )
			self->inc[v_i] = inc[self->vertMap[v_i]];
		if( Simplex_Search3D( mesh->verts, self->inc,
                                      self->num_simplexes[1], self->tetInds,
                                      point, bc, &inside ) ) {
			*dim = MT_VOLUME;
			*ind = elInd;
			return True;
		}
	}
	else {
		if( Simplex_Search3D( mesh->verts, (unsigned*)inc,
                                      self->num_simplexes[1], self->tetInds,
                                      point, bc, &inside ) ) {
			*dim = MT_VOLUME;
			*ind = elInd;
			return True;
		}
	}

	return False;
}

Bool Mesh_HexType_ElementHasPoint3DWithIncidence( Mesh_HexType* self, unsigned elInd, double* point, 
						  MeshTopology_Dim* dim, unsigned* ind )
{
  Mesh*		mesh;
  unsigned	nInc;
  Bool		fnd;
  double		bc[4];
  IGraph*		topo;
  unsigned	inside;
  const int*	inc;

  assert( self && Stg_CheckType( self, Mesh_HexType ) );
  assert( Mesh_GetDimSize( self->mesh ) == 3 );
  assert( elInd < Mesh_GetDomainSize( self->mesh, Mesh_GetDimSize( self->mesh ) ) );
  assert( point );
  assert( dim );
  assert( ind );

  /* Shortcuts. */
  mesh = self->mesh;
  topo = (IGraph*)mesh->topo;

  /* Get element to vertex incidence. */
  Mesh_GetIncidence( mesh, Mesh_GetDimSize( mesh ), elInd, MT_VERTEX, self->incArray );
  nInc = IArray_GetSize( self->incArray );
  inc = IArray_GetPtr( self->incArray );

  /* Search for tetrahedra. */
  if( self->mapSize ) {
    unsigned	v_i;

    for( v_i = 0; v_i < self->mapSize; v_i++ )
      self->inc[v_i] = inc[self->vertMap[v_i]];
    fnd = Simplex_Search3D( mesh->verts, self->inc,
                            self->num_simplexes[1], self->tetInds,
                            point, bc, &inside );
  }
  else
    fnd = Simplex_Search3D( mesh->verts, (unsigned*)inc,
                            self->num_simplexes[1], self->tetInds,
                            point, bc, &inside );
  if( fnd ) {
    unsigned*	inds = self->tetInds[inside];

    /* Check boundary ownership. */
    if( bc[0] == 0.0 || bc[0] == -0.0 ) {
      if( bc[1] == 0.0 || bc[1] == -0.0 ) {
        if( bc[2] == 0.0 || bc[2] == -0.0 ) {
          *dim = MT_VERTEX;
          *ind = topo->incEls[MT_VOLUME][MT_VERTEX][elInd][inds[3]];
        }
        else if( bc[3] == 0.0 || bc[3] == -0.0 ) {
          *dim = MT_VERTEX;
          *ind = topo->incEls[MT_VOLUME][MT_VERTEX][elInd][inds[2]];
        }
        else {
          switch(self->num_simplexes[1])
            {
            case 10:
              if( inside == 0 ) {
                *dim = MT_FACE;
                *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][4];
              }
              else if( inside == 1 ) {
                *dim = MT_EDGE;
                *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][11];
              }
              else if( inside == 2 ) {
                *dim = MT_EDGE;
                *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][7];
              }
              else if( inside == 3 ) {
                *dim = MT_EDGE;
                *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][5];
              }
              else if( inside == 4 ) {
                *dim = MT_FACE;
                *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][1];
              }
              else if( inside == 5 ) {
                *dim = MT_FACE;
                *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][5];
              }
              else if( inside == 6 ) {
                *dim = MT_FACE;
                *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][1];
              }
              else if( inside == 7 ) {
                *dim = MT_FACE;
                *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][3];
              }
              else if( inside == 8 ) {
                *dim = MT_EDGE;
                *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][5];
              }
              else {
                *dim = MT_FACE;
                *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][1];
              }
              break;
            case 80:
              abort();
              break;
            default:
              abort;
            }
        }
      }
      else if( bc[2] == 0.0 || bc[2] == -0.0 ) {
        if( bc[3] == 0.0 || bc[3] == -0.0 ) {
          *dim = MT_VERTEX;
          *ind = topo->incEls[MT_VOLUME][MT_VERTEX][elInd][inds[1]];
        }
        else {
          switch(self->num_simplexes[1])
            {
            case 10:
              if( inside == 0 ) {
                *dim = MT_FACE;
                *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][2];
              }
              else if( inside == 1 ) {
                *dim = MT_FACE;
                *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][3];
              }
              else if( inside == 2 ) {
                *dim = MT_FACE;
                *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][1];
              }
              else if( inside == 3 ) {
                *dim = MT_FACE;
                *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][1];
              }
              else if( inside == 4 ) {
                *dim = MT_FACE;
                *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][3];
              }
              else if( inside == 5 ) {
                *dim = MT_EDGE;
                *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][9];
              }
              else if( inside == 6 ) {
                *dim = MT_EDGE;
                *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][6];
              }
              else if( inside == 7 ) {
                *dim = MT_EDGE;
                *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][10];
              }
              else if( inside == 8 ) {
                *dim = MT_EDGE;
                *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][7];
              }
              else {
                *dim = MT_FACE;
                *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][3];
              }
              break;
            case 80:
              abort();
              break;
            default:
              abort;
            }
        }
      }
      else if( bc[3] == 0.0 || bc[3] == -0.0 ) {
        switch(self->num_simplexes[1])
          {
          case 10:
            if( inside == 0 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][0];
            }
            else if( inside == 1 ) {
              *dim = MT_EDGE;
              *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][1];
            }
            else if( inside == 2 ) {
              *dim = MT_EDGE;
              *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][4];
            }
            else if( inside == 3 ) {
              *dim = MT_EDGE;
              *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][6];
            }
            else if( inside == 4 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][4];
            }
            else if( inside == 5 ) {
              *dim = MT_EDGE;
              *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][3];
            }
            else if( inside == 6 ) {
              *dim = MT_EDGE;
              *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][4];
            }
            else if( inside == 7 ) {
              *dim = MT_EDGE;
              *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][1];
            }
            else if( inside == 8 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][1];
            }
            else {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][5];
            }
            break;
          case 80:
            abort();
            break;
          default:
            abort;
          }
      }
      else {
        switch(self->num_simplexes[1])
          {
          case 10:
            if( inside == 0 ) {
              *dim = MT_VOLUME;
              *ind = elInd;
            }
            else if( inside == 1 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][3];
            }
            else if( inside == 2 ) {
              *dim = 	MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][1];
            }
            else if( inside == 3 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][1];
            }
            else if( inside == 4 ) {
              *dim = MT_VOLUME;
              *ind = elInd;
            }
            else if( inside == 5 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][5];
            }
            else if( inside == 6 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][1];
            }
            else if( inside == 7 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][3];
            }
            else if( inside == 8 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][1];
            }
            else {
              *dim = MT_VOLUME;
              *ind = elInd;
            }
            break;
          case 80:
            abort();
            break;
          default:
            abort;
          }
      }
    }
    else if( bc[1] == 0.0 || bc[1] == -0.0 ) {
      if( bc[2] == 0.0 || bc[2] == -0.0 ) {
        if( bc[3] == 0.0 || bc[3] == -0.0 ) {
          *dim = MT_VERTEX;
          *ind = topo->incEls[MT_VOLUME][MT_VERTEX][elInd][inds[0]];
        }
        else {
          switch(self->num_simplexes[1])
            {
            case 10:
              if( inside == 0 ) {
                *dim = MT_EDGE;
                *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][8];
              }
              else if( inside == 1 ) {
                *dim = MT_FACE;
                *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][5];
              }
              else if( inside == 2 ) {
                *dim = MT_FACE;
                *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][5];
              }
              else if( inside == 3 ) {
                *dim = MT_FACE;
                *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][3];
              }
              else if( inside == 4 ) {
                *dim = MT_FACE;
                *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][5];
              }
              else if( inside == 5 ) {
                *dim = MT_FACE;
                *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][2];
              }
              else if( inside == 6 ) {
                *dim = MT_FACE;
                *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][4];
              }
              else if( inside == 7 ) {
                *dim = MT_FACE;
                *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][4];
              }
              else if( inside == 8 ) {
                *dim = MT_EDGE;
                *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][11];
              }
              else {
                *dim = MT_FACE;
                *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][4];
              }
              break;
            case 80:
              abort();
              break;
            default:
              abort;
            }
        }
      }
      else if( bc[3] == 0.0 || bc[3] == -0.0 ) {
        switch(self->num_simplexes[1])
          {
          case 10:
            if( inside == 0 ) {
              *dim = MT_EDGE;
              *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][2];
            }
            else if( inside == 1 ) {
              *dim = MT_EDGE;
              *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][3];
            }
            else if( inside == 2 ) {
              *dim = MT_EDGE;
              *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][9];
            }
            else if( inside == 3 ) {
              *dim = MT_EDGE;
              *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][10];
            }
            else if( inside == 4 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][2];
            }
            else if( inside == 5 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][0];
            }
            else if( inside == 6 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][2];
            }
            else if( inside == 7 ) {
              *dim = 	MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][0];
            }
            else if( inside == 8 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][3];
            }
            else {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][2];
            }
            break;
          case 80:
            abort();
            break;
          default:
            abort;
          }
      }
      else {
        switch(self->num_simplexes[1])
          {
          case 10:
            if( inside == 0 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][4];
            }
            else if( inside == 1 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][5];
            }
            else if( inside == 2 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][5];
            }
            else if( inside == 3 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][3];
            }
            else if( inside == 4 ) {
              *dim = MT_VOLUME;
              *ind = elInd;
            }
            else if( inside == 5 ) {
              *dim = MT_VOLUME;
              *ind = elInd;
            }
            else if( inside == 6 ) {
              *dim = MT_VOLUME;
              *ind = elInd;
            }
            else if( inside == 7 ) {
              *dim = MT_VOLUME;
              *ind = elInd;
            }
            else if( inside == 8 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][3];
            }
            else {
              *dim = MT_VOLUME;
              *ind = elInd;
            }
            break;
          case 80:
            abort();
            break;
          default:
            abort;
          }
      }
    }
    else if( bc[2] == 0.0 || bc[2] == -0.0 ) {
      if( bc[3] == 0.0 || bc[3] == -0.0 ) {
        switch(self->num_simplexes[1])
          {
          case 10:
            if( inside == 0 ) {
              *dim = MT_EDGE;
              *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][0];
            }
            else if( inside == 1 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][0];
            }
            else if( inside == 2 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][2];
            }
            else if( inside == 3 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][4];
            }
            else if( inside == 4 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][0];
            }
            else if( inside == 5 ) {
              *dim = MT_EDGE;
              *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][0];
            }
            else if( inside == 6 ) {
              *dim = MT_EDGE;
              *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][8];
            }
            else if( inside == 7 ) {
              *dim = MT_EDGE;
              *ind = topo->incEls[MT_VOLUME][MT_EDGE][elInd][2];
            }
            else if( inside == 8 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][5];
            }
            else {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][0];
            }
            break;
          case 80:
            abort();
            break;
          default:
            abort;
          }
      }
      else {
        switch(self->num_simplexes[1])
          {
          case 10:
            if( inside == 0 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][2];
            }
            else if( inside == 1 ) {
              *dim = MT_VOLUME;
              *ind = elInd;
            }
            else if( inside == 2 ) {
              *dim = MT_VOLUME;
              *ind = elInd;
            }
            else if( inside == 3 ) {
              *dim = MT_VOLUME;
              *ind = elInd;
            }
            else if( inside == 4 ) {
              *dim = MT_VOLUME;
              *ind = elInd;
            }
            else if( inside == 5 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][2];
            }
            else if( inside == 6 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][4];
            }
            else if( inside == 7 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][4];
            }
            else if( inside == 8 ) {
              *dim = MT_FACE;
              *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][5];
            }
            else {
              *dim = MT_VOLUME;
              *ind = elInd;
            }
            break;
          case 80:
            abort();
            break;
          default:
            abort;
          }
      }
    }
    else if( bc[3] == 0.0 || bc[3] == -0.0 ) {
      switch(self->num_simplexes[1])
        {
        case 10:
          if( inside == 0 ) {
            *dim = MT_FACE;
            *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][0];
          }
          else if( inside == 1 ) {
            *dim = MT_FACE;
            *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][0];
          }
          else if( inside == 2 ) {
            *dim = MT_FACE;
            *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][2];
          }
          else if( inside == 3 ) {
            *dim = MT_FACE;
            *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][4];
          }
          else if( inside == 4 ) {
            *dim = MT_VOLUME;
            *ind = elInd;
          }
          else if( inside == 5 ) {
            *dim = MT_FACE;
            *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][0];
          }
          else if( inside == 6 ) {
            *dim = MT_FACE;
            *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][2];
          }
          else if( inside == 7 ) {
            *dim = MT_FACE;
            *ind = topo->incEls[MT_VOLUME][MT_FACE][elInd][0];
          }
          else if( inside == 8 ) {
            *dim = MT_VOLUME;
            *ind = elInd;
          }
          else {
            *dim = MT_VOLUME;
            *ind = elInd;
          }
          break;
        case 80:
          abort();
          break;
        default:
          abort;
        }
    }
    else {
      *dim = MT_VOLUME;
      *ind = elInd;
    }
    return True;
  }
  return False;
}

Bool Mesh_HexType_ElementHasPoint2DGeneral( Mesh_HexType* self, unsigned elInd, double* point, 
					    MeshTopology_Dim* dim, unsigned* ind )
{
	Mesh*		mesh;
	unsigned	nInc;
	Bool		fnd;
	double		bc[3];
	unsigned	inside;
	const int*	inc;

	assert( self && Stg_CheckType( self, Mesh_HexType ) );
	assert( Mesh_GetDimSize( self->mesh ) == 2 );
	assert( elInd < Mesh_GetDomainSize( self->mesh, Mesh_GetDimSize( self->mesh ) ) );
	assert( point );
	assert( dim );
	assert( ind );

	/* Shortcuts. */
	mesh = self->mesh;

	/* Get element to vertex incidence. */
	Mesh_GetIncidence( mesh, Mesh_GetDimSize( mesh ), elInd, MT_VERTEX, self->incArray );
	nInc = IArray_GetSize( self->incArray );
	inc = IArray_GetPtr( self->incArray );

	/* Search for triangle. */
	if( self->mapSize ) {
		unsigned	v_i;

		for( v_i = 0; v_i < self->mapSize; v_i++ )
			self->inc[v_i] = inc[self->vertMap[v_i]];
		fnd = Simplex_Search2D( mesh->verts, self->inc,
                                        self->num_simplexes[0], self->triInds,
                                        point, bc, &inside );
	}
	else
		fnd = Simplex_Search2D( mesh->verts, (unsigned*)inc,
                                        self->num_simplexes[0], self->triInds,
                                        point, bc, &inside );
	if( fnd ) {
		*dim = MT_FACE;
		*ind = elInd;
		return True;
	}

	return False;
}

Bool Mesh_HexType_ElementHasPoint2DWithIncidence( Mesh_HexType* self, unsigned elInd, double* point, 
						  MeshTopology_Dim* dim, unsigned* ind )
{
  Mesh*		mesh;
  unsigned	nInc;
  Bool		fnd;
  double		bc[3];
  IGraph*		topo;
  unsigned	inside;
  const int*	inc;

  assert( self && Stg_CheckType( self, Mesh_HexType ) );
  assert( Mesh_GetDimSize( self->mesh ) == 2 );
  assert( elInd < Mesh_GetDomainSize( self->mesh, Mesh_GetDimSize( self->mesh ) ) );
  assert( point );
  assert( dim );
  assert( ind );

  /* Shortcuts. */
  mesh = self->mesh;
  topo = (IGraph*)mesh->topo;

  /* Get element to vertex incidence. */
  Mesh_GetIncidence( mesh, Mesh_GetDimSize( mesh ), elInd, MT_VERTEX, self->incArray );
  nInc = IArray_GetSize( self->incArray );
  inc = IArray_GetPtr( self->incArray );

  /* Search for triangle. */
  if( self->mapSize ) {
    unsigned	v_i;

    for( v_i = 0; v_i < self->mapSize; v_i++ )
      self->inc[v_i] = inc[self->vertMap[v_i]];
    fnd = Simplex_Search2D( mesh->verts, self->inc,
                            self->num_simplexes[0],
                            self->triInds, point, bc,
                            &inside );
  }
  else
    fnd = Simplex_Search2D( mesh->verts, (unsigned*)inc,
                            self->num_simplexes[0],
                            self->triInds, point, bc, &inside );
  if( fnd ) {
    unsigned	*inds = self->triInds[inside];

    /* Check boundary ownership. */
    if( bc[0] == 0.0 || bc[0] == -0.0 ) {
      if( bc[1] == 0.0 || bc[1] == -0.0 ) {
        *dim = MT_VERTEX;
        *ind = topo->incEls[MT_FACE][MT_VERTEX][elInd][inds[2]];
      }
      else if( bc[2] == 0.0 || bc[2] == -0.0 ) {
        *dim = MT_VERTEX;
        *ind = topo->incEls[MT_FACE][MT_VERTEX][elInd][inds[1]];
      }
      else {
        switch(self->num_simplexes[0])
          {
          case 2:
            switch(inside)
              {
              case 0:
                *dim = MT_FACE;
                *ind = elInd;
                break;
              case 1:
                *dim = MT_EDGE;
                *ind = topo->incEls[MT_FACE][MT_EDGE][elInd][1];
                break;
              }
            break;
          case 8:
            switch(inside)
              {
              case 0:
              case 1:
              case 2:
              case 3:
              case 4:
              case 6:
                *dim = MT_FACE;
                *ind = elInd;
                break;
              case 5:
              case 7:
                *dim = MT_EDGE;
                *ind = topo->incEls[MT_FACE][MT_EDGE][elInd][1];
                break;
              }
            break;
          default:
            abort();
          }
      }
    }
    else if( bc[1] == 0.0 || bc[1] == -0.0 ) {
      if( bc[2] == 0.0 || bc[2] == -0.0 ) {
        *dim = MT_VERTEX;
        *ind = topo->incEls[MT_FACE][MT_VERTEX][elInd][inds[0]];
      }
      else {
        switch(self->num_simplexes[0])
          {
          case 2:
            switch(inside)
              {
              case 0:
                *dim = MT_EDGE;
                *ind = topo->incEls[MT_FACE][MT_EDGE][elInd][2];
                break;
              case 1:
                *dim = MT_FACE;
                *ind = elInd;
                break;
              }
            break;
          case 8:
            switch(inside)
              {
              case 0:
              case 4:
                *dim = MT_EDGE;
                *ind = topo->incEls[MT_FACE][MT_EDGE][elInd][2];
                break;
              case 1:
              case 2:
              case 3:
              case 5:
              case 6:
              case 7:
                *dim = MT_FACE;
                *ind = elInd;
                break;
              }
          default:
            abort;
          }
      }
    }
    else if( bc[2] == 0.0 || bc[2] == -0.0 ) {
      switch(self->num_simplexes[0])
        {
        case 2:
          switch(inside)
            {
            case 0:
              *dim = MT_EDGE;
              *ind = topo->incEls[MT_FACE][MT_EDGE][elInd][0];
              break;
            case 1:
              *dim = MT_EDGE;
              *ind = topo->incEls[MT_FACE][MT_EDGE][elInd][3];
              break;
            }
          break;
        case 8:
          switch(inside)
            {
            case 0:
            case 2:
              *dim = MT_EDGE;
              *ind = topo->incEls[MT_FACE][MT_EDGE][elInd][0];
              break;
            case 1:
            case 4:
            case 5:
            case 6:
              *dim = MT_FACE;
              *ind = elInd;
              break;
            case 3:
            case 7:
              *dim = MT_EDGE;
              *ind = topo->incEls[MT_FACE][MT_EDGE][elInd][3];
              break;
            }
          break;
        default:
          abort;
        }
    }
    else {
      *dim = MT_FACE;
      *ind = elInd;
    }
    return True;
  }
  return False;
}

Bool Mesh_HexType_ElementHasPoint1DGeneral( Mesh_HexType* self, unsigned elInd, double* point, 
					    MeshTopology_Dim* dim, unsigned* ind )
{
	Mesh*		mesh;
	unsigned	nInc;
	const int*	inc;

	assert( self && Stg_CheckType( self, Mesh_HexType ) );
	assert( Mesh_GetDimSize( self->mesh ) == 1 );
	assert( elInd < Mesh_GetDomainSize( self->mesh, Mesh_GetDimSize( self->mesh ) ) );
	assert( point );
	assert( dim );
	assert( ind );

	mesh = self->mesh;
	Mesh_GetIncidence( mesh, MT_EDGE, elInd, MT_VERTEX, self->incArray );
	nInc = IArray_GetSize( self->incArray );
	inc = IArray_GetPtr( self->incArray );

	if( point[0] > *Mesh_GetVertex( mesh, inc[self->vertMap[0]] ) && 
	    point[0] < *Mesh_GetVertex( mesh, inc[self->vertMap[1]] ) )
	{
		*dim = MT_EDGE;
		*ind = elInd;
		return True;
	}

	return False;
}

Bool Mesh_HexType_ElementHasPoint1DWithIncidence( Mesh_HexType* self, unsigned elInd, double* point, 
						  MeshTopology_Dim* dim, unsigned* ind )
{
	Mesh*		mesh;
	unsigned	nInc;
	const int*	inc;

	assert( self && Stg_CheckType( self, Mesh_HexType ) );
	assert( Mesh_GetDimSize( self->mesh ) == 1 );
	assert( elInd < Mesh_GetDomainSize( self->mesh, Mesh_GetDimSize( self->mesh ) ) );
	assert( point );
	assert( dim );
	assert( ind );

	mesh = self->mesh;
	Mesh_GetIncidence( mesh, MT_EDGE, elInd, MT_VERTEX, self->incArray );
	nInc = IArray_GetSize( self->incArray );
	inc = IArray_GetPtr( self->incArray );
	assert( nInc == 2 );

	if( point[0] > *Mesh_GetVertex( mesh, inc[self->vertMap[0]] ) && 
	    point[0] < *Mesh_GetVertex( mesh, inc[self->vertMap[1]] ) )
	{
		*dim = MT_EDGE;
		*ind = elInd;
		return True;
	}
	else if( point[0] == *Mesh_GetVertex( mesh, inc[self->vertMap[0]] ) ) {
		*dim = MT_VERTEX;
		*ind = inc[0];
		return True;
	}
	else if( point[0] == *Mesh_GetVertex( mesh, inc[self->vertMap[1]] ) ) {
		*dim = MT_VERTEX;
		*ind = inc[1];
		return True;
	}

	return False;
}


