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
const int cube_vertices(8);

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

	self->vertMap = AllocArray( unsigned, cube_vertices );
	self->inc = AllocArray( unsigned, max_vertices );
	Mesh_HexType_SetVertexMap( self );

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
	unsigned		e_i;
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

void Mesh_HexType_SetVertexMap( void* hexType ) {
  Mesh_HexType*	self = (Mesh_HexType*)hexType;
  int v_i;

  assert( self && Stg_CheckType( self, Mesh_HexType ) );

  for( v_i = 0; v_i < cube_vertices; v_i++ )
    self->vertMap[v_i] = v_i;
}

/* Modify triInds and tetInds for Q2 elements */

void Mesh_HexType_SetQ2Inds( void* hexType) {
  Mesh_HexType*	self = (Mesh_HexType*)hexType;

  unsigned index_map[]={0,1,3,4,9,10,12,13};

  /* Set vertmap so that MinimumSeparation will work */
  for(int i=0; i<cube_vertices; i++)
    self->vertMap[i] = index_map[i];
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

namespace {
  bool Face_Search3D(const int *inc, double **verts, 
                     const double* point, int on_face[6]);
}

bool Mesh_HexType_ElementHasPoint3DGeneral( Mesh_HexType* self, unsigned elInd, double* point, 
					    MeshTopology_Dim* dim, unsigned* ind )
{
  Mesh*		mesh;
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
  inc = IArray_GetPtr( self->incArray );

  bool fnd;
  /* If quad elements, use face search instead of simplexes.*/
  if(IArray_GetSize(self->incArray)==27)
    {
      int on_face[6];
      fnd=Face_Search3D(inc,mesh->verts,point,on_face);
    }
  /* Search for tetrahedra. */
  else
    {
      fnd=Simplex_Search3D( mesh->verts, (unsigned*)inc, 10, self->tetInds,
                            point, bc, &inside);
    }
  if(fnd)
    {
      *dim = MT_VOLUME;
      *ind = elInd;
    }
  return fnd;
}

namespace {
  void set_index(IGraph* topo,
                 const MeshTopology_Dim &volume, const unsigned elInd,
                 const MeshTopology_Dim incidence_dim[],
                 const int incidence_index[], const int &inside,
                 MeshTopology_Dim* dim, unsigned* ind)
  {
    *dim=incidence_dim[inside];
    if(*dim==volume)
      {
        *ind=elInd;
      }
    else
      {
        *ind=topo->incEls[volume][*dim][elInd][incidence_index[inside]];
      }
  }
}


bool Mesh_HexType_ElementHasPoint3DWithIncidence
( Mesh_HexType* self, unsigned elInd, double* point, 
  MeshTopology_Dim* dim, unsigned* ind )
{
  Mesh*		mesh;
  bool		fnd;
  double		bc[4];
  IGraph*		topo;
  unsigned	inside;
  const int*	inc;

  assert( self && Stg_CheckType( self, Mesh_HexType ) );
  assert( Mesh_GetDimSize( self->mesh ) == 3 );
  assert( elInd < Mesh_GetDomainSize( self->mesh,
                                      Mesh_GetDimSize( self->mesh ) ) );
  assert( point );
  assert( dim );
  assert( ind );

  /* Shortcuts. */
  mesh = self->mesh;
  topo = (IGraph*)mesh->topo;

  /* Get element to vertex incidence. */
  Mesh_GetIncidence( mesh, Mesh_GetDimSize( mesh ), elInd, MT_VERTEX,
                     self->incArray );
  inc = IArray_GetPtr( self->incArray );

  /* If quad elements, use face search instead of simplexes.*/
  if(IArray_GetSize(self->incArray)==27)
    {
      int on_face[6];
      fnd=Face_Search3D(inc,mesh->verts,point,on_face);

      if(on_face[0]==0)
        {
          if(on_face[2]==0)
            {
              if(on_face[4]==0)
                {
                  *dim=MT_VERTEX;
                  *ind=topo->incEls[MT_VOLUME][*dim][elInd][0];
                }
              else if(on_face[5]==0)
                {
                  *dim=MT_VERTEX;
                  *ind=topo->incEls[MT_VOLUME][*dim][elInd][2];
                }
              else
                {
                  *dim=MT_EDGE;
                  *ind=topo->incEls[MT_VOLUME][*dim][elInd][0];
                }
            }
          else if(on_face[3]==0)
            {
              if(on_face[4]==0)
                {
                  *dim=MT_VERTEX;
                  *ind=topo->incEls[MT_VOLUME][*dim][elInd][6];
                }
              else if(on_face[5]==0)
                {
                  *dim=MT_VERTEX;
                  *ind=topo->incEls[MT_VOLUME][*dim][elInd][8];
                }
              else
                {
                  *dim=MT_EDGE;
                  *ind=topo->incEls[MT_VOLUME][*dim][elInd][1];
                }
            }
          else if(on_face[4]==0)
            {
              *dim=MT_EDGE;
              *ind=topo->incEls[MT_VOLUME][*dim][elInd][2];
              }
          else if(on_face[5]==0)
            {
              *dim=MT_EDGE;
              *ind=topo->incEls[MT_VOLUME][*dim][elInd][3];
            }
          else
            {
              *dim=MT_FACE;
              *ind=topo->incEls[MT_VOLUME][*dim][elInd][0];
            }
        }
      else if(on_face[1]==0)
        {
          if(on_face[2]==0)
            {
              if(on_face[4]==0)
                {
                  *dim=MT_VERTEX;
                  *ind=topo->incEls[MT_VOLUME][*dim][elInd][18];
                }
              else if(on_face[5]==0)
                {
                  *dim=MT_VERTEX;
                  *ind=topo->incEls[MT_VOLUME][*dim][elInd][20];
                }
              else
                {
                  *dim=MT_EDGE;
                  *ind=topo->incEls[MT_VOLUME][*dim][elInd][4];
                }
            }
          else if(on_face[3]==0)
            {
              if(on_face[4]==0)
                {
                  *dim=MT_VERTEX;
                  *ind=topo->incEls[MT_VOLUME][*dim][elInd][24];
                }
              else if(on_face[5]==0)
                {
                  *dim=MT_VERTEX;
                  *ind=topo->incEls[MT_VOLUME][*dim][elInd][26];
                }
              else
                {
                  *dim=MT_EDGE;
                  *ind=topo->incEls[MT_VOLUME][*dim][elInd][5];
                }
            }
          else if(on_face[4]==0)
            {
              *dim=MT_EDGE;
              *ind=topo->incEls[MT_VOLUME][*dim][elInd][6];
              }
          else if(on_face[5]==0)
            {
              *dim=MT_EDGE;
              *ind=topo->incEls[MT_VOLUME][*dim][elInd][7];
            }
          else
            {
              *dim=MT_FACE;
              *ind=topo->incEls[MT_VOLUME][*dim][elInd][1];
            }
        }
      else if(on_face[2]==0)
        {
          if(on_face[4]==0)
            {
              *dim=MT_EDGE;
              *ind=topo->incEls[MT_VOLUME][*dim][elInd][8];
            }
          else if(on_face[5]==0)
            {
              *dim=MT_EDGE;
              *ind=topo->incEls[MT_VOLUME][*dim][elInd][9];
            }
          else
            {
              *dim=MT_FACE;
              *ind=topo->incEls[MT_VOLUME][*dim][elInd][2];
            }
        }
      else if(on_face[3]==0)
        {
          if(on_face[4]==0)
            {
              *dim=MT_EDGE;
              *ind=topo->incEls[MT_VOLUME][*dim][elInd][10];
            }
          else if(on_face[5]==0)
            {
              *dim=MT_EDGE;
              *ind=topo->incEls[MT_VOLUME][*dim][elInd][11];
            }
          else
            {
              *dim=MT_FACE;
              *ind=topo->incEls[MT_VOLUME][*dim][elInd][3];
            }
        }
      else if(on_face[4]==0)
        {
          *dim=MT_FACE;
          *ind=topo->incEls[MT_VOLUME][*dim][elInd][4];
        }
      else if(on_face[5]==0)
        {
          *dim=MT_FACE;
          *ind=topo->incEls[MT_VOLUME][*dim][elInd][5];
        }
      else
        {
          *dim=MT_VOLUME;
          *ind=elInd;
        }
    }
  /* Search for tetrahedra. */
  else
    {
      fnd = Simplex_Search3D( mesh->verts, (unsigned*)inc, 10, self->tetInds,
                              point, bc, &inside );
      if( fnd ) {
        unsigned*	inds = self->tetInds[inside];

        /* Check boundary ownership. */
        if( bc[0] == 0.0 || bc[0] == -0.0 )
          {
            if( bc[1] == 0.0 || bc[1] == -0.0 )
              {
                if( bc[2] == 0.0 || bc[2] == -0.0 )
                  {
                    *dim = MT_VERTEX;
                    *ind = topo->incEls[MT_VOLUME][MT_VERTEX][elInd][inds[3]];
                  }
                else if( bc[3] == 0.0 || bc[3] == -0.0 )
                  {
                    *dim = MT_VERTEX;
                    *ind = topo->incEls[MT_VOLUME][MT_VERTEX][elInd][inds[2]];
                  }
                else
                  {
                    const MeshTopology_Dim incidence_dim[]=
                      { MT_FACE, MT_EDGE, MT_EDGE, MT_EDGE, MT_FACE, MT_FACE,
                        MT_FACE, MT_FACE, MT_EDGE, MT_FACE };
                    const int incidence_index[]={ 4, 11, 7, 5, 1, 5, 1, 3, 5, 1 };
                    set_index(topo,MT_VOLUME,elInd,
                              incidence_dim,incidence_index,inside,dim,ind);
                  }
              }
            else if( bc[2] == 0.0 || bc[2] == -0.0 )
              {
                if( bc[3] == 0.0 || bc[3] == -0.0 )
                  {
                    *dim = MT_VERTEX;
                    *ind = topo->incEls[MT_VOLUME][MT_VERTEX][elInd][inds[1]];
                  }
                else
                  {
                    MeshTopology_Dim incidence_dim[]=
                      { MT_FACE, MT_FACE, MT_FACE, MT_FACE, MT_FACE, MT_EDGE,
                        MT_EDGE, MT_EDGE, MT_EDGE, MT_FACE };
                    int incidence_index[]={ 2, 3, 1, 1, 3, 9, 6, 10, 7, 3 };
                    set_index(topo,MT_VOLUME,elInd,
                              incidence_dim,incidence_index,inside,dim,ind);
                  }
              }
            else if( bc[3] == 0.0 || bc[3] == -0.0 )
              {
                MeshTopology_Dim incidence_dim[]=
                  { MT_FACE, MT_EDGE, MT_EDGE, MT_EDGE, MT_FACE, MT_EDGE,
                    MT_EDGE, MT_EDGE, MT_FACE, MT_FACE };
                int incidence_index[]={ 0, 1, 4, 6, 4, 3, 4, 1, 1, 5 };
                set_index(topo,MT_VOLUME,elInd,
                          incidence_dim,incidence_index,inside,dim,ind);
              }
            else
              {
                MeshTopology_Dim incidence_dim[]=
                  { MT_VOLUME, MT_FACE, MT_FACE, MT_FACE, MT_VOLUME, MT_FACE,
                    MT_FACE, MT_FACE, MT_FACE, MT_VOLUME };
                int incidence_index[]={ 0, 3, 1, 1, 0, 5, 1, 3, 1, 0 };
                set_index(topo,MT_VOLUME,elInd,
                          incidence_dim,incidence_index,inside,dim,ind);
              }
          }
        else if( bc[1] == 0.0 || bc[1] == -0.0 )
          {
            if( bc[2] == 0.0 || bc[2] == -0.0 )
              {
                if( bc[3] == 0.0 || bc[3] == -0.0 )
                  {
                    *dim = MT_VERTEX;
                    *ind = topo->incEls[MT_VOLUME][MT_VERTEX][elInd][inds[0]];
                  }
                else
                  {
                    MeshTopology_Dim incidence_dim[]=
                      { MT_EDGE, MT_FACE, MT_FACE, MT_FACE, MT_FACE, MT_FACE,
                        MT_FACE, MT_FACE, MT_EDGE, MT_FACE };
                    int incidence_index[]={ 8, 5, 5, 3, 5, 2, 4, 4, 11, 4 };
                    set_index(topo,MT_VOLUME,elInd,
                              incidence_dim,incidence_index,inside,dim,ind);
                  }
              }
            else if( bc[3] == 0.0 || bc[3] == -0.0 )
              {
                MeshTopology_Dim incidence_dim[]=
                  { MT_EDGE, MT_EDGE, MT_EDGE, MT_EDGE, MT_FACE, MT_FACE,
                    MT_FACE, MT_FACE, MT_FACE, MT_FACE };
                int incidence_index[]={ 2, 3, 9, 10, 2, 0, 2, 0, 3, 2 };
                set_index(topo,MT_VOLUME,elInd,
                          incidence_dim,incidence_index,inside,dim,ind);
              }
            else
              {
                MeshTopology_Dim incidence_dim[]=
                  { MT_FACE, MT_FACE, MT_FACE, MT_FACE, MT_VOLUME, MT_VOLUME,
                    MT_VOLUME, MT_VOLUME, MT_FACE, MT_VOLUME };
                int incidence_index[]={ 4, 5, 5, 3, 0, 0, 0, 0, 3, 0};
                set_index(topo,MT_VOLUME,elInd,
                          incidence_dim,incidence_index,inside,dim,ind);
              }
          }
        else if( bc[2] == 0.0 || bc[2] == -0.0 )
          {
            if( bc[3] == 0.0 || bc[3] == -0.0 )
              {
                MeshTopology_Dim incidence_dim[]=
                  { MT_EDGE, MT_FACE, MT_FACE, MT_FACE, MT_FACE, MT_EDGE,
                    MT_EDGE, MT_EDGE, MT_FACE, MT_FACE };
                int incidence_index[]={ 0, 0, 2, 4, 0, 0, 8, 2, 5, 0 };
                set_index(topo,MT_VOLUME,elInd,
                          incidence_dim,incidence_index,inside,dim,ind);
              }
            else
              {
                MeshTopology_Dim incidence_dim[]=
                  { MT_FACE, MT_VOLUME, MT_VOLUME, MT_VOLUME, MT_VOLUME,
                    MT_FACE, MT_FACE, MT_FACE, MT_FACE, MT_VOLUME };
                int incidence_index[]={ 2, 0, 0, 0, 0, 2, 4, 4, 5, 0 };
                set_index(topo,MT_VOLUME,elInd,
                          incidence_dim,incidence_index,inside,dim,ind);
              }
          }
        else if( bc[3] == 0.0 || bc[3] == -0.0 )
          {
            MeshTopology_Dim incidence_dim[]=
              { MT_FACE, MT_FACE, MT_FACE, MT_FACE, MT_VOLUME,
                MT_FACE, MT_FACE, MT_FACE, MT_VOLUME, MT_VOLUME };
            int incidence_index[]={ 0, 0, 2, 4, 0, 0, 2, 0, 0, 0 };
            set_index(topo,MT_VOLUME,elInd,
                      incidence_dim,incidence_index,inside,dim,ind);
          }
        else
          {
            *dim = MT_VOLUME;
            *ind = elInd;
          }
      }
    }
  return fnd;
}

namespace {
  bool Side_Search2D(const int *inc, double **verts, 
                     const double* point, int on_edge[4]);
}

bool Mesh_HexType_ElementHasPoint2DGeneral
( Mesh_HexType* self, unsigned elInd, double* point, 
  MeshTopology_Dim* dim, unsigned* ind )
{
  Mesh*		mesh;
  bool		fnd;
  double		bc[3];
  unsigned	inside;
  const int*	inc;

  assert( self && Stg_CheckType( self, Mesh_HexType ) );
  assert( Mesh_GetDimSize( self->mesh ) == 2 );
  assert( elInd < Mesh_GetDomainSize( self->mesh,
                                      Mesh_GetDimSize( self->mesh ) ) );
  assert( point );
  assert( dim );
  assert( ind );

  /* Shortcuts. */
  mesh = self->mesh;

  /* Get element to vertex incidence. */
  Mesh_GetIncidence( mesh, Mesh_GetDimSize( mesh ),
                     elInd, MT_VERTEX, self->incArray );
  inc = IArray_GetPtr( self->incArray );

  /* If a quadratic element, use side search instead of simplexes. */
  if(IArray_GetSize(self->incArray)==9)
    {
      int on_edge[4];
      fnd=Side_Search2D(inc,mesh->verts,point,on_edge);
    }
  /* Search for triangle. */
  else
    fnd = Simplex_Search2D( mesh->verts, (unsigned*)inc, 2, self->triInds,
                            point, bc, &inside );
  if( fnd ) {
    *dim = MT_FACE;
    *ind = elInd;
  }
  return fnd;
}

bool Mesh_HexType_ElementHasPoint2DWithIncidence
( Mesh_HexType* self, unsigned elInd, double* point, 
  MeshTopology_Dim* dim, unsigned* ind )
{
  Mesh*		mesh;
  bool		fnd;
  double		bc[3];
  IGraph*		topo;
  unsigned	inside;
  const int*	inc;

  assert( self && Stg_CheckType( self, Mesh_HexType ) );
  assert( Mesh_GetDimSize( self->mesh ) == 2 );
  assert( elInd < Mesh_GetDomainSize(self->mesh,Mesh_GetDimSize(self->mesh)));
  assert( point );
  assert( dim );
  assert( ind );

  /* Shortcuts. */
  mesh = self->mesh;
  topo = (IGraph*)mesh->topo;

  /* Get element to vertex incidence. */
  Mesh_GetIncidence( mesh, Mesh_GetDimSize( mesh ), elInd, MT_VERTEX,
                     self->incArray );
  inc = IArray_GetPtr( self->incArray );

  /* If a quadratic element, use side search instead of simplexes. */
  if(IArray_GetSize(self->incArray)==9)
    {
      int on_edge[4];
      fnd=Side_Search2D(inc,mesh->verts,point,on_edge);

      if(fnd)
        {
          if(on_edge[0]==0)
            {
              if(on_edge[2]==0)
                {
                  *dim=MT_VERTEX;
                  *ind=topo->incEls[MT_FACE][MT_VERTEX][elInd][0];
                }
              else if(on_edge[3]==0)
                {
                  *dim=MT_VERTEX;
                  *ind=topo->incEls[MT_FACE][MT_VERTEX][elInd][2];
                }
              else
                {
                  *dim=MT_EDGE;
                  *ind = topo->incEls[MT_FACE][MT_EDGE][elInd][0];
                }
            }
          else if(on_edge[1]==0)
            {
              if(on_edge[2]==0)
                {
                  *dim=MT_VERTEX;
                  *ind=topo->incEls[MT_FACE][MT_VERTEX][elInd][6];
                }
              else if(on_edge[3]==0)
                {
                  *dim=MT_VERTEX;
                  *ind=topo->incEls[MT_FACE][MT_VERTEX][elInd][8];
                }
              else
                {
                  *dim=MT_EDGE;
                  *ind = topo->incEls[MT_FACE][MT_EDGE][elInd][1];
                }
            }
          else if(on_edge[2]==0)
            {
              *dim=MT_EDGE;
              *ind = topo->incEls[MT_FACE][MT_EDGE][elInd][2];
            }
          else if(on_edge[3]==0)
            {
              *dim=MT_EDGE;
              *ind = topo->incEls[MT_FACE][MT_EDGE][elInd][3];
            }
          else
            {
              *dim=MT_FACE;
              *ind=elInd;
            }
        }
    }
  else
    {
      /* Search for triangle. */
      fnd = Simplex_Search2D( mesh->verts, (unsigned*)inc, 2, self->triInds, 
                              point, bc, &inside );
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
          }
        }
        else if( bc[1] == 0.0 || bc[1] == -0.0 ) {
          if( bc[2] == 0.0 || bc[2] == -0.0 ) {
            *dim = MT_VERTEX;
            *ind = topo->incEls[MT_FACE][MT_VERTEX][elInd][inds[0]];
          }
          else {
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
          }
        }
        else if( bc[2] == 0.0 || bc[2] == -0.0 ) {
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
        }
        else {
          *dim = MT_FACE;
          *ind = elInd;
        }
      }
    }
  return fnd;
}

Bool Mesh_HexType_ElementHasPoint1DGeneral( Mesh_HexType* self, unsigned elInd, double* point, 
					    MeshTopology_Dim* dim, unsigned* ind )
{
	Mesh*		mesh;
	const int*	inc;

	assert( self && Stg_CheckType( self, Mesh_HexType ) );
	assert( Mesh_GetDimSize( self->mesh ) == 1 );
	assert( elInd < Mesh_GetDomainSize( self->mesh, Mesh_GetDimSize( self->mesh ) ) );
	assert( point );
	assert( dim );
	assert( ind );

	mesh = self->mesh;
	Mesh_GetIncidence( mesh, MT_EDGE, elInd, MT_VERTEX, self->incArray );
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
	const int*	inc;

	assert( self && Stg_CheckType( self, Mesh_HexType ) );
	assert( Mesh_GetDimSize( self->mesh ) == 1 );
	assert( elInd < Mesh_GetDomainSize( self->mesh, Mesh_GetDimSize( self->mesh ) ) );
	assert( point );
	assert( dim );
	assert( ind );

	mesh = self->mesh;
	Mesh_GetIncidence( mesh, MT_EDGE, elInd, MT_VERTEX, self->incArray );
	assert(2 == IArray_GetSize( self->incArray ));
	inc = IArray_GetPtr( self->incArray );

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


namespace {

  int is_face_over_point(const int mm, const int mz, const int mp,
                         const int zm, const int zz, const int zp,
                         const int pm, const int pz, const int pp,
                         double **verts, const int z, const double *point);

  bool about_equal(const double &a, const double &b)
  {
    return (fabs(a)>1 && Num_Approx(1,b/a))
      || Num_Approx(a,b);
  }

  bool about_equal(const double &a, const double &b, const double &c)
  {
    return (fabs(a)>1 && Num_Approx(1,b/a) && Num_Approx(1,c/a))
      || (Num_Approx(a,b) && Num_Approx(a,c));
  }

  bool Face_Search3D(const int *inc, double **verts, 
                     const double* point, int on_face[6])
  {
    on_face[0]=-1*is_face_over_point(inc[0],inc[1],inc[2],inc[3],inc[4],inc[5],
                                     inc[6],inc[7],inc[8],verts,2,point);
    on_face[1]=is_face_over_point(inc[18],inc[19],inc[20],inc[21],inc[22],
                                  inc[23],inc[24],inc[25],inc[26],
                                  verts,2,point);
    on_face[2]=-1*is_face_over_point(inc[0],inc[9],inc[18],inc[1],inc[10],
                                     inc[19],inc[2],inc[11],inc[20],
                                     verts,1,point);
    on_face[3]=is_face_over_point(inc[6],inc[15],inc[24],inc[7],inc[16],
                                  inc[25],inc[8],inc[17],inc[26],
                                  verts,1,point);
    on_face[4]=-1*is_face_over_point(inc[0],inc[3],inc[6],inc[9],inc[12],
                                     inc[15],inc[18],inc[21],inc[24],
                                     verts,0,point);
    on_face[5]=is_face_over_point(inc[2],inc[5],inc[8],inc[11],inc[14],
                                  inc[17],inc[20],inc[23],inc[26],
                                  verts,0,point);
    return (on_face[0]>=0) && (on_face[1]>=0) && (on_face[2]>=0)
      && (on_face[3]>=0) && (on_face[4]>=0) && (on_face[5]>=0);
  }

  void set_N(const double &xi, double N[3])
  {
    N[0]=xi*(xi-1)/2;
    N[1]=(1-xi*xi);
    N[2]=xi*(xi+1)/2;
  }

  template <typename T> int sgn(T val)
  {
    return (val > T(0)) - (val < T(0));
  }

  int is_face_over_point(const int mm, const int zm, const int pm,
                         const int mz, const int zz, const int pz,
                         const int mp, const int zp, const int pp,
                         double **verts, const int z, const double *point)
  {
    int x((z+1)%3), y((z+2)%3);
    /* First, find which walls are vertical and evenly spaced. */

    bool x_equal=about_equal(verts[mz][x],verts[mm][x],verts[mp][x])
      && about_equal(verts[zz][x],verts[zm][x],verts[zp][x])
      && about_equal(verts[pz][x],verts[pm][x],verts[pp][x])
      && about_equal(verts[pz][x]-verts[zz][x],verts[zz][x]-verts[mz][x]);

    bool y_equal=about_equal(verts[zm][y],verts[mm][y],verts[pm][y])
      && about_equal(verts[zz][y],verts[mz][y],verts[pz][y])
      && about_equal(verts[zp][y],verts[mp][y],verts[pp][y])
      && about_equal(verts[zp][y]-verts[zz][y],verts[zz][y]-verts[zm][y]);

    Journal_Firewall(x_equal || y_equal,
                     Journal_Register( Error_Type,"Mesh_HexType"),
                     "The coordinates in either direction %d or %d must be lined up and evenly spaced.\n\t(%g %g %g)\n\t(%g %g %g)\n\t(%g %g %g)\n\t(%g %g %g)\n\t(%g %g %g)\n\t(%g %g %g)\n\t(%g %g %g)\n\t(%g %g %g)\n\t(%g %g %g)\nDid you forget to enable EulerDeform?\n",
                     x,y,verts[mm][0],verts[mm][1],verts[mm][2],
                     verts[mz][0],verts[mz][1],verts[mz][2],
                     verts[mp][0],verts[mp][1],verts[mp][2],
                     verts[zm][0],verts[zm][1],verts[zm][2],
                     verts[zz][0],verts[zz][1],verts[zz][2],
                     verts[zp][0],verts[zp][1],verts[zp][2],
                     verts[pm][0],verts[pm][1],verts[pm][2],
                     verts[pz][0],verts[pz][1],verts[pz][2],
                     verts[pp][0],verts[pp][1],verts[pp][2]);

    /* Compute xi and eta */

    const double dx[]={ (verts[pm][x]-verts[mm][x])/2,
                        (verts[pz][x]-verts[mz][x])/2,
                        (verts[pp][x]-verts[mp][x])/2 };

    const double dy[]={ (verts[mp][y]-verts[mm][y])/2,
                        (verts[zp][y]-verts[zm][y])/2,
                        (verts[pp][y]-verts[pm][y])/2 };

    double xi, eta, Nx[3], Ny[3];

    if(x_equal)
      {
        xi=(point[x]-verts[zz][x])/dx[1];
        set_N(xi,Nx);

        eta=(point[y]
             - (Nx[0]*verts[mz][y] + Nx[1]*verts[zz][y] + Nx[2]*verts[pz][y]))/
        (Nx[0]*dy[0] + Nx[1]*dy[1] + Nx[2]*dy[2]);
        set_N(eta,Ny);
      }
    else
      {
        eta=(point[y]-verts[zz][y])/dy[1];
        set_N(eta,Ny);

        xi=(point[x]
             - (Ny[0]*verts[zm][y] + Ny[1]*verts[zz][y] + Ny[2]*verts[zp][y]))/
        (Ny[0]*dx[0] + Ny[1]*dx[1] + Ny[2]*dx[2]);
        set_N(xi,Nx);
      }

    /* Compute z */
    double computed_z=
      Nx[0]*(Ny[0]*verts[mm][z] + Ny[1]*verts[mz][z] + Ny[2]*verts[mp][z])
      + Nx[1]*(Ny[0]*verts[zm][z] + Ny[1]*verts[zz][z] + Ny[2]*verts[zp][z])
      + Nx[2]*(Ny[0]*verts[pm][z] + Ny[1]*verts[pz][z] + Ny[2]*verts[pp][z]);

    if(about_equal(computed_z,point[z]))
      return 0;
    else
      return sgn(computed_z-point[z]);
  }



  /* Check whether a point is inside an quadratic element by checking
     whether it is above or below the 4 edges of the element.  For now,
     only works when the spacing on each edge is even.  Currently, this
     is enforced by EulerDeform. */

  int is_line_over_point(const int minus, const int zero, const int plus,
                         double **verts, const int x, const double *point);

  bool Side_Search2D(const int *inc, double **verts, 
                     const double* point, int on_edge[4])
  {
    on_edge[0]=-1*is_line_over_point(inc[0],inc[1],inc[2],verts,0,point);
    on_edge[1]=is_line_over_point(inc[6],inc[7],inc[8],verts,0,point);
    on_edge[2]=-1*is_line_over_point(inc[0],inc[3],inc[6],verts,1,point);
    on_edge[3]=is_line_over_point(inc[8],inc[5],inc[2],verts,1,point);

    return (on_edge[0]>=0) && (on_edge[1]>=0) && (on_edge[2]>=0)
      && (on_edge[3]>=0);
  }

  /* Checks whether a point is above, below, or on a line.  To reverse
     the direction, reverse the order of indces (plus, zero, minus).
     Only works for when the spacing between the "x" points within the
     element is the same.  This is enforced by EulerDeform.  This could
     be generalized to any spacing, but there is no good solution in
     3D.  */

  int is_line_over_point(const int minus, const int zero, const int plus,
                         double **verts, const int x, const double *point)
  {
    int y=(x+1)%2;
  
    double x_0, y_0, dy_p, dy_m, dx, b, c, computed_y;

    x_0=verts[zero][x];
    y_0=verts[zero][y];
    dy_p=verts[plus][y]-y_0;
    dy_m=verts[minus][y]-y_0;

    dx=verts[plus][x]-x_0;
    Journal_Firewall(about_equal(dx,x_0-verts[minus][x]),
                     Journal_Register( Error_Type,"Mesh_HexType"),
                     "The coordinates in direction %d must be evenly spaced.\n\t(%g %g)\n\t(%g %g)\n\t(%g %g)\nDid you forget to enable EulerDeform?\n",
                     x,verts[minus][0],verts[minus][1],
                     verts[zero][0],verts[zero][1],
                     verts[plus][0],verts[plus][1]);
  
    c=(dy_p + dy_m)/(2*dx*dx);
    b=(dy_p - dy_m)/(2*dx);

    computed_y=y_0 + b*(point[x]-x_0) + c*(point[x]-x_0)*(point[x]-x_0);

    if(about_equal(computed_y,point[y]))
      return 0;
    else
      return sgn(computed_y-point[y]);
  }
}
