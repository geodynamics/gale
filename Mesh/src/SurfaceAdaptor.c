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
** $Id: SurfaceAdaptor.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>
#include <StGermain/StGermain.h>

#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>

#include "types.h"
#include "shortcuts.h"
#include "Grid.h"
#include "Decomp.h"
#include "Sync.h"
#include "MeshTopology.h"
#include "IGraph.h"
#include "MeshClass.h"
#include "MeshGenerator.h"
#include "MeshAdaptor.h"
#include "SurfaceAdaptor.h"
#include "Remesher.h"


typedef double (SurfaceAdaptor_DeformFunc)( SurfaceAdaptor* self, Mesh* mesh,
					    unsigned* globalSize, unsigned vertex, unsigned* vertexInds);


/* Textual name of this class */
const Type SurfaceAdaptor_Type = "SurfaceAdaptor";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

SurfaceAdaptor* SurfaceAdaptor_New( Name name, AbstractContext* context ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(SurfaceAdaptor);
	Type                                                      type = SurfaceAdaptor_Type;
	Stg_Class_DeleteFunction*                              _delete = _SurfaceAdaptor_Delete;
	Stg_Class_PrintFunction*                                _print = _SurfaceAdaptor_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = (void* (*)(Name))_SurfaceAdaptor_New;
	Stg_Component_ConstructFunction*                    _construct = _SurfaceAdaptor_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _SurfaceAdaptor_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _SurfaceAdaptor_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _SurfaceAdaptor_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _SurfaceAdaptor_Destroy;
	AllocationType                              nameAllocationType = NON_GLOBAL;
	MeshGenerator_SetDimSizeFunc*                   setDimSizeFunc = _MeshGenerator_SetDimSize;
	MeshGenerator_GenerateFunc*                       generateFunc = SurfaceAdaptor_Generate;

	SurfaceAdaptor* self = _SurfaceAdaptor_New(  SURFACEADAPTOR_PASSARGS  );

   _MeshGenerator_Init( (MeshGenerator*)self, context );
   _MeshAdaptor_Init( (MeshAdaptor*)self );
	_SurfaceAdaptor_Init( self );

   return self;
}

SurfaceAdaptor* _SurfaceAdaptor_New(  SURFACEADAPTOR_DEFARGS  ) {
	SurfaceAdaptor* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(SurfaceAdaptor) );
	self = (SurfaceAdaptor*)_MeshAdaptor_New(  MESHADAPTOR_PASSARGS  );

	/* Virtual info */
	return self;
}

void _SurfaceAdaptor_Init( SurfaceAdaptor* self ) {
	self->surfaceType = SurfaceAdaptor_SurfaceType_Invalid;
	memset( &self->info, 0, sizeof(SurfaceAdaptor_SurfaceInfo) );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _SurfaceAdaptor_Delete( void* adaptor ) {
	SurfaceAdaptor*	self = (SurfaceAdaptor*)adaptor;

	/* Delete the parent. */
	_MeshAdaptor_Delete( self );
}

void _SurfaceAdaptor_Print( void* adaptor, Stream* stream ) {
	SurfaceAdaptor*	self = (SurfaceAdaptor*)adaptor;
	
	/* Set the Journal for printing informations */
	Stream* adaptorStream;
	adaptorStream = Journal_Register( InfoStream_Type, (Name)"SurfaceAdaptorStream"  );

	/* Print parent */
	Journal_Printf( stream, "SurfaceAdaptor (ptr): (%p)\n", self );
	_MeshAdaptor_Print( self, stream );
}

void _SurfaceAdaptor_AssignFromXML( void* adaptor, Stg_ComponentFactory* cf, void* data ) {
	SurfaceAdaptor*	self = (SurfaceAdaptor*)adaptor;
	Dictionary*	dict;
	char*		surfaceType;

	assert( self );
	assert( cf );

	/* Call parent construct. */
	_MeshAdaptor_AssignFromXML( self, cf, data );

	/* Rip out the components structure as a dictionary. */
	dict = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( cf->componentDict, (Dictionary_Entry_Key)self->name )  );

        /* Check if we want to keep a certain depth at the bottom reserved for
           contact elements. */
        self->contactDepth = Stg_ComponentFactory_GetInt( cf, self->name, (Dictionary_Entry_Key)"contactDepth", 0  );

	/* What kind of surface do we want? */
	surfaceType = Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"surfaceType", ""  );
	if( !strcmp( surfaceType, "wedge" ) ) {
		self->surfaceType = SurfaceAdaptor_SurfaceType_Wedge;
		self->info.wedge.offs[0] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"beginOffset", 0.0  );
		self->info.wedge.endOffs[0] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"endOffset", 1.0  );
		self->info.wedge.grad[0] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"gradient", 0.5  );
		/* get the parameters for the z-axis */
		self->info.wedge.offs[1] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"beginOffsetZ", 0.0  );
		self->info.wedge.endOffs[1] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"endOffsetZ", 1.0  );
		self->info.wedge.grad[1] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"gradientZ", 0.5  );
	}
	else if( !strcmp( surfaceType, "plateau" ) ) {
		self->surfaceType = SurfaceAdaptor_SurfaceType_Plateau;
		self->info.plateau.x1 = Stg_ComponentFactory_GetDouble( cf, self->name, "x1", 0.0 );
		self->info.plateau.x2 = Stg_ComponentFactory_GetDouble( cf, self->name, "x2", 0.0 );
		self->info.plateau.x3 = Stg_ComponentFactory_GetDouble( cf, self->name, "x3", 0.0 );
		self->info.plateau.x4 = Stg_ComponentFactory_GetDouble( cf, self->name, "x4", 0.0 );
		self->info.plateau.z1 = Stg_ComponentFactory_GetDouble( cf, self->name, "z1", 0.0 );
		self->info.plateau.z2 = Stg_ComponentFactory_GetDouble( cf, self->name, "z2", 0.0 );
		self->info.plateau.z3 = Stg_ComponentFactory_GetDouble( cf, self->name, "z3", 0.0 );
		self->info.plateau.z4 = Stg_ComponentFactory_GetDouble( cf, self->name, "z4", 0.0 );
		self->info.plateau.height = Stg_ComponentFactory_GetDouble( cf, self->name, "height", 0.0 );
	}
	else if( !strcmp( surfaceType, "topo_data" ) ) {
                FILE *fp;
                char* surfaceName;
                int i,j,ii,jj;
                surfaceName = Stg_ComponentFactory_GetString( cf, self->name, "surfaceName", "ascii_topo" );
		self->info.topo_data.nx = Stg_ComponentFactory_GetInt( cf, self->name, "nx", 0 );
		self->info.topo_data.nz = Stg_ComponentFactory_GetInt( cf, self->name, "nz", 0 );
	        self->info.topo_data.minX = Stg_ComponentFactory_GetDouble( cf, self->name, "minX", 0 );
	        self->info.topo_data.minZ = Stg_ComponentFactory_GetDouble( cf, self->name, "minZ", 0 );
	        self->info.topo_data.maxX = Stg_ComponentFactory_GetDouble( cf, self->name, "maxX", 0 );
	        self->info.topo_data.maxZ = Stg_ComponentFactory_GetDouble( cf, self->name, "maxZ", 0 );
                self->info.topo_data.dx=
                  (self->info.topo_data.maxX-self->info.topo_data.minX)
                  /(self->info.topo_data.nx-1);
                self->info.topo_data.dz=
                  (self->info.topo_data.maxZ-self->info.topo_data.minZ)
                  /(self->info.topo_data.nz-1);
                self->info.topo_data.heights=
                  malloc(sizeof(double)*self->info.topo_data.nx
                         *self->info.topo_data.nz);
		self->surfaceType = SurfaceAdaptor_SurfaceType_Topo_Data;
                fp=fopen(surfaceName,"r");
                if(!fp)
                  {
                    printf("Can not open the file %s\n",surfaceName);
                    abort();
                  }
                for(i=0;i<self->info.topo_data.nx;++i)
                  for(j=0;j<self->info.topo_data.nz;++j)
                    {
                      float h;
                      fscanf(fp,"%d %d %f",&ii,&jj,&h);
                      self->info.topo_data.heights[ii+self->info.topo_data.nx*jj]=h;
                    }
                fclose(fp);
	}
	else if( !strcmp( surfaceType, "sine" ) || !strcmp( surfaceType, "cosine" ) ) {
		Dictionary_Entry_Value*	originList;

		if( !strcmp( surfaceType, "sine" ) )
			self->surfaceType = SurfaceAdaptor_SurfaceType_Sine;
		else
			self->surfaceType = SurfaceAdaptor_SurfaceType_Cosine;

		originList = Dictionary_Get( dict, (Dictionary_Entry_Key)"origin" );
		if( originList ) {
			unsigned	nDims;
			unsigned	d_i;

			nDims = Dictionary_Entry_Value_GetCount( originList );
			for( d_i = 0; d_i < nDims; d_i++  ) {
				Dictionary_Entry_Value*	val;

				val = Dictionary_Entry_Value_GetElement( originList, d_i );
				self->info.trig.origin[d_i] = Dictionary_Entry_Value_AsDouble( val );
			}
		}
		else
			memset( self->info.trig.origin, 0, sizeof(double) * 2 );

		self->info.trig.amp = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"amplitude", 1.0  );
		self->info.trig.freq = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"frequency", 1.0 );
	}
	else
		_SurfaceAdaptor_Init( self  );
}

void _SurfaceAdaptor_Build( void* adaptor, void* data ) {
   SurfaceAdaptor* self = (SurfaceAdaptor*)adaptor;

   _MeshAdaptor_Build( self, data );
}

void _SurfaceAdaptor_Initialise( void* adaptor, void* data ) {
   SurfaceAdaptor* self = (SurfaceAdaptor*)adaptor;

   _MeshAdaptor_Initialise( self, data );
}

void _SurfaceAdaptor_Execute( void* adaptor, void* data ) {
   SurfaceAdaptor* self = (SurfaceAdaptor*)adaptor;

   _MeshAdaptor_Execute( self, data );
}

void _SurfaceAdaptor_Destroy( void* adaptor, void* data ) {
   SurfaceAdaptor* self = (SurfaceAdaptor*)adaptor;

   _MeshAdaptor_Destroy( self, data );
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void SurfaceAdaptor_Generate( void* adaptor, void* _mesh, void* data ) {
	SurfaceAdaptor* self = (SurfaceAdaptor*)adaptor;
	Mesh* mesh = (Mesh*)_mesh;
	const Sync* sync;
	SurfaceAdaptor_DeformFunc* deformFunc;
	Grid *grid;
	unsigned* inds;
	unsigned n_i;

	/* Build base mesh, which is assumed to be cartesian. */
	MeshGenerator_Generate( self->generator, mesh, data );

	/* if loading from checkpoint then forget about this step */
	if( ((Context*)data)->loadFromCheckPoint == True )
		return;
	/* If we're not 2D or 3D, forget about it. */
	if( mesh->topo->nDims != 2 && mesh->topo->nDims != 3 )
		return;

	/* What kind of surface do we want? */
	switch( self->surfaceType ) {
	case SurfaceAdaptor_SurfaceType_Wedge:
		if ( mesh->topo->nDims == 3 ) { deformFunc = SurfaceAdaptor_Wedge3D ; }
		else { deformFunc = SurfaceAdaptor_Wedge2D; }
		break;
	case SurfaceAdaptor_SurfaceType_Plateau:
		deformFunc = SurfaceAdaptor_Plateau;
		break;
	case SurfaceAdaptor_SurfaceType_Topo_Data:
		deformFunc = SurfaceAdaptor_Topo_Data;
		break;
	case SurfaceAdaptor_SurfaceType_Sine:
		deformFunc = SurfaceAdaptor_Sine;
		break;
	case SurfaceAdaptor_SurfaceType_Cosine:
		deformFunc = SurfaceAdaptor_Cosine;
		break;
	default:
		break;
	};

	/* Extract the cartesian information. */
	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, (Name)"vertexGrid" )  );
	inds = AllocArray( unsigned, Mesh_GetDimSize( mesh ) );

	/* Loop over domain nodes. */
	sync = IGraph_GetDomain( mesh->topo, MT_VERTEX );
	for( n_i = 0; n_i < Sync_GetNumDomains( sync ); n_i++ ) {
		unsigned	gNode;
		double		height;
                double deform;

		gNode = Sync_DomainToGlobal( sync, n_i );
		Grid_Lift( grid, gNode, inds );

		/* Check if we're inside the contact depth. */
		if( inds[1] <= self->contactDepth )
		   continue;

		/* Calculate a height percentage. */
		height = (double)(inds[1] - self->contactDepth) / (double)(grid->sizes[1] - 1);

		/* Deform this node. */
                deform = deformFunc( self, mesh, grid->sizes, n_i, inds);
		mesh->verts[n_i][1] += height * deform;
	}

	/* Free resources. */
	FreeArray( inds );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

double SurfaceAdaptor_Wedge2D( SurfaceAdaptor* self, Mesh* mesh, 
			     unsigned* globalSize, unsigned vertex, unsigned* vertexInds )
{
   if( mesh->verts[vertex][0] >= self->info.wedge.offs[0] ) {
      if( mesh->verts[vertex][0] >= self->info.wedge.endOffs[0] )
         return (self->info.wedge.endOffs[0] - self->info.wedge.offs[0]) * self->info.wedge.grad[0];
      else
         return (mesh->verts[vertex][0] - self->info.wedge.offs[0]) * self->info.wedge.grad[0];
   }
   else 
      return 0.0;
}

double SurfaceAdaptor_Wedge3D( SurfaceAdaptor* self, Mesh* mesh, 
			     unsigned* globalSize, unsigned vertex, unsigned* vertexInds )
{
   if( mesh->verts[vertex][0] >= self->info.wedge.offs[0] ) {
      if( mesh->verts[vertex][0] >= self->info.wedge.endOffs[0] ) {
         return (self->info.wedge.endOffs[0] - self->info.wedge.offs[0]) * self->info.wedge.grad[0] + 
					 (mesh->verts[vertex][2] - self->info.wedge.offs[1]) * self->info.wedge.grad[1];
			} else {
         return (mesh->verts[vertex][0] - self->info.wedge.offs[0]) * self->info.wedge.grad[0] +
					 (mesh->verts[vertex][2] - self->info.wedge.offs[1]) * self->info.wedge.grad[1];
			}
   }
   else 
      return 0.0;
}
double SurfaceAdaptor_Plateau( SurfaceAdaptor* self, Mesh* mesh, 
                               unsigned* globalSize, unsigned vertex,
                               unsigned* vertexInds )
{
  double x_factor, z_factor;
  x_factor =1;
  z_factor=1;
  if( mesh->verts[vertex][0] < self->info.plateau.x1
      || mesh->verts[vertex][0] > self->info.plateau.x4)
    {
      x_factor=0;
    }
  else if( mesh->verts[vertex][0] <= self->info.plateau.x2)
    {
      x_factor=(mesh->verts[vertex][0] - self->info.plateau.x1)
        /(self->info.plateau.x2 - self->info.plateau.x1);
    }
  else if( mesh->verts[vertex][0] <= self->info.plateau.x3)
    {
      x_factor=1;
    }
  else if( mesh->verts[vertex][0] <= self->info.plateau.x4)
    {
      x_factor=(self->info.plateau.x4 - mesh->verts[vertex][0])
        /(self->info.plateau.x4 - self->info.plateau.x3);
    }

  if(mesh->topo->nDims==3)
    {
      if( mesh->verts[vertex][2] < self->info.plateau.z1
          || mesh->verts[vertex][2] > self->info.plateau.z4)
        {
          z_factor=0;
        }
      else if( mesh->verts[vertex][2] <= self->info.plateau.z2)
        {
          z_factor=(mesh->verts[vertex][2] - self->info.plateau.z1)
            /(self->info.plateau.z2 - self->info.plateau.z1);
        }
      else if( mesh->verts[vertex][2] <= self->info.plateau.z3)
        {
          z_factor=1;
        }
      else if( mesh->verts[vertex][2] <= self->info.plateau.z4)
        {
          z_factor=(self->info.plateau.z4 - mesh->verts[vertex][2])
            /(self->info.plateau.z4 - self->info.plateau.z3);
        }
    }

  return x_factor*z_factor*self->info.plateau.height;
}

double SurfaceAdaptor_Topo_Data( SurfaceAdaptor* self, Mesh* mesh, 
                                 unsigned* globalSize, unsigned vertex,
                                 unsigned* vertexInds )
{
  int i,k,ip,kp;
  double dx,dz;

  i=floor((mesh->verts[vertex][0] - self->info.topo_data.minX)
          /self->info.topo_data.dx + 0.5);
  k=floor((mesh->verts[vertex][2] - self->info.topo_data.minZ)
          /self->info.topo_data.dz + 0.5);

  if(i<0 || i>self->info.topo_data.nx-1
     || k<0 || k>self->info.topo_data.nz-1)
    {
      printf("Coordinate not covered by the topography file: %g %g\n\tminX: %g\n\tmaxX: %g\n\tminZ: %g\n\tmaxZ: %g\n\tnx: %d\n\tnz: %d\n",
             mesh->verts[vertex][0],
             mesh->verts[vertex][2],
             self->info.topo_data.minX,
             self->info.topo_data.maxX,
             self->info.topo_data.minZ,
             self->info.topo_data.maxZ,
             self->info.topo_data.nx,
             self->info.topo_data.nz);
      abort();
    }

  /* Interpolate the height */
  ip=i+1;
  kp=k+1;
  if(ip>self->info.topo_data.nx-1)
    ip=i;
  if(kp>self->info.topo_data.nz-1)
    kp=k;

  dx=(mesh->verts[vertex][0]
      - (i*self->info.topo_data.dx+self->info.topo_data.minX))
    /self->info.topo_data.dx;
  dz=(mesh->verts[vertex][2]
      - (k*self->info.topo_data.dz+self->info.topo_data.minZ))
    /self->info.topo_data.dz;

  return self->info.topo_data.heights[i+self->info.topo_data.nx*k]*(1-dx)*(1-dz)
    + self->info.topo_data.heights[i+self->info.topo_data.nx*kp]*(1-dx)*dz
    + self->info.topo_data.heights[ip+self->info.topo_data.nx*k]*dx*(1-dz)
    + self->info.topo_data.heights[ip+self->info.topo_data.nx*kp]*dx*dz;
}

double SurfaceAdaptor_Sine( SurfaceAdaptor* self, Mesh* mesh, 
			    unsigned* globalSize, unsigned vertex, unsigned* vertexInds )
{
	double	dx, dy;
	double	rad;

	dx = mesh->verts[vertex][0] - self->info.trig.origin[0];
	rad = dx * dx;
	if( mesh->topo->nDims == 3 ) {
		dy = mesh->verts[vertex][1] - self->info.trig.origin[1];
		rad += dy * dy;
	}
	rad = sqrt( rad );

	return self->info.trig.amp * sin( self->info.trig.freq * rad );
}

double SurfaceAdaptor_Cosine( SurfaceAdaptor* self, Mesh* mesh, 
			      unsigned* globalSize, unsigned vertex, unsigned* vertexInds )
{
	double	dx, dz;
	double	rad;

	dx = mesh->verts[vertex][0] - self->info.trig.origin[0];
	rad = dx * dx;
	if( mesh->topo->nDims == 3 ) {
		dz = mesh->verts[vertex][2] - self->info.trig.origin[2];
		rad += dz * dz;
	}
	rad = sqrt( rad );

	return self->info.trig.amp * cos( self->info.trig.freq * rad );
}


