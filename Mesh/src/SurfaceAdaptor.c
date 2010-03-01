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
#include <StgDomain/Mesh/Mesh.h>

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
	self->topSurfaceType = SurfaceAdaptor_SurfaceType_Invalid;
	self->bottomSurfaceType = SurfaceAdaptor_SurfaceType_Invalid;
	memset( &self->top_info, 0, sizeof(SurfaceAdaptor_SurfaceInfo ));
	memset( &self->bottom_info, 0, sizeof(SurfaceAdaptor_SurfaceInfo));
        self->topDeformFunc=self->bottomDeformFunc=NULL;
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

void _SurfaceAdaptor_AssignFromXML_Surface(Stg_ComponentFactory* cf,
                                           Name name,
                                           SurfaceAdaptor_SurfaceType *surfaceType,
                                           Dictionary* dict,
                                           SurfaceAdaptor_SurfaceInfo *info,
                                           SurfaceAdaptor_DeformFunc **deformFunc,
                                           char *surface);

void _SurfaceAdaptor_AssignFromXML( void* adaptor, Stg_ComponentFactory* cf, void* data ) {
	SurfaceAdaptor*	self = (SurfaceAdaptor*)adaptor;
	Dictionary*	dict;

	assert( self );
	assert( cf );

	/* Call parent construct. */
	_MeshAdaptor_AssignFromXML( self, cf, data );

	/* Rip out the components structure as a dictionary. */
	dict = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( cf->componentDict, (Dictionary_Entry_Key)self->name )  );

        self->topDeformFunc=self->bottomDeformFunc=NULL;
        char surface[100];
        _SurfaceAdaptor_AssignFromXML_Surface(cf,self->name,
                                              &(self->topSurfaceType),
                                              dict,
                                              &(self->top_info),
                                              &(self->topDeformFunc),
                                              strcpy(surface,"top"));
        _SurfaceAdaptor_AssignFromXML_Surface(cf,self->name,
                                              &(self->bottomSurfaceType),
                                              dict,
                                              &(self->bottom_info),
                                              &(self->bottomDeformFunc),
                                              strcpy(surface,"bottom"));
}

void _SurfaceAdaptor_AssignFromXML_Surface(Stg_ComponentFactory* cf,
                                           Name name,
                                           SurfaceAdaptor_SurfaceType *surfaceType,
                                           Dictionary* dict,
                                           SurfaceAdaptor_SurfaceInfo *info,
                                           SurfaceAdaptor_DeformFunc **deformFunc,
                                           char *surface)
{
  char temp[100];
  char *surfaceName;
  strcpy(temp,surface);
  /* What kind of surface do we want? */
  surfaceName = 
    Stg_ComponentFactory_GetString( cf, name,
                                    (Dictionary_Entry_Key)strcat(temp,"SurfaceType"), ""  );
  if( !strcmp( surfaceName, "wedge" ) ) {
    *surfaceType = SurfaceAdaptor_SurfaceType_Wedge;
    *deformFunc = SurfaceAdaptor_Wedge ;
    strcpy(temp,surface);
    info->wedge.offs[0] =
      Stg_ComponentFactory_GetDouble( cf, name,
                                      (Dictionary_Entry_Key)strcat(temp,"BeginOffset"), 0.0  );
    strcpy(temp,surface);
    info->wedge.endOffs[0] =
      Stg_ComponentFactory_GetDouble( cf, name,
                                      (Dictionary_Entry_Key)strcat(temp,"EndOffset"), 1.0  );
    strcpy(temp,surface);
    info->wedge.grad[0] = 
      Stg_ComponentFactory_GetDouble( cf, name,
                                      (Dictionary_Entry_Key)strcat(temp,"Gradient"), 0.5  );
    /* get the parameters for the z-axis */
    strcpy(temp,surface);
    info->wedge.offs[1] = 
      Stg_ComponentFactory_GetDouble( cf, name,
                                      (Dictionary_Entry_Key)strcat(temp,"BeginOffsetZ"), 0.0  );
    strcpy(temp,surface);
    info->wedge.endOffs[1] = 
      Stg_ComponentFactory_GetDouble( cf, name,
                                      (Dictionary_Entry_Key)strcat(temp,"EndOffsetZ"), 1.0  );
    strcpy(temp,surface);
    info->wedge.grad[1] = 
      Stg_ComponentFactory_GetDouble( cf, name,
                                      (Dictionary_Entry_Key)strcat(temp,"GradientZ"), 0.5  );
  }
  else if( !strcmp( surfaceName, "plateau" ) ) {
    *surfaceType = SurfaceAdaptor_SurfaceType_Plateau;
    *deformFunc = SurfaceAdaptor_Plateau;
    strcpy(temp,surface);
    info->plateau.x1 =
      Stg_ComponentFactory_GetDouble( cf, name,
                                      strcat(temp,"X1"), 0.0 );
    strcpy(temp,surface);
    info->plateau.x2 = 
      Stg_ComponentFactory_GetDouble( cf, name,
                                      strcat(temp,"X2"), 0.0 );
    strcpy(temp,surface);
    info->plateau.x3 = 
      Stg_ComponentFactory_GetDouble( cf, name,
                                      strcat(temp,"X3"), 0.0 );
    strcpy(temp,surface);
    info->plateau.x4 = 
      Stg_ComponentFactory_GetDouble( cf, name,
                                      strcat(temp,"X4"), 0.0 );
    strcpy(temp,surface);
    info->plateau.z1 = 
      Stg_ComponentFactory_GetDouble( cf, name,
                                      strcat(temp,"Z1"), 0.0 );
    strcpy(temp,surface);
    info->plateau.z2 = 
      Stg_ComponentFactory_GetDouble( cf, name,
                                      strcat(temp,"Z2"), 0.0 );
    strcpy(temp,surface);
    info->plateau.z3 = 
      Stg_ComponentFactory_GetDouble( cf, name,
                                      strcat(temp,"Z3"), 0.0 );
    strcpy(temp,surface);
    info->plateau.z4 = 
      Stg_ComponentFactory_GetDouble( cf, name,
                                      strcat(temp,"Z4"), 0.0 );
    strcpy(temp,surface);
    info->plateau.height =
      Stg_ComponentFactory_GetDouble( cf, name,
                                      strcat(temp,"Height"), 0.0 );
  }
  else if( !strcmp( surfaceName, "topo_data" ) ) {
    FILE *fp;
    char* surfaceFileName;
    int i,j,ii,jj;
    *surfaceType = SurfaceAdaptor_SurfaceType_Topo_Data;
    *deformFunc = SurfaceAdaptor_Topo_Data;
    strcpy(temp,surface);
    surfaceFileName =
      Stg_ComponentFactory_GetString( cf, name,
                                      strcat(temp,"SurfaceName"),
                                      "ascii_topo" );
    strcpy(temp,surface);
    info->topo_data.nx =
      Stg_ComponentFactory_GetInt( cf, name,
                                   strcat(temp,"Nx"), 0 );
    strcpy(temp,surface);
    info->topo_data.nz = 
      Stg_ComponentFactory_GetInt( cf, name,
                                   strcat(temp,"Nz"), 0 );
    strcpy(temp,surface);
    info->topo_data.minX = 
      Stg_ComponentFactory_GetDouble( cf, name,
                                      strcat(temp,"MinX"), 0 );
    strcpy(temp,surface);
    info->topo_data.minZ = 
      Stg_ComponentFactory_GetDouble( cf, name,
                                      strcat(temp,"MinZ"), 0 );
    strcpy(temp,surface);
    info->topo_data.maxX = 
      Stg_ComponentFactory_GetDouble( cf, name,
                                      strcat(temp,"MaxX"), 0 );
    strcpy(temp,surface);
    info->topo_data.maxZ =
      Stg_ComponentFactory_GetDouble( cf, name,
                                      strcat(temp,"MaxZ"), 0 );
    info->topo_data.dx=
      (info->topo_data.maxX-info->topo_data.minX)
      /(info->topo_data.nx-1);
    info->topo_data.dz=
      (info->topo_data.maxZ-info->topo_data.minZ)
      /(info->topo_data.nz-1);
    info->topo_data.heights=
      malloc(sizeof(double)*info->topo_data.nx
             *info->topo_data.nz);
    fp=fopen(surfaceFileName,"r");
    if(!fp)
      {
        printf("Can not open the file %s\n",surfaceFileName);
        abort();
      }
    for(i=0;i<info->topo_data.nx;++i)
      for(j=0;j<info->topo_data.nz;++j)
        {
          float h;
          fscanf(fp,"%d %d %f",&ii,&jj,&h);
          info->topo_data.heights[ii+info->topo_data.nx*jj]=h;
        }
    fclose(fp);
  }
  else if( !strcmp( surfaceName, "sine" )
           || !strcmp( surfaceName, "cosine" ) ) {
    Dictionary_Entry_Value*	originList;

    if( !strcmp( surfaceName, "sine" ) )
      {
        *surfaceType = SurfaceAdaptor_SurfaceType_Sine;
        *deformFunc = SurfaceAdaptor_Sine;
      }
    else
      {
        *surfaceType = SurfaceAdaptor_SurfaceType_Cosine;
        *deformFunc = SurfaceAdaptor_Cosine;
      }
    strcpy(temp,surface);
    originList =
      Dictionary_Get( dict, (Dictionary_Entry_Key)strcat(temp,"Origin") );
    if( originList ) {
      unsigned	nDims;
      unsigned	d_i;

      nDims = Dictionary_Entry_Value_GetCount( originList );
      for( d_i = 0; d_i < nDims; d_i++  ) {
        Dictionary_Entry_Value*	val;

        val = Dictionary_Entry_Value_GetElement( originList, d_i );
        info->trig.origin[d_i] = Dictionary_Entry_Value_AsDouble( val );
      }
    }
    else
      memset( info->trig.origin, 0, sizeof(double) * 2 );

    strcpy(temp,surface);
    info->trig.amp =
      Stg_ComponentFactory_GetDouble( cf, name,
                                      (Dictionary_Entry_Key)strcat(temp,"Amplitude"), 1.0  );
    strcpy(temp,surface);
    info->trig.freq =
      Stg_ComponentFactory_GetDouble( cf, name,
                                      (Dictionary_Entry_Key)strcat(temp,"Frequency"), 1.0 );
  } else if( !strcmp( surfaceName, "cylinder" ) ) {
    *surfaceType = SurfaceAdaptor_SurfaceType_Cylinder;
    *deformFunc = SurfaceAdaptor_Cylinder ;
    strcpy(temp,surface);
    info->cylinder.origin[0] =
      Stg_ComponentFactory_GetDouble( cf, name,
                                      (Dictionary_Entry_Key)strcat(temp,"X0"), 0.0  );
    strcpy(temp,surface);
    info->cylinder.origin[1] =
      Stg_ComponentFactory_GetDouble( cf, name,
                                      (Dictionary_Entry_Key)strcat(temp,"Y0"), 0.0  );
    strcpy(temp,surface);
    info->cylinder.r = 
      Stg_ComponentFactory_GetDouble( cf, name,
                                      (Dictionary_Entry_Key)strcat(temp,"Radius"), 0.0  );
    strcpy(temp,surface);
    info->cylinder.minX = 
      Stg_ComponentFactory_GetDouble( cf, name,
                                      (Dictionary_Entry_Key)strcat(temp,"MinX"),
                                      info->cylinder.origin[0] - info->cylinder.r);
    strcpy(temp,surface);
    info->cylinder.maxX = 
      Stg_ComponentFactory_GetDouble( cf, name,
                                      (Dictionary_Entry_Key)strcat(temp,"MaxX"),
                                      info->cylinder.origin[0] + info->cylinder.r);
    strcpy(temp,surface);
    info->cylinder.sign = 
      Stg_ComponentFactory_GetBool( cf, name,
                                      (Dictionary_Entry_Key)strcat(temp,"Sign"), True  );
  }
  else
    Journal_Firewall(!strcmp(surfaceName,""),Journal_Register( Error_Type, name ),
                     "Unknown type of surface for SurfaceAdaptor: %s\n",
                     surfaceName);
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

	/* Extract the cartesian information. */
	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, (Name)"vertexGrid" )  );
	inds = AllocArray( unsigned, Mesh_GetDimSize( mesh ) );

	/* Loop over domain nodes. */
	sync = IGraph_GetDomain( mesh->topo, MT_VERTEX );
	for( n_i = 0; n_i < Sync_GetNumDomains( sync ); n_i++ ) {
		unsigned gNode;
		double percentage, min_height, max_height;
                double topDeform=0.0;
                double bottomDeform=0.0;

		gNode = Sync_DomainToGlobal( sync, n_i );
		Grid_Lift( grid, gNode, inds );

		/* Calculate a height percentage. */
		percentage = (double)(inds[1]) / (double)(grid->sizes[1] - 1);

		/* Deform this node. */
                if(self->topDeformFunc)
                  topDeform = self->topDeformFunc( &(self->top_info), mesh,
                                                   grid->sizes, n_i, inds);
                if(self->bottomDeformFunc)
                  bottomDeform =
                    self->bottomDeformFunc(&(self->bottom_info),mesh,
                                           grid->sizes, n_i, inds);
                min_height=((CartesianGenerator*)self->generator)->crdMin[1] + bottomDeform;
                max_height=((CartesianGenerator*)self->generator)->crdMax[1] + topDeform;

		mesh->verts[n_i][1] = percentage * (max_height - min_height)
                  + min_height;
	}

	/* Free resources. */
	FreeArray( inds );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

double SurfaceAdaptor_Wedge( SurfaceAdaptor_SurfaceInfo *info, Mesh* mesh, 
                             unsigned* globalSize, unsigned vertex,
                             unsigned* vertexInds )
{
  if ( mesh->topo->nDims != 3 )
    {
      if( mesh->verts[vertex][0] >= info->wedge.offs[0] )
        {
          if( mesh->verts[vertex][0] >= info->wedge.endOffs[0] )
            return (info->wedge.endOffs[0] - info->wedge.offs[0])
              * info->wedge.grad[0];
          else
            return (mesh->verts[vertex][0] - info->wedge.offs[0])
              * info->wedge.grad[0];
        }
      else 
        return 0.0;
    }
  else
    {
      if( mesh->verts[vertex][0] >= info->wedge.offs[0] )
        {
          if( mesh->verts[vertex][0] >= info->wedge.endOffs[0] )
            {
              return (info->wedge.endOffs[0] - info->wedge.offs[0])
                * info->wedge.grad[0]
                + (mesh->verts[vertex][2] - info->wedge.offs[1])
                * info->wedge.grad[1];
            }
          else
            {
              return (mesh->verts[vertex][0] - info->wedge.offs[0])
                * info->wedge.grad[0]
                + (mesh->verts[vertex][2] - info->wedge.offs[1])
                * info->wedge.grad[1];
            }
        }
      else 
        return 0.0;
    }
}

double SurfaceAdaptor_Plateau( SurfaceAdaptor_SurfaceInfo *info, Mesh* mesh, 
                               unsigned* globalSize, unsigned vertex,
                               unsigned* vertexInds )
{
  double x_factor, z_factor;
  x_factor =1;
  z_factor=1;
  if( mesh->verts[vertex][0] < info->plateau.x1
      || mesh->verts[vertex][0] > info->plateau.x4)
    {
      x_factor=0;
    }
  else if( mesh->verts[vertex][0] <= info->plateau.x2)
    {
      x_factor=(mesh->verts[vertex][0] - info->plateau.x1)
        /(info->plateau.x2 - info->plateau.x1);
    }
  else if( mesh->verts[vertex][0] <= info->plateau.x3)
    {
      x_factor=1;
    }
  else if( mesh->verts[vertex][0] <= info->plateau.x4)
    {
      x_factor=(info->plateau.x4 - mesh->verts[vertex][0])
        /(info->plateau.x4 - info->plateau.x3);
    }

  if(mesh->topo->nDims==3)
    {
      if( mesh->verts[vertex][2] < info->plateau.z1
          || mesh->verts[vertex][2] > info->plateau.z4)
        {
          z_factor=0;
        }
      else if( mesh->verts[vertex][2] <= info->plateau.z2)
        {
          z_factor=(mesh->verts[vertex][2] - info->plateau.z1)
            /(info->plateau.z2 - info->plateau.z1);
        }
      else if( mesh->verts[vertex][2] <= info->plateau.z3)
        {
          z_factor=1;
        }
      else if( mesh->verts[vertex][2] <= info->plateau.z4)
        {
          z_factor=(info->plateau.z4 - mesh->verts[vertex][2])
            /(info->plateau.z4 - info->plateau.z3);
        }
    }

  return x_factor*z_factor*info->plateau.height;
}

double SurfaceAdaptor_Topo_Data( SurfaceAdaptor_SurfaceInfo *info, Mesh* mesh, 
                                 unsigned* globalSize, unsigned vertex,
                                 unsigned* vertexInds )
{
  int i,k,ip,kp;
  double dx,dz;

  i=floor((mesh->verts[vertex][0] - info->topo_data.minX)
          /info->topo_data.dx + 0.5);
  k=floor((mesh->verts[vertex][2] - info->topo_data.minZ)
          /info->topo_data.dz + 0.5);

  if(i<0 || i>info->topo_data.nx-1
     || k<0 || k>info->topo_data.nz-1)
    {
      printf("Coordinate not covered by the topography file: %g %g\n\tminX: %g\n\tmaxX: %g\n\tminZ: %g\n\tmaxZ: %g\n\tnx: %d\n\tnz: %d\n",
             mesh->verts[vertex][0],
             mesh->verts[vertex][2],
             info->topo_data.minX,
             info->topo_data.maxX,
             info->topo_data.minZ,
             info->topo_data.maxZ,
             info->topo_data.nx,
             info->topo_data.nz);
      abort();
    }

  /* Interpolate the height */
  ip=i+1;
  kp=k+1;
  if(ip>info->topo_data.nx-1)
    ip=i;
  if(kp>info->topo_data.nz-1)
    kp=k;

  dx=(mesh->verts[vertex][0]
      - (i*info->topo_data.dx+info->topo_data.minX))
    /info->topo_data.dx;
  dz=(mesh->verts[vertex][2]
      - (k*info->topo_data.dz+info->topo_data.minZ))
    /info->topo_data.dz;

  return info->topo_data.heights[i+info->topo_data.nx*k]*(1-dx)*(1-dz)
    + info->topo_data.heights[i+info->topo_data.nx*kp]*(1-dx)*dz
    + info->topo_data.heights[ip+info->topo_data.nx*k]*dx*(1-dz)
    + info->topo_data.heights[ip+info->topo_data.nx*kp]*dx*dz;
}

double SurfaceAdaptor_Sine( SurfaceAdaptor_SurfaceInfo *info, Mesh* mesh, 
			    unsigned* globalSize, unsigned vertex, unsigned* vertexInds )
{
	double	dx, dy;
	double	rad;

	dx = mesh->verts[vertex][0] - info->trig.origin[0];
	rad = dx * dx;
	if( mesh->topo->nDims == 3 ) {
		dy = mesh->verts[vertex][1] - info->trig.origin[1];
		rad += dy * dy;
	}
	rad = sqrt( rad );

	return info->trig.amp * sin( info->trig.freq * rad );
}

double SurfaceAdaptor_Cosine( SurfaceAdaptor_SurfaceInfo *info, Mesh* mesh, 
			      unsigned* globalSize, unsigned vertex, unsigned* vertexInds )
{
	double	dx, dz;
	double	rad;

	dx = mesh->verts[vertex][0] - info->trig.origin[0];
	rad = dx * dx;
	if( mesh->topo->nDims == 3 ) {
		dz = mesh->verts[vertex][2] - info->trig.origin[2];
		rad += dz * dz;
	}
	rad = sqrt( rad );

	return info->trig.amp * cos( info->trig.freq * rad );
}

double SurfaceAdaptor_Cylinder( SurfaceAdaptor_SurfaceInfo *info, Mesh* mesh, 
                                unsigned* globalSize, unsigned vertex,
                                unsigned* vertexInds )
{
  double x, x0, minX, maxX, y0, r;

  x = mesh->verts[vertex][0];
  x0= info->cylinder.origin[0];
  minX=info->cylinder.minX;
  maxX=info->cylinder.maxX;

  y0= info->cylinder.origin[1];
  r = info->cylinder.r;
  if(x<minX)
    {
      x=minX;
    }
  else if(x>maxX)
    {
      x=maxX;
    }
  return info->cylinder.sign ? y0+sqrt(r*r-(x-x0)*(x-x0))
    : y0-sqrt(r*r-(x-x0)*(x-x0));
}


