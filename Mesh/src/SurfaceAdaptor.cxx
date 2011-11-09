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
#include <StgDomain/Utils/Utils.h>

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
#include <list>

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
	
	/* Print parent */
	Journal_Printf( stream, "SurfaceAdaptor (ptr): (%p)\n", self );
	_MeshAdaptor_Print( self, stream );
}

void _SurfaceAdaptor_AssignFromXML_Surface(Stg_ComponentFactory* cf,
                                           Name name,
                                           Dictionary* dict,
                                           SurfaceAdaptor_SurfaceInfo *info,
                                           SurfaceAdaptor_DeformFunc **deformFunc,
                                           const std::string &surface);

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
        _SurfaceAdaptor_AssignFromXML_Surface(cf,self->name,
                                              dict,
                                              &(self->top_info),
                                              &(self->topDeformFunc),
                                              "top");
        _SurfaceAdaptor_AssignFromXML_Surface(cf,self->name,
                                              dict,
                                              &(self->bottom_info),
                                              &(self->bottomDeformFunc),
                                              "bottom");
}

void _SurfaceAdaptor_AssignFromXML_Surface(Stg_ComponentFactory* cf,
                                           Name name,
                                           Dictionary* dict,
                                           SurfaceAdaptor_SurfaceInfo *info,
                                           SurfaceAdaptor_DeformFunc **deformFunc,
                                           const std::string &surface)
{
  std::string surfaceName;
  /* What kind of surface do we want? */
  surfaceName = 
    Stg_ComponentFactory_GetString(cf, name,
                                   (surface+"SurfaceType").c_str(),"");
  if(surfaceName.empty() || surfaceName=="equation")
    {
      char *equation=Stg_ComponentFactory_GetString(cf,name,(surface+"Equation").c_str(),"");
      if(strlen(equation)!=0)
        {
          *deformFunc = SurfaceAdaptor_Equation;
          /* This will never get free'd */
          info->equation=StG_Strdup(equation);
        }
    }
  else if(surfaceName=="wedge")
    {
      *deformFunc = SurfaceAdaptor_Wedge ;
      info->wedge.offs[0] =
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"BeginOffset").c_str(),
                                       0.0);
      info->wedge.endOffs[0] =
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"EndOffset").c_str(),
                                       1.0);
      info->wedge.grad[0] = 
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"Gradient").c_str(),0.5);
      /* get the parameters for the z-axis */
      info->wedge.offs[1] = 
        Stg_ComponentFactory_GetDouble(cf,name,
                                       (surface+"BeginOffsetZ").c_str(),0.0);
      info->wedge.endOffs[1] = 
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"EndOffsetZ").c_str(),
                                       1.0);
      info->wedge.grad[1] = 
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"GradientZ").c_str(),
                                       0.5);
                                     
    }
  else if(surfaceName=="plateau")
    {
      *deformFunc = SurfaceAdaptor_Plateau;
      info->plateau.x1 =
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"X1").c_str(), 0.0 );
      info->plateau.x2 = 
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"X2").c_str(), 0.0 );
      info->plateau.x3 = 
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"X3").c_str(), 0.0 );
      info->plateau.x4 = 
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"X4").c_str(), 0.0 );
      info->plateau.z1 = 
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"Z1").c_str(), 0.0 );
      info->plateau.z2 = 
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"Z2").c_str(), 0.0 );
      info->plateau.z3 = 
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"Z3").c_str(), 0.0 );
      info->plateau.z4 = 
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"Z4").c_str(), 0.0 );
      info->plateau.height =
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"Height").c_str(), 0.0 );
    }
  else if(surfaceName=="topo_data")
    {
      FILE *fp;
      char* surfaceFileName;
      int i,j,ii,jj;
      *deformFunc = SurfaceAdaptor_Topo_Data;
      surfaceFileName =
        Stg_ComponentFactory_GetString(cf,name,(surface+"SurfaceName").c_str(),
                                       "ascii_topo" );
      info->topo_data.nx =
        Stg_ComponentFactory_GetInt(cf,name,(surface+"Nx").c_str(),0);
      info->topo_data.nz = 
        Stg_ComponentFactory_GetInt(cf,name,(surface+"Nz").c_str(),1);
      info->topo_data.minX = 
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"MinX").c_str(),0);
      info->topo_data.minZ = 
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"MinZ").c_str(),0);
      info->topo_data.maxX = 
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"MaxX").c_str(),0);
      info->topo_data.maxZ =
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"MaxZ").c_str(),0);
      info->topo_data.dx=
        (info->topo_data.maxX-info->topo_data.minX)/(info->topo_data.nx-1);
      info->topo_data.dz=
        (info->topo_data.maxZ-info->topo_data.minZ)/(info->topo_data.nz-1);
      info->topo_data.heights=
        (double*)malloc(sizeof(double)*info->topo_data.nx*info->topo_data.nz);
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
  else if(surfaceName=="sine" || surfaceName=="cosine")
    {
      Dictionary_Entry_Value*	originList;

      if(surfaceName=="sine")
        {
          *deformFunc = SurfaceAdaptor_Sine;
        }
      else
        {
          *deformFunc = SurfaceAdaptor_Cosine;
        }
      originList=Dictionary_Get(dict,(surface+"Origin").c_str());
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
        memset( info->trig.origin, 0, sizeof(double) * 3 );

      info->trig.amp =
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"Amplitude").c_str(),1.0);
      info->trig.freq =
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"Frequency").c_str(),1.0);
    }
  else if(surfaceName=="cylinder")
    {
      *deformFunc = SurfaceAdaptor_Cylinder ;
      info->cylinder.origin[0] =
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"X0").c_str(),0.0);
      info->cylinder.origin[1] =
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"Y0").c_str(),0.0);
      info->cylinder.r = 
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"Radius").c_str(),0.0);
      info->cylinder.minX = 
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"MinX").c_str(),
                                       info->cylinder.origin[0]-info->cylinder.r);
      info->cylinder.maxX = 
        Stg_ComponentFactory_GetDouble(cf,name,(surface+"MaxX").c_str(),
                                       info->cylinder.origin[0]+info->cylinder.r);
      info->cylinder.sign = 
        Stg_ComponentFactory_GetBool(cf,name,(surface+"Sign").c_str(),True);
    }
  else
    Journal_Firewall(surfaceName.empty(),Journal_Register( Error_Type, name ),
                     "Unknown type of surface for SurfaceAdaptor: %s\n",
                     surfaceName.c_str());
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
	for( n_i = 0; n_i < (unsigned)Sync_GetNumDomains( sync ); n_i++ ) {
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

double SurfaceAdaptor_Equation(SurfaceAdaptor_SurfaceInfo *info, Mesh* mesh, 
                               unsigned* globalSize, unsigned vertex,
                               unsigned* vertexInds )
{
  return Equation_eval(mesh->verts[vertex],(DomainContext*)(mesh->context),
                       info->equation);
}

double SurfaceAdaptor_Wedge(SurfaceAdaptor_SurfaceInfo *info, Mesh* mesh, 
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

  if(mesh->topo->nDims==3)
    {
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
  else
    {
      i=floor((mesh->verts[vertex][0] - info->topo_data.minX)
              /info->topo_data.dx + 0.5);
      if(i<0 || i>info->topo_data.nx-1)
        {
          printf("Coordinate not covered by the topography file: %g\n\tminX: %g\n\tmaxX: %g\n\tnx: %d\n\n",
                 mesh->verts[vertex][0],
                 info->topo_data.minX,
                 info->topo_data.maxX,
                 info->topo_data.nx);
          abort();
        }

      /* Interpolate the height */
      ip=i+1;
      if(ip>info->topo_data.nx-1)
        ip=i;

      dx=(mesh->verts[vertex][0]
          - (i*info->topo_data.dx+info->topo_data.minX))
        /info->topo_data.dx;

      return info->topo_data.heights[i]*(1-dx)
        + info->topo_data.heights[ip]*dx;
    }
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


