/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**	Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
** Contributors:
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
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
** $Id: /cig/src/StgFEM/plugins/Output/PrintFeVariableDiscreteValues/Plugin.c 21 2006-04-07T21:29:57.251207Z walter  $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include "VTKOutput.h"

#include <stdlib.h>
#include <string.h>

const Type Underworld_VTKOutput_Type = "Underworld_VTKOutput";

void _Underworld_VTKOutput_Construct( void* component, Stg_ComponentFactory* cf, void *data ) {
	UnderworldContext* context;

	context = (UnderworldContext*)Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data );

	ContextEP_Append( context, AbstractContext_EP_Dump,
                          VTKOutput );
}

void* _Underworld_VTKOutput_DefaultNew( Name name ) {
	return Codelet_New(
			Underworld_VTKOutput_Type,
			_Underworld_VTKOutput_DefaultNew,
			_Underworld_VTKOutput_Construct,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}

Index Underworld_VTKOutput_Register( PluginsManager* pluginsManager ) {
	Journal_DPrintf( StgFEM_Debug, "In: %s( void* )\n", __func__ );
	
	return PluginsManager_Submit( pluginsManager, Underworld_VTKOutput_Type, "0", _Underworld_VTKOutput_DefaultNew );
}

void VTKOutput_particles(IntegrationPointsSwarm*  picswarm, 
                         double defaultDiffusivity,
                         int stepping,
                         char *outputPath, int timeStep, int dim, int myRank,
                         int nprocs);
void VTKOutput_fields(void *context, int myRank, int nprocs);

void VTKOutput( void* _context ) {
	UnderworldContext*	context = (UnderworldContext*)_context;
	Dictionary*             dictionary         = context->dictionary;

        int myRank, nprocs;
        MPI_Comm comm;

	if(context->picIntegrationPoints)
	    comm = Comm_GetMPIComm( Mesh_GetCommTopology( context->picIntegrationPoints->mesh, MT_VERTEX ) );
	else
	    comm = MPI_COMM_WORLD;

	MPI_Comm_rank( comm, (int*)&myRank );
        MPI_Comm_size( comm, (int*)&nprocs );

        /* Only dump if at the right time step. */
        if(context->timeStep % context->dumpEvery != 0)
          return;
	
        /* Write the particles and then all of the fields. */

	if(context->picIntegrationPoints) {
	    VTKOutput_particles(context->picIntegrationPoints,
				Dictionary_GetDouble_WithDefault
				(dictionary,"defaultDiffusivity",1.0),
				Dictionary_GetInt_WithDefault
				(dictionary,"particleStepping",1),
				context->outputPath, context->timeStep,
				context->dim,myRank,nprocs);
	}
        VTKOutput_fields(context,myRank,nprocs);
}

void VTKOutput_particles(IntegrationPointsSwarm*  picswarm,
                         double defaultDiffusivity,
                         int stepping, char *outputPath,
                         int timeStep, int dim, int myRank, int nprocs) {
  double *coord;
  int iteration, i;
  Particle_Index          num_particles = picswarm->particleLocalCount;
  Particle_Index          lParticle_I;
  
  RheologyMaterial*       material;
  MaterialPointsSwarm*    materialSwarm;
  MaterialPoint*          materialparticle;
 
  Rheology_Register*      rheology_register;
  Rheology_Index      rheology_I; 
  Rheology_Index      rheologyCount;
  YieldRheology*      rheology; 
 
  FILE *fp, *pfp;
  Name filename;
  
  /* Open the processor specific output file */
  Stg_asprintf(&filename,"%s/particles.%d.%05d.vtu",outputPath,myRank,timeStep);
  fp=fopen(filename,"w");
  Memory_Free( filename );

  /* Open the parallel control file */
  if(myRank==0)
    {
      Stg_asprintf( &filename, "%s/particles.%05d.pvtu", outputPath, timeStep);
      pfp=fopen(filename,"w");
      Memory_Free( filename );
    }
      
  /* Write a header */
  fprintf(fp,"<?xml version=\"1.0\"?>\n");
  fprintf(fp,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
  fprintf(fp,"  <UnstructuredGrid>\n");
  fprintf(fp,"    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"1\">\n",
          (num_particles-1)/stepping+1);
  fprintf(fp,"      <Points>\n");
  fprintf(fp,"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");

  if(myRank==0)
    {
      fprintf(pfp,"<?xml version=\"1.0\"?>\n\
<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n\
  <PUnstructuredGrid GhostLevel=\"0\">\n\
    <PPoints>\n\
      <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"/>\n\
    </PPoints>\n\
    <PCellData></PCellData>\n");
    }
  
  /* We need many iterations, because the values are written
     separately from the coordinates. */
  for(iteration=0; iteration<9; ++iteration)
    {
      /* Loop over all of the particles */
      for ( lParticle_I = 0 ; lParticle_I < num_particles ;
            lParticle_I+=stepping ){
        double yielding, viscosity, density, alpha, diffusivity;
        Material_Index material_index;
        SymmetricTensor stress;
        BuoyancyForceTerm_MaterialExt*   materialExt;
        Material *extension_info;
        double normal[3];
        
        IntegrationPoint* integrationparticle = (IntegrationPoint*)Swarm_ParticleAt( picswarm, lParticle_I );
        
        material = (RheologyMaterial*) IntegrationPointsSwarm_GetMaterialOn( picswarm, integrationparticle );
        materialparticle = OneToOneMapper_GetMaterialPoint( picswarm->mapper, integrationparticle, &materialSwarm );
        
        density=Dictionary_GetDouble_WithDefault( material->dictionary,
                                                  "density", 0.0 );
        alpha=Dictionary_GetDouble_WithDefault( material->dictionary,
                                                "alpha", 0.0 );
	material_index=material->index;
        diffusivity=Dictionary_GetDouble_WithDefault
          ( material->dictionary, "diffusivity", defaultDiffusivity );
        rheology_register=(Rheology_Register*)material->rheology_Register;

        if( rheology_register && strcmp( material->name, Material_Type ) )
           rheologyCount = Rheology_Register_GetCount( rheology_register );
        else
           rheologyCount = 0;
        
        coord = materialparticle->coord;
        yielding=0;
        viscosity=0;
        SymmetricTensor_Zero(stress);
        
        /* First print out only the coordinates. */
        if(iteration==0)
          {
            if (dim == 2) {
              fprintf(fp,"%lf %lf 0.0 ",(double)coord[0],
                      (double)coord[1]);
            } else {
              fprintf(fp,"%lf %lf %lf ", (double)coord[0],
                      (double)coord[1], (double)coord[2]);
            }
          }
        else
          {
            /* Loop over all of the rheologies for a particle. */
             memset( normal, 0, 3 * sizeof(double) );
            for( rheology_I = 0; rheology_I < rheologyCount ; rheology_I++ ) { 
              rheology = (YieldRheology*)Rheology_Register_GetByIndex( rheology_register, rheology_I ); 

              /* Get yielding information */
/*
              if(!strcmp(rheology->type, "DruckerPrager") || 
                 !strcmp(rheology->type, "VonMises") ||
                 !strcmp(rheology->type, "FaultingMoresiMuhlhaus2006") ||
                 !strcmp(rheology->type, "MohrCoulomb"))
                {
                  yielding=StrainWeakening_CalcRatio(rheology->strainWeakening,
                                                     materialparticle);
                }
*/
              /* Get viscosity */
              if(!strcmp(rheology->type,"StoreVisc"))
                {
                  StoreVisc* self = (StoreVisc*) rheology;
                  StoreVisc_ParticleExt* particleExt;
                  particleExt=
                    ExtensionManager_Get( materialSwarm->particleExtensionMgr, materialparticle, self->particleExtHandle );
                  viscosity=particleExt->effVisc;
                }
              /* Get stress */
              if(!strcmp(rheology->type,"StoreStress"))
                {
                  StoreStress* self = (StoreStress*) rheology;
                  StoreStress_ParticleExt* particleExt;
                  particleExt=
                    ExtensionManager_Get( materialSwarm->particleExtensionMgr, materialparticle, self->particleExtHandle );
                  stress[0]=particleExt->stress[0];
                  stress[1]=particleExt->stress[1];
                  stress[2]=particleExt->stress[2];
                  stress[3]=particleExt->stress[3];
                  if(dim==3)
                    {
                      stress[4]=particleExt->stress[4];
                      stress[5]=particleExt->stress[5];
                    }
                }
              if(!strcmp(rheology->type, "FaultingMoresiMuhlhaus2006")) {
                 Director* director = ((FaultingMoresiMuhlhaus2006*)rheology)->director;
                 Director_GetNormal( director, materialparticle, normal );
              }
            }
            switch(iteration)
              {
              case 1:
                fprintf(fp,"%lf ",viscosity);
                break;
              case 2:
                fprintf(fp,"%lf ",yielding);
                break;
              case 3:
                if(dim==2)
                  {
                    fprintf(fp,"%lf %lf 0.0 %lf %lf 0.0 0.0 0.0 0.0 ",
                            stress[0],stress[1],
                            stress[1],stress[2]);
                  }
                else
                  {
                    fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf ",
                            stress[0],stress[1],stress[2],
                            stress[1],stress[3],stress[4],
                            stress[2],stress[4],stress[5]);
                  }
                break;
              case 4:
                fprintf(fp,"%lf ",density);
                break;
              case 5:
                fprintf(fp,"%d ",material_index);
                break;
              case 6:
                fprintf(fp,"%lf ",alpha);
                break;
              case 7:
                fprintf(fp,"%lf ",diffusivity);
                break;
                 case 8:
                    fprintf(fp, "%lf %lf %lf ", normal[0], normal[1], normal[2]);
              }
          }
      }
      switch(iteration)
        {
        case 0:
          fprintf(fp,"\n        </DataArray>\n");
          fprintf(fp,"      </Points>\n");
          fprintf(fp,"      <PointData Scalars=\"Viscosity\" Tensors=\"Stress\">\n");
          fprintf(fp,"        <DataArray type=\"Float64\" Name=\"Viscosity\" format=\"ascii\">\n");
          if(myRank==0)
            {
              fprintf(pfp,"      <PPointData Scalars=\"Viscosity\" Tensors=\"Stress\">\n");
              fprintf(pfp,"        <PDataArray type=\"Float64\" Name=\"Viscosity\" format=\"ascii\"/>\n");
            }
          break;
        case 1:
          fprintf(fp,"\n        </DataArray>\n");
          fprintf(fp,"        <DataArray type=\"Float64\" Name=\"Yielding_fraction\" format=\"ascii\">\n");
          if(myRank==0)
            fprintf(pfp,"        <PDataArray type=\"Float64\" Name=\"Yielding_fraction\" format=\"ascii\"/>\n");
          break;
        case 2:
          fprintf(fp,"\n        </DataArray>\n");
          fprintf(fp,"        <DataArray type=\"Float64\" Name=\"Stress\" format=\"ascii\" NumberOfComponents=\"9\">\n");
          if(myRank==0)
            fprintf(pfp,"        <PDataArray type=\"Float64\" Name=\"Stress\" format=\"ascii\" NumberOfComponents=\"9\"/>\n");
          break;
        case 3:
          fprintf(fp,"\n        </DataArray>\n");
          fprintf(fp,"        <DataArray type=\"Float64\" Name=\"Density\" format=\"ascii\">\n");
          if(myRank==0)
            fprintf(pfp,"        <PDataArray type=\"Float64\" Name=\"Density\" format=\"ascii\"/>\n");
          break;
        case 4:
          fprintf(fp,"\n        </DataArray>\n");
          fprintf(fp,"        <DataArray type=\"Int32\" Name=\"Material_Index\" format=\"ascii\">\n");
          if(myRank==0)
            fprintf(pfp,"        <PDataArray type=\"Int32\" Name=\"Material_Index\" format=\"ascii\"/>\n");
          break;
        case 5:
          fprintf(fp,"\n        </DataArray>\n");
          fprintf(fp,"        <DataArray type=\"Float64\" Name=\"Alpha\" format=\"ascii\">\n");
          if(myRank==0)
            fprintf(pfp,"        <PDataArray type=\"Float64\" Name=\"Alpha\" format=\"ascii\"/>\n");
          break;
        case 6:
          fprintf(fp,"\n        </DataArray>\n");
          fprintf(fp,"        <DataArray type=\"Float64\" Name=\"Thermal_Diffusivity\" format=\"ascii\">\n");
          if(myRank==0)
            fprintf(pfp,"        <PDataArray type=\"Float64\" Name=\"Thermal_Diffusivity\" format=\"ascii\"/>\n");
          break;
        case 7:
          fprintf(fp,"\n        </DataArray>\n");
          fprintf(fp,"        <DataArray type=\"Float64\" Name=\"Director_Normals\" format=\"ascii\" NumberOfComponents=\"3\">\n");
          if(myRank==0)
            fprintf(pfp,"        <PDataArray type=\"Float64\" Name=\"Director_Normals\" format=\"ascii\" NumberOfComponents=\"3\"/>\n");
          break;
        case 8:
          fprintf(fp,"\n        </DataArray>\n");
          fprintf(fp,"      </PointData>\n");
          fprintf(fp,"      <CellData>\n");
          fprintf(fp,"      </CellData>\n");
          fprintf(fp,"      <Cells>\n");
          fprintf(fp,"        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
          for(i=0;i<(num_particles-1)/stepping+1;++i)
            fprintf(fp,"%d ",i);
          fprintf(fp,"\n        </DataArray>\n");
          fprintf(fp,"        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
          fprintf(fp,"          %d\n",(num_particles-1)/stepping+1);
          fprintf(fp,"        </DataArray>\n");
          fprintf(fp,"        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
          fprintf(fp,"          2\n");
          fprintf(fp,"        </DataArray>\n");
          fprintf(fp,"      </Cells>\n");
          fprintf(fp,"    </Piece>\n");
          fprintf(fp,"  </UnstructuredGrid>\n");
          fprintf(fp,"</VTKFile>\n");
        }
    }
  fclose(fp);
  if(myRank==0)
    {
      fprintf(pfp,"        </PPointData>\n");
      for(i=0;i<nprocs;++i)
        fprintf(pfp,"    <Piece Source=\"particles.%d.%05d.vtu\"/>\n",
                i,timeStep);
      fprintf(pfp,"  </PUnstructuredGrid>\n\
</VTKFile>\n");
      fclose(pfp);
    }
}

/* Print out the coordinates of the mesh. */

void VTKOutput_print_coords(FILE *fp, FeMesh *feMesh, Grid *grid, int nDims,
                            int lower[3], int upper[3]) {
  IJK ijk;

  fprintf(fp,"      <Points>\n");
  fprintf(fp,"        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  
  for(ijk[2]=lower[2];ijk[2]<upper[2];++ijk[2])
    for(ijk[1]=lower[1];ijk[1]<upper[1];++ijk[1])
      for(ijk[0]=lower[0];ijk[0]<upper[0];++ijk[0])
        {
          double *coord;
          unsigned local;
          Mesh_GlobalToDomain(feMesh,MT_VERTEX,Grid_Project(grid,ijk),&local);
          coord=Mesh_GetVertex(feMesh,local);
          if(nDims==2)
            fprintf(fp, "%.15g %.15g 0\n", coord[0], coord[1]);
          else
            fprintf(fp, "%.15g %.15g %.15g\n", coord[0], coord[1], coord[2]);
        }
  fprintf(fp,"        </DataArray>\n      </Points>\n");
}


/* Pressure is stored on cell centers, while everything else is stored
   on cell vertices.  So we have to make two different files, one for
   pressure, and one for everything else. */
void VTKOutput_fields(void *context, int myRank, int nprocs) {
  
  FiniteElementContext*     self = (FiniteElementContext*) context;
  Index var_I;
  int header_printed=0;
  int nDims, pn, i;
  int lower[3], upper[3], p_lower[3], p_upper[3];

  Name field_filename, pressure_filename;
  FILE *fp, *field_fp, *pressure_fp, *pfp, *pfield_fp, *ppressure_fp;

  /* We need to save the grids to be used later when printing out the
     pressure coordinates and to map between 1D and 3D indices for the
     values. */
  Grid *elGrid, *vertGrid;

  /* Open the file */

  Stg_asprintf( &pressure_filename, "%s/pressure.%d.%05d.vts",
                self->outputPath,myRank,
                self->timeStep);
  Stg_asprintf( &field_filename, "%s/fields.%d.%05d.vts", self->outputPath,
                myRank, self->timeStep);

  field_fp=fopen(field_filename,"w");
  pressure_fp=fopen(pressure_filename,"w");

  Memory_Free( pressure_filename );
  Memory_Free( field_filename );

  /* Print out the parallel control files if rank==0 */
  if(myRank==0)
    {
      Stg_asprintf( &pressure_filename, "%s/pressure.%05d.pvts",
                    self->outputPath,
                    self->timeStep);
      Stg_asprintf( &field_filename, "%s/fields.%05d.pvts", self->outputPath,
                    self->timeStep);
      pfield_fp=fopen(field_filename,"w");
      ppressure_fp=fopen(pressure_filename,"w");
      Memory_Free( pressure_filename );
      Memory_Free( field_filename );
    }

  /* First, output the coordinates.  We have to do a huge song and
     dance just to get the extents of the mesh. */

  for ( var_I = 0; var_I < self->fieldVariable_Register->objects->count;
        var_I++ ) {
    FieldVariable* fieldVar;
    fieldVar = FieldVariable_Register_GetByIndex( self->fieldVariable_Register,
                                                  var_I );
    if (Stg_Class_IsInstance( fieldVar, FeVariable_Type )) {
      FeVariable* feVar;
      Node_LocalIndex    lNode_I;
      Dof_Index          dofAtEachNodeCount;
      int *low, *up;
      Grid *grid;
      IJK ijk;

      feVar=(FeVariable*)fieldVar;

      if(!strcmp(feVar->name,"HeightField"))
        continue;

      if(!strcmp(feVar->name,"VelocityField") && self->timeStep==0)
        FeVariable_SyncShadowValues(feVar);
      if(!header_printed)
        {
          Mesh *mesh;
          CartesianGenerator *gen;

          mesh=(Mesh*)(feVar->feMesh);
          gen=((CartesianGenerator *)(mesh->generator));
          /* If we got the surface adaptor instead of the cartesian
             mesh generator, go to the mesh generator.  */
          if(!strcmp(gen->type,"SurfaceAdaptor") || !strcmp(gen->type,"FieldVariableSurfaceAdaptor"))
            gen=(CartesianGenerator *)((SurfaceAdaptor *)(gen))->generator;

          elGrid=gen->elGrid;
          vertGrid=gen->vertGrid;
          nDims=gen->nDims;

          p_lower[2]=lower[2]=0;
          p_upper[2]=upper[2]=1;
          for(i=0;i<nDims;++i)
            {
              p_lower[i]=gen->origin[i];
              p_upper[i]=p_lower[i]+gen->range[i];

              /* The grid is split up by elements, so for the pressure
                 grid, we simply add the ghost zones.  Ghost zones are
                 only added at the top, not the bottom. */
              if(p_upper[i]!=gen->elGrid->sizes[i])
                ++p_upper[i];

              /* Then we get the vertices by adding appropriate
                 offsets to the element grid. */
              upper[i]=p_upper[i]+1;
              lower[i]=p_lower[i];
              if(p_lower[i]!=0)
                --lower[i];
            }

          /* Now that we have the extents, write the header. */

          fprintf(field_fp,"<?xml version=\"1.0\"?>\n\
<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n\
  <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n\
    <Piece Extent=\"%d %d %d %d %d %d\">\n\
      <CellData></CellData>\n",
                  lower[0],upper[0]-1,lower[1],upper[1]-1,lower[2],upper[2]-1,
                  lower[0],upper[0]-1,lower[1],upper[1]-1,lower[2],upper[2]-1);
          fprintf(pressure_fp,"<?xml version=\"1.0\"?>\n\
<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n\
  <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n\
    <Piece Extent=\"%d %d %d %d %d %d\">\n\
      <CellData></CellData>\n",
                  p_lower[0],p_upper[0]-1,p_lower[1],p_upper[1]-1,
                  p_lower[2],p_upper[2]-1,
                  p_lower[0],p_upper[0]-1,p_lower[1],p_upper[1]-1,
                  p_lower[2],p_upper[2]-1);

          /* Write the coordinates for the fields, but not the pressure */
          VTKOutput_print_coords(field_fp, feVar->feMesh, gen->vertGrid, nDims,
                                 lower, upper);
          fprintf(field_fp,"      <PointData Scalars=\"StrainRateInvariantField\" Vectors=\"VelocityField\" Tensors=\"Stress\">\n");
          
          /* Write out parallel control file headers. */
          if(myRank==0)
            {
              int global[3], p_global[3];
              global[0]=gen->vertGrid->sizes[0];
              global[1]=gen->vertGrid->sizes[1];
              if(nDims==3)
                {
                  global[2]=gen->vertGrid->sizes[2];
                  p_global[2]=global[2]-1;
                }
              else
                {
                  global[2]=p_global[2]=1;
                }
              fprintf(pfield_fp,"<?xml version=\"1.0\"?>\n\
<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n\
  <PStructuredGrid GhostLevel=\"0\" WholeExtent=\"0 %d 0 %d 0 %d\">\n\
    <PCellData></PCellData>\n\
    <PPoints>\n\
      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n\
    </PPoints>\n\
    <PPointData Scalars=\"StrainRateInvariantField\" Vectors=\"VelocityField\" Tensors=\"Stress\">\n",global[0]-1,global[1]-1,global[2]-1);
              fprintf(ppressure_fp,"<?xml version=\"1.0\"?>\n\
<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n \
  <PStructuredGrid GhostLevel=\"0\" WholeExtent=\"0 %d 0 %d 0 %d\">\n\
    <PCellData></PCellData>\n\
    <PPoints>\n\
      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n\
    </PPoints>\n\
    <PPointData Scalars=\"PressureField\">\n",
                      global[0]-2,global[1]-2,p_global[2]-1);
            }
          header_printed=1;
        }

      /* Write the coordinates for the pressure, and set the
         appropriate file pointer to output for this variable. */
      if(!strcmp(feVar->name,"PressureField") && !strcmp(feVar->feMesh->name, "constantMesh")){
        VTKOutput_print_coords(pressure_fp, feVar->feMesh, elGrid, nDims,
                               p_lower, p_upper);
        fprintf(pressure_fp,"      <PointData Scalars=\"PressureField\">\n");
        fp=pressure_fp;
        pfp=ppressure_fp;
        low=p_lower;
        up=p_upper;
        grid=elGrid;
      } else {
        fp=field_fp;
        pfp=pfield_fp;
        low=lower;
        up=upper;
        grid=vertGrid;
      }
      
      /* Finally, output the fields.  For now, just output every field */

      dofAtEachNodeCount = feVar->fieldComponentCount;

      switch(dofAtEachNodeCount*nDims)
        {
          /* Scalars */
        case 2:
        case 3:
          fprintf(fp,"        <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n",feVar->name);
          if(myRank==0)
            fprintf(pfp,"      <PDataArray type=\"Float64\" Name=\"%s\" format=\"ascii\"/>\n",feVar->name);
          break;
          /* Vectors */
        case 4:
        case 9:
          fprintf(fp,"        <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\" NumberOfComponents=\"3\">\n",feVar->name);
          if(myRank==0)
            fprintf(pfp,"      <PDataArray type=\"Float64\" Name=\"%s\" format=\"ascii\" NumberOfComponents=\"3\"/>\n",feVar->name);
          break;
          /* Rank 2 Tensors */
        case 6:
        case 8:
        case 18:
        case 27:
          fprintf(fp,"        <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\" NumberOfComponents=\"9\">\n",feVar->name);
          if(myRank==0)
            fprintf(pfp,"      <PDataArray type=\"Float64\" Name=\"%s\" format=\"ascii\" NumberOfComponents=\"9\"/>\n",feVar->name);
          break;
          /* Unknown */
        default:
          fprintf(stderr,
                  "Bad number of degrees of freedom for variable %s: %d\n",
                 feVar->name, dofAtEachNodeCount);
          abort();
        }

      for(ijk[2]=low[2];ijk[2]<up[2];++ijk[2])
        for(ijk[1]=low[1];ijk[1]<up[1];++ijk[1])
          for(ijk[0]=low[0];ijk[0]<up[0];++ijk[0])
            {
              double variableValues[MAX_FIELD_COMPONENTS];	
              Dof_Index          dof_I;
              unsigned local;
              Mesh_GlobalToDomain(feVar->feMesh,MT_VERTEX,
                                  Grid_Project(grid,ijk),&local);
              FeVariable_GetValueAtNode(feVar,local,variableValues);
                                        
              switch(dofAtEachNodeCount*nDims)
                {
                  /* If writing scalars or 3D objects, then just write
                     the values.  Otherwise, we have to fill in the 3D
                     components with zeros. */
                case 2:
                case 3:
                case 9:
                case 27:
                  for ( dof_I = 0; dof_I < dofAtEachNodeCount; dof_I++ ) {
                    fprintf(fp, "%.15g ", variableValues[dof_I] );
                  }
                  break;
                case 4:
                  fprintf(fp, "%.15g %.15g 0", variableValues[0],
                          variableValues[1] );
                  break;
                case 6:
                  /* Ordering is xx, yy, xy */
                  fprintf(fp, "%.15g %.15g 0 %.15g %.15g 0 0 0 0 ",
                          variableValues[0], variableValues[2],
                          variableValues[2], variableValues[1]);
                  break;
                case 8:
                  fprintf(fp, "%.15g %.15g 0 %.15g %.15g 0 0 0 0 ",
                          variableValues[0], variableValues[1],
                          variableValues[2], variableValues[3]);
                  break;
                case 18:
                  /* Ordering is xx, yy, zz, xy, xz, yz */
                  fprintf(fp, "%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g ",
                          variableValues[0], variableValues[3],
                          variableValues[4], variableValues[3],
                          variableValues[1], variableValues[5],
                          variableValues[4], variableValues[5],
                          variableValues[2]);
                  break;
                }
              fprintf(fp, "\n" );
            }
      fprintf(fp,"        </DataArray>\n");
    }
  }
  fprintf(pressure_fp,"      </PointData>\n    </Piece>\n\
  </StructuredGrid>\n\
</VTKFile>\n");
  fprintf(field_fp,"      </PointData>\n    </Piece>\n\
  </StructuredGrid>\n\
</VTKFile>\n");
  fclose(pressure_fp);
  fclose(field_fp);

  if(myRank==0)
    {
      fprintf(ppressure_fp,"      </PPointData>\n");
      fprintf(pfield_fp,"      </PPointData>\n");
      for(i=0;i<nprocs;++i)
        {
          fprintf(pfield_fp,"    <Piece Extent=\"%d %d %d %d %d %d\"\n\
             Source=\"fields.%d.%05d.vts\"/>\n",lower[0],upper[0]-1,
                  lower[1],upper[1]-1,lower[2],upper[2]-1,
                  i,self->timeStep);
          fprintf(ppressure_fp,"    <Piece Extent=\"%d %d %d %d %d %d\"\n\
             Source=\"pressure.%d.%05d.vts\"/>\n",p_lower[0],p_upper[0]-1,
                  p_lower[1],p_upper[1]-1,p_lower[2],p_upper[2]-1,
                  i,self->timeStep);
        }
      fprintf(ppressure_fp,"  </PStructuredGrid>\n\
</VTKFile>\n");
      fprintf(pfield_fp,"  </PStructuredGrid>\n\
</VTKFile>\n");
      fclose(ppressure_fp);
      fclose(pfield_fp);
     }

}
     
