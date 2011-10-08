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
#include <PICellerator/Utils/HydrostaticTerm.h>

#include <stdlib.h>
#include <string.h>
#include <string>

const Type Underworld_VTKOutput_Type = "Underworld_VTKOutput";

void _Underworld_VTKOutput_AssignFromXML( void* component, Stg_ComponentFactory* cf, void *data ) {
	UnderworldContext* context;

	context = (UnderworldContext*)Stg_ComponentFactory_ConstructByName( cf, (Name)"context", UnderworldContext, True, data  );

        /* Run after the solve, but before time integration
           (e.g. advection).  This gives us data that is always valid.
           It also simplifies intepretation because we do not have to
           separate mesh deformation from what the solver gave us. */

	ContextEP_Append( context, AbstractContext_EP_Solve,
                          VTKOutput );
}

void* _Underworld_VTKOutput_DefaultNew( Name name ) {
	return Codelet_New(
			Underworld_VTKOutput_Type,
			_Underworld_VTKOutput_DefaultNew,
			_Underworld_VTKOutput_AssignFromXML,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}

Index Underworld_VTKOutput_Register( PluginsManager* pluginsManager ) {
	Journal_DPrintf( StgFEM_Debug, "In: %s( void* )\n", __func__ );
	
	return PluginsManager_Submit( pluginsManager, Underworld_VTKOutput_Type, (Name)"0", _Underworld_VTKOutput_DefaultNew  );
}

void VTKOutput_particles(Swarm* swarm, 
                         double defaultDiffusivity,
                         int stepping,
                         char *outputPath, const int timeStep,
                         int dim, int myRank, int nprocs);
void VTKOutput_fields(void *context, int myRank, int nprocs,
                      const int timeStep);

void VTKOutput( void* _context ) {
  UnderworldContext* context = (UnderworldContext*)_context;
  Dictionary* dictionary = context->dictionary;

  /* Only dump if at the right time step.  We use timeStep-1, because
     we are outputing after a solve, but before advection.  So
     timeStep-1 makes more sense in terms of when the simulation looks
     like this. */
  if(context->timeStep % context->dumpEvery != 0)
    return;
	

  /* Write the particles and then all of the fields. */
  
  if(Dictionary_GetBool_WithDefault(dictionary,"VTKOutput_Particles",
                                    True))
    {
      Name swarmName;
      Swarm* swarm;
      Dictionary_Entry_Value* swarm_list;
      int swarm_list_size, s;
      swarm_list=Dictionary_Get(dictionary, "VTKOutput_SwarmList" );

      if(!swarm_list)
        {
          swarm_list=Dictionary_Entry_Value_NewList();
          Dictionary_Entry_Value_AddElement(swarm_list,Dictionary_Entry_Value_FromString("picIntegrationPoints"));
        }
      swarm_list_size=Dictionary_Entry_Value_GetCount(swarm_list);
      for(s=0;s<swarm_list_size;++s)
        {
          swarmName=Dictionary_Entry_Value_AsString(Dictionary_Entry_Value_GetElement(swarm_list,s));
          swarm=(Swarm*)LiveComponentRegister_Get(context->CF->LCRegister,
                                                  swarmName);

          VTKOutput_particles(swarm,
                              Dictionary_GetDouble_WithDefault
                              (dictionary,"defaultDiffusivity",1.0),
                              Dictionary_GetInt_WithDefault
                              (dictionary,"particleStepping",1),
                              context->outputPath, context->timeStep,
                              context->dim,context->rank,context->nproc);
        }
    }
  if(Dictionary_GetBool_WithDefault(dictionary,"VTKOutput_Fields",
                                    True))
    VTKOutput_fields(context,context->rank,context->nproc,
                     context->timeStep);
}

void VTKOutput_particles(Swarm* swarm,
                         double defaultDiffusivity,
                         int stepping, char *outputPath,
                         const int timeStep, int dim, int myRank, int nprocs) {
  double *coord;
  int iteration;
  Particle_Index          num_particles = swarm->particleLocalCount;
  Particle_Index          lParticle_I;
  
  RheologyMaterial*       material;
  MaterialPointsSwarm*    materialSwarm;
  MaterialPoint*          materialparticle;
 
  Rheology_Register*      rheology_register;
  Rheology_Index      rheology_I; 
  Rheology_Index      rheologyCount;
  YieldRheology*      rheology; 
 
  FILE *fp, *pfp;
  char* filename;
  
  /* Open the processor specific output file */
  Stg_asprintf(&filename,"%s/%s.%d.%05d.vtu",outputPath,swarm->name,myRank,
               timeStep);
  fp=fopen(filename,"w");
  Memory_Free( filename );

  if(fp==NULL)
    {
      fprintf(stderr,"WARNING: Can not open particle file for rank %d step %d\n",myRank,
              timeStep);
      return;
    }
  /* Open the parallel control file */
  if(myRank==0)
    {
      Stg_asprintf( &filename, "%s/%s.%05d.pvtu", outputPath, swarm->name, timeStep);
      pfp=fopen(filename,"w");
      Memory_Free( filename );

      if(pfp==NULL)
        {
          fprintf(stderr,"FATAL ERROR: Can not open master particle file for step %d",timeStep);
          abort();
        }
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
  for(iteration=0; iteration<10; ++iteration)
    {
      /* Loop over all of the particles */
      for ( lParticle_I = 0 ; lParticle_I < num_particles ;
            lParticle_I+=stepping ){
        double postFailureStrain, viscosity, density, alpha, diffusivity;
        int currently_yielding;
        Material_Index material_index;
        SymmetricTensor stress;
        XYZ normal;
        
        if(Stg_Class_IsInstance(swarm,IntegrationPointsSwarm_Type))
          {
            IntegrationPoint* integrationparticle = (IntegrationPoint*)Swarm_ParticleAt( ((IntegrationPointsSwarm*)swarm), lParticle_I );
        
            material = (RheologyMaterial*) IntegrationPointsSwarm_GetMaterialOn( ((IntegrationPointsSwarm*)swarm), integrationparticle );
            materialparticle = OneToOneMapper_GetMaterialPoint( ((IntegrationPointsSwarm*)swarm)->mapper, integrationparticle, &materialSwarm );
        
            density=Dictionary_GetDouble_WithDefault( material->dictionary, (Dictionary_Entry_Key)"density", 0.0  );
            alpha=Dictionary_GetDouble_WithDefault( material->dictionary, (Dictionary_Entry_Key)"alpha", 0.0  );
            diffusivity=Dictionary_GetDouble_WithDefault( material->dictionary, (Dictionary_Entry_Key)"diffusivity", defaultDiffusivity );

            material_index=material->index;
            rheology_register=(Rheology_Register* )material->rheology_Register;
        
            if( rheology_register && strcmp( material->name, Material_Type ) )
              rheologyCount = Rheology_Register_GetCount( rheology_register );
            else
              rheologyCount = 0;
            coord = materialparticle->coord;
          }
        else
          {
            GlobalParticle *particle;
            particle=(GlobalParticle *)Swarm_ParticleAt(swarm,lParticle_I);
            coord = particle->coord;
            material_index=0;
            density=alpha=0;
            diffusivity=defaultDiffusivity;
            rheologyCount=0;
          }
        postFailureStrain=0;
        viscosity=0;
        currently_yielding=0;
        SymmetricTensor_Zero(stress);
        normal[0]=normal[1]=normal[2]=0;
        
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
                  
            for( rheology_I = 0; rheology_I < rheologyCount ; rheology_I++ ) { 
              rheology = (YieldRheology*)Rheology_Register_GetByIndex( rheology_register, rheology_I ); 
                
              /* Get postFailureStrain and current_yielding information */
              if(Stg_Class_IsInstance(rheology,YieldRheology_Type)
                 && rheology->strainWeakening)
                {
                  postFailureStrain=
                    StrainWeakening_GetPostFailureWeakening
                    (rheology->strainWeakening,materialparticle);
                  currently_yielding=YieldRheology_GetParticleFlag
                    (rheology,materialSwarm,materialparticle);
                }
              /* Get viscosity */
              if(Stg_Class_IsInstance(rheology,StoreVisc_Type))
                {
                  StoreVisc* self = (StoreVisc*) rheology;
                  StoreVisc_ParticleExt* particleExt;
                  particleExt=
                    (StoreVisc_ParticleExt*)ExtensionManager_Get( materialSwarm->particleExtensionMgr, materialparticle, self->particleExtHandle );
                  viscosity=particleExt->effVisc;
                }
              /* Get stress */
              if(Stg_Class_IsInstance(rheology,StoreStress_Type))
                {
                  StoreStress* self = (StoreStress*) rheology;
                  StoreStress_ParticleExt* particleExt;
                  particleExt=
                    (StoreStress_ParticleExt*)ExtensionManager_Get( materialSwarm->particleExtensionMgr, materialparticle, self->particleExtHandle );
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

              if(Stg_Class_IsInstance(rheology,Anisotropic_Type))
                {
                  Anisotropic* self = (Anisotropic*) rheology;
                  XYZ normal;
                  Director_GetNormal(self->director,materialparticle,normal);
                }


            }
            switch(iteration)
              {
              case 1:
                fprintf(fp,"%lf ",viscosity);
                break;
              case 2:
                fprintf(fp,"%lf ",postFailureStrain);
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
                fprintf(fp,"%d ",currently_yielding);
                break;
              case 9:
                if(dim==2)
                  {
                    fprintf(fp,"%lf %lf 0.0 ",normal[0],normal[1]);
                  }
                else
                  {
                    fprintf(fp,"%lf %lf %lf ",normal[0],normal[1],normal[2]);
                  }
                break;
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
          fprintf(fp,"        <DataArray type=\"Float64\" Name=\"Post_Failure_Strain\" format=\"ascii\">\n");
          if(myRank==0)
            fprintf(pfp,"        <PDataArray type=\"Float64\" Name=\"Post_Failure_Strain\" format=\"ascii\"/>\n");
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
          fprintf(fp,"        <DataArray type=\"Int32\" Name=\"Currently_Yielding\" format=\"ascii\">\n");
          if(myRank==0)
            fprintf(pfp,"        <PDataArray type=\"Int32\" Name=\"Currently_Yielding\" format=\"ascii\"/>\n");
          break;
        case 8:
          fprintf(fp,"\n        </DataArray>\n");
          fprintf(fp,"        <DataArray type=\"Float64\" Name=\"Orientation\" format=\"ascii\" NumberOfComponents=\"3\">\n");
          if(myRank==0)
            fprintf(pfp,"        <PDataArray type=\"Float64\" Name=\"Orientation\" format=\"ascii\" NumberOfComponents=\"3\"/>\n");
          break;
        case 9:
          fprintf(fp,"\n        </DataArray>\n");
          fprintf(fp,"      </PointData>\n");
          fprintf(fp,"      <CellData>\n");
          fprintf(fp,"      </CellData>\n");
          fprintf(fp,"      <Cells>\n");
          fprintf(fp,"        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
          for(uint i=0;i<(num_particles-1)/stepping+1;++i)
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
      for(int i=0;i<nprocs;++i)
        fprintf(pfp,"    <Piece Source=\"%s.%d.%05d.vtu\"/>\n",
                swarm->name,i,timeStep);
      fprintf(pfp,"  </PUnstructuredGrid>\n\
</VTKFile>\n");
      fclose(pfp);
    }
}

/* Print out the coordinates of the mesh. */

void VTKOutput_print_coords(FILE *fp, FeMesh *feMesh, Grid *grid, int nDims,
                            uint lower[3], uint upper[3]) {
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


/* Everything is stored on cell vertices. */
void VTKOutput_fields(void *context, int myRank, int nprocs,
                      const int timeStep) {
  
  FiniteElementContext*     self = (FiniteElementContext*) context;
  Index var_I;
  int header_printed=0;
  int nDims, i;
  uint lower[3], upper[3];

  HydrostaticTerm *hydrostaticTerm;
  char* field_filename;
  FILE *field_fp, *pfield_fp;

  Dictionary_Entry_Value* field_list;
  int field_list_size;

  /* We need to save the grids to map between 1D and 3D indices for
     the values. */
  Grid *vertGrid;

  hydrostaticTerm =
    (HydrostaticTerm*)LiveComponentRegister_Get(self->CF->LCRegister,
                                                "hydrostaticTerm" );
  /* Open the file */

  field_list=Dictionary_Get(self->dictionary, "VTKOutput_FieldList" );
  if(field_list)
    {
      field_list_size=Dictionary_Entry_Value_GetCount(field_list);
      Journal_Firewall(field_list_size>=0,
                       Journal_Register( Error_Type, self->type ),
                       "The list of fields in VTKOutput_FieldList is %d but must not be negative.\n",
                       field_list_size );
    }

  Stg_asprintf( &field_filename, "%s/fields.%d.%05d.vts", self->outputPath,
                myRank, timeStep);
  field_fp=fopen(field_filename,"w");
  Memory_Free( field_filename );

  if(field_fp==NULL)
    {
      fprintf(stderr,"WARNING: Can not open fields file for rank %d step %d\n",
              myRank,timeStep);
      return;
    }

  /* Print out the parallel control files if rank==0 */
  if(myRank==0)
    {
      Stg_asprintf( &field_filename, "%s/fields.%05d.pvts", self->outputPath,
                    timeStep);
      pfield_fp=fopen(field_filename,"w");
      Memory_Free( field_filename );
      if(pfield_fp==NULL)
        {
          fprintf(stderr,
                  "FATAL ERROR: Can not open master fields file for step %d\n",
                  timeStep);
          abort();
        }
    }

  /* First, update the fields that are derived from particles.  This
     is needed to get an accurate pressure, because it needs the
     stresses. */

  for ( var_I = 0; var_I < self->fieldVariable_Register->objects->count;
        var_I++ ) {
    FieldVariable* fieldVar;
    fieldVar = FieldVariable_Register_GetByIndex( self->fieldVariable_Register,
                                                  var_I );
    if (Stg_Class_IsInstance( fieldVar, ParticleFeVariable_Type ))
      {
        ParticleFeVariable_Update( fieldVar );
      }
  }

  /* First, output the coordinates.  We have to do a huge song and
     dance just to get the extents of the mesh. */

  for ( var_I = 0; var_I < self->fieldVariable_Register->objects->count;
        var_I++ ) {
    int fields;
    FieldVariable* fieldVar;
    fieldVar = FieldVariable_Register_GetByIndex( self->fieldVariable_Register,
                                                  var_I );
    if (Stg_Class_IsInstance( fieldVar, FeVariable_Type )
        && ((FeVariable*)fieldVar)->feMesh->feElFamily!=std::string("linear-inner")) {
      FeVariable* feVar;
      Dof_Index          dofAtEachNodeCount;
      IJK ijk;

      feVar=(FeVariable*)fieldVar;

      if(!header_printed)
        {
          Mesh *mesh;
          CartesianGenerator *gen;

          mesh=(Mesh*)(feVar->feMesh);
          gen=((CartesianGenerator *)(mesh->generator));
          /* If we got the surface adaptor instead of the cartesian
             mesh generator, go to the mesh generator.  */
          if(!strcmp(gen->type,"SurfaceAdaptor"))
            gen=(CartesianGenerator *)((SurfaceAdaptor *)(gen))->generator;

          vertGrid=gen->vertGrid;
          nDims=gen->nDims;

          lower[2]=0;
          upper[2]=1;
          for(i=0;i<nDims;++i)
            {
              lower[i]=gen->vertOrigin[i];
              upper[i]=lower[i]+gen->vertRange[i];

              /* The grid is split up by elements, so for the element
                 grid, we simply add the ghost zones.  Ghost zones are
                 only added at the top, not the bottom. */
              if(upper[i]!=gen->vertGrid->sizes[i])
                ++upper[i];

              if(lower[i]!=0)
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

          /* Write the coordinates for the fields */
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
            }
          header_printed=1;
        }
      
      /* Finally, output the fields. */
      
      /* Check whether the field is in the field list */
      if(field_list)
        {
          for(fields=0; fields<field_list_size; ++fields)
            {
              if(!strcmp(Dictionary_Entry_Value_AsString(Dictionary_Entry_Value_GetElement(field_list,fields)),
                         fieldVar->name))
                {
                  break;
                }
            }
          if(fields==field_list_size)
            continue;
        }

      dofAtEachNodeCount = feVar->fieldComponentCount;

      switch(dofAtEachNodeCount*nDims)
        {
          /* Scalars */
        case 2:
        case 3:
          fprintf(field_fp,"        <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n",feVar->name);
          if(myRank==0)
            fprintf(pfield_fp,"      <PDataArray type=\"Float64\" Name=\"%s\" format=\"ascii\"/>\n",feVar->name);
          break;
          /* Vectors */
        case 4:
        case 9:
          fprintf(field_fp,"        <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\" NumberOfComponents=\"3\">\n",feVar->name);
          if(myRank==0)
            fprintf(pfield_fp,"      <PDataArray type=\"Float64\" Name=\"%s\" format=\"ascii\" NumberOfComponents=\"3\"/>\n",feVar->name);
          break;
          /* Rank 2 Tensors */
        case 6:
        case 8:
        case 18:
        case 27:
          fprintf(field_fp,"        <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\" NumberOfComponents=\"9\">\n",feVar->name);
          if(myRank==0)
            fprintf(pfield_fp,"      <PDataArray type=\"Float64\" Name=\"%s\" format=\"ascii\" NumberOfComponents=\"9\"/>\n",feVar->name);
          break;
          /* Unknown */
        default:
          fprintf(stderr,
                  "Bad number of degrees of freedom for variable %s: %d\n",
                 feVar->name, dofAtEachNodeCount);
          abort();
        }

      for(ijk[2]=lower[2];ijk[2]<upper[2];++ijk[2])
        for(ijk[1]=lower[1];ijk[1]<upper[1];++ijk[1])
          for(ijk[0]=lower[0];ijk[0]<upper[0];++ijk[0])
            {
              double variableValues[MAX_FIELD_COMPONENTS];	
              Dof_Index          dof_I;
              unsigned local;
              Mesh_GlobalToDomain(feVar->feMesh,MT_VERTEX,
                                  Grid_Project(vertGrid,ijk),&local);
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
                  /* Special case the pressure */
                  if(!strcmp(feVar->name,"PressureField"))
                    {
                      double p;
                      /* First add the trace of the stress */

                      for(uint stress_I = 0;
                          stress_I<self->fieldVariable_Register->objects->count;
                          stress_I++ ) {
                        FieldVariable* stressVar;
                        stressVar =
                          FieldVariable_Register_GetByIndex(self->fieldVariable_Register,
                                                            stress_I );
                        if(!strcmp(stressVar->name,"StressField"))
                          {
                            if (Stg_Class_IsInstance( stressVar,
                                                      FeVariable_Type ))
                              {
                                FeVariable* sVar;
                                double stressValues[MAX_FIELD_COMPONENTS];	
                                sVar=(FeVariable*)stressVar;
                                
                                FeVariable_GetValueAtNode(sVar,local,stressValues);
                                if(nDims==2)
                                  {
                                    p=(stressValues[0]+stressValues[1])/2;
                                  }
                                else
                                  {
                                    p=(stressValues[0]+stressValues[1]
                                       +stressValues[2])/3;
                                  }
                              }
                            break;
                          }
                      }
                      /* Next add the hydrostatic term */
                      if(hydrostaticTerm)
                        {
                          double *coord;
                          coord=Mesh_GetVertex(feVar->feMesh,local);
                          p+=HydrostaticTerm_Pressure(hydrostaticTerm,coord);
                        }
                      p+=variableValues[0];
                      fprintf(field_fp, "%.15g ", p );
                    }
                  else
                    {
                      for ( dof_I = 0; dof_I < dofAtEachNodeCount; dof_I++ ) {
                        fprintf(field_fp, "%.15g ", variableValues[dof_I] );
                      }
                    }
                  break;
                case 4:
                  fprintf(field_fp, "%.15g %.15g 0", variableValues[0],
                          variableValues[1] );
                  break;
                case 6:
                  /* Ordering is xx, yy, xy */
                  fprintf(field_fp, "%.15g %.15g 0 %.15g %.15g 0 0 0 0 ",
                          variableValues[0], variableValues[2],
                          variableValues[2], variableValues[1]);
                  break;
                case 8:
                  fprintf(field_fp, "%.15g %.15g 0 %.15g %.15g 0 0 0 0 ",
                          variableValues[0], variableValues[1],
                          variableValues[2], variableValues[3]);
                  break;
                case 18:
                  /* Ordering is xx, yy, zz, xy, xz, yz */
                  fprintf(field_fp, "%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g ",
                          variableValues[0], variableValues[3],
                          variableValues[4], variableValues[3],
                          variableValues[1], variableValues[5],
                          variableValues[4], variableValues[5],
                          variableValues[2]);
                  break;
                }
              fprintf(field_fp, "\n" );
            }
      fprintf(field_fp,"        </DataArray>\n");
    }
  }
  fprintf(field_fp,"      </PointData>\n    </Piece>\n\
  </StructuredGrid>\n\
</VTKFile>\n");
  fclose(field_fp);

  if(myRank==0)
    {
      fprintf(pfield_fp,"      </PPointData>\n");
      for(i=0;i<nprocs;++i)
        {
          fprintf(pfield_fp,"    <Piece Extent=\"%d %d %d %d %d %d\"\n\
             Source=\"fields.%d.%05d.vts\"/>\n",lower[0],upper[0]-1,
                  lower[1],upper[1]-1,lower[2],upper[2]-1,
                  i,timeStep);
        }
      fprintf(pfield_fp,"  </PStructuredGrid>\n\
</VTKFile>\n");
      fclose(pfield_fp);
     }

}
