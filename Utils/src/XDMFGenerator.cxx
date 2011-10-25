/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
** Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
** Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
** Australian Computational Earth Systems Simulator - http://www.access.edu.au
** Monash Cluster Computing - http://www.mcc.monash.edu.au
** Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
** Contributors:
** Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
** Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
** Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
** David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
** Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
** Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
** Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
** Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
** Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
** Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
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
** $Id: $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include "XDMFGenerator.h"

#include <stdlib.h>
#include <string.h>
#include <string>
#include <sstream>

#ifndef MASTER
	#define MASTER 0
#endif

const Type XDMFGenerator_Type = "XDMFGenerator";

void XDMFGenerator_GenerateAll( void* _context ) {
   UnderworldContext*   context    = (UnderworldContext*)_context;
   Stream*              stream;

   /** only the MASTER process writes to the file.  other processes send information to the MASTER where required **/   
   if(context->rank == MASTER) {
      Bool                 fileOpened;
      Stream*              errorStream  = Journal_Register( Error_Type, (Name)CURR_MODULE_NAME  );
      char*                 filename;
      char*                outputPathString;

      /** Create Stream **/
      stream = Journal_Register( InfoStream_Type, (Name)"XDMFOutputFile"  );
   
      /** Set auto flush on stream **/
      Stream_SetAutoFlush( stream, True );

		outputPathString = Context_GetCheckPointWritePrefixString( context );
   
      /** Get name of XDMF schema file **/
      Stg_asprintf( &filename, "%s/XDMF.%05d.xmf", outputPathString, context->timeStep );
   
      /** Init file, always overwriting any existing **/
      fileOpened = Stream_RedirectFile( stream, filename );
      Journal_Firewall( fileOpened, errorStream, 
            "Could not open file %s. Possibly directory %s does not exist or is not writable.\n"
            "Check 'checkpointWritePath' in input file.\n", filename, outputPathString );
      Memory_Free( filename );
      Memory_Free( outputPathString );
      
      /** write header information **/
      _XDMFGenerator_WriteHeader( context, stream );

      /** Write all (checkpointed) FeVariable field information  **/
      _XDMFGenerator_WriteFieldSchema( context, stream );

      /** Write all (checkpointed) FeVariable field information  **/
      if ( (context->timeStep % context->checkpointEvery == 0) || 
           (context->checkpointAtTimeInc && (context->currentTime >= context->nextCheckpointTime)) ) _XDMFGenerator_WriteSwarmSchema( context, stream);

      /** writes footer information and close file/stream **/
      _XDMFGenerator_WriteFooter( context, stream );

      /** close the file **/
      Stream_CloseFile( stream );	

   } else {
      /** other process send information about swarms populations to MASTER **/
      _XDMFGenerator_SendInfo( context );
   }

}

namespace {
  void print_Pm1_corner(Stream* stream,
                        const int &num_columns,
                        const int &num_rows,
                        const int &e,
                        const std::string &mesh_name,
                        const std::string &x,
                        const std::string &y)
  {
    Journal_Printf(stream,"    <DataItem ItemType=\"Function\"  Dimensions=\"%u\" Function=\"JOIN(",
                   num_rows*num_columns);
    for(int i=0;i<num_columns-1;++i)
      Journal_Printf(stream,"$%u; ",i);
    Journal_Printf(stream,"$%u)\" Name=\"C%u%s\">\n",
                   num_columns-1,e,(x+y).c_str());

    for(int i=0;i<num_columns;++i)
      {
        Journal_Printf(stream,"      <DataItem Reference=\"XML\">/Xdmf/Domain/Grid[@Name=\"FEM_Grid_%s\"]/DataItem[@Name=\"C%u%s%u\"]</DataItem>\n",
                       mesh_name.c_str(),e,x.c_str(),i);
      }
    Journal_Printf(stream,"    </DataItem>\n");
  }
}

namespace {
  void write_continuous_geometry(Stream *stream, FeMesh *feMesh,
                                 const std::string &element_family,
                                 const int &timestep,
                                 const int &nDims,
                                 const std::string &topologyType);
  void write_continuous_FeVariable(Stream *stream,
                                   FeMesh *feMesh,
                                   FeVariable *feVar,
                                   const int &timestep,
                                   const int &nDims,
                                   Bool saveCoords);
  void write_discontinuous_geometry(Stream *stream, Mesh *mesh,
                                    const int &timestep,
                                    const int &nDims);
  void write_discontinuous_FeVariable(Stream *stream,
                                      const std::string &mesh_name,
                                      const std::string &var_name,
                                      const int &timestep,
                                      const int &num_elements);
}

void _XDMFGenerator_WriteFieldSchema( UnderworldContext* context, Stream* stream ) {
  FieldVariable*       fieldVar = NULL;
  FeVariable*          feVar    = NULL;
  Mesh*                mesh     = NULL;
  FeMesh*              feMesh   = NULL;
  Index                var_I = 0;
  std::string          topologyType;
  unsigned             componentCount = LiveComponentRegister_GetCount(stgLiveComponentRegister);
  unsigned             compI;
  Stg_Component*       stgComp;

  /** search for entire live component register for feMesh types  **/
  for( compI = 0 ; compI < componentCount ; compI++ ){
    stgComp = LiveComponentRegister_At( stgLiveComponentRegister, compI );
    /* check that component is of type FeMesh, and that its element
       family is supported */

    std::string element_family(Stg_Class_IsInstance(stgComp,FeMesh_Type) ?
                               ((FeMesh*)stgComp)->feElFamily : "");

    if(element_family=="linear" || element_family=="linear-inner"
       || element_family=="quadratic") {
      mesh   = (  Mesh*)stgComp;
      feMesh = (FeMesh*)stgComp;
      const unsigned nDims=Mesh_GetDimSize( mesh );
      int elementGlobalSize = FeMesh_GetElementGlobalSize(feMesh);

      /** now write all the xdmf geometry info **/
      /**----------------------- START GEOMETRY   --------------------------- **/
      if(element_family=="linear" || element_family=="quadratic")
        {
          topologyType= (nDims==2) ? "Quadrilateral" : "Hexahedron";
        }
      else if(element_family=="linear-inner")
        {
          topologyType= (nDims==2) ? "Triangle" : "Tetrahedron";
        }
      /** first create the grid which is applicable to the
          checkpointed fevariables **/
      Journal_Printf( stream, "  <Grid Name=\"FEM_Grid_%s\">\n\n", feMesh->name);
      Journal_Printf( stream, "    <Time Value=\"%f\" />\n\n", context->currentTime );
      /** now print out topology info **/

      /* Discontinuous elements are written separately. */
      if(element_family=="linear-inner" && nDims==2)
        {
          write_discontinuous_geometry(stream,mesh,context->timeStep,nDims);
        }
      else
        {
          write_continuous_geometry(stream,feMesh,element_family,
                                    context->timeStep,nDims,topologyType);
        }

      /** now write FeVariable data **/   
      for ( var_I = 0; var_I < context->fieldVariable_Register->objects->count;
            var_I++ ) {
        fieldVar =
          FieldVariable_Register_GetByIndex(context->fieldVariable_Register,
                                            var_I);

        if ( Stg_Class_IsInstance( fieldVar, FeVariable_Type ) ) {
          feVar = (FeVariable*)fieldVar;
          if ( (feVar->isCheckpointedAndReloaded
                && (context->timeStep % context->checkpointEvery == 0))
               || (feVar->isCheckpointedAndReloaded
                   && (context->checkpointAtTimeInc
                       && (context->currentTime >= context->nextCheckpointTime)))
               || (feVar->isSavedData
                   && (context->timeStep % context->saveDataEvery   == 0)) ){
            FeMesh* feVarMesh = NULL;
            /** check what type of generator was used to know where
                elementMesh is **/
            if(Stg_Class_IsInstance(feVar->feMesh->generator,C0Generator_Type))
              feVarMesh=
                (FeMesh*)(((C0Generator*)feVar->feMesh->generator)->elMesh);
            if(Stg_Class_IsInstance(feVar->feMesh->generator,
                                    CartesianGenerator_Type)
               || Stg_Class_IsInstance(feVar->feMesh->generator,MeshAdaptor_Type)
               || Stg_Class_IsInstance(feVar->feMesh->generator,
                                       InnerGenerator_Type))
              feVarMesh = feVar->feMesh;

            /** make sure that the fevariable femesh is the same as
                that used above for the geometry definition, if so
                proceed **/
            if( feVarMesh == feMesh ){
              if(element_family=="linear-inner" && nDims==2)
                {
                  write_discontinuous_FeVariable(stream,feMesh->name,feVar->name,
                                                 context->timeStep,
                                                 elementGlobalSize);
                }
              else
                {
                  write_continuous_FeVariable
                    (stream,feMesh,feVar,context->timeStep,nDims,
                     Dictionary_GetBool_WithDefault
                     (context->dictionary,
                      (Dictionary_Entry_Key)"saveCoordsWithFields",False));
                      
                }
            }
          }
        }
      }
      Journal_Printf( stream, "  </Grid>\n\n" );
    }
  }
}

void _XDMFGenerator_WriteSwarmSchema( UnderworldContext* context, Stream* stream ) {
   Swarm_Register* swarmRegister = Swarm_Register_GetSwarm_Register();
   Index           swarmCount;
   Index           swarmcountindex;
   Index           variablecount;
   Index           dofCountIndex;
   Index           swarmParticleLocalCount;
   Index           countindex;
   Swarm*          currentSwarm;
   SwarmVariable*  swarmVar;
   char*           swarmVarName;
   char*           variableType = NULL;
   char*           filename_part = NULL;
   Stream*         errorStream  = Journal_Register( Error_Type, (Name)CURR_MODULE_NAME  );
	const int       FINISHED_WRITING_TAG = 100;
	MPI_Status      status;
   
   /** get total number of different swarms **/
   swarmCount  = swarmRegister->swarmList->count;
   
   /** parse swarm list, checking which are actually stored to HDF5.  **/
   /** We assume that all processes have the same number of swarms on them, with the same swarms checkpointed**/
   for(swarmcountindex = 0; swarmcountindex < swarmCount; ++swarmcountindex){
      currentSwarm  = Swarm_Register_At( swarmRegister, swarmcountindex );

      if ( currentSwarm->isSwarmTypeToCheckPointAndReload != True ) 
         continue;

      swarmParticleLocalCount = currentSwarm->particleLocalCount;
      /** first create a grid collection which will contain the collection of swarms from each process. 
          there will be one of these collections for each swarm (not swarmvariable) that is checkpointed **/
                              Journal_Printf( stream, "   <Grid Name=\"%s\" GridType=\"Collection\">\n\n", currentSwarm->name );
                              Journal_Printf( stream, "      <Time Value=\"%f\" />\n\n", context->currentTime );

      Stg_asprintf( &filename_part, "" );
      for (int ii = 0 ; ii < context->nproc ; ++ii) {

         /** get the number of particles in each swarm from each process **/
         if (ii != MASTER        ) MPI_Recv( &swarmParticleLocalCount, 1, MPI_INT, ii, FINISHED_WRITING_TAG, context->communicator, &status );
         if (context->nproc != 1 ) Stg_asprintf( &filename_part, ".%uof%u", (ii + 1), context->nproc );
         
         /** first write all the MASTER procs swarm info **/
         if (swarmParticleLocalCount != 0) {
                              Journal_Printf( stream, "      <Grid Name=\"%s_proc_%u\">\n\n", currentSwarm->name, ii );
            
            /** now write all the xdmf geometry info **/
            /**----------------------- START GEOMETRY   ------------------------------------------------------------------------------------------------------------------- **/
                              Journal_Printf( stream, "         <Topology Type=\"POLYVERTEX\" NodesPerElement=\"%u\"> </Topology>\n", swarmParticleLocalCount );
                              Journal_Printf( stream, "         <Geometry Type=\"XYZ\">\n" );
            Stg_asprintf( &swarmVarName, "%s-Position", currentSwarm->name );
            swarmVar = SwarmVariable_Register_GetByName( currentSwarm->swarmVariable_Register, swarmVarName );
            if (!swarmVar) 
               Journal_DPrintf( errorStream, "\n\n Error: Could not find required Position SwarmVariable. \n\n" );
            Memory_Free( swarmVarName );
            
            /** check what precision it Position variable is stored at **/
            if(        swarmVar->variable->dataTypes[0] == Variable_DataType_Int ){
               Journal_DPrintf( errorStream, "\n\n Error: Position variable can not be of type Int. \n\n" );
            } else if ( swarmVar->variable->dataTypes[0] == Variable_DataType_Char){
               Journal_DPrintf( errorStream, "\n\n Error: Position variable can not be of type Char. \n\n" );
            } else if ( swarmVar->variable->dataTypes[0] == Variable_DataType_Float ){
               Stg_asprintf( &variableType, "NumberType=\"Float\" Precision=\"4\"" );
            } else {
               Stg_asprintf( &variableType, "NumberType=\"Float\" Precision=\"8\"" );
            }
           
            if(         swarmVar->dofCount == 2 ){
               /** note that for 2d, we feed back a quasi 3d array, with the 3rd Dof zeroed.  so in effect we always work in 3d.
                   this is done because paraview in particular seems to do everything in 3d, and if you try and give it a 2d vector 
                   or array, it complains.... and hence the verbosity of the following 2d definitions**/
                              Journal_Printf( stream, "            <DataItem ItemType=\"Function\"  Dimensions=\"%u 3\" Function=\"JOIN($0, $1, 0*$1)\">\n", swarmParticleLocalCount );
                              Journal_Printf( stream, "               <DataItem ItemType=\"HyperSlab\" Dimensions=\"%u 1\" Name=\"XCoords\">\n", swarmParticleLocalCount );
                              Journal_Printf( stream, "                  <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 0 1 1 %u 1 </DataItem>\n", swarmParticleLocalCount );
                              Journal_Printf( stream, "                  <DataItem Format=\"HDF\" %s Dimensions=\"%u 2\">%s.%05d%s.h5:/Position</DataItem>\n", variableType, swarmParticleLocalCount, currentSwarm->name, context->timeStep, filename_part );
                              Journal_Printf( stream, "               </DataItem>\n" );
                              Journal_Printf( stream, "               <DataItem ItemType=\"HyperSlab\" Dimensions=\"%u 1\" Name=\"YCoords\">\n", swarmParticleLocalCount );
                              Journal_Printf( stream, "                  <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 1 1 1 %u 1 </DataItem>\n", swarmParticleLocalCount );
                              Journal_Printf( stream, "                  <DataItem Format=\"HDF\" %s Dimensions=\"%u 2\">%s.%05d%s.h5:/Position</DataItem>\n", variableType, swarmParticleLocalCount, currentSwarm->name, context->timeStep, filename_part );
                              Journal_Printf( stream, "               </DataItem>\n" );
                              Journal_Printf( stream, "            </DataItem>\n" );
            } else if ( swarmVar->dofCount == 3 ) {
               /** in 3d we simply feed back the 3d hdf5 array, nice and easy **/
                              Journal_Printf( stream, "            <DataItem Format=\"HDF\" %s Dimensions=\"%u 3\">%s.%05d%s.h5:/Position</DataItem>\n", variableType, swarmParticleLocalCount, currentSwarm->name, context->timeStep, filename_part );
            } else {
               Journal_DPrintf( errorStream, "\n\n Error: Position SwarmVariable is not of dofCount 2 or 3.\n\n" );
            }
                              Journal_Printf( stream, "         </Geometry>\n\n" );
            /**----------------------- FINISH GEOMETRY  ------------------------------------------------------------------------------------------------------------------- **/
   
            
            /** now write all the swarm attributes.. ie all the checkpointed swarmVariables **/
            
            variablecount = currentSwarm->swarmVariable_Register->objects->count;
            for(countindex = 0; countindex < variablecount; ++countindex){
               swarmVar = SwarmVariable_Register_GetByIndex( currentSwarm->swarmVariable_Register, countindex );
               if( swarmVar->isCheckpointedAndReloaded ) {
            /**----------------------- START ATTRIBUTES ------------------------------------------------------------------------------------------------------------------- **/
                  if(         swarmVar->variable->dataTypes[0] == Variable_DataType_Int ){
                     Stg_asprintf( &variableType, "NumberType=\"Int\"" );
                  } else if ( swarmVar->variable->dataTypes[0] == Variable_DataType_Char){
                     Stg_asprintf( &variableType, "NumberType=\"Char\"" );
                  } else if ( swarmVar->variable->dataTypes[0] == Variable_DataType_Float ){
                     Stg_asprintf( &variableType, "NumberType=\"Float\" Precision=\"4\"" );
                  } else {
                     Stg_asprintf( &variableType, "NumberType=\"Float\" Precision=\"8\"" );
                  }
                  if (        swarmVar->dofCount == 1 ) {
                              Journal_Printf( stream, "         <Attribute Type=\"Scalar\" Center=\"Node\" Name=\"%s\">\n", swarmVar->name);
                              Journal_Printf( stream, "            <DataItem Format=\"HDF\" %s Dimensions=\"%u 1\">%s.%05d%s.h5:/%s</DataItem>\n", variableType, swarmParticleLocalCount, currentSwarm->name, context->timeStep, filename_part, swarmVar->name + strlen(currentSwarm->name)+1 );
                              Journal_Printf( stream, "         </Attribute>\n\n" );
                  } else if ( swarmVar->dofCount == 2 ){
                     /** note that for 2d, we feed back a quasi 3d array, with the 3rd Dof zeroed.  so in effect we always work in 3d.
                         this is done because paraview in particular seems to do everything in 3d, and if you try and give it a 2d vector 
                         or array, it complains.... and hence the verbosity of the following 2d definitions **/
                              Journal_Printf( stream, "         <Attribute Type=\"Vector\" Center=\"Node\" Name=\"%s\">\n", swarmVar->name);
                              Journal_Printf( stream, "            <DataItem ItemType=\"Function\"  Dimensions=\"%u 3\" Function=\"JOIN($0, $1, 0*$1)\">\n", swarmParticleLocalCount );
                              Journal_Printf( stream, "               <DataItem ItemType=\"HyperSlab\" Dimensions=\"%u 1\" Name=\"XValue\">\n", swarmParticleLocalCount );
                              Journal_Printf( stream, "                  <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 0 1 1 %u 1 </DataItem>\n", swarmParticleLocalCount );
                              Journal_Printf( stream, "                  <DataItem Format=\"HDF\" %s Dimensions=\"%u 2\">%s.%05d%s.h5:/%s</DataItem>\n", variableType, swarmParticleLocalCount, currentSwarm->name, context->timeStep, filename_part, swarmVar->name + strlen(currentSwarm->name)+1 );
                              Journal_Printf( stream, "               </DataItem>\n" );
                              Journal_Printf( stream, "               <DataItem ItemType=\"HyperSlab\" Dimensions=\"%u 1\" Name=\"YValue\">\n", swarmParticleLocalCount );
                              Journal_Printf( stream, "                  <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 1 1 1 %u 1 </DataItem>\n", swarmParticleLocalCount );
                              Journal_Printf( stream, "                  <DataItem Format=\"HDF\" %s Dimensions=\"%u 2\">%s.%05d%s.h5:/%s</DataItem>\n", variableType, swarmParticleLocalCount, currentSwarm->name, context->timeStep, filename_part, swarmVar->name + strlen(currentSwarm->name)+1 );
                              Journal_Printf( stream, "               </DataItem>\n" );
                              Journal_Printf( stream, "            </DataItem>\n" );
                              Journal_Printf( stream, "         </Attribute>\n\n" );
                  } else if ( swarmVar->dofCount == 3 ) {
                     /** in 3d we simply feed back the 3d hdf5 array, nice and easy **/
                              Journal_Printf( stream, "         <Attribute Type=\"Vector\" Center=\"Node\" Name=\"%s\">\n", swarmVar->name);
                              Journal_Printf( stream, "            <DataItem Format=\"HDF\" %s Dimensions=\"%u 3\">%s.%05d%s.h5:/%s</DataItem>\n", variableType, swarmParticleLocalCount, currentSwarm->name, context->timeStep, filename_part, swarmVar->name + strlen(currentSwarm->name)+1 );
                              Journal_Printf( stream, "         </Attribute>\n\n" );
                  } else {
                     /** where there are more than 3 components, we write each one out as a scalar **/
                     for(dofCountIndex = 0 ; dofCountIndex < swarmVar->dofCount ; ++dofCountIndex){
                              Journal_Printf( stream, "         <Attribute Type=\"Scalar\" Center=\"Node\" Name=\"%s-Component-%u\">\n", swarmVar->name, dofCountIndex);
                              Journal_Printf( stream, "            <DataItem ItemType=\"HyperSlab\" Dimensions=\"%u 1\" >\n", swarmParticleLocalCount );
                              Journal_Printf( stream, "               <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 %u 1 1 %u 1 </DataItem>\n", dofCountIndex, swarmParticleLocalCount );
                              Journal_Printf( stream, "               <DataItem Format=\"HDF\" %s Dimensions=\"%u %u\">%s.%05d%s.h5:/%s</DataItem>\n", variableType, swarmParticleLocalCount, swarmVar->dofCount, currentSwarm->name, context->timeStep, filename_part, swarmVar->name + strlen(currentSwarm->name)+1 );
                              Journal_Printf( stream, "            </DataItem>\n" );
                              Journal_Printf( stream, "         </Attribute>\n" );
                     }
                              Journal_Printf( stream, "\n" );
                  }
            /**----------------------- END ATTRIBUTES   ------------------------------------------------------------------------------------------------------------------- **/
               }
         
            }
   
                              Journal_Printf( stream, "      </Grid>\n\n" );
   
         }
      
      }

                              Journal_Printf( stream, "  </Grid>\n\n" );

   if(variableType)  Memory_Free( variableType );
   if(filename_part) Memory_Free( filename_part );
      
   }


}

void  _XDMFGenerator_WriteHeader( UnderworldContext* context, Stream* stream ) {
   
	/** Print XDMF header info **/
                              Journal_Printf( stream, "<?xml version=\"1.0\" ?>\n" );
                              Journal_Printf( stream, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n" );
                              Journal_Printf( stream, "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n" );
                              Journal_Printf( stream, "\n" );
                              Journal_Printf( stream, "<Domain>\n" );
                              Journal_Printf( stream, "\n" );

}


void _XDMFGenerator_WriteFooter( UnderworldContext* context, Stream* stream ) {

                              Journal_Printf( stream, "</Domain>\n" );
                              Journal_Printf( stream, "\n" );
                              Journal_Printf( stream, "</Xdmf>\n" );
                              Journal_Printf( stream, "\n" );

}

void _XDMFGenerator_SendInfo( UnderworldContext* context ) {
   Swarm_Register* swarmRegister = Swarm_Register_GetSwarm_Register();
   Index           swarmCount;
   Index           swarmcountindex;
   Index           swarmParticleLocalCount;
   Swarm*          currentSwarm;
	const int       FINISHED_WRITING_TAG = 100;
   
   /** get total number of different swarms **/
   swarmCount  = swarmRegister->swarmList->count;
   
   /** parse swarm list, checking which are actually stored to HDF5.  **/
   /** We assume that all processes have the same number of swarms on them, with the same swarms checkpointed**/
   for(swarmcountindex = 0; swarmcountindex < swarmCount; ++swarmcountindex){
      currentSwarm  = Swarm_Register_At( swarmRegister, swarmcountindex );

      if ( currentSwarm->isSwarmTypeToCheckPointAndReload != True ) 
         continue;

      swarmParticleLocalCount = currentSwarm->particleLocalCount;
		MPI_Ssend( &swarmParticleLocalCount, 1, MPI_INT, MASTER, FINISHED_WRITING_TAG, context->communicator );
   }
}


namespace {

  void write_continuous_geometry(Stream *stream, FeMesh *feMesh,
                                 const std::string &element_family,
                                 const int &timestep,
                                 const int &nDims,
                                 const std::string &topologyType)
  {
    int totalVerts        = Mesh_GetGlobalSize( feMesh, (MeshTopology_Dim)0 );
    int elementGlobalSize = FeMesh_GetElementGlobalSize(feMesh);

    unsigned int maxNodes;
    /* get connectivity array size */
    if (feMesh->nElTypes == 1)
      maxNodes = FeMesh_GetElementNodeSize( feMesh, 0);
    else {
      /* determine the maximum number of nodes each element has */
      maxNodes = 0;
      for (unsigned int gElement_I = 0; gElement_I < FeMesh_GetElementGlobalSize(feMesh);
           gElement_I++ ) {
        unsigned int numNodes;
        numNodes = FeMesh_GetElementNodeSize( feMesh, gElement_I);
        if( maxNodes < numNodes ) maxNodes = numNodes;
      }
    }
    /* Quadratic elements */
    if(element_family=="quadratic")
      {
        /* First separate each node of each element into its own
           variable C0-C8 or C0-C26 */
        for(unsigned int i=0;i<maxNodes;++i)
          {
            Journal_Printf(stream,"    <DataItem ItemType=\"HyperSlab\" Dimensions=\"%u 1\" Name=\"C%u\">\n",
                           elementGlobalSize, i);
            Journal_Printf(stream,"      <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 %u 1 1 %u 1 </DataItem>\n",
                           i, elementGlobalSize);
            Journal_Printf(stream,"      <DataItem Format=\"HDF\" DataType=\"Int\"  Dimensions=\"%u %u\">Mesh.%s.%05d.h5:/connectivity</DataItem>\n",
                           elementGlobalSize, maxNodes, feMesh->name, timestep );
            Journal_Printf(stream,"    </DataItem>\n");
          }
        /* Print out the topology */
        int num_nodes=(maxNodes==9) ? 4 : 8;
        Journal_Printf(stream,"    <Topology Type=\"%s\" NumberOfElements=\"%u\"> \n",
                       topologyType.c_str(),elementGlobalSize*num_nodes);
        Journal_Printf(stream,"      <DataItem ItemType=\"Function\"  Dimensions=\"%u %u\" Function=\"JOIN(",
                       elementGlobalSize*num_nodes, num_nodes);
        for(int j=0;j<num_nodes-1;++j)
          {
            Journal_Printf(stream,"$%u; ",j);
          }
        Journal_Printf(stream,"$%u)\">\n",num_nodes-1);

        int first_mapping[8]={0,1,4,3,9,10,13,12};
        int mapping[8][8];
        for(int i=0;i<8;++i)
          {
            mapping[0][i]=first_mapping[i];
            mapping[1][i]=first_mapping[i]+1;
            mapping[2][i]=first_mapping[i]+3;
            mapping[3][i]=first_mapping[i]+4;
            mapping[4][i]=first_mapping[i]+9;
            mapping[5][i]=first_mapping[i]+10;
            mapping[6][i]=first_mapping[i]+12;
            mapping[7][i]=first_mapping[i]+13;
          }

        /* Print out the quadrilaterals made up of only some of the
           nodes */
        for(int j=0;j<num_nodes;++j)
          {
            Journal_Printf(stream,"        <DataItem ItemType=\"Function\"  Dimensions=\"%u %u\" Function=\"JOIN(",
                           elementGlobalSize, num_nodes);
            for(int k=0;k<num_nodes-1;++k)
              {
                Journal_Printf(stream,"$%u, ",k);
              }
            Journal_Printf(stream,"$%u)\">\n",num_nodes-1);
            for(int k=0;k<num_nodes;++k)
              {
                Journal_Printf(stream,"          <DataItem Reference=\"XML\">/Xdmf/Domain/Grid[@Name=\"FEM_Grid_%s\"]/DataItem[@Name=\"C%u\"] </DataItem>\n",
                               feMesh->name, mapping[j][k]);
              }
            Journal_Printf(stream,"        </DataItem>\n");
          }
        Journal_Printf(stream,"      </DataItem>\n");
        Journal_Printf(stream,"    </Topology>\n\n" );
      }

    /* Continuous linear elements */
    else
      {
        Journal_Printf( stream, "         <Topology Type=\"%s\" NumberOfElements=\"%u\"> \n",
                        topologyType.c_str(), elementGlobalSize );
        Journal_Printf( stream, "            <DataItem Format=\"HDF\" DataType=\"Int\"  Dimensions=\"%u %u\">Mesh.%s.%05d.h5:/connectivity</DataItem>\n",
                        elementGlobalSize, maxNodes, feMesh->name,
                        timestep );
        Journal_Printf( stream, "         </Topology>\n\n" );
      }

    /* Print out coordinates of the vertices */
    Journal_Printf( stream, "    <Geometry Type=\"XYZ\">\n" );
    std::string variableType("NumberType=\"Float\" Precision=\"8\"");
    if(nDims==2){
      /** note that for 2d, we feed back a quasi 3d array, with the
          3rd Dof zeroed.  so in effect we always work in 3d.  this
          is done because paraview in particular seems to do
          everything in 3d, and if you try and give it a 2d vector
          or array, it complains.... and hence the verbosity of the
          following 2d definitions**/
      Journal_Printf( stream, "      <DataItem ItemType=\"Function\" Dimensions=\"%u 3\" Function=\"JOIN($0, $1, 0*$1)\">\n", totalVerts );
      Journal_Printf( stream, "        <DataItem ItemType=\"HyperSlab\" Dimensions=\"%u 1\" Name=\"XCoords\">\n", totalVerts );
      Journal_Printf( stream, "          <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 0 1 1 %u 1 </DataItem>\n", totalVerts );
      Journal_Printf( stream, "          <DataItem Format=\"HDF\" %s Dimensions=\"%u 2\">Mesh.%s.%05d.h5:/vertices</DataItem>\n", variableType.c_str(), totalVerts, feMesh->name,  timestep );
      Journal_Printf( stream, "        </DataItem>\n" );
      Journal_Printf( stream, "        <DataItem ItemType=\"HyperSlab\" Dimensions=\"%u 1\" Name=\"YCoords\">\n", totalVerts );
      Journal_Printf( stream, "          <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 1 1 1 %u 1 </DataItem>\n", totalVerts );
      Journal_Printf( stream, "          <DataItem Format=\"HDF\" %s Dimensions=\"%u 2\">Mesh.%s.%05d.h5:/vertices</DataItem>\n", variableType.c_str(), totalVerts, feMesh->name, timestep );
      Journal_Printf( stream, "        </DataItem>\n" );
      Journal_Printf( stream, "      </DataItem>\n" );
    } else if ( nDims == 3 ) {
      /** in 3d we simply feed back the 3d hdf5 array, nice and easy **/
      Journal_Printf( stream, "            <DataItem Format=\"HDF\" %s Dimensions=\"%u 3\">Mesh.%s.%05d.h5:/vertices</DataItem>\n", variableType.c_str(), totalVerts, feMesh->name, timestep );
    } else {
      Stream* errorStream=Journal_Register(Error_Type,(Name)CURR_MODULE_NAME);
      Journal_DPrintf( errorStream, "\n\n Error: Mesh vertex location is not of dofCount 2 or 3.\n\n" );
    }
    Journal_Printf( stream, "    </Geometry>\n\n" );
  }

  void write_continuous_FeVariable(Stream *stream, FeMesh *feMesh,
                                   FeVariable *feVar,
                                   const int &timestep,
                                   const int &nDims,
                                   Bool saveCoords)
  {
    std::string centering;
    Index  offset = 0;
    Index  meshSize = Mesh_GetGlobalSize( feVar->feMesh, (MeshTopology_Dim)0 );

    /**----------------------- START ATTRIBUTES ----------------- **/
    /** if coordinates are being stored with feVariable,
        account for this **/
    if( saveCoords) offset = nDims; 
    /** all feVariables are currently stored as doubles **/
    std::string variableType="NumberType=\"Float\" Precision=\"8\"";
    /** determine whether feVariable data is cell centered (like
        Pressure with P0 elements), or on the nodes (like
        Velocity) **/
    unsigned int totalVerts(Mesh_GetGlobalSize(feMesh,(MeshTopology_Dim)0));
    unsigned int elementGlobalSize(FeMesh_GetElementGlobalSize(feMesh));
    if(meshSize == elementGlobalSize){
      centering="Cell";
    } else if(meshSize == totalVerts) {
      centering="Node";
    } else {
      /* unknown/unsupported type */
      centering="UNKNOWN_POSSIBLY_ERROR";
    }

    /** how many degrees of freedom does the fevariable have? **/
    unsigned int dofAtEachNodeCount = feVar->fieldComponentCount;
    if (dofAtEachNodeCount == 1) {
      Journal_Printf( stream, "    <Attribute Type=\"Scalar\" Center=\"%s\" Name=\"%s\">\n", centering.c_str(),  feVar->name);
      Journal_Printf( stream, "       <DataItem ItemType=\"HyperSlab\" Dimensions=\"%u 1\" >\n", meshSize );
      Journal_Printf( stream, "          <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 %u 1 1 %u 1 </DataItem>\n", offset, meshSize );
      Journal_Printf( stream, "          <DataItem Format=\"HDF\" %s Dimensions=\"%u %u\">%s.%05d.h5:/data</DataItem>\n", variableType.c_str(), meshSize, (offset + dofAtEachNodeCount), feVar->name, timestep);
      Journal_Printf( stream, "       </DataItem>\n" );
      Journal_Printf( stream, "    </Attribute>\n\n" );
    } else if (dofAtEachNodeCount == 2){
      /** note that for 2d, we feed back a quasi 3d array, with the 3rd Dof zeroed.  so in effect we always work in 3d.
          this is done because paraview in particular seems to do everything in 3d, and if you try and give it a 2d vector 
          or array, it complains.... and hence the verbosity of the following 2d definitions **/
      Journal_Printf( stream, "    <Attribute Type=\"Vector\" Center=\"%s\" Name=\"%s\">\n", centering.c_str(),  feVar->name);
      Journal_Printf( stream, "      <DataItem ItemType=\"Function\"  Dimensions=\"%u 3\" Function=\"JOIN($0, $1, 0*$1)\">\n", meshSize );
      Journal_Printf( stream, "        <DataItem ItemType=\"HyperSlab\" Dimensions=\"%u 1\" Name=\"XValue\">\n", meshSize );
      Journal_Printf( stream, "          <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 %u 1 1 %u 1 </DataItem>\n", offset, meshSize );
      Journal_Printf( stream, "          <DataItem Format=\"HDF\" %s Dimensions=\"%u %u\">%s.%05d.h5:/data</DataItem>\n", variableType.c_str(), meshSize, (offset + dofAtEachNodeCount), feVar->name, timestep);
      Journal_Printf( stream, "        </DataItem>\n" );
      Journal_Printf( stream, "        <DataItem ItemType=\"HyperSlab\" Dimensions=\"%u 1\" Name=\"YValue\">\n", meshSize );
      Journal_Printf( stream, "          <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 %u 1 1 %u 1 </DataItem>\n", (offset+1), meshSize );
      Journal_Printf( stream, "          <DataItem Format=\"HDF\" %s Dimensions=\"%u %u\">%s.%05d.h5:/data</DataItem>\n", variableType.c_str(), meshSize, (offset + dofAtEachNodeCount), feVar->name, timestep);
      Journal_Printf( stream, "        </DataItem>\n" );
      Journal_Printf( stream, "      </DataItem>\n" );
      Journal_Printf( stream, "    </Attribute>\n\n" );
    } else if (dofAtEachNodeCount == 3) {
      /** in 3d we simply feed back the 3d hdf5 array, nice and easy **/
      Journal_Printf( stream, "         <Attribute Type=\"Vector\" Center=\"%s\" Name=\"%s\">\n", centering.c_str(),  feVar->name);
      Journal_Printf( stream, "            <DataItem ItemType=\"HyperSlab\" Dimensions=\"%u 3\" >\n", meshSize );
      Journal_Printf( stream, "               <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 %u 1 1 %u 3 </DataItem>\n", offset, meshSize );
      Journal_Printf( stream, "               <DataItem Format=\"HDF\" %s Dimensions=\"%u %u\">%s.%05d.h5:/data</DataItem>\n", variableType.c_str(), meshSize, (offset + dofAtEachNodeCount), feVar->name, timestep);
      Journal_Printf( stream, "            </DataItem>\n" );
      Journal_Printf( stream, "         </Attribute>\n\n" );
    } else {
      /** where there are more than 3 components, we write each one out as a scalar **/
      for(unsigned int dofCountIndex = 0 ; dofCountIndex < dofAtEachNodeCount ; ++dofCountIndex){
        Journal_Printf( stream, "         <Attribute Type=\"Scalar\" Center=\"%s\" Name=\"%s-Component-%u\">\n", centering.c_str(),  feVar->name, dofCountIndex);
        Journal_Printf( stream, "            <DataItem ItemType=\"HyperSlab\" Dimensions=\"%u 1\" >\n", meshSize );
        Journal_Printf( stream, "               <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 %u 1 1 %u 1 </DataItem>\n", (offset+dofCountIndex), meshSize );
        Journal_Printf( stream, "               <DataItem Format=\"HDF\" %s Dimensions=\"%u %u\">%s.%05d.h5:/data</DataItem>\n", variableType.c_str(), meshSize, (offset + dofAtEachNodeCount), feVar->name, timestep);
        Journal_Printf( stream, "            </DataItem>\n" );
        Journal_Printf( stream, "         </Attribute>\n" );
      }
      Journal_Printf( stream, "\n" );
    }
  }

  void write_discontinuous_geometry(Stream *stream, Mesh *mesh,
                                    const int &timestep,
                                    const int &nDims)
  {
    if(nDims==2)
      {
        Grid** grid=
          (Grid** )ExtensionManager_Get
          (mesh->info, mesh, 
           ExtensionManager_GetHandle(mesh->info, (Name)"elementGrid"));

        unsigned *sizes=Grid_GetSizes(*grid);
        unsigned num_elements(sizes[0]*sizes[1]);
        unsigned vert_sizes[]={2*sizes[0]+1,2*sizes[1]+1};
        unsigned num_vertices=vert_sizes[0] * vert_sizes[1];

        char
          *vert_mesh_name(((InnerGenerator*)(mesh->generator))->elMesh->name);

        /* We need a number that goes from 0 to num_elements.  We hack
           it with a WHERE clause that is always true */
        Journal_Printf(stream,"    <DataItem ItemType=\"Function\" Dimensions=\"%u 1\" Function=\"WHERE( $0 != -1)\" Name=\"N\">\n",
                       num_elements);
        Journal_Printf(stream,"      <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"%u 3\">Mesh.%s.%05d.h5:/connectivity</DataItem>\n",
                       num_elements,mesh->name,timestep);
        Journal_Printf(stream,"    </DataItem>\n");

        /* Write out the Quadrilaterals connectivity.  There are 4 quads per element. */
        Journal_Printf(stream,"    <Topology Type=\"Quadrilateral\" NumberOfElements=\"%u\">\n",
                       num_elements*4);
        Journal_Printf(stream,"      <DataItem ItemType=\"Function\" Dimensions=\"%u 4\" Function=\"($0; $1; $2; $3)\">\n",
                       num_elements*4);
        /* Quad 0 */
        Journal_Printf(stream,"        <DataItem ItemType=\"Function\" Dimensions=\"%u 4\" Function=\"(9*$0), (1 + 9*$0), (4 + 9*$0), (3 + 9*$0)\">\n",
                       num_elements);
        Journal_Printf(stream,"          <DataItem Reference=\"XML\">/Xdmf/Domain/Grid[@Name=\"FEM_Grid_%s\"]/DataItem[@Name=\"N\"] </DataItem>\n",
                       mesh->name);
        Journal_Printf(stream,"        </DataItem>\n");
        /* Quad 1 */
        Journal_Printf(stream,"        <DataItem ItemType=\"Function\" Dimensions=\"%u 4\" Function=\"(1 + 9*$0), (2 + 9*$0), (5 + 9*$0), (4 + 9*$0)\">\n",
                       num_elements);
        Journal_Printf(stream,"          <DataItem Reference=\"XML\">/Xdmf/Domain/Grid[@Name=\"FEM_Grid_%s\"]/DataItem[@Name=\"N\"] </DataItem>\n",
                       mesh->name);
        Journal_Printf(stream,"        </DataItem>\n");
        /* Quad 2 */
        Journal_Printf(stream,"        <DataItem ItemType=\"Function\" Dimensions=\"%u 4\" Function=\"(3 + 9*$0), (4 + 9*$0), (7 + 9*$0), (6 + 9*$0)\">\n",
                       num_elements);
        Journal_Printf(stream,"          <DataItem Reference=\"XML\">/Xdmf/Domain/Grid[@Name=\"FEM_Grid_%s\"]/DataItem[@Name=\"N\"] </DataItem>\n",
                       mesh->name);
        Journal_Printf(stream,"        </DataItem>\n");
        /* Quad 3 */
        Journal_Printf(stream,"        <DataItem ItemType=\"Function\" Dimensions=\"%u 4\" Function=\"(4 + 9*$0), (5 + 9*$0), (8 + 9*$0), (7 + 9*$0)\">\n",
                       num_elements);
        Journal_Printf(stream,"          <DataItem Reference=\"XML\">/Xdmf/Domain/Grid[@Name=\"FEM_Grid_%s\"]/DataItem[@Name=\"N\"] </DataItem>\n",
                       mesh->name);
        Journal_Printf(stream,"        </DataItem>\n");
        Journal_Printf(stream,"      </DataItem>\n");
        Journal_Printf(stream,"    </Topology>\n");


        /* Write out the coordinates.  This is a bit complicated,
           because we want separate vertices for each element, so the
           sides of elements can not share.  */

        /* First get the coordinates from the vertex mesh */
        Journal_Printf(stream,"    <DataItem ItemType=\"HyperSlab\" Dimensions=\"%u %u\" Name=\"XCoords\">\n",
                       sizes[0]*2+1,sizes[1]*2+1);
        Journal_Printf(stream,"      <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 0 1 1 %u 1 </DataItem>\n",
                       num_vertices);
        Journal_Printf(stream,"      <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"%u 2\">Mesh.%s.%05d.h5:/vertices</DataItem>\n",
                       num_vertices,vert_mesh_name,timestep);
        Journal_Printf(stream,"    </DataItem>\n");

        Journal_Printf(stream,"    <DataItem ItemType=\"HyperSlab\" Dimensions=\"%u %u\" Name=\"YCoords\">\n",
                       sizes[0]*2+1,sizes[1]*2+1);
        Journal_Printf(stream,"      <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 1 1 1 %u 1 </DataItem>\n",
                       num_vertices);
        Journal_Printf(stream,"      <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"%u 2\">Mesh.%s.%05d.h5:/vertices</DataItem>\n",
                       num_vertices,vert_mesh_name,timestep);
        Journal_Printf(stream,"    </DataItem>\n");

        std::string direction[]={"X", "Y"};
        for(int d=0;d<2;++d)
          {
            for(int n=0;n<9;++n)
              {
                Journal_Printf(stream,"    <DataItem ItemType=\"HyperSlab\" Dimensions=\"%u 1\" Name=\"%s%u\">\n",
                               sizes[0]*sizes[1],direction[d].c_str(),n);
                Journal_Printf(stream,"      <DataItem Dimensions=\"3 2\" Format=\"XML\"> %u %u 2 2 %u %u </DataItem>\n",
                               n/3,n%3,sizes[0],sizes[1]);
                Journal_Printf(stream,"      <DataItem Reference=\"XML\">/Xdmf/Domain/Grid[@Name=\"FEM_Grid_%s\"]/DataItem[@Name=\"%sCoords\"] </DataItem>\n",
                               mesh->name,direction[d].c_str());
                Journal_Printf(stream,"    </DataItem>\n");
              }
            Journal_Printf(stream,"    <DataItem ItemType=\"Function\" Dimensions=\"%u 1\" Function=\"$0, $1, $2, $3, $4, $5, $6, $7, $8\" Name=\"%s\">\n",
                           9*num_elements,direction[d].c_str());
            for(int i=0;i<9;++i)
              {
                Journal_Printf(stream,"      <DataItem Reference=\"XML\">/Xdmf/Domain/Grid[@Name=\"FEM_Grid_%s\"]/DataItem[@Name=\"%s%d\"] </DataItem>\n",
                               mesh->name,direction[d].c_str(),i);
              }
            Journal_Printf(stream,"    </DataItem>\n");
          }
        Journal_Printf(stream,"    <Geometry Type=\"XYZ\">\n");
        Journal_Printf(stream,"      <DataItem ItemType=\"Function\" Dimensions=\"%u 3\" Function=\"JOIN($0, $1, 0*$1)\">\n",
                       9*num_elements);
        Journal_Printf(stream,"        <DataItem Reference=\"XML\">/Xdmf/Domain/Grid[@Name=\"FEM_Grid_%s\"]/DataItem[@Name=\"X\"] </DataItem>\n",
                       mesh->name);
        Journal_Printf(stream,"        <DataItem Reference=\"XML\">/Xdmf/Domain/Grid[@Name=\"FEM_Grid_%s\"]/DataItem[@Name=\"Y\"] </DataItem>\n",
                       mesh->name);
        Journal_Printf(stream,"      </DataItem>\n");
        Journal_Printf(stream,"    </Geometry>\n");
      }
    /* nDim==3 */
    else
      {
        abort();
      }
  }

  void write_three_component_function(Stream *stream, 
                                      const std::string &start,
                                      const std::string &end,
                                      const int &num_elements,
                                      const std::string &equation,
                                      const std::string &prefix,
                                      const std::string &result,
                                      const std::string &first,
                                      const std::string &second,
                                      const std::string &third)
  {
    std::stringstream ss;
    ss << "    <DataItem ItemType=\"Function\"  Dimensions=\""
       << num_elements
       << " 1\" Function=\""
       << equation
       << "\" Name=\""
       << prefix
       << result
       << "\">\n";
    Journal_Printf(stream, ss.str().c_str());

    Journal_Printf(stream,(start + first + end).c_str());
    Journal_Printf(stream,(start + second + end).c_str());
    Journal_Printf(stream,(start + third + end).c_str());

    Journal_Printf(stream,"    </DataItem>\n");
  }

  void write_discontinuous_FeVariable(Stream *stream,
                                      const std::string &mesh_name,
                                      const std::string &var_name,
                                      const int &timestep,
                                      const int &num_elements)
  {
    for(int i=0;i<3;++i)
      {
        Journal_Printf(stream,"    <DataItem ItemType=\"HyperSlab\" Dimensions=\"%u 1\" Name=\"%s%d\">\n",
                       num_elements,var_name.c_str(),i);
        Journal_Printf(stream,"      <DataItem Dimensions=\"3 2\" Format=\"XML\"> %u 0 3 1 %u 1 </DataItem>\n",
                       i,num_elements);
        Journal_Printf(stream,"      <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"%u 1\">%s.%05d.h5:/data</DataItem>\n",
                       3*num_elements,var_name.c_str(),timestep);
        Journal_Printf(stream,"    </DataItem>\n");
      }

    /* Write P_avg, dP/dx, dP/dy */
    std::string start("      <DataItem Reference=\"XML\">/Xdmf/Domain/Grid[@Name=\"FEM_Grid_"
                      + mesh_name
                      + "\"]/DataItem[@Name=\""
                      + var_name);
    std::string end("\"] </DataItem>\n");
    write_three_component_function(stream,start,end,num_elements,
                                   "($2 + ($1 + $0) / 2) / 2",
                                   var_name,"_avg","0","1","2");
    write_three_component_function(stream,start,end,num_elements,
                                   "($1 - $0)",var_name,"_dx","0","1","2");
    write_three_component_function(stream,start,end,num_elements,
                                   "($2 - ($1 + $0) / 2)",
                                   var_name,"_dy","0","1","2");

    /* Write P at each vertex */
    write_three_component_function(stream,start,end,num_elements,
                                   "($0 - $1 - $2)",var_name,
                                   "_0","_avg","_dx","_dy");
    write_three_component_function(stream,start,end,num_elements,
                                   "$0 - $2",var_name,"_1","_avg","_dx","_dy");
    write_three_component_function(stream,start,end,num_elements,
                                   "$0 + $1 - $2",var_name,
                                   "_2","_avg","_dx","_dy");
    write_three_component_function(stream,start,end,num_elements,
                                   "$0 - $1",var_name,"_3","_avg","_dx","_dy");
    write_three_component_function(stream,start,end,num_elements,
                                   "$0",var_name,"_4","_avg","_dx","_dy");
    write_three_component_function(stream,start,end,num_elements,
                                   "$0 + $1",var_name,"_5","_avg","_dx","_dy");
    write_three_component_function(stream,start,end,num_elements,
                                   "$0 - $1 + $2",var_name,
                                   "_6","_avg","_dx","_dy");
    write_three_component_function(stream,start,end,num_elements,
                                   "$0 + $2",var_name,"_7","_avg","_dx","_dy");
    write_three_component_function(stream,start,end,num_elements,
                                   "$0 + $1 + $2",var_name,
                                   "_8","_avg","_dx","_dy");

    /* Bundle them all together for the field */
    Journal_Printf(stream,"    <Attribute Type=\"Scalar\" Center=\"Node\" Name=\"%s\">\n",
                   var_name.c_str());
    Journal_Printf(stream,"      <DataItem ItemType=\"Function\"  Dimensions=\"%u 1\" Function=\"$0, $1, $2, $3, $4, $5, $6, $7, $8\">\n",
                   num_elements*9);
    for(int i=0;i<9;++i)
      {
        std::stringstream ss;
        ss << "  "
           << start
           << "_"
           << i
           << end;
        Journal_Printf(stream,ss.str().c_str());
      }
    Journal_Printf(stream,"      </DataItem>\n");
    Journal_Printf(stream,"    </Attribute>\n");
  }
}
