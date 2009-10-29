/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
** Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
** Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
** Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
** Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
** Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
** Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
** $Id: FileParticleLayout.c 4081 2007-04-27 06:20:07Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef READ_HDF5
#include <hdf5.h>
#endif

#include <mpi.h>
#include <StGermain/StGermain.h>

#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>
#include <StgDomain/Utils/Utils.h>

#include "types.h"
#include "shortcuts.h"
#include "ParticleLayout.h"
#include "GlobalParticleLayout.h"
#include "FileParticleLayout.h"

#include "SwarmClass.h"
#include "StandardParticle.h"
#include "ShadowInfo.h"
#include "CellLayout.h"
#include "ElementCellLayout.h"
#include "IntegrationPoint.h"
#include "SwarmVariable.h"
#include "SwarmVariable_Register.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

const Type FileParticleLayout_Type = "FileParticleLayout";

FileParticleLayout* FileParticleLayout_New( Name name, Name filename, Index checkpointnfiles )
{
   FileParticleLayout* self = (FileParticleLayout*) _FileParticleLayout_DefaultNew( name );
   _FileParticleLayout_Init( self, filename, checkpointnfiles );
   return self;
}


FileParticleLayout* _FileParticleLayout_New( 
                SizeT                                            _sizeOfSelf,
                Type                                             type,
                Stg_Class_DeleteFunction*                        _delete,
                Stg_Class_PrintFunction*                         _print,
                Stg_Class_CopyFunction*                          _copy,
                Stg_Component_DefaultConstructorFunction*        _defaultConstructor,
                Stg_Component_ConstructFunction*                 _construct,
                Stg_Component_BuildFunction*                     _build,
                Stg_Component_InitialiseFunction*                _initialise,
                Stg_Component_ExecuteFunction*                   _execute,
                Stg_Component_DestroyFunction*                   _destroy,
                ParticleLayout_SetInitialCountsFunction*         _setInitialCounts,
                ParticleLayout_InitialiseParticlesFunction*      _initialiseParticles,
                GlobalParticleLayout_InitialiseParticleFunction* _initialiseParticle,
                Name                                             name,
                Bool                                             initFlag,
                Index                                            checkpointnfiles,
                Name                                             filename )
{
   FileParticleLayout* self;
   
   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof( FileParticleLayout ) );
   self = (FileParticleLayout*)_GlobalParticleLayout_New( 
         _sizeOfSelf, 
         type,
         _delete,
         _print,
         _copy, 
         _defaultConstructor,
         _construct,
         _build,
         _initialise,
         _execute,
         _destroy,
         _setInitialCounts,
         _initialiseParticles,
         _initialiseParticle,
         name,
         initFlag,
         GlobalCoordSystem,
         False,
         0,
         0.0 );

   if ( initFlag ) {
      _FileParticleLayout_Init( self, filename, checkpointnfiles );
   }

   return self;
}


void _FileParticleLayout_Init( void* particleLayout, Name filename, Index checkpointnfiles )
{
   FileParticleLayout* self = (FileParticleLayout*) particleLayout;

   self->filename = StG_Strdup( filename );
   self->checkpointnfiles = checkpointnfiles;
   self->errorStream = Journal_MyStream( Error_Type, self );
   _GlobalParticleLayout_Init( self, GlobalCoordSystem, False, 0, 0.0 );
}


void _FileParticleLayout_Delete( void* particleLayout ) {
   FileParticleLayout* self = (FileParticleLayout*)particleLayout;

   Memory_Free( self->filename );

   /* Stg_Class_Delete parent class */
   _GlobalParticleLayout_Delete( self );
}

void _FileParticleLayout_Print( void* particleLayout, Stream* stream ) {
   FileParticleLayout* self = (FileParticleLayout*)particleLayout;
   
   /* General info */
   Journal_Printf( stream, "FileParticleLayout (ptr): %p:\n", self );
   Stream_Indent( stream );
   
   /* Parent class info */
   _GlobalParticleLayout_Print( self, stream );
   
   /* FileParticleLayout */
   Journal_Printf( stream, "filename: %s\n", self->filename );
   
   Stream_UnIndent( stream );
}


void* _FileParticleLayout_Copy( void* particleLayout, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
   FileParticleLayout*     self                    = (FileParticleLayout*)particleLayout;
   FileParticleLayout*     newFileParticleLayout;
   
   newFileParticleLayout = (FileParticleLayout*)_GlobalParticleLayout_Copy( self, dest, deep, nameExt, ptrMap );
   
   newFileParticleLayout->filename = self->filename;

   return (void*)newFileParticleLayout;
}

void* _FileParticleLayout_DefaultNew( Name name ) {
   return (void*)_FileParticleLayout_New( 
         sizeof(FileParticleLayout),
         FileParticleLayout_Type,
         _FileParticleLayout_Delete,
         _FileParticleLayout_Print,
         _FileParticleLayout_Copy,
         _FileParticleLayout_DefaultNew,
         _FileParticleLayout_AssignFromXML,
         _FileParticleLayout_Build,
         _FileParticleLayout_Initialise,
         _FileParticleLayout_Execute,
         _FileParticleLayout_Destroy,
         _FileParticleLayout_SetInitialCounts,
         _FileParticleLayout_InitialiseParticles,
         _FileParticleLayout_InitialiseParticle,
         name,
         False,
         0,/* checkpointnfiles*/
         NULL /* filename */ );
}

void _FileParticleLayout_AssignFromXML( void* particleLayout, Stg_ComponentFactory *cf, void* data ) {
   FileParticleLayout* self     = (FileParticleLayout*) particleLayout;
   Name                filename;
   Index               checkpointnfiles;

   self->context = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Context", AbstractContext, False, data );
   if( !self->context )
      self->context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data );

   filename         = Stg_ComponentFactory_GetString( cf, self->name, "filename", "Swarm" );
   
   checkpointnfiles = Stg_ComponentFactory_GetInt( cf, self->name, "checkpointnfiles", 1 );
#ifdef READ_HDF5
   /* if doing checkpoint restart, grab number of particles swarmdump previously stored against */
   if( self->context->loadFromCheckPoint ) checkpointnfiles = _FileParticleLayout_GetFileCountFromTimeInfoFile( self->context );
   Journal_Firewall( checkpointnfiles > 0, self->errorStream,
              "Error in %s for %s '%s' - determined number of fileParticleLayout checkpoint files (%d) for reload is not valid.\n",
              __func__,
              self->type,
              self->name,
              checkpointnfiles );
#endif

   _FileParticleLayout_Init( self, filename, checkpointnfiles );
}
   
void _FileParticleLayout_Build( void* particleLayout, void* data ) {
}
void _FileParticleLayout_Initialise( void* particleLayout, void* data ) {
}
void _FileParticleLayout_Execute( void* particleLayout, void* data ) {
}
void _FileParticleLayout_Destroy( void* particleLayout, void* data ) {
}

void _FileParticleLayout_SetInitialCounts( void* particleLayout, void* _swarm ) {
   FileParticleLayout*  self         = (FileParticleLayout*)particleLayout;
   Swarm*               swarm        = (Swarm*)_swarm;
   Name                 filename     = self->filename;
   
#ifdef READ_HDF5
   hid_t                file, fileData, fileSpace;
   hid_t                group_id, attrib_id;
   hsize_t              size[2];
   char                 dataSpaceName[1024];
   SwarmVariable*       swarmVar;
   Index                ii;
   Index                swarmVar_I;
   int                  nParticles;
   herr_t               status;
#else
   MPI_File             mpiFile;
   int                  openResult;
   MPI_Offset           bytesCount;
   SizeT                particleSize = swarm->particleExtensionMgr->finalSize;
   ldiv_t               division;
#endif

   Journal_DPrintf( self->debug, "In %s(): for ParticleLayout \"%s\", of type %s\n",
      __func__, self->name, self->type );
   Stream_IndentBranch( Swarm_Debug ); 

#ifdef READ_HDF5
   self->lastParticleIndex  = Memory_Alloc_Array( Index, self->checkpointnfiles, "lastParticleIndex" );
   self->totalInitialParticles = 0;
   for( ii = 1 ; ii <= self->checkpointnfiles ; ii++ ){
      char* filenameTemp = NULL;
      /* Open the swarm checkpointing file */
      if(self->checkpointnfiles == 1)
         Stg_asprintf( &filenameTemp, "%s.h5", filename );
      else 
         Stg_asprintf( &filenameTemp, "%s.%dof%d.h5", filename, ii, self->checkpointnfiles );
      
      file = H5Fopen( filenameTemp, H5F_ACC_RDONLY, H5P_DEFAULT );
      Journal_Firewall( file >= 0,
              self->errorStream,
              "Error in %s for %s '%s' - Cannot open file %s.\n",
              __func__,
              self->type,
              self->name,
              filenameTemp );
      
      /* get the file attributes to determine if this file contains particles */
      #if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
         group_id  = H5Gopen(file, "/");
         attrib_id = H5Aopen_name(group_id, "Swarm Particle Count");
      #else
         group_id  = H5Gopen(file, "/", H5P_DEFAULT);
         attrib_id = H5Aopen(group_id, "Swarm Particle Count", H5P_DEFAULT);
      #endif
       Journal_Firewall( attrib_id > 0,
              self->errorStream,
              "\nError in %s for %s '%s' - Swarm Particle Count group not present in checkpoint file %s.\n  Perhaps you are trying to restart from an old checkpoint file.",
              __func__,
              self->type,
              self->name,
              filenameTemp );
      status = H5Aread(attrib_id, H5T_NATIVE_INT, &nParticles);
      H5Aclose(attrib_id);
      H5Gclose(group_id);

      self->totalInitialParticles += nParticles;
      self->lastParticleIndex[ii-1] = self->totalInitialParticles;
      
      /* Close the dataspace and file */
      H5Fclose( file );
      Memory_Free( filenameTemp );
   }
#else
   Journal_DPrintf( self->debug, "Finding number of bytes in checkpoint file \"%s\":\n",
      self->filename );
      
   openResult = MPI_File_open( swarm->comm, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpiFile );
   
   Journal_Firewall( 
      openResult == 0, 
      self->errorStream,
      "Error in %s for %s '%s' - Cannot open file %s.\n", 
      __func__, 
      self->type, 
      self->name, 
      filename );
   
   MPI_File_get_size( mpiFile, &bytesCount );
   MPI_File_close( &mpiFile );

   Journal_DPrintf( self->debug, "...calculated bytes total of %u.\n", bytesCount );
   
   /* Divide by particle size to get number of particles */
   division = ldiv( bytesCount, (long) particleSize );
   self->totalInitialParticles = (unsigned int) division.quot;
   
   Journal_DPrintf( self->debug, "given bytes total %u / particle size %u ->\n"
      "\ttotalInitialParticles = %u.\n", bytesCount, (unsigned int)particleSize,
      self->totalInitialParticles );

   Journal_Firewall( 
      division.rem == 0,
      self->errorStream,
      "Error in func %s for %s '%s' - Trying to read particle information from %s which stores %u bytes.\n"
      "This doesn't produce an integer number of particles of size %u - It gives remainder = %u\n", 
      __func__, 
      self->type, 
      self->name, 
      filename, 
      bytesCount, 
      (unsigned int)particleSize, 
      division.rem ); 

   Journal_DPrintf( self->debug, "calling parent func to set cell counts:\n", bytesCount );
#endif

   _GlobalParticleLayout_SetInitialCounts( self, swarm );

   Stream_UnIndentBranch( Swarm_Debug );  
   Journal_DPrintf( self->debug, "...finished %s() for ParticleLayout \"%s\".\n",
      __func__, self->name );
}

void _FileParticleLayout_InitialiseParticles( void* particleLayout, void* _swarm ) {
   FileParticleLayout*    self  = (FileParticleLayout*)particleLayout;
   Swarm                 *swarm = (Swarm*)_swarm;
   
#ifdef READ_HDF5
   SwarmVariable*         swarmVar;
   Index                  swarmVar_I;
   char                   dataSpaceName[1024];
   hid_t                  file[self->checkpointnfiles];
   Index                  ii, jj, kk;
   hid_t                  group_id, attrib_id;
   int                    nParticles;
   herr_t                 status;
     
   /* Allocate space to store arrays of dataspaces */   
   assert( swarm->swarmVariable_Register );  
   self->fileData  = Memory_Alloc_2DArray( hid_t, swarm->swarmVariable_Register->objects->count, self->checkpointnfiles, "fileData" );
   self->fileSpace = Memory_Alloc_2DArray( hid_t, swarm->swarmVariable_Register->objects->count, self->checkpointnfiles, "fileSpace" );
   /* set these spaces to null initially */
   for( jj = 0 ; jj < swarm->swarmVariable_Register->objects->count ; jj++)
      for( kk = 0 ; kk < self->checkpointnfiles ; kk++){
         self->fileData [jj][kk] = NULL;
         self->fileSpace[jj][kk] = NULL;
      }
      
   /* Open the files */
   for( ii = 1 ; ii <= self->checkpointnfiles ; ii++ ){
      char*  filenameTemp = NULL;
      /* Open the swarm checkpointing file */
      if(self->checkpointnfiles == 1)
         Stg_asprintf( &filenameTemp, "%s.h5", self->filename );
      else 
         Stg_asprintf( &filenameTemp, "%s.%dof%d.h5", self->filename, ii, self->checkpointnfiles );

      file[ii-1] = H5Fopen( filenameTemp, H5F_ACC_RDONLY, H5P_DEFAULT );
      Journal_Firewall( 
            file[ii-1] >= 0, 
            self->errorStream,
            "Error in %s for %s '%s' - Cannot open file %s.\n", 
            __func__, 
            self->type, 
            self->name, 
            self->filename );

      /* get the file attributes to determine if this file contains particles */
      #if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
         group_id  = H5Gopen(file[ii-1], "/");
         attrib_id = H5Aopen_name(group_id, "Swarm Particle Count");
      #else
         group_id  = H5Gopen(file[ii-1], "/", H5P_DEFAULT);
         attrib_id = H5Aopen(group_id, "Swarm Particle Count", H5P_DEFAULT);
      #endif
      status = H5Aread(attrib_id, H5T_NATIVE_INT, &nParticles);

      H5Aclose(attrib_id);
      H5Gclose(group_id);

      if(nParticles > 0){
         /* Open a dataspace for each swarmVariable */
         for( swarmVar_I = 0; swarmVar_I < swarm->swarmVariable_Register->objects->count; swarmVar_I++ ) {
            swarmVar = SwarmVariable_Register_GetByIndex( swarm->swarmVariable_Register, swarmVar_I );
            
            if( swarmVar->isCheckpointedAndReloaded ) {
               sprintf( dataSpaceName, "/%s", swarmVar->name + strlen(swarm->name)+1 );
               
               #if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
                  self->fileData[swarmVar_I][ii-1]  = H5Dopen( file[ii-1], dataSpaceName );
               #else
                  self->fileData[swarmVar_I][ii-1]  = H5Dopen( file[ii-1], dataSpaceName, H5P_DEFAULT );
               #endif
                  self->fileSpace[swarmVar_I][ii-1] = H5Dget_space( self->fileData[swarmVar_I][ii-1] );
               
               Variable_Update( swarmVar->variable );
            }
         }
      }
      Memory_Free( filenameTemp );
   }
       
   self->start[1] = 0;
   self->count[0] = 1; 
   
   _GlobalParticleLayout_InitialiseParticles( self, _swarm );

   /* Close dataspaces and the file */
   for( ii = 1 ; ii <= self->checkpointnfiles ; ii++ ){
      for( swarmVar_I = 0; swarmVar_I < swarm->swarmVariable_Register->objects->count; swarmVar_I++ ) {
         swarmVar = SwarmVariable_Register_GetByIndex( swarm->swarmVariable_Register, swarmVar_I );
         if( swarmVar->isCheckpointedAndReloaded ) {
            if ( self->fileSpace[swarmVar_I][ii-1] ) H5Sclose( self->fileSpace[swarmVar_I][ii-1] );
            if ( self->fileData [swarmVar_I][ii-1] ) H5Dclose( self->fileData [swarmVar_I][ii-1] );
         }
      }
      H5Fclose( file[ii-1] );
   }
   
   Memory_Free( self->fileData );
   Memory_Free( self->fileSpace );
   Memory_Free( self->lastParticleIndex );
   
#else
   self->file = fopen( self->filename, "rb" );
   Journal_Firewall( 
      self->file != NULL, 
      self->errorStream,
      "Error in %s for %s '%s' - Cannot open file %s.\n", 
      __func__, 
      self->type, 
      self->name, 
      self->filename );

   _GlobalParticleLayout_InitialiseParticles( self, _swarm );
   
   fclose( self->file );
#endif
}  
	
void _FileParticleLayout_InitialiseParticle( 
      void*              particleLayout, 
      void*              _swarm, 
      Particle_Index     newParticle_I,
      void*              particle )
{
   FileParticleLayout*     self         = (FileParticleLayout*)particleLayout;
   Swarm*                  swarm        = (Swarm*)_swarm;
   SizeT                   particleSize = swarm->particleExtensionMgr->finalSize;
   int                     result;
   IntegrationPoint*       newParticle  = (IntegrationPoint*)particle;

#ifdef READ_HDF5
   SwarmVariable*    swarmVar;
   Index             swarmVar_I;
   Index             ii;
   hid_t             memSpace; 

   /* find out which file particle is contained within */
   for( ii = 1 ; ii <= self->checkpointnfiles ; ii++ ){
      if( newParticle_I < self->lastParticleIndex[ii-1]) break;
   }
   self->start[0] = newParticle_I;
   if(ii != 1) self->start[0] -= self->lastParticleIndex[ii-2];
   
   for( swarmVar_I = 0; swarmVar_I < swarm->nSwarmVars; swarmVar_I++ ) {
      swarmVar = SwarmVariable_Register_GetByIndex( swarm->swarmVariable_Register, swarmVar_I );
      /* only retrieve variable if it does not have a parent, as these
         are not stored, and the data is retrieved when the parent is. 
         Also do make sure variable is of type SwarmVariable (as opposed 
         to materialSwarmVariable, which is not checkpointed). */
      if( swarmVar->isCheckpointedAndReloaded ) {          
         /* Update the hyperslab. */   
         self->count[1] = swarmVar->dofCount;
         memSpace = H5Screate_simple( 2, self->count, NULL );
         H5Sselect_hyperslab( self->fileSpace[swarmVar_I][ii-1], H5S_SELECT_SET, self->start, NULL, self->count, NULL );
         H5Sselect_all( memSpace );
         
         /* Treat the data differently depending on its type */
         if( swarmVar->variable->dataTypes[0] == Variable_DataType_Int ) {
            int* particleInfo = Memory_Alloc_Array( int, swarmVar->dofCount, "particleCheckpointInfo" );
            
            /* Read particle data. */
            H5Dread( self->fileData[swarmVar_I][ii-1], H5T_NATIVE_INT, memSpace, self->fileSpace[swarmVar_I][ii-1], H5P_DEFAULT, particleInfo );
            
            Variable_SetValue( swarmVar->variable, swarm->particleLocalCount, particleInfo );
            
            Memory_Free( particleInfo );
         }
        
         else if( swarmVar->variable->dataTypes[0] == Variable_DataType_Char) {
            char* particleInfo = Memory_Alloc_Array( char, swarmVar->dofCount, "particleCheckpointInfo" );
            
            /* Read particle data. */
            H5Dread( self->fileData[swarmVar_I][ii-1], H5T_NATIVE_CHAR, memSpace, self->fileSpace[swarmVar_I][ii-1], H5P_DEFAULT, particleInfo );
            
            Variable_SetValue( swarmVar->variable, swarm->particleLocalCount, particleInfo );
            
            Memory_Free( particleInfo );
         }
            
         else if( swarmVar->variable->dataTypes[0] == Variable_DataType_Float ) {
            float* particleInfo = Memory_Alloc_Array( float, swarmVar->dofCount, "particleCheckpointInfo" );
            
            /* Read particle data. */
            H5Dread( self->fileData[swarmVar_I][ii-1], H5T_NATIVE_FLOAT, memSpace, self->fileSpace[swarmVar_I][ii-1], H5P_DEFAULT, particleInfo );
            
            Variable_SetValue( swarmVar->variable, swarm->particleLocalCount, particleInfo );
            
            Memory_Free( particleInfo );
         }
            
         else {
            double* particleInfo = Memory_Alloc_Array( double, swarmVar->dofCount, "particleCheckpointInfo" );
            
            /* Read particle data. */
            H5Dread( self->fileData[swarmVar_I][ii-1], H5T_NATIVE_DOUBLE, memSpace, self->fileSpace[swarmVar_I][ii-1], H5P_DEFAULT, particleInfo );
            
            Variable_SetValue( swarmVar->variable, swarm->particleLocalCount, particleInfo );
            
            Memory_Free( particleInfo );
         }   
         
         H5Sclose( memSpace );
      }
   }
      
#else
   result = fread( particle, particleSize, 1, self->file );

   Journal_Firewall( 
      result == 1,
      self->errorStream,
      "Error in func %s for %s '%s':\n"
      "\tCouldn't read in particle %u - May have reached end-of-file.\n",
      __func__, 
      self->type, 
      self->name, 
      newParticle_I );
#endif
}

#ifdef READ_HDF5
Index _FileParticleLayout_GetFileCountFromTimeInfoFile( void* _context ){
   AbstractContext*     context = (AbstractContext*)_context;
   char*                timeInfoFileName = NULL;
   char*                timeInfoFileNamePart = NULL;
   FILE*                timeInfoFile;     
   Stream*              errorStr = Journal_Register( Error_Type, FileParticleLayout_Type );
   hid_t                file, fileSpace, fileData;
   Index                checkpointnproc;
   
   timeInfoFileNamePart = Context_GetCheckPointReadPrefixString( context );
   Stg_asprintf( &timeInfoFileName, "%stimeInfo.%.5u.h5", timeInfoFileNamePart, context->restartTimestep );
    
   /* Open the file and data set. */
   file = H5Fopen( timeInfoFileName, H5F_ACC_RDONLY, H5P_DEFAULT );
   Journal_Firewall( 
      file >= 0, 
      errorStr, "\n\nError- in %s(), Couldn't find checkpoint time info file with "
      "filename \"%s\" - aborting.\n", __func__, timeInfoFileName );
   
   /* Read previous nproc from file */
   #if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
   fileData = H5Dopen( file, "/nproc" );
   #else
   fileData = H5Dopen( file, "/nproc", H5P_DEFAULT );
   #endif
   fileSpace = H5Dget_space( fileData );
      
   H5Dread( fileData, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &checkpointnproc );
      
   H5Sclose( fileSpace );
   H5Dclose( fileData );
   
   H5Fclose( file );
      
   Memory_Free( timeInfoFileName );
   Memory_Free( timeInfoFileNamePart );
   return checkpointnproc;
}
#endif


