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

FileParticleLayout* FileParticleLayout_New( Name name, Name filename )
{
	FileParticleLayout* self = (FileParticleLayout*) _FileParticleLayout_DefaultNew( name );
	_FileParticleLayout_Init( self, filename );
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
		_FileParticleLayout_Init( self, filename );
	}

	return self;
}


void _FileParticleLayout_Init( void* particleLayout, Name filename )
{
	FileParticleLayout* self = (FileParticleLayout*) particleLayout;

	self->filename = StG_Strdup( filename );
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
	FileParticleLayout* self  = (FileParticleLayout*)particleLayout;
	
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
	FileParticleLayout*		self                    = (FileParticleLayout*)particleLayout;
	FileParticleLayout*		newFileParticleLayout;
	
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
			_FileParticleLayout_Construct,
			_FileParticleLayout_Build,
			_FileParticleLayout_Initialise,
			_FileParticleLayout_Execute,
			_FileParticleLayout_Destroy,
			_FileParticleLayout_SetInitialCounts,
			_FileParticleLayout_InitialiseParticles,
			_FileParticleLayout_InitialiseParticle,
			name,
			False, 
			NULL /* filename */ );
}

void _FileParticleLayout_Construct( void* particleLayout, Stg_ComponentFactory *cf, void* data ) {
	FileParticleLayout* self     = (FileParticleLayout*) particleLayout;
	Name                filename;

	filename = Stg_ComponentFactory_GetString( cf, self->name, "filename", "Swarm.dat" );
	
	_FileParticleLayout_Init( self, filename );
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
	hsize_t              size[2];
	char                 dataSpaceName[1024];
	SwarmVariable*       swarmVar;
#else
	MPI_File             mpiFile;
	int                  openResult;
	MPI_Offset           bytesCount;
	SizeT                particleSize = swarm->particleExtensionMgr->finalSize;
	ldiv_t                division;
#endif

	Journal_DPrintf( self->debug, "In %s(): for ParticleLayout \"%s\", of type %s\n",
		__func__, self->name, self->type );
	Stream_IndentBranch( Swarm_Debug );	

#ifdef READ_HDF5
   /* Open the swarm checkpointing file */
	file = H5Fopen( filename, H5F_ACC_RDONLY, H5P_DEFAULT );
	Journal_Firewall( 
		file >= 0, 
		self->errorStream,
		"Error in %s for %s '%s' - Cannot open file %s.\n", 
		__func__, 
		self->type, 
		self->name, 
		filename );

   /* Open a dataspace */
   swarmVar = SwarmVariable_Register_GetByIndex( swarm->swarmVariable_Register, 0 );
   sprintf( dataSpaceName, "/%s", swarmVar->name );
      
   #if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
	fileData = H5Dopen( file, dataSpaceName );
   #else
	fileData = H5Dopen( file, dataSpaceName, H5P_DEFAULT );
   #endif
	fileSpace = H5Dget_space( fileData );
	
	/* Get the dimensions of the open dataspace */
   H5Sget_simple_extent_dims( fileSpace, size, NULL ); 
   
   self->totalInitialParticles = size[0];
   
   /* Close the dataspace and file */
   H5Sclose( fileSpace );
	H5Dclose( fileData );	 
	H5Fclose( file );	
   
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
	FileParticleLayout*        self             = (FileParticleLayout*)particleLayout;
	Swarm *swarm = (Swarm*)_swarm;
	
#ifdef READ_HDF5
	SwarmVariable*          swarmVar;
   Index                   swarmVar_I, dof_I;
   char                    dataSpaceName[1024];
   hid_t                   file;
     
   /* Allocate space to store arrays of dataspaces */   
   assert( swarm->swarmVariable_Register );  
   self->fileData = Memory_Alloc_Array( hid_t, swarm->swarmVariable_Register->objects->count, "fileData" );
   self->fileSpace = Memory_Alloc_Array( hid_t, swarm->swarmVariable_Register->objects->count, "fileSpace" );
	 
	/* Open the file */
	file = H5Fopen( self->filename, H5F_ACC_RDONLY, H5P_DEFAULT );
	Journal_Firewall( 
		file >= 0, 
		self->errorStream,
		"Error in %s for %s '%s' - Cannot open file %s.\n", 
		__func__, 
		self->type, 
		self->name, 
		self->filename );

   /* Open a dataspace for each swarmVariable */
   for( swarmVar_I = 0; swarmVar_I < swarm->swarmVariable_Register->objects->count; swarmVar_I++ ) {
      swarmVar = SwarmVariable_Register_GetByIndex( swarm->swarmVariable_Register, swarmVar_I );
      sprintf( dataSpaceName, "/%s", swarmVar->name );
      
      #if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
	   self->fileData[swarmVar_I] = H5Dopen( file, dataSpaceName );
      #else
	   self->fileData[swarmVar_I] = H5Dopen( file, dataSpaceName, H5P_DEFAULT );
      #endif
	   self->fileSpace[swarmVar_I] = H5Dget_space( self->fileData[swarmVar_I] );
   
      Variable_Update( swarmVar->variable );
   }
       
	self->start[1] = 0;
	self->count[0] = 1; 
	
	_GlobalParticleLayout_InitialiseParticles( self, _swarm );

   /* Close dataspaces and the file */
   for( swarmVar_I = 0; swarmVar_I < swarm->swarmVariable_Register->objects->count; swarmVar_I++ ) {
	   H5Sclose( self->fileSpace[swarmVar_I] );
	   H5Dclose( self->fileData[swarmVar_I] );
	}
	H5Fclose( file );	
	
	Memory_Free( self->fileData );
	Memory_Free( self->fileSpace );
	
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
	FileParticleLayout*        self                 = (FileParticleLayout*)particleLayout;
	Swarm*                     swarm                = (Swarm*)_swarm;
	SizeT                      particleSize         = swarm->particleExtensionMgr->finalSize;
	int                        result;
	IntegrationPoint*          newParticle          = (IntegrationPoint*)particle;

#ifdef READ_HDF5
   SwarmVariable*          swarmVar;
   Index                   swarmVar_I;
   hid_t                   memSpace; 
    
	self->start[0] = newParticle_I;
	
	for( swarmVar_I = 0; swarmVar_I < swarm->nSwarmVars; swarmVar_I++ ) {
      swarmVar = SwarmVariable_Register_GetByIndex( swarm->swarmVariable_Register, swarmVar_I );
       
      /* Update the hyperslab. */   
      self->count[1] = swarmVar->dofCount;
	   memSpace = H5Screate_simple( 2, self->count, NULL );
	   H5Sselect_hyperslab( self->fileSpace[swarmVar_I], H5S_SELECT_SET, self->start, NULL, self->count, NULL );
      H5Sselect_all( memSpace );
      
      /* Treat the data differently depending on its type */
      if( swarmVar->variable->dataTypes[0] == Variable_DataType_Int ) {
            int* particleInfo = Memory_Alloc_Array( int, swarmVar->dofCount, "particleCheckpointInfo" );
            
            /* Read particle data. */
	         H5Dread( self->fileData[swarmVar_I], H5T_NATIVE_INT, memSpace, self->fileSpace[swarmVar_I], H5P_DEFAULT, particleInfo );
	         
	         Variable_SetValue( swarmVar->variable, swarm->particleLocalCount, particleInfo );
	         
	         Memory_Free( particleInfo );
	    }
	  
	    else if( swarmVar->variable->dataTypes[0] == Variable_DataType_Char) {
	         char* particleInfo = Memory_Alloc_Array( char, swarmVar->dofCount, "particleCheckpointInfo" );
            
            /* Read particle data. */
	         H5Dread( self->fileData[swarmVar_I], H5T_NATIVE_CHAR, memSpace, self->fileSpace[swarmVar_I], H5P_DEFAULT, particleInfo );
	         
	         Variable_SetValue( swarmVar->variable, swarm->particleLocalCount, particleInfo );
	        
	         Memory_Free( particleInfo );
	   }
	           
      else if( swarmVar->variable->dataTypes[0] == Variable_DataType_Float ) {
            float* particleInfo = Memory_Alloc_Array( float, swarmVar->dofCount, "particleCheckpointInfo" );
               
            /* Read particle data. */
	         H5Dread( self->fileData[swarmVar_I], H5T_NATIVE_FLOAT, memSpace, self->fileSpace[swarmVar_I], H5P_DEFAULT, particleInfo );
	         
	         Variable_SetValue( swarmVar->variable, swarm->particleLocalCount, particleInfo );
	         
	         Memory_Free( particleInfo );
	   }
	      
      else {
            double* particleInfo = Memory_Alloc_Array( double, swarmVar->dofCount, "particleCheckpointInfo" );
            
            /* Read particle data. */
	         H5Dread( self->fileData[swarmVar_I], H5T_NATIVE_DOUBLE, memSpace, self->fileSpace[swarmVar_I], H5P_DEFAULT, particleInfo );
	         
	         Variable_SetValue( swarmVar->variable, swarm->particleLocalCount, particleInfo );
	         
	         Memory_Free( particleInfo );
	   }   
	   
	   H5Sclose( memSpace );
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
		



