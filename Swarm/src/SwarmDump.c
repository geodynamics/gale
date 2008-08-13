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
** $Id: SwarmDump.c 3851 2006-10-12 08:57:22Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>

#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>
#include <StgDomain/Utils/Utils.h>

#include "types.h"
#include "SwarmClass.h"
#include "StandardParticle.h"
#include "SwarmDump.h"
#include "SwarmVariable.h"
#include "SwarmVariable_Register.h"

#include <assert.h>
#include <string.h>

#ifdef WRITE_HDF5
#include <hdf5.h>
#endif

const Type SwarmDump_Type = "SwarmDump";


SwarmDump* SwarmDump_New(		
		Name                                               name,
		void*                                              context,
		Swarm**                                            swarmList,
		Index                                              swarmCount,
		Bool                                               newFileEachTime )
{
	SwarmDump* self = _SwarmDump_DefaultNew( name );

	_SwarmDump_Init( self, context, swarmList, swarmCount, newFileEachTime );
	return self;
}

SwarmDump* _SwarmDump_New(
		SizeT                                              _sizeOfSelf, 
		Type                                               type,
		Stg_Class_DeleteFunction*	                       _delete,
		Stg_Class_PrintFunction*	                       _print, 
		Stg_Class_CopyFunction*	                           _copy, 
		Stg_Component_DefaultConstructorFunction*          _defaultConstructor,
		Stg_Component_ConstructFunction*                   _construct,
		Stg_Component_BuildFunction*                       _build,
		Stg_Component_InitialiseFunction*                  _initialise,
		Stg_Component_ExecuteFunction*                     _execute,
		Stg_Component_DestroyFunction*                     _destroy,
		Name                                               name ) 
{
	SwarmDump*		self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(SwarmDump) );
	self = (SwarmDump*)_Stg_Component_New( 
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
			name, 
			NON_GLOBAL );
	
	/* Virtual functions */

	return self;
}

void _SwarmDump_Init( 
		SwarmDump*                                         self,
		void*                                              context,
		Swarm**                                            swarmList,
		Index                                              swarmCount,
		Bool                                               newFileEachTime )
{
	self->isConstructed = True;

	self->swarmList = Memory_Alloc_Array( Swarm*, swarmCount, "swarmList" );
	memcpy( self->swarmList, swarmList, swarmCount * sizeof(Swarm*) );
	self->swarmCount = swarmCount;

	self->newFileEachTime = newFileEachTime;
		
	/* Only append hook to context's save EP if context is given */
	if ( context ) {
		EP_AppendClassHook( Context_GetEntryPoint( context, AbstractContext_EP_SaveClass ), SwarmDump_Execute, self );
	}
}


void _SwarmDump_Delete( void* swarmDump ) {
	SwarmDump* self = (SwarmDump*) swarmDump;
	
	Memory_Free( self->swarmList );
	_Stg_Component_Delete( self );
}

void _SwarmDump_Print( void* _swarmDump, Stream* stream ) {
	SwarmDump* self = (SwarmDump*) _swarmDump;
	Index      swarm_I;

	Journal_Printf( stream, "SwarmDump - '%s'\n", self->name );
	Stream_Indent( stream );
	_Stg_Component_Print( self, stream );

	for ( swarm_I = 0 ; swarm_I < self->swarmCount ; swarm_I++ ) {
		Journal_Printf( stream, "Swarm - '%s'\n", self->swarmList[ swarm_I ]->name );
	}

	Stream_UnIndent( stream );
}

void* _SwarmDump_Copy( void* swarmDump, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	SwarmDump*	self = (SwarmDump*)swarmDump;
	SwarmDump*	newSwarmDump;
	PtrMap*			map = ptrMap;
	Bool			ownMap = False;
	
	if( !map ) {
		map = PtrMap_New( 10 );
		ownMap = True;
	}
	
	newSwarmDump = _Stg_Component_Copy( self, dest, deep, nameExt, map );
	memcpy( newSwarmDump->swarmList, self->swarmList, self->swarmCount * sizeof(Swarm*) );
	newSwarmDump->swarmCount = self->swarmCount;

	if( ownMap ) {
		Stg_Class_Delete( map );
	}
				
	return (void*)newSwarmDump;
}


void* _SwarmDump_DefaultNew( Name name ) {
		return (void*) _SwarmDump_New( 
			sizeof(SwarmDump), 
			SwarmDump_Type, 
			_SwarmDump_Delete, 
			_SwarmDump_Print,
			_SwarmDump_Copy, 
			_SwarmDump_DefaultNew,
			_SwarmDump_Construct,
			_SwarmDump_Build, 
			_SwarmDump_Initialise, 
			_SwarmDump_Execute, 
			_SwarmDump_Destroy, 
			name );
}
void _SwarmDump_Construct( void* swarmDump, Stg_ComponentFactory* cf, void* data ) {
	SwarmDump*	            self         = (SwarmDump*)swarmDump;
	Swarm**                 swarmList;
	AbstractContext*        context;
	Bool                    newFileEachTime;
	Index                   swarmCount;

	context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data ) ;
	swarmList = Stg_ComponentFactory_ConstructByList( 
		cf, 
		self->name, 
		"Swarm", 
		Stg_ComponentFactory_Unlimited, 
		Swarm, 
		True, 
		&swarmCount,
		data ) ;
	newFileEachTime = Stg_ComponentFactory_GetBool( cf, self->name, "newFileEachTime", True );

	_SwarmDump_Init( 
			self,
			context,
			swarmList,
			swarmCount,
			newFileEachTime );

	Memory_Free( swarmList );
}

void _SwarmDump_Build( void* swarmDump, void* data ) {
	SwarmDump*	 self                = (SwarmDump*)     swarmDump;
	Index        swarm_I;

	for ( swarm_I = 0 ; swarm_I < self->swarmCount ; swarm_I++ ) {
		Stg_Component_Build( self->swarmList[ swarm_I ], data, False );
	}
}

void _SwarmDump_Initialise( void* swarmDump, void* data ) {
	SwarmDump*	 self                = (SwarmDump*)     swarmDump;
	Index        swarm_I;

	for ( swarm_I = 0 ; swarm_I < self->swarmCount ; swarm_I++ ) {
		Stg_Component_Initialise( self->swarmList[ swarm_I ], data, False );
	}
}

void _SwarmDump_Execute( void* swarmDump, void* data ) {
	SwarmDump*	      self                = (SwarmDump*)     swarmDump;
	AbstractContext*  context             = Stg_CheckType( data, AbstractContext );
	Stream*           stream              = Journal_Register( MPIStream_Type, Swarm_Type );
	Particle_Index    particleLocalCount;
	SizeT             particleSize;
	Name              filename;
	Index             swarm_I;
	Swarm*            swarm;
	Stream*           info = Journal_Register( Info_Type, self->type );
	Processor_Index   rank_I;

	Journal_DPrintf( info, "Proc %d: beginning Swarm binary checkpoint in %s():\n", self->swarmList[0]->myRank, __func__ );
	Stream_Indent( info );
	
	for ( swarm_I = 0 ; swarm_I < self->swarmCount ; swarm_I++ ) {
		swarm = self->swarmList[ swarm_I ];
		particleLocalCount = swarm->particleLocalCount;
		particleSize = (SizeT) swarm->particleExtensionMgr->finalSize;

		if ( self->newFileEachTime ) {
			if ( strlen(context->checkPointPrefixString) > 0 ) {
				Stg_asprintf( &filename, "%s/%s.%s.%05d", context->checkpointWritePath,
					context->checkPointPrefixString, swarm->name, context->timeStep );
			}
			else {
				Stg_asprintf( &filename, "%s/%s.%05d", context->checkpointWritePath,
					swarm->name, context->timeStep );
			}
		}	
		else { 
			if ( strlen(context->checkPointPrefixString) > 0 ) {
				Stg_asprintf( &filename, "%s/%s.%s", context->checkpointWritePath,
					context->checkPointPrefixString, swarm->name );
			}
			else {
				Stg_asprintf( &filename, "%s/%s", context->checkpointWritePath, swarm->name );
			}
		}	

		for ( rank_I = 0; rank_I < swarm->nProc; rank_I++ ) {
			if ( swarm->myRank == rank_I ) {
				Journal_DPrintf( info, "Proc %d: for swarm \"%s\", dumping its %u particles of size %u bytes "
					"each (= %g bytes total) to file %s\n", swarm->myRank, swarm->name, particleLocalCount,
					particleSize, (float)(particleLocalCount * particleSize), filename );
			}	
			MPI_Barrier( swarm->comm );
		}

#ifdef WRITE_HDF5
      Stg_asprintf( &filename, "%s.h5", filename );
		SwarmDump_DumpToHDF5( self, swarm, filename );
#else
      Stg_asprintf( &filename, "%s.dat", filename );
      Stream_RedirectFile( stream, filename );
		MPIStream_WriteAllProcessors( stream, swarm->particles, particleSize, (SizeT) particleLocalCount, swarm->comm );
		Stream_CloseFile( stream );
#endif

		Memory_Free( filename );
	}
	Stream_UnIndent( info );
	Journal_DPrintf( info, "Proc %d: finished Swarm binary checkpoint.\n", self->swarmList[0]->myRank );
}

void _SwarmDump_Destroy( void* swarmDump, void* data ) {
}

/** Virtual Function Wrappers */
void SwarmDump_Execute( void* swarmDump, void* context ) {
	SwarmDump*	  self                = (SwarmDump*)     swarmDump;

	self->_execute( self, context );
}

#ifdef WRITE_HDF5
void SwarmDump_DumpToHDF5( SwarmDump* self, Swarm* swarm, const char* filename ) {
   hid_t                   file, fileSpace, fileData;
   hid_t                   memSpace;
   hid_t                   props;
   hsize_t                 size[2];
   int                     intSize;
   int                     rank, nRanks, offset;
   hsize_t                 start[2], count[2], slab[2];
   Particle_Index          lParticle_I = 0;
   Particle_Index          particleGlobalCount;
   Stream*                 errorStr = Journal_Register( Error_Type, self->type );
   SwarmVariable*          swarmVar;
   Index                   swarmVar_I, dof_I;
   double*                 value;
   unsigned                maxDofs = 0;
   char                    dataSpaceName[1024];
    
   /* Caluclate the maximum dofCount of all the swarmVariables */   
   for( swarmVar_I = 0; swarmVar_I < swarm->swarmVariable_Register->objects->count; swarmVar_I++ ) {
      swarmVar = SwarmVariable_Register_GetByIndex( swarm->swarmVariable_Register, swarmVar_I );
      if( swarmVar->dofCount > maxDofs )
         maxDofs = swarmVar->dofCount;
   }
   
   /* Allocate space for the values to be written to file */
   value = Memory_Alloc_Array( double, maxDofs, "swarmVariableValue" );
         
   /* Create parallel file property list. */
   props = H5Pcreate( H5P_FILE_ACCESS );
   H5Pset_fapl_mpio( props, MPI_COMM_WORLD, MPI_INFO_NULL );

   /* Open the HDF5 output file. */
   file = H5Fcreate( filename, H5F_ACC_TRUNC, H5P_DEFAULT, props );
   Journal_Firewall( 
		file >= 0, 
		errorStr,
		"Error in %s for %s '%s' - Cannot create file %s.\n", 
		__func__, 
		self->type, 
		self->name, 
		filename );
		
   H5Pclose( props );

   /* Create our output space and data objects. */
   MPI_Allreduce( &swarm->particleLocalCount, &particleGlobalCount, 1, MPI_UNSIGNED, MPI_SUM, swarm->comm );
   size[0] = particleGlobalCount; 

   /* Calculate our file offset. */
   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   MPI_Comm_size( MPI_COMM_WORLD, &nRanks );
   offset = 0;
   if( rank > 0 ) {
      MPI_Status status;
      MPI_Recv( &offset, 1, MPI_INT, rank - 1, 455, MPI_COMM_WORLD, &status );
   }
   start[0] = offset;                     start[1] = 0;
   count[0] = 1;                          count[1] = maxDofs;
   offset += swarm->particleLocalCount;
   if( rank < nRanks - 1 )
      MPI_Send( &offset, 1, MPI_INT, rank + 1, 455, MPI_COMM_WORLD );

   props = H5Pcreate( H5P_DATASET_XFER );
   H5Pset_dxpl_mpio( props, H5FD_MPIO_INDEPENDENT );
   
   /* Loop through the swarmVariable_Register */
   for( swarmVar_I = 0; swarmVar_I < swarm->swarmVariable_Register->objects->count; swarmVar_I++ ) {
      swarmVar = SwarmVariable_Register_GetByIndex( swarm->swarmVariable_Register, swarmVar_I );
       
      size[1] = swarmVar->dofCount;   
      fileSpace = H5Screate_simple( 2, size, NULL );
      count[1] = size[1];
      slab[0] = start[0];
      slab[1] = start[1];
      
      /* Create our memory space. */
      memSpace = H5Screate_simple( 2, count, NULL );
      
      /* Create a new dataspace */
      sprintf( dataSpaceName, "/%s", swarmVar->name );
                  
      #if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
      fileData = H5Dcreate( file, dataSpaceName, H5T_NATIVE_DOUBLE, fileSpace, H5P_DEFAULT );
      #else
      fileData = H5Dcreate( file, dataSpaceName, H5T_NATIVE_DOUBLE, fileSpace,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      #endif
            
      /* Loop through local particles */      
      for( lParticle_I=0; lParticle_I < swarm->particleLocalCount; lParticle_I++ ) {       
         H5Sselect_hyperslab( fileSpace, H5S_SELECT_SET, slab, NULL, count, NULL );
         H5Sselect_all( memSpace );       
   
         /* Write the value of the current swarmVariable at the current particle to the dataspace */
         SwarmVariable_ValueAt( swarmVar, lParticle_I, value );
         H5Dwrite( fileData, H5T_NATIVE_DOUBLE, memSpace, fileSpace, props, value );
         
         slab[0] += 1;
      }
      
      /* Close the dataspace */           
      H5Dclose( fileData );
      H5Sclose( fileSpace );
   }  

   /* Close off all our handles. */
   H5Pclose( props );
   H5Fclose( file );
   
   Memory_Free( value );
}
#endif
