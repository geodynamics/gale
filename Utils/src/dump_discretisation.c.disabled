#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_HDF5
#include <hdf5.h>
#endif
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>


#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
#define H5DCREATE( file, path, type, space, props1, props2, props3 )	\
   H5Dcreate( file, path, type, space, props1 )
#else
#define H5DCREATE( file, path, type, space, props1, props2, props3 )	\
   H5Dcreate( file, path, type, space, props1, props2, props3 )
#endif


void dump_discretisation( Mesh *mesh, Swarm *swarm, const char *filename ) {
#ifdef HAVE_HDF5
   hid_t file, fileSpace, fileData;
   hid_t memSpace;
   hid_t props;
   hid_t particleType;
   hsize_t size[2];
   int intSize[2];
   int rank, nRanks, offset;
   hsize_t start[2], count[2];
   int nDims, nLocals, nGlobals;
   char dataName[13];
   int nTrashBytes;
   int nParticles;
   int ii;

   /* Create parallel file property list. */
   props = H5Pcreate( H5P_FILE_ACCESS );
   H5Pset_fapl_mpio( props, MPI_COMM_WORLD, MPI_INFO_NULL );

   /* Open the HDF5 output file. */
   file = H5Fcreate( filename, H5F_ACC_TRUNC, H5P_DEFAULT, props );
   assert( file );
   H5Pclose( props );

   /* Dump the size so we don't have to do any divisions later on. */
   props = H5Pcreate( H5P_DATASET_XFER );
   H5Pset_dxpl_mpio( props, H5FD_MPIO_COLLECTIVE );
   size[0] = 1;
   fileSpace = H5Screate_simple( 1, size, NULL );
   fileData = H5Dcreate( file, "/numDims", H5T_NATIVE_INT, fileSpace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
   nDims = Mesh_GetDimSize( mesh );
   H5Dwrite( fileData, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, props, &nDims );
   H5Dclose( fileData );
   nLocals = Mesh_GetLocalSize( mesh, 0 );
   MPI_Allreduce( &nLocals, &nGlobals, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
   fileData = H5Dcreate( file, "/numMeshGlobals", H5T_NATIVE_INT, fileSpace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
   H5Dwrite( fileData, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, props, &nGlobals );
   H5Dclose( fileData );
   MPI_Allreduce( &swarm->particleLocalCount, &nParticles, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
   fileData = H5Dcreate( file, "/numParticles", H5T_NATIVE_INT, fileSpace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
   H5Dwrite( fileData, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, props, &nParticles );
   H5Dclose( fileData );
   H5Sclose( fileSpace );
   H5Pclose( props );

   /* Dump mesh vertices */
   size[0] = nGlobals;
   size[1] = nDims;
   fileSpace = H5Screate_simple( 2, size, NULL );
   fileData = H5Dcreate( file, "/meshVertices", H5T_NATIVE_CHAR, fileSpace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   MPI_Comm_size( MPI_COMM_WORLD, &nRanks );
   offset = 0;
   if( rank > 0 ) {
      MPI_Status status;
      MPI_Recv( &offset, 1, MPI_INT, rank - 1, 455, MPI_COMM_WORLD, &status );
   }
   start[0] = offset;  start[1] = 0;
   count[0] = nLocals; count[1] = nDims;
   offset += nLocals;
   if( rank < nRanks - 1 )
      MPI_Send( &offset, 1, MPI_INT, rank + 1, 455, MPI_COMM_WORLD );
   H5Sselect_hyperslab( fileSpace, H5S_SELECT_SET, start, NULL, count, NULL );

   memSpace = H5Screate_simple( 2, count, NULL );

   props = H5Pcreate( H5P_DATASET_XFER );
   H5Pset_dxpl_mpio( props, H5FD_MPIO_INDEPENDENT );
   H5Dwrite( fileData, H5T_NATIVE_CHAR, memSpace, fileSpace, props, mesh->verts );

   H5Pclose( props );
   H5Sclose( memSpace );
   H5Dclose( fileData );
   H5Sclose( fileSpace );

   /* Dump particle coordinates. */
   particleType = H5Tcreate( H5T_COMPOUND, swarm->particleExtensionMgr->finalSize );
   offset = 0;
   H5Tinsert( particleType, "owning cell", 0, H5T_NATIVE_INT );
   offset += sizeof(int);
   for( ii = 0; ii < nDims; ii++ ) {
      sprintf( dataName, "coordinate %d", ii );
      H5Tinsert( particleType, dataName, offset, H5T_NATIVE_DOUBLE );
      offset += sizeof(double);
   }
   nTrashBytes = swarm->particleExtensionMgr->finalSize - sizeof(int) - nDims * sizeof(double);
   for( ii = 0; ii < nTrashBytes; ii++ ) {
      sprintf( dataName, "trash %d", ii );
      H5Tinsert( particleType, dataName, offset, H5T_NATIVE_CHAR );
      offset++;
   }

   size[0] = nParticles;
   fileSpace = H5Screate_simple( 1, size, NULL );
   fileData = H5Dcreate( file, "/particles", particleType, fileSpace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   MPI_Comm_size( MPI_COMM_WORLD, &nRanks );
   offset = 0;
   if( rank > 0 ) {
      MPI_Status status;
      MPI_Recv( &offset, 1, MPI_INT, rank - 1, 455, MPI_COMM_WORLD, &status );
   }
   start[0] = offset;                    start[1] = 0;
   count[0] = swarm->particleLocalCount; count[1] = 1;
   offset += swarm->particleLocalCount;
   if( rank < nRanks - 1 )
      MPI_Send( &offset, 1, MPI_INT, rank + 1, 455, MPI_COMM_WORLD );
   H5Sselect_hyperslab( fileSpace, H5S_SELECT_SET, start, NULL, count, NULL );

   size[0] = swarm->particleLocalCount;
   memSpace = H5Screate_simple( 1, size, NULL );

   props = H5Pcreate( H5P_DATASET_XFER );
   H5Pset_dxpl_mpio( props, H5FD_MPIO_INDEPENDENT );
   H5Dwrite( fileData, particleType, memSpace, fileSpace, props, mesh->verts );

   H5Pclose( props );
   H5Sclose( memSpace );
   H5Dclose( fileData );
   H5Sclose( fileSpace );

   /* Close off all our handles. */
   H5Fclose( file );
#else
   printf( "*** Warning: cannot dump discretisation, not configured with HDF5.\n" );
#endif
}
