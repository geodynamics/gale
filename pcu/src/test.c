#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include "types.h"
#include "test.h"
#include "source.h"
#include "suite.h"
#include "listener.h"

void pcu_test_gathersources( pcu_test_t* test );

void pcu_test_run( pcu_test_t* test, pcu_listener_t* lsnr ) {
   int rank;
   int passed;
   pcu_source_t* src;

   MPI_Comm_rank( MPI_COMM_WORLD, &rank );

   /* Must have a listener. */
   assert( lsnr );

   /* Begin this test. */
   lsnr->testbegin( lsnr, test );

   /* Run the test. */
   assert( test->func );
   assert( test->suite );
   test->func( test->suite->data );

   /* Need to collect information from all ranks to
      determine if the test passed. */
   passed = 1;
   src = test->srcs;
   while( src ) {
      if( !src->result ) {
	 passed = 0;
	 break;
      }
      src = src->next;
   }
   MPI_Reduce( &passed, &test->globalresult, 1, MPI_INT, MPI_LAND, 
	       0, MPI_COMM_WORLD ); /* did anyone fail? */

   /* Update the suite's passed count. */
   if( rank == 0 && test->globalresult )
      test->suite->npassed++;

   /* Gather up all the sources onto the master rank. */
   pcu_test_gathersources( test );

   /* End the test. */
   lsnr->testend( lsnr, test );
}

pcu_source_t* pcu_test_addSource( pcu_test_t* test, pcu_source_t* src ) {
   assert( test );
   assert( src );
   if( test->srcs ) {
      pcu_source_t* cur = test->srcs;

      while( cur->next )
	 cur = cur->next;
      cur->next = src;
   }
   else
      test->srcs = src;
   test->nsrcs++;

   return src;
}

void pcu_test_gathersources( pcu_test_t* test ) {
   int rank, nranks;
   int buflen;
   pcu_source_t* cur;
   void* buf;
   void* ptr;
   int totalsize;
   void* totalbuf;
   int* alllens;
   int* disps;
   int ii;

   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   MPI_Comm_size( MPI_COMM_WORLD, &nranks );

   /* Pack all our sources. */
   buflen = 0;
   cur = test->srcs;
   while( cur ) {
      buflen += pcu_source_getPackLen( cur );
      cur = cur->next;
   }
   if( buflen )
      buf = malloc( buflen );
   else
      buf = NULL;
   ptr = buf;
   while( test->srcs ) {
      pcu_source_pack( test->srcs, ptr );
      ptr += pcu_source_getPackLen( test->srcs );
      cur = test->srcs->next;
      pcu_source_clear( test->srcs );
      free( test->srcs );
      test->srcs = cur;
   }
   test->nsrcs = 0;
   test->srcs = NULL;

   /* Gather them all up. */
   totalsize = 0;
   MPI_Reduce( &buflen, &totalsize, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
   if( rank == 0 && totalsize ) {
      totalbuf = malloc( totalsize );
      alllens = (int*)malloc( nranks * sizeof(int) );
      disps = (int*)malloc( nranks * sizeof(int) );
   }
   MPI_Gather( &buflen, 1, MPI_INT, alllens, 1, MPI_INT, 0, MPI_COMM_WORLD );
   if( rank == 0 ) {
      disps[0] = 0;
      for( ii = 1; ii < nranks; ii++ )
	 disps[ii] = disps[ii - 1] + alllens[ii - 1];
   }
   MPI_Gatherv( buf, buflen, MPI_BYTE, totalbuf, alllens, disps, 
		MPI_BYTE, 0, MPI_COMM_WORLD );
   if( buf )
      free( buf );

   if( rank == 0 ) {
      /* Free arrays. */
      if( alllens )
	 free( alllens );
      if( disps )
	 free( disps );

      /* Unpack sources into the list. */
      ptr = totalbuf;
      while( ptr < totalbuf + totalsize ) {
	 cur = (pcu_source_t*)malloc( sizeof(pcu_source_t) );
	 pcu_source_init( cur );
	 pcu_source_unpack( cur, ptr );
	 cur->test = test;
	 ptr += pcu_source_getPackLen( cur );
	 pcu_test_addSource( test, cur );
      }

      /* Free global buffer. */
      if( totalbuf )
	 free( totalbuf );
   }
}
