#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "types.h"
#include "listener.h"
#include "test.h"
#include "suite.h"
#include "textoutput.h"
#include "test.h"
#include "source.h"

typedef struct {
      int rank;
} textoutputdata_t;

void printstatus( pcu_listener_t* lsnr, pcu_suite_t* suite, int final );
void printsources( pcu_listener_t* lsnr, pcu_suite_t* suite );

void pcu_textoutput_suitebegin( pcu_listener_t* lsnr, pcu_suite_t* suite ) {
   printstatus( lsnr, suite, 0 );
}

void pcu_textoutput_suiteend( pcu_listener_t* lsnr, pcu_suite_t* suite ) {
   printstatus( lsnr, suite, 1 );
   printsources( lsnr, suite );
}

void pcu_textoutput_testbegin( pcu_listener_t* lsnr, pcu_test_t* test ) {
}

void pcu_textoutput_testend( pcu_listener_t* lsnr, pcu_test_t* test ) {
   printstatus( lsnr, test->suite, 0 );
}

void pcu_textoutput_assertdone( pcu_listener_t* lsnr, pcu_source_t* src ) {
}

pcu_listener_t* pcu_textoutput_create() {
   pcu_listener_t* lsnr;

   lsnr = (pcu_listener_t*)malloc( sizeof(pcu_listener_t) );
   lsnr->suitebegin = pcu_textoutput_suitebegin;
   lsnr->suiteend = pcu_textoutput_suiteend;
   lsnr->testbegin = pcu_textoutput_testbegin;
   lsnr->testend = pcu_textoutput_testend;
   lsnr->assertdone = pcu_textoutput_assertdone;
   lsnr->data = malloc( sizeof(textoutputdata_t) );
   MPI_Comm_rank( MPI_COMM_WORLD, 
		  &((textoutputdata_t*)lsnr->data)->rank );

   return lsnr;
}

void pcu_textoutput_destroy( pcu_listener_t* lsnr ) {
   if( lsnr->data )
      free( lsnr->data );
   free( lsnr );
}

void printstatus( pcu_listener_t* lsnr, pcu_suite_t* suite, int final ) {
   if( ((textoutputdata_t*)lsnr->data)->rank )
      return;

   printf( "Running suite '%s', passes: %d/%d%s", 
	   suite->name, suite->npassed, suite->ntests, 
	   final ? "\n" : "\r" );
}

void printsources( pcu_listener_t* lsnr, pcu_suite_t* suite ) {
   pcu_test_t* test;
   pcu_source_t* src;
   int nfails;

   if( ((textoutputdata_t*)lsnr->data)->rank )
      return;

   nfails = 0;
   test = suite->tests;
   while( test ) {
      src = test->srcs;
      while( src ) {
	 if( !src->result ) {
	    printf( "\n\tAssert '%s' failed:\n", src->type );
	    printf( "\t\tLocation: \t%s:%d\n", src->file, src->line );
	    printf( "\t\tTest name: \t%s\n", src->test->name );
	    printf( "\t\tExpression: \t%s\n", src->expr );
	    printf( "\t\tRank: \t\t%d\n", src->rank );
	    nfails++;
	 }
	 src = src->next;
      }
      test = test->next;
   }

   if( nfails )
     printf( "\n" );
}
