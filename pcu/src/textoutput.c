/*
**  Copyright 2008 Luke Hodkinson
**
**  This file is part of pcu.
**
**  pcu is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  Foobar is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <mpi.h>
#include "types.h"
#include "listener.h"
#include "test.h"
#include "suite.h"
#include "textoutput.h"
#include "test.h"
#include "source.h"

extern int PCU_PRINT_DOCS;

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

void pcu_textoutput_checkdone( pcu_listener_t* lsnr, pcu_source_t* src ) {
}

pcu_listener_t* pcu_textoutput_create( int printdocs ) {
   pcu_listener_t* lsnr;

   lsnr = (pcu_listener_t*)malloc( sizeof(pcu_listener_t) );
   lsnr->suitebegin = pcu_textoutput_suitebegin;
   lsnr->suiteend = pcu_textoutput_suiteend;
   lsnr->testbegin = pcu_textoutput_testbegin;
   lsnr->testend = pcu_textoutput_testend;
   lsnr->checkdone = pcu_textoutput_checkdone;
   lsnr->data = malloc( sizeof(textoutputdata_t) );
   assert( printdocs == 1 || printdocs == 0 );
   lsnr->printdocs = printdocs;
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
      if( lsnr->printdocs ) {
         if (test->globalresult ) {
            printf( " * (P) Test %s: ", test->name );
         }
         else {
            printf( " * (F) Test %s: ", test->name );
         }
         if ( test->docString ) {
            printf( "%s\n", test->docString );
         }
         else {
            printf( "(undocumented)\n" );
         }
      }
      src = test->srcs;
      while( src ) {
         if( !src->result ) {
            printf( "\n\tCheck '%s' failed:\n", src->type );
            printf( "\t\tLocation: \t%s:%d\n", src->file, src->line );
            printf( "\t\tTest name: \t%s\n", src->test->name );
            printf( "\t\tExpression: \t%s\n", src->expr );
            if( src->msg )
               printf( "\t\tMessage: \t%s\n", src->msg );
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
