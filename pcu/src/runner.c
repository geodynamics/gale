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
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "types.h"
#include "runner.h"
#include "suite.h"
#include "test.h"
#include "source.h"

int pcu_nsuites;
pcu_suite_t* pcu_suites;
pcu_suite_t* pcu_cursuite;

int pcu_nnames;
char** pcu_names;

int PCU_PRINT_DOCS=0;

void pcu_runner_searchHierarchy( pcu_suite_t* suite );

void pcu_runner_init( int argc, char* argv[] ) {
   /* Clear global values. */
   pcu_nsuites = 0;
   pcu_suites = NULL;
   pcu_cursuite = NULL;

   /* Extract suite names from the command line. */
   if ( argc > 1 && 0 == strcmp( argv[1], "--withDocs" ) ) {
      pcu_nnames = argc - 2;
      pcu_names = argv + 2;
      PCU_PRINT_DOCS = 1;
   }
   else {
      pcu_nnames = argc - 1;
      pcu_names = argv + 1;
   }

}

void pcu_runner_finalise() {
   pcu_suite_t* tmp;

   while( pcu_suites ) {
      pcu_suite_clear( pcu_suites );
      tmp = pcu_suites->next;
      free( pcu_suites );
      pcu_suites = tmp;
   }
}

PCU_Runner_Status pcu_runner_run( pcu_listener_t* lsnr ) {
   pcu_suite_t* cur;
   unsigned int totalPasses=0; 
   unsigned int totalTests=0; 
   PCU_Runner_Status returnStatus;

   cur = pcu_suites;
   while( cur ) {
      pcu_suite_run( cur, lsnr );
      totalPasses += cur->npassed;
      totalTests += cur->ntests;
      cur = cur->next;
   }

   if ( pcu_nsuites >= 1 ) {
      printf( "-----------------------------------------------------------\n" );
      printf( "[PCU] Total Passes: (%d/%d)\n", totalPasses, totalTests );
      printf( "-----------------------------------------------------------\n" );
   }

   if ( totalPasses == totalTests ) {
      returnStatus = PCU_RUNNER_ALLPASS;
   }
   else {
      returnStatus = PCU_RUNNER_FAILS;
   }
   
   return returnStatus;
}

void _pcu_runner_addSuite( const char* name, void (initfunc)( pcu_suite_t* ), const char* moduleDir ) {
   pcu_suite_t* suite;

   assert( initfunc );
   assert( name );
   assert( moduleDir );

   /* Setup the new suite. */
   suite = (pcu_suite_t*)malloc( sizeof(pcu_suite_t) );
   suite->name = strdup( name );
   suite->moduleDir = strdup( moduleDir );
   suite->ntests = 0;
   suite->tests = NULL;
   suite->npassed = 0;
   suite->curtest = NULL;
   suite->lsnr = NULL;
   suite->next = NULL;
   suite->nsubsuites = 0;
   suite->subsuites = NULL;
   suite->setup = NULL;
   suite->teardown = NULL;
   suite->data = NULL;
   initfunc( suite );

   /* Don't add the suite if it's not in our list of names. */
   pcu_runner_searchHierarchy( suite );
}

void pcu_runner_searchHierarchy( pcu_suite_t* suite ) {
   if( pcu_nnames ) {
      int ii;

      for( ii = 0; ii < pcu_nnames; ii++ ) {
	 if( !strcmp( pcu_names[ii], suite->name ) )
	    break;
      }
      if( ii == pcu_nnames ) {
         pcu_suite_t* cur;

         cur = suite->subsuites;
         while( cur ) {
            pcu_runner_searchHierarchy( cur );
            cur = cur->next;
         }
         pcu_suite_clear( suite );
         free( suite );
         return;
      }
   }

   /* Add to the list of current suites. */
   if( pcu_suites ) {
      pcu_suite_t* cur;

      cur = pcu_suites;
      while( cur->next )
	 cur = cur->next;
      cur->next = suite;
   }
   else
      pcu_suites = suite;
   pcu_nsuites++;
}


