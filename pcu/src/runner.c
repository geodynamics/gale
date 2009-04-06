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

void pcu_runner_searchHierarchy( pcu_suite_t* suite );

void pcu_runner_init( int argc, char* argv[] ) {
   /* Clear global values. */
   pcu_nsuites = 0;
   pcu_suites = NULL;
   pcu_cursuite = NULL;

   /* Extract suite names from the command line. */
   pcu_nnames = argc - 1;
   pcu_names = argv + 1;
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

void pcu_runner_run( pcu_listener_t* lsnr ) {
   pcu_suite_t* cur;

   cur = pcu_suites;
   while( cur ) {
      pcu_suite_run( cur, lsnr );
      cur = cur->next;
   }
}

void _pcu_runner_addSuite( const char* name, 
			   void (initfunc)( pcu_suite_t* ) )
{
   pcu_suite_t* suite;

   assert( initfunc );

   /* Setup the new suite. */
   suite = (pcu_suite_t*)malloc( sizeof(pcu_suite_t) );
   suite->name = name;
   suite->ntests = 0;
   suite->tests = NULL;
   suite->npassed = 0;
   suite->curtest = NULL;
   suite->lsnr = NULL;
   suite->next = NULL;
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
