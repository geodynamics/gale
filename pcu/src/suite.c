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
#include "listener.h"
#include "test.h"
#include "source.h"
#include "suite.h"

extern pcu_suite_t* pcu_cursuite;

void pcu_suite_run( pcu_suite_t* suite, pcu_listener_t* lsnr ) {
   /* Must have a listener. */
   assert( lsnr );

   /* Temporarily set the listener, we need this so assert
      macros can work properly. */
   assert( suite );
   suite->lsnr = lsnr;
   suite->npassed = 0;

   /* Run all sub-suites first. */
   if( suite->subsuites ) {
      pcu_suite_t* sub;

      sub = suite->subsuites;
      while( sub ) {
         pcu_suite_run( sub, lsnr );
         sub = sub->next;
      }
   }

   /* Set this as the active suite. */
   pcu_cursuite = suite;

   /* Begin this suite. */
   lsnr->suitebegin( lsnr, suite );

   /* Temporarily set the current test, once again needed
      for assert macros. */
   suite->curtest = suite->tests;

   /* Loop over all the tests. */
   while( suite->curtest ) {
      /* Run the test. */
      if( suite->setup )
	 suite->setup( suite->data );
      pcu_test_run( suite->curtest, lsnr );
      if( suite->teardown )
	 suite->teardown( suite->data );

      /* Move to the next test. */
      suite->curtest = suite->curtest->next;
   }

   /* End the suite. */
   lsnr->suiteend( lsnr, suite );

   /* Clear temporary settings. */
   suite->curtest = NULL;
   suite->lsnr = NULL;
   pcu_cursuite = NULL;
}

void _pcu_suite_setFixtures( pcu_suite_t* suite, 
			     pcu_fixture_t* setup, pcu_fixture_t* teardown )
{
   assert( suite );
   suite->setup = setup;
   suite->teardown = teardown;
}

void _pcu_suite_setData( pcu_suite_t* suite, int size ) {
   if( suite->data )
      free( suite->data );
   if( size )
      suite->data = malloc( size );
   else
      suite->data = NULL;
}

void _pcu_suite_addTest( pcu_suite_t* suite, pcu_testfunc_t* func, const char* name ) {
   pcu_test_t* test;

   /* Extract test name. */
   /* TODO */

   /* Create the new test. */
   test = (pcu_test_t*)malloc( sizeof(pcu_test_t) );
   test->name = name;
   test->suite = suite;
   test->func = func;
   test->next = NULL;
   test->nsrcs = 0;
   test->srcs = NULL;

   /* Add the new test. */
   if( suite->tests ) {
      pcu_test_t* cur = suite->tests;

      while( cur->next )
	 cur = cur->next;
      cur->next = test;
   }
   else
      suite->tests = test;
   suite->ntests++;
}

void _pcu_suite_addSubSuite( pcu_suite_t* suite, const char* name,
                             void (initfunc)( pcu_suite_t* ) )
{
   pcu_suite_t* subsuite;

   assert( initfunc );

   /* Setup the new suite. */
   subsuite = (pcu_suite_t*)malloc( sizeof(pcu_suite_t) );
   subsuite->name = name;
   subsuite->ntests = 0;
   subsuite->tests = NULL;
   subsuite->npassed = 0;
   subsuite->curtest = NULL;
   subsuite->lsnr = NULL;
   subsuite->next = NULL;
   subsuite->subsuites = NULL;
   subsuite->setup = NULL;
   subsuite->teardown = NULL;
   subsuite->data = NULL;
   initfunc( subsuite );

   /* Add to our list. */
   if( suite->subsuites ) {
      pcu_suite_t* cur;

      cur = suite->subsuites;
      while( cur->next )
	 cur = cur->next;
      cur->next = subsuite;
   }
   else
     suite->subsuites = subsuite;
}

void pcu_suite_clear( pcu_suite_t* suite ) {
   pcu_test_t* tst;
   pcu_source_t* src;
   pcu_suite_t* sub;

   while( suite->tests ) {
      while( suite->tests->srcs ) {
         src = suite->tests->srcs->next;
         pcu_source_clear( suite->tests->srcs );
         free( suite->tests->srcs );
         suite->tests->srcs = src;
      }

      tst = suite->tests->next;
      free( suite->tests );
      suite->tests = tst;
   }
   suite->ntests = 0;
   suite->curtest = NULL;

   while( suite->subsuites ) {
      pcu_suite_clear( suite->subsuites );
      sub = suite->subsuites->next;
      free( suite->subsuites );
      suite->subsuites = sub;
   }

   if( suite->data ) {
      free( suite->data );
      suite->data = NULL;
   }
}
