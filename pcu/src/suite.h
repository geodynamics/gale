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

#ifndef pcu_suite_h
#define pcu_suite_h

struct pcu_suite_t {
      const char* name;
      int ntests;
      pcu_test_t* tests;
      int npassed;
      pcu_test_t* curtest;
      pcu_listener_t* lsnr;
      pcu_suite_t* next;
      pcu_suite_t* subsuites;

      pcu_fixture_t* setup;
      pcu_fixture_t* teardown;
      void* data;
};

void pcu_suite_run( pcu_suite_t* suite, pcu_listener_t* lsnr );
void _pcu_suite_setFixtures( pcu_suite_t* suite, 
			     pcu_fixture_t* setup, pcu_fixture_t* teardown );
void _pcu_suite_addTest( pcu_suite_t* suite, pcu_testfunc_t* func, const char* name );
void _pcu_suite_addSubSuite( pcu_suite_t* suite, const char* name,
                             void (initfunc)( pcu_suite_t* ) );
void _pcu_suite_setData( pcu_suite_t* suite, int size );
void pcu_suite_clear( pcu_suite_t* suite );

#define pcu_suite_setFixtures( suite, setup, teardown )		\
   _pcu_suite_setFixtures( suite, (pcu_fixture_t*)setup,	\
			   (pcu_fixture_t*)teardown )

#define pcu_suite_addTest( suite, func )		\
   _pcu_suite_addTest( suite, (pcu_testfunc_t*)func,	\
		       #func )

#define pcu_suite_addSubSuite( suite, subsuite )                \
   _pcu_suite_addSubSuite( suite, #subsuite, subsuite##_init )

#define pcu_suite_setData( suite, type )	\
   _pcu_suite_setData( suite, sizeof(type) )

#endif
