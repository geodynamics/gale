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

#ifndef pcu_test_h
#define pcu_test_h

struct pcu_test_t {
      const char* name;
      pcu_suite_t* suite;
      pcu_testfunc_t* func;
      pcu_test_t* next;

      int globalresult;
      int nsrcs;
      pcu_source_t* srcs;
};

void pcu_test_run( pcu_test_t* test, pcu_listener_t* lsnr );
pcu_source_t* pcu_test_addSource( pcu_test_t* test, pcu_source_t* src );

#endif
