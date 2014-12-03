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

#ifndef pcu_types_h
#define pcu_types_h

typedef struct pcu_source_t     pcu_source_t;
typedef struct pcu_test_t       pcu_test_t;
typedef struct pcu_suite_t      pcu_suite_t;
typedef struct pcu_listener_t   pcu_listener_t;
typedef struct pcu_textoutput_t pcu_textoutput_t;
typedef struct pcu_runner_t     pcu_runner_t;

typedef void (pcu_testfunc_t)   ( void* data );
typedef void (pcu_fixture_t)    ( void* data );
typedef void (pcu_suiteentry_t) ( pcu_listener_t* lsnr, pcu_suite_t* suite );
typedef void (pcu_testentry_t)  ( pcu_listener_t* lsnr, pcu_test_t* test );
typedef void (pcu_checkentry_t) ( pcu_listener_t* lsnr, pcu_source_t* src );

#endif
