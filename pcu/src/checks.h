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
**  pcu is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with pcu.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef pcu_checks_h
#define pcu_checks_h

#include "types.h"
#include "listener.h"
#include "suite.h"
#include "source.h"
#include "assert.h"

extern pcu_suite_t* pcu_cursuite;

#define _pcu_check_eval( expr, strexpr, msg, type )     \
   pcu_cursuite->lsnr->checkdone(                       \
      pcu_cursuite->lsnr,                               \
      pcu_test_addSource(                               \
         pcu_cursuite->curtest,                         \
         pcu_source_create( (expr) ? 1 : 0,             \
                            type,                       \
                            __FILE__,                   \
                            __LINE__,                   \
                            #strexpr,                   \
                            msg,                        \
                            pcu_cursuite->curtest )     \
         )                                              \
      )

#define pcu_check_true( expr )			\
   _pcu_check_eval( expr, expr, NULL, "true" )

#define pcu_check_gt( a, b )                                            \
   _pcu_check_eval( (a) > (b), (a) > (b), NULL, "greater than" )

#define pcu_check_lt( a, b )                                    \
   _pcu_check_eval( (a) < (b), (a) < (b), NULL, "less than" )

#define pcu_check_ge( a, b )                                            \
   _pcu_check_eval( (a) >= (b), (a) >= (b), NULL, "greater than or equal" )

#define pcu_check_le( a, b )                                            \
   _pcu_check_eval( (a) <= (b), (a) <= (b), NULL, "less than or equal" )

#ifndef NDEBUG

#define pcu_check_noassert( stmnt )                                     \
   do {                                                                 \
      pcu_jump_ready = 1;                                               \
      if( setjmp( pcu_jump_env ) ) {                                    \
         _pcu_check_eval( 0, stmnt, pcu_assert_cur, "assertion free" ); \
         pcu_assert_cur = NULL;                                         \
      }                                                                 \
      else {                                                            \
         stmnt;                                                         \
      }                                                                 \
      pcu_jump_ready = 0;                                               \
   } while( 0 )

#define pcu_check_assert( stmnt )                                       \
   do {                                                                 \
      pcu_jump_ready = 1;                                               \
      if( setjmp( pcu_jump_env ) )                                      \
         pcu_jump_ready = 0;                                            \
      else {                                                            \
         stmnt;                                                         \
      }                                                                 \
      if( pcu_jump_ready ) {                                            \
         _pcu_check_eval( 0, stmnt, NULL, "assertion generated" );      \
         pcu_jump_ready = 0;                                            \
      }                                                                 \
   } while( 0 )

#else

#define pcu_check_noassert( stmnt ) stmnt
#define pcu_check_assert( stmnt ) assert( 0 )

#endif

#endif
