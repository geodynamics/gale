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

#ifndef pcu_assert_h
#define pcu_assert_h

#include <setjmp.h>
#include <assert.h>
#include "checks.h"

#ifndef NDEBUG

extern int pcu_jump_ready;
extern jmp_buf pcu_jump_env;
extern const char* pcu_assert_cur;

#define pcu_assert( expr )                                      \
   (pcu_jump_ready ?                                            \
    (!(expr) ?                                                  \
     (pcu_assert_cur = #expr, longjmp( pcu_jump_env, 1 )) :     \
      (void)0) :                                                \
     assert( expr ))

#else

#define pcu_assert( expr ) assert( expr )

#endif

#endif
