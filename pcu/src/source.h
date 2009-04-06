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

#ifndef pcu_source_h
#define pcu_source_h

struct pcu_source_t {
      int result;
      char* type;
      char* file;
      int line;
      char* expr;
      char* msg;
      pcu_test_t* test;
      int rank;
      pcu_source_t* next;
};

void pcu_source_init( pcu_source_t* src );
pcu_source_t* pcu_source_create( int result, const char* type,
				 const char* file, int line,
				 const char* expr, const char* msg,
				 pcu_test_t* test );
int pcu_source_getPackLen( pcu_source_t* src );
void pcu_source_pack( pcu_source_t* src, void* buf );
void pcu_source_unpack( pcu_source_t* src, void* buf );
void pcu_source_clear( pcu_source_t* src );

#endif
