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
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include "types.h"
#include "utils.h"
#include "source.h"
#include "test.h"
#include "suite.h"

void pcu_source_init( pcu_source_t* src ) {
   assert( src );
   memset( src, 0, sizeof(pcu_source_t) );
}

pcu_source_t* pcu_source_create( int result, const char* type,
				 const char* file, int line,
				 const char* expr, const char* msg,
				 pcu_test_t* test )
{
   pcu_source_t* src;

   src = (pcu_source_t*)malloc( sizeof(pcu_source_t) );
   src->result = result;
   src->type = pcu_strdup( type );
   src->file = pcu_strdup( file );
   src->line = line;
   src->expr = pcu_strdup( expr );
   src->msg = pcu_strdup( msg );
   src->test = test;
   src->next = NULL;
   MPI_Comm_rank( MPI_COMM_WORLD, &src->rank );

   return src;
}

int pcu_source_getPackLen( pcu_source_t* src ) {
   int len = 0;

   len += sizeof(int);

   len += sizeof(int);
   len += strlen( src->type ) + 1;

   len += sizeof(int);
   len += strlen( src->file ) + 1;

   len += sizeof(int);

   len += sizeof(int);
   len += strlen( src->expr ) + 1;

   len += sizeof(int);
   if( src->msg )
      len += strlen( src->msg ) + 1;

   len += sizeof(int);

   return len;
}

void pcu_source_pack( pcu_source_t* src, void* buf ) {
   char* tmp = (char*)buf;
   int len;

   assert( tmp );

   *((int*)tmp) = src->result;
   tmp += sizeof(int);

   len = strlen( src->type ) + 1;
   *((int*)tmp) = len;
   tmp += sizeof(int);
   memcpy( tmp, src->type, len * sizeof(char) );
   tmp += len;

   len = strlen( src->file ) + 1;
   *((int*)tmp) = len;
   tmp += sizeof(int);
   memcpy( tmp, src->file, len * sizeof(char) );
   tmp += len;

   *((int*)tmp) = src->line;
   tmp += sizeof(int);

   len = strlen( src->expr ) + 1;
   *((int*)tmp) = len;
   tmp += sizeof(int);
   memcpy( tmp, src->expr, len * sizeof(char) );
   tmp += len;

   if( src->msg ) {
      len = strlen( src->msg ) + 1;
      *((int*)tmp) = len;
      tmp += sizeof(int);
      memcpy( tmp, src->msg, len * sizeof(char) );
      tmp += len;
   }
   else {
      *((int*)tmp) = 0;
      tmp += sizeof(int);
   }

   *((int*)tmp) = src->rank;
   tmp += sizeof(int);
}

void pcu_source_unpack( pcu_source_t* src, void* buf ) {
   char* tmp = (char*)buf;
   int len;

   pcu_source_clear( src );

   src->result = *(int*)tmp;
   tmp += sizeof(int);

   len = *(int*)tmp;
   tmp += sizeof(int);
   src->type = (char*)pcu_memdup( tmp, len );
   tmp += len;

   len = *(int*)tmp;
   tmp += sizeof(int);
   src->file = (char*)pcu_memdup( tmp, len );

   tmp += len;
   src->line = *(int*)tmp;
   tmp += sizeof(int);

   len = *(int*)tmp;
   tmp += sizeof(int);
   src->expr = (char*)pcu_memdup( tmp, len );
   tmp += len;

   len = *(int*)tmp;
   tmp += sizeof(int);
   if( len ) {
      src->msg = (char*)pcu_memdup( tmp, len );
      tmp += len;
   }

   src->rank = *(int*)tmp;
   tmp += sizeof(int);

   src->next = NULL;
}

void pcu_source_clear( pcu_source_t* src ) {
   assert( src );
   if( src->type ) {
      free( src->type );
      src->type = NULL;
   }
   if( src->file ) {
      free( src->file );
      src->file = NULL;
   }
   if( src->expr ) {
      free( src->expr );
      src->expr = NULL;
   }
   if( src->msg ) {
      free( src->msg );
      src->msg = NULL;
   }
}


