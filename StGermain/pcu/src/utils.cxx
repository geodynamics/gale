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
#include "utils.h"

void* pcu_memdup( const void* buf, int len ) {
   void* tmp;

   tmp = malloc( len );
   memcpy( tmp, buf, len );
   return tmp;
}

char* pcu_strdup( const char* str ) {
   return str ? (char*)pcu_memdup( str, strlen( str ) + 1 ) : NULL;
}


