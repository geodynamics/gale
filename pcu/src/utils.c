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
   return str ? pcu_memdup( str, strlen( str ) + 1 ) : NULL;
}
