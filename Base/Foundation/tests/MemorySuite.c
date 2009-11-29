/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**   Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**   Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**   Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**   Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**   Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**   Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Lesser General Public
**  License as published by the Free Software Foundation; either
**  version 2.1 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License along with this library; if not, write to the Free Software
**  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**
** Role:
**   Tests accuracy of memory statistics generation.
**
** $Id: testMemory0.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/Foundation/forwardDecl.h" /* For Journal stuff */
#include "MemorySuite.h"

/* A basic tolerance, not for serious arithmetic calculations but just since we are using doubles 
 * in increments of 0.1 as part of read-write tests */
#define DOUBLE_TOLERANCE 0.00001

typedef struct {
} MemorySuiteData;

void MemorySuite_Setup( MemorySuiteData* data ) {
   #ifdef MEMORY_STATS
   stgMemory->enable = True;
   #endif
}

void MemorySuite_Teardown( MemorySuiteData* data ) {
}


void MemorySuite_Test2DArray( MemorySuiteData* data ) {
   double**       array2d = NULL;
   Index          ii = 0;
   Index          jj = 0;
   MemoryPointer* memoryPtr;

   array2d = Memory_Alloc_2DArray_Unnamed( double, 2, 3 );
   pcu_check_true( NULL != array2d );
   #ifdef MEMORY_STATS
   pcu_check_true( Memory_IsAllocated( array2d ) == True );
   memoryPtr = (MemoryPointer*)Memory_Find_Pointer( array2d );
   pcu_check_true( NULL != memoryPtr );
   pcu_check_true( memoryPtr->allocType == MEMORY_2DARRAY );
   pcu_check_streq( memoryPtr->type->value, "double" );
   pcu_check_streq( memoryPtr->file->value, "StGermain/Base/Foundation/tests/MemorySuite.c" );
   pcu_check_streq( memoryPtr->func->value, __func__ );
   pcu_check_true( memoryPtr->itemSize == sizeof(double));
   pcu_check_true( memoryPtr->totalSize == sizeof(double)*2*3 + sizeof(double*)*2 );
   pcu_check_true( memoryPtr->length.twoD[0] == 2 );
   pcu_check_true( memoryPtr->length.twoD[1] == 3 );
   #endif

   /* Simply write values to specified array entries, and check they're able to be read back properly */
   for ( ii = 0; ii < 2; ii++ ) {
      for ( jj = 0; jj < 3; jj++ ) {
         array2d[ii][jj] = ii + (jj / 10.0);
      }
   }
   for ( ii = 0; ii < 2; ii++ ) {
      for ( jj = 0; jj < 3; jj++ ) {
         pcu_check_true( fabs(array2d[ii][jj] - (ii + (jj / 10.0))) <= DOUBLE_TOLERANCE );
      }
   }

   array2d = Memory_Realloc_2DArray( array2d, double, 4, 4 );
   /* Check that the realloc hasn't corrupted existing values, then re-write and re-check */
   for ( ii = 0; ii < 2; ii++ ) {
      for ( jj = 0; jj < 3; jj++ ) {
         pcu_check_true( fabs(array2d[ii][jj] - (ii + (jj / 10.0))) <= DOUBLE_TOLERANCE );
      }
   }
   for ( ii = 0; ii < 4; ii++ ) {
      for ( jj = 0; jj < 4; jj++ ) {
         array2d[ii][jj] = ii + (jj / 20.0);
      }
   }
   for ( ii = 0; ii < 4; ii++ ) {
      for ( jj = 0; jj < 4; jj++ ) {
         pcu_check_true( fabs(array2d[ii][jj] - (ii + (jj / 20.0))) <= DOUBLE_TOLERANCE );
      }
   }

   #ifdef MEMORY_STATS
   pcu_check_true( Memory_CountGet( array2d ) == 0 );
   Memory_CountInc( array2d );
   pcu_check_true( Memory_CountGet( array2d ) == 1 );
   Memory_CountDec( array2d );
   Memory_CountGet( array2d );
   #endif
   Memory_Free( array2d );
   #ifdef MEMORY_STATS
   pcu_check_true( Memory_IsAllocated( array2d ) == False );
   #endif

}


void MemorySuite_Test3DArray( MemorySuiteData* data ) {
   double***   array3d=NULL;
   Index       ii = 0;
   Index       jj = 0;
   Index       kk = 0;

   array3d = Memory_Alloc_3DArray_Unnamed( double, 4, 3, 2 );
   pcu_check_true( NULL != array3d );
   #ifdef MEMORY_STATS
   pcu_check_true( Memory_IsAllocated( array3d ) == True );
   #endif

   for ( ii = 0; ii < 4; ii++ ) {
      for ( jj = 0; jj < 3; jj++ ) {
         for ( kk = 0; kk < 2; kk++ ) {
            array3d[ii][jj][kk] = ii + (jj / 10.0) + (kk / 100.0);
         }
      }
   }
   for ( ii = 0; ii < 4; ii++ ) {
      for ( jj = 0; jj < 3; jj++ ) {
         for ( kk = 0; kk < 2; kk++ ) {
            pcu_check_true( fabs(array3d[ii][jj][kk] - (ii + (jj / 10.0) + (kk / 100.0))) <= DOUBLE_TOLERANCE );
         }
      }
   }

   array3d = Memory_Realloc_3DArray( array3d, double, 5, 5, 5 );
   /* Check that the realloc hasn't corrupted existing values, then re-write and re-check */
   for ( ii = 0; ii < 4; ii++ ) {
      for ( jj = 0; jj < 3; jj++ ) {
         for ( kk = 0; kk < 2; kk++ ) {
            pcu_check_true( fabs(array3d[ii][jj][kk] - (ii + (jj / 10.0) + (kk / 100.0))) <= DOUBLE_TOLERANCE );
         }
      }
   }
   for ( ii = 0; ii < 5; ii++ ) {
      for ( jj = 0; jj < 5; jj++ ) {
         for ( kk = 0; kk < 5; kk++ ) {
            array3d[ii][jj][kk] = ii + (jj / 20.0) + (kk / 200.0);
         }
      }
   }
   for ( ii = 0; ii < 5; ii++ ) {
      for ( jj = 0; jj < 5; jj++ ) {
         for ( kk = 0; kk < 5; kk++ ) {
            pcu_check_true( fabs(array3d[ii][jj][kk] - (ii + (jj / 20.0) + (kk / 200.0))) <= DOUBLE_TOLERANCE );
         }
      }
   }

   /* Note: this test would work better compiled with MEMORY_STATS enabled. Not sure how to
    * do this in scons yet. -- PatrickSunter, 7 Apr 2009 */
   #ifdef MEMORY_STATS
   pcu_check_true( Memory_CountGet( array3d ) == 0 );
   Memory_CountInc( array3d );
   pcu_check_true( Memory_CountGet( array3d ) == 1 );
   Memory_CountDec( array3d );
   Memory_CountGet( array3d );
   #endif
   Memory_Free( array3d );
   #ifdef MEMORY_STATS
   pcu_check_true( Memory_IsAllocated( array3d ) == False );
   #endif
}


void MemorySuite_Test4DArray( MemorySuiteData* data ) {
   double****   array4d=NULL;
   Index       ii = 0;
   Index       jj = 0;
   Index       kk = 0;
   Index       ll = 0;

   array4d = Memory_Alloc_4DArray_Unnamed( double, 4, 3, 2, 3 );
   pcu_check_true( NULL != array4d );
   #ifdef MEMORY_STATS
   pcu_check_true( Memory_IsAllocated( array4d ) == True );
   #endif

   for ( ii = 0; ii < 4; ii++ ) {
      for ( jj = 0; jj < 3; jj++ ) {
         for ( kk = 0; kk < 2; kk++ ) {
            for ( ll = 0; ll < 3; ll++ ) {
               array4d[ii][jj][kk][ll] = ii + (jj / 10.0) + (kk / 100.0) + (ll / 1000.0);
            }
         }
      }
   }
   for ( ii = 0; ii < 4; ii++ ) {
      for ( jj = 0; jj < 3; jj++ ) {
         for ( kk = 0; kk < 2; kk++ ) {
            for ( ll = 0; ll < 3; ll++ ) {
               pcu_check_true( fabs( array4d[ii][jj][kk][ll] - (ii + (jj / 10.0) + (kk / 100.0)
                               + (ll / 1000.0))) <= DOUBLE_TOLERANCE );
            }
         }
      }
   }

   /* There doesn't yet seem to be a function to Re-Alloc a 4D array, so no test for it.
    *  -- Patrick Sunter, 7 Apr 2009 */

   #ifdef MEMORY_STATS
   pcu_check_true( Memory_CountGet( array4d ) == 0 );
   Memory_CountInc( array4d );
   pcu_check_true( Memory_CountGet( array4d ) == 1 );
   Memory_CountDec( array4d );
   Memory_CountGet( array4d );
   #endif
   Memory_Free( array4d );
   #ifdef MEMORY_STATS
   pcu_check_true( Memory_IsAllocated( array4d ) == False );
   #endif
}


void MemorySuite_Test2DArrayAs1D( MemorySuiteData* data ) {
   double*    one2d = NULL;
   Index      ii = 0;
   Index      jj = 0;

   one2d = Memory_Alloc_2DArrayAs1D_Unnamed( double, 3, 2 );
   pcu_check_true( NULL != one2d );

   for ( ii = 0; ii < 3; ii++ ) {
      for ( jj = 0; jj < 2; jj++ ) {
         Memory_Access2D( one2d, ii, jj, 2 ) = ii + (jj / 10.0);
      }
   }
   for ( ii = 0; ii < 3; ii++ ) {
      for ( jj = 0; jj < 2; jj++ ) {
         pcu_check_true( fabs( Memory_Access2D( one2d, ii, jj, 2 ) - (ii + (jj / 10.0))) <= DOUBLE_TOLERANCE );
      }
   }

   one2d = Memory_Realloc_2DArrayAs1D( one2d, double, 3, 2, 4, 4 );
   /* Check that the realloc hasn't corrupted existing values, then re-write and re-check */
   for ( ii = 0; ii < 3; ii++ ) {
      for ( jj = 0; jj < 2; jj++ ) {
         pcu_check_true( fabs( Memory_Access2D( one2d, ii, jj, 4 ) - (ii + (jj / 10.0))) <= DOUBLE_TOLERANCE );
      }
   }
   for ( ii = 0; ii < 4; ii++ ) {
      for ( jj = 0; jj < 4; jj++ ) {
         Memory_Access2D( one2d, ii, jj, 4 ) = ii + (jj / 20.0);
      }
   }
   for ( ii = 0; ii < 4; ii++ ) {
      for ( jj = 0; jj < 4; jj++ ) {
         pcu_check_true( fabs( Memory_Access2D( one2d, ii, jj, 4 ) - (ii + (jj / 20.0))) <= DOUBLE_TOLERANCE );
      }
   }

   Memory_Free( one2d );
}


void MemorySuite_Test3DArrayAs1D( MemorySuiteData* data ) {
   double*    one3d = NULL;
   Index      ii = 0;
   Index      jj = 0;
   Index      kk = 0;

   one3d = Memory_Alloc_3DArrayAs1D_Unnamed( double, 4, 3, 2 );
   pcu_check_true( NULL != one3d );

   for ( ii = 0; ii < 4; ii++ ) {
      for ( jj = 0; jj < 3; jj++ ) {
         for ( kk = 0; kk < 2; kk++ ) {
            Memory_Access3D( one3d, ii, jj, kk, 3, 2 ) = ii + (jj / 10.0) + (kk / 100.0);
         }
      }
   }
   for ( ii = 0; ii < 4; ii++ ) {
      for ( jj = 0; jj < 3; jj++ ) {
         for ( kk = 0; kk < 2; kk++ ) {
            pcu_check_true( fabs( Memory_Access3D( one3d, ii, jj, kk, 3, 2 ) - 
               (ii + (jj / 10.0) + (kk / 100.0))) <= DOUBLE_TOLERANCE );
         }
      }
   }

   one3d = Memory_Realloc_3DArrayAs1D( one3d, double, 4, 3, 2, 5, 5, 5 );
   /* Check that the realloc hasn't corrupted existing values, then re-write and re-check */
   for ( ii = 0; ii < 4; ii++ ) {
      for ( jj = 0; jj < 3; jj++ ) {
         for ( kk = 0; kk < 2; kk++ ) {
            pcu_check_true( fabs( Memory_Access3D( one3d, ii, jj, kk, 5, 5 ) - 
               (ii + (jj / 10.0) + (kk / 100.0))) <= DOUBLE_TOLERANCE );
         }
      }
   }
   for ( ii = 0; ii < 5; ii++ ) {
      for ( jj = 0; jj < 5; jj++ ) {
         for ( kk = 0; kk < 5; kk++ ) {
            Memory_Access3D( one3d, ii, jj, kk, 5, 5 ) = ii + (jj / 20.0) + (kk / 200.0);
         }
      }
   }
   for ( ii = 0; ii < 5; ii++ ) {
      for ( jj = 0; jj < 5; jj++ ) {
         for ( kk = 0; kk < 5; kk++ ) {
            pcu_check_true( fabs( Memory_Access3D( one3d, ii, jj, kk, 5, 5 ) -
               (ii + (jj / 20.0) + (kk / 200.0))) <= DOUBLE_TOLERANCE );
         }
      }
   }

   Memory_Free( one3d );
}


void MemorySuite_Test4DArrayAs1D( MemorySuiteData* data ) {
   double*    one4d = NULL;
   Index      ii = 0;
   Index      jj = 0;
   Index      kk = 0;
   Index      ll = 0;

   one4d = Memory_Alloc_4DArrayAs1D_Unnamed( double, 4, 3, 2, 3 );
   pcu_check_true( NULL != one4d );

   for ( ii = 0; ii < 4; ii++ ) {
      for ( jj = 0; jj < 3; jj++ ) {
         for ( kk = 0; kk < 2; kk++ ) {
            for ( ll = 0; ll < 3; ll++ ) {
               Memory_Access4D( one4d, ii, jj, kk, ll, 3, 2, 3 ) = 
                  ii + (jj / 10.0) + (kk / 100.0) + (ll / 1000.0);
            }
         }
      }
   }

   for ( ii = 0; ii < 4; ii++ ) {
      for ( jj = 0; jj < 3; jj++ ) {
         for ( kk = 0; kk < 2; kk++ ) {
            for ( ll = 0; ll < 3; ll++ ) {
               pcu_check_true( fabs(Memory_Access4D( one4d, ii, jj, kk, ll, 3, 2, 3 ) - 
                  (ii + (jj / 10.0) + (kk / 100.0) + (ll / 1000.0))) <= DOUBLE_TOLERANCE );
            }
         }
      }
   }

   Memory_Free( one4d );
}


void MemorySuite_Test2DComplexArray( MemorySuiteData* data ) {
   double**    complex2d;
   Index       x1 = 5;
   Index       y1[] = { 1, 2, 3, 4, 5 };
   Index      ii = 0;
   Index      jj = 0;

   complex2d = Memory_Alloc_2DComplex_Unnamed( double, x1, y1 );
   pcu_check_true( NULL != complex2d );

   for (ii = 0; ii < x1; ii++) {
      for (jj = 0; jj < y1[ii]; jj++) {
         complex2d[ii][jj] = ii + (double)(jj / 10.0);
      }
   }

   for (ii = 0; ii < x1; ii++) {
      for (jj = 0; jj < y1[ii]; jj++) {
         pcu_check_true( (complex2d[ii][jj] - (ii + (double)(jj / 10.0))) <= DOUBLE_TOLERANCE );
      }
   }

   Memory_Free( complex2d );
}

   
void MemorySuite_Test3DComplexArray( MemorySuiteData* data ) {
   Index**     setup;
   double***   complex3d;
   Index       x2 = 3;
   Index       y2[] = { 4, 2, 3 };
   Index       ii = 0;
   Index       jj = 0;
   Index       kk = 0;

   setup = Memory_Alloc_3DSetup( x2, y2 );
   pcu_check_true( NULL != setup );
   setup[0][0] = 2;
   setup[0][1] = 3;
   setup[0][2] = 2;
   setup[0][3] = 3;
   setup[1][0] = 1;
   setup[1][1] = 5;
   setup[2][0] = 2;
   setup[2][1] = 4;
   setup[2][2] = 2;

   complex3d = Memory_Alloc_3DComplex_Unnamed( double, x2, y2, setup );   
   pcu_check_true( NULL != complex3d );
   for (ii = 0; ii < x2; ii++) {
      for (jj = 0; jj < y2[ii]; jj++) {
         for (kk = 0; kk < setup[ii][jj]; kk++) {
            complex3d[ii][jj][kk] = ii + (jj / 10.0) + (kk / 100.0);
         }
      }
   }
   for (ii = 0; ii < x2; ii++) {
      for (jj = 0; jj < y2[ii]; jj++) {
         for (kk = 0; kk < setup[ii][jj]; kk++) {
            pcu_check_true( (complex3d[ii][jj][kk] - (ii + (jj / 10.0) + (kk / 100.0))) <= DOUBLE_TOLERANCE);
         }
      }
   }

   Memory_Free( setup );
   Memory_Free( complex3d );
}


void MemorySuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, MemorySuiteData );
   pcu_suite_setFixtures( suite, MemorySuite_Setup, MemorySuite_Teardown );
   pcu_suite_addTest( suite, MemorySuite_Test2DArray );
   pcu_suite_addTest( suite, MemorySuite_Test3DArray );
   pcu_suite_addTest( suite, MemorySuite_Test4DArray );
   pcu_suite_addTest( suite, MemorySuite_Test2DArrayAs1D );
   pcu_suite_addTest( suite, MemorySuite_Test3DArrayAs1D );
   pcu_suite_addTest( suite, MemorySuite_Test4DArrayAs1D );
   pcu_suite_addTest( suite, MemorySuite_Test2DComplexArray );
   pcu_suite_addTest( suite, MemorySuite_Test3DComplexArray );
}


