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
** $Id: testIndexSet.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "StGermain/Base/Container/Container.h"
#include "IndexSetSuite.h"


#define IS1_SIZE   25
#define IS2_SIZE   30

typedef struct {
      IndexSet* is;
      IndexSet* is2;
} IndexSetSuiteData;


void IndexSetSuite_Setup( IndexSetSuiteData* data ) {
   // Index sets are deliberately different sizes, to test all aspects of merge functionality
   data->is = IndexSet_New(IS1_SIZE);
   data->is2 = IndexSet_New(IS2_SIZE);
}

void IndexSetSuite_Teardown( IndexSetSuiteData* data ) {
   Stg_Class_Delete( data->is );
   Stg_Class_Delete( data->is2 );
}


/* Start with this test, so we can be confident to returen an IS to a zero state, for other tests */
void IndexSetSuite_TestRemoveAll( IndexSetSuiteData* data ) {
   Index ii;

  /* a couple of additions, so we know there's something to remove */
   IndexSet_Add( data->is, 0 );
   IndexSet_Add( data->is, 24 );
   IndexSet_RemoveAll( data->is );
   for ( ii=0; ii< IS1_SIZE; ii++ ) {
      pcu_assert_true( 0 == IndexSet_IsMember( data->is, ii ) );
   }
}


void IndexSetSuite_TestInsertion( IndexSetSuiteData* data ) {
   Index ii=0;

   IndexSet_RemoveAll( data->is );
   IndexSet_Add( data->is, 0 );
   IndexSet_Add( data->is, 7 );
   IndexSet_Add( data->is, 8 );
   IndexSet_Add( data->is, 22 );
   IndexSet_Add( data->is, 24 );
   for ( ii=0; ii< IS1_SIZE; ii++ ) {
      if ( (ii==0) || (ii==7) || (ii==8) || (ii==22) || (ii==24) ) {
         pcu_assert_true( True == IndexSet_IsMember( data->is, ii ) );
      }
      else {
         pcu_assert_true( False == IndexSet_IsMember( data->is, ii ) );
      }
   }
}


void IndexSetSuite_TestRemoval( IndexSetSuiteData* data ) {
   Index ii=0;

   IndexSet_RemoveAll( data->is );
   IndexSet_Add( data->is, 0 );
   IndexSet_Add( data->is, 7 );
   IndexSet_Add( data->is, 8 );
   IndexSet_Add( data->is, 22 );
   IndexSet_Add( data->is, 24 );

   //Now remove a couple of these again
   IndexSet_Remove( data->is, 7 );
   IndexSet_Remove( data->is, 24 );

   //Thus list of members should be (0,7,8,22,24) - (7,24) = (0,8,22)

   for ( ii=0; ii< IS1_SIZE; ii++ ) {
      //printf( "index %u: value %u\n", ii, IndexSet_IsMember( data->is, ii ) );
      if ( (ii==0) || (ii==8) || (ii==22) ) {
         pcu_assert_true( True == IndexSet_IsMember( data->is, ii ) );
      }
      else {
         pcu_assert_true( False == IndexSet_IsMember( data->is, ii ) );
      }
   }
}

void IndexSetSuite_TestUpdateMembersCount( IndexSetSuiteData* data ) {
   Index ii=0;

   IndexSet_RemoveAll( data->is );

   IndexSet_UpdateMembersCount( data->is );
   pcu_assert_true( 0 == data->is->membersCount );

   // Add some members, to generate a count
   for ( ii=0; ii < 5; ii++ ) {
      IndexSet_Add( data->is, ii );
   }
   for ( ii=20; ii < 25; ii++ ) {
      IndexSet_Add( data->is, ii );
   }
   IndexSet_UpdateMembersCount( data->is );
   pcu_assert_true( 10 == data->is->membersCount );

}


#if 0
void IndexSetSuite_TestGetIndexOfNthMember( IndexSetSuiteData* data ) {
   for( i = 0; i <= data->is->membersCount; i++ ){
      printf( "* Test GetIndexOfNthMember %d*\n", i );
      index = IndexSet_GetIndexOfNthMember( data->is, i );
      printf( "Index of member %d=%d", i, index );
      if ( IndexSet_Invalid(data->is) == index ) {
         printf(" (invalid)");
      }
      printf("\n");
   }
}
#endif


void IndexSetSuite_TestGetMembers( IndexSetSuiteData* data ) {
   Index*	setArray;
   unsigned int	setArraySize;

   IndexSet_RemoveAll( data->is );

   // add some test set members
   IndexSet_Add( data->is, 0 );
   IndexSet_Add( data->is, 7 );
   IndexSet_Add( data->is, 8 );
   IndexSet_Add( data->is, 22 );
   IndexSet_Add( data->is, 24 );

   IndexSet_GetMembers( data->is, &setArraySize, &setArray );
   pcu_assert_true( 5 == setArraySize );
   pcu_assert_true( 0 == setArray[0] );
   pcu_assert_true( 7 == setArray[1] );
   pcu_assert_true( 8 == setArray[2] );
   pcu_assert_true( 22 == setArray[3] );
   pcu_assert_true( 24 == setArray[4] );

   Memory_Free( setArray );
}


void IndexSetSuite_TestGetVacancies( IndexSetSuiteData* data ) {
   Index        ii;
   Index*	setArray;
   unsigned int	setArraySize;

   IndexSet_RemoveAll( data->is );

   // set all indices to be included, except those exactly divisible by 3
   for ( ii=0; ii< IS1_SIZE; ii++ ) {
      if ( 0 != (ii % 5) ) {
         IndexSet_Add( data->is, ii );
      } 
   }

   IndexSet_GetVacancies( data->is, &setArraySize, &setArray );
   // All should be enabled except for 0, 5, 10, 15, 20 - thus count should be 5 
   pcu_assert_true( 5 == setArraySize );
   pcu_assert_true( 0 == setArray[0] );
   pcu_assert_true( 5 == setArray[1] );
   pcu_assert_true( 10 == setArray[2] );
   pcu_assert_true( 15 == setArray[3] );
   pcu_assert_true( 20 == setArray[4] );

   Memory_Free( setArray );
}


#if 0
void IndexSetSuite_TestAddAll( IndexSetSuiteData* data ) {
   printf( "* Test AddAll *\n" );
   IndexSet_AddAll( data->is );
   Print_Container( data->is );
}


void IndexSetSuite_TestDuplicate( IndexSetSuiteData* data ) {
   printf( "* Test Duplicate *\n" );
   is2 = IndexSet_Duplicate( data->is );
}


void IndexSetSuite_TestMerge_OR( IndexSetSuiteData* data ) {
   printf( "* Test Merge_OR *\n" );
   IndexSet_Add( data->is, 0 );
   IndexSet_Add( data->is, 7 );
   IndexSet_Add( data->is2, 7 );
   IndexSet_Add( data->is2, 8 );
   IndexSet_Add( data->is2, 22 );
   printf( "Pre_Merge\n" );
   Print_Container( data->is );
   Print_Container( data->is2 );
   IndexSet_Merge_OR( data->is, data->is2 );
   printf( "Post_Merge\n" );
   Print_Container( data->is );
   Print_Container( data->is2 );
}
#endif

void IndexSetSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, IndexSetSuiteData );
   pcu_suite_setFixtures( suite, IndexSetSuite_Setup, IndexSetSuite_Teardown );
   pcu_suite_addTest( suite, IndexSetSuite_TestRemoveAll );
   pcu_suite_addTest( suite, IndexSetSuite_TestInsertion );
   pcu_suite_addTest( suite, IndexSetSuite_TestRemoval );
   pcu_suite_addTest( suite, IndexSetSuite_TestUpdateMembersCount );
   pcu_suite_addTest( suite, IndexSetSuite_TestGetMembers );
   pcu_suite_addTest( suite, IndexSetSuite_TestGetVacancies );
}
