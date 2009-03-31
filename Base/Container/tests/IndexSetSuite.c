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
#define IS2_SIZE   40

typedef struct {
      IndexSet* is;
      IndexSet* is2;
} IndexSetSuiteData;


void IndexSetSuite_Setup( IndexSetSuiteData* data ) {
   /* Index sets are deliberately different sizes, to test all aspects of merge functionality */
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


void IndexSetSuite_TestAddAll( IndexSetSuiteData* data ) {
   Index ii;

   /* Index set will be starting in empty state */
   IndexSet_AddAll( data->is );
   for ( ii=0; ii< IS1_SIZE; ii++ ) {
      pcu_assert_true( True == IndexSet_IsMember( data->is, ii ) );
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

   /*Now remove a couple of these again */
   IndexSet_Remove( data->is, 7 );
   IndexSet_Remove( data->is, 24 );

   /*Thus list of members should be (0,7,8,22,24) - (7,24) = (0,8,22) */

   for ( ii=0; ii< IS1_SIZE; ii++ ) {
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

   /* Add some members, to generate a count */
   for ( ii=0; ii < 5; ii++ ) {
      IndexSet_Add( data->is, ii );
   }
   for ( ii=20; ii < 25; ii++ ) {
      IndexSet_Add( data->is, ii );
   }
   IndexSet_UpdateMembersCount( data->is );
   pcu_assert_true( 10 == data->is->membersCount );

}


void IndexSetSuite_TestGetIndexOfNthMember( IndexSetSuiteData* data ) {
   Index ii;
   Index retreivedIndex;
   #define TEST_INDICES_LENGTH 5
   Index testIndices[TEST_INDICES_LENGTH] = {0,3,9,11,22};

   /* add all the test indices */
   for( ii = 0; ii < TEST_INDICES_LENGTH; ii++ ) {
      IndexSet_Add( data->is, testIndices[ii] );
   }

   IndexSet_UpdateMembersCount( data->is );
   /* now check we are able to recover their indices correctly */
   for( ii = 0; ii < data->is->membersCount; ii++ ) {
      retreivedIndex = IndexSet_GetIndexOfNthMember( data->is, ii );
      pcu_assert_true( testIndices[ii] == retreivedIndex );
   }
}


void IndexSetSuite_TestGetMembers( IndexSetSuiteData* data ) {
   Index*	setArray;
   unsigned int	setArraySize;

   IndexSet_RemoveAll( data->is );

   /* add some test set members */
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

   /* set all indices to be included, except those exactly divisible by an example number */
   for ( ii=0; ii< IS1_SIZE; ii++ ) {
      if ( 0 != (ii % 5) ) {
         IndexSet_Add( data->is, ii );
      } 
   }

   IndexSet_GetVacancies( data->is, &setArraySize, &setArray );

   pcu_assert_true( (IS1_SIZE/5) == setArraySize );
   for ( ii=0; ii < IS1_SIZE/5; ii++ ) {
      pcu_assert_true( ii*5 == setArray[ii] );
   }

   Memory_Free( setArray );
}


void IndexSetSuite_TestDuplicate( IndexSetSuiteData* data ) {
   IndexSet* isDup;
   Index     ii;

   /* add some test set members to is1 for the duplicate */
   IndexSet_RemoveAll( data->is );
   IndexSet_Add( data->is, 0 );
   IndexSet_Add( data->is, 7 );
   IndexSet_Add( data->is, 8 );
   IndexSet_Add( data->is, 22 );
   IndexSet_UpdateMembersCount( data->is );
   
   isDup = IndexSet_Duplicate( data->is );

   /* check it's exactly the same */
   pcu_assert_true( isDup->size == data->is->size );
   pcu_assert_true( isDup->membersCount == data->is->membersCount );

   for ( ii=0; ii< IS1_SIZE; ii++ ) {
      pcu_assert_true( IndexSet_IsMember( isDup, ii ) == IndexSet_IsMember( data->is, ii ) );
   }
}


void IndexSetSuite_TestMerge_OR( IndexSetSuiteData* data ) {
   Index     ii;
   #define   TEST_INDICES1_LENGTH 5
   Index     testIndices1[TEST_INDICES1_LENGTH] = {0,3,9,11,   22};
   /* Note is2 is slightly larger, so we can add some test data outside range of is1 */
   #define   TEST_INDICES2_LENGTH 6
   Index     testIndices2[TEST_INDICES2_LENGTH] = {0,  9,   14,22,28,38};
   IndexSet* is1_preMerge = NULL;  
   IndexSet* is2_preMerge = NULL; 

   /* Create two copies of the index sets to use for comparison later */
   is1_preMerge = IndexSet_New(IS1_SIZE);
   is2_preMerge = IndexSet_New(IS2_SIZE);

   /* set up some test data. */
   for ( ii=0; ii< TEST_INDICES1_LENGTH; ii++ ) {
      IndexSet_Add( data->is, testIndices1[ii] );
      IndexSet_Add( is1_preMerge, testIndices1[ii] );
   }
   for ( ii=0; ii< TEST_INDICES2_LENGTH; ii++ ) {
      IndexSet_Add( data->is2, testIndices2[ii] );
      IndexSet_Add( is2_preMerge, testIndices2[ii] );
   }

   /* at this point, data->is and is1_preMerge should be identical, and likewise for is2_preMerge
   ** and testIndices2
   */
   
   IndexSet_Merge_OR( data->is, data->is2 );

   /* after the merge, we now expect data->is to hold the merged-or list */

   for ( ii=0; ii< IS1_SIZE; ii++ ) {
      /* check that an OR of the original list gives the same result as in merged list */
      if ( IndexSet_IsMember( is1_preMerge, ii ) || IndexSet_IsMember( is2_preMerge, ii ) ) {
         pcu_assert_true( IndexSet_IsMember( data->is, ii ) );
      }
   }
   /* merge function isn't supposed to modify its second argument. Check this. */
   for ( ii=0; ii< IS2_SIZE; ii++ ) {
      pcu_assert_true( IndexSet_IsMember( data->is2, ii ) == IndexSet_IsMember( is2_preMerge, ii ) );
   }

   /* ok - now try reverse operation, merging into is2. Result up to IS1_SIZE should be
      identical,
      but after that is2 list should be unchanged. */

   /* first re-set is1 */
   IndexSet_RemoveAll( data->is );
   for ( ii=0; ii< TEST_INDICES1_LENGTH; ii++ ) {
      IndexSet_Add( data->is, testIndices1[ii] );
   }

   /* then merge again, arguments reversed */
   IndexSet_Merge_OR( data->is2, data->is );
   /* and check is2 is merged for common length */
   for ( ii=0; ii< IS1_SIZE; ii++ ) {
      /* check that an OR of the original list gives the same result as in merged list */
      if ( IndexSet_IsMember( is1_preMerge, ii ) || IndexSet_IsMember( is2_preMerge, ii ) ) {
         pcu_assert_true( IndexSet_IsMember( data->is2, ii ) );
      }
   }
   /* check unchanged for last part of is2 */
   for ( ii=IS1_SIZE; ii< IS2_SIZE; ii++ ) {
      pcu_assert_true( IndexSet_IsMember( data->is2, ii ) == IndexSet_IsMember( is2_preMerge, ii ) );
   }
   /* merge function isn't supposed to modify its second argument. Check this. */
   for ( ii=0; ii< IS1_SIZE; ii++ ) {
      pcu_assert_true( IndexSet_IsMember( data->is, ii ) == IndexSet_IsMember( is1_preMerge, ii ) );
   }
}


void IndexSetSuite_TestMerge_AND( IndexSetSuiteData* data ) {
   Index     ii;
   #define   TEST_INDICES1_LENGTH 5
   Index     testIndices1[TEST_INDICES1_LENGTH] = {0,3,9,11,   22};
   /* Note is2 is slightly larger, so we can add some test data outside range of is1 */
   #define   TEST_INDICES2_LENGTH 6
   Index     testIndices2[TEST_INDICES2_LENGTH] = {0,  9,   14,22,28,38};
   IndexSet* is1_preMerge = NULL;  
   IndexSet* is2_preMerge = NULL; 

   /* Create two copies of the index sets to use for comparison later */
   is1_preMerge = IndexSet_New(IS1_SIZE);
   is2_preMerge = IndexSet_New(IS2_SIZE);

   /* set up some test data.  */
   for ( ii=0; ii< TEST_INDICES1_LENGTH; ii++ ) {
      IndexSet_Add( data->is, testIndices1[ii] );
      IndexSet_Add( is1_preMerge, testIndices1[ii] );
   }
   for ( ii=0; ii< TEST_INDICES2_LENGTH; ii++ ) {
      IndexSet_Add( data->is2, testIndices2[ii] );
      IndexSet_Add( is2_preMerge, testIndices2[ii] );
   }

   /* at this point, data->is and is1_preMerge should be identical, and likewise
    *  for is2_preMerge and testIndices2 */
   
   IndexSet_Merge_AND( data->is, data->is2 );

   for ( ii=0; ii< IS1_SIZE; ii++ ) {
      /* check that an AND of the original list gives the same result as in merged list */
      if ( IndexSet_IsMember( is1_preMerge, ii ) && IndexSet_IsMember( is2_preMerge, ii ) ) {
         pcu_assert_true( IndexSet_IsMember( data->is, ii ) );
      }
   }
   /* merge function isn't supposed to modify its second argument. Check this. */
   for ( ii=0; ii< IS2_SIZE; ii++ ) {
      pcu_assert_true( IndexSet_IsMember( data->is2, ii ) == IndexSet_IsMember( is2_preMerge, ii ) );
   }

   /* ok - now try reverse operation, merging into is2. Result up to IS1_SIZE should be
    * identical, but after that is2 list should be unchanged.
    */

   /* first re-set is1 */
   IndexSet_RemoveAll( data->is );
   for ( ii=0; ii< TEST_INDICES1_LENGTH; ii++ ) {
      IndexSet_Add( data->is, testIndices1[ii] );
   }

   IndexSet_Merge_AND( data->is2, data->is );

   /* and check is2 is merged for common length */
   for ( ii=0; ii< IS1_SIZE; ii++ ) {
      /* check that an AND of the original list gives the same result as in merged list */
      if ( IndexSet_IsMember( is1_preMerge, ii ) && IndexSet_IsMember( is2_preMerge, ii ) ) {
         pcu_assert_true( IndexSet_IsMember( data->is2, ii ) );
      }
   }
   /* check last part of is2, beyond length of is1, has all entries set to 0 */
   for ( ii=IS1_SIZE; ii< IS2_SIZE; ii++ ) {
      pcu_assert_true( False == IndexSet_IsMember( data->is2, ii ) );
   }
   /* merge function isn't supposed to modify its second argument. Check this. */
   for ( ii=0; ii< IS1_SIZE; ii++ ) {
      pcu_assert_true( IndexSet_IsMember( data->is, ii ) == IndexSet_IsMember( is1_preMerge, ii ) );
   }
}


void IndexSetSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, IndexSetSuiteData );
   pcu_suite_setFixtures( suite, IndexSetSuite_Setup, IndexSetSuite_Teardown );
   pcu_suite_addTest( suite, IndexSetSuite_TestRemoveAll );
   pcu_suite_addTest( suite, IndexSetSuite_TestAddAll );
   pcu_suite_addTest( suite, IndexSetSuite_TestInsertion );
   pcu_suite_addTest( suite, IndexSetSuite_TestRemoval );
   pcu_suite_addTest( suite, IndexSetSuite_TestUpdateMembersCount );
   pcu_suite_addTest( suite, IndexSetSuite_TestGetIndexOfNthMember );
   pcu_suite_addTest( suite, IndexSetSuite_TestGetMembers );
   pcu_suite_addTest( suite, IndexSetSuite_TestGetVacancies );
   pcu_suite_addTest( suite, IndexSetSuite_TestDuplicate );
   pcu_suite_addTest( suite, IndexSetSuite_TestMerge_OR );
   pcu_suite_addTest( suite, IndexSetSuite_TestMerge_AND );
}
