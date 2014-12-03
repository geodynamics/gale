#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include "StgDomain/StgDomain.h"

#include "LibDiscretisationSuite.h"

typedef struct {
} LibDiscretisationSuiteData;

void LibDiscretisationSuite_Setup( LibDiscretisationSuiteData* data ) {
}

void LibDiscretisationSuite_Teardown( LibDiscretisationSuiteData* data ) {
}

void LibDiscretisationSuite_DirectoryStGermain( LibDiscretisationSuiteData* data ) {
	Stg_Object* testDirectoryStGermain;

	testDirectoryStGermain = (Stg_Object*)Stg_ObjectList_Get( Project_XMLSearchPaths, (Name)"StGermain" );
	pcu_check_true( testDirectoryStGermain != NULL );
}

void LibDiscretisationSuite_DirectoryDiscretisation( LibDiscretisationSuiteData * data  ) {
	Stg_Object* testDirectoryStGermain;
	Stg_Object* testDirectoryDiscretisation;

	testDirectoryStGermain = (Stg_Object*)Stg_ObjectList_Get( Project_XMLSearchPaths, (Name)"StGermain"  );
	testDirectoryDiscretisation = (Stg_Object*)Stg_ObjectList_Get( Project_XMLSearchPaths, (Name)"StgDomain" );

	pcu_check_true( ( strcmp((char* )LIB_DIR, (char*)testDirectoryStGermain) ) || ( testDirectoryDiscretisation != NULL ) );
}

void LibDiscretisationSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, LibDiscretisationSuiteData );
	pcu_suite_setFixtures( suite, LibDiscretisationSuite_Setup, LibDiscretisationSuite_Teardown);

	pcu_suite_addTest( suite, LibDiscretisationSuite_DirectoryStGermain );
	pcu_suite_addTest( suite, LibDiscretisationSuite_DirectoryDiscretisation );
}


