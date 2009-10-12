#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include "StgDomain/StgDomain.h"
#include <StgFEM/StgFEM.h>

#include "PICellerator/PICellerator.h"
#include "PICellerator/Init.h"
#include "PICellerator/Finalise.h"

#include "LibPICelleratorSuite.h"

typedef struct {
} LibPICelleratorSuiteData;

void LibPICelleratorSuite_Setup( LibPICelleratorSuiteData* data ) {
}

void LibPICelleratorSuite_Teardown( LibPICelleratorSuiteData* data ) {
}

void LibPICelleratorSuite_DirectoryStGermain( LibPICelleratorSuiteData* data ) {
    Stg_Object* testDirectoryStGermain;
    testDirectoryStGermain = Stg_ObjectList_Get( Project_XMLSearchPaths, "StGermain" );
    pcu_check_true( testDirectoryStGermain != NULL );
}

void LibPICelleratorSuite_DirectoryStgFEM( LibPICelleratorSuiteData * data ) {
    Stg_Object* testDirectoryStGermain;
    Stg_Object* testDirectoryStgFEM;

    testDirectoryStGermain = Stg_ObjectList_Get( Project_XMLSearchPaths, "StGermain" );
    testDirectoryStgFEM= Stg_ObjectList_Get( Project_XMLSearchPaths, "StgFEM" );

    pcu_check_true( ( strcmp((char*)LIB_DIR, (char*)testDirectoryStGermain) ) || ( testDirectoryStgFEM != NULL ) );
}


void LibPICelleratorSuite_DirectoryPICellerator( LibPICelleratorSuiteData * data ) {
    Stg_Object* testDirectoryStGermain;
    Stg_Object* testDirectoryPICellerator;

    testDirectoryStGermain = Stg_ObjectList_Get( Project_XMLSearchPaths, "StGermain" );
    testDirectoryPICellerator= Stg_ObjectList_Get( Project_XMLSearchPaths, "PICellerator" );

    pcu_check_true( ( strcmp((char*)LIB_DIR, (char*)testDirectoryStGermain) ) || ( testDirectoryPICellerator != NULL ) );
}

void LibPICelleratorSuite( pcu_suite_t* suite ) {

    pcu_suite_setData( suite, LibPICelleratorSuiteData );
    pcu_suite_setFixtures( suite, LibPICelleratorSuite_Setup, LibPICelleratorSuite_Teardown);

    pcu_suite_addTest( suite, LibPICelleratorSuite_DirectoryStGermain );
    pcu_suite_addTest( suite, LibPICelleratorSuite_DirectoryStgFEM);
    pcu_suite_addTest( suite, LibPICelleratorSuite_DirectoryPICellerator);
}
