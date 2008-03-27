#ifndef StGermain_ProgressSuite_h
#define StGermain_ProgressSuite_h

typedef struct ProgressSuiteData ProgressSuiteData;

void ProgressSuite_TestSetStream( struct ProgressSuiteData* data );
void ProgressSuite_TestSetTitle( struct ProgressSuiteData* data );
void ProgressSuite_TestSetPrefix( struct ProgressSuiteData* data );
void ProgressSuite_TestSetRange( struct ProgressSuiteData* data );
void ProgressSuite_TestCalcStatus( struct ProgressSuiteData* data );

void ProgressSuite( pcu_suite_t* suite );
void ProgressSuite_Setup( struct ProgressSuiteData* data );
void ProgressSuite_Teardown( struct ProgressSuiteData* data );

#endif
