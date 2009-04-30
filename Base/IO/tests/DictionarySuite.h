#ifndef StGermain_DictionarySuite_h
#define StGermain_DictionarySuite_h

#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"

void DictionarySuite( pcu_suite_t* suite );

/* The following are useful for several other testing suites in this dir */

typedef struct {
   double startx;
   double starty;
   double startz;
} Geom;

typedef struct {
   double   height;
   Bool     anisotropic;
   char*    person;
   Geom     geom;
} TestStruct;

typedef struct {
   Index                      testEntriesCount;
   char**                     testKeys;
   Dictionary_Entry_Value**   testValues;
   char*                      testString;
   double                     testDouble;
   unsigned int               testUint;
   int                        testInt;
   unsigned long              testUnsignedlong;
   Bool                       testBool;
   double*                    testList;
   Index                      testListCount;
   TestStruct*                testStruct;
} DictionarySuite_TestDictData;

void DictionarySuite_SetupTestDictData( DictionarySuite_TestDictData* testDD );
void DictionarySuite_DictionaryData_Free( DictionarySuite_TestDictData* testDD );
void DictionarySuite_PopulateDictWithTestValues( Dictionary* dict, DictionarySuite_TestDictData* testDD );

#endif
