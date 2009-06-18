#ifndef StGermain_LinkedListSuite_h
#define StGermain_LinkedListSuite_h

void LinkedListSuite( pcu_suite_t* suite );

int LinkedListSuite_CompareFunction(void *data1, void *data2);
void LinkedListSuite_DataPrintFunction( void *nodeData, void *args );
void LinkedListSuite_DataCopyFunction( void **nodeData, void *newData, SizeT dataSize);

#endif
