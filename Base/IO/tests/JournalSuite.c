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
**   Tests the journal functionality
**
** $Id: testJournal.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "JournalSuite.h"


void JournalSuite_Setup( JournalSuiteData* data ) {
   /* This is where we'll save the Journal, to be restored after the test, and create one for testing */
   data->savedJournal = stJournal;
  	stJournal = Journal_New();
   Journal_SetupDefaultTypedStreams();

   /* For testing, we want our custom Journal to output to saved files */
   Stg_asprintf( &data->testStdOutFilename, "./testStdOut.txt" );
   Stg_asprintf( &data->testStdErrFilename, "./testStdErr.txt" );

   Stg_Class_Delete( stJournal->stdOut );
   Stg_Class_Delete( stJournal->stdErr );
   stJournal->stdOut = CFile_New();
   stJournal->stdErr = CFile_New();
   JournalFile_Open( stJournal->stdOut, data->testStdOutFilename );
   JournalFile_Open( stJournal->stdErr, data->testStdErrFilename );

   Stream_SetFile( Journal_GetTypedStream(Info_Type), stJournal->stdOut );
   Stream_SetFile( Journal_GetTypedStream(Debug_Type), stJournal->stdOut );
   Stream_SetFile( Journal_GetTypedStream(Dump_Type), stJournal->stdOut );
   Stream_SetFile( Journal_GetTypedStream(Error_Type), stJournal->stdErr );

   data->testStdOutFile = fopen( data->testStdOutFilename, "r" );
   data->testStdErrFile = fopen( data->testStdErrFilename, "r" );
}

void JournalSuite_Teardown( JournalSuiteData* data ) {
   /* Delete temporary Journal, then restore the regular one */
   Journal_Delete();
   stJournal = data->savedJournal;

   fclose( data->testStdOutFile );
   fclose( data->testStdErrFile );
   remove( data->testStdOutFilename );
   remove( data->testStdErrFilename );
}


void JournalSuite_TestRegister( JournalSuiteData* data ) {
   Stream* myInfo;
   Stream* myDebug;
   Stream* myDump;
   Stream* myError;
   Stream* allNew;   /* Will use for testing a non-standard type stream */
   
   myInfo = Journal_Register( Info_Type, "MyInfo" );
   myDebug = Journal_Register( Debug_Type, "MyDebug" );
   myDump = Journal_Register( Dump_Type, "MyDump" );
   myError = Journal_Register( Error_Type, "MyError" );
   allNew = Journal_Register( "New_Type", "allNew" );

   /* Check the streams themselves were created properly */
   pcu_check_true( 0 == strcmp( myInfo->name, "MyInfo" ) );
   pcu_check_true( myInfo->_parent == Journal_GetTypedStream( Info_Type ) );
   pcu_check_true( myInfo->_children->count == 0 );
   pcu_check_true( 0 == strcmp( myDebug->name, "MyDebug" ) );
   pcu_check_true( myDebug->_parent == Journal_GetTypedStream( Debug_Type ) );
   pcu_check_true( myDebug->_children->count == 0 );
   pcu_check_true( 0 == strcmp( myDump->name, "MyDump" ) );
   pcu_check_true( myDump->_parent == Journal_GetTypedStream( Dump_Type ) );
   pcu_check_true( myDump->_children->count == 0 );
   pcu_check_true( 0 == strcmp( myError->name, "MyError" ) );
   pcu_check_true( myError->_parent == Journal_GetTypedStream( Error_Type ) );
   pcu_check_true( myError->_children->count == 0 );
   pcu_check_true( 0 == strcmp( allNew->name, "allNew" ) );
   pcu_check_true( allNew->_parent == Journal_GetTypedStream( "New_Type" ) );
   pcu_check_true( allNew->_children->count == 0 );

   /* Now check they were inserted in Journal hierarchy correctly */
   pcu_check_true( Stg_ObjectList_Get( Journal_GetTypedStream(Info_Type)->_children, "MyInfo" ) == myInfo );
   pcu_check_true( Stg_ObjectList_Get( Journal_GetTypedStream(Debug_Type)->_children, "MyDebug" ) == myDebug );
   pcu_check_true( Stg_ObjectList_Get( Journal_GetTypedStream(Dump_Type)->_children, "MyDump" ) == myDump );
   pcu_check_true( Stg_ObjectList_Get( Journal_GetTypedStream(Error_Type)->_children, "MyError" ) == myError );
   pcu_check_true( Stg_ObjectList_Get( Journal_GetTypedStream("New_Type")->_children, "allNew" ) == allNew );
}


void JournalSuite_TestRegister2( JournalSuiteData* data ) {
   Stream* register2Stream;
   Stream* register2Test;
   
   register2Stream = Journal_Register2( Info_Type, "Component", "Instance" );
   register2Test   = Journal_Register(  Info_Type, "Component.Instance"    );

   pcu_check_true( register2Stream == register2Test );
}


void JournalSuite_TestPrintBasics( JournalSuiteData* data ) {
   Stream*     myInfo;
   Stream*     myDebug;
   Stream*     myDump;
   Stream*     myError;
   #define     MAXLINE 1000
   char        outLine[MAXLINE];

   /* Check as is expected - see Base/IO/src/Init.c . Important for later tests */
   pcu_check_true( Stream_IsEnable( Journal_GetTypedStream(Info_Type)) == True );
   pcu_check_true( Stream_IsEnable( Journal_GetTypedStream(Debug_Type)) == False );
   pcu_check_true( Stream_IsEnable( Journal_GetTypedStream(Dump_Type)) == False );
   pcu_check_true( Stream_IsEnable( Journal_GetTypedStream(Error_Type)) == True ) ;

   myInfo = Journal_Register( InfoStream_Type, "MyInfo");
   myDebug = Journal_Register( DebugStream_Type, "MyDebug" );
   myDump = Journal_Register( Dump_Type, "MyDump" );
   myError = Journal_Register( ErrorStream_Type, "MyError" );

   Journal_Printf( myInfo, "%s\n", "HELLOInfo" );
   pcu_check_true(         fgets( outLine, MAXLINE, data->testStdOutFile ));
   pcu_check_true( 0 == strcmp( outLine, "HELLOInfo\n" ));
   Journal_Printf( myDebug, "%s\n", "WORLDDebug" );
   pcu_check_true( NULL == fgets( outLine, MAXLINE, data->testStdErrFile ));
   Journal_Printf( myDump, "%s\n", "HELLODump" );
   pcu_check_true( NULL == fgets( outLine, MAXLINE, data->testStdOutFile ));
   Journal_Printf( myError, "%s\n", "WORLDError" );
   pcu_check_true(         fgets( outLine, MAXLINE, data->testStdErrFile ));
   pcu_check_true( 0 == strcmp( outLine, "WORLDError\n" ));

   Journal_Enable_NamedStream( Info_Type, "MyInfo", False );
   Journal_Printf( myInfo, "%s\n", "HELLOInfo2" );
   pcu_check_true( NULL == fgets( outLine, MAXLINE, data->testStdOutFile ));

   Journal_Enable_TypedStream( Dump_Type, True );
   Journal_Enable_NamedStream( Dump_Type, "MyDump", True );
   Journal_Printf( myDump, "%s\n", "HELLODump2" );
   /* This stream should have auto-flush set to false. Check first, then flush and check again */
   pcu_check_true( Journal_GetTypedStream(Dump_Type)->_autoFlush == False );
   pcu_check_true( NULL == fgets( outLine, MAXLINE, data->testStdOutFile ));
   Stream_Flush( Journal_GetTypedStream(Dump_Type) );
   pcu_check_true(         fgets( outLine, MAXLINE, data->testStdOutFile ));
   pcu_check_true( 0 == strcmp( outLine, "HELLODump2\n" ));
   
   stJournal->enable = False;
   Journal_Printf( myDump, "%s\n", "HELLODump3" );
   pcu_check_true( NULL == fgets( outLine, MAXLINE, data->testStdOutFile ));
}

/* TODO */
void JournalSuite_TestPrintfL( JournalSuiteData* data ) {
   Stream* myStream;

   myStream = Journal_Register( InfoStream_Type, "myComponent");
   Journal_PrintfL( myStream, 1, "Hello\n" );
   Journal_PrintfL( myStream, 2, "Hello\n" );
}


/* TODO */
void JournalSuite_TestPrintfD( JournalSuiteData* data ) {
   Stream* myInfo;

   Journal_DPrintf( myInfo, "DPrintf\n" );
}


/* TODO */
void JournalSuite_TestChildStreams( JournalSuiteData* data ) {
   Stream* myStream;

   Stream* childStream1;
   Stream* childStream2;

   myStream = Journal_Register( InfoStream_Type, "myComponent");

  /* Make sure the hierarchy works*/
   childStream1 = Stream_RegisterChild( myStream, "child1" );
   childStream2 = Stream_RegisterChild( childStream1, "child2" );

   Journal_Printf( myStream, "0 no indent\n" );
   Stream_IndentBranch( myStream );
   Journal_Printf( childStream1, "1 with 1 indent\n" );
   Stream_IndentBranch( myStream );
   Journal_Printf( childStream2, "2 with 2 indent\n" );
   Stream_UnIndentBranch( myStream );
   Journal_Printf( childStream2, "2 with 1 indent\n" );
   Stream_UnIndentBranch( myStream );
   Journal_Printf( childStream1, "1 with no indent\n" );
   Journal_Printf( childStream2, "2 with no indent\n" );
}


void JournalSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, JournalSuiteData );
   pcu_suite_setFixtures( suite, JournalSuite_Setup, JournalSuite_Teardown );
   pcu_suite_addTest( suite, JournalSuite_TestRegister );
   pcu_suite_addTest( suite, JournalSuite_TestRegister2 );
   pcu_suite_addTest( suite, JournalSuite_TestPrintBasics );
}
