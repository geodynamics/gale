/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
** $Id: MemoryReport.c 3803 2006-09-27 03:17:12Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include "types.h"
#include "forwardDecl.h"

#include "MemoryField.h"
#include "MemoryPointer.h"
#include "MemoryReport.h"
#include "Memory.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stddef.h>

const Type MemoryReport_Type = "MemoryReport";

const int MEMORYREPORT_DELTA = 1;	/**< Number of items to grow by when array resizes. */
const int MEMORYREPORT_SIZE = 2;	/**< Number of items an array begins with. */


/** Returns the index of a given group in an array and -1 if not found. */
int MemoryReport_Find_Group( int numGroups, MemoryReportGroup* groups, MemoryReportGroup search );


MemoryReport* MemoryReport_New( ) {
	MemoryReport* result = (MemoryReport*) malloc( sizeof(MemoryReport) );
	
	_MemoryReport_Init( result );
	return result;
}
	

void _MemoryReport_Init( MemoryReport* memoryReport )
{
	char     reportQueryName[1000];
	Index    ii=0;

	memoryReport->groupCount = 0;
	memoryReport->groupSize = MEMORYREPORT_SIZE;
	memoryReport->groups = (MemoryReportGroup*) malloc( sizeof(MemoryReportGroup) * MEMORYREPORT_SIZE );
	memoryReport->conditionCount = 0;
	memoryReport->conditionSize = MEMORYREPORT_SIZE;
	memoryReport->conditionGroups = (MemoryReportGroup*) malloc( sizeof(MemoryReportGroup) * MEMORYREPORT_SIZE );
	memoryReport->conditionValues = (char**) malloc( sizeof(char*) * MEMORYREPORT_SIZE );
	for ( ii=0; ii < MEMORYREPORT_SIZE; ii++ ) {
		memoryReport->conditionValues[ii] = NULL;
	}

	
	memoryReport->reportField = MemoryField_New( "Report Query:" );
}
	

void MemoryReport_Delete( MemoryReport* memoryReport )
{
	Index i;
	
	MemoryField_Delete( memoryReport->reportField );

	free( memoryReport->groups );
	free( memoryReport->conditionGroups );
	
	for ( i = 0; i < memoryReport->conditionCount; ++i )
	{
		if ( memoryReport->conditionValues[i] != NULL )
		{
	 		free( memoryReport->conditionValues[i] );
	 	}
	}
	free( memoryReport->conditionValues );
	free( memoryReport );
}


void MemoryReport_AddGroup( MemoryReport* memoryReport, MemoryReportGroup group )
{
	if ( MemoryReport_Find_Group( memoryReport->groupCount, memoryReport->groups, group ) >= 0 ) {
		return;
	}

	/* Extend the groups array if needed. */
	if ( memoryReport->groupCount == memoryReport->groupSize ) {
		memoryReport->groupSize += MEMORYREPORT_DELTA;
		memoryReport->groups = (MemoryReportGroup*)
			realloc( memoryReport->groups, sizeof(MemoryReportGroup) * memoryReport->groupSize );	
	}
	
	memoryReport->groups[memoryReport->groupCount] = group;
	memoryReport->groupCount++;
}

void MemoryReport_AddCondition( MemoryReport* memoryReport, MemoryReportGroup group, const char* condition )
{
	/* Add this group if it does not already exist. */
	if ( MemoryReport_Find_Group( memoryReport->groupCount, memoryReport->groups, group ) < 0 ) {
		MemoryReport_AddGroup( memoryReport, group );	
	}
	
	/* Extend the condition arrays if needed. */
	if ( memoryReport->conditionCount == memoryReport->conditionSize ) {
		memoryReport->conditionSize += MEMORYREPORT_DELTA;
		memoryReport->conditionGroups = (MemoryReportGroup*)
			realloc( memoryReport->conditionGroups, sizeof(MemoryReportGroup) * memoryReport->conditionSize );
		memoryReport->conditionValues = (char**)
			realloc( memoryReport->conditionValues, sizeof(char*) * memoryReport->conditionSize );		
	}
	
	memoryReport->conditionGroups[memoryReport->conditionCount] = group;
	
	if ( condition ) {
		char*	ptr = (char*)malloc( (strlen(condition) + 1) * sizeof(char) );
		strcpy( ptr, condition );
		memoryReport->conditionValues[memoryReport->conditionCount] = ptr;
	}
	else {
		/* NULL is a condition as well, such as Type_Invalid and Name_Invalid. */
		memoryReport->conditionValues[memoryReport->conditionCount] = NULL;
	}
	
	memoryReport->conditionCount++;
}


void MemoryReport_Print( void* memoryReport )
{
	MemoryReport*  self = (MemoryReport*) memoryReport;

	if ( self->groupCount == 0 ) {
		return;
	}
	

	/* Algorithm:
	 * - Iterate through all MemoryPointers recorded.
	 * - Tuples matching the condition are added to the results.
	 * - Statistics are derived from tuples.
	 *
	 * Reason:
	 * Allows flexibility to produce any report required. The down side is that the peak bytes used cannot be derived this way.
	 *
	 * The alternative is to always record stats for all combinations (useful ones) but that will have a large impact on run
	 * time as well as memory space.
	 */
		
	/* Derive the statistics, using a BTree parse */
	BTree_ParseTree( stgMemory->pointers, MemoryReport_Print_Helper, self );
	
	//prevField = rootField;
	//while ( prevField->subCount == 1 ) {
	//	Journal_Printf( stgMemory->infoStream, "%s \n", prevField->value );
	//	prevField = prevField->subFields[0];
	//}
	
	// TODO: replace reportField->value with a name representative of conditions
	MemoryField_PrintSummary( self->reportField, "~Report~", (MEMORYFIELD_ALL-MEMORYFIELD_PEAK) );
}


/* Used for a BTree parse to gather statistics */
void MemoryReport_Print_Helper( void *memoryPointer, void* memoryReport ) {
	MemoryPointer* memPtr = (MemoryPointer*) memoryPointer;
	MemoryReport*  memReport = (MemoryReport*)memoryReport;
	MemoryField*   subField = NULL;
	Bool           valid;             /* Whether a memory pointer record matches the conditions. */
	Index          iGroup, iCondition;/* Iterators. */
	const char*    valueStr = NULL;

	assert ( memPtr );

	/* check condition */
	valid = True;
	for ( iCondition = 0; iCondition < memReport->conditionCount && valid; ++iCondition ) {
		valueStr = _MemoryReport_GetValue( memReport, memReport->conditionGroups[iCondition], memPtr );
		if ( MemoryField_StringCompare( valueStr, memReport->conditionValues[iCondition] ) != 0 ) {
			valid = False;
		}
	}
		
	if ( valid ) {
		/* Add this entry, sorted by the groups of the report. */
		/* Start at the root field */
		subField = memReport->reportField;
		/*  The way MemoryReport is designed, keep adding "sub-fields" of more specialised info. */
		for ( iGroup = 0; iGroup < memReport->groupCount; ++iGroup ) {
			valueStr = _MemoryReport_GetValue( memReport, memReport->groups[iGroup], memPtr );
			subField = MemoryField_Register( subField, valueStr );
		}

		/* Add statistics for this entry - "subField" will now be most specialised info */
		subField->allocCount++;
		if ( memPtr->ptr == NULL ) {
			subField->freeCount++;
		}
		else {
		       subField->currentAllocation += memPtr->totalSize;
		}
		subField->totalAllocation += memPtr->totalSize;
	}
}


int MemoryReport_Find_Group( int numGroups, MemoryReportGroup* groups, MemoryReportGroup search )
{
	int result;
	
	for ( result = 0; result < numGroups; ++result )
	{
		if ( groups[result] == search )
		{
			return result;
		}
	}
	
	return -1;
}


const char* _MemoryReport_GetValue( MemoryReport* memoryReport, MemoryReportGroup reportGroup, MemoryPointer* memPtr ) {
	const char* valueString = NULL;

	switch ( reportGroup ) {
		case MEMORYREPORT_TYPE:
			valueString = memPtr->type->value;
			break;
		case MEMORYREPORT_NAME:
			valueString = memPtr->name->value;
			break;
		case MEMORYREPORT_FILE:
			valueString = memPtr->file->value;
			break;
		case MEMORYREPORT_FUNC:
			valueString = memPtr->func->value;
			break;
	}
	return valueString;
}
