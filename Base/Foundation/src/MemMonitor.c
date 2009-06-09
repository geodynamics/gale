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
** $Id: MemMonitor.c 3157 2005-08-07 23:43:05Z AlanLo $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdarg.h>
#include <mpi.h>

#include "pcu/pcu.h"
#include "types.h"
#include "shortcuts.h"
#include "forwardDecl.h"
#include "MemoryPointer.h"
#include "MemoryField.h"
#include "Memory.h"
#include "MemMonitor.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

double Stg_MemoryWatchCriteria = -1;

const Type Stg_MemMonitor_Type = "Stg_MemMonitor";
const Type Stg_MemMonitor_InfoStreamName = "MemMonitor";
const Type Stg_MemMonitor_TagType = "Stg_MemMonitor_Tag";

void _Memory_Print_AllocsAboveThreshold_Helper( void* memoryPointer, void* args ) {
	MemoryPointer* memPtr;
	void** arguments;
	SizeT threshold;
	MemoryOpStamp begin;
	MemoryOpStamp end;

	assert( memoryPointer );
	assert( args );

	memPtr = (MemoryPointer*)memoryPointer;
	arguments = (void**)args;
	threshold = *((SizeT*)arguments[0]);
	begin = *((MemoryOpStamp*)arguments[1]);
	end = *((MemoryOpStamp*)arguments[2]);

	if ( memPtr->ptr != NULL && memPtr->status != MEMORY_POINTER_RELEASED ) {
		if ( memPtr->totalSize > threshold && memPtr->stamp >= begin && memPtr->stamp <= end ) {
			MemoryPointer_Print( memPtr, MEMORYPOINTER_NAME | MEMORYPOINTER_TOTALSIZE );
		}
	}
}

void Stg_MemMonitor_Initialise() {
	Stg_MemoryWatchCriteria = 0.2;
}
void Stg_MemMonitor_Finalise() {

}
void Stg_MemMonitor_SetMemoryWatchCriteria( double ratioOfTotalMemory ) {
	Stg_MemoryWatchCriteria = ratioOfTotalMemory;
}

Stg_MemMonitor* Stg_MemMonitor_New( char* tag, Bool criteria, Bool print, int comm ) {
	Stg_MemMonitor* mm;
	
	mm = Memory_Alloc_Unnamed( Stg_MemMonitor );
	mm->tag = Memory_Alloc_Bytes_Unnamed( strlen( tag ) + 1, Stg_MemMonitor_TagType );
	strcpy( mm->tag, tag );
	mm->criteria = criteria;
	mm->print = print;
	mm->comm = comm;
	
	return mm;
}

void Stg_MemMonitor_Delete( Stg_MemMonitor* mm ) {
	if( mm->tag ) {
		Memory_Free( mm->tag );
	}
	
	Memory_Free( mm );
}


void Stg_MemMonitor_Begin( Stg_MemMonitor* mm ) {
#ifdef MEMORY_STATS
	mm->m1 = stgMemory->stamp;
	mm->m2 = mm->m1;
	MemoryField_UpdateAsSumOfSubFields( stgMemory->types );
	mm->totalMem1 = stgMemory->types->currentAllocation;
	mm->totalMem2 = mm->totalMem1;
#endif
}

double Stg_MemMonitor_End( Stg_MemMonitor* mm, MemMonitorData* mmData ) {
#ifdef MEMORY_STATS
	long memSumDiff;
	int rank;
	int size;
#endif

	mmData->avgProcMemDiff = 0;

#ifdef MEMORY_STATS
	
	mm->m2 = stgMemory->stamp;
	MemoryField_UpdateAsSumOfSubFields( stgMemory->types );

	mm->totalMem2 = stgMemory->types->currentAllocation;
	
	mmData->memFinal = mm->totalMem2;
	/* Casts are necessary below because both totalMems are SizeT, and the answer may be negative (IE, have net
	 *  free'd memory */
	mmData->memDiff = (int)((int)mm->totalMem2 - (int)mm->totalMem1);
	
	MPI_Comm_size( mm->comm, &size );

	/*
	MPI_Reduce( &mmData->memDiff, &maxProcMemDiff, 1, MPI_LONG, MPI_MAX, 0, mm->comm );
	MPI_Reduce( &mmData->memDiff, &minProcMemDiff, 1, MPI_LONG, MPI_MIN, 0, mm->comm );
	MPI_Allreduce( &mmData->memDiff, &memSumDiff, 1, MPI_LONG, MPI_SUM, mm->comm );
	
	avgProcMemDiff = (double)memSumDiff / size;
	*/
	/* Above is commented and replaced with below. See TimeMonitor.c for reason */
	memSumDiff = mmData->memDiff;

	mmData->minProcMemDiff = mmData->memDiff;
	mmData->maxProcMemDiff = mmData->memDiff;
	mmData->avgProcMemDiff = mmData->memDiff;
	
	/* Note: maybe Stg_Components should store rank and comm??? how do the find their comm? */
	
	MPI_Comm_rank( mm->comm, &rank );
	
	mmData->percentChange = (mmData->avgProcMemDiff / (double)mm->totalMem1) * 100.0;

	mmData->criterionPassed = False;
	if ( mm->criteria ) {
		mmData->criterionPassed = fabs( mmData->percentChange/100.0 ) >= Stg_MemoryWatchCriteria;
	}

	if( (rank == 0) && mm->print && ( (mm->criteria==False) || mmData->criterionPassed ) ) {
		void*    args[3];
		SizeT    threshold = (SizeT)(Stg_MemoryWatchCriteria * mmData->memFinal);

		if( size == 1 ) {
			double   memFinalPrint;
			char     memFinalUnit[100];
			double   memDiffPrint;
			char     memDiffUnit[100];

			Stg_MemMonitor_ConvertBytesToPrintingUnit( mmData->memFinal, &memFinalPrint, memFinalUnit );
			Stg_MemMonitor_ConvertBytesToPrintingUnit( mmData->memDiff, &memDiffPrint, memDiffUnit );

			Journal_Printf( 
				Journal_Register( Info_Type, Stg_MemMonitor_InfoStreamName ),
				"\t%s(%s): memory allocated at end= %.2f%s, of which %.2f%s during monitoring (%.2f%% change)\n", 
				Stg_MemMonitor_InfoStreamName,
				mm->tag,
				memFinalPrint, memFinalUnit,
				memDiffPrint, memDiffUnit,
				mmData->percentChange );
		}
		else {
			double   memFinalPrint;
			char     memFinalUnit[100];
			double   avgProcMemDiffPrint;
			char     avgProcMemDiffUnit[100];
			double   minProcMemDiffPrint;
			char     minProcMemDiffUnit[100];
			double   maxProcMemDiffPrint;
			char     maxProcMemDiffUnit[100];

			Stg_MemMonitor_ConvertBytesToPrintingUnit( mmData->memFinal, &memFinalPrint, memFinalUnit );
			Stg_MemMonitor_ConvertBytesToPrintingUnit( mmData->avgProcMemDiff, &avgProcMemDiffPrint, avgProcMemDiffUnit );
			Stg_MemMonitor_ConvertBytesToPrintingUnit( mmData->minProcMemDiff, &minProcMemDiffPrint, minProcMemDiffUnit );
			Stg_MemMonitor_ConvertBytesToPrintingUnit( mmData->maxProcMemDiff, &maxProcMemDiffPrint, maxProcMemDiffUnit );

			Journal_Printf( 
				Journal_Register( Info_Type, Stg_MemMonitor_InfoStreamName ),
				"\t%s(%s): memory allocated at end= %.2f%s, of which ave %.2f%s/proc during monitoring (%.2f%% change)\n"
				"\t\t(individual proc usage during monitoring min/max = %.2f%s/%.2f%s)\n",
				Stg_MemMonitor_InfoStreamName,
				mm->tag,
				memFinalPrint, memFinalUnit,
				avgProcMemDiffPrint, avgProcMemDiffUnit,
				mmData->percentChange,
				minProcMemDiffPrint, minProcMemDiffUnit,
				maxProcMemDiffPrint, maxProcMemDiffUnit );
		}
		args[0] = &threshold;
		args[1] = &mm->m1;
		args[2] = &mm->m2;
		Journal_Printf( Journal_Register( Info_Type, Stg_MemMonitor_InfoStreamName ), "Allocations after threshold reached were:\n\n" );
		BTree_ParseTree( stgMemory->pointers, _Memory_Print_AllocsAboveThreshold_Helper, (void*)args );
	}
#endif
	return mmData->avgProcMemDiff;
}


void Stg_MemMonitor_ConvertBytesToPrintingUnit( int bytesInput, double* convertedAmt, char* unitString ) {
	double    one_kb = 1024.0;
	double    one_mb = 1024.0*one_kb;
	double    one_gb = 1024.0*one_mb;
	double    one_tb = 1024.0*one_gb;
	unsigned  absBytesInput = abs(bytesInput);

	pcu_assert( unitString != NULL );

	/* Remember, bytesInput may be negative representing a free'd amount */

	if ( absBytesInput < one_kb) {
		/* leave as bytes */
		*convertedAmt = bytesInput;
		sprintf( unitString, "b" );
	}
	else if ( absBytesInput < one_mb ) {
		/* convert to KB */
		*convertedAmt = bytesInput / one_kb; 
		sprintf( unitString, "kb" );
	}
	else if ( absBytesInput < one_gb ) {
		/* convert to MB */
		*convertedAmt = bytesInput / one_mb;
		sprintf( unitString, "mb" );
	}
	else if ( absBytesInput < one_tb ) {
		/* convert to GB */
		*convertedAmt = bytesInput / one_gb;
		sprintf( unitString, "gb" );
	}
	else {
		/* convert to TB - don't think we need to worry about PetaBytes! */
		*convertedAmt = bytesInput / one_tb;
		sprintf( unitString, "tb" );
	}
}
