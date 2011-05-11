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

#ifndef __StGermain_Base_Foundation_MemMonitor_h__
#define __StGermain_Base_Foundation_MemMonitor_h__

extern double Stg_MemoryWatchCriteria;

extern const Type Stg_MemMonitor_Type;
extern const Type Stg_MemMonitor_InfoStreamName;
extern const Type Stg_MemMonitor_TagType;

typedef struct {
	SizeT    memFinal;          /* memory used at end of monitoring period - note this is different from Peak memory */
	int      memDiff;           /* Is an int rather than SizeT, since it may be negative */
	double   avgProcMemDiff;	
	double   minProcMemDiff;
	double   maxProcMemDiff;
	double   percentChange;    /* Unlike the TimeMonitor, decided a percent change over the monitored period was better than
                                 percent of monitored period over total, as allocated memory can go up _and_ down */
	Bool     criterionPassed;
} MemMonitorData;


typedef struct {
	MemoryOpStamp m1;
	MemoryOpStamp m2;
	SizeT totalMem1;
	SizeT totalMem2;
	char* tag;
	Bool criteria;
	Bool print;
	MPI_Comm comm;
} Stg_MemMonitor;

void Stg_MemMonitor_Initialise();
void Stg_MemMonitor_Finalise();
/* Note: it's assumed this criteria should be considered passed for both +ve and -ve changes of greater magnitude than
 * ratio specified (thus same semantics as the percentChange calculated). */
void Stg_MemMonitor_SetMemoryWatchCriteria( double ratioOfTotalMemory );

Stg_MemMonitor* Stg_MemMonitor_New(Name tag, Bool criteria, Bool print, MPI_Comm comm );
void Stg_MemMonitor_Delete( Stg_MemMonitor* mm );

void Stg_MemMonitor_Begin( Stg_MemMonitor* mm );
double Stg_MemMonitor_End( Stg_MemMonitor* mm, MemMonitorData* mmData );

void Stg_MemMonitor_ConvertBytesToPrintingUnit( int bytesInput, double* convertedAmt, char* unitString );

#endif
