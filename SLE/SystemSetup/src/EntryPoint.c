/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**	Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
** Contributors:
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
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
** $Id: EntryPoint.c 964 2007-10-11 08:03:06Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "units.h"
#include "types.h"
#include "shortcuts.h"
#include "EntryPoint.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdarg.h>

/* Textual name of this class */
const Type FeEntryPoint_Type = "FeEntryPoint";


FeEntryPoint* FeEntryPoint_New( const Name name, unsigned int castType ) {
	/* Variables set in this function */
	SizeT                       _sizeOfSelf = sizeof(FeEntryPoint);
	Type                               type = FeEntryPoint_Type;
	Stg_Class_DeleteFunction*       _delete = _EntryPoint_Delete;
	Stg_Class_PrintFunction*         _print = _EntryPoint_Print;
	Stg_Class_CopyFunction*           _copy = _EntryPoint_Copy;
	EntryPoint_GetRunFunction*      _getRun = _FeEntryPoint_GetRun;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return _FeEntryPoint_New(  FEENTRYPOINT_PASSARGS  );
}

void FeEntryPoint_Init( void* feEntryPoint, Name name, unsigned int castType ) {
	FeEntryPoint* self = (FeEntryPoint*)feEntryPoint;
	
	/* General info */
	self->type = FeEntryPoint_Type;
	self->_sizeOfSelf = sizeof(FeEntryPoint);
	self->_deleteSelf = False;
	
	/* Virtual info */
	self->_delete = _EntryPoint_Delete;
	self->_print = _EntryPoint_Print;
	self->_copy = _EntryPoint_Copy;
	self->_getRun = _FeEntryPoint_GetRun;
	_Stg_Class_Init( (Stg_Class*)self );
	_Stg_Object_Init( (Stg_Object*)self, name, GLOBAL );
	_EntryPoint_Init( (EntryPoint*)self, castType );
	
	/* FeEntryPoint info */
	_FeEntryPoint_Init( self );
}

FeEntryPoint* _FeEntryPoint_New(  FEENTRYPOINT_DEFARGS  )
{
	FeEntryPoint* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(FeEntryPoint) );
	self = (FeEntryPoint*)_EntryPoint_New(  ENTRYPOINT_PASSARGS  );
	
	/* General info */
	
	/* Virtual info */
	
	/* FeEntryPoint info */
	_FeEntryPoint_Init( self );
	
	return self;
}

void _FeEntryPoint_Init( FeEntryPoint* self ) {
	/* General and Virtual info should already be set */
	
	/* FeEntryPoint info */
}


Func_Ptr _FeEntryPoint_GetRun( void* feEntryPoint ) {
	FeEntryPoint* self = (FeEntryPoint*)feEntryPoint;
	
	/* Most frequently called EPs are put first in the switch statement */
	switch( self->castType ) {
		case FeEntryPoint_AssembleStiffnessMatrix_CastType:
			return (Func_Ptr)_FeEntryPoint_Run_AssembleStiffnessMatrix;
		
		case FeEntryPoint_AssembleForceVector_CastType:
			return (Func_Ptr)_FeEntryPoint_Run_AssembleForceVector;
		
		default:
			return (Func_Ptr)_EntryPoint_GetRun( self );
	}
}


void _FeEntryPoint_Run_AssembleStiffnessMatrix( 
		void* feEntryPoint, 
		void* stiffnessMatrix, 
		Bool bcRemoveQuery,
		void* _sle,
		void* _context )
{
	FeEntryPoint* self = (FeEntryPoint*)feEntryPoint;
	Hook_Index hookIndex;
  double wallTime;
	
	#ifdef USE_PROFILE
		Stg_CallGraph_Push( stgCallGraph, _FeEntryPoint_Run_AssembleStiffnessMatrix, self->name );
	#endif

	for( hookIndex = 0; hookIndex < self->hooks->count; hookIndex++ ) {
    wallTime = MPI_Wtime();

		((FeEntryPoint_AssembleStiffnessMatrix_Function*)((Hook*)self->hooks->data[hookIndex])->funcPtr) (
			stiffnessMatrix,
			bcRemoveQuery,
			_sle,
			_context );

    stg_profile_EntryPoint( self->name, self->hooks->data[hookIndex]->name, MPI_Wtime() - wallTime );
	}

	#ifdef USE_PROFILE
		Stg_CallGraph_Pop( stgCallGraph );
	#endif
}

void _FeEntryPoint_Run_AssembleForceVector( 
		void* feEntryPoint, 
		void* forceVector ) 
{
  FeEntryPoint* self = (FeEntryPoint*)feEntryPoint;
  Hook_Index hookIndex;
  double wallTime;
	
	#ifdef USE_PROFILE
		Stg_CallGraph_Push( stgCallGraph, _FeEntryPoint_Run_AssembleForceVector, self->name );
	#endif

	for( hookIndex = 0; hookIndex < self->hooks->count; hookIndex++ ) {
    wallTime = MPI_Wtime();

    ((FeEntryPoint_AssembleForceVector_Function*)((Hook*)self->hooks->data[hookIndex])->funcPtr) (
      forceVector );

    stg_profile_EntryPoint( self->name, self->hooks->data[hookIndex]->name, MPI_Wtime() - wallTime );
	}

	#ifdef USE_PROFILE
		Stg_CallGraph_Pop( stgCallGraph );
	#endif
}




