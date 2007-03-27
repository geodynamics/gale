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
** $Id:  $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <mpi.h>
#include <StGermain/StGermain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/LinearAlgebra/LinearAlgebra.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"

#include "units.h"
#include "types.h"
#include "shortcuts.h"
#include "MGCoarsener.h"


/* Textual name of this class */
const Type MGCoarsener_Type = "MGCoarsener";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

MGCoarsener* _MGCoarsener_New(
		SizeT					_sizeOfSelf, 
		Type						type,
		Stg_Class_DeleteFunction*	_delete,
		Stg_Class_PrintFunction*		_print, 
		Stg_Class_CopyFunction*		_copy, 
		MGCoarsener_SetMeshFunc*		_setMesh, 
		MGCoarsener_CoarsenFunc*		_coarsen )
{
	MGCoarsener*	self;
	
	/* Allocate memory. */
	self = (MGCoarsener*)_Stg_Class_New(
		_sizeOfSelf,
		type,
		_delete,
		_print, 
		_copy );
	
	/* General info */
	
	/* Virtual info */
	self->_setMesh = _setMesh;
	self->_coarsen = _coarsen;
	
	/* MGCoarsener info */
	_MGCoarsener_Init( self );
	
	return self;
}


void _MGCoarsener_Init( MGCoarsener* self ) {
	/* General and Virtual info should already be set */
	
	/* MGCoarsener info */
	self->mesh = NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _MGCoarsener_Delete( void* coarsener ) {
	MGCoarsener*	self = (MGCoarsener*)coarsener;
	
	/* Delete the class itself */
	
	/* Delete parent */
	_Stg_Class_Delete( self );
}


void _MGCoarsener_Print( void* coarsener, Stream* stream ) {
	MGCoarsener*	self = (MGCoarsener*)coarsener;
	Stream*	myStream;
	
	/* Set the Journal for printing informations */
	myStream = Journal_Register( InfoStream_Type, "MGCoarsenerStream" );
	
	/* Print parent */
	_Stg_Class_Print( self, stream );
	
	/* General info */
	Journal_Printf( myStream, "MGCoarsener (ptr): (%p)\n", self );
	
	/* Virtual info */
	
	/* MGCoarsener info */
}


void* _MGCoarsener_Copy( void* coarsener, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
#if 0
	MGCoarsener*	self = (MGCoarsener*)coarsener;
	MGCoarsener*	newMGCoarsener;
	PtrMap*	map = ptrMap;
	Bool		ownMap = False;
	
	if( !map ) {
		map = PtrMap_New( 10 );
		ownMap = True;
	}
	
	/* Copy goes here. */
	
	if( ownMap ) {
		Stg_Class_Delete( map );
	}
	
	return (void*)newMGCoarsener;
#endif
	return NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void _MGCoarsener_SetMesh( void* coarsener, void* mesh ) {
	((MGCoarsener*)coarsener)->mesh = (Mesh*)mesh;
}

