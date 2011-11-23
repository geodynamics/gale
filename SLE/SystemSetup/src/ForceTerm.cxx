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
** $Id: ForceTerm.c 964 2007-10-11 08:03:06Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "StgFEM/Discretisation/Discretisation.h"

#include "units.h"
#include "types.h"
#include "shortcuts.h"
#include "FiniteElementContext.h"
#include "ForceTerm.h"
#include "SolutionVector.h"
#include "ForceVector.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "EntryPoint.h"

/* Textual name of this class */
const Type ForceTerm_Type = "ForceTerm";

ForceTerm* ForceTerm_New(
	Name							name,
	FiniteElementContext*	context,
	ForceVector*				forceVector,
	Swarm*						integrationSwarm,
	Stg_Component*				extraInfo )		
{
	ForceTerm* self = (ForceTerm*) _ForceTerm_DefaultNew( name );

	self->isConstructed = True;
	_ForceTerm_Init( self, context, forceVector, integrationSwarm, extraInfo );

	return self;
}

void* _ForceTerm_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(ForceTerm);
	Type                                                      type = ForceTerm_Type;
	Stg_Class_DeleteFunction*                              _delete = _ForceTerm_Delete;
	Stg_Class_PrintFunction*                                _print = _ForceTerm_Print;
	Stg_Class_CopyFunction*                                  _copy = _ForceTerm_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _ForceTerm_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _ForceTerm_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _ForceTerm_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _ForceTerm_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _ForceTerm_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _ForceTerm_Destroy;
	AllocationType                              nameAllocationType = NON_GLOBAL;
	ForceTerm_AssembleElementFunction*            _assembleElement = _ForceTerm_AssembleElement;

	return _ForceTerm_New(  FORCETERM_PASSARGS  );
}

ForceTerm* _ForceTerm_New(  FORCETERM_DEFARGS  ) {
	ForceTerm* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(ForceTerm) );
	self = (ForceTerm*)_Stg_Component_New(  STG_COMPONENT_PASSARGS  );

	self->_assembleElement = _assembleElement;
	
	return self;
}


void _ForceTerm_Init(
	void*							forceTerm,
	FiniteElementContext*	context,
	ForceVector*				forceVector,
	Swarm*						integrationSwarm,
	Stg_Component*				extraInfo )
{
	ForceTerm* self			= (ForceTerm*) forceTerm;
	self->context 				= context;	
	self->debug					= Stream_RegisterChild( StgFEM_SLE_SystemSetup_Debug, self->type );
	self->extraInfo			= extraInfo;
	self->integrationSwarm	= integrationSwarm;	
   self->forceVector       = forceVector;

	ForceVector_AddForceTerm( forceVector, self );
}


void _ForceTerm_Delete( void* forceTerm ) {
	ForceTerm* self = (ForceTerm*)forceTerm;
	
	Journal_DPrintf( self->debug, "In %s - for %s\n", __func__, self->name );

	/* Stg_Class_Delete parent*/
	_Stg_Component_Delete( self );
}


void _ForceTerm_Print( void* forceTerm, Stream* stream ) {
	ForceTerm* self = (ForceTerm*)forceTerm;

	/* General info */
	Journal_Printf( stream, "ForceTerm (ptr): %p\n", self );
	
	/* Print parent */
	_Stg_Component_Print( self, stream );
	
	/* ForceTerm info */
        if(self->integrationSwarm!=NULL)
          Journal_Printf( stream, "\tintegrationSwarm (ptr): %p\n", self->integrationSwarm );
	Journal_Printf( stream, "\textraInfo (ptr): %p\n", self->extraInfo );
}


void* _ForceTerm_Copy( const void* forceTerm, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	ForceTerm*	self = (ForceTerm*)forceTerm;
	ForceTerm*	newForceTerm;
	PtrMap*		map = ptrMap;
	Bool			ownMap = False;
	
	if( !map ) {
		map = PtrMap_New( 10 );
		ownMap = True;
	}
	
	newForceTerm = (ForceTerm*)_Stg_Component_Copy( self, dest, deep, nameExt, map );
	
	newForceTerm->extraInfo = self->extraInfo;
	if(deep && self->integrationSwarm!=NULL) {
		newForceTerm->integrationSwarm = (Swarm*)Stg_Class_Copy( self->integrationSwarm, NULL, deep, nameExt, map );
	}
	else {
		newForceTerm->integrationSwarm = self->integrationSwarm;
	}
	
	if( ownMap ) {
		Stg_Class_Delete( map );
	}
	
	return (void*)newForceTerm;
}

void _ForceTerm_AssignFromXML( void* forceTerm, Stg_ComponentFactory* cf, void* data ) {
	FiniteElementContext*	context;
	ForceTerm*					self = (ForceTerm*)forceTerm;
	Swarm*						swarm = NULL;
	Stg_Component*				extraInfo;
	ForceVector*				forceVector;

	context = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Context", FiniteElementContext, False, data );

	if( !context  )
		context = Stg_ComponentFactory_ConstructByName( cf, (Name)"context", FiniteElementContext, True, data  );

	forceVector = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"ForceVector", ForceVector, True, data  ) ;
	swarm       = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Swarm", Swarm, False, data  ) ;
	extraInfo   = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"ExtraInfo", Stg_Component, False, data  ) ;

	_ForceTerm_Init( self, context, forceVector, swarm, extraInfo );
}

void _ForceTerm_Build( void* forceTerm, void* data ) {
	ForceTerm* self = (ForceTerm*)forceTerm;
	
	Journal_DPrintf( self->debug, "In %s - for %s\n", __func__, self->name );
	Stream_IndentBranch( StgFEM_Debug );
	
	/* ensure integrationSwarm is built */
        if(self->integrationSwarm!=NULL)
          Stg_Component_Build( self->integrationSwarm, data, False );
	Stg_Component_Build( self->forceVector, data, False );

	if ( self->extraInfo ) 
		Stg_Component_Build( self->extraInfo, data, False );
		
	Stream_UnIndentBranch( StgFEM_Debug );
}


void _ForceTerm_Initialise( void* forceTerm, void* data ) {
	ForceTerm* self = (ForceTerm*)forceTerm;
	
	Journal_DPrintf( self->debug, "In %s - for %s\n", __func__, self->name );
	Stream_IndentBranch( StgFEM_Debug );

        if(self->integrationSwarm!=NULL)
          Stg_Component_Initialise( self->integrationSwarm, data, False );
	Stg_Component_Initialise( self->forceVector, data, False );
	if ( self->extraInfo ) 
		Stg_Component_Initialise( self->extraInfo, data, False );
	
	Stream_UnIndentBranch( StgFEM_Debug );
}

void _ForceTerm_Execute( void* forceTerm, void* data ) {
}

void _ForceTerm_Destroy( void* forceTerm, void* data ) {
  ForceTerm* self = (ForceTerm*)forceTerm;

  if ( self->extraInfo ) 
    Stg_Component_Destroy( self->extraInfo, data, False );

  if(self->integrationSwarm!=NULL)
    Stg_Component_Destroy( self->integrationSwarm, data, False );
  Stg_Component_Destroy( self->forceVector, data, False );
}

void ForceTerm_AssembleElement( 
			void*                             forceTerm, 
			ForceVector*                      forceVector, 
			Element_LocalIndex                lElement_I,
			double*                           elForceVecToAdd ) 
{
	ForceTerm* self = (ForceTerm*)forceTerm;

	self->_assembleElement( self, forceVector, lElement_I, elForceVecToAdd );
}
	
void _ForceTerm_AssembleElement( 
			void*                             forceTerm, 
			ForceVector*                      forceVector, 
			Element_LocalIndex                lElement_I,
			double*                           elForceVecToAdd ) 
{
	ForceTerm* self        = (ForceTerm*)forceTerm;
	Stream*    errorStream = Journal_Register( Error_Type, (Name)self->type  );

	Journal_Printf( errorStream, "Error in func %s for %s '%s' - "
			"This function is the default function which should never be called - "
			"Please set this virtual function with appropriate application dependent function.\n",
			__func__, self->type, self->name );
	abort();
}	

void ForceTerm_SetAssembleElementFunction( void* forceTerm, ForceTerm_AssembleElementFunction* assembleElementFunction ) {
	ForceTerm* self        = (ForceTerm*)forceTerm;

	self->_assembleElement = assembleElementFunction;
}


