/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005, Monash Cluster Computing 
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
** 		* Redistributions of source code must retain the above copyright notice, 
** 			this list of conditions and the following disclaimer.
** 		* Redistributions in binary form must reproduce the above copyright 
**			notice, this list of conditions and the following disclaimer in the 
**			documentation and/or other materials provided with the distribution.
** 		* Neither the name of the Monash University nor the names of its contributors 
**			may be used to endorse or promote products derived from this software 
**			without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
** THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
** PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS 
** BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
** CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
** SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
** HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
** LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
** OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
**
** Contact:
*%		Louis Moresi - Louis.Moresi@sci.monash.edu.au
*%
** Contributors:
*+		Robert Turnbull
*+		Vincent Lemiale
*+		Louis Moresi
*+		David May
*+		David Stegman
*+		Mirko Velic
*+		Patrick Sunter
*+		Julian Giordani
*+
** $Id: Context.c 725 2008-05-08 05:15:45Z WendySharples $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Rheology/Rheology.h>
#include <Underworld/Rheology/Rheology.h>

#include "types.h"
#include "UnderworldContext.h"
#include "XDMFGenerator.h"

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type UnderworldContext_Type = "UnderworldContext";

/* Constructors ------------------------------------------------------------------------------------------------*/

UnderworldContext* UnderworldContext_New( 
	Name			name,
	double		start,
	double		stop,
	MPI_Comm		communicator,
	Dictionary*	dictionary )
{
	UnderworldContext* self = _UnderworldContext_DefaultNew( name );

	self->isConstructed = True;
	_AbstractContext_Init( (AbstractContext*) self );
	_DomainContext_Init( (DomainContext*) self );	
	_FiniteElementContext_Init( (FiniteElementContext*) self );
	_PICelleratorContext_Init( (PICelleratorContext*) self );
	_UnderworldContext_Init( self );

	return self;
}	

void* _UnderworldContext_DefaultNew( Name name ) {
	return (void*) _UnderworldContext_New(
		sizeof(UnderworldContext),
		UnderworldContext_Type,
		_UnderworldContext_Delete,
		_UnderworldContext_Print,
		NULL, 
		_UnderworldContext_DefaultNew,
		_UnderworldContext_AssignFromXML,
		_FiniteElementContext_Build,
		_FiniteElementContext_Initialise,
		_AbstractContext_Execute,
		_UnderworldContext_Destroy,
		name,
		NON_GLOBAL,
		_UnderworldContext_SetDt,
		0,
		0,
		MPI_COMM_WORLD,
		NULL );
}

UnderworldContext* _UnderworldContext_New( UNDERWORLDCONTEXT_DEFARGS ) {
	UnderworldContext* self;
	
	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. 
		At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	self = (UnderworldContext*)_PICelleratorContext_New( PICELLERATORCONTEXT_PASSARGS );
	
	/* General info */
	
	/* Function pointers for this class that are not on the parent class should be set here */
	
	return self;
}

void _UnderworldContext_Init( UnderworldContext* self ) {
	self->isConstructed = True;
	self->Vrms = 0.0;
   self->timeIntegrator = NULL;
   self->stokesSLE = NULL;
   self->energySLE = NULL;
   self->compositionSLE = NULL;
   self->constitutiveMatrix = NULL;
   /* turn off because the context no longer holds things
	EntryPoint_Append_AlwaysLast( Context_GetEntryPoint( self, AbstractContext_EP_AssignFromXML ),
		"Underworld App Assign Pointers",
		UnderworldContext_AssignPointers,
		"Underworld_App_Construct" );
      */
}

/* Virtual Functions -------------------------------------------------------------------------------------------------------------*/

void _UnderworldContext_AssignFromXML( void* context, Stg_ComponentFactory* cf, void* data ) {
	UnderworldContext* self = (UnderworldContext*)context;

	_PICelleratorContext_AssignFromXML( context, cf, data );

	/* always generate XDMF files when we generate HDF5 checkpoints */
#ifdef WRITE_HDF5
	if( Dictionary_Entry_Value_AsBool( Dictionary_GetDefault( self->dictionary, "generateXDMF", Dictionary_Entry_Value_FromBool( True ) ) ) ){
		ContextEP_Append( self, AbstractContext_EP_Save, XDMFGenerator_GenerateAll );
		ContextEP_Append( self, AbstractContext_EP_DataSave, XDMFGenerator_GenerateAll );
	}
#endif
	_UnderworldContext_Init( self );
}

void _UnderworldContext_Delete( void* context ) {
	UnderworldContext* self = (UnderworldContext*)context;
	
	Journal_DPrintf( self->debug, "In: %s()\n", __func__ );

	/* Stg_Class_Delete parent */
	_PICelleratorContext_Delete( self );
}

void _UnderworldContext_Destroy( void* context ) {
	UnderworldContext* self = (UnderworldContext*)context;
	
	_PICelleratorContext_Destroy( self );
}

void _UnderworldContext_Print( void* context, Stream* stream ) {
	UnderworldContext* self = (UnderworldContext*)context;
	
	/* General info */
	Journal_Printf( (void*) stream, "UnderworldContext (ptr): %p\n", self );
	
	/* Print parent */
	_PICelleratorContext_Print( self, stream );

	Journal_PrintPointer( stream, self->stokesSLE );
	Journal_PrintPointer( stream, self->energySLE );
	Journal_PrintPointer( stream, self->compositionSLE );
	Journal_PrintPointer( stream, self->constitutiveMatrix );
}


void _UnderworldContext_SetDt( void* context, double dt ) {
	UnderworldContext* self = (UnderworldContext*)context;
	
	self->dt = dt;
}


/* Public Functions ----------------------------------------------------------------------------------------------------*/


/* EntryPoint Hooks ----------------------------------------------------------------------------------------------------*/

void UnderworldContext_AssignPointers( void* context, void* ptrToContext ) {
	UnderworldContext*      self         = (UnderworldContext*) context;
	FieldVariable_Register* fV_Register  = self->fieldVariable_Register;
	
	Stream_IndentBranch( StgFEM_Debug );
	Journal_DPrintf( self->debug, "In: %s()\n", __func__ );

	if ( !self->CF )
		return;
	
	self->timeIntegrator = (TimeIntegrator*)  LiveComponentRegister_Get( self->CF->LCRegister, "timeIntegrator" );

	/* Get SLEs */
	self->stokesSLE = (Stokes_SLE*)            LiveComponentRegister_Get( self->CF->LCRegister, "stokesEqn" );
	self->energySLE = (AdvectionDiffusionSLE*) LiveComponentRegister_Get( self->CF->LCRegister, "EnergyEqn" );
	self->compositionSLE = (AdvectionDiffusionSLE*) LiveComponentRegister_Get( self->CF->LCRegister, "CompositionEqn" );
	self->constitutiveMatrix = (ConstitutiveMatrix*) LiveComponentRegister_Get( self->CF->LCRegister, "constitutiveMatrix" );
	
	Stream_UnIndentBranch( StgFEM_Debug );
}
