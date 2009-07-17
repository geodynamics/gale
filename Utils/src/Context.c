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
#include "Context.h"
#include "XDMFGenerator.h"

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type UnderworldContext_Type = "UnderworldContext";

/* Constructors ------------------------------------------------------------------------------------------------*/
void* _UnderworldContext_DefaultNew( Name name ) {
	return (void*) _UnderworldContext_New(
		sizeof(UnderworldContext),
		UnderworldContext_Type,
		_UnderworldContext_Delete,
		_UnderworldContext_Print,
		NULL, 
		_UnderworldContext_DefaultNew,
		_AbstractContext_Construct,
		_AbstractContext_Build,
		_AbstractContext_Initialise,
		_AbstractContext_Execute,
		_AbstractContext_Destroy,
		_UnderworldContext_SetDt,
		name,
		False,
		0,
		0,
		0,
		NULL );
}

UnderworldContext* UnderworldContext_New( 
		Name                        name,
		double						start,
		double						stop,
		MPI_Comm					communicator,
		Dictionary*					dictionary )
{
	return _UnderworldContext_New(
		sizeof(UnderworldContext),
		UnderworldContext_Type,
		_UnderworldContext_Delete,
		_UnderworldContext_Print,
		NULL, 
		_UnderworldContext_DefaultNew,
		_AbstractContext_Construct,
		_AbstractContext_Build,
		_AbstractContext_Initialise,
		_AbstractContext_Execute,
		_AbstractContext_Destroy,
		_UnderworldContext_SetDt,
		name,
		True,
		start,
		stop,
		communicator,
		dictionary );
}	


UnderworldContext* _UnderworldContext_New( 
		SizeT                                  sizeOfSelf,
		Type                                   type,
		Stg_Class_DeleteFunction*              _delete,
		Stg_Class_PrintFunction*               _print,
		Stg_Class_CopyFunction*                _copy, 
		Stg_Component_DefaultConstructorFunction*  _defaultConstructor,
		Stg_Component_ConstructFunction*       _construct,
		Stg_Component_BuildFunction*           _build,
		Stg_Component_InitialiseFunction*      _initialise,
		Stg_Component_ExecuteFunction*         _execute,
		Stg_Component_DestroyFunction*         _destroy,
		AbstractContext_SetDt*                 _setDt,
		Name                                   name,
		Bool                                   initFlag,
		double                                 start,
		double                                 stop,
		MPI_Comm                               communicator,
		Dictionary*                            dictionary )		
{
	UnderworldContext* self;
	
	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	self = (UnderworldContext*)_PICelleratorContext_New( 
		sizeOfSelf, 
		type, 
		_delete, 
		_print, 
		_copy,
		_defaultConstructor,
		_construct,
		_build,
		_initialise,
		_execute,
		_destroy,
		_setDt, 
		name,
		initFlag,
		start, 
		stop, 
		communicator, 
		dictionary );
	
	/* General info */
	
	/* Function pointers for this class that are not on the parent class should be set here */
	
	if( initFlag ){
		_UnderworldContext_Init( self );
	}
	
	return self;
}


void _UnderworldContext_Init( UnderworldContext* self ) {
	self->isConstructed = True;
	self->Vrms = 0.0;
	EntryPoint_Append_AlwaysLast( Context_GetEntryPoint( self, AbstractContext_EP_Construct ),
			   "Underworld App Assign Pointers",
			   UnderworldContext_AssignPointers,
			   "Underworld_App_Construct" );

   /* always generate XDMF files when we generate HDF5 checkpoints */
#ifdef WRITE_HDF5
	if( Dictionary_Entry_Value_AsBool( Dictionary_GetDefault( self->dictionary, "generateXDMF", Dictionary_Entry_Value_FromBool( True ) ) ) )
      ContextEP_Append( self, AbstractContext_EP_Save, XDMFGenerator_GenerateAll );
#endif
   

}


/* Virtual Functions -------------------------------------------------------------------------------------------------------------*/

void _UnderworldContext_Delete( void* context ) {
	UnderworldContext* self = (UnderworldContext*)context;
	
	Journal_DPrintf( self->debug, "In: %s()\n", __func__ );

	/* Stg_Class_Delete parent */
	_PICelleratorContext_Delete( self );
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
	
	Journal_PrintPointer( stream, self->gaussSwarm );
	
	Journal_PrintPointer( stream, self->velocityField );
	Journal_PrintPointer( stream, self->velocityGradientsField );
	Journal_PrintPointer( stream, self->strainRateField );
	Journal_PrintPointer( stream, self->strainRateInvField );
	Journal_PrintPointer( stream, self->pressureField );
	Journal_PrintPointer( stream, self->temperatureField );
	Journal_PrintPointer( stream, self->temperatureGradientsField );
	Journal_PrintPointer( stream, self->viscosityField );
	Journal_PrintPointer( stream, self->stressField );
	Journal_PrintPointer( stream, self->stressInvField );
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
	
	/* Swarms */
	self->gaussSwarm     = (Swarm*) LiveComponentRegister_Get( self->CF->LCRegister, "gaussSwarm" );
	self->picIntegrationPoints = (IntegrationPointsSwarm*) LiveComponentRegister_Get( self->CF->LCRegister, "picIntegrationPoints" );

	/* Get copy of fields from register */
	self->velocityField             = (FeVariable*) FieldVariable_Register_GetByName( fV_Register, "VelocityField" );
	self->velocityGradientsField    = (FeVariable*) FieldVariable_Register_GetByName( fV_Register, "VelocityGradientsField" );
	self->strainRateField           = (FeVariable*) FieldVariable_Register_GetByName( fV_Register, "StrainRateField" );
	self->strainRateInvField        = (FeVariable*) FieldVariable_Register_GetByName( fV_Register, "StrainRateInvariantField" );
	
	self->pressureField             = (FeVariable*) FieldVariable_Register_GetByName( fV_Register, "PressureField" );
	
	self->temperatureField          = (FeVariable*) FieldVariable_Register_GetByName( fV_Register, "TemperatureField" );
	self->temperatureGradientsField = (FeVariable*) FieldVariable_Register_GetByName( fV_Register, "TemperatureGradientsField" );
	
	self->viscosityField            = (FeVariable*) FieldVariable_Register_GetByName( fV_Register, "ViscosityField" );
	self->stressField               = (FeVariable*) FieldVariable_Register_GetByName( fV_Register, "StressField" );
	self->stressInvField            = (FeVariable*) FieldVariable_Register_GetByName( fV_Register, "StressInvField" );

	Stream_UnIndentBranch( StgFEM_Debug );
}
