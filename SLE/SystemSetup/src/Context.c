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
** $Id: Context.c 833 2007-05-16 01:12:22Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/LinearAlgebra/LinearAlgebra.h"
#include "units.h"
#include "types.h"
#include "shortcuts.h"
#include "Context.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "SystemLinearEquations.h"
#include "SolutionVector.h"

/* Textual name of this class */
const Type FiniteElementContext_Type = "FiniteElementContext";
const Name defaultFiniteElementContextETypeRegisterName = "defaultFiniteElementContextETypeRegisterName";
const Name FiniteElementContext_EP_CalcDt = "FiniteElementContext_EP_CalcDt";


/* Constructors ------------------------------------------------------------------------------------------------------------------*/

void* FiniteElementContext_DefaultNew( Name name )
{
	return _FiniteElementContext_New(
		sizeof(FiniteElementContext),
		FiniteElementContext_Type,
		_FiniteElementContext_Delete,
		_FiniteElementContext_Print,
		NULL,
		FiniteElementContext_DefaultNew,
		_AbstractContext_Construct,
		_AbstractContext_Build,
		_AbstractContext_Initialise,
		_AbstractContext_Execute,
		_AbstractContext_Destroy,
		name,
		False,
		_FiniteElementContext_SetDt,
		0,
		0,
		MPI_COMM_WORLD,
		NULL );
}

FiniteElementContext*				FiniteElementContext_New(
		Name						name,
		double						start,
		double						stop,
		MPI_Comm					communicator,
		Dictionary*					dictionary )
{
	return _FiniteElementContext_New(
		sizeof(FiniteElementContext),
		FiniteElementContext_Type,
		_FiniteElementContext_Delete,
		_FiniteElementContext_Print,
		NULL,
		FiniteElementContext_DefaultNew,
		_AbstractContext_Construct,
		_AbstractContext_Build,
		_AbstractContext_Initialise,
		_AbstractContext_Execute,
		_AbstractContext_Destroy,
		name,
		True,
		_FiniteElementContext_SetDt,
		start,
		stop,
		communicator,
		dictionary );
}	


FiniteElementContext* _FiniteElementContext_New( 
		SizeT						sizeOfSelf,
		Type						type,
		Stg_Class_DeleteFunction*				_delete,
		Stg_Class_PrintFunction*				_print,
		Stg_Class_CopyFunction*				_copy, 
		Stg_Component_DefaultConstructorFunction*    _defaultConstructor,
		Stg_Component_ConstructFunction*_construct,
		Stg_Component_BuildFunction*    _build,
		Stg_Component_InitialiseFunction*   _initialise,
		Stg_Component_ExecuteFunction*  _execute,
		Stg_Component_DestroyFunction*  _destroy,
		Name 							name,
		Bool							initFlag,
		AbstractContext_SetDt*				_setDt,
		double						start,
		double						stop,
		MPI_Comm					communicator,
		Dictionary*					dictionary )
{
	FiniteElementContext* self;
	
	/* Allocate memory */
	self = (FiniteElementContext*)_DiscretisationContext_New( 
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
		name,
		initFlag,
		_setDt, 
		start, 
		stop, 
		communicator, 
		dictionary );
	
	/* General info */
	
	/* Virtual info */
	
	if( initFlag ){
		_FiniteElementContext_Init( self );
	}
	
	return self;
}


void _FiniteElementContext_Init( FiniteElementContext* self ) {
	Stream*  errorStream = Journal_Register( Error_Type, self->type );

	/* Set up stream preferences */
	Journal_Enable_NamedStream( InfoStream_Type, "StgFEM_VerboseConfig", False );


		

	/* register this current stream (the context) as a child of the FE stream */
	/* TODO: Want to be able to recombine this with the Abs context's debug stream at some stage */
	self->isConstructed = True;
	self->debug = Stream_RegisterChild( StgFEM_Debug, FiniteElementContext_Type );
	
	self->dt = 0.0f;
	self->prevTimestepDt = 0.0;
	self->limitTimeStepIncreaseRate = Dictionary_GetBool_WithDefault( self->dictionary, "limitTimeStepIncreaseRate", False );
	self->maxTimeStepIncreasePercentage = Dictionary_GetDouble_WithDefault( self->dictionary,
		"maxTimeStepIncreasePercentage", 10.0 );
	Journal_Firewall( self->maxTimeStepIncreasePercentage >= 0, errorStream,
		"Error - in %s(): maxTimeStepIncreasePercentage must be >= 0\n", __func__ );
	
	/* set up s.l.e list */
	self->slEquations = Stg_ObjectList_New();

	/* Set up the element type register and add basic types */
	self->elementType_Register = ElementType_Register_New(defaultFiniteElementContextETypeRegisterName);
	ElementType_Register_Add( self->elementType_Register, (ElementType*)ConstantElementType_New("constantElementType") );
	ElementType_Register_Add( self->elementType_Register, (ElementType*)BilinearElementType_New("bilinearElementType") );
	ElementType_Register_Add( self->elementType_Register, (ElementType*)TrilinearElementType_New("trilinearElementType") );
	ElementType_Register_Add( self->elementType_Register, (ElementType*)LinearTriangleElementType_New("linearElementType") );
	Stg_ObjectList_ClassAppend( self->register_Register, (void*)self->elementType_Register, "ElementType_Register" );

	/* Create Entry Point for Calculating timestep */
	self->calcDtEP = EntryPoint_New( FiniteElementContext_EP_CalcDt, EntryPoint_Minimum_VoidPtr_CastType );
	EntryPoint_Register_Add( self->entryPoint_Register, self->calcDtEP );
	
	/* Add hooks to existing entry points... use name "default" so that plugin, etc can exert same behaviour on other contexts*/
	EntryPoint_Prepend( 
		Context_GetEntryPoint( self, AbstractContext_EP_Build ),
		"default", 
		_FiniteElementContext_Build, 
		FiniteElementContext_Type );
	EntryPoint_Prepend( 
		Context_GetEntryPoint( self, AbstractContext_EP_Initialise ),
		"default", 
		_FiniteElementContext_Initialise, 
		FiniteElementContext_Type );
	EntryPoint_Append( 
		Context_GetEntryPoint( self, AbstractContext_EP_Solve ),
		"default", 
		_FiniteElementContext_Solve, 
		FiniteElementContext_Type );
	EntryPoint_Append( 
		Context_GetEntryPoint( self, AbstractContext_EP_Solve ),
		"postSolve", 
		_FiniteElementContext_PostSolve, 
		FiniteElementContext_Type );
	EntryPoint_Append( 
		Context_GetEntryPoint( self, AbstractContext_EP_Dt ),
		"default", 
		_FiniteElementContext_GetDt, 
		FiniteElementContext_Type );

	EntryPoint_Append(
		Context_GetEntryPoint( self, AbstractContext_EP_Save ),
		"saveFeVariables",
		_FiniteElementContext_SaveFeVariables,
		FiniteElementContext_Type );
	/* The FEM context needs to save gauss swarms so they can be re-loaded for restart later.
	   This will automatically save material point swarms too if PICellerator is used.
	 */
	EntryPoint_Append(
		Context_GetEntryPoint( self, AbstractContext_EP_Save ),
		"saveSwarms",
		_FiniteElementContext_SaveSwarms,
		FiniteElementContext_Type );

}


/* Virtual Functions -------------------------------------------------------------------------------------------------------------*/

void _FiniteElementContext_Delete( void* context ) {
	FiniteElementContext* self = (FiniteElementContext*)context;
	
	Journal_DPrintf( self->debug, "In: %s()\n", __func__ );

	Stream_IndentBranch( StgFEM_Debug );
	Journal_DPrintfL( self->debug, 2, "Deleting the element type register (and hence all element types).\n" );
	Stg_Class_Delete( self->elementType_Register );
	Journal_DPrintfL( self->debug, 2, "Deleting all SLEs and the SLE list.\n" );
	Stg_ObjectList_DeleteAllObjects( self->slEquations ); 
	Stg_Class_Delete( self->slEquations ); 
	Stream_UnIndentBranch( StgFEM_Debug );

	/* Stg_Class_Delete parent */
	_DiscretisationContext_Delete( self );
}


void _FiniteElementContext_Print( void* context, Stream* stream ) {
	FiniteElementContext* self = (FiniteElementContext*)context;
	
	/* General info */
	Journal_Printf( (void*) stream, "FiniteElementContext (ptr): %p\n", self );
	
	/* Print parent */
	_DiscretisationContext_Print( self, stream );

	Journal_Printf( (void*) stream, "\tslEquations (ptr): %p\n", self->slEquations );
	Stg_Class_Print( self->slEquations, stream );

	Journal_Printf( (void*) stream, "\telementType_Register (ptr): %p\n", self->elementType_Register );
	Stg_Class_Print( self->elementType_Register, stream );
}


void _FiniteElementContext_SetDt( void* context, double dt ) {
	FiniteElementContext* self = (FiniteElementContext*)context;
	
	self->dt = dt;
}


/* Public Functions --------------------------------------------------------------------------------------------------------------*/

void FiniteElementContext_AddSLE_Func( void* context, void* sle ) {
	FiniteElementContext* self = (FiniteElementContext*)context;
	
	FiniteElementContext_AddSLE_Macro( self, sle );
}


SystemLinearEquations* FiniteElementContext_GetSLE_Func( void* context, Name sleName ) {
	FiniteElementContext* self = (FiniteElementContext*)context;
	
	return FiniteElementContext_GetSLE_Macro( self, sleName );
}


/* EntryPoint Hooks --------------------------------------------------------------------------------------------------------------*/

void _FiniteElementContext_Construct( void* context, Stg_ComponentFactory* cf, void* data ){
	FiniteElementContext *self = (FiniteElementContext*) context;
	int componentIndex = 0;

	self->dictionary = cf->rootDict;

/* I don't think this is needed any more...
	_AbstractContext_Init( (AbstractContext*)self, 0, 0, MPI_COMM_WORLD );
	_DiscretisationContext_Init( (DiscretisationContext*)self );*/
	_FiniteElementContext_Init( self );
	
	componentIndex = Stg_ObjectList_GetIndex( cf->LCRegister->componentList, self->name );
	
	Journal_Firewall( 
		!componentIndex, 
		Journal_Register( Error_Type, self->type ), 
		"Context should be the first component in the 'components' list..!" );

	self->CF = cf;
	self->CF->registerRegister = self->register_Register;
}


void _FiniteElementContext_Build( void* context ) {
	FiniteElementContext* self = (FiniteElementContext*)context;
	SystemLinearEquations_Index sle_I;
	
	Stream_IndentBranch( StgFEM_Debug );
	Journal_DPrintf( self->debug, "In: %s()\n", __func__ );

	/* build all the systems of linear equations */
	for ( sle_I = 0; sle_I < self->slEquations->count; sle_I++ ) {
		Build( self->slEquations->data[sle_I], self, False );
	}

	/* TODO:
	This call shouldn't really be necessary - each variable used should be built
	by the FeVariable that needs it, via its dofLayout. However, unfortunately with
	"Vector" variables that use it, the app fails without this call - since otherwise
	the "velocity" variable doesn't get build, only "vx", "vy" and "vz".
	Need to debug this properly later.
	*/
	Variable_Register_BuildAll( self->variable_Register );

	Stream_UnIndentBranch( StgFEM_Debug );
}


void _FiniteElementContext_Initialise( void* context ) {
	FiniteElementContext* self = (FiniteElementContext*)context;
	SystemLinearEquations_Index sle_I;
	
	Stream_IndentBranch( StgFEM_Debug );
	Journal_DPrintf( self->debug, "In: %s()\n", __func__ );

	/* initialise all the systems of linear equations */
	for ( sle_I = 0; sle_I < self->slEquations->count; sle_I++ ) {
		Initialise( self->slEquations->data[sle_I], self, False );
	}

	Stream_UnIndentBranch( StgFEM_Debug );
}


void _FiniteElementContext_Solve( void* context ) {
	FiniteElementContext* self = (FiniteElementContext*)context;
	SystemLinearEquations_Index sle_I;
	
	Journal_DPrintf( self->debug, "In: %s()\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );

	/* solve all the systems of linear equations */
	for ( sle_I = 0; sle_I < self->slEquations->count; sle_I++ ) {
		SystemLinearEquations*	currentSLE = (SystemLinearEquations*)self->slEquations->data[sle_I];
		Journal_DPrintf( self->debug, "Solving for this timestep the %s SLE:\n", self->slEquations->data[sle_I]->name );
		/* TODO: FeVariable should have the option of rebuilding ID and LM, based on sim.
		loop if geometry or BCs change...need to improve interface. */
		/* We set the "force" flag to True here - want the SLE to be re-solved every timestep */
		Execute( currentSLE, self, True );
	}
	
	Stream_UnIndentBranch( StgFEM_Debug );
}

void _FiniteElementContext_PostSolve( void* context ) {
	FiniteElementContext* self = (FiniteElementContext*)context;

	Journal_DPrintf( self->debug, "In: %s()\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );

	FiniteElementContext_CalcNewDt( self ) ;

	Stream_UnIndentBranch( StgFEM_Debug );
}

double _FiniteElementContext_GetDt( void* context ) {
	FiniteElementContext* self = (FiniteElementContext*)context;

	return self->dt;
}

double FiniteElementContext_CalcNewDt( void* context ) {
	FiniteElementContext* self = (FiniteElementContext*)context;

	self->prevTimestepDt = self->dt;
	
	if ( self->calcDtEP->hooks->count == 0 ) {
		self->dt = 0.0;
	}
	else {
		self->dt = _EntryPoint_Run_Class_Minimum_VoidPtr( self->calcDtEP, self );
	}	
		
	if ( ( self->timeStep > 1 ) && ( self->limitTimeStepIncreaseRate == True ) ) {
		double  maxAllowedDt = self->prevTimestepDt * ( 1 + self->maxTimeStepIncreasePercentage / 100 );

		if ( self->dt > maxAllowedDt ) {
			int prevContextPrintingRank = Stream_GetPrintingRank( self->info );
			/* We assume the dt calculation will be the same across all procs since its a global
			  operation, so only print this once */
			Stream_SetPrintingRank( self->info, 0 );
			Journal_Printf( 
				self->info, 
				"In %s(): dt calculated was %g (time), but prev timestep's dt\n"
				"was %g (time) and max allowed increase percentage is %.2f\n, thus limiting current\n"
				"dt to %g (time).\n", 
				__func__, 
				self->dt, 
				self->prevTimestepDt,
				self->maxTimeStepIncreasePercentage, 
				maxAllowedDt );
			self->dt = maxAllowedDt;
			Stream_SetPrintingRank( self->info, prevContextPrintingRank );
		}
	}
	
	return self->dt;
}


void _FiniteElementContext_SaveFeVariables( void* context ) {
	FiniteElementContext*     self = (FiniteElementContext*) context;
	Index                     var_I = 0;
	FieldVariable*            fieldVar = NULL;
	FeVariable*               feVar = NULL;
	char*                     outputPathString = NULL;
	Index                     outputStrLen = 0;
	Dictionary*               dictionary = self->dictionary;
	Dictionary_Entry_Value*   feVarsList = NULL;

	outputStrLen = strlen(self->outputPath) + 1 + 1;
	if ( strlen(self->checkPointPrefixString) > 0 ) {
		outputStrLen += strlen(self->checkPointPrefixString) + 1;
	}
	outputPathString = Memory_Alloc_Array( char, outputStrLen, "outputPathString" );

	if ( strlen(self->checkPointPrefixString) > 0 ) {
		sprintf( outputPathString, "%s/%s.", self->outputPath, self->checkPointPrefixString );
	}
	else {
		sprintf( outputPathString, "%s/", self->outputPath );
	}

	/* Save the variables that have had their "isCheckpointedAndReloaded" flag enabled - 
	 *  default is true, but the user may restrict the list by specifying the "FieldVariablesToCheckpoint"
	 *  flag in their constructor - see _FeVariable_Construct().
	 */ 	
	for ( var_I = 0; var_I < self->fieldVariable_Register->objects->count; var_I++ ) {
		fieldVar = FieldVariable_Register_GetByIndex( self->fieldVariable_Register, var_I );

		if ( Stg_Class_IsInstance( fieldVar, FeVariable_Type ) ) {
			feVar = (FeVariable*)fieldVar;

			if ( feVar->isCheckpointedAndReloaded ) {
				FeVariable_SaveToFile( feVar, outputPathString, self->timeStep );
			}
		}
	}

	Memory_Free( outputPathString );
}


void _FiniteElementContext_SaveSwarms( void* context ) {
	
	Swarm_Register_SaveAllRegisteredSwarms( 
		Swarm_Register_GetSwarm_Register(), context );
}
