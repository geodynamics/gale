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
** $Id: AdvectionDiffusionSLE.c 999 2008-01-09 04:13:42Z DavidLee $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <assert.h>
#include <string.h>

#include "mpi.h"
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"

#include "types.h"
#include "AdvectionDiffusionSLE.h"
#include "UpwindParameter.h"
#include "Residual.h"
#include "Multicorrector.h"
#include "Timestep.h"

const Type AdvectionDiffusionSLE_Type = "AdvectionDiffusionSLE";

AdvectionDiffusionSLE* AdvectionDiffusionSLE_New( 
		Name                                               name,
		SLE_Solver*                                        solver,
		FiniteElementContext*                              context,
		Bool                                               isNonLinear,
		double                                             nonLinearTolerance,
		Iteration_Index                                    nonLinearMaxIterations,
		Bool                                               killNonConvergent,		
		EntryPoint_Register*                               entryPoint_Register,
		MPI_Comm                                           comm,
		FeVariable*                                        phiField,
		ForceVector*                                       residual,
		Stg_Component*                                     massMatrix,
		Dimension_Index                                    dim,
		double                                             courantFactor,
		Variable_Register*                                 variable_Register,
		FieldVariable_Register*                            fieldVariable_Register ) 
{	
	AdvectionDiffusionSLE* self = (AdvectionDiffusionSLE*)	 _AdvectionDiffusionSLE_DefaultNew( name );
	
	AdvectionDiffusionSLE_InitAll( 		
			self,
			solver,
			context,
			isNonLinear,
			nonLinearTolerance, 
			nonLinearMaxIterations,
			killNonConvergent, 			
			entryPoint_Register,
			comm,
			phiField,
			residual,
			massMatrix,
			dim,
			courantFactor,
			variable_Register,
			fieldVariable_Register );

	 return self;
}

/* Creation implementation / Virtual constructor */
AdvectionDiffusionSLE* _AdvectionDiffusionSLE_New( 
		SizeT                                              sizeOfSelf,
		Type                                               type,
		Stg_Class_DeleteFunction*                          _delete,
		Stg_Class_PrintFunction*                           _print,
		Stg_Class_CopyFunction*                            _copy, 
		Stg_Component_DefaultConstructorFunction*          _defaultConstructor,
		Stg_Component_ConstructFunction*                   _construct,
		Stg_Component_BuildFunction*                       _build,
		Stg_Component_InitialiseFunction*                  _initialise,
		Stg_Component_ExecuteFunction*                     _execute,
		Stg_Component_DestroyFunction*                     _destroy,
		SystemLinearEquations_LM_SetupFunction*            _LM_Setup,
		SystemLinearEquations_MatrixSetupFunction*         _matrixSetup,
		SystemLinearEquations_VectorSetupFunction*         _vectorSetup,		
		SystemLinearEquations_UpdateSolutionOntoNodesFunc* _updateOntoNodes,
		SystemLinearEquations_MG_SelectStiffMatsFunc*      _mgSelectStiffMats, 
		Name                                               name )
{
	AdvectionDiffusionSLE* self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(AdvectionDiffusionSLE) );
	self = (AdvectionDiffusionSLE*) _SystemLinearEquations_New( 
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
		_LM_Setup,
		_matrixSetup,
		_vectorSetup, 
		_updateOntoNodes,
		_mgSelectStiffMats, 
		name );

	return self;
}

void _AdvectionDiffusionSLE_Init(
		void*                                              sle,
		FeVariable*                                        phiField,
		ForceVector*                                   residual,
		Stg_Component*                                     massMatrix,
		Dimension_Index                                    dim,
		double                                             courantFactor,
		Variable_Register*                                 variable_Register,  
		FieldVariable_Register*                            fieldVariable_Register )		
{
	AdvectionDiffusionSLE* self          = (AdvectionDiffusionSLE*)sle;

	self->isConstructed = True;

	/* Assign values */
	self->phiField                 = phiField;
	self->residual                 = residual;
	self->massMatrix               = massMatrix;
	self->dim                      = dim;
	self->courantFactor            = courantFactor;

	/* Solution Vectors are loaded up as part of the algorithm so we can remove this one */
	EP_Remove( self->executeEP, "UpdateSolutionOntoNodes" );
	EP_Remove( self->executeEP, "MatrixSetup" );
	EP_Remove( self->executeEP, "VectorSetup" );

	/* Put Pointer of Solver onto vectors */
	if (residual) {
		residual->applicationDepExtraInfo   = (Stg_Component*) self;
		SystemLinearEquations_AddForceVector( self, residual );
	}
	if (massMatrix && Stg_Class_IsInstance( massMatrix, ForceVector_Type ) ) {
		((ForceVector*) massMatrix)->applicationDepExtraInfo = (Stg_Component*) self;
		SystemLinearEquations_AddForceVector( self, massMatrix );
	}
	else if (massMatrix && Stg_Class_IsInstance( massMatrix, StiffnessMatrix_Type ) ) {
		((StiffnessMatrix*) massMatrix)->applicationDepInfo = (Stg_Component*) self;
		SystemLinearEquations_AddStiffnessMatrix( self, massMatrix );
	}

	self->variableReg = variable_Register;
	self->fieldVariableReg = fieldVariable_Register;

	if ( self->context ) {
    /* Create a specific name for the calcDt hook */
    char* tmpName = Memory_Alloc_Array_Unnamed( char, strlen(self->name) + 7 + 1 );
    sprintf( tmpName, "%s_CalcDt", self->name );
		EntryPoint_AppendClassHook( self->context->calcDtEP, tmpName, AdvectionDiffusionSLE_CalculateDt, self->type, self );
		//EP_AppendClassHook( self->context->calcDtEP, AdvectionDiffusionSLE_CalculateDt, self );
    Memory_Free( tmpName );
  }
}	

void AdvectionDiffusionSLE_InitAll(
		void*                                              sle,
		SLE_Solver*                                        solver,
		FiniteElementContext*                              context,
		Bool                                               isNonLinear,
		double                                             nonLinearTolerance,
		Iteration_Index                                    nonLinearMaxIterations,
		Bool                                               killNonConvergent,		
		EntryPoint_Register*                               entryPoint_Register,
		MPI_Comm                                           comm,		
		FeVariable*                                        phiField,
		ForceVector*                                   residual,
		Stg_Component*                                     massMatrix,
		Dimension_Index                                    dim,
		double                                             courantFactor,
		Variable_Register*                                 variable_Register,  
		FieldVariable_Register*                            fieldVariable_Register )
{
	AdvectionDiffusionSLE* self          = (AdvectionDiffusionSLE*)sle;

	SystemLinearEquations_InitAll( 
			self,
			solver,
			NULL, /* not sure if a non linear solver is req'd. for the advection diffusion SLE */
			context,
			isNonLinear,
			nonLinearTolerance, 
			nonLinearMaxIterations,
			killNonConvergent,
			entryPoint_Register, 
			comm );
	_AdvectionDiffusionSLE_Init( 		
			self,
			phiField,
			residual,
			massMatrix,
			dim,
			courantFactor,
			variable_Register,
			fieldVariable_Register );
}

/** Virtual Functions from "Class" Class */
void _AdvectionDiffusionSLE_Delete( void* sle ) {
	AdvectionDiffusionSLE* self = (AdvectionDiffusionSLE*)sle;

	Stg_Class_Delete( self->massMatrix );
	Stg_Class_Delete( self->residual );

	Memory_Free( self->phiDotArray );
	_SystemLinearEquations_Delete( self );
}

void _AdvectionDiffusionSLE_Print( void* sle, Stream* stream ) {
	AdvectionDiffusionSLE* self = (AdvectionDiffusionSLE*) sle;
	
	Journal_Printf( stream, "Printing contents of SLE Advection Diffusion '%s':\n", self->name );
	Stream_Indent( stream );
		Journal_PrintPointer( stream, self );
		_Stg_Component_Print( self, stream );

		Journal_Printf( stream, "Items from Constructor:\n");
		Stream_Indent( stream );
			Journal_PrintPointer( stream, self->phiField );
			Journal_PrintPointer( stream, self->residual );
			Journal_PrintPointer( stream, self->massMatrix );
			Journal_PrintUnsignedInt( stream, self->dim );
		Stream_UnIndent( stream );

		Journal_Printf( stream, "Items from created by SLE self:\n");
		Stream_Indent( stream );
			Journal_PrintPointer( stream, self->phiDotArray );
			Journal_PrintPointer( stream, self->phiDotField );
			Journal_PrintPointer( stream, self->phiDotDofLayout );
			Journal_PrintPointer( stream, self->phiVector );
			Journal_PrintPointer( stream, self->phiDotVector );
		Stream_UnIndent( stream );

		Journal_Printf( stream, "Parameters set from dictionary:\n");
		Stream_Indent( stream );
			Journal_PrintDouble( stream, self->courantFactor );
		Stream_UnIndent( stream );
					
		Journal_Printf( stream, "Stored Values:\n");
		Stream_Indent( stream );
			Journal_PrintDouble( stream, self->maxDiffusivity );
		Stream_UnIndent( stream );

	Stream_UnIndent( stream );
}


void* _AdvectionDiffusionSLE_Copy( void* _sle, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	AdvectionDiffusionSLE*	self = (AdvectionDiffusionSLE*) _sle;
	AdvectionDiffusionSLE*	newSLE;
	
	newSLE = _SystemLinearEquations_Copy( self, dest, deep, nameExt, ptrMap );
	
	/* TODO: Copy Method */
	abort();
	
	return (void*)newSLE;
}

void* _AdvectionDiffusionSLE_DefaultNew( Name name ) {
	return (void*) _AdvectionDiffusionSLE_New( 
		sizeof(AdvectionDiffusionSLE), 
		AdvectionDiffusionSLE_Type, 
		_AdvectionDiffusionSLE_Delete, 
		_AdvectionDiffusionSLE_Print, 
		_AdvectionDiffusionSLE_Copy,
		_AdvectionDiffusionSLE_DefaultNew,
		_AdvectionDiffusionSLE_Construct,
		_AdvectionDiffusionSLE_Build,
		_AdvectionDiffusionSLE_Initialise,
		_AdvectionDiffusionSLE_Execute,
		_AdvectionDiffusionSLE_Destroy,
		_SystemLinearEquations_LM_Setup,
		_SystemLinearEquations_MatrixSetup,
		_SystemLinearEquations_VectorSetup,
		_SystemLinearEquations_UpdateSolutionOntoNodes, 
		_SystemLinearEquations_MG_SelectStiffMats, 
		name );
}

void _AdvectionDiffusionSLE_Construct( void* sle, Stg_ComponentFactory* cf, void* data ) {
	AdvectionDiffusionSLE*            self                   = (AdvectionDiffusionSLE*) sle;
	Stream*                           error                  = Journal_Register( Error_Type, self->type );
	FeVariable*                       phiField;
	ForceVector*                      residual;
	Stg_Component*                    massMatrix;
	Dimension_Index                   dim;
	double                            courantFactor;
	FieldVariable_Register*           fieldVariable_Register;
	Variable_Register*                variable_Register;

	/* Construct Parent */
	_SystemLinearEquations_Construct( self, cf, data );

	/* Get Registers */
	variable_Register = self->context->variable_Register; 
	assert( variable_Register );
	fieldVariable_Register = self->context->fieldVariable_Register; 
	assert( fieldVariable_Register );

	/* Get Dependency Stg_Components */
	phiField   =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "PhiField",   FeVariable,    True, data );
	residual   =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "Residual",   ForceVector,   True, data );
	massMatrix =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "MassMatrix", Stg_Component, True, data );

	dim = Stg_ComponentFactory_GetRootDictUnsignedInt( cf, "dim", 0 );

	courantFactor = Stg_ComponentFactory_GetDouble( cf, self->name, "courantFactor", 0.5 );
	Journal_Firewall( 0.0 < courantFactor && courantFactor <= 1.0, 
		error, "In func %s: CourantFactor read in from dictionary = %2.4f - This must be from 0 - 1.\n", 
		__func__, courantFactor );

	_AdvectionDiffusionSLE_Init(
			self,
			phiField,
			residual,
			massMatrix,
			dim,
			courantFactor,
			variable_Register, 
			fieldVariable_Register );
}

void _AdvectionDiffusionSLE_Destroy( void* sle, void* data ) {
	AdvectionDiffusionSLE*            self           = (AdvectionDiffusionSLE*) sle;

	_SystemLinearEquations_Destroy( self, data );
}


/** Virtual Functions from "Stg_Component" Class */
void _AdvectionDiffusionSLE_Build( void* sle, void* data ) {
	AdvectionDiffusionSLE*	self = (AdvectionDiffusionSLE*) sle;
	Stream*                 errorStream = Journal_MyStream( Error_Type, self );
	Index                   forceTerm_I;
	Index                   forceTermCount = Stg_ObjectList_Count( self->residual->forceTermList );
	ForceTerm*              forceTerm;
	unsigned int           *nodeDomainCountPtr;
	Variable*              variable;
	Node_DomainIndex       node_I;

	Journal_DPrintf( self->debug, "In %s()\n", __func__ );

	/* Create New FeVariable for Phi Dot */
	if (self->phiField) {
		Variable_Register*	variable_Register;
		FieldVariable_Register*	fieldVariable_Register;
		char* fieldName, *dofName, *fieldDotName, *phiVecName, *phiDotVecName;

		variable_Register = self->variableReg;
		fieldVariable_Register = self->fieldVariableReg;

		Stg_Component_Build( self->phiField->feMesh, NULL, False );

		assert( Class_IsSuper( self->phiField->feMesh->topo, IGraph ) );
		nodeDomainCountPtr = &((IGraph*)self->phiField->feMesh->topo)->remotes[MT_VERTEX]->nDomains;

    /* must create unique names otherwise multiple instances of this component
     * will index incorrect instances of this component's data */
    fieldName = Memory_Alloc_Array_Unnamed( char, strlen(self->name)+8 );
    sprintf( fieldName, "%s-phiDot", self->name );

    dofName = Memory_Alloc_Array_Unnamed( char, strlen(self->name)+11 );
    sprintf( dofName, "%s-dofLayout", self->name );

    fieldDotName  = Memory_Alloc_Array_Unnamed( char, strlen(self->name)+13 );
    sprintf( fieldDotName, "%s-phiDotField", self->name );

    phiVecName = Memory_Alloc_Array_Unnamed( char, strlen(self->name)+11 );
    sprintf( phiVecName, "%s-phiVector", self->name );

    phiDotVecName = Memory_Alloc_Array_Unnamed( char, strlen(self->name)+14 );
    sprintf( phiDotVecName, "%s-phiDotVector", self->name );

		variable = Variable_NewScalar( 
			fieldName, 
			Variable_DataType_Double, 
			nodeDomainCountPtr, 
			NULL,
			(void**)&self->phiDotArray, 
			variable_Register );

    self->phiDotDofLayout = DofLayout_New( dofName, variable_Register, *nodeDomainCountPtr, NULL );
		//self->phiDotDofLayout = DofLayout_New( "dofLayout1", variable_Register, *nodeDomainCountPtr, NULL );
		for( node_I = 0; node_I < *nodeDomainCountPtr ; node_I++ ) 
			DofLayout_AddDof_ByVarName( self->phiDotDofLayout, variable->name, node_I );

		self->phiDotField = FeVariable_New_FromTemplate(
			fieldDotName,
			self->phiField,
			self->phiDotDofLayout,
			NULL,
			False, False,
			fieldVariable_Register );

		/* Construct Solution Vectors */
		self->phiVector    = SolutionVector_New( phiVecName, self->phiField->communicator, self->phiField );
		self->phiDotVector = SolutionVector_New( phiDotVecName, self->phiField->communicator, self->phiDotField );

    /* free original name variables */
    Memory_Free(fieldName);
    Memory_Free(dofName);
    Memory_Free(fieldDotName);
    Memory_Free(phiVecName);
    Memory_Free(phiDotVecName);
	}

	_SystemLinearEquations_Build( self, data );

	/* Get pointer to residual force term 
	 * All Advection Diffusion SLE's need a force term of type AdvDiffResidualForceTerm for the algorithm to work
	 * this chunk of code is making sure that one and only one is registered to the force vector */
	for ( forceTerm_I = 0 ; forceTerm_I < forceTermCount ; forceTerm_I++ ) {
		forceTerm = (ForceTerm*) Stg_ObjectList_At( self->residual->forceTermList, forceTerm_I );

		if ( Stg_Class_IsInstance( forceTerm, AdvDiffResidualForceTerm_Type ) ) {
			/* Check to make sure this force term is unique */
			if (self->advDiffResidualForceTerm != NULL) {
				Journal_Firewall( self->advDiffResidualForceTerm == NULL, errorStream, 
						"Error - More than one force term of type '%s' registered on %s '%s'. \n\tThey are %s and %s.\n", 
						AdvDiffResidualForceTerm_Type, self->type, self->name, 
						self->advDiffResidualForceTerm->name, forceTerm->name);
			}

			/* Store pointer to force term */
			self->advDiffResidualForceTerm = (AdvDiffResidualForceTerm*) forceTerm;
                       
			/* HACK */
			Stg_Component_Build( self->advDiffResidualForceTerm->velocityField, data, False );

		}
	}
	/* Ensure that there is at least one force term with proper residual type */
	Journal_Firewall( self->advDiffResidualForceTerm != NULL, errorStream,
			"Error - No force terms of type '%s' registered on %s '%s'.\n",
			AdvDiffResidualForceTerm_Type, self->type, self->name );

	if ( self->phiDotField ) {
		self->phiDotArray = Memory_Alloc_Array( 
			double, Mesh_GetDomainSize( self->phiDotField->feMesh, MT_VERTEX ), "phiDotArray" );
	}

	/* Force Vectors */
	if ( self->residual )
		Stg_Component_Build( self->residual, data, False );
	if ( self->massMatrix )
		Stg_Component_Build( self->massMatrix, data, False );
	
	/* Solution Vectors */
	if ( self->phiVector )
		Stg_Component_Build( self->phiVector, data, False );
	if ( self->phiDotVector )
		Stg_Component_Build( self->phiDotVector, data, False );
}

void _AdvectionDiffusionSLE_Initialise( void* sle, void* data ) {
	AdvectionDiffusionSLE* self = (AdvectionDiffusionSLE*) sle;
	FiniteElementContext* context = (FiniteElementContext*) data;

	Journal_DPrintf( self->debug, "In %s()\n", __func__ );

	_SystemLinearEquations_Initialise( self, data );
	
	Stg_Component_Initialise( self->phiDotField, data, False );

/* 	Stream* stream = Journal_Register( Info_Type, self->type ); */
/* 	FeVariable_PrintLocalDiscreteValues( self->phiDotField, stream ); */
	
	if ( False == context->loadFromCheckPoint ) {
		DofLayout_SetAllToZero( self->phiDotField->dofLayout );
	}

	/* Force Vectors */
	Stg_Component_Initialise( self->residual, data, False );
	Stg_Component_Initialise( self->massMatrix, data, False );
	
	/* Solution Vectors */
	Stg_Component_Initialise( self->phiVector, data, False );
	Stg_Component_Initialise( self->phiDotVector, data, False );

	AdvectionDiffusionSLE_ResetStoredValues( self );
}

void _AdvectionDiffusionSLE_Execute( void* sle, void* _context ) {
	AdvectionDiffusionSLE*     self  = (AdvectionDiffusionSLE*) sle;
	FiniteElementContext*      context = (FiniteElementContext*) _context;
	double                     dt      = context->dt;
	
	AdvectionDiffusionSLE_ResetStoredValues( self );
	self->currentDt = dt;
	
	_SystemLinearEquations_Execute( self, context );
}


//Vector* _AdvectionDiffusionSLE_GetResidual( void* sle, Index fv_I ) {
Vec _AdvectionDiffusionSLE_GetResidual( void* sle, Index fv_I ) {
	AdvectionDiffusionSLE* self  = (AdvectionDiffusionSLE*) sle;
	return self->residual->vector;
}

#define SMALL_VALUE 1.0e-99

void AdvectionDiffusionSLE_ResetStoredValues( void* sle ) {
	AdvectionDiffusionSLE* self  = (AdvectionDiffusionSLE*) sle;
	
	self->maxDiffusivity = SMALL_VALUE;
}

