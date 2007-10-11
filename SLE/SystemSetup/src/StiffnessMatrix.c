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
** $Id: StiffnessMatrix.c 964 2007-10-11 08:03:06Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <mpi.h>

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/LinearAlgebra/LinearAlgebra.h"

#include "units.h"
#include "types.h"
#include "shortcuts.h"
#include "StiffnessMatrix.h"
#include "StiffnessMatrixTerm.h"
#include "SystemLinearEquations.h"
#include "PETScShellMatrix.h"
#include "EntryPoint.h"
#include "SolutionVector.h"
#include "ForceVector.h"
#include "Assembler.h"


void StiffnessMatrix_NewAssemble( void* stiffnessMatrix, Bool removeBCs, void* _sle, void* _context );
Bool StiffnessMatrix_ZeroBCsAsm_RowR( void* stiffMat, Assembler* assm );
Bool StiffnessMatrix_ZeroBCsAsm_ColR( void* stiffMat, Assembler* assm );
Bool StiffnessMatrix_BCAsm_ColR( void* stiffMat, Assembler* assm );
Bool StiffnessMatrix_TransBCAsm_ColR( void* stiffMat, Assembler* assm );
Bool StiffnessMatrix_DiagBCsAsm_RowR( void* stiffMat, Assembler* assm );


/* Textual name of this class */
const Type StiffnessMatrix_Type = "StiffnessMatrix";

/** First part of name for build entry point */
static const char	StiffnessMatrix_assembleStiffnessMatrixStr[] = "assembleStiffnessMatrix";


void* StiffnessMatrix_DefaultNew( Name name )
{
	return _StiffnessMatrix_New( 
		sizeof(StiffnessMatrix), 
		StiffnessMatrix_Type, 
		_StiffnessMatrix_Delete,
		_StiffnessMatrix_Print, 
		_StiffnessMatrix_Copy,
		StiffnessMatrix_DefaultNew,
		_StiffnessMatrix_Construct,
		_StiffnessMatrix_Build, 
		_StiffnessMatrix_Initialise,
		_StiffnessMatrix_Execute,
		_StiffnessMatrix_Destroy,
		name,
		False,
		_StiffnessMatrix_CalculateNonZeroEntries, 
		NULL, 
		NULL, 
		NULL, 
		NULL,
		0,
		False,
		False,
		NULL,
		0 );
}


StiffnessMatrix* StiffnessMatrix_New(
	Name                                             name,
	void*                                            rowVariable,
	void*                                            columnVariable,
	void*                                            rhs,
	Stg_Component*                                   applicationDepInfo,
	Dimension_Index                                  dim,
	Bool                                             isNonLinear,
	Bool                                             allowZeroElementContributions,
	void*                                            entryPoint_Register,
	MPI_Comm                                         comm )
{
	return _StiffnessMatrix_New( 
		sizeof(StiffnessMatrix), 
		StiffnessMatrix_Type, 
		_StiffnessMatrix_Delete,
		_StiffnessMatrix_Print, 
		_StiffnessMatrix_Copy,
		StiffnessMatrix_DefaultNew,
		_StiffnessMatrix_Construct,
		_StiffnessMatrix_Build, 
		_StiffnessMatrix_Initialise,
		_StiffnessMatrix_Execute,
		_StiffnessMatrix_Destroy,
		name, 
		True,
		_StiffnessMatrix_CalculateNonZeroEntries,
		rowVariable, 
		columnVariable, 
		rhs,
		applicationDepInfo,
		dim,
		isNonLinear,
		allowZeroElementContributions,
		entryPoint_Register,
		comm );
}


StiffnessMatrix* _StiffnessMatrix_New(
	SizeT                                            _sizeOfSelf,
	Type                                             type,
	Stg_Class_DeleteFunction*                        _delete,
	Stg_Class_PrintFunction*                         _print,
	Stg_Class_CopyFunction*                          _copy, 
	Stg_Component_DefaultConstructorFunction*        _defaultConstructor,
	Stg_Component_ConstructFunction*                 _construct,
	Stg_Component_BuildFunction*                     _build,
	Stg_Component_InitialiseFunction*                _initialise,
	Stg_Component_ExecuteFunction*                   _execute,
	Stg_Component_DestroyFunction*                   _destroy,
	Name                                             name,
	Bool                                             initFlag,
	StiffnessMatrix_CalculateNonZeroEntriesFunction* _calculateNonZeroEntries,
	void*                                            rowVariable,
	void*                                            columnVariable,
	void*                                            rhs,
	Stg_Component*                                   applicationDepInfo,
	Dimension_Index                                  dim,
	Bool                                             isNonLinear,
	Bool                                             allowZeroElementContributions,
	void*                                            entryPoint_Register,
	MPI_Comm                                         comm )
{
	StiffnessMatrix*	self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(StiffnessMatrix) );
	self = (StiffnessMatrix*)_Stg_Component_New(
		_sizeOfSelf,
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
		NON_GLOBAL );
	
	/* General info */
	
	/* Virtual functions */
	self->_calculateNonZeroEntries = _calculateNonZeroEntries;
	
	if( initFlag ){
		_StiffnessMatrix_Init( self, rowVariable, columnVariable, rhs, applicationDepInfo, dim,
				       isNonLinear, allowZeroElementContributions, entryPoint_Register, comm );
	}
	
	return self;
}

void _StiffnessMatrix_Init(
	StiffnessMatrix*                                 self,
	void*                                            rowVariable,
	void*                                            columnVariable,
	void*                                            rhs,
	Stg_Component*                                   applicationDepInfo,
	Dimension_Index                                  dim,
	Bool                                             isNonLinear,
	Bool                                             allowZeroElementContributions,
	void*                                            entryPoint_Register,
	MPI_Comm                                         comm )
{
	Stream*		error = Journal_Register( ErrorStream_Type, self->type );
	Stream*		stream;

	/* General and Virtual info should already be set */
	stream = Journal_Register( Info_Type, self->type );
	Stream_SetPrintingRank( stream, 0 );
	
	/* StiffnessMatrix info */
	self->isConstructed = True;
	self->debug = Stream_RegisterChild( StgFEM_SLE_SystemSetup_Debug, self->type );
	Journal_Firewall( (rowVariable != NULL), error, "Error: NULL row FeVariable provided to \"%s\" %s.\n",
			  self->name, self->type );
	self->rowVariable = (FeVariable*)rowVariable;
	Journal_Firewall( (columnVariable != NULL), error, "Error: NULL column FeVariable provided to \"%s\" %s.\n",
			  self->name, self->type );
	self->columnVariable = (FeVariable*)columnVariable;
	Journal_Firewall( (rhs != NULL), error, "Error: NULL rhs ForceVector provided to \"%s\" %s.\n",
			  self->name, self->type );
	self->rhs = (ForceVector*)rhs;
	self->applicationDepInfo = applicationDepInfo;
	self->comm = comm;
	self->dim = dim;
	self->isNonLinear = isNonLinear;
	self->allowZeroElementContributions = allowZeroElementContributions;
	
	self->rowLocalSize = 0;
	self->colLocalSize = 0;
	self->nonZeroCount = 0;
	self->diagonalNonZeroCount = 0;
	self->offDiagonalNonZeroCount = 0;
	self->diagonalNonZeroIndices = NULL;
	self->offDiagonalNonZeroIndices = NULL;
	
	self->entryPoint_Register = (EntryPoint_Register*)entryPoint_Register;

	Stg_asprintf( &self->_assembleStiffnessMatrixEPName, "%s-%s", self->name, StiffnessMatrix_assembleStiffnessMatrixStr );
	self->assembleStiffnessMatrix = FeEntryPoint_New( self->_assembleStiffnessMatrixEPName,
							  FeEntryPoint_AssembleStiffnessMatrix_CastType );
	EntryPoint_Register_Add( self->entryPoint_Register, self->assembleStiffnessMatrix );

	self->stiffnessMatrixTermList = Stg_ObjectList_New();

	/* Set default function for Global Stiffness Matrix Assembly */
	/*EP_ReplaceAll( self->assembleStiffnessMatrix, StiffnessMatrix_GlobalAssembly_General );*/
	EP_ReplaceAll( self->assembleStiffnessMatrix, StiffnessMatrix_NewAssemble );

	/* We need some assembler contexts. */
	self->zeroBCsAsm = Assembler_New();
	self->bcAsm = Assembler_New();
	self->transBCAsm = Assembler_New();
	if( rowVariable == columnVariable )
		self->diagBCsAsm = Assembler_New();

	self->elStiffMat = NULL;
	self->bcVals = NULL;
	self->nRowDofs = 0;
	self->nColDofs = 0;
	self->transRHS = NULL;

	self->rowInc = IArray_New();
	self->colInc = IArray_New();
}


void _StiffnessMatrix_Delete( void* stiffnessMatrix ) {
	StiffnessMatrix* self = (StiffnessMatrix*)stiffnessMatrix;
	
	Journal_DPrintf( self->debug, "In %s - for matrix %s\n", __func__, self->name );
	FreeObject( self->stiffnessMatrixTermList );
	FreeArray( self->_assembleStiffnessMatrixEPName );
	FreeArray( self->diagonalNonZeroIndices );
	FreeArray( self->offDiagonalNonZeroIndices );
	FreeObject( self->zeroBCsAsm );
	FreeObject( self->bcAsm );
	FreeObject( self->transBCAsm );
	FreeObject( self->diagBCsAsm );
	/* Don't delete entry points: E.P. register will delete them automatically */

	NewClass_Delete( self->rowInc );
	NewClass_Delete( self->colInc );

	/* Stg_Class_Delete parent*/
	_Stg_Component_Delete( self );
	
}

void _StiffnessMatrix_Print( void* stiffnessMatrix, Stream* stream ) {
	StiffnessMatrix* self = (StiffnessMatrix*)stiffnessMatrix;
	
	/* Set the Journal for printing informations */
	Stream* stiffnessMatrixStream = stream;
	
	/* General info */
	Journal_Printf( stiffnessMatrixStream, "StiffnessMatrix (ptr): %p\n", self );
	
	/* Print parent */
	_Stg_Component_Print( self, stiffnessMatrixStream );
	
	/* Virtual info */
	Journal_Printf( stiffnessMatrixStream, "\t_build (func ptr): %p\n", self->_build );
	
	Journal_Printf( stiffnessMatrixStream, "\tassembleStiffnessMatrix e.p. (ptr): %p\n", self->assembleStiffnessMatrix );
	EntryPoint_PrintConcise( self->assembleStiffnessMatrix, stream );
	
	/* StiffnessMatrix info */
	Journal_Printf( stiffnessMatrixStream, "\trowVariable (ptr): %p\n", self->rowVariable );
	Journal_Printf( stiffnessMatrixStream, "\t\tvariable name: %s\n", self->rowVariable->name );
	Journal_Printf( stiffnessMatrixStream, "\tcolumnVariable (ptr): %p\n", self->columnVariable );
	Journal_Printf( stiffnessMatrixStream, "\t\tvariable name: %s\n", self->columnVariable->name );
	Journal_Printf( stiffnessMatrixStream, "\tMatrix (ptr): %p\n", self->matrix );
	Journal_Printf( stiffnessMatrixStream, "\tComm: %u\n", self->comm );
	Journal_Printf( stiffnessMatrixStream, "\trowLocalSize: %u\n", self->rowLocalSize );
	Journal_Printf( stiffnessMatrixStream, "\tcolLocalSize: %u\n", self->colLocalSize );
	Journal_Printf( stiffnessMatrixStream, "\tnonZeroCount: %u\n", self->nonZeroCount );
	Journal_Printf( stiffnessMatrixStream, "\tisNonLinear: %s\n", StG_BoolToStringMap[self->isNonLinear] );
	Journal_Printf( stiffnessMatrixStream, "\tallowZeroElementContributions: %s\n",
			StG_BoolToStringMap[self->allowZeroElementContributions] );
}


void* _StiffnessMatrix_Copy( void* stiffnessMatrix, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	StiffnessMatrix*	self = (StiffnessMatrix*)stiffnessMatrix;
	StiffnessMatrix*	newStiffnessMatrix;
	PtrMap*			map = ptrMap;
	Bool			ownMap = False;
	
	if( !map ) {
		map = PtrMap_New( 10 );
		ownMap = True;
	}
	
	newStiffnessMatrix = _Stg_Component_Copy( self, dest, deep, nameExt, map );
	
	/* Virtual functions */
	newStiffnessMatrix->_calculateNonZeroEntries = self->_calculateNonZeroEntries;
	
	/* TODO: copy matrix */
	newStiffnessMatrix->matrix = self->matrix;
	newStiffnessMatrix->entryPoint_Register = self->entryPoint_Register;
	newStiffnessMatrix->comm = self->comm;
	newStiffnessMatrix->rowLocalSize = self->rowLocalSize;
	newStiffnessMatrix->colLocalSize = self->colLocalSize;
	newStiffnessMatrix->dim = self->dim;
	newStiffnessMatrix->nonZeroCount = self->nonZeroCount;
	newStiffnessMatrix->diagonalNonZeroCount = self->diagonalNonZeroCount;
	newStiffnessMatrix->offDiagonalNonZeroCount = self->offDiagonalNonZeroCount;
	
	if( deep ) {
		newStiffnessMatrix->debug = (Stream*)Stg_Class_Copy( self->debug, NULL, deep, nameExt, map );
		newStiffnessMatrix->rowVariable = (FeVariable*)Stg_Class_Copy( self->rowVariable, NULL, deep, nameExt, map );
		newStiffnessMatrix->columnVariable = (FeVariable*)Stg_Class_Copy( self->columnVariable, NULL, deep, nameExt, map );
		newStiffnessMatrix->rhs =(ForceVector*)Stg_Class_Copy( self->rhs, NULL, deep, nameExt, map );
		newStiffnessMatrix->assembleStiffnessMatrix = (FeEntryPoint*)Stg_Class_Copy( self->assembleStiffnessMatrix, NULL, deep, nameExt, map );
		
		if( self->_assembleStiffnessMatrixEPName ) {
			if( nameExt ) {
				Stg_asprintf( &newStiffnessMatrix->_assembleStiffnessMatrixEPName, "%s%s", 
					      self->_assembleStiffnessMatrixEPName, nameExt );
			}
			else {
				newStiffnessMatrix->_assembleStiffnessMatrixEPName = StG_Strdup( self->_assembleStiffnessMatrixEPName );
			}
		}
		else {
			newStiffnessMatrix->_assembleStiffnessMatrixEPName = NULL;
		}
		
		/* Arrays */
		if( (newStiffnessMatrix->diagonalNonZeroIndices = PtrMap_Find( map, self->diagonalNonZeroIndices )) == NULL ) {
			if( self->diagonalNonZeroIndices ) {
				newStiffnessMatrix->diagonalNonZeroIndices = Memory_Alloc_Array( Index, 
												 newStiffnessMatrix->rowLocalSize, "diagonalNonZeroIndices" );
				memcpy( newStiffnessMatrix->diagonalNonZeroIndices, self->diagonalNonZeroIndices, 
					newStiffnessMatrix->rowLocalSize * sizeof( Index ) );
				PtrMap_Append( map, self->diagonalNonZeroIndices, newStiffnessMatrix->diagonalNonZeroIndices );
			}
			else {
				newStiffnessMatrix->diagonalNonZeroIndices = NULL;
			}
		}
		
		if( (newStiffnessMatrix->offDiagonalNonZeroIndices = PtrMap_Find( map, self->offDiagonalNonZeroIndices )) == NULL ) {
			if( self->offDiagonalNonZeroIndices ) {
				newStiffnessMatrix->offDiagonalNonZeroIndices = Memory_Alloc_Array( Index, 
												    newStiffnessMatrix->rowLocalSize, "diagonalNonZeroIndices" );
				memcpy( newStiffnessMatrix->offDiagonalNonZeroIndices, self->offDiagonalNonZeroIndices, 
					newStiffnessMatrix->rowLocalSize * sizeof( Index ) );
				PtrMap_Append( map, self->offDiagonalNonZeroIndices, newStiffnessMatrix->offDiagonalNonZeroIndices );
			}
			else {
				newStiffnessMatrix->offDiagonalNonZeroIndices = NULL;
			}
		}
	}
	else {
		newStiffnessMatrix->debug = self->debug;
		newStiffnessMatrix->rowVariable = self->rowVariable;
		newStiffnessMatrix->columnVariable = self->columnVariable;
		newStiffnessMatrix->rhs = self->rhs;
		newStiffnessMatrix->diagonalNonZeroIndices = self->diagonalNonZeroIndices;
		newStiffnessMatrix->offDiagonalNonZeroIndices = self->offDiagonalNonZeroIndices;
	}
	
	if( ownMap ) {
		Stg_Class_Delete( map );
	}
	
	return (void*)newStiffnessMatrix;
}

void _StiffnessMatrix_Construct( void* stiffnessMatrix, Stg_ComponentFactory* cf, void* data ) {
	StiffnessMatrix* self               = (StiffnessMatrix*)stiffnessMatrix;
	Stream*		 stream;
	FeVariable*      rowVar             = NULL;
	FeVariable*      colVar             = NULL;
	ForceVector*     fVector            = NULL;
	Stg_Component*   applicationDepInfo = NULL;
	void*            entryPointRegister = NULL;
	Dimension_Index  dim                = 0;
	Bool             isNonLinear;
	Bool             allowZeroElementContributions;
	
	rowVar             = Stg_ComponentFactory_ConstructByKey( cf, self->name, "RowVariable",        FeVariable,    True, data );
	colVar             = Stg_ComponentFactory_ConstructByKey( cf, self->name, "ColumnVariable",     FeVariable,    True, data );
	fVector            = Stg_ComponentFactory_ConstructByKey( cf, self->name, "RHS",                ForceVector,   False, data );
	applicationDepInfo = Stg_ComponentFactory_ConstructByKey( cf, self->name, "ApplicationDepInfo", Stg_Component, False, data);
	
	entryPointRegister = Stg_ObjectList_Get( cf->registerRegister, "EntryPoint_Register" );
	assert( entryPointRegister );
	
	dim = Stg_ComponentFactory_GetRootDictUnsignedInt( cf, "dim", 0 );
	assert( dim );

	isNonLinear = Stg_ComponentFactory_GetBool( cf, self->name, "isNonLinear", False );

	/* Default is to allow zero element contributions - to allow backward compatibility */
	allowZeroElementContributions = Stg_ComponentFactory_GetBool( cf, self->name, "allowZeroElementContributions", True );

	_StiffnessMatrix_Init( 
		self, 
		rowVar, 
		colVar, 
		fVector, 
		applicationDepInfo, 
		dim,
		isNonLinear,
		allowZeroElementContributions,
		entryPointRegister, 
		0 );

	/* Do we need to use the transpose? */
	self->transRHS = Stg_ComponentFactory_ConstructByKey( cf, self->name, "transposeRHS", ForceVector, False, data );

	/* Read the matrix type from the dictionary. */
	self->matrix = Stg_ComponentFactory_ConstructByKey( cf, self->name, "matrix", Matrix, False, data );
	if( !self->matrix ) {
#ifdef HAVE_PETSC
		self->matrix = (Matrix*)PETScMatrix_New( "" );
#else
#error *** No linear algebra package found.
#endif
	}
	else if( !strcmp( self->matrix->type, PETScShellMatrix_Type ) ) {
		EP_ReplaceAll( self->assembleStiffnessMatrix, StiffnessMatrix_ShellAssembly );
	}

	/* Setup the stream. */
	stream = Journal_Register( Info_Type, self->type );
	if( Dictionary_GetBool_WithDefault( cf->rootDict, "watchAll", False ) == True )
		Stream_SetPrintingRank( stream, STREAM_ALL_RANKS );
	else {
		unsigned	rankToWatch;

		rankToWatch = Dictionary_GetUnsignedInt_WithDefault( cf->rootDict, "rankToWatch", 0 );
		Stream_SetPrintingRank( stream, rankToWatch );
	}
}

void _StiffnessMatrix_Build( void* stiffnessMatrix, void* data ) {
	StiffnessMatrix* self = (StiffnessMatrix*)stiffnessMatrix;

	Journal_DPrintf( self->debug, "In %s - for matrix %s\n", __func__, self->name );
	Stream_IndentBranch( StgFEM_Debug );

	/* ensure variables are built */
	if( self->rowVariable )
		Stg_Component_Build( self->rowVariable, data, False );

	/* If we don't have a communicator, grab one off the mesh. */
	if( !self->comm ) {
		self->comm = Mesh_GetCommTopology( self->rowVariable->feMesh, MT_VERTEX )->mpiComm;
		Journal_Firewall( (self->comm != 0), self->debug, "Error: NULL Comm provided to \"%s\" %s.\n",
				  self->name, self->type );
	}
	
	if( self->columnVariable )
		Stg_Component_Build( self->columnVariable, data, False );
	
	
#if DEBUG
	if ( Stream_IsPrintableLevel( self->debug, 3 ) ) {
		Journal_DPrintf( self->debug, "Row variable(%s) I.D. array calc. as:\n", self->rowVariable->name );
		FeEquationNumber_PrintDestinationArray( self->rowVariable->eqNum, self->debug );
		Journal_DPrintf( self->debug, "Column variable(%s) I.D. array calc. as:\n", self->columnVariable->name );
		FeEquationNumber_PrintDestinationArray( self->columnVariable->eqNum, self->debug );
	}
#endif
	
	/* update the row and column sizes for the variables */	
	self->rowLocalSize = self->rowVariable->eqNum->localEqNumsOwnedCount;
	assert( self->rowLocalSize );
	self->colLocalSize = self->columnVariable->eqNum->localEqNumsOwnedCount;
	assert( self->colLocalSize );
	
	/* update the number of non zero entries from the finite element variables */
	StiffnessMatrix_CalculateNonZeroEntries( self );
	
	Journal_DPrintf( self->debug, "row(%s) localSize = %d : col(%s) localSize = %d \n", self->rowVariable->name,
			 self->rowLocalSize, self->columnVariable->name, self->colLocalSize );
	
/* 	assert( self->nonZeroCount ); */

	StiffnessMatrix_RefreshMatrix( self );

	Journal_DPrintf( self->debug, "Matrix allocated.\n" );

	Assembler_SetVariables( self->zeroBCsAsm, self->rowVariable, self->columnVariable );
	Assembler_SetCallbacks( self->zeroBCsAsm, 
				NULL, 
				StiffnessMatrix_ZeroBCsAsm_RowR, NULL, 
				StiffnessMatrix_ZeroBCsAsm_ColR, NULL, 
				self );
	Assembler_SetVariables( self->bcAsm, self->rowVariable, self->columnVariable );
	Assembler_SetCallbacks( self->bcAsm, 
				NULL, 
				NULL, NULL, 
				StiffnessMatrix_BCAsm_ColR, NULL, 
				self );
	Assembler_SetVariables( self->transBCAsm, self->columnVariable, self->rowVariable );
	Assembler_SetCallbacks( self->transBCAsm, 
				NULL, 
				NULL, NULL, 
				StiffnessMatrix_TransBCAsm_ColR, NULL, 
				self );
	if( self->rowVariable == self->columnVariable ) {
		Assembler_SetVariables( self->diagBCsAsm, self->rowVariable, self->columnVariable );
		Assembler_SetCallbacks( self->diagBCsAsm, 
					NULL, 
					StiffnessMatrix_DiagBCsAsm_RowR, NULL, 
					NULL, NULL, 
					self );
	}

	Stream_UnIndentBranch( StgFEM_Debug );
}


void _StiffnessMatrix_Initialise( void* stiffnessMatrix, void* data ) {
	StiffnessMatrix* self = (StiffnessMatrix*)stiffnessMatrix;
	
	Journal_DPrintf( self->debug, "In %s - for matrix %s\n", __func__, self->name );
	/* ensure variables are initialised */
	if( self->rowVariable )
		Stg_Component_Initialise( self->rowVariable, data, False );
	
	if( self->columnVariable )
		Stg_Component_Initialise( self->columnVariable, data, False );
}


void _StiffnessMatrix_Execute( void* stiffnessMatrix, void* data ) {
}

void _StiffnessMatrix_Destroy( void* stiffnessMatrix, void* data ) {
}

void StiffnessMatrix_CalculateNonZeroEntries( void* stiffnessMatrix ) {
	StiffnessMatrix* self = (StiffnessMatrix*)stiffnessMatrix;
	
	self->_calculateNonZeroEntries( self );
}

void _StiffnessMatrix_CalculateNonZeroEntries( void* stiffnessMatrix ) {
	StiffnessMatrix*	self = (StiffnessMatrix*)stiffnessMatrix;
	FeMesh*			rFeMesh = self->rowVariable->feMesh;
	FeMesh*			cFeMesh = self->columnVariable->feMesh;
	FeEquationNumber*	rowEqNum = self->rowVariable->eqNum;
	Dof_EquationNumber	currMatrixRow = 0;
	Node_LocalIndex		rowNode_lI = 0;
	Dof_EquationNumber	lowestActiveEqNumAtCurrRowNode = 0;
	Index			activeEqsAtCurrRowNodeCount = 0;
#ifdef DEBUG
	unsigned int            totalNonZeroEntries = 0;
#endif

	Journal_DPrintf( self->debug, "In %s - for matrix %s\n", __func__, self->name );
	Stream_IndentBranch( StgFEM_Debug );

	assert ( self->rowVariable );
	assert ( self->columnVariable );

	Journal_DPrintfL( self->debug, 1, "row nodeLocalCount %d\n", FeMesh_GetNodeLocalSize( rFeMesh ) );
	Journal_DPrintfL( self->debug, 1, "column nodeLocalCount %d\n", FeMesh_GetNodeLocalSize( cFeMesh ) );
	
	self->diagonalNonZeroIndices = Memory_Alloc_Array  ( Index, self->rowLocalSize, "diagonalNonZeroIndices" );
	self->offDiagonalNonZeroIndices = Memory_Alloc_Array  ( Index, self->rowLocalSize, "offDiagonalNonZeroIndices" );
	for ( rowNode_lI = 0; rowNode_lI < self->rowLocalSize; rowNode_lI++ ) { 
		self->diagonalNonZeroIndices[rowNode_lI] = 0;
		self->offDiagonalNonZeroIndices[rowNode_lI] = 0;
	}

	for( rowNode_lI = 0; rowNode_lI < FeMesh_GetNodeLocalSize( rFeMesh ); rowNode_lI++ ) {
		activeEqsAtCurrRowNodeCount = FeEquationNumber_CalculateActiveEqCountAtNode( rowEqNum, rowNode_lI,
											     &lowestActiveEqNumAtCurrRowNode );

		if ( activeEqsAtCurrRowNodeCount > 0 ) {
			currMatrixRow = lowestActiveEqNumAtCurrRowNode - rowEqNum->firstOwnedEqNum;
			
			if ( currMatrixRow >= rowEqNum->localEqNumsOwnedCount ) {
				Journal_DPrintfL( self->debug, 1, "row node %d - eq nums not owned by this proc -> continue\n",
						  rowNode_lI );
				continue;
			}
			
			_StiffnessMatrix_CalcAndUpdateNonZeroEntriesAtRowNode( self, rowNode_lI, currMatrixRow,
									       activeEqsAtCurrRowNodeCount );
		}
	}

#ifdef DEBUG
	for ( rowNode_lI = 0; rowNode_lI < self->rowLocalSize; rowNode_lI++ ) { 
		totalNonZeroEntries += self->diagonalNonZeroIndices[rowNode_lI];
		totalNonZeroEntries += self->offDiagonalNonZeroIndices[rowNode_lI];
	}	

	Journal_DPrintfL( self->debug, 1, "Calculated %d non-zero entries in Matrix (results in %d bytes storage)\n",
			  totalNonZeroEntries, totalNonZeroEntries * sizeof(double) );
#endif	

	Stream_UnIndentBranch( StgFEM_Debug );
}


void _StiffnessMatrix_CalcAndUpdateNonZeroEntriesAtRowNode(
	StiffnessMatrix*	self,
	Node_LocalIndex		rowNode_lI,
	Dof_EquationNumber	currMatrixRow,
	Index			activeEqsAtCurrRowNode )
{
	FeMesh*			rFeMesh = self->rowVariable->feMesh;
	FeMesh*			cFeMesh = self->columnVariable->feMesh;
	DofLayout*		colDofLayout = self->columnVariable->dofLayout;
	FeEquationNumber*	rowEqNum = self->rowVariable->eqNum;
	FeEquationNumber*	colEqNum = self->columnVariable->eqNum;
	Element_Index		rowNodeElement_I = 0;
	Element_DomainIndex	element_dI = 0;
	Node_DomainIndex	colNode_dI = 0;
	Dof_Index		colNodeDof_I = 0;
	Index*			countTableToAdjust = 0;
	Dof_EquationNumber	currColEqNum = 0;
	Node_DomainIndex*	uniqueRelatedColNodes = NULL;
	Node_Index		uniqueRelatedColNodesCount = 0;
	Node_Index		uniqueRelatedColNodes_AllocCount = 0;
	Node_Index		uniqueRelatedColNode_I = 0;
	Dof_Index		currNodeDof_I = 0;
	Dof_EquationNumber	currDofMatrixRow = 0;
	unsigned		nNodeInc, *nodeInc;

	Journal_DPrintfL( self->debug, 3, "In %s - for row local node %d\n", __func__, rowNode_lI );
	Stream_Indent( self->debug );

	FeMesh_GetNodeElements( rFeMesh, rowNode_lI, self->rowInc );
	nNodeInc = IArray_GetSize( self->rowInc );
	nodeInc = IArray_GetPtr( self->rowInc );
	for ( rowNodeElement_I = 0; rowNodeElement_I < nNodeInc; rowNodeElement_I++ ) {
		unsigned	nElInc;

		/* note the dI (domain index) - some of these elements may be in shadow space */
		element_dI = nodeInc[rowNodeElement_I];

		nElInc = FeMesh_GetElementNodeSize( cFeMesh, element_dI );
		uniqueRelatedColNodes_AllocCount += nElInc;
	}
	
	Journal_DPrintfL( self->debug, 3, "Calculated the max possible number of unique related nodes as %d\n",
			  uniqueRelatedColNodes_AllocCount );
	uniqueRelatedColNodes = Memory_Alloc_Array( Node_DomainIndex, uniqueRelatedColNodes_AllocCount, "uniqueRelatedColNodes" );

	_StiffnessMatrix_CalculatedListOfUniqueRelatedColNodes( self, rowNode_lI, uniqueRelatedColNodes, &uniqueRelatedColNodesCount);
	
	Journal_DPrintfL( self->debug, 3, "Searching the %d unique related col nodes for active dofs\n",
			  uniqueRelatedColNodesCount );
	Stream_Indent( self->debug );
	for ( uniqueRelatedColNode_I = 0; uniqueRelatedColNode_I < uniqueRelatedColNodesCount; uniqueRelatedColNode_I++ ) {
		colNode_dI = uniqueRelatedColNodes[uniqueRelatedColNode_I];

		Journal_DPrintfL( self->debug, 3, "col node_dI %d: has %d dofs\n", colNode_dI, colDofLayout->dofCounts[colNode_dI] );
		Stream_Indent( self->debug );
		for ( colNodeDof_I = 0; colNodeDof_I < colDofLayout->dofCounts[colNode_dI]; colNodeDof_I++ ) {

			currColEqNum = colEqNum->destinationArray[colNode_dI][colNodeDof_I];
			Journal_DPrintfL( self->debug, 3, "dof %d: ", colNodeDof_I );
			if ( currColEqNum == -1 ) {
				Journal_DPrintfL( self->debug, 3, "is a BC.\n" );
			}
			else {
				if ( (currColEqNum >= colEqNum->firstOwnedEqNum ) && 
				     (currColEqNum <= colEqNum->lastOwnedEqNum ) )
				{
					Journal_DPrintfL( self->debug, 3, "is diagonal (eq %d)\n", currColEqNum );
					countTableToAdjust = self->diagonalNonZeroIndices;
				}
				else {
					Journal_DPrintfL( self->debug, 3, "is off-diagonal (eq %d)\n", currColEqNum );
					countTableToAdjust = self->offDiagonalNonZeroIndices;
				}

				for ( currNodeDof_I = 0; currNodeDof_I < self->rowVariable->dofLayout->dofCounts[rowNode_lI]; currNodeDof_I++) {
					if ( -1 != rowEqNum->destinationArray[rowNode_lI][currNodeDof_I] ) {
						currDofMatrixRow = rowEqNum->destinationArray[rowNode_lI][currNodeDof_I] -
							rowEqNum->firstOwnedEqNum;
						/* Because of periodic BCs, the eq num may be lower than the normal
						 * lowest held on this processor, so we need to check this */
						if ( currDofMatrixRow >= self->rowLocalSize ) {	
							Journal_DPrintfL( self->debug, 3, "Found currDofMatRow(=%d) >= self->rowLocalSize(=%d) : for "
									  "rowNode_lI=%d, currMatRow=%d, colNode_dI=%d, colNodeDof_I = %d, "
									  "currNodeDof_I = %d\n", currDofMatrixRow,
									  self->rowLocalSize, rowNode_lI, currMatrixRow,
									  colNode_dI, colNodeDof_I, currNodeDof_I ); 
						}		
						else if ( currDofMatrixRow < 0 ) {	
							Journal_DPrintfL( self->debug, 3, "Found currDofMatRow(=%d) < 0 : for "
									  "rowNode_lI=%d, currMatRow=%d, colNode_dI=%d, colNodeDof_I = %d, "
									  "currNodeDof_I = %d\n", currDofMatrixRow,
									  rowNode_lI, currMatrixRow,
									  colNode_dI, colNodeDof_I, currNodeDof_I ); 
						}		
						else {	
							Journal_DPrintfL( self->debug, 3, "(incrementing app. count at row %d)\n",
									  currDofMatrixRow );

							countTableToAdjust[ currDofMatrixRow ] += 1;
						}
					}
				}
			}	
		}	
		Stream_UnIndent( self->debug );
	}		
	Stream_UnIndent( self->debug );

	Journal_DPrintfL( self->debug, 3, "diagonal count\t%d off diagonal count\t%d\n", 
			  self->diagonalNonZeroIndices[currMatrixRow], self->offDiagonalNonZeroIndices[currMatrixRow]);

	Memory_Free( uniqueRelatedColNodes );

/* 	 TODO: do we need to check that diag is set to at least 1, as the PETSc webpage on MatCreate_MPIAIJ suggests? */

	Stream_UnIndent( self->debug );
}			


void _StiffnessMatrix_CalculatedListOfUniqueRelatedColNodes(
	StiffnessMatrix*	self,
	Node_LocalIndex		rowNode_lI,
	Node_DomainIndex*	uniqueRelatedColNodes,
	Node_Index*		uniqueRelatedColNodesCountPtr )
{
	FeMesh*			rFeMesh = self->rowVariable->feMesh;
	FeMesh*			cFeMesh = self->columnVariable->feMesh;
	Element_Index		rowNodeElement_I = 0;
	Element_DomainIndex	element_dI = 0;
	Node_Index		colElLocalNode_I = 0;
	Node_DomainIndex	colNode_dI = 0;
	Node_Index		uniqueRelatedColNode_I = 0;
	unsigned		nNodeInc, *nodeInc;

	FeMesh_GetNodeElements( rFeMesh, rowNode_lI, self->rowInc );
	nNodeInc = IArray_GetSize( self->rowInc );
	nodeInc = IArray_GetPtr( self->rowInc );
	Journal_DPrintfL( self->debug, 3, "Searching the %d elements this node belongs to for unique related col nodes:\n",
			  nNodeInc );
	
	Stream_Indent( self->debug );
	for ( rowNodeElement_I = 0; rowNodeElement_I < nNodeInc; rowNodeElement_I++ ) {
		unsigned	nElInc, *elInc;

		/* note the dI (domain index) - some of these elements may be in shadow space */
		element_dI = nodeInc[rowNodeElement_I];

		Journal_DPrintfL( self->debug, 3, "rowNodeElement_I: ", rowNodeElement_I );
		Journal_DPrintfL( self->debug, 3, "domain element %d\n", element_dI );
		
		Stream_Indent( self->debug );
		FeMesh_GetElementNodes( cFeMesh, element_dI, self->colInc );
		nElInc = IArray_GetSize( self->colInc );
		elInc = IArray_GetPtr( self->colInc );
		Journal_DPrintfL( self->debug, 3, "Searching the %d column var nodes in this el:\n", nElInc );
		for ( colElLocalNode_I =0; colElLocalNode_I < nElInc; colElLocalNode_I++ ) {
			colNode_dI = elInc[colElLocalNode_I];
			
			Journal_DPrintfL( self->debug, 3, "Col domain node %d: ", colNode_dI );
			for ( uniqueRelatedColNode_I = 0; uniqueRelatedColNode_I < (*uniqueRelatedColNodesCountPtr); uniqueRelatedColNode_I++ )
			{
				if ( colNode_dI == uniqueRelatedColNodes[uniqueRelatedColNode_I] ) {
					Journal_DPrintfL( self->debug, 3, "already in list -> skip to next.\n" );
					break;
				}
			}
			if ( uniqueRelatedColNode_I == (*uniqueRelatedColNodesCountPtr) ) {
				Journal_DPrintfL( self->debug, 3, "is unique so far -> add to list.\n" );
				uniqueRelatedColNodes[uniqueRelatedColNode_I] = colNode_dI;
				(*uniqueRelatedColNodesCountPtr)++;
			}
		}
		Stream_UnIndent( self->debug );
	}
	Stream_UnIndent( self->debug );
}


void StiffnessMatrix_Assemble( void* stiffnessMatrix, Bool bcRemoveQuery, void* _sle, void* _context ) {
	StiffnessMatrix* self = (StiffnessMatrix*)stiffnessMatrix;

	StiffnessMatrix_RefreshMatrix( self );

	Journal_DPrintf( self->debug, "In %s - for matrix \"%s\" - calling the \"%s\" E.P.\n", __func__, self->name,
			 self->assembleStiffnessMatrix->name );
	/* Call the Entry point directly from the base class */
	/* Note that it may be empty: this is deliberate. */
	((FeEntryPoint_AssembleStiffnessMatrix_CallFunction*)EntryPoint_GetRun( self->assembleStiffnessMatrix ))(
		self->assembleStiffnessMatrix,
		self,
		bcRemoveQuery,
		_sle,
		_context );
}


void StiffnessMatrix_GlobalAssembly_General( void* stiffnessMatrix, Bool bcRemoveQuery, void* _sle, void* _context )
{
	StiffnessMatrix*        self = (StiffnessMatrix*)stiffnessMatrix;
	SystemLinearEquations*  sle  = (SystemLinearEquations*)_sle;
	FiniteElementContext*   context  = (FiniteElementContext*)_context;
	
	FeVariable*             feVars[MAX_FE_VARS]; /* Set later */
	Index                   numFeVars = 2;
	Index                   feVar_I;

	/* Location Matrix for a single element - may be a ptr into full LM array */
	Dof_EquationNumber**    elementLM[MAX_FE_VARS];
	/* Number of dofs at an element: required for matrix assembly */
	Dof_Index*              totalDofsThisElement[MAX_FE_VARS];
	Dof_Index*              totalDofsPrevElement[MAX_FE_VARS];
	/* The actual element stiffness matrix values, per element: filled in by this class's entry point */
	double**                elStiffMatToAdd = NULL;

	/* counts and indices used in building per-element locationMatrix and BC info */
	Element_LocalIndex      nLocalElements;
	Element_LocalIndex      element_lI;
	Node_ElementLocalIndex	nodeCountThisEl;
	/* Shortcut to node IDs at a particular element */
	Element_Nodes           nodeIdsThisEl;

	/* related to correction of BCs */
	Bool                    makeBC_Corrections_row = False;
	Bool                    modifiedRHS_Vec_cont = False;
	Dof_Index*              bcLM_Id[ MAX_FE_VARS ] = { NULL, NULL };
	double*                 bcValues[ MAX_FE_VARS ] = { NULL, NULL };
	int                     nBC_NodalDof[ MAX_FE_VARS ] = { 0, 0 };
	double*                 h2Add = NULL;
	
	/* For output printing */
	double                  outputPercentage      = 10;	/* Controls how often to give a status update of assembly progress*/
	int                     outputInterval;
	double                  startTime, totalTime;
	double                  matAddingStart;
	double                  matAddingTime         = 0;
	double                  elStiffMatBuildStart;
	double                  elStiffMatBuildTime   = 0;
	Bool			updateRHS;
	MPI_Comm		comm;

	/* Do some type checking */
	assert( !sle || Stg_CheckType( sle, SystemLinearEquations ) );
	assert( self->rhs || self->transRHS );

	feVars[0] = self->rowVariable;
	feVars[1] = self->columnVariable;

	/* Get communicator. */
	comm = Mesh_GetCommTopology( feVars[0]->feMesh, MT_VERTEX )->mpiComm;
	
	startTime = MPI_Wtime();

	Journal_DPrintf( self->debug, "In %s - for matrix \"%s\"\n", __func__, self->name );
	Stream_IndentBranch( StgFEM_Debug );

	Journal_Firewall( Stg_ObjectList_Count( self->stiffnessMatrixTermList ) != 0,
			  Journal_Register(Error_Type, self->type),
			  "Error in func %s for %s '%s' - No StiffnessMatrixTerms registered.\n", 
			  __func__, self->type, self->name );

	totalDofsThisElement[ROW_VAR] = Memory_Alloc( Dof_Index, "el nodal dofs" );
	totalDofsPrevElement[ROW_VAR] = Memory_Alloc( Dof_Index, "el nodal dofs (prev element)" );
	*totalDofsPrevElement[ROW_VAR] = 0;
	/* check if row variable and col variable are the same */
	if ( self->rowVariable == self->columnVariable ) {
		numFeVars = 1;
		Journal_DPrintfL( self->debug, 2, "Detected both row and column FeVariable to assemble over are \"%s\", "
				  "so only processing once per element.\n", self->rowVariable->name );
		/* since Row and Col FeVars are same, set LM_col and totalDofsThisElement to be same as row ones */
		totalDofsThisElement[COL_VAR] = totalDofsThisElement[ROW_VAR];	
		totalDofsPrevElement[COL_VAR] = totalDofsPrevElement[ROW_VAR];	
	}
	else {
		Journal_DPrintfL( self->debug, 2, "Since row FeVariable \"%s\" and column FeVariable \"%s\" to assemble over "
				  "are different, processing both per element.\n", self->rowVariable->name, self->columnVariable->name );
		totalDofsThisElement[COL_VAR] = Memory_Alloc( Dof_Index, "el nodal dofs (col)" );
		totalDofsPrevElement[COL_VAR] = Memory_Alloc( Dof_Index, "el nodal dofs (col) (prev element)" );
		*totalDofsPrevElement[COL_VAR] = 0;
	}
	
	/* Assumes that both row and col variables have same number of variables */
	nLocalElements = FeMesh_GetElementLocalSize( feVars[ROW_VAR]->feMesh );
	
	/* Initialise matrix */
	Matrix_Zero( self->matrix );
	outputInterval = (int)( (outputPercentage/100.0)*(double)(nLocalElements) );
	if( outputInterval == 0 ) { outputInterval = nLocalElements; }

	for( element_lI = 0; element_lI < nLocalElements; element_lI++ ) {  
		
		/* This loop is how we handle the possiblity of different row-column variables: if both variables are
		 * the same, then numFeVars is set to 1 and the loop only progresses through once.
		 -- PatrickSunter 13 September 2004 */
		for ( feVar_I = 0; feVar_I < numFeVars; feVar_I++ ) {
			FeEquationNumber* 		feEqNum = feVars[feVar_I]->eqNum;
			DofLayout*			dofLayout = feVars[feVar_I]->dofLayout;
			Dof_Index			dofCountLastNode;
			unsigned			nDofsThisEl;
			unsigned			n_i;

			/* Get the local node ids */
			FeMesh_GetElementNodes( feVars[feVar_I]->feMesh, element_lI, 
						self->rowInc );
			nodeCountThisEl = IArray_GetSize( self->rowInc );
			nodeIdsThisEl = IArray_GetPtr( self->rowInc );

			/* Set value of elementLM: will automatically use large one if built */
			elementLM[feVar_I] = FeEquationNumber_BuildOneElementLocationMatrix( feEqNum, element_lI );

			/* work out number of dofs at the node, based on LM */
			dofCountLastNode = dofLayout->dofCounts[nodeIdsThisEl[nodeCountThisEl-1]];
			nDofsThisEl = dofLayout->dofCounts[nodeIdsThisEl[0]];
			for( n_i = 1; n_i < nodeCountThisEl; n_i++ )
				nDofsThisEl += dofLayout->dofCounts[nodeIdsThisEl[n_i]];
			*totalDofsThisElement[feVar_I] = nDofsThisEl;
/*
 *totalDofsThisElement[feVar_I] = &elementLM[feVar_I][nodeCountThisEl-1][dofCountLastNode-1]
 - &elementLM[feVar_I][0][0] + 1;
*/

			if ( *totalDofsThisElement[feVar_I] > *totalDofsPrevElement[feVar_I] ) {
				Journal_DPrintfL( self->debug, 2, "Reallocating bcLM_Id and bcValues for fe var \"%s\" "
						  "to size %d\n", feVars[feVar_I]->name, *totalDofsThisElement[feVar_I] );

				if ( bcLM_Id[feVar_I] ) Memory_Free( bcLM_Id[feVar_I] );
				bcLM_Id[feVar_I] = Memory_Alloc_Array( Dof_Index, *totalDofsThisElement[feVar_I], "bcLM_Id" );

				if ( bcValues[feVar_I] ) Memory_Free( bcValues[feVar_I] );
				bcValues[feVar_I] = Memory_Alloc_Array( double, *totalDofsThisElement[feVar_I], "bcValues" );
			}
		}

		/* Reallocate el stiff mat and other arrays if necessary */
		if ( (*totalDofsThisElement[ROW_VAR] != *totalDofsPrevElement[ROW_VAR]) ||
		     (*totalDofsThisElement[COL_VAR] != *totalDofsPrevElement[COL_VAR]) )
		{
			if (h2Add) Memory_Free( h2Add );
			Journal_DPrintfL( self->debug, 2, "Reallocating h2Add to size %d\n",
					  *totalDofsThisElement[COL_VAR] ); 
			h2Add = Memory_Alloc_Array( double, *totalDofsThisElement[COL_VAR], "h2Add" );
			
			if (elStiffMatToAdd) Memory_Free( elStiffMatToAdd );
			Journal_DPrintfL( self->debug, 2, "Reallocating elStiffMatToAdd to size %d*%d\n",
					  *totalDofsThisElement[ROW_VAR], *totalDofsThisElement[COL_VAR] ); 
			elStiffMatToAdd = Memory_Alloc_2DArray( double, *totalDofsThisElement[ROW_VAR], *totalDofsThisElement[COL_VAR], "elStiffMatToAdd" );
		}

		/* Initialise the elStiffMat to zero */
		/* Note we have to dereference the ptr once ... so we don't clobber the ptrs, just the values */
		memset( elStiffMatToAdd[0], 0, sizeof(double) * *totalDofsThisElement[ROW_VAR] * *totalDofsThisElement[COL_VAR] );

		/* Assemble this element's element stiffness matrix: call the entry point */
		elStiffMatBuildStart = MPI_Wtime();
		StiffnessMatrix_AssembleElement( self, element_lI, sle, context, elStiffMatToAdd );
		elStiffMatBuildTime += MPI_Wtime() - elStiffMatBuildStart;
		if ( False == self->allowZeroElementContributions ) {
			StiffnessMatrix_CheckElementAssembly( self, element_lI, elStiffMatToAdd, *totalDofsThisElement[ROW_VAR],
							      *totalDofsThisElement[COL_VAR] );
		}	

		/* This loop is how we handle the possiblity of different row-column variables: if both variables are
		 * the same, then numFeVars is set to 1 and the loop only progresses through once.
		 -- PatrickSunter 13 September 2004 */
		for ( feVar_I = 0; feVar_I < numFeVars; feVar_I++ ) {
			/* Reset from previous calculation */
			nBC_NodalDof[feVar_I] = 0;
			makeBC_Corrections_row = False;

			/* Check if we want to build table of corrected BC info for the matrix */
			if( bcRemoveQuery == True ) {
				unsigned	nNodeInc, *nodeInc;

				FeMesh_GetElementNodes( feVars[feVar_I]->feMesh, element_lI, 
							self->rowInc );
				nNodeInc = IArray_GetSize( self->rowInc );
				nodeInc = IArray_GetPtr( self->rowInc );
				_StiffnessMatrix_UpdateBC_CorrectionTables(
					self,
					feVars[feVar_I]->eqNum, 
					feVars[feVar_I]->dofLayout,
					elementLM[feVar_I],
					nNodeInc, 
					nodeInc, 
					bcLM_Id[feVar_I],
					bcValues[feVar_I],
					&(nBC_NodalDof[feVar_I]) );
			}
		}

		/* If there is only one feVar, use the same LM for both row and col insertion */
		if ( numFeVars == 1 ) elementLM[COL_VAR] = elementLM[ROW_VAR];

/*
  #if DEBUG
  if ( Stream_IsPrintableLevel( self->debug, 3 ) ) {
  Journal_DPrintf( self->debug, "Handling Element %d:\n", element_lI );
  for ( feVar_I = 0; feVar_I < numFeVars; feVar_I++ ) {
  FeEquationNumber_PrintElementLocationMatrix( feVars[feVar_I]->eqNum,
  elementLM[feVar_I], element_lI, self->debug );
				
  }
			
  Journal_DPrintf( self->debug, "El stiff Mat about to be added:\n" );
  _StiffnessMatrix_PrintElementStiffnessMatrix( self, element_lI, elementLM[ROW_VAR],
  elementLM[COL_VAR], elStiffMatToAdd );
  }	
  #endif
*/


		/*
		** Need to make corrections to the RHS before filling in the matrix.  This is because we may modify the
		** element matrix values before they get there.
		*/

		/* Set flag for making corrcetions to rhs from bc's */
		if( nBC_NodalDof[ROW_VAR] != 0 ) {
			makeBC_Corrections_row = True;
			modifiedRHS_Vec_cont = True;	/* this guy only needs to be toggled once */
		}
		
		if( makeBC_Corrections_row == True && bcRemoveQuery ) {
			_StiffnessMatrix_CorrectForceVectorWithOneElementsBoundaryConditions(
				self,
				elementLM,
				h2Add,
				elStiffMatToAdd,
				totalDofsThisElement,
				bcLM_Id,
				bcValues,
				nBC_NodalDof[ROW_VAR], 
				element_lI );
		}

		/*
		** If not keeping BCs in, we may need to zero some of these element values, or set them to one.
		*/

		if( (!self->rowVariable || !self->rowVariable->eqNum->removeBCs) && 
		    (!self->columnVariable || !self->columnVariable->eqNum->removeBCs) )
		{
			FeEquationNumber*	rowEqNum = self->rowVariable->eqNum;
			FeEquationNumber*	colEqNum = self->columnVariable->eqNum;
			unsigned		nRows = *totalDofsThisElement[ROW_VAR];
			unsigned		row_i;

			/* Loop over elementLM rows */
			for( row_i = 0; row_i < nRows; row_i++ ) {
				unsigned	rowEqInd = elementLM[ROW_VAR][0][row_i] - rowEqNum->firstOwnedEqNum;

				/* If row is bc */
				if( IndexSet_IsMember( rowEqNum->bcEqNums, rowEqInd ) ) {
					unsigned	nCols = *totalDofsThisElement[COL_VAR];
					unsigned	col_i;

					/* Loop over elementLM cols */
					for( col_i = 0; col_i < nCols; col_i++ ) {
						unsigned	colEqInd = elementLM[COL_VAR][0][col_i] - colEqNum->firstOwnedEqNum;

						/* If col is bc */
						if( IndexSet_IsMember( colEqNum->bcEqNums, colEqInd ) ) {
							elStiffMatToAdd[0][row_i * nCols + col_i] = 0.0;
						}
						else {
							elStiffMatToAdd[0][row_i * nCols + col_i] = 0.0;
						}
					}
				}
				else {
					unsigned	nCols = *totalDofsThisElement[COL_VAR];
					unsigned	col_i;

					for( col_i = 0; col_i < nCols; col_i++ ) {
						unsigned	colEqInd = elementLM[COL_VAR][0][col_i];

						/* If col is bc */
						if( IndexSet_IsMember( colEqNum->bcEqNums, colEqInd ) ) {
							elStiffMatToAdd[0][row_i * nCols + col_i] = 0.0;
						}
					}
				}
			}
		}

		/* Add to the global matrix. */
		matAddingStart = MPI_Wtime();
		Matrix_AddEntries( self->matrix, *totalDofsThisElement[ROW_VAR], (Index *)(elementLM[ROW_VAR][0]), *totalDofsThisElement[COL_VAR],
				   (Index *)(elementLM[COL_VAR][0]), elStiffMatToAdd[0] );
		matAddingTime += MPI_Wtime() - matAddingStart;	
		
		
#if DEBUG
		if( element_lI % outputInterval == 0 ) {
			Journal_DPrintfL( self->debug, 2, "done %d percent of global element stiffness assembly (general) \n",
					  (int)(100.0*((double)element_lI/(double)nLocalElements)) );
		}
#endif

		for ( feVar_I = 0; feVar_I < numFeVars; feVar_I++ ) {
			/* Save the number of dofs in this element */
			*totalDofsPrevElement[feVar_I] = *totalDofsThisElement[feVar_I];
			/* If we haven't built the big LM for all elements, free the temporary one */
			if ( False == feVars[feVar_I]->eqNum->locationMatrixBuilt ) {
				Memory_Free( elementLM[feVar_I] );
			}	
		}	
	}

/* 	//////////////////////////////////////////////////////// */
		
	Matrix_AssemblyBegin( self->matrix );
	Matrix_AssemblyEnd( self->matrix );
	
	MPI_Allreduce( &modifiedRHS_Vec_cont, &updateRHS, 1, MPI_UNSIGNED, MPI_LOR, comm );
	if( updateRHS == True && bcRemoveQuery ) {
		/* stuff was submitted to vec (at least once) and needs to be assembled */
		Vector_AssemblyBegin( self->rhs->vector );
		Vector_AssemblyEnd( self->rhs->vector );
	}

/* 	//////////////////////////////////////////////////////// */

	/*
	** If keeping BCs in, modify the RHS to reflect BCs.
	*/

	if( (!self->rowVariable || !self->rowVariable->eqNum->removeBCs) && 
	    (!self->columnVariable || !self->columnVariable->eqNum->removeBCs) )
	{
		FeEquationNumber*	colEqNum = self->columnVariable->eqNum;
		FeMesh*			feMesh = self->columnVariable->feMesh;
		DofLayout*		dofLayout = self->columnVariable->dofLayout;
		unsigned		nNodes = FeMesh_GetNodeLocalSize( feMesh );
		unsigned		node_i;

		/* What to do if they aren't the same? */
		assert( self->rowVariable == self->columnVariable );

		for( node_i = 0; node_i < nNodes; node_i++ ) {
			unsigned	nDofs = dofLayout->dofCounts[node_i];
			unsigned	dof_i;

			for( dof_i = 0; dof_i < nDofs; dof_i++ ) {
				unsigned	eqInd = colEqNum->destinationArray[node_i][dof_i];

				if( IndexSet_IsMember( colEqNum->bcEqNums, eqInd - colEqNum->firstOwnedEqNum ) ) {
					double	bcVal = DofLayout_GetValueDouble( dofLayout, node_i, dof_i );
					double	one = 1.0;

					Vector_InsertEntries( self->rhs->vector, 1, &eqInd, &bcVal );
					Matrix_InsertEntries( self->matrix, 1, &eqInd, 1, &eqInd, &one );
				}
			}
		}

		/* Need to reassemble. */
		Matrix_AssemblyBegin( self->matrix );
		Matrix_AssemblyEnd( self->matrix );	
		Vector_AssemblyBegin( self->rhs->vector );
		Vector_AssemblyEnd( self->rhs->vector ); 
	}

/*
  #if DEBUG
  if ( Stream_IsPrintableLevel( self->debug, 3 ) ) {
  Journal_DPrintf( self->debug, "Completely built matrix %s is:\n", self->name );
  Matrix_View( self->matrix, self->debug );
  Journal_DPrintf( self->debug, "Updated RHS vector %s is now:\n", self->rhs->name );
  Vector_View( self->rhs->vector, self->debug );
  }
  #endif
*/

	for ( feVar_I = 0; feVar_I < numFeVars; feVar_I++ ) {
		Memory_Free( totalDofsThisElement[feVar_I] );
		Memory_Free( totalDofsPrevElement[feVar_I] );
		Memory_Free( bcLM_Id[feVar_I] );
		Memory_Free( bcValues[feVar_I] );
	}
	Memory_Free( elStiffMatToAdd );
	Memory_Free( h2Add );

	totalTime = MPI_Wtime() - startTime;
	Journal_DPrintfL( self->debug, 2, "Total time used by %s: %.5g s\n", __func__, totalTime );
	Journal_DPrintfL( self->debug, 2, "Of which %.5g s spend building elStiffMats\n", elStiffMatBuildTime );
	Journal_DPrintfL( self->debug, 2, "And %.5g s spend adding to global matrix\n", matAddingTime );
	Stream_UnIndentBranch( StgFEM_Debug );
}

void StiffnessMatrix_NewAssemble( void* stiffnessMatrix, Bool removeBCs, void* _sle, void* _context ) {
	StiffnessMatrix*		self = (StiffnessMatrix*)stiffnessMatrix;
	SystemLinearEquations*		sle = (SystemLinearEquations*)_sle;
	FeVariable			*rowVar, *colVar;
	FeMesh				*rowMesh, *colMesh;
	FeEquationNumber		*rowEqNum, *colEqNum;
	DofLayout			*rowDofs, *colDofs;
	unsigned			nRowEls;
	unsigned			nRowNodes, *rowNodes;
	unsigned			nColNodes, *colNodes;
	unsigned			maxDofs, maxRCDofs, nDofs, nRowDofs, nColDofs;
	double**			elStiffMat;
	double*				bcVals;
	Matrix*				matrix;
	Vector				*vector, *transVector;
	unsigned			e_i, n_i;

	assert( self && Stg_CheckType( self, StiffnessMatrix ) );

	rowVar = self->rowVariable;
	colVar = self->columnVariable ? self->columnVariable : rowVar;
	rowEqNum = rowVar->eqNum;
	colEqNum = colVar->eqNum;
	rowMesh = rowVar->feMesh;
	colMesh = colVar->feMesh;
	rowDofs = rowVar->dofLayout;
	colDofs = colVar->dofLayout;
	nRowEls = FeMesh_GetElementLocalSize( rowMesh );
	assert( (rowVar == colVar) ? !self->transRHS : 1 );

	matrix = self->matrix;
	vector = self->rhs ? self->rhs->vector : NULL;
	transVector = self->transRHS ? self->transRHS->vector : NULL;
	elStiffMat = NULL;
	bcVals = NULL;
	maxDofs = 0;

	/* Begin assembling each element. */
	for( e_i = 0; e_i < nRowEls; e_i++ ) {
		FeMesh_GetElementNodes( rowMesh, e_i, self->rowInc );
		nRowNodes = IArray_GetSize( self->rowInc );
		rowNodes = IArray_GetPtr( self->rowInc );
		FeMesh_GetElementNodes( colMesh, e_i, self->colInc );
		nColNodes = IArray_GetSize( self->colInc );
		colNodes = IArray_GetPtr( self->colInc );

		/* Do we need more space to assemble this element? */
		nRowDofs = 0;
		for( n_i = 0; n_i < nRowNodes; n_i++ )
			nRowDofs += rowDofs->dofCounts[rowNodes[n_i]];
		nColDofs = 0;
		for( n_i = 0; n_i < nColNodes; n_i++ )
			nColDofs += colDofs->dofCounts[colNodes[n_i]];
		nDofs = nRowDofs * nColDofs;
		self->nRowDofs = nRowDofs;
		self->nColDofs = nColDofs;
		if( nDofs > maxDofs ) {
			maxRCDofs = (nRowDofs > nColDofs) ? nRowDofs : nColDofs;
			elStiffMat = ReallocArray2D( elStiffMat, double, nRowDofs, nColDofs );
			bcVals = ReallocArray( bcVals, double, maxRCDofs );
			maxDofs = nDofs;
			self->elStiffMat = elStiffMat;
			self->bcVals = bcVals;
		}

		/* Assemble the element. */
		memset( elStiffMat[0], 0, nDofs * sizeof(double) );
		StiffnessMatrix_AssembleElement( self, e_i, sle, _context, elStiffMat );

		/* Correct for BCs providing I'm not keeping them in. */
		if( vector ) {
			memset( bcVals, 0, nRowDofs * sizeof(double) );
			Assembler_LoopMatrixElement( self->bcAsm, e_i );
			Vector_AddEntries( vector, nRowDofs, (unsigned*)rowEqNum->locationMatrix[e_i][0], bcVals );
		}
		if( transVector ) {
			memset( bcVals, 0, nColDofs * sizeof(double) );
			Assembler_LoopMatrixElement( self->transBCAsm, e_i );
			Vector_AddEntries( transVector, nColDofs, (unsigned*)colEqNum->locationMatrix[e_i][0], bcVals );
		}

		/* If keeping BCs in, zero corresponding entries in the element stiffness matrix. */
		if( !rowEqNum->removeBCs || !colEqNum->removeBCs )
			Assembler_LoopMatrixElement( self->zeroBCsAsm, e_i );

		/* Add to stiffness matrix. */
		Matrix_AddEntries( matrix, 
				   nRowDofs, (unsigned*)rowEqNum->locationMatrix[e_i][0], 
				   nColDofs, (unsigned*)colEqNum->locationMatrix[e_i][0], 
				   elStiffMat[0] );
	}

	FreeArray( elStiffMat );
	FreeArray( bcVals );

	/* If keeping BCs in and rows and columnns use the same variable, put ones in all BC'd diagonals. */
	if( !colEqNum->removeBCs && rowVar == colVar )
		Assembler_LoopMatrixDiagonal( self->diagBCsAsm );

	/* Reassemble the matrix and vectors. */
	Matrix_AssemblyBegin( matrix );
	Matrix_AssemblyEnd( matrix );
	if( vector ) {
		Vector_AssemblyBegin( vector );
		Vector_AssemblyEnd( vector );
	}
	if( transVector) {
		Vector_AssemblyBegin( transVector );
		Vector_AssemblyEnd( transVector );
	}

	Matrix_AssemblyBegin( matrix );
	Matrix_AssemblyEnd( matrix );
}

void StiffnessMatrix_ShellAssembly( void* stiffnessMatrix, Bool removeBCs, void* data ) {
	StiffnessMatrix*	self = (StiffnessMatrix*)stiffnessMatrix;
	Vector*			rhs;
	FeVariable		*rowVar, *colVar;
	FeMesh			*rowMesh, *colMesh;
	FeEquationNumber	*rowEqNum, *colEqNum;
	DofLayout		*rowDofs, *colDofs;
	unsigned		nRowEls;
	unsigned		nRowNodes, *rowNodes;
	unsigned		nColNodes, *colNodes;
	unsigned		maxDofs, nDofs, nRowDofs, nColDofs;
	double**		elStiffMat;
	double*			values;
	unsigned*		indices;
	unsigned		curRow, curCol;
	unsigned		colEq;
	double			bc;
	unsigned		e_i, n_i, n_j, dof_i, dof_j;

	assert( self && Stg_CheckType( self, StiffnessMatrix ) );

	/* The whole point of this routine is to remove the BCs. */
	if( !removeBCs )
		return;

	rhs = self->rhs->vector;
	rowVar = self->rowVariable;
	colVar = self->columnVariable ? self->columnVariable : rowVar;
	if( rowVar != colVar ) {
		FeVariable* 	tmp = rowVar;
		rowVar = colVar;
		colVar = tmp;
	}
	rowEqNum = rowVar->eqNum;
	colEqNum = colVar->eqNum;
	rowMesh = rowVar->feMesh;
	colMesh = colVar->feMesh;
	rowDofs = rowVar->dofLayout;
	colDofs = colVar->dofLayout;
	nRowEls = FeMesh_GetElementLocalSize( rowMesh );
	elStiffMat = NULL;
	values = NULL;
	indices = NULL;
	maxDofs = 0;

	for( e_i = 0; e_i < nRowEls; e_i++ ) {
		FeMesh_GetElementNodes( rowMesh, e_i, self->rowInc );
		nRowNodes = IArray_GetSize( self->rowInc );
		rowNodes = IArray_GetPtr( self->rowInc );
		FeMesh_GetElementNodes( colMesh, e_i, self->colInc );
		nColNodes = IArray_GetSize( self->colInc );
		colNodes = IArray_GetPtr( self->colInc );

		/* If none of the column equations on this element have BCs then skip it. */
		for( n_i = 0; n_i < nRowNodes; n_i++ ) {
			for( dof_i = 0; dof_i < rowDofs->dofCounts[rowNodes[n_i]]; dof_i++ ) {
				if( rowEqNum->locationMatrix[e_i][n_i][dof_i] == (unsigned)-1 )
					continue;
				for( n_j = 0; n_j < nColNodes; n_j++ ) {
					for( dof_j = 0; dof_j < colDofs->dofCounts[colNodes[n_j]]; dof_j++ ) {
						if( colEqNum->locationMatrix[e_i][n_j][dof_j] == (unsigned)-1 )
							break;
					}
					if( dof_j < colDofs->dofCounts[colNodes[n_j]] )
						break;
				}
				if( n_j < nColNodes )
					break;
			}
			if( dof_i < rowDofs->dofCounts[rowNodes[n_i]] )
				break;
		}
		if( n_i == nRowNodes )
			continue;

		/* Do we need more space to assemble this element? */
		nRowDofs = 0;
		for( n_i = 0; n_i < nRowNodes; n_i++ )
			nRowDofs += rowDofs->dofCounts[rowNodes[n_i]];
		nColDofs = 0;
		for( n_i = 0; n_i < nColNodes; n_i++ )
			nColDofs += colDofs->dofCounts[colNodes[n_i]];
		nDofs = nRowDofs * nColDofs;
		if( nDofs > maxDofs ) {
#if 0
			elStiffMat = ReallocArray2D( elStiffMat, double, nRowDofs, nColDofs );
#endif
			values = ReallocArray( values, double, nDofs );
			indices = ReallocArray( indices, unsigned, nDofs );
			maxDofs = nDofs;
		}

#if 0
		/* Assemble the element. */
		memset( &elStiffMat[0][0], 0, nDofs * sizeof(double) );
		StiffnessMatrix_AssembleElement( self, e_i, sle, elStiffMat );
#endif

		elStiffMat = ((PETScShellMatrix*)self->matrix)->elStiffMat;	assert( elStiffMat );

		/* Update the force vector with BCs. */
		curRow = 0;
		memset( values, 0, nDofs * sizeof(double) );
		for( n_i = 0; n_i < nRowNodes; n_i++ ) {
			for( dof_i = 0; dof_i < rowDofs->dofCounts[rowNodes[n_i]]; dof_i++ ) {
				indices[curRow] = rowEqNum->locationMatrix[e_i][n_i][dof_i];
				if( indices[curRow] == (unsigned)-1 ) {
					curRow++;
					continue;
				}

				curCol = 0;
				for( n_j = 0; n_j < nColNodes; n_j++ ) {
					for( dof_j = 0; dof_j < colDofs->dofCounts[colNodes[n_j]]; dof_j++ ) {
						colEq = colEqNum->locationMatrix[e_i][n_j][dof_j];
						if( colEq != (unsigned)-1 ) {
							curCol++;
							continue;
						}

						bc = DofLayout_GetValueDouble( colDofs, colNodes[n_j], dof_j );
						values[curRow] -= elStiffMat[curRow][curCol] * bc;

						curCol++;
					}
				}

				curRow++;
			}
		}

		Vector_AddEntries( rhs, curRow, indices, values );
	}

	FreeArray( values );
	FreeArray( indices );

	Vector_AssemblyBegin( rhs );
	Vector_AssemblyEnd( rhs );
}



/* +++ PRIVATE FUNCTIONS +++ */
	
void _StiffnessMatrix_UpdateBC_CorrectionTables(
	StiffnessMatrix*	self,
	FeEquationNumber*	eqNum, 
	DofLayout*		dofLayout,
	Dof_EquationNumber**	elementLM,
	Node_ElementLocalIndex	nodeCountThisEl,
	Element_Nodes		nodeIdsThisEl,
	Dof_Index*		bcLM_Id,
	double*			bcValues,
	int*			nBC_NodalDofPtr )
{
	Node_ElementLocalIndex		node_elLocalI = 0;
	Node_LocalIndex			node_lI = 0;
	Dof_Index*			dofCounts = dofLayout->dofCounts;
	Dof_Index			dofCountThisNode=0;
	Dof_Index			dof_nodeLocalI=0;
	unsigned			pos = 0;

	for( node_elLocalI = 0; node_elLocalI < nodeCountThisEl; node_elLocalI++ ) {
		node_lI = nodeIdsThisEl[node_elLocalI];
		dofCountThisNode = dofCounts[node_lI];
		
		for( dof_nodeLocalI = 0; dof_nodeLocalI < dofCountThisNode; dof_nodeLocalI++ ) {
			Bool	isBC = False;

			/* Can only use 'elementLM' if FeEquationNumber has been told to remove BCs.  Otherwise
			   we'll need to determine if the VariableCondition has a value specified for this 
			   node/dof. - Luke */
			if( elementLM[node_elLocalI][dof_nodeLocalI] != (unsigned)-1 ) {
				unsigned	lEqNum = elementLM[node_elLocalI][dof_nodeLocalI] - eqNum->_lowestLocalEqNum;

				if( eqNum->bcEqNums && IndexSet_IsMember( eqNum->bcEqNums, lEqNum ) ) {
					isBC = True;
				}
			}
			else {
				isBC = True;
			}

			if ( isBC ) {
				/* offset into the elementStiffness matrix */
				bcLM_Id[ *nBC_NodalDofPtr ] = pos;
/*
  bcLM_Id[ *nBC_NodalDofPtr ] = &elementLM[node_elLocalI][dof_nodeLocalI] - &elementLM[0][0];
*/

/* 				 get bc values from the bc_layout  */
				bcValues[ *nBC_NodalDofPtr ]  = DofLayout_GetValueDouble( dofLayout, node_lI, dof_nodeLocalI );

				Journal_DPrintfL( self->debug, 3, "bcValues[%d]: at &LM[0][0] + %d=%d, is %f\n",
						  *nBC_NodalDofPtr, bcLM_Id[ *nBC_NodalDofPtr ],
						  elementLM[0][ bcLM_Id[ *nBC_NodalDofPtr ] ],
						  bcValues[ *nBC_NodalDofPtr ] );

				(*nBC_NodalDofPtr)++;
			}

			/* Move to next element stiffness matrix entry. */
			pos++;
		}
	}
}	


void _StiffnessMatrix_CorrectForceVectorWithOneElementsBoundaryConditions(
	StiffnessMatrix*	self,
	Dof_EquationNumber** 	elementLM[MAX_FE_VARS],
	double*			h2Add,
	double** 		elStiffMatToAdd,
	Dof_Index*		totalDofsThisElement[MAX_FE_VARS],
	Dof_Index*		bcLM_Id[MAX_FE_VARS],
	double*			bcValues[MAX_FE_VARS],
	int			nBC_NodalDof_Row, 
	unsigned		elementInd /* NEW ONE */ )
{
#if 0
	int 		rowEqId;
	double		bc_value;
	Dof_Index	colDof_elLocalI;
#endif

	memset( h2Add, 0, (*totalDofsThisElement[COL_VAR]) * sizeof(double) );

	/*
	** Something fishy is up with BC corrections, adding this for now.
	*/

	if( self->rowVariable != self->columnVariable )	{
		Mesh			*rowMesh, *colMesh;
		FeEquationNumber	*rowEqNum, *colEqNum;
		DofLayout		*rowDofs, *colDofs;
		unsigned		nRowElNodes, *rowElNodes, nColElNodes, *colElNodes;
		unsigned		nRowDofs, nColDofs;
		unsigned		dofI, dofJ, elIndI, elIndJ;
		double			bcValue;
		unsigned		n_i, n_j, d_i, d_j;

		rowMesh = (Mesh*)self->rowVariable->feMesh;
		colMesh = (Mesh*)self->columnVariable->feMesh;
		rowEqNum = self->rowVariable->eqNum;
		colEqNum = self->columnVariable->eqNum;
		rowDofs = rowEqNum->dofLayout;
		colDofs = colEqNum->dofLayout;
		Mesh_GetIncidence( rowMesh, Mesh_GetDimSize( rowMesh ), elementInd, MT_VERTEX, self->rowInc );
		nRowElNodes = IArray_GetSize( self->rowInc );
		rowElNodes = IArray_GetPtr( self->rowInc );
		Mesh_GetIncidence( colMesh, Mesh_GetDimSize( colMesh ), elementInd, MT_VERTEX, self->colInc );
		nColElNodes = IArray_GetSize( self->colInc );
		colElNodes = IArray_GetPtr( self->colInc );
		nRowDofs = rowDofs->dofCounts[0];
		nColDofs = colDofs->dofCounts[0];

		for( n_i = 0; n_i < nColElNodes; n_i++ ) {
			for( d_i = 0; d_i < nColDofs; d_i++ ) {
				dofI = colEqNum->locationMatrix[elementInd][n_i][d_i];
				if( dofI == -1 )
					continue;

				elIndI = n_i * nColDofs + d_i;
				for( n_j = 0; n_j < nRowElNodes; n_j++ ) {
					for( d_j = 0; d_j < nRowDofs; d_j++ ) {
						dofJ = rowEqNum->locationMatrix[elementInd][n_j][d_j];
						if( dofJ != -1 )
							continue;

						elIndJ = n_j * nRowDofs + d_j;
						bcValue = DofLayout_GetValueDouble( rowDofs, rowElNodes[n_j], d_j );
						h2Add[elIndI] -= elStiffMatToAdd[elIndJ][elIndI] * bcValue;
					}
				}
			}
		}
	}
	else {
		Mesh			*rowMesh, *colMesh;
		FeEquationNumber	*rowEqNum, *colEqNum;
		DofLayout		*rowDofs, *colDofs;
		unsigned		nRowElNodes, *rowElNodes, nColElNodes, *colElNodes;
		unsigned		nRowDofs, nColDofs;
		unsigned		dofI, dofJ, elIndI, elIndJ;
		double			bcValue;
		unsigned		n_i, n_j, d_i, d_j;

		rowMesh = (Mesh*)self->rowVariable->feMesh;
		colMesh = (Mesh*)self->columnVariable->feMesh;
		rowEqNum = self->rowVariable->eqNum;
		colEqNum = self->columnVariable->eqNum;
		rowDofs = rowEqNum->dofLayout;
		colDofs = colEqNum->dofLayout;
		Mesh_GetIncidence( rowMesh, Mesh_GetDimSize( rowMesh ), elementInd, MT_VERTEX, self->rowInc );
		nRowElNodes = IArray_GetSize( self->rowInc );
		rowElNodes = IArray_GetPtr( self->rowInc );
		Mesh_GetIncidence( colMesh, Mesh_GetDimSize( colMesh ), elementInd, MT_VERTEX, self->colInc );
		nColElNodes = IArray_GetSize( self->colInc );
		colElNodes = IArray_GetPtr( self->colInc );
		nRowDofs = rowDofs->dofCounts[0];
		nColDofs = colDofs->dofCounts[0];

		for( n_i = 0; n_i < nRowElNodes; n_i++ ) {
			for( d_i = 0; d_i < nRowDofs; d_i++ ) {
				dofI = rowEqNum->locationMatrix[elementInd][n_i][d_i];
				if( dofI == -1 )
					continue;

				elIndI = n_i * nRowDofs + d_i;
				for( n_j = 0; n_j < nColElNodes; n_j++ ) {
					for( d_j = 0; d_j < nColDofs; d_j++ ) {
						dofJ = colEqNum->locationMatrix[elementInd][n_j][d_j];
						if( dofJ != -1 )
							continue;

						elIndJ = n_j * nColDofs + d_j;
						bcValue = DofLayout_GetValueDouble( colDofs, colElNodes[n_j], d_j );
						h2Add[elIndI] -= elStiffMatToAdd[elIndI][elIndJ] * bcValue;
					}
				}
			}
		}
	}

#if 0
	for( colDof_elLocalI=0; colDof_elLocalI < *totalDofsThisElement[COL_VAR]; colDof_elLocalI++ ) {
		double		sum = 0.0;

		for( bcDof_I=0; bcDof_I < nBC_NodalDof_Row; bcDof_I++ ) {
			rowEqId = bcLM_Id[ROW_VAR][bcDof_I];
			bc_value = bcValues[ROW_VAR][bcDof_I];
/* 			printf("bc_value = %f \n",bc_value );  */
			
			/* this index is gets us to the right */
			sum = sum - ( elStiffMatToAdd[rowEqId][colDof_elLocalI] * bc_value );
		}
		h2Add[colDof_elLocalI] = sum;
	}
#endif

	Vector_AddEntries( self->rhs->vector, *totalDofsThisElement[COL_VAR], (Index *)(elementLM[COL_VAR][0]), h2Add );
	
	/* assume that K is symetric, so corrections are made with KTrans.
	   this allows us to use this func with G, when we want velocity
	   corrections from GTrans to appear in H.
	*/

	/*
	  for( bcDof_I=0; bcDof_I < nBC_NodalDof_row; bcDof_I++ ) {
	  rowEqId = bcLM_Id[ROW_VAR][bcDof_I];
	  bc_value = bcValues[ROW_VAR][bcDof_I];
	  printf("bc_value = %f \n",bc_value ); 
		
	  printf("bc value = %f \n",bc_value );
	  correct = 0;
	  for( colDof_elLocalI=0; colDof_elLocalI< *totalDofsThisElement[COL_VAR]; colDof_elLocalI++ ) {
	  h2Add[correct] = h2Add[correct]
	  - elStiffMatToAdd[rowEqId* (*totalDofsThisElement[COL_VAR]) + colDof_elLocalI] * bc_value;
	  hIdx[correct] = elementLM[COL_VAR][0][ colDof_elLocalI ];
			
	  correct = correct + 1;
	  }
	  }
	  Vector_AddTo( self->rhs->vector, correct, hIdx, h2Add );
	*/
	
/* 	 does not work */
	/*
	  for( rowDof_elLocalI=0; rowDof_elLocalI < *totalDofsThisElement[ROW_VAR]; rowDof_elLocalI++ ) {
	  double sum = 0.0;
	  for( rowDof_elLocalI=0; rowDof_elLocalI<nBC_NodalDof[ROW_VAR]; rowDof_elLocalI++ ) {
	  colEqId = bcLM_Id[COL_VAR][rowDof_elLocalI];
	  bc_value = bcValues[COL_VAR][rowDof_elLocalI];
			
	  sum = sum + elStiffMatToAdd[ i* (*totalDofsThisElement[COL_VAR]) + colEqId] * bc_value;
	  }
	  h2Add[rowDof_elLocalI] = -sum;
	  }
	  Vector_AddTo( self->rhs->vector, *totalDofsThisElement[ROW_VAR], elementLM[ROW_VAR][0], h2Add );
	*/
}


void _StiffnessMatrix_PrintElementStiffnessMatrix(
	StiffnessMatrix* self,
	Element_LocalIndex element_lI,
	Dof_EquationNumber** rowElementLM,
	Dof_EquationNumber** colElementLM,
	double** elStiffMatToAdd )
{
	FeMesh*			rFeMesh = self->rowVariable->feMesh;
	FeMesh*			cFeMesh = self->columnVariable->feMesh;
	Dof_Index		rowDofsPerNode;
	Dof_Index		colDofsPerNode;
	Node_LocalIndex		rowNodesThisEl;
	Node_LocalIndex		colNodesThisEl;
	Node_LocalIndex		rowNode_I, colNode_I;
	Dof_Index		rowDof_I, colDof_I;
	Index			rowIndex, colIndex;
	unsigned		nRowElInc, *rowElInc;
	unsigned		nColElInc, *colElInc;

	FeMesh_GetElementNodes( rFeMesh, element_lI, self->rowInc );
	nRowElInc = IArray_GetSize( self->rowInc );
	rowElInc = IArray_GetPtr( self->rowInc );
	FeMesh_GetElementNodes( cFeMesh, element_lI, self->colInc );
	nColElInc = IArray_GetSize( self->colInc );
	colElInc = IArray_GetPtr( self->colInc );

	rowDofsPerNode = self->rowVariable->dofLayout->dofCounts[rowElInc[0]];
	colDofsPerNode = self->columnVariable->dofLayout->dofCounts[colElInc[0]];
	rowNodesThisEl = nRowElInc;
	colNodesThisEl = nColElInc;

	for ( rowNode_I=0; rowNode_I < rowNodesThisEl; rowNode_I++ ) {
		for ( rowDof_I = 0; rowDof_I < rowDofsPerNode; rowDof_I++ ) {
			for ( colNode_I=0; colNode_I < colNodesThisEl; colNode_I++ ) {
				for ( colDof_I = 0; colDof_I < colDofsPerNode; colDof_I++ ) {
					rowIndex = rowNode_I*rowDofsPerNode + rowDof_I;
					colIndex = colNode_I*colDofsPerNode + colDof_I;

					Journal_DPrintf( self->debug, "Row [%d][%d], Col [%d][%d] (LM (%4d,%4d)) = %.3f\n",
							 rowNode_I, rowDof_I,
							 colNode_I, colDof_I,
							 rowElementLM[rowNode_I][rowDof_I],
							 colElementLM[colNode_I][colDof_I],
							 elStiffMatToAdd[rowIndex][colIndex] ); 
				}			
			}
		}
	}
}

void StiffnessMatrix_AssembleElement(
	void* stiffnessMatrix,
	Element_LocalIndex element_lI,
	SystemLinearEquations* sle,
	FiniteElementContext* context,
	double** elStiffMatToAdd
	)
{
	StiffnessMatrix*        self                      = (StiffnessMatrix*) stiffnessMatrix;
	Index                   stiffnessMatrixTermCount  = Stg_ObjectList_Count( self->stiffnessMatrixTermList );
	Index                   stiffnessMatrixTerm_I;
	StiffnessMatrixTerm*    stiffnessMatrixTerm;

	for ( stiffnessMatrixTerm_I = 0 ; stiffnessMatrixTerm_I < stiffnessMatrixTermCount ; stiffnessMatrixTerm_I++ ) {
		stiffnessMatrixTerm = (StiffnessMatrixTerm*)
			Stg_ObjectList_At( self->stiffnessMatrixTermList, stiffnessMatrixTerm_I );

		StiffnessMatrixTerm_AssembleElement( stiffnessMatrixTerm, self, element_lI, sle, context, elStiffMatToAdd );
	}
}


void StiffnessMatrix_CheckElementAssembly( 
	void* stiffnessMatrix,
	Element_LocalIndex element_lI,
	double** elStiffMatToAdd,
	Index elStiffMatToAddRowSize,
	Index elStiffMatToAddColSize )
{
	StiffnessMatrix*  self = (StiffnessMatrix*)stiffnessMatrix;
	Bool              atLeastOneNonZeroElementContributionEntry = False;
	Index             elStiffMat_rowI = 0;
	Index             elStiffMat_colI = 0;
	Stream*           errorStream = Journal_Register( Error_Type, self->type );
	
	for ( elStiffMat_colI = 0; elStiffMat_colI < elStiffMatToAddColSize; elStiffMat_colI++ ) {
		for ( elStiffMat_rowI = 0; elStiffMat_rowI < elStiffMatToAddColSize; elStiffMat_rowI++ ) {
			if ( elStiffMatToAdd[elStiffMat_rowI][elStiffMat_colI] != 0.0 ) {
				atLeastOneNonZeroElementContributionEntry = True;
				break;
			}	
		}	
		if ( atLeastOneNonZeroElementContributionEntry == True ) {
			break;
		}	
	}

	Journal_Firewall( atLeastOneNonZeroElementContributionEntry == True, errorStream,
			  "Error - in %s(): while assembling matrix \"%s\", for element %u - elStiffMatToAdd assembled at this "
			  "element is all zeros."
			  "Did you register a stiffnessMatrixTerm? Is there at least one integration point in this "
			  "element?\n", __func__, self->name, element_lI  );
}


void StiffnessMatrix_AddStiffnessMatrixTerm( void* stiffnessMatrix, StiffnessMatrixTerm* stiffnessMatrixTerm ) {
	StiffnessMatrix*        self                 = (StiffnessMatrix*) stiffnessMatrix;

	stiffnessMatrixTerm = Stg_CheckType( stiffnessMatrixTerm, StiffnessMatrixTerm );
	Stg_ObjectList_Append( self->stiffnessMatrixTermList, stiffnessMatrixTerm );
}


void StiffnessMatrix_RefreshMatrix( StiffnessMatrix* self ) {
	assert( self && Stg_CheckType( self, StiffnessMatrix ) );

	/* Note: I'd like to make this a dereference, just in case there is another class still using 
	   the old matrix, but that'd require two matrices to exist at one time; i.e. lots of memory. 
	   Keeping it as a free just means that other classes need to assume they never own the matrix. */
	Matrix_SetComm( self->matrix, self->comm );
	Matrix_SetLocalSize( self->matrix, self->rowLocalSize, self->colLocalSize );

#ifdef HAVE_PETSC
	if( Stg_CheckType( self->matrix, PETScMatrix ) ) {
		PETScMatrix_SetNonZeroStructure( (PETScMatrix*)self->matrix, self->nonZeroCount, 
						 self->diagonalNonZeroIndices, self->offDiagonalNonZeroIndices );
	}
#endif
}


#if 0
void StiffnessMatrix_CalcNonZeros( void* stiffnessMatrix ) {
	StiffnessMatrix*	self = (StiffnessMatrix*)stiffnessMatrix;
	Stream*			stream;
	FeVariable		*rowVar, *colVar;
	FeMesh			*rowMesh, *colMesh;
	FeEquationNumber	*rowEqNum, *colEqNum;
	DofLayout		*rowDofs, *colDofs;
	unsigned		nRowEqs, nColEqs;
	unsigned		nRowNodes, *rowNodes;
	unsigned		nColNodes, *colNodes;
	unsigned		*nMaxNonZeros, maxNZ;
	unsigned		*nNonZeros, **nonZeros;
	unsigned		*nDiagNZs, *nOffDiagNZs;
	unsigned		rowEq, colEq;
	unsigned		netNonZeros;
	unsigned		e_i, nz_i, eq_i;
	unsigned		n_i, dof_i;
	unsigned		n_j, dof_j;

	assert( self && Stg_CheckType( self, StiffnessMatrix ) );
	assert( self->rowVariable );

	stream = Journal_Register( Info_Type, self->type );
	Journal_Printf( stream, "Stiffness matrix: '%s'\n", self->name );
	Stream_Indent( stream );
	Journal_Printf( stream, "Calculating number of nonzero entries...\n" );
	Stream_Indent( stream );

	rowVar = self->rowVariable;
	colVar = self->columnVariable ? self->columnVariable : rowVar;
	rowMesh = rowVar->feMesh;
	colMesh = colVar->feMesh;
	rowEqNum = rowVar->eqNum;
	colEqNum = colVar->eqNum;
	nRowEqs = rowEqNum->localEqNumsOwnedCount;
	nColEqs = colEqNum->localEqNumsOwnedCount;
	rowDofs = rowVar->dofLayout;
	colDofs = colVar->dofLayout;

	nMaxNonZeros = AllocArray( unsigned, nRowEqs );
	memset( nMaxNonZeros, 0, nRowEqs * sizeof(unsigned) );

	for( e_i = 0; e_i < FeMesh_GetElementDomainSize( rowMesh ); e_i++ ) {
		FeMesh_GetElementNodes( rowMesh, e_i, &nRowNodes, &rowNodes );
		FeMesh_GetElementNodes( colMesh, e_i, &nColNodes, &colNodes );

		for( n_i = 0; n_i < nRowNodes; n_i++ ) {
			if( rowNodes[n_i] >= FeMesh_GetNodeLocalSize( rowMesh ) )
				continue;

			for( dof_i = 0; dof_i < rowDofs->dofCounts[rowNodes[n_i]]; dof_i++ ) {
				rowEq = rowEqNum->locationMatrix[e_i][n_i][dof_i];
				if( rowEq == (unsigned)-1 )
					continue;

				for( n_j = 0; n_j < nColNodes; n_j++ ) {
					for( dof_j = 0; dof_j < colDofs->dofCounts[colNodes[n_j]]; dof_j++ ) {
						colEq = colEqNum->locationMatrix[e_i][n_j][dof_j];
						if( colEq == (unsigned)-1 )
							continue;

						nMaxNonZeros[rowEq - rowEqNum->firstOwnedEqNum]++;
					}
				}
			}
		}
	}

	maxNZ = 0;
	for( eq_i = 0; eq_i < nRowEqs; eq_i++ ) {
		if( nMaxNonZeros[eq_i] > maxNZ )
			maxNZ = nMaxNonZeros[eq_i];
	}

	FreeArray( nMaxNonZeros );

	nNonZeros = AllocArray( unsigned, nRowEqs );
	memset( nNonZeros, 0, nRowEqs * sizeof(unsigned) );
	nonZeros = AllocArray2D( unsigned, nRowEqs, maxNZ );

	for( e_i = 0; e_i < FeMesh_GetElementDomainSize( rowMesh ); e_i++ ) {
		FeMesh_GetElementNodes( rowMesh, e_i, &nRowNodes, &rowNodes );
		FeMesh_GetElementNodes( colMesh, e_i, &nColNodes, &colNodes );

		for( n_i = 0; n_i < nRowNodes; n_i++ ) {
			if( rowNodes[n_i] >= FeMesh_GetNodeLocalSize( rowMesh ) )
				continue;

			for( dof_i = 0; dof_i < rowDofs->dofCounts[rowNodes[n_i]]; dof_i++ ) {
				rowEq = rowEqNum->locationMatrix[e_i][n_i][dof_i];
				if( rowEq == (unsigned)-1 )
					continue;

				for( n_j = 0; n_j < nColNodes; n_j++ ) {
					for( dof_j = 0; dof_j < colDofs->dofCounts[colNodes[n_j]]; dof_j++ ) {
						colEq = colEqNum->locationMatrix[e_i][n_j][dof_j];
						if( colEq == (unsigned)-1 )
							continue;

						StiffnessMatrix_TrackUniqueEqs( self, rowEq - rowEqNum->firstOwnedEqNum, colEq, 
										nNonZeros, nonZeros );
					}
				}
			}
		}
	}

	nDiagNZs = AllocArray( unsigned, nRowEqs );
	memset( nDiagNZs, 0, nRowEqs * sizeof(unsigned) );
	nOffDiagNZs = AllocArray( unsigned, nRowEqs );
	memset( nOffDiagNZs, 0, nRowEqs * sizeof(unsigned) );

	netNonZeros = 0;
	for( eq_i = 0; eq_i < nRowEqs; eq_i++ ) {
		netNonZeros += nNonZeros[eq_i];

		for( nz_i = 0; nz_i < nNonZeros[eq_i]; nz_i++ ) {
			if( nonZeros[eq_i][nz_i] < nRowEqs )
				nDiagNZs[eq_i]++;
			else
				nOffDiagNZs[eq_i]++;
		}
	}

	FreeArray( nNonZeros );
	FreeArray( nonZeros );

	self->diagonalNonZeroIndices = nDiagNZs;
	self->offDiagonalNonZeroIndices = nOffDiagNZs;

	Journal_Printf( stream, "Found %d nonzero entries.\n", netNonZeros );
	Journal_Printf( stream, "Done.\n" );
	Stream_UnIndent( stream );
	Stream_UnIndent( stream );
}


void StiffnessMatrix_TrackUniqueEqs( StiffnessMatrix* self, unsigned rowEq, unsigned colEq, 
				     unsigned* nNonZeros, unsigned** nonZeros )
{
	unsigned	nz_i;

	assert( self );
	assert( nNonZeros && nonZeros );

	for( nz_i = 0; nz_i < nNonZeros[rowEq]; nz_i++ ) {
		if( nonZeros[rowEq][nz_i] == colEq )
			break;
	}
	if( nz_i < nNonZeros[rowEq] )
		return;

	nNonZeros[rowEq]++;
	nonZeros[rowEq][nNonZeros[rowEq] - 1] = colEq;
}
#endif


Bool StiffnessMatrix_ZeroBCsAsm_RowR( void* stiffMat, Assembler* assm ) {
	memset( ((StiffnessMatrix*)stiffMat)->elStiffMat[assm->rowInd], 0, 
		((StiffnessMatrix*)stiffMat)->nColDofs * sizeof(double) );
	return False;
}

Bool StiffnessMatrix_ZeroBCsAsm_ColR( void* stiffMat, Assembler* assm ) {
	((StiffnessMatrix*)stiffMat)->elStiffMat[assm->rowInd][assm->colInd] = 0.0;
	return True;
}

Bool StiffnessMatrix_BCAsm_ColR( void* stiffMat, Assembler* assm ) {
	double	bc;

	bc = DofLayout_GetValueDouble( assm->colVar->dofLayout, 
				       assm->colNodeInd, 
				       assm->colDofInd );
	((StiffnessMatrix*)stiffMat)->bcVals[assm->rowInd] -= bc * 
		((StiffnessMatrix*)stiffMat)->elStiffMat[assm->rowInd][assm->colInd];
	return True;
}

Bool StiffnessMatrix_TransBCAsm_ColR( void* stiffMat, Assembler* assm ) {
	double	bc;

	bc = DofLayout_GetValueDouble( assm->colVar->dofLayout, 
				       assm->colNodeInd, 
				       assm->colDofInd );
	((StiffnessMatrix*)stiffMat)->bcVals[assm->rowInd] -= bc * 
		((StiffnessMatrix*)stiffMat)->elStiffMat[assm->colInd][assm->rowInd];
	return True;
}

Bool StiffnessMatrix_DiagBCsAsm_RowR( void* stiffMat, Assembler* assm ) {
	static const double	one = 1.0;

	Matrix_AddEntries( ((StiffnessMatrix*)stiffMat)->matrix, 
			   1, &assm->rowEq, 
			   1, &assm->rowEq, 
			   (double*)&one );
	return True;
}
