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
** $Id: Stokes_SLE_PenaltySolver.c 964 2007-10-11 08:03:06Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/LinearAlgebra/LinearAlgebra.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"
#include "types.h"
#include "Stokes_SLE_PenaltySolver.h"

#include <assert.h>
#include <string.h>

#include "Stokes_SLE.h"

const Type Stokes_SLE_PenaltySolver_Type = "Stokes_SLE_PenaltySolver";

void* Stokes_SLE_PenaltySolver_DefaultNew( Name name ) {
	return _Stokes_SLE_PenaltySolver_New( 
		sizeof(Stokes_SLE_PenaltySolver), 
		Stokes_SLE_PenaltySolver_Type, 
		_Stokes_SLE_PenaltySolver_Delete, 
		_Stokes_SLE_PenaltySolver_Print, 
		_Stokes_SLE_PenaltySolver_Copy,
		Stokes_SLE_PenaltySolver_DefaultNew,
		_Stokes_SLE_PenaltySolver_Construct,
		_Stokes_SLE_PenaltySolver_Build,
		_SLE_Solver_Initialise,
		_SLE_Solver_Execute,
		_SLE_Solver_Destroy,
		_Stokes_SLE_PenaltySolver_SolverSetup, 
		_Stokes_SLE_PenaltySolver_Solve, 
		_Stokes_SLE_PenaltySolver_GetResidual,
		name );
}

Stokes_SLE_PenaltySolver* Stokes_SLE_PenaltySolver_New(
		Name                                        name,
		Bool                                        useStatSolve, 
		int                                         statReps )
{
	Stokes_SLE_PenaltySolver* self = (Stokes_SLE_PenaltySolver*) Stokes_SLE_PenaltySolver_DefaultNew( name );

	Stokes_SLE_PenaltySolver_InitAll( self, useStatSolve, statReps );

	return self;
}


/* Creation implementation / Virtual constructor */
Stokes_SLE_PenaltySolver* _Stokes_SLE_PenaltySolver_New( 
		SizeT                                       sizeOfSelf,
		Type                                        type,
		Stg_Class_DeleteFunction*                   _delete,
		Stg_Class_PrintFunction*                    _print,
		Stg_Class_CopyFunction*                     _copy, 
		Stg_Component_DefaultConstructorFunction*   _defaultConstructor,
		Stg_Component_ConstructFunction*            _construct,
		Stg_Component_BuildFunction*                _build,
		Stg_Component_InitialiseFunction*           _initialise,
		Stg_Component_ExecuteFunction*              _execute,
		Stg_Component_DestroyFunction*              _destroy,
		SLE_Solver_SolverSetupFunction*             _solverSetup,
		SLE_Solver_SolveFunction*                   _solve,
		SLE_Solver_GetResidualFunc*                 _getResidual, 
		Name                                        name )
{
	Stokes_SLE_PenaltySolver* self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(Stokes_SLE_PenaltySolver) );
	self = (Stokes_SLE_PenaltySolver*) _SLE_Solver_New( 
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
		_solverSetup,
		_solve,
		_getResidual, 
		name );

	/* Virtual info */
	return self;
}


void _Stokes_SLE_PenaltySolver_Init( void* solver ) {
	Stokes_SLE_PenaltySolver* self = (Stokes_SLE_PenaltySolver*)solver;

	self->isConstructed       = True;
}

void Stokes_SLE_PenaltySolver_InitAll( 
		void*                        solver,
		Bool                         useStatSolve,
		int                          statReps )
{
	Stokes_SLE_PenaltySolver* self = (Stokes_SLE_PenaltySolver*)solver;

	SLE_Solver_InitAll( self, useStatSolve, statReps );
	_Stokes_SLE_PenaltySolver_Init( self );
}


void _Stokes_SLE_PenaltySolver_Delete( void* solver ) {
	Stokes_SLE_PenaltySolver* self = (Stokes_SLE_PenaltySolver*)solver;
		
	Journal_DPrintf( self->debug, "In: %s \n", __func__);

	Stream_IndentBranch( StgFEM_Debug );

	Stream_UnIndentBranch( StgFEM_Debug );
}       


void _Stokes_SLE_PenaltySolver_Print( void* solver, Stream* stream ) {
	Stokes_SLE_PenaltySolver* self = (Stokes_SLE_PenaltySolver*)solver;

	_SLE_Solver_Print( self, stream );
}


void* _Stokes_SLE_PenaltySolver_Copy( void* stokesSlePenaltySolver, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	Stokes_SLE_PenaltySolver* self = (Stokes_SLE_PenaltySolver*)stokesSlePenaltySolver;
	Stokes_SLE_PenaltySolver*	newStokesSlePenaltySolver;
	
	newStokesSlePenaltySolver = _SLE_Solver_Copy( self, dest, deep, nameExt, ptrMap );
	
	return (void*) newStokesSlePenaltySolver;
}


void _Stokes_SLE_PenaltySolver_Build( void* solver, void* stokesSLE ) {
	Stokes_SLE_PenaltySolver*	self  = (Stokes_SLE_PenaltySolver*)solver;

 	Journal_DPrintf( self->debug, "In %s\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );
	Stream_UnIndentBranch( StgFEM_Debug );
}

void _Stokes_SLE_PenaltySolver_Construct( void* solver, Stg_ComponentFactory* cf, void* data ) {
	Stokes_SLE_PenaltySolver* self         = (Stokes_SLE_PenaltySolver*) solver;

	_SLE_Solver_Construct( self, cf, data );
	
	_Stokes_SLE_PenaltySolver_Init( self );
}

void _Stokes_SLE_PenaltySolver_Execute( void* solver, void* data ) {
}

void _Stokes_SLE_PenaltySolver_Destroy( void* solver, void* data ) {
}

void _Stokes_SLE_PenaltySolver_Initialise( void* solver, void* stokesSLE ) {
	Stokes_SLE_PenaltySolver* self = (Stokes_SLE_PenaltySolver*) solver;
	Stokes_SLE*             sle  = (Stokes_SLE*)             stokesSLE;
	
	/* Initialise Parent */
	_SLE_Solver_Initialise( self, sle );
}

/* SolverSetup */
void _Stokes_SLE_PenaltySolver_SolverSetup( void* solver, void* stokesSLE ) {
	Stokes_SLE_PenaltySolver* self = (Stokes_SLE_PenaltySolver*) solver;
	
 	Journal_DPrintf( self->debug, "In %s:\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );
	Stream_UnIndentBranch( StgFEM_Debug );
}


void _Stokes_SLE_PenaltySolver_Solve( void* solver,void* stokesSLE ) {
	Stokes_SLE_PenaltySolver* self            = (Stokes_SLE_PenaltySolver*)solver;	
	Stokes_SLE*             sle             = (Stokes_SLE*)stokesSLE;
	/* Create shortcuts to stuff needed on sle */
	Matrix*                 kMatrix         = sle->kStiffMat->matrix;
	Matrix*                 gradMat         = sle->gStiffMat->matrix;
	Matrix*                 divMat          = NULL;
	Matrix*                 C_Mat           = sle->cStiffMat->matrix;
	Vector*                 uVec            = sle->uSolnVec->vector;
	Vector*                 pVec            = sle->pSolnVec->vector;
	Vector*                 fVec            = sle->fForceVec->vector;
	Vector*                 hVec            = sle->hForceVec->vector;
	Vector* 		hTempVec;
	Vector*			fTempVec;

	Matrix*			GTrans;
	MatrixSolver*		sles_v;
	double	 		negOne=-1.0;
	double	 		one=1.0;
	Matrix*			kHat;
	Matrix*			C_InvMat;
	Vector*			diagC;
	
	Journal_DPrintf( self->debug, "In %s():\n", __func__ );
	
	Vector_Duplicate( hVec, (void**)&hTempVec );
	Vector_SetLocalSize( hTempVec, Vector_GetLocalSize( hVec ) );
	Vector_Duplicate( fVec, (void**)&fTempVec );
	Vector_SetLocalSize( fTempVec, Vector_GetLocalSize( fVec ) );
	Vector_Duplicate( pVec, (void**)&diagC );
	Vector_SetLocalSize( diagC, Vector_GetLocalSize( pVec ) );
	
	if( sle->dStiffMat == NULL ) {
		Journal_DPrintf( self->debug, "Div matrix == NULL : Problem is assumed to be symmetric. ie Div = GTrans \n");
		Matrix_Duplicate( gradMat, (void**)&GTrans );
		Matrix_Transpose( gradMat, GTrans );
		divMat = GTrans;
	}
	else {
		/* make a copy we can play with */
		Matrix_Duplicate( sle->dStiffMat->matrix, (void**)&GTrans );
		Matrix_CopyEntries( sle->dStiffMat->matrix, GTrans );
		divMat = GTrans;
	}
	
	/* Create CInv */
	Matrix_GetDiagonal( C_Mat, diagC );
	Vector_Reciprocal( diagC );
	Matrix_DiagonalInsertEntries( C_Mat, diagC );
	C_InvMat = C_Mat;				/* Use pointer CInv since C has been inverted */
	
	/* Build RHS : rhs = f - GCInv h */
	Matrix_Multiply( C_InvMat, hVec, hTempVec );		/* hTemp = CInv h	*/
	Vector_Scale( hTempVec, negOne );		/* hTemp = -hTemp	: -CInv h	*/
	Matrix_MultiplyAdd( gradMat, hTempVec, fVec, fTempVec );	/* fTemp = F + G hTemp	: fTemp = F - G CInv h */
	
	/* Build G CInv GTrans */
/* 	MatTranspose( gradMat, &GTrans ); */
/* 	 since CInv is diagonal we can just scale mat entries by the diag vector */
	Matrix_DiagonalScale( divMat, diagC, NULL );	/*  Div = CInve Div */
        /* 	MatMatMult_any( &tmpMat, C, *divMat );	*/ /* tmpMat = CInv Div */
	
	
	Journal_DPrintf( self->debug, "UpdivMat mat mat mult \n");
	Matrix_Duplicate( gradMat, (void**)&kHat );
	Matrix_MatrixMultiply( gradMat, divMat, kHat );		/* tmpMat2 = G CInv Div */
	Journal_DPrintf( self->debug, "done mult \n");
	Matrix_Scale( kHat, -1 );			/* tmpMat2 = - G CInv Div */
	Matrix_AddScaled( kHat, one, kMatrix );	/* kHat = kHat + kMatrix */
	
	/* Setup solver context and make sure that it uses a direct solver */
#ifndef HAVE_PETSC
#error Need PETSc!
#endif
	sles_v = (MatrixSolver*)PETScMatrixSolver_New( "" );
	MatrixSolver_SetComm( sles_v, sle->comm );
	MatrixSolver_SetMatrix( sles_v, kHat );
	PETScMatrixSolver_SetKSPType( sles_v, PETScMatrixSolver_KSPType_PreOnly );
	PETScMatrixSolver_SetPCType( sles_v, PETScMatrixSolver_PCType_LU );

	MatrixSolver_Solve( sles_v, fTempVec, uVec );
	
	/* Recover p */
	if( sle->dStiffMat == NULL ) {
/* 		 since Div was modified when C is diagonal, re build the transpose */
		Matrix_Transpose( gradMat, GTrans );
		divMat = GTrans;
	}
	else {
/* 		 never modified Div_null so set divMat to point back to it */
		divMat = sle->dStiffMat->matrix;
	}
	Matrix_Multiply( divMat, uVec, hTempVec );		/* hTemp = Div v */
	Vector_AddScaled( hVec, negOne, hTempVec );	/* hTemp = H - hTemp	: hTemp = H - Div v */
	Matrix_Multiply( C_InvMat, hTempVec, pVec );		/* p = CInv hTemp	: p = CInv ( H - Div v ) */
	
	
	FreeObject( kHat );
	FreeObject( fTempVec );
	FreeObject( hTempVec );
	FreeObject( diagC );
	FreeObject( sles_v );
	FreeObject( GTrans );
}


Vector* _Stokes_SLE_PenaltySolver_GetResidual( void* solver, Index fv_I ) {
/* 	 TODO */
	return NULL;
}
