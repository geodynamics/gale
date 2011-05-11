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
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"
#include "types.h"
#include "Stokes_SLE_PenaltySolver.h"

#include <assert.h>
#include <string.h>

#include "Stokes_SLE.h"

const Type Stokes_SLE_PenaltySolver_Type = "Stokes_SLE_PenaltySolver";

void* Stokes_SLE_PenaltySolver_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(Stokes_SLE_PenaltySolver);
	Type                                                      type = Stokes_SLE_PenaltySolver_Type;
	Stg_Class_DeleteFunction*                              _delete = _Stokes_SLE_PenaltySolver_Delete;
	Stg_Class_PrintFunction*                                _print = _Stokes_SLE_PenaltySolver_Print;
	Stg_Class_CopyFunction*                                  _copy = _Stokes_SLE_PenaltySolver_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = Stokes_SLE_PenaltySolver_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _Stokes_SLE_PenaltySolver_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Stokes_SLE_PenaltySolver_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _SLE_Solver_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _SLE_Solver_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _SLE_Solver_Destroy;
	SLE_Solver_SolverSetupFunction*                   _solverSetup = _Stokes_SLE_PenaltySolver_SolverSetup;
	SLE_Solver_SolveFunction*                               _solve = _Stokes_SLE_PenaltySolver_Solve;
	SLE_Solver_GetResidualFunc*                       _getResidual = _Stokes_SLE_PenaltySolver_GetResidual;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return _Stokes_SLE_PenaltySolver_New(  STOKES_SLE_PENALTYSOLVER_PASSARGS  );
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
Stokes_SLE_PenaltySolver* _Stokes_SLE_PenaltySolver_New(  STOKES_SLE_PENALTYSOLVER_DEFARGS  )
{
	Stokes_SLE_PenaltySolver* self;

	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(Stokes_SLE_PenaltySolver) );
	self = (Stokes_SLE_PenaltySolver*) _SLE_Solver_New(  SLE_SOLVER_PASSARGS  );

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


void* _Stokes_SLE_PenaltySolver_Copy( const void* stokesSlePenaltySolver, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	Stokes_SLE_PenaltySolver* self = (Stokes_SLE_PenaltySolver*)stokesSlePenaltySolver;
	Stokes_SLE_PenaltySolver*	newStokesSlePenaltySolver;
	
	newStokesSlePenaltySolver = (Stokes_SLE_PenaltySolver*)_SLE_Solver_Copy( self, dest, deep, nameExt, ptrMap );
	
	return (void*) newStokesSlePenaltySolver;
}


void _Stokes_SLE_PenaltySolver_Build( void* solver, void* stokesSLE ) {
	Stokes_SLE_PenaltySolver*	self  = (Stokes_SLE_PenaltySolver*)solver;

 	Journal_DPrintf( self->debug, "In %s\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );
	Stream_UnIndentBranch( StgFEM_Debug );
}

void _Stokes_SLE_PenaltySolver_AssignFromXML( void* solver, Stg_ComponentFactory* cf, void* data ) {
	Stokes_SLE_PenaltySolver* self         = (Stokes_SLE_PenaltySolver*) solver;

	_SLE_Solver_AssignFromXML( self, cf, data );
	
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
	Mat                     kMatrix         = sle->kStiffMat->matrix;
	Mat                     gradMat         = sle->gStiffMat->matrix;
	Mat                     divMat          = NULL;
	Mat                     C_Mat           = sle->cStiffMat->matrix;
	Vec                     uVec            = sle->uSolnVec->vector;
	Vec                     pVec            = sle->pSolnVec->vector;
	Vec                     fVec            = sle->fForceVec->vector;
	Vec                     hVec            = sle->hForceVec->vector;
	Vec     		hTempVec;
	Vec    			fTempVec;
	Mat    			GTrans;
	KSP			sles_v;
	double	 		negOne=-1.0;
	double	 		one=1.0;
	Mat    			kHat;
	Mat    			C_InvMat;
	Vec    			diagC;
	PC			pc;

	Journal_DPrintf( self->debug, "In %s():\n", __func__ );

	VecDuplicate( hVec, &hTempVec );
	VecDuplicate( fVec, &fTempVec );
	VecDuplicate( pVec, &diagC );
	
	if( sle->dStiffMat == NULL ) {
		Journal_DPrintf( self->debug, "Div matrix == NULL : Problem is assumed to be symmetric. ie Div = GTrans \n");
#if( PETSC_VERSION_MAJOR <= 2 )
		MatTranspose( gradMat, &GTrans );
#else
		MatTranspose( gradMat, MAT_INITIAL_MATRIX, &GTrans );
#endif
		divMat = GTrans;
	}
	else {
		/* make a copy we can play with */
		MatCreate( sle->comm, &GTrans );
		MatCopy( sle->dStiffMat->matrix, GTrans, DIFFERENT_NONZERO_PATTERN );
		divMat = GTrans;
	}
	
	/* Create CInv */
	MatGetDiagonal( C_Mat, diagC );
	VecReciprocal( diagC );
	MatDiagonalSet( C_Mat, diagC, INSERT_VALUES );
	C_InvMat = C_Mat;				/* Use pointer CInv since C has been inverted */
	
	/* Build RHS : rhs = f - GCInv h */
	MatMult( C_InvMat, hVec, hTempVec );
	VecScale( hTempVec, negOne );
	MatMultAdd( gradMat, hTempVec, fVec, fTempVec );
	
	/* Build G CInv GTrans */
/* 	MatTranspose( gradMat, &GTrans ); */
/* 	 since CInv is diagonal we can just scale mat entries by the diag vector */
	MatDiagonalScale( divMat, diagC, PETSC_NULL );  /*  Div = CInve Div */
        /* 	MatMatMult_any( &tmpMat, C, *divMat );	*/ /* tmpMat = CInv Div */
	
	
	Journal_DPrintf( self->debug, "UpdivMat mat mat mult \n");
	MatPtAP( gradMat, divMat, MAT_INITIAL_MATRIX, 1.0, &kHat );
	Journal_DPrintf( self->debug, "done mult \n");
	MatScale( kHat, -1 );
	MatAXPY( kHat, one, kMatrix, SAME_NONZERO_PATTERN );
	
	/* Setup solver context and make sure that it uses a direct solver */
	KSPCreate( sle->comm, &sles_v );
	KSPSetOperators( sles_v, kHat, kHat, DIFFERENT_NONZERO_PATTERN );
	KSPSetType( sles_v, KSPPREONLY );
	KSPGetPC( sles_v, &pc );
	PCSetType( pc, PCLU );

	KSPSolve( sles_v, fTempVec, uVec );
	
	/* Recover p */
	if( sle->dStiffMat == NULL ) {
/* 		 since Div was modified when C is diagonal, re build the transpose */
#if( PETSC_VERSION_MAJOR <= 2 )
		MatTranspose( gradMat, &GTrans );
#else
		MatTranspose( gradMat, MAT_INITIAL_MATRIX, &GTrans );
#endif
		divMat = GTrans;
	}
	else {
/* 		 never modified Div_null so set divMat to point back to it */
		divMat = sle->dStiffMat->matrix;
	}
	MatMult( divMat, uVec, hTempVec );    /* hTemp = Div v */
	VecAXPY( hVec, negOne, hTempVec );    /* hTemp = H - hTemp   : hTemp = H - Div v */
	MatMult( C_InvMat, hTempVec, pVec );  /* p = CInv hTemp      : p = CInv ( H - Div v ) */
	
	if( kHat != PETSC_NULL )     MatDestroy( kHat );
	if( fTempVec != PETSC_NULL ) VecDestroy( fTempVec );
	if( hTempVec != PETSC_NULL ) VecDestroy( hTempVec );
	if( diagC != PETSC_NULL )    VecDestroy( diagC );
	if( sles_v != PETSC_NULL )   KSPDestroy( sles_v );
	if( GTrans != PETSC_NULL )   MatDestroy( GTrans );
}


Vec _Stokes_SLE_PenaltySolver_GetResidual( void* solver, Index fv_I ) {
/* 	 TODO */
	return NULL;
}



