/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
** $Id: PETScVector.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>
#include "StGermain/StGermain.h"
#include "Discretisation/Discretisation.h"
#include "LinearAlgebra.h"


/* Textual name of this class */
const Type PETScVector_Type = "PETScVector";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

PETScVector* PETScVector_New( Name name ) {
	return _PETScVector_New( sizeof(PETScVector), 
				 PETScVector_Type, 
				 _PETScVector_Delete, 
				 _PETScVector_Print, 
				 NULL, 
				 (void* (*)(Name))PETScVector_New, 
				 _PETScVector_Construct, 
				 _PETScVector_Build, 
				 _PETScVector_Initialise, 
				 _PETScVector_Execute, 
				 _PETScVector_Destroy, 
				 name, 
				 NON_GLOBAL, 
				 PETScVector_SetComm, 
				 PETScVector_SetGlobalSize, 
				 PETScVector_SetLocalSize, 
				 PETScVector_AddEntries, 
				 PETScVector_InsertEntries, 
				 PETScVector_SetScalar, 
				 PETScVector_Zero, 
				 PETScVector_AssemblyBegin, 
				 PETScVector_AssemblyEnd, 
				 PETScVector_Add, 
				 PETScVector_AddScaled, 
				 PETScVector_ScaleAdd, 
				 PETScVector_Subtract, 
				 PETScVector_Scale, 
				 PETScVector_DotProduct, 
				 PETScVector_L2Norm, 
				 PETScVector_LInfNorm, 
				 PETScVector_PointwiseMultiply, 
				 PETScVector_PointwiseDivide, 
				 PETScVector_Reciprocal, 
				 PETScVector_GetGlobalSize, 
				 PETScVector_GetLocalSize, 
				 PETScVector_GetArray, 
				 PETScVector_RestoreArray, 
				 PETScVector_Duplicate, 
				 PETScVector_CopyEntries, 
				 _Vector_View );
}

PETScVector* _PETScVector_New( PETSCVECTOR_DEFARGS ) {
	PETScVector*	self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(PETScVector) );
	self = (PETScVector*)_Vector_New( VECTOR_PASSARGS );

	/* Virtual info */

	/* PETScVector info */
	_PETScVector_Init( self );

	return self;
}

void _PETScVector_Init( PETScVector* self ) {
	assert( self && Stg_CheckType( self, PETScVector ) );

	self->petscVec = PETSC_NULL;
	PETScVector_SetComm( self, MPI_COMM_WORLD );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _PETScVector_Delete( void* vector ) {
	PETScVector*	self = (PETScVector*)vector;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScVector ) );

	ec = VecDestroy( self->petscVec );
	CheckPETScError( ec );

	/* Delete the parent. */
	_Vector_Delete( self );
}

void _PETScVector_Print( void* vector, Stream* stream ) {
	PETScVector*	self = (PETScVector*)vector;
	
	/* Set the Journal for printing informations */
	Stream* vectorStream;
	vectorStream = Journal_Register( InfoStream_Type, "PETScVectorStream" );

	assert( self && Stg_CheckType( self, PETScVector ) );

	/* Print parent */
	Journal_Printf( stream, "PETScVector (ptr): (%p)\n", self );
	_Vector_Print( self, stream );
}

void _PETScVector_Construct( void* vector, Stg_ComponentFactory* cf, void* data ) {
	PETScVector*	self = (PETScVector*)vector;

	assert( self && Stg_CheckType( self, PETScVector ) );
	assert( cf );
}

void _PETScVector_Build( void* vector, void* data ) {
}

void _PETScVector_Initialise( void* vector, void* data ) {
}

void _PETScVector_Execute( void* vector, void* data ) {
}

void _PETScVector_Destroy( void* vector, void* data ) {
}

void PETScVector_SetComm( void* vector, MPI_Comm comm ) {
	PETScVector*	self = (PETScVector*)vector;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScVector ) );

	_Vector_SetComm( self, comm );

	if( self->petscVec != PETSC_NULL )
		VecDestroy( self->petscVec );
	ec = VecCreate( comm, &self->petscVec );
	CheckPETScError( ec );

}

void PETScVector_SetGlobalSize( void* vector, unsigned size ) {
	PETScVector*	self = (PETScVector*)vector;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScVector ) );

	ec = VecSetSizes( self->petscVec, PETSC_DECIDE, (PetscInt)size );
	CheckPETScError( ec );
	ec = VecSetFromOptions( self->petscVec );
	CheckPETScError( ec );

        #if( (PETSC_VERSION_MINOR == 3) && (PETSC_VERSION_SUBMINOR > 2) )
        VecSetOption( self->petscVec, VEC_IGNORE_NEGATIVE_INDICES );
        #endif

}

void PETScVector_SetLocalSize( void* vector, unsigned size ) {
	PETScVector*	self = (PETScVector*)vector;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScVector ) );

	ec = VecSetSizes( self->petscVec, (PetscInt)size, PETSC_DECIDE );
	CheckPETScError( ec );
	ec = VecSetFromOptions( self->petscVec );
	CheckPETScError( ec );

        #if( (PETSC_VERSION_MINOR == 3) && (PETSC_VERSION_SUBMINOR > 2) )
        VecSetOption( self->petscVec, VEC_IGNORE_NEGATIVE_INDICES );
        #endif

}

void PETScVector_AddEntries( void* vector, unsigned nEntries, unsigned* indices, double* values ) {
	PETScVector*	self = (PETScVector*)vector;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScVector ) );
	assert( !nEntries || (indices && values) );

	ec = VecSetValues( self->petscVec, nEntries, (PetscInt*)indices, (PetscScalar*)values, ADD_VALUES );
	CheckPETScError( ec );
}

void PETScVector_InsertEntries( void* vector, unsigned nEntries, unsigned* indices, double* values ) {
	PETScVector*	self = (PETScVector*)vector;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScVector ) );
	assert( !nEntries || (indices && values) );

	ec = VecSetValues( self->petscVec, nEntries, (PetscInt*)indices, (PetscScalar*)values, INSERT_VALUES );
	CheckPETScError( ec );
}

void PETScVector_SetScalar( void* vector, double scalar ) {
	PETScVector*	self = (PETScVector*)vector;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScVector ) );

	ec = VecSet( self->petscVec, scalar );
	CheckPETScError( ec );
}

void PETScVector_Zero( void* vector ) {
	PETScVector_SetScalar( vector, 0.0 );
}

void PETScVector_AssemblyBegin( void* vector ) {
	PETScVector*	self = (PETScVector*)vector;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScVector ) );

	ec = VecAssemblyBegin( self->petscVec );
	CheckPETScError( ec );
}

void PETScVector_AssemblyEnd( void* vector ) {
	PETScVector*	self = (PETScVector*)vector;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScVector ) );

	ec = VecAssemblyEnd( self->petscVec );
	CheckPETScError( ec );
}

void PETScVector_Add( void* vector, void* vector0 ) {
	PETScVector_AddScaled( vector, 1.0, vector0 );
}

void PETScVector_AddScaled( void* vector, double scalar, void* _vector0 ) {
	PETScVector*	self = (PETScVector*)vector;
	PETScVector*	vector0 = (PETScVector*)_vector0;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScVector ) );
	assert( vector0 && Stg_CheckType( vector0, PETScVector ) );

	ec = VecAXPY( self->petscVec, (PetscScalar)scalar, vector0->petscVec );
	CheckPETScError( ec );
}

void PETScVector_ScaleAdd( void* vector, double scalar, void* _vector0 ) {
	PETScVector*	self = (PETScVector*)vector;
	PETScVector*	vector0 = (PETScVector*)_vector0;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScVector ) );
	assert( vector0 && Stg_CheckType( vector0, PETScVector ) );

	ec = VecAYPX( self->petscVec, (PetscScalar)scalar, vector0->petscVec );
	CheckPETScError( ec );
}

void PETScVector_Subtract( void* vector, void* vector0 ) {
	PETScVector_AddScaled( vector, -1.0, vector0 );
}

void PETScVector_Scale( void* vector, double factor ) {
	PETScVector*	self = (PETScVector*)vector;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScVector ) );

	ec = VecScale( self->petscVec, (PetscScalar)factor );
	CheckPETScError(ec);
}

double PETScVector_DotProduct( void* vector, void* _vector0 ) {
	PETScVector*	self = (PETScVector*)vector;
	PETScVector*	vector0 = (PETScVector*)_vector0;
	PetscErrorCode	ec;
	double		dot;

	assert( self && Stg_CheckType( self, PETScVector ) );
	assert( vector0 && Stg_CheckType( vector0, PETScVector ) );

	ec = VecDot( self->petscVec, vector0->petscVec, (PetscScalar*)&dot );
	CheckPETScError(ec);

	return dot;
}

double PETScVector_L2Norm( void* vector ) {
	PETScVector*	self = (PETScVector*)vector;
	PetscErrorCode	ec;
	double         norm;

	ec = VecNorm( self->petscVec, NORM_2, (PetscScalar*)&norm );
	CheckPETScError(ec);

	return norm;
}

double PETScVector_LInfNorm( void* vector ) {
	PETScVector*	self = (PETScVector*)vector;
	PetscErrorCode	ec;
	double         norm;

	ec = VecNorm( self->petscVec, NORM_INFINITY, (PetscScalar*)&norm );
	CheckPETScError(ec);

	return norm;
}

void PETScVector_PointwiseMultiply( void* vector, void* _enumVec, void* _denomVec ) {
	PETScVector*	self = (PETScVector*)vector;
	PETScVector*	enumVec = (PETScVector*)_enumVec;
	PETScVector*	denomVec = (PETScVector*)_denomVec;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScVector ) );
	assert( enumVec && Stg_CheckType( enumVec, PETScVector ) );
	assert( denomVec && Stg_CheckType( denomVec, PETScVector ) );

	ec = VecPointwiseMult( self->petscVec, enumVec->petscVec, denomVec->petscVec ); 
	CheckPETScError(ec);
}

void PETScVector_PointwiseDivide( void* vector, void* _enumVec, void* _denomVec ) {
	PETScVector*	self = (PETScVector*)vector;
	PETScVector*	enumVec = (PETScVector*)_enumVec;
	PETScVector*	denomVec = (PETScVector*)_denomVec;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScVector ) );
	assert( enumVec && Stg_CheckType( enumVec, PETScVector ) );
	assert( denomVec && Stg_CheckType( denomVec, PETScVector ) );
	
	ec = VecPointwiseDivide( self->petscVec, enumVec->petscVec, denomVec->petscVec ); 
	CheckPETScError(ec);
}

void PETScVector_Reciprocal( void* vector ) {
	PETScVector*	self = (PETScVector*)vector;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScVector ) );

	ec = VecReciprocal( self->petscVec ); 
	CheckPETScError(ec);
}

unsigned PETScVector_GetGlobalSize( void* vector ) {
	PETScVector*	self = (PETScVector*)vector;
	PetscErrorCode	ec;
	unsigned	size;

	assert( self && Stg_CheckType( self, PETScVector ) );

	ec = VecGetSize( self->petscVec, (PetscInt*)&size );
	CheckPETScError( ec );

	return size;
}

unsigned PETScVector_GetLocalSize( void* vector ) {
	PETScVector*	self = (PETScVector*)vector;
	PetscErrorCode	ec;
	unsigned	size;

	assert( self && Stg_CheckType( self, PETScVector ) );

	ec = VecGetLocalSize( self->petscVec, (PetscInt*)&size );
	CheckPETScError( ec );

	return size;
}

void PETScVector_GetArray( void* vector, double** array ) {
	PETScVector*	self = (PETScVector*)vector;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScVector ) );
	assert( array );

	ec = VecGetArray( self->petscVec, (PetscScalar**)array );
	CheckPETScError( ec );
}

void PETScVector_RestoreArray( void* vector, double** array ) {
	PETScVector*	self = (PETScVector*)vector;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScVector ) );
	assert( array );

	ec = VecRestoreArray( self->petscVec, (PetscScalar**)array );
	CheckPETScError( ec );
}

/* Why isn't VecDuplicate() being used Luke ?? */
void PETScVector_Duplicate( void* vector, void** dstVector ) {
	PETScVector*	self = (PETScVector*)vector;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScVector ) );
	assert( dstVector );

	*dstVector = (void*)self->_defaultConstructor( "" );
	ec = VecDestroy( ((PETScVector*)*dstVector)->petscVec );
	CheckPETScError( ec );

	ec = VecCreate( self->comm, &((PETScVector*)*dstVector)->petscVec );


	CheckPETScError( ec );
}

void PETScVector_CopyEntries( void* vector, void* _dstVector ) {
	PETScVector*	self = (PETScVector*)vector;
	PETScVector*	dstVector = (PETScVector*)_dstVector;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScVector ) );
	assert( dstVector && Stg_CheckType( dstVector, PETScVector ) );

	ec = VecCopy( self->petscVec, dstVector->petscVec );
	CheckPETScError( ec );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
