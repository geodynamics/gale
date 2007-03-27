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
** $Id: Vector.c 744 2007-02-15 00:34:22Z DavidMay $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include "SLE/LinearAlgebra/LinearAlgebra.h"

#include <petsc.h>
#include <petscvec.h>
#include <StGermain/compatibility/petsccompat.h>

#include "ErrorChecking.h"

/* Use/assume PETSc as the linear algebra technology
	Warning: PETSc interface changes are handled by macros in the StGermain compatability layer
	so that outdated PETSc versions remain supported. Please bear this in mind when adding new wrappers.
 */

Vector* Vector_New_SpecifyLocalSize( MPI_Comm comm, Index size ) {
	PetscErrorCode errorFlag;
	Vec	    vector;
		
	errorFlag = VecCreateMPI( comm, size, PETSC_DETERMINE, &vector ); 
	CheckPETScError( errorFlag );

	return (Vector*)vector;
}

Vector* Vector_New_SpecifyGlobalSize( MPI_Comm comm, Index size ) {
	PetscErrorCode errorFlag;
	Vec	    vector;
		
	errorFlag = VecCreateMPI( comm, PETSC_DECIDE, size, &vector ); 
	CheckPETScError( errorFlag );
		
	return (Vector*)vector;
}


Vector* Vector_New_FromArray( MPI_Comm comm, Index size, double* array ) {
	PetscErrorCode errorFlag;
	Vec	    vector;

	errorFlag = VecCreateMPIWithArray( comm, PETSC_DECIDE, size, (PetscScalar*)array, &vector ); 
	CheckPETScError( errorFlag );
	
	return (Vector*)vector;
}


Vector *Vector_New_Ghost(MPI_Comm comm, Index size, Index ghostCount, Index *ghosts){
	PetscErrorCode errorFlag;
	Vec	    vector;
	
	errorFlag = VecCreateGhost( comm, size - ghostCount, PETSC_DECIDE, ghostCount, (int*)ghosts, &vector ); 
	CheckPETScError( errorFlag );
	
	return (Vector *)vector;
}


Vector *Vector_New_GhostFromArray(MPI_Comm comm, Index size, Index ghostCount, Index *ghosts, double *array){
	PetscErrorCode errorFlag;
	Vec	vector;
	
	errorFlag = VecCreateGhostWithArray( comm, size - ghostCount, PETSC_DECIDE, ghostCount, (int*)ghosts, (PetscScalar*)array, &vector );
	CheckPETScError( errorFlag );
	
	return (Vector *)vector;
}


Vector* Vector_New_Seq( Index size ) {
	PetscErrorCode errorFlag;
	Vec	    v;

	errorFlag = VecCreateSeq( PETSC_COMM_SELF, (int)size, &v ); CheckPETScError( errorFlag );
	CheckPETScError(errorFlag);

	return (Vector*)v;
}


void Vector_Destroy(Vector *vector) {
	PetscErrorCode errorFlag;

	errorFlag = VecDestroy((Vec)vector); CheckPETScError( errorFlag );
	CheckPETScError(errorFlag);
}

void Vector_View( Vector* vector, Stream* stream ) {
	PetscErrorCode     errorFlag;

	if ( stream == NULL ) {
		errorFlag = PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INDEX ); 
		CheckPETScError( errorFlag );
		
		errorFlag = VecView( (Vec)vector, PETSC_VIEWER_STDOUT_WORLD ); 
		CheckPETScError( errorFlag );
	}
	else {
		char*        name;
		Index        localSize = Vector_LocalSize( vector );
		double       *array;
		Index        array_I;

		Vector_Get( vector, &array );

		errorFlag = PetscObjectGetName((PetscObject)vector,&name);
		CheckPETScError( errorFlag );

		/* Print Vector */
		Journal_Printf( stream, "%s = [\n", name );
		for ( array_I = 0 ; array_I < localSize ; array_I++ ) 
			Journal_Printf( stream , "\t%u: \t %2.4g\n", array_I, array[ array_I ] );
		Journal_Printf( stream, "];\n");
		
		Vector_Restore( vector, &array );
	}
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#ifdef PETSC_USE_COMPLEX

void Vector_AddTo( Vector* vector, Index count, Index indices[], double* values ) 
{
	PetscInt i;
	PetscScalar *_values;
	
	PetscMalloc( count*sizeof(PetscScalar), &_values );
	PetscMemzero( _values, count*sizeof(PetscScalar) );
	for( i=0; i<count; i++ ) {
		PetscRealPart(_values[i]) = (PetscReal)values[i];
	}
	
	VecSetValues( (Vec)vector, count, (int*)indices, _values, ADD_VALUES ); 
	PetscFree( _values );
}

#else

void Vector_AddTo( Vector* vector, Index count, Index indices[], double* values ) {
	PetscErrorCode errorFlag;
	
	errorFlag = VecSetValues( (Vec)vector, count, (int*)indices, (PetscScalar*)values, ADD_VALUES ); 
	
	#ifdef DEBUG
		CheckPETScError( errorFlag );
	#endif
}

#endif

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#ifdef PETSC_USE_COMPLEX

void Vector_Insert( Vector* vector, Index count, Index indices[], double* values ) 
{
	PetscInt i;
	PetscScalar *_values;
	
	
	PetscMalloc( count*sizeof(PetscScalar), &_values );
	PetscMemzero( _values, count*sizeof(PetscScalar) );
	
	for( i=0; i<count; i++ ) {
		PetscRealPart(_values[i]) = (PetscReal)values[i];
	}
/*	PetscScalarView( count, _values, PETSC_VIEWER_STDOUT_WORLD ); */
	VecSetValues( (Vec)vector, count, (int*)indices, _values, INSERT_VALUES ); 
	
	PetscFree( _values );
}

#else

void Vector_Insert( Vector* vector, Index count, Index indices[], double* values ) {
	PetscErrorCode errorFlag;
	
	errorFlag = VecSetValues( (Vec)vector, count, (int*)indices, (PetscScalar*)values, INSERT_VALUES ); 

	#ifdef DEBUG
		CheckPETScError( errorFlag );
	#endif
}

#endif 

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#ifdef PETSC_USE_COMPLEX

void Vector_SetEntry( Vector* vector, Index ind, double val ) 
{
	PetscScalar _val;
	
	PetscRealPart(_val) = (PetscReal)val;
	PetscImaginaryPart(_val) = 0.0;
/*	PetscScalarView( 1, &_val, PETSC_VIEWER_STDOUT_WORLD ); */
	VecSetValue( (Vec)vector, ind, _val, INSERT_VALUES );
}

#else

void Vector_SetEntry( Vector* vector, Index ind, double val ) {
#ifdef DEBUG
	PetscErrorCode	ec;

	ec = VecSetValue( (Vec)vector, ind, (PetscScalar)val, INSERT_VALUES );
	CheckPETScError( ec );
#else
	VecSetValue( (Vec)vector, ind, (PetscScalar)val, INSERT_VALUES );
#endif
}

#endif 



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#ifdef PETSC_USE_COMPLEX

void Vector_AddEntry( Vector* vector, Index ind, double val ) 
{
	PetscScalar _val;
	
	PetscRealPart(_val) = (PetscReal)val;
	PetscImaginaryPart(_val) = 0.0;
/*	PetscScalarView( 1, &_val, PETSC_VIEWER_STDOUT_WORLD ); */
	VecSetValue( (Vec)vector, ind, _val, ADD_VALUES );
}

#else

void Vector_AddEntry( Vector* vector, Index ind, double val ) {
#ifdef DEBUG
	PetscErrorCode	ec;

	ec = VecSetValue( (Vec)vector, ind, (PetscScalar)val, ADD_VALUES );
	CheckPETScError( ec );
#else
	VecSetValue( (Vec)vector, ind, (PetscScalar)val, ADD_VALUES );
#endif
}

#endif 

void Vector_AssemblyBegin( Vector* vector ) {
	PetscErrorCode errorFlag;

	errorFlag = VecAssemblyBegin( (Vec)vector );
	CheckPETScError( errorFlag );
}

void Vector_AssemblyEnd( Vector* vector ) {
	PetscErrorCode errorFlag;

	errorFlag = VecAssemblyEnd( (Vec)vector );
	CheckPETScError( errorFlag );
}

void Vector_Zero( Vector* vector ) {
	PetscErrorCode errorFlag;

	PetscScalar zero = 0.0;

	errorFlag = VecSet( &zero, (Vec)vector );
	CheckPETScError(errorFlag);
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#ifdef PETSC_USE_COMPLEX

void Vector_Get(Vector *vector, double **array )
{
	Vec x = (Vec)vector;
	double *a;
	PetscInt M;
	PetscScalar *xv;
	PetscInt i;
	
	VecGetSize( x, &M );
	
	a = Memory_Alloc_Array( double, M, "Vector_Get_cmplx" );
	
	VecGetArray( (Vec)vector, &xv );
	for( i=0; i<M; i++ ) {
		a[i]  = (double)PetscRealPart(xv[i]);
	}
	VecRestoreArray( (Vec)vector, &xv );
	
	(*array) = a;
}

#else

void Vector_Get(Vector *vector, double **array){
	PetscErrorCode errorFlag;

	errorFlag = VecGetArray((Vec)vector, array);
	CheckPETScError(errorFlag);
}

#endif


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#ifdef PETSC_USE_COMPLEX

void Vector_GetValues( Vector* vector, unsigned nInds, unsigned* inds, double* array ) 
{
	Vec x = (Vec)vector;
	PetscScalar *a;
	PetscInt M;
	PetscInt i;
	
	VecGetSize( x, &M );
	
	PetscMalloc( M*sizeof(PetscScalar), &a );
	PetscMemzero( a, M*sizeof(PetscScalar) );
	VecGetValues( (Vec)vector, (PetscInt)nInds, (PetscInt*)inds, a );
	
	for( i=0; i<M; i++ ) {
		array[i]  = (double)PetscRealPart(a[i]);
	}
	
	PetscFree( a );
}

#else

void Vector_GetValues( Vector* vector, unsigned nInds, unsigned* inds, double* array ) {
	PetscErrorCode	ec;

	ec = VecGetValues( (Vec)vector, (PetscInt)nInds, (PetscInt*)inds, (PetscScalar*)array );
	CheckPETScError( ec );
}

#endif



void Vector_GhostGetLocal(Vector *vector, Vector **local) {
	PetscErrorCode errorFlag;

	errorFlag = VecGhostGetLocalForm((Vec)vector, (Vec *)local);
	CheckPETScError(errorFlag);
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#ifdef PETSC_USE_COMPLEX

void Vector_Restore(Vector *vector, double **array)
{
	Memory_Free( (*array) );
	(*array) = NULL;
}

#else

void Vector_Restore(Vector *vector, double **array){
	PetscErrorCode errorFlag;

	errorFlag = VecRestoreArray((Vec)vector, array);
	CheckPETScError(errorFlag);
}

#endif


void Vector_GhostRestoreLocal(Vector *vector, Vector **local) {
	PetscErrorCode errorFlag;

	errorFlag = VecGhostGetLocalForm((Vec)vector, (Vec *)local);
	CheckPETScError(errorFlag);
}

void Vector_Duplicate(Vector *vector, Vector **newVector) {
	PetscErrorCode errorFlag;

	errorFlag = VecDuplicate((Vec)vector, (Vec *)newVector);
	CheckPETScError(errorFlag);
}


Index Vector_LocalSize(Vector *vector) {
	PetscErrorCode errorFlag;
	int	    size;

	errorFlag = VecGetLocalSize((Vec)vector, &size);
	CheckPETScError(errorFlag);
	
	return (Index)size;
}

Index Vector_GlobalSize(Vector *vector) {
	PetscErrorCode errorFlag;
	int	    size;

	errorFlag = VecGetSize((Vec)vector, &size);
	CheckPETScError(errorFlag);
	
	return (Index)size;
}


void Vector_SetLocalSize( Vector* vec, unsigned size ) {
	PetscErrorCode ef;
	
	ef = VecSetSizes( (Vec)vec, (PetscInt)size, PETSC_DECIDE );
	CheckPETScError( ef );
}


void Vector_SetContents( Vector* x, double a ) {
	PetscErrorCode errorFlag;

	errorFlag = VecSet( (PetscScalar*)(&a), (Vec)x );
	CheckPETScError(errorFlag);
}       

void Vector_ScaleContents( Vector* x, double a ) {
	PetscErrorCode errorFlag;
	
	errorFlag = VecScale( (PetscScalar*)(&a), (Vec)x );
	CheckPETScError(errorFlag);
}       
 
void Vector_CopyContents( Vector* x, Vector* y ) {
	PetscErrorCode errorFlag;
	
	errorFlag = VecCopy( (Vec)x, (Vec)y );
	CheckPETScError(errorFlag);
}       
 
void Vector_DotProduct( Vector* x, Vector* y, double* dot_prod )
{ 
	PetscErrorCode errorFlag;
	PetscScalar x_dot_y;
	
	errorFlag = VecDot( (Vec)x, (Vec)y, &x_dot_y );
	*dot_prod = PetscRealPart(x_dot_y);
	
	CheckPETScError(errorFlag);
}
 
void Vector_AddScaledVector( Vector* vector, double scaleFactor, Vector* vectorToAdd ) 
{
	PetscErrorCode errorFlag;

	errorFlag = VecAXPY( (PetscScalar)scaleFactor, (Vec)vectorToAdd, (Vec)vector );
	CheckPETScError(errorFlag);
}       
        
void Vector_ScaleAndAddVector( Vector* vector, double scaleFactor, Vector* vectorToAdd ) {
	PetscErrorCode errorFlag;
        
	errorFlag = VecAYPX( (PetscScalar)scaleFactor, (Vec)vectorToAdd, (Vec)vector );

	/* Error Checking */
	CheckPETScError(errorFlag);
}
        
double Vector_L2_Norm( Vector* x ) {       
	PetscErrorCode errorFlag;
	double         norm;
    
	errorFlag = VecNorm( (Vec)x, NORM_2, (PetscReal*) &norm );
	CheckPETScError(errorFlag);

	return norm;
}

double Vector_Linf_Norm( Vector* x ) {       
	PetscErrorCode errorFlag;
	double         norm;
    
	errorFlag = VecNorm( (Vec)x, NORM_INFINITY, (PetscReal*) &norm );
	CheckPETScError(errorFlag);

	return norm;
}


void Vector_MaxComponent ( Vector* x, int* position, double* maxValue) {
	PetscErrorCode errorFlag;
	
	errorFlag = VecMax( (Vec) x, (PetscInt *) position, (PetscReal *) maxValue );
	CheckPETScError(errorFlag);	
		
}



void Vector_PointwiseDivide( Vector* x, Vector* y, Vector* result ) {
	PetscErrorCode errorFlag;
	
	errorFlag = VecPointwiseDivide( (Vec) x, (Vec) y, (Vec) result ); 
	CheckPETScError(errorFlag);	
}	

void Vector_PointwiseMultiply( Vector* x, Vector* y, Vector* result ) {
		PetscErrorCode errorFlag;

		errorFlag = VecPointwiseMult( (Vec) x, (Vec) y, (Vec) result ); 		
		CheckPETScError(errorFlag);
	}

void Vector_Transfer( Vector*		dst, 
		      Vector*	src, 
		      Index		idxCnt, 
		      Index		dstIndices[], 
		      Index		srcIndices[], 
		      MPI_Comm		comm )
{
	PetscErrorCode errorFlag;
	
	VecScatter	scatter;
	IS		isFrom;
	IS		isTo;
	
	ISCreateGeneral( comm, idxCnt, (int*)srcIndices, &isFrom );
	ISCreateGeneral( comm, idxCnt, (int*)dstIndices, &isTo );

	errorFlag = VecScatterCreate( (Vec)src, isFrom, (Vec)dst, isTo, &scatter ); CheckPETScError(errorFlag);
	errorFlag = VecScatterBegin( (Vec)src, (Vec)dst, INSERT_VALUES, SCATTER_FORWARD, scatter ); CheckPETScError(errorFlag);
	errorFlag = VecScatterEnd( (Vec)src, (Vec)dst, INSERT_VALUES, SCATTER_FORWARD, scatter ); CheckPETScError(errorFlag);

	ISDestroy( isFrom );
	ISDestroy( isTo );

	errorFlag = VecScatterDestroy( scatter ); CheckPETScError(errorFlag);
}


void Vector_DupNewSize( MPI_Comm comm, Vector* src, unsigned localSize, Vector** dst ) {
	PetscErrorCode	ef;
	VecType		type;
	
	ef = VecCreate( comm, (Vec*)dst );						CheckPETScError( ef );
	ef = VecSetSizes( (Vec)*dst, localSize, PETSC_DECIDE );	CheckPETScError( ef );
	ef = VecGetType( (Vec)src, &type );					CheckPETScError( ef );
	ef = VecSetType( (Vec)*dst, type );					CheckPETScError( ef );
}

void Vector_Reciprocal( Vector* vec ) {
	PetscErrorCode errorFlag;
	
	errorFlag = VecReciprocal( (Vec) vec ); 
	CheckPETScError(errorFlag);
}
 


