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
** $Id: Vector.c 656 2006-10-18 06:45:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include "SLE/LinearAlgebra/LinearAlgebra.h"

Vector* Vector_New( MPI_Comm comm, Index size ) {
	return (Vector*)0;
}


Vector* Vector_New_FromArray( MPI_Comm comm, Index size, double* array ) {
	return (Vector*)0;
}


Vector *Vector_New_Ghost( MPI_Comm comm, Index size, Index ghostCount, Index *ghosts ) {
	return 0;
}


Vector *Vector_New_GhostFromArray( MPI_Comm comm, Index size, Index ghostCount, Index *ghosts, double *array ) {
	return 0;
}


Vector* Vector_New_Seq( Index size ) {
	return (Vector*)0;
}


void Vector_Destroy( Vector *vector ) {
}


void Vector_View( Vector* vector, Stream* stream ) {
}

void Vector_AddTo( Vector* vector, Index count, Index indices[], double* values ) {
}

void Vector_Insert( Vector* vector, Index count, Index indices[], double* values ) {
}

void Vector_SetEntry( Vector* vector, Index ind, double val ) {
}

void Vector_AddEntry( Vector* vector, Index ind, double val ) {
}

void Vector_AssemblyBegin( Vector* vector ) {
}

void Vector_AssemblyEnd( Vector* vector ) {
}

void Vector_Zero( Vector* vector ) {
}

void Vector_Get( Vector *vector, double **array ) {
}

void Vector_GetValues( Vector* vector, unsigned nInds, unsigned* inds, double* array ) {
}


void Vector_GhostGetLocal( Vector *vector, Vector **local ) {
}


void Vector_Restore( Vector *vector, double **array ) {
}


void Vector_GhostRestoreLocal( Vector *vector, Vector **local ) {
}

void Vector_Duplicate(Vector *vector, Vector **newVector) {
}

Index Vector_LocalSize( Vector *vector ) {
	return 0;
}	

Index Vector_GlobalSize( Vector *vector ) {
		return 0;
	}

void Vector_SetContents( Vector* x, double a ) {
}
                                                                                               
void Vector_ScaleContents( Vector* x, double a ) {
}
                                                                                               
void Vector_CopyContents( Vector* x, Vector* y ) {
}
                                                                                               
void Vector_DotProduct( Vector* x, Vector* y, double* dot_prod ) {
}

void Vector_AddScaledVector( Vector* vector, double scaleFactor, Vector* vectorToAdd ) {
}

void Vector_ScaleAndAddVector( Vector* vector, double scaleFactor, Vector* vectorToAdd ) {
}
                                                                                               
double Vector_L2_Norm( Vector* x ) {
	return 0.0;
}

double Vector_Linf_Norm( Vector* x ) {
	return 0.0;
}

void Vector_PointwiseDivide( Vector* x, Vector* y, Vector* result ) {
}

void Vector_PointwiseMultiply( Vector* x, Vector* y, Vector* result ) {
}


void Vector_Transfer( Vector*		dst, 
		      Vector*		src, 
		      Index		idxCnt, 
		      Index		dstIndices[], 
		      Index		srcIndices[], 
		      MPI_Comm		comm )
{
}

void Vector_SetLocalSize( Vector* vec, unsigned size ) {
}

void Vector_Reciprocal( Vector* vec ) {
}

void Vector_DupNewSize( MPI_Comm comm, Vector* src, unsigned localSize, Vector** dst ) {
}
