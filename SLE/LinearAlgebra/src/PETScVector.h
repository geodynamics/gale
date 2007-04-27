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
*/
/** \file
**  Role:
**
** Assumptions:
**
** Invariants:
**
** Comments:
**
** $Id: PETScVector.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_LinearAlgebra_PETScVector_h__
#define __StgFEM_SLE_LinearAlgebra_PETScVector_h__

	/** Textual name of this class */
	extern const Type PETScVector_Type;

	/** Virtual function types */

	/** PETScVector class contents */
	#define __PETScVector			\
		/* General info */		\
		__Vector			\
						\
		/* Virtual info */		\
						\
		/* PETScVector info */		\
		double*		array;		\
		Vec		petscVec;

	struct PETScVector { __PETScVector };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define PETSCVECTOR_DEFARGS	\
		VECTOR_DEFARGS

	#define PETSCVECTOR_PASSARGS	\
		VECTOR_PASSARGS

	PETScVector* PETScVector_New( Name name );
	PETScVector* _PETScVector_New( VECTOR_DEFARGS );
	void _PETScVector_Init( PETScVector* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _PETScVector_Delete( void* vector );
	void _PETScVector_Print( void* vector, Stream* stream );
	void _PETScVector_Construct( void* vector, Stg_ComponentFactory* cf, void* data );
	void _PETScVector_Build( void* vector, void* data );
	void _PETScVector_Initialise( void* vector, void* data );
	void _PETScVector_Execute( void* vector, void* data );
	void _PETScVector_Destroy( void* vector, void* data );

	void PETScVector_SetComm( void* vector, MPI_Comm comm );
	void PETScVector_SetGlobalSize( void* vector, unsigned size );
	void PETScVector_SetLocalSize( void* vector, unsigned size );
	void PETScVector_AddEntries( void* vector, unsigned nEntries, unsigned* indices, double* values );
	void PETScVector_InsertEntries( void* vector, unsigned nEntries, unsigned* indices, double* values );
	void PETScVector_SetScalar( void* vector, double scalar );
	void PETScVector_Zero( void* vector );
	void PETScVector_AssemblyBegin( void* vector );
	void PETScVector_AssemblyEnd( void* vector );

	void PETScVector_Add( void* vector, void* vector0 );
	void PETScVector_AddScaled( void* vector, double scalar, void* vector0 );
	void PETScVector_ScaleAdd( void* vector, double scalar, void* vector0 );
	void PETScVector_Subtract( void* vector, void* vector0 );
	void PETScVector_Scale( void* vector, double factor );
	double PETScVector_DotProduct( void* vector, void* vector0 );
	double PETScVector_L2Norm( void* vector );
	double PETScVector_LInfNorm( void* vector );
	void PETScVector_PointwiseMultiply( void* vector, void* enumVec, void* denomVec );
	void PETScVector_PointwiseDivide( void* vector, void* enumVec, void* denomVec );
	void PETScVector_Reciprocal( void* vector );

	unsigned PETScVector_GetGlobalSize( void* vector );
	unsigned PETScVector_GetLocalSize( void* vector );
	void PETScVector_GetArray( void* vector, double** array );
	void PETScVector_RestoreArray( void* vector, double** array );
	void PETScVector_Duplicate( void* vector, void** dstVector );
	void PETScVector_CopyEntries( void* vector, void* dstVector );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void PETScVector_SetComm( void* vector, MPI_Comm comm );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

#endif /* __StgFEM_SLE_LinearAlgebra_PETScVector_h__ */
