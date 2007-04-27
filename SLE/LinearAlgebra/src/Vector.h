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
** $Id: Vector.h 822 2007-04-27 06:20:35Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_LinearAlgebra_Vector_h__
#define __StgFEM_SLE_LinearAlgebra_Vector_h__

	/** Textual name of this class */
	extern const Type Vector_Type;

	/** Virtual function types */
	typedef void (Vector_SetCommFunc)( void* vector, MPI_Comm comm );
	typedef void (Vector_SetGlobalSizeFunc)( void* vector, unsigned size );
	typedef void (Vector_SetLocalSizeFunc)( void* vector, unsigned size );
	typedef void (Vector_AddEntriesFunc)( void* vector, unsigned nEntries, unsigned* indices, double* values );
	typedef void (Vector_InsertEntriesFunc)( void* vector, unsigned nEntries, unsigned* indices, double* values );
	typedef void (Vector_SetScalarFunc)( void* vector, double scalar );
	typedef void (Vector_ZeroFunc)( void* vector );
	typedef void (Vector_AssemblyBeginFunc)( void* vector );
	typedef void (Vector_AssemblyEndFunc)( void* vector );

	typedef void (Vector_AddFunc)( void* vector, void* vector0 );
	typedef void (Vector_AddScaledFunc)( void* vector, double scalar, void* vector0 );
	typedef void (Vector_ScaleAddFunc)( void* vector, double scalar, void* vector0 );
	typedef void (Vector_SubtractFunc)( void* vector, void* vector0 );
	typedef void (Vector_ScaleFunc)( void* vector, double factor );
	typedef double (Vector_DotProductFunc)( void* vector, void* vector0 );
	typedef double (Vector_L2NormFunc)( void* vector );
	typedef double (Vector_LInfNormFunc)( void* vector );
	typedef void (Vector_PointwiseMultiplyFunc)( void* vector, void* enumVec, void* denomVec );
	typedef void (Vector_PointwiseDivideFunc)( void* vector, void* enumVec, void* denomVec );
	typedef void (Vector_ReciprocalFunc)( void* vector );

	typedef unsigned (Vector_GetGlobalSizeFunc)( void* vector );
	typedef unsigned (Vector_GetLocalSizeFunc)( void* vector );
	typedef void (Vector_GetArrayFunc)( void* vector, double** array );
	typedef void (Vector_RestoreArrayFunc)( void* vector, double** array );
	typedef void (Vector_DuplicateFunc)( void* vector, void** dstVector );
	typedef void (Vector_CopyEntriesFunc)( void* vector, void* dstVector );
	typedef void (Vector_ViewFunc)( void* vector, void* stream );

	/** Vector class contents */
	#define __Vector						\
		/* General info */					\
		__Stg_Component						\
									\
		/* Virtual info */					\
		Vector_SetCommFunc*		setCommFunc;		\
		Vector_SetGlobalSizeFunc*	setGlobalSizeFunc;	\
		Vector_SetLocalSizeFunc*	setLocalSizeFunc;	\
		Vector_AddEntriesFunc*		addEntriesFunc;		\
		Vector_InsertEntriesFunc*	insertEntriesFunc;	\
		Vector_SetScalarFunc*		setScalarFunc;		\
		Vector_ZeroFunc*		zeroFunc;		\
		Vector_AssemblyBeginFunc*	assemblyBeginFunc;	\
		Vector_AssemblyEndFunc*		assemblyEndFunc;	\
									\
		Vector_AddFunc*			addFunc;		\
		Vector_AddScaledFunc*		addScaledFunc;		\
		Vector_ScaleAddFunc*		scaleAddFunc;		\
		Vector_SubtractFunc*		subtractFunc;		\
		Vector_ScaleFunc*		scaleFunc;		\
		Vector_DotProductFunc*		dotProductFunc;		\
		Vector_L2NormFunc*		l2NormFunc;		\
		Vector_LInfNormFunc*		lInfNormFunc;		\
		Vector_PointwiseMultiplyFunc*	pointwiseMultiplyFunc;	\
		Vector_PointwiseDivideFunc*	pointwiseDivideFunc;	\
		Vector_ReciprocalFunc*		reciprocalFunc;		\
									\
		Vector_GetGlobalSizeFunc*	getSizeFunc;		\
		Vector_GetLocalSizeFunc*	getLocalSizeFunc;	\
		Vector_GetArrayFunc*		getArrayFunc;		\
		Vector_RestoreArrayFunc*	restoreArrayFunc;	\
		Vector_DuplicateFunc*		duplicateFunc;		\
		Vector_CopyEntriesFunc*		copyEntriesFunc;	\
		Vector_ViewFunc*		viewFunc;		\
									\
		/* Vector info */					\
		MPI_Comm		comm;

	struct Vector { __Vector };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define VECTOR_DEFARGS						\
		STG_COMPONENT_DEFARGS,					\
		Vector_SetCommFunc*		setCommFunc, 		\
		Vector_SetGlobalSizeFunc*	setGlobalSizeFunc,	\
		Vector_SetLocalSizeFunc*	setLocalSizeFunc,	\
		Vector_AddEntriesFunc*		addEntriesFunc,		\
		Vector_InsertEntriesFunc*	insertEntriesFunc,	\
		Vector_SetScalarFunc*		setScalarFunc,		\
		Vector_ZeroFunc*		zeroFunc,		\
		Vector_AssemblyBeginFunc*	assemblyBeginFunc,	\
		Vector_AssemblyEndFunc*		assemblyEndFunc,	\
									\
		Vector_AddFunc*			addFunc,		\
		Vector_AddScaledFunc*		addScaledFunc,		\
		Vector_ScaleAddFunc*		scaleAddFunc,		\
		Vector_SubtractFunc*		subtractFunc,		\
		Vector_ScaleFunc*		scaleFunc,		\
		Vector_DotProductFunc*		dotProductFunc,		\
		Vector_L2NormFunc*		l2NormFunc,		\
		Vector_LInfNormFunc*		lInfNormFunc,		\
		Vector_PointwiseMultiplyFunc*	pointwiseMultiplyFunc,	\
		Vector_PointwiseDivideFunc*	pointwiseDivideFunc,	\
		Vector_ReciprocalFunc*		reciprocalFunc,		\
									\
		Vector_GetGlobalSizeFunc*	getSizeFunc,		\
		Vector_GetLocalSizeFunc*	getLocalSizeFunc,	\
		Vector_GetArrayFunc*		getArrayFunc,		\
		Vector_RestoreArrayFunc*	restoreArrayFunc,	\
		Vector_DuplicateFunc*		duplicateFunc,		\
		Vector_CopyEntriesFunc*		copyEntriesFunc,	\
		Vector_ViewFunc*		viewFunc

	#define VECTOR_PASSARGS			\
		STG_COMPONENT_PASSARGS, 	\
		setCommFunc, 			\
		setGlobalSizeFunc,		\
		setLocalSizeFunc,		\
		addEntriesFunc,			\
		insertEntriesFunc,		\
		setScalarFunc,			\
		zeroFunc,			\
		assemblyBeginFunc,		\
		assemblyEndFunc,		\
						\
		addFunc,			\
		addScaledFunc,			\
		scaleAddFunc,			\
		subtractFunc,			\
		scaleFunc,			\
		dotProductFunc,			\
		l2NormFunc,			\
		lInfNormFunc,			\
		pointwiseMultiplyFunc,		\
		pointwiseDivideFunc,		\
		reciprocalFunc,			\
						\
		getSizeFunc,			\
		getLocalSizeFunc,		\
		getArrayFunc,			\
		restoreArrayFunc,		\
		duplicateFunc,			\
		copyEntriesFunc,		\
		viewFunc

	Vector* _Vector_New( VECTOR_DEFARGS );
	void _Vector_Init( Vector* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _Vector_Delete( void* vector );
	void _Vector_Print( void* vector, Stream* stream );
	void _Vector_Construct( void* vector, Stg_ComponentFactory* cf, void* data );
	void _Vector_Build( void* vector, void* data );
	void _Vector_Initialise( void* vector, void* data );
	void _Vector_Execute( void* vector, void* data );
	void _Vector_Destroy( void* vector, void* data );

	void _Vector_SetComm( void* vector, MPI_Comm comm );
	void _Vector_View( void* vector, void* stream );

	#define Vector_SetComm( self, comm )			\
		VirtualCall( self, setCommFunc, self, comm )

	#define Vector_SetGlobalSize( self, size )	\
		VirtualCall( self, setGlobalSizeFunc, self, size )

	#define Vector_SetLocalSize( self, size )	\
		VirtualCall( self, setLocalSizeFunc, self, size )

	#define Vector_AddEntries( self, nEntries, indices, values )	\
		VirtualCall( self, addEntriesFunc, self, nEntries, indices, values )

	#define Vector_InsertEntries( self, nEntries, indices, values )	\
		VirtualCall( self, insertEntriesFunc, self, nEntries, indices, values )

	#define Vector_SetScalar( self, scalar )	\
		VirtualCall( self, setScalarFunc, self, scalar )

	#define Vector_Zero( self )	\
		VirtualCall( self, zeroFunc, self )

	#define Vector_AssemblyBegin( self )	\
		VirtualCall( self, assemblyBeginFunc, self )

	#define Vector_AssemblyEnd( self )	\
		VirtualCall( self, assemblyEndFunc, self )

	#define Vector_Add( self, vector0 )		\
		VirtualCall( self, addFunc, self, vector0 )

	#define Vector_AddScaled( self, scalar, vector0 )	\
		VirtualCall( self, addScaledFunc, self, scalar, vector0 )

	#define Vector_ScaleAdd( self, scalar, vector0 )	\
		VirtualCall( self, scaleAddFunc, self, scalar, vector0 )

	#define Vector_Subtract( self, vector0 )	\
		VirtualCall( self, subtractFunc, self, vector0 )

	#define Vector_Scale( self, factor )	\
		VirtualCall( self, scaleFunc, self, factor )

	#define Vector_DotProduct( self, vector0 )	\
		VirtualCall( self, dotProductFunc, self, vector0 )

	#define Vector_L2Norm( self )	\
		VirtualCall( self, l2NormFunc, self )

	#define Vector_LInfNorm( self )	\
		VirtualCall( self, lInfNormFunc, self )

	#define Vector_PointwiseMultiply( self, enumVec, denomVec )			\
		VirtualCall( self, pointwiseMultiplyFunc, self, enumVec, denomVec )

	#define Vector_PointwiseDivide( self, enumVec, denomVec )			\
		VirtualCall( self, pointwiseDivideFunc, self, enumVec, denomVec )

	#define Vector_Reciprocal( self )	\
		VirtualCall( self, reciprocalFunc, self )

	#define Vector_GetGlobalSize( self )	\
		VirtualCall( self, getSizeFunc, self )

	#define Vector_GetLocalSize( self )	\
		VirtualCall( self, getLocalSizeFunc, self )

	#define Vector_GetArray( self, array )			\
		VirtualCall( self, getArrayFunc, self, array )

	#define Vector_RestoreArray( self, array )	\
		VirtualCall( self, restoreArrayFunc, self, array )

	#define Vector_Duplicate( self, dstVector )	\
		VirtualCall( self, duplicateFunc, self, dstVector )

	#define Vector_CopyEntries( self, dstVector )	\
		VirtualCall( self, copyEntriesFunc, self, dstVector )

	#define Vector_View( self, stream )			\
		VirtualCall( self, viewFunc, self, stream )

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void Vector_Dump( void* vector, const char* filename );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

#endif /* __StgFEM_SLE_LinearAlgebra_Vector_h__ */
