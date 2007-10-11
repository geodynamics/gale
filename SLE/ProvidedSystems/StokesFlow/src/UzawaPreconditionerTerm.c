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
** $Id: UzawaPreconditionerTerm.c 964 2007-10-11 08:03:06Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include "StgFEM/SLE/LinearAlgebra/LinearAlgebra.h"
#include <StgFEM/SLE/SystemSetup/SystemSetup.h>

#include "types.h"
#include "UzawaPreconditionerTerm.h"
#include "Stokes_SLE.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class */
const Type UzawaPreconditionerTerm_Type = "UzawaPreconditionerTerm";

UzawaPreconditionerTerm* UzawaPreconditionerTerm_New( 
		Name                                                name,
		StiffnessMatrix*                                    stiffnessMatrix,
		Swarm*                                              integrationSwarm )
{
	UzawaPreconditionerTerm* self = (UzawaPreconditionerTerm*) _UzawaPreconditionerTerm_DefaultNew( name );

	UzawaPreconditionerTerm_InitAll( 
			self,
			stiffnessMatrix,
			integrationSwarm );

	return self;
}

/* Creation implementation / Virtual constructor */
UzawaPreconditionerTerm* _UzawaPreconditionerTerm_New( 
		SizeT                                               sizeOfSelf,  
		Type                                                type,
		Stg_Class_DeleteFunction*                           _delete,
		Stg_Class_PrintFunction*                            _print,
		Stg_Class_CopyFunction*                             _copy, 
		Stg_Component_DefaultConstructorFunction*           _defaultConstructor,
		Stg_Component_ConstructFunction*                    _construct,
		Stg_Component_BuildFunction*                        _build,
		Stg_Component_InitialiseFunction*                   _initialise,
		Stg_Component_ExecuteFunction*                      _execute,
		Stg_Component_DestroyFunction*                      _destroy,
		StiffnessMatrixTerm_AssembleElementFunction*        _assembleElement,
		Name                                                name )
{
	UzawaPreconditionerTerm* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(UzawaPreconditionerTerm) );
	self = (UzawaPreconditionerTerm*) _StiffnessMatrixTerm_New( 
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
		_assembleElement,
		name );
	
	/* Virtual info */
	
	return self;
}

void _UzawaPreconditionerTerm_Init( 
		UzawaPreconditionerTerm*                                    self )
{
}

void UzawaPreconditionerTerm_InitAll( 
		void*                                               matrixTerm,
		StiffnessMatrix*                                    stiffnessMatrix,
		Swarm*                                              integrationSwarm )
{
	UzawaPreconditionerTerm* self = (UzawaPreconditionerTerm*) matrixTerm;

	StiffnessMatrixTerm_InitAll( self, stiffnessMatrix, integrationSwarm, NULL );
	_UzawaPreconditionerTerm_Init( self );
}

void _UzawaPreconditionerTerm_Delete( void* matrixTerm ) {
	UzawaPreconditionerTerm* self = (UzawaPreconditionerTerm*)matrixTerm;

	_StiffnessMatrixTerm_Delete( self );
}

void _UzawaPreconditionerTerm_Print( void* matrixTerm, Stream* stream ) {
	UzawaPreconditionerTerm* self = (UzawaPreconditionerTerm*)matrixTerm;
	
	_StiffnessMatrixTerm_Print( self, stream );

	/* General info */
}

void* _UzawaPreconditionerTerm_DefaultNew( Name name ) {
	return (void*)_UzawaPreconditionerTerm_New( 
		sizeof(UzawaPreconditionerTerm), 
		UzawaPreconditionerTerm_Type,
		_UzawaPreconditionerTerm_Delete,
		_UzawaPreconditionerTerm_Print,
		NULL,
		_UzawaPreconditionerTerm_DefaultNew,
		_UzawaPreconditionerTerm_Construct,
		_UzawaPreconditionerTerm_Build,
		_UzawaPreconditionerTerm_Initialise,
		_UzawaPreconditionerTerm_Execute,
		_UzawaPreconditionerTerm_Destroy,
		_UzawaPreconditionerTerm_AssembleElement,
		name );
}

void _UzawaPreconditionerTerm_Construct( void* matrixTerm, Stg_ComponentFactory* cf, void* data ) {
	UzawaPreconditionerTerm*            self             = (UzawaPreconditionerTerm*)matrixTerm;

	/* Construct Parent */
	_StiffnessMatrixTerm_Construct( self, cf, data );

	_UzawaPreconditionerTerm_Init( self );
}

void _UzawaPreconditionerTerm_Build( void* matrixTerm, void* data ) {
	UzawaPreconditionerTerm*             self             = (UzawaPreconditionerTerm*)matrixTerm;

	_StiffnessMatrixTerm_Build( self, data );
}

void _UzawaPreconditionerTerm_Initialise( void* matrixTerm, void* data ) {
	UzawaPreconditionerTerm*             self             = (UzawaPreconditionerTerm*)matrixTerm;

	_StiffnessMatrixTerm_Initialise( self, data );
}

void _UzawaPreconditionerTerm_Execute( void* matrixTerm, void* data ) {
	_StiffnessMatrixTerm_Execute( matrixTerm, data );
}

void _UzawaPreconditionerTerm_Destroy( void* matrixTerm, void* data ) {
	_StiffnessMatrixTerm_Destroy( matrixTerm, data );
}


void _UzawaPreconditionerTerm_AssembleElement( 
		void*                                              matrixTerm,
		StiffnessMatrix*                                   stiffnessMatrix, 
		Element_LocalIndex                                 lElement_I, 
		SystemLinearEquations*                             _sle,
		FiniteElementContext*                              context,
		double**                                           elPreconditioner ) 
{
	Stokes_SLE*            sle              = Stg_DCheckType( _sle, Stokes_SLE );
	StiffnessMatrix*       gMatrix          = sle->gStiffMat;
	FeVariable*            gFeVariable_col  = gMatrix->columnVariable;
	ElementType*           gElementType_col = FeMesh_GetElementType( gFeVariable_col->feMesh, lElement_I );
	Node_ElementLocalIndex gColCount        = gElementType_col->nodeCount;
	double**               gElementMatrix;

	StiffnessMatrix*       kMatrix          = sle->kStiffMat;
	FeVariable*            kFeVariable_row  = kMatrix->rowVariable;
	ElementType*           kElementType_row = FeMesh_GetElementType( kFeVariable_row->feMesh, lElement_I );
	Node_ElementLocalIndex kRowCount        = kElementType_row->nodeCount;

	Index                  velocityDofCount;
	Index                  pressureDofCount;
	double**               kElementMatrix;

	double**               kInverseG;
	Index                  row_I;
	Index                  col_I;
	Index                  vel_I;

	pressureDofCount = gColCount;
	velocityDofCount = kRowCount * kMatrix->dim;

	/* Allocate Memory */
	gElementMatrix = Memory_Alloc_2DArray( double, velocityDofCount, pressureDofCount, "g element matrix" );
	memset( gElementMatrix[0], 0, sizeof(double) * velocityDofCount * pressureDofCount );
	kElementMatrix = Memory_Alloc_2DArray( double, velocityDofCount, velocityDofCount, "k element matrix" );
	memset( kElementMatrix[0], 0, sizeof(double) * velocityDofCount * velocityDofCount );
	kInverseG      = Memory_Alloc_2DArray( double, velocityDofCount, pressureDofCount, "[ diag(K) ]^-1 G" );

	/* Assemble Element G */
	StiffnessMatrix_AssembleElement( gMatrix, lElement_I, _sle, context, gElementMatrix );

	/* Assemble Element K */
	StiffnessMatrix_AssembleElement( kMatrix, lElement_I, _sle, context, kElementMatrix );

	/* Assemble [ diag(K) ]^-1 G */
	for ( row_I = 0 ; row_I < velocityDofCount ; row_I++ ) {
		for ( col_I = 0 ; col_I < pressureDofCount ; col_I++ ) {
			kInverseG[ row_I ][ col_I ] = gElementMatrix[ row_I ][ col_I ] / kElementMatrix[ row_I ][ row_I ];
		}
	}

	/* Assemble Gt [ diag(K) ]^-1 G */
	for ( row_I = 0 ; row_I < pressureDofCount ; row_I++ ) {
		for ( col_I = 0 ; col_I < pressureDofCount ; col_I++ ) {

			for ( vel_I = 0 ; vel_I < velocityDofCount ; vel_I++ ) 
				elPreconditioner[ row_I ][ col_I ] += gElementMatrix[ vel_I ][ col_I ] * kInverseG[ vel_I ][ row_I ];
		}
	}
	Memory_Free( kInverseG );
	Memory_Free( gElementMatrix );
	Memory_Free( kElementMatrix );


	/* Correct for Compressibility */
	if ( sle->cStiffMat ) {
		double**               mElementMatrix;
		StiffnessMatrix*       cMatrix        = sle->cStiffMat;
		
		mElementMatrix = Memory_Alloc_2DArray( double, pressureDofCount, pressureDofCount, "m element matrix" );
		memset( mElementMatrix[0], 0, sizeof(double) * pressureDofCount * pressureDofCount );

		StiffnessMatrix_AssembleElement( cMatrix, lElement_I, _sle, context, mElementMatrix );

		for ( row_I = 0 ; row_I < pressureDofCount ; row_I++ ) {
			for ( col_I = 0 ; col_I < pressureDofCount ; col_I++ ) {
				elPreconditioner[ row_I ][ col_I ] -= mElementMatrix[ row_I ][ col_I ];
			}
		}

		Memory_Free( mElementMatrix );
	}
}
