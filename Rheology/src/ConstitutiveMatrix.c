/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005, Monash Cluster Computing 
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
** 		* Redistributions of source code must retain the above copyright notice, 
** 			this list of conditions and the following disclaimer.
** 		* Redistributions in binary form must reproduce the above copyright 
**			notice, this list of conditions and the following disclaimer in the 
**			documentation and/or other materials provided with the distribution.
** 		* Neither the name of the Monash University nor the names of its contributors 
**			may be used to endorse or promote products derived from this software 
**			without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
** THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
** PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS 
** BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
** CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
** SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
** HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
** LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
** OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
**
** Contact:
*%		Louis Moresi - Louis.Moresi@sci.monash.edu.au
*%
** Contributors:
*+		Robert Turnbull
*+		Vincent Lemiale
*+		Louis Moresi
*+		David May
*+		David Stegman
*+		Mirko Velic
*+		Patrick Sunter
*+		Julian Giordani
*+
** $Id: ConstitutiveMatrix.c 803 2008-09-11 05:22:20Z LukeHodkinson $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "Rheology_Register.h"
#include "RheologyClass.h"
#include "ConstitutiveMatrix.h"
#include "RheologyMaterial.h"

#include <assert.h>
#include <string.h>

const Type ConstitutiveMatrix_Type = "ConstitutiveMatrix";

ConstitutiveMatrix* _ConstitutiveMatrix_New( 
		SizeT                                        _sizeOfSelf,
		Type                                         type,
		Stg_Class_DeleteFunction*                    _delete,
		Stg_Class_PrintFunction*                     _print,
		Stg_Class_CopyFunction*                      _copy, 
		Stg_Component_DefaultConstructorFunction*    _defaultConstructor,
		Stg_Component_ConstructFunction*             _construct,
		Stg_Component_BuildFunction*                 _build,
		Stg_Component_InitialiseFunction*            _initialise,
		Stg_Component_ExecuteFunction*               _execute,
		Stg_Component_DestroyFunction*               _destroy,
		StiffnessMatrixTerm_AssembleElementFunction* _assembleElement,
		ConstitutiveMatrix_SetValueFunc*             _setValue,
		ConstitutiveMatrix_GetValueFunc*             _getViscosity,
		ConstitutiveMatrix_SetValueFunc*             _isotropicCorrection,
		ConstitutiveMatrix_SetSecondViscosityFunc*   _setSecondViscosity,
		ConstitutiveMatrix_Assemble_D_B_Func*        _assemble_D_B,
		ConstitutiveMatrix_CalculateStressFunc*      _calculateStress,
		Name                                         name )
{
	ConstitutiveMatrix*	self;
	
	assert( _sizeOfSelf >= sizeof(ConstitutiveMatrix) );
	
	/* General info */
	self = (ConstitutiveMatrix*)_StiffnessMatrixTerm_New( 
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
			_assembleElement,
			name );
	
	/* Virtual functions */
	self->_setValue            = _setValue;
	self->_getViscosity        = _getViscosity;
	self->_isotropicCorrection = _isotropicCorrection;
	self->_setSecondViscosity  = _setSecondViscosity;
	self->_assemble_D_B        = _assemble_D_B;
	self->_calculateStress     = _calculateStress;
	
	return self;
}

void _ConstitutiveMatrix_Init(
      ConstitutiveMatrix*	                   self,
      Dimension_Index                        dim,	
      Bool                                   storeConstitutiveMatrix,
      FiniteElementContext*                  context,	
      Materials_Register*                    materials_Register ) 
{
	/* General and Function pointers for this class that are not on the parent class should be set here should already be set */
	
	/* ConstitutiveMatrix info */
	self->isConstructed = True;
  self->storeConstitutiveMatrix = storeConstitutiveMatrix;

	self->matrixData = NULL;
	self->dim                 = dim;
	self->isSwarmTypeIntegrationPointsSwarm = Stg_Class_IsInstance( self->integrationSwarm, IntegrationPointsSwarm_Type );
	Journal_Firewall(
		self->isSwarmTypeIntegrationPointsSwarm,
		Journal_MyStream( Error_Type, self ),
		"Error In %s - ConstitutiveMatrix %s cannot use %s. An instance of IntegrationPointsSwarm is required.\n",
		__func__,
		self->name,
		self->integrationSwarm->name );

	Journal_Firewall(
		Stg_Class_IsInstance( ((IntegrationPointsSwarm*)self->integrationSwarm)->mapper, OneToOneMapper_Type ),
		Journal_MyStream( Error_Type, self ),
		"Error In %s - ConstitutiveMatrix %s cannot use %s. ConstitutiveMatrix only works with IntegrationPointsSwarms"
		" which uses one-to-one mapping\n",
		__func__,
		self->name,
		self->integrationSwarm->name );

	self->materials_Register  = materials_Register;
	self->isDiagonal          = False;
	self->columnSize          = 0;
	self->rowSize             = 0;

	/* If we are restarting, there will be an existing valid solution for the velocity, pressure
	etc fields - thus we record this so any yield rheologies will behave correctly */
	if ( True == context->loadFromCheckPoint ) {
		self->previousSolutionExists = True;
	}
	else {
		/* Otherwise, we don't want to set this as true till we've done at least one iteration of the
		first solve */
		self->previousSolutionExists = False;
	}

	self->sle = NULL;
	self->sleNonLinearIteration_I = 0;

}


void ConstitutiveMatrix_InitAll( 
		void*	                                     constitutiveMatrix,
		StiffnessMatrix*                             stiffnessMatrix,
		Swarm*                                       swarm,
		Dimension_Index                              dim,
		FiniteElementContext*                        context,
		Materials_Register*                          materials_Register ) 
{
	ConstitutiveMatrix* self = (ConstitutiveMatrix*)constitutiveMatrix;

	_StiffnessMatrixTerm_Init( self, stiffnessMatrix, swarm, NULL );
	_ConstitutiveMatrix_Init( self, dim, False, context, materials_Register );
}

void _ConstitutiveMatrix_Delete( void* constitutiveMatrix ) {
	ConstitutiveMatrix* self = (ConstitutiveMatrix*)constitutiveMatrix;
	
	Journal_DPrintf( self->debug, "In %s - for matrix %s\n", __func__, self->name );
	Memory_Free( self->matrixData );

	/* Stg_Class_Delete parent*/
	_StiffnessMatrixTerm_Delete( self );
	
}

void _ConstitutiveMatrix_Print( void* constitutiveMatrix, Stream* stream ) {
	ConstitutiveMatrix* self = (ConstitutiveMatrix*)constitutiveMatrix;

	/* General info */
	Journal_PrintPointer( stream, constitutiveMatrix );
	Stream_Indent( stream );
	
	/* Print parent */
	_StiffnessMatrixTerm_Print( self, stream );
	
	/* Function pointers for this class that are not on the parent class should be set here */
	Journal_PrintPointer( stream, self->_setValue );
	Journal_PrintPointer( stream, self->_getViscosity );
	Journal_PrintPointer( stream, self->_isotropicCorrection );
	Journal_PrintPointer( stream, self->_setSecondViscosity );
	Journal_PrintPointer( stream, self->_assemble_D_B );
	Journal_PrintPointer( stream, self->_calculateStress );

	/* Regular Info */
	Journal_PrintPointer( stream, self->debug );
	ConstitutiveMatrix_PrintContents( self, stream );
				
	Journal_PrintBool( stream, self->isDiagonal );
	Journal_PrintValue( stream, self->dim );
	Journal_PrintValue( stream, self->columnSize );
	Journal_PrintValue( stream, self->rowSize );
	Journal_PrintBool( stream, self->previousSolutionExists );
	
	Stream_UnIndent( stream );
}


void* _ConstitutiveMatrix_Copy( void* constitutiveMatrix, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	ConstitutiveMatrix*	self = (ConstitutiveMatrix*)constitutiveMatrix;
	ConstitutiveMatrix*	newConstitutiveMatrix;
	
	/* TODO */ abort();
	if (deep) {
		newConstitutiveMatrix->matrixData = Memory_Alloc_2DArray( double, self->columnSize, self->rowSize, self->name );
	}
	return (void*)newConstitutiveMatrix;
}

void _ConstitutiveMatrix_Construct( void* constitutiveMatrix, Stg_ComponentFactory* cf, void* data ) {
	ConstitutiveMatrix*         self          = (ConstitutiveMatrix*)constitutiveMatrix;
	Dimension_Index             dim;
	Materials_Register*         materialsRegister;
	Bool                        storeConstitutiveMatrix;
	PICelleratorContext*	    context;

	_StiffnessMatrixTerm_Construct( self, cf, data );
	
	context = (PICelleratorContext*)self->context;
	assert( Stg_CheckType( context, PICelleratorContext ) );
	materialsRegister = context->materials_Register; 
	assert( materialsRegister );

	dim = Stg_ComponentFactory_GetRootDictUnsignedInt( cf, "dim", 0 );

	storeConstitutiveMatrix = Stg_ComponentFactory_GetBool( cf, self->name, "storeConstitutiveMatrix", False );

	_ConstitutiveMatrix_Init( self, dim, storeConstitutiveMatrix, (FiniteElementContext*)context, materialsRegister );
}

void _ConstitutiveMatrix_Build( void* constitutiveMatrix, void* data ) {
	ConstitutiveMatrix* self = (ConstitutiveMatrix*)constitutiveMatrix;
	Material_Index      material_I;
	Material_Index      materialCount = Materials_Register_GetCount( self->materials_Register );
	RheologyMaterial*   material;
	Stream*             errorStream = Journal_Register( Error_Type, self->type );

	_StiffnessMatrixTerm_Build( self, data );

	Journal_DPrintf( self->debug, "In %s - for matrix %s\n", __func__, self->name );

	self->matrixData = Memory_Alloc_2DArray( double, self->columnSize, self->rowSize, self->name );
/*
        self->deriv = AllocArray2D( double, (self->dim + 1), (self->dim + 1) );
*/

	for ( material_I = 0 ; material_I < materialCount ; material_I++ ) {
		material = (RheologyMaterial*) Materials_Register_GetByIndex( self->materials_Register, material_I );
		Journal_Firewall( Stg_Class_IsInstance( material, RheologyMaterial_Type ),
			errorStream,
			"Error - in %s(): while checking if each "
				"RheologyMaterial is non-linear, found a material of type %s, which is "
				"not a subclass of RheologyMaterial. Currently, if you wish to use the "
				"underworld constitutive matrix to assemble over a swarm, _every_ material "
				"must be a RheologyMaterial.\n", __func__, material->type );

		if ( RheologyMaterial_IsNonLinear( material ) ) 
			ConstitutiveMatrix_SetToNonLinear( self );
	}
}


void _ConstitutiveMatrix_Initialise( void* constitutiveMatrix, void* data ) {
	ConstitutiveMatrix* self          = (ConstitutiveMatrix*)constitutiveMatrix;

	Journal_DPrintf( self->debug, "In %s - for matrix %s\n", __func__, self->name );
	
  
	_StiffnessMatrixTerm_Initialise( self, data );

	ConstitutiveMatrix_ZeroMatrix( self ) ;
}


void _ConstitutiveMatrix_Execute( void* constitutiveMatrix, void* data ) {
	_StiffnessMatrixTerm_Execute( constitutiveMatrix, data );
}

void _ConstitutiveMatrix_Destroy( void* constitutiveMatrix, void* data ) {
	_StiffnessMatrixTerm_Destroy( constitutiveMatrix, data );
}

/* +++ Private Functions +++ */

/* +++ Public Functions +++ */

void ConstitutiveMatrix_Assemble( 
		void*                                              constitutiveMatrix,
		Element_LocalIndex                                 lElement_I,
		int                                                particleIndex,
		IntegrationPoint*                                  particle )
{
	ConstitutiveMatrix*     self          = (ConstitutiveMatrix*)constitutiveMatrix;
	IntegrationPointsSwarm* swarm         = (IntegrationPointsSwarm*)self->integrationSwarm;
	RheologyMaterial*       material;
	MaterialPointsSwarm*    materialSwarm;
	MaterialPoint*          materialPoint;

	/* Big fat assumption!
	 * Because of Rheology framework vs PIC IP mapping change, ConstitutiveMatrix assumes that
	 * we are using a one-to-one mapping.
	 * This is because the Rheology stuff was made based on operating on MaterialPoints
	 * rather than IntegrationPoints even though its assembling. However it does require xi
	 * which we are passing from the integration point for speed (rather than re-computing it)
	 *
	 * _Init() firewalls this assumption.
	 *
	 * Proper fix for this is to reassess the Rheology/ConstitutiveMatrix design
	 *
	 * -- Alan 20060421
	 */
	material = (RheologyMaterial*) IntegrationPointsSwarm_GetMaterialOn( swarm, particle );
	materialPoint = OneToOneMapper_GetMaterialPoint( swarm->mapper, particle, &materialSwarm );
	self->currentParticleIndex = particleIndex;

	RheologyMaterial_RunRheologies( material, self, materialSwarm, lElement_I, materialPoint, particle->xi );

  if( self->storeConstitutiveMatrix ) {
    /* copy the recently calculated self->matrixData, the constitutive matrix, onto the particle extension */
    double* cMatrix = ExtensionManager_Get( materialSwarm->particleExtensionMgr, materialPoint, self->storedConstHandle );
    Index row_I, rowSize = self->rowSize;
    Index columnSize = self->columnSize;

    /* flatten the matrix into a 1D array */
    for( row_I = 0 ; row_I < rowSize ; row_I++ ) 
      memcpy( &cMatrix[columnSize*row_I], self->matrixData[row_I], columnSize*sizeof(double) );
  }
	Journal_DPrintfL( self->debug, 3, "Viscosity = %g\n", ConstitutiveMatrix_GetIsotropicViscosity( self ) );
}

void ConstitutiveMatrix_ZeroMatrix( void* constitutiveMatrix ) {
	ConstitutiveMatrix* self   = (ConstitutiveMatrix*)constitutiveMatrix;

	memset( self->matrixData[0], 0, (self->columnSize * self->rowSize)*sizeof(double) );
        memset( self->derivs, 0, 3 * 3 * sizeof(double) );
	self->isDiagonal = True;
}

void ConstitutiveMatrix_SetIsotropicViscosity( void* constitutiveMatrix, double viscosity ) {
	ConstitutiveMatrix* self = (ConstitutiveMatrix*)constitutiveMatrix;

	ConstitutiveMatrix_ZeroMatrix( self );
	ConstitutiveMatrix_IsotropicCorrection( self, viscosity );
	self->isDiagonal = True;
}	

void ConstitutiveMatrix_MultiplyByValue( void* constitutiveMatrix, double factor ) {
	ConstitutiveMatrix* self       = (ConstitutiveMatrix*)constitutiveMatrix;
	Index               row_I;
	Index               col_I;
	Index               columnSize = self->columnSize;
	Index               rowSize    = self->rowSize;
	double*             columnValue;

	for ( col_I = 0 ; col_I < columnSize ; col_I++ ) {
		columnValue = self->matrixData[ col_I ];

		for ( row_I = 0 ; row_I < rowSize ; row_I++ ) {
			columnValue[ row_I ] *= factor;
		}

	}
}
	
void ConstitutiveMatrix_PrintContents( void* constitutiveMatrix, Stream* stream ) {
	ConstitutiveMatrix* self   = (ConstitutiveMatrix*)constitutiveMatrix;
	Index               row_I;
	Index               col_I;

	for ( col_I = 0 ; col_I < self->columnSize ; col_I++ ) {
		for ( row_I = 0 ; row_I < self->rowSize ; row_I++ ) {
			Journal_Printf( stream, "matrixData[ %u ][ %u ] = %.4g; \t", col_I, row_I, self->matrixData[ col_I ][ row_I ] );
		}
		Journal_Printf( stream, "\n" );
	}
}

Index ConstitutiveMatrix_GetParticleConstExtHandle( void* constitutiveMatrix ) {
  /* Function definition:
   * gets the handle that defines the "constitutive matrix on particles" extension */
	ConstitutiveMatrix* self   = (ConstitutiveMatrix*)constitutiveMatrix;
  return self->storedConstHandle;
}

void ConstitutiveMatrix_GetStoredMatrixOnParticle( 
    void* constitutiveMatrix,
    IntegrationPoint* particle,
    double** cm ) {
  /* Function definition:
   * given a integration point the function returns the stored constitutiveMatrix,
   * if it's defined on the materialPoints, which maps to the int. particle
   */
	ConstitutiveMatrix* self   = (ConstitutiveMatrix*)constitutiveMatrix;
  ExtensionInfo_Index handle = self->storedConstHandle;

  double* ext = _OneToOneMapper_GetExtensionOn( 
      ((IntegrationPointsSwarm*)self->integrationSwarm)->mapper, 
      particle, handle );

	if( self->dim == 2 ) {
		cm[0][0] = ext[0]; cm[0][1] = ext[1]; cm[0][2] = ext[2]; 
		cm[1][0] = ext[3]; cm[1][1] = ext[4]; cm[1][2] = ext[5]; 
		cm[2][0] = ext[6]; cm[2][1] = ext[7]; cm[2][2] = ext[8];
	} else {
		cm[0][0] = ext[0]; cm[0][1] = ext[1]; cm[0][2] = ext[2]; cm[0][3] = ext[3] ; cm[0][4] = ext[4]; cm[0][5] = ext[5]; 
		cm[1][0] = ext[6]; cm[1][1] = ext[7]; cm[1][2] = ext[8]; cm[1][3] = ext[9] ; cm[1][4] = ext[10]; cm[1][5] = ext[11]; 
		cm[2][0] = ext[12]; cm[2][1] = ext[13]; cm[2][2] = ext[14]; cm[2][3] = ext[15] ; cm[2][4] = ext[16]; cm[2][5] = ext[17]; 
		cm[3][0] = ext[18]; cm[3][1] = ext[19]; cm[3][2] = ext[20]; cm[3][3] = ext[21] ; cm[3][4] = ext[22]; cm[3][5] = ext[23]; 
		cm[4][0] = ext[24]; cm[4][1] = ext[25]; cm[4][2] = ext[26]; cm[4][3] = ext[27] ; cm[4][4] = ext[28]; cm[4][5] = ext[29]; 
		cm[5][0] = ext[30]; cm[5][1] = ext[31]; cm[5][2] = ext[32]; cm[5][3] = ext[33] ; cm[5][4] = ext[34]; cm[5][5] = ext[35]; 
	}
}

