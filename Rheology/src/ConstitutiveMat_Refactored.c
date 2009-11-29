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
** $Id: ConstitutiveMat_Refactored.c 803 2008-09-11 05:22:20Z LukeHodkinson $
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
#include "ConstitutiveMat_Refactored.h"
#include "RheologyMaterial.h"

#include <assert.h>
#include <string.h>

const Type ConstitutiveMat_Refactored_Type = "ConstitutiveMat_Refactored";

ConstitutiveMat_Refactored* _ConstitutiveMat_Refactored_New(  CONSTITUTIVEMAT_REFACTORED_DEFARGS  ) {
	ConstitutiveMat_Refactored* self;
	
	assert( _sizeOfSelf >= sizeof(ConstitutiveMat_Refactored) );
	
	/* General info */
	self = (ConstitutiveMat_Refactored*)_Stg_Component_New(  STG_COMPONENT_PASSARGS  );
	
	/* Virtual functions */
	self->_setValue            = _setValue;
	self->_getViscosity        = _getViscosity;
	self->_isotropicCorrection = _isotropicCorrection;
	self->_setSecondViscosity  = _setSecondViscosity;
	self->_assemble_D_B        = _assemble_D_B;
	self->_calculateStress     = _calculateStress;
	
	return self;
}

void _ConstitutiveMat_Refactored_Init(
	void*						constitutiveMatrix,
	Dimension_Index		dim,	
	PICelleratorContext*	context,	
	Materials_Register*	materials_Register ) 
{
	ConstitutiveMat_Refactored* self = (ConstitutiveMat_Refactored*)constitutiveMatrix;
	/* General and Function pointers for this class that are not on the parent class should be set here should already be set */
	
	/* ConstitutiveMat_Refactored info */
	self->context = context;
	self->debug = Journal_MyStream( Debug_Type, self );

	self->matrixData = NULL;
	self->dim = dim;

	self->materials_Register = materials_Register;
	self->isDiagonal = False;
	self->columnSize = 0;
	self->rowSize = 0;

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


void _ConstitutiveMat_Refactored_Delete( void* constitutiveMatrix ) {
	ConstitutiveMat_Refactored* self = (ConstitutiveMat_Refactored*)constitutiveMatrix;
	
	Journal_DPrintf( self->debug, "In %s - for matrix %s\n", __func__, self->name );
	Memory_Free( self->matrixData );

	/* Stg_Class_Delete parent*/
	_Stg_Component_Delete( self );
	
}

void _ConstitutiveMat_Refactored_Print( void* constitutiveMatrix, Stream* stream ) {
	ConstitutiveMat_Refactored* self = (ConstitutiveMat_Refactored*)constitutiveMatrix;

	/* General info */
	Journal_PrintPointer( stream, constitutiveMatrix );
	Stream_Indent( stream );
	
	/* Print parent */
	_Stg_Component_Print( self, stream );
	
	/* Function pointers for this class that are not on the parent class should be set here */
	Journal_PrintPointer( stream, self->_setValue );
	Journal_PrintPointer( stream, self->_getViscosity );
	Journal_PrintPointer( stream, self->_isotropicCorrection );
	Journal_PrintPointer( stream, self->_setSecondViscosity );
	Journal_PrintPointer( stream, self->_assemble_D_B );
	Journal_PrintPointer( stream, self->_calculateStress );

	/* Regular Info */
	Journal_PrintPointer( stream, self->debug );
	ConstitutiveMat_Refactored_PrintContents( self, stream );
				
	Journal_PrintBool( stream, self->isDiagonal );
	Journal_PrintValue( stream, self->dim );
	Journal_PrintValue( stream, self->columnSize );
	Journal_PrintValue( stream, self->rowSize );
	Journal_PrintBool( stream, self->previousSolutionExists );
	
	Stream_UnIndent( stream );
}


void* _ConstitutiveMat_Refactored_Copy( void* constitutiveMatrix, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	ConstitutiveMat_Refactored*	self = (ConstitutiveMat_Refactored*)constitutiveMatrix;
	ConstitutiveMat_Refactored*	newConstitutiveMat_Refactored;
	
	/* TODO */ abort();
	if (deep) {
		newConstitutiveMat_Refactored->matrixData = Memory_Alloc_2DArray( double, self->columnSize, self->rowSize, self->name );
	}
	return (void*)newConstitutiveMat_Refactored;
}

void _ConstitutiveMat_Refactored_AssignFromXML( void* constitutiveMatrix, Stg_ComponentFactory* cf, void* data ) {
	ConstitutiveMat_Refactored*	self = (ConstitutiveMat_Refactored*)constitutiveMatrix;
	Dimension_Index					dim;
	Materials_Register*				materialsRegister;
	PICelleratorContext*				context;

	context = (PICelleratorContext*)Stg_ComponentFactory_ConstructByName( cf, "context", PICelleratorContext, True, data );

	materialsRegister = context->materials_Register;
	assert( materialsRegister );

	dim = Stg_ComponentFactory_GetRootDictUnsignedInt( cf, "dim", 0 );

	_ConstitutiveMat_Refactored_Init( self, dim, context, materialsRegister );
}

void _ConstitutiveMat_Refactored_Build( void* constitutiveMatrix, void* data ) {
	ConstitutiveMat_Refactored* self = (ConstitutiveMat_Refactored*)constitutiveMatrix;
	Material_Index      material_I;
	Material_Index      materialCount = Materials_Register_GetCount( self->materials_Register );
	RheologyMaterial*   material;
	Stream*             errorStream = Journal_Register( Error_Type, self->type );

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

		/* calls the StiffnessMatrix_SetNonLinear macro which just sets the isNonLinear flag for the 
		 * stiffnessMatrix. No longer using stiffnessMatrices, and the isNonLinear flag doesn't seem to
		 * be used anywhere, so its probably ok not to set it (I HOPE!!!) - dave. 27.02.09 */
		//if ( RheologyMaterial_IsNonLinear( material ) ) 
		//	ConstitutiveMat_Refactored_SetToNonLinear( self );
	}
}


void _ConstitutiveMat_Refactored_Initialise( void* constitutiveMatrix, void* data ) {
	ConstitutiveMat_Refactored* self          = (ConstitutiveMat_Refactored*)constitutiveMatrix;

	Journal_DPrintf( self->debug, "In %s - for matrix %s\n", __func__, self->name );
	
	ConstitutiveMat_Refactored_ZeroMatrix( self ) ;
}


void _ConstitutiveMat_Refactored_Execute( void* constitutiveMatrix, void* data ) {
	Stg_Component_Execute( constitutiveMatrix, data, False );
}

void _ConstitutiveMat_Refactored_Destroy( void* constitutiveMatrix, void* data ) {
	Stg_Component_Destroy( constitutiveMatrix, data, False );
}

/* +++ Private Functions +++ */

/* +++ Public Functions +++ */

void ConstitutiveMat_Refactored_Assemble( 
		void*                                              constitutiveMatrix,
		Element_LocalIndex                                 lElement_I,
		IntegrationPointsSwarm*				   swarm,
		int                                                particleIndex,
		IntegrationPoint*                                  particle )
{
	ConstitutiveMat_Refactored*     self          = (ConstitutiveMat_Refactored*)constitutiveMatrix;
	RheologyMaterial*       material;
	MaterialPointsSwarm*    materialSwarm;
	MaterialPoint*          materialPoint;

	/* Big fat assumption!
	 * Because of Rheology framework vs PIC IP mapping change, ConstitutiveMat_Refactored assumes that
	 * we are using a one-to-one mapping.
	 * This is because the Rheology stuff was made based on operating on MaterialPoints
	 * rather than IntegrationPoints even though its assembling. However it does require xi
	 * which we are passing from the integration point for speed (rather than re-computing it)
	 *
	 * _Init() firewalls this assumption.
	 *
	 * Proper fix for this is to reassess the Rheology/ConstitutiveMat_Refactored design
	 *
	 * -- Alan 20060421
	 */
	material = (RheologyMaterial*) IntegrationPointsSwarm_GetMaterialOn( swarm, particle );
	materialPoint = OneToOneMapper_GetMaterialPoint( swarm->mapper, particle, &materialSwarm );
	self->currentParticleIndex = particleIndex;

	/* need to change the interface for this function, so it takes in a new constitutive matrix */
	/* TODO!!! not sure if passing a ConstitutiveMat_Refactored to this function is going to be a problem !!!TODO */
	RheologyMaterial_RunRheologies( material, (ConstitutiveMatrix*)self, materialSwarm, lElement_I, materialPoint, particle->xi );
	Journal_DPrintfL( self->debug, 3, "Viscosity = %g\n", ConstitutiveMat_Refactored_GetIsotropicViscosity( self ) );
}

void ConstitutiveMat_Refactored_ZeroMatrix( void* constitutiveMatrix ) {
	ConstitutiveMat_Refactored* self   = (ConstitutiveMat_Refactored*)constitutiveMatrix;

	memset( self->matrixData[0], 0, (self->columnSize * self->rowSize)*sizeof(double) );
        memset( self->derivs, 0, 3 * 3 * sizeof(double) );
	self->isDiagonal = True;
}

void ConstitutiveMat_Refactored_SetIsotropicViscosity( void* constitutiveMatrix, double viscosity ) {
	ConstitutiveMat_Refactored* self = (ConstitutiveMat_Refactored*)constitutiveMatrix;

	ConstitutiveMat_Refactored_ZeroMatrix( self );
	ConstitutiveMat_Refactored_IsotropicCorrection( self, viscosity );
	self->isDiagonal = True;
}	

void ConstitutiveMat_Refactored_MultiplyByValue( void* constitutiveMatrix, double factor ) {
	ConstitutiveMat_Refactored* self       = (ConstitutiveMat_Refactored*)constitutiveMatrix;
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
	
void ConstitutiveMat_Refactored_PrintContents( void* constitutiveMatrix, Stream* stream ) {
	ConstitutiveMat_Refactored* self   = (ConstitutiveMat_Refactored*)constitutiveMatrix;
	Index               row_I;
	Index               col_I;

	for ( col_I = 0 ; col_I < self->columnSize ; col_I++ ) {
		for ( row_I = 0 ; row_I < self->rowSize ; row_I++ ) {
			Journal_Printf( stream, "matrixData[ %u ][ %u ] = %.4g; \t", col_I, row_I, self->matrixData[ col_I ][ row_I ] );
		}
		Journal_Printf( stream, "\n" );
	}
}



