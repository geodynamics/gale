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
** $Id: ConstitutiveMatCartesian_Refactored.c 803 2008-09-11 05:22:20Z LukeHodkinson $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/



#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "Rheology_Register.h"
#include "ConstitutiveMat_Refactored.h"
#include "ConstitutiveMatCartesian_Refactored.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>


/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type ConstitutiveMatCartesian_Refactored_Type = "ConstitutiveMatCartesian_Refactored";

ConstitutiveMatCartesian_Refactored* ConstitutiveMatCartesian_Refactored_New( 
		Name                                                name,
		Dimension_Index                                     dim,
		FiniteElementContext*                               context,
		Materials_Register*                                 materials_Register )
{
	ConstitutiveMatCartesian_Refactored* self = (ConstitutiveMatCartesian_Refactored*) _ConstitutiveMatCartesian_Refactored_DefaultNew( name );

	ConstitutiveMatCartesian_Refactored_InitAll( self, dim, context, materials_Register );

	return self;
}

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
ConstitutiveMatCartesian_Refactored* _ConstitutiveMatCartesian_Refactored_New( 
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
		ConstitutiveMat_Refactored_SetValueFunc*            _setValue,
		ConstitutiveMat_Refactored_GetValueFunc*            _getViscosity,
		ConstitutiveMat_Refactored_SetValueFunc*            _isotropicCorrection,
		ConstitutiveMat_Refactored_SetSecondViscosityFunc*  _setSecondViscosity,
		ConstitutiveMat_Refactored_Assemble_D_B_Func*       _assemble_D_B,
		ConstitutiveMat_Refactored_CalculateStressFunc*     _calculateStress,
		Name                                                name )
{
	ConstitutiveMatCartesian_Refactored* self;
	
	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(ConstitutiveMatCartesian_Refactored) );
	self = (ConstitutiveMatCartesian_Refactored*) _ConstitutiveMat_Refactored_New( 
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
		_setValue,
		_getViscosity,
		_isotropicCorrection,
		_setSecondViscosity,
		_assemble_D_B,
		_calculateStress,
		name );
	
	/* Function pointers for this class that are not on the parent class should be set here */
	
	return self;
}

void _ConstitutiveMatCartesian_Refactored_Init( 
		ConstitutiveMatCartesian_Refactored*                 self )
{
	self->rowSize    = StGermain_nSymmetricTensorVectorComponents( self->dim );
	self->columnSize = StGermain_nSymmetricTensorVectorComponents( self->dim );
}

void ConstitutiveMatCartesian_Refactored_InitAll( 
		void*                                        constitutiveMatrix,
		Dimension_Index                              dim,
		FiniteElementContext*                        context,
		Materials_Register*                          materials_Register )
{
	ConstitutiveMatCartesian_Refactored* self = (ConstitutiveMatCartesian_Refactored*) constitutiveMatrix;

	ConstitutiveMat_Refactored_InitAll( self, dim, context, materials_Register );
	_ConstitutiveMatCartesian_Refactored_Init( self );
}

void _ConstitutiveMatCartesian_Refactored_Delete( void* constitutiveMatrix ) {
	ConstitutiveMatCartesian_Refactored* self = (ConstitutiveMatCartesian_Refactored*)constitutiveMatrix;

	_ConstitutiveMat_Refactored_Delete( self );
}

void _ConstitutiveMatCartesian_Refactored_Print( void* constitutiveMatrix, Stream* stream ) {
	ConstitutiveMatCartesian_Refactored* self = (ConstitutiveMatCartesian_Refactored*)constitutiveMatrix;
	
	_ConstitutiveMat_Refactored_Print( self, stream );

	/* General info */
}

void* _ConstitutiveMatCartesian_Refactored_DefaultNew( Name name ) {
	return (void*)_ConstitutiveMatCartesian_Refactored_New( 
		sizeof(ConstitutiveMatCartesian_Refactored), 
		ConstitutiveMatCartesian_Refactored_Type,
		_ConstitutiveMatCartesian_Refactored_Delete,
		_ConstitutiveMatCartesian_Refactored_Print,
		NULL,
		_ConstitutiveMatCartesian_Refactored_DefaultNew,
		_ConstitutiveMatCartesian_Refactored_Construct,
		_ConstitutiveMatCartesian_Refactored_Build,
		_ConstitutiveMatCartesian_Refactored_Initialise,
		_ConstitutiveMatCartesian_Refactored_Execute,
		_ConstitutiveMatCartesian_Refactored_Destroy,
		_ConstitutiveMatCartesian_Refactored_SetValueInAllEntries,
		_ConstitutiveMatCartesian_Refactored_GetIsotropicViscosity,
		_ConstitutiveMatCartesian_Refactored_IsotropicCorrection,
		_ConstitutiveMatCartesian_Refactored_SetSecondViscosity,
		_ConstitutiveMatCartesian_Refactored_Assemble_D_B,
		_ConstitutiveMatCartesian_Refactored_CalculateStress,
		name );
}

void _ConstitutiveMatCartesian_Refactored_Construct( void* constitutiveMatrix, Stg_ComponentFactory* cf, void* data ) {
	ConstitutiveMatCartesian_Refactored*            self             = (ConstitutiveMatCartesian_Refactored*)constitutiveMatrix;

	/* Construct Parent */
	_ConstitutiveMat_Refactored_Construct( self, cf, data );

	_ConstitutiveMatCartesian_Refactored_Init( self );
}

void _ConstitutiveMatCartesian_Refactored_Build( void* constitutiveMatrix, void* data ) {
	ConstitutiveMatCartesian_Refactored*             self             = (ConstitutiveMatCartesian_Refactored*)constitutiveMatrix;

	_ConstitutiveMat_Refactored_Build( self, data );
}

void _ConstitutiveMatCartesian_Refactored_Initialise( void* constitutiveMatrix, void* data ) {
	ConstitutiveMatCartesian_Refactored*             self             = (ConstitutiveMatCartesian_Refactored*)constitutiveMatrix;

	_ConstitutiveMat_Refactored_Initialise( self, data );
}

void _ConstitutiveMatCartesian_Refactored_Execute( void* constitutiveMatrix, void* data ) {
	_ConstitutiveMat_Refactored_Execute( constitutiveMatrix, data );
}

void _ConstitutiveMatCartesian_Refactored_Destroy( void* constitutiveMatrix, void* data ) {
	_ConstitutiveMat_Refactored_Destroy( constitutiveMatrix, data );
}

void _ConstitutiveMatCartesian_Refactored_SetValueInAllEntries( void* constitutiveMatrix, double value ) {
	ConstitutiveMat_Refactored* self   = (ConstitutiveMat_Refactored*) constitutiveMatrix;

	self->_setValue = ( self->dim == 2 ? 
			_ConstitutiveMatCartesian_Refactored2D_SetValueInAllEntries :
			_ConstitutiveMatCartesian_Refactored3D_SetValueInAllEntries );

	ConstitutiveMat_Refactored_SetValueInAllEntries( self, value );
}

void _ConstitutiveMatCartesian_Refactored2D_SetValueInAllEntries( void* constitutiveMatrix, double value ) {
	ConstitutiveMat_Refactored* self   = (ConstitutiveMat_Refactored*) constitutiveMatrix;

	if ( fabs( value ) < 1.0e-20 ) 
		ConstitutiveMat_Refactored_ZeroMatrix( self );
	else {
		double**            D      = self->matrixData;

		D[0][0] = D[0][1] = D[0][2] = value;
		D[1][0] = D[1][1] = D[1][2] = value;
		D[2][0] = D[2][1] = D[2][2] = value;
	
		self->isDiagonal = False;
	}
}

void _ConstitutiveMatCartesian_Refactored3D_SetValueInAllEntries( void* _constitutiveMatrix, double value ) {
	ConstitutiveMat_Refactored* self   = (ConstitutiveMat_Refactored*)_constitutiveMatrix;

	if ( fabs( value ) < 1.0e-20 ) 
		ConstitutiveMat_Refactored_ZeroMatrix( self );
	else {
		double**            D      = self->matrixData;

		D[0][0] = D[0][1] = D[0][2] = D[0][3] = D[0][4] = D[0][5] = value;
		D[1][0] = D[1][1] = D[1][2] = D[1][3] = D[1][4] = D[1][5] = value;
		D[2][0] = D[2][1] = D[2][2] = D[2][3] = D[2][4] = D[2][5] = value;
		D[3][0] = D[3][1] = D[3][2] = D[3][3] = D[3][4] = D[3][5] = value;
		D[4][0] = D[4][1] = D[4][2] = D[4][3] = D[4][4] = D[4][5] = value;
		D[5][0] = D[5][1] = D[5][2] = D[5][3] = D[5][4] = D[5][5] = value;
	
		self->isDiagonal = False;
	}
}

double _ConstitutiveMatCartesian_Refactored_GetIsotropicViscosity( void* constitutiveMatrix ) {
	ConstitutiveMat_Refactored* self   = (ConstitutiveMat_Refactored*) constitutiveMatrix;

	self->_getViscosity = ( self->dim == 2 ? 
			_ConstitutiveMatCartesian_Refactored2D_GetIsotropicViscosity :
			_ConstitutiveMatCartesian_Refactored3D_GetIsotropicViscosity );

	return ConstitutiveMat_Refactored_GetIsotropicViscosity( self );
}

double _ConstitutiveMatCartesian_Refactored2D_GetIsotropicViscosity( void* constitutiveMatrix ) {
	ConstitutiveMat_Refactored* self = (ConstitutiveMat_Refactored*) constitutiveMatrix;

	return self->matrixData[2][2];
}

double _ConstitutiveMatCartesian_Refactored3D_GetIsotropicViscosity( void* constitutiveMatrix ) {
	ConstitutiveMat_Refactored* self = (ConstitutiveMat_Refactored*) constitutiveMatrix;

	return self->matrixData[3][3];
}

void _ConstitutiveMatCartesian_Refactored_IsotropicCorrection( void* constitutiveMatrix, double isotropicCorrection ) {
	ConstitutiveMat_Refactored* self   = (ConstitutiveMat_Refactored*) constitutiveMatrix;

	self->_isotropicCorrection = ( self->dim == 2 ? 
			_ConstitutiveMatCartesian_Refactored2D_IsotropicCorrection :
			_ConstitutiveMatCartesian_Refactored3D_IsotropicCorrection );

	ConstitutiveMat_Refactored_IsotropicCorrection( self, isotropicCorrection );
}

void _ConstitutiveMatCartesian_Refactored2D_IsotropicCorrection( void* constitutiveMatrix, double isotropicCorrection ) {
	ConstitutiveMat_Refactored* self   = (ConstitutiveMat_Refactored*) constitutiveMatrix;
	double**            D      = self->matrixData;
		
	D[0][0] += 2.0 * isotropicCorrection;
	D[1][1] += 2.0 * isotropicCorrection;
	D[2][2] += isotropicCorrection;
}

void _ConstitutiveMatCartesian_Refactored3D_IsotropicCorrection( void* constitutiveMatrix, double isotropicCorrection ) {
	ConstitutiveMat_Refactored* self   = (ConstitutiveMat_Refactored*) constitutiveMatrix;
	double**            D      = self->matrixData;

	D[0][0] += 2.0 * isotropicCorrection;
	D[1][1] += 2.0 * isotropicCorrection;
	D[2][2] += 2.0 * isotropicCorrection;
	
	D[3][3] += isotropicCorrection;
	D[4][4] += isotropicCorrection;
	D[5][5] += isotropicCorrection;
}

void _ConstitutiveMatCartesian_Refactored_SetSecondViscosity( void* constitutiveMatrix, double deltaViscosity, XYZ director ) {
	ConstitutiveMat_Refactored* self   = (ConstitutiveMat_Refactored*) constitutiveMatrix;

	self->_setSecondViscosity = ( self->dim == 2 ? 
			_ConstitutiveMatCartesian_Refactored2D_SetSecondViscosity :
			_ConstitutiveMatCartesian_Refactored3D_SetSecondViscosity );

	ConstitutiveMat_Refactored_SetSecondViscosity( self, deltaViscosity, director );
}

void _ConstitutiveMatCartesian_Refactored2D_SetSecondViscosity( void* constitutiveMatrix, double deltaViscosity, XYZ director ) {
	ConstitutiveMat_Refactored* self      = (ConstitutiveMat_Refactored*) constitutiveMatrix;
	double**            D         = self->matrixData;
	double              n1        = director[ I_AXIS ];
	double              n2        = director[ J_AXIS ];
	double              a0;
	double              a1;

	a0 = 4.0 * deltaViscosity * n1 * n1 * n2 * n2;
	a1 = 2.0 * deltaViscosity * n1 * n2 * (n2*n2 - n1*n1);

	D[0][0] += -a0 ;	D[0][1] +=  a0 ;	D[0][2] += -a1 ;
	D[1][0] +=  a0 ;	D[1][1] += -a0 ;	D[1][2] +=  a1 ;
	D[2][0] += -a1 ;	D[2][1] +=  a1 ;	D[2][2] +=  a0 - deltaViscosity ;

	self->isDiagonal = False;
}

void _ConstitutiveMatCartesian_Refactored3D_SetSecondViscosity( void* constitutiveMatrix, double deltaViscosity, XYZ director ) {
	ConstitutiveMat_Refactored* self      = (ConstitutiveMat_Refactored*) constitutiveMatrix;
	double**            D         = self->matrixData;
	double              n1        = director[ I_AXIS ];
	double              n2        = director[ J_AXIS ];
	double              n3        = director[ K_AXIS ];
	double              a00,a01,a02,a03,a04,a05;
	double                  a11,a12,a13,a14,a15;
	double                      a22,a23,a24,a25;
	double                          a33,a34,a35;
	double                              a44,a45;
	double                                  a55;	

	a00 = -4 * n1*n1 * ( 1 - n1*n1 ) * deltaViscosity; 
	a01 =  4 * n1*n1 * n2*n2 * deltaViscosity; 
	a02 =  4 * n1*n1 * n3*n3 * deltaViscosity; 
	a03 =  2 * n1*n2 * (2*n1*n1-1) * deltaViscosity;
	a04 =  2 * n1*n3 * (2*n1*n1-1) * deltaViscosity; 
	a05 =  4 * n1*n1 * n2*n3 * deltaViscosity;
		 
	a11= 4 * n2*n2 * (n2*n2-1) * deltaViscosity; 
	a12= 4 * n2*n2 * n3*n3 * deltaViscosity; 
	a13= 2 * n1*n2 * (2*n2*n2-1) * deltaViscosity; 
	a14= 4 * n1*n2 * n2*n3 * deltaViscosity;  
	a15= 2 * n2*n3 * (2*n2*n2-1) * deltaViscosity;
		
	a22 = 4 * n3*n3 * (n3*n3-1) * deltaViscosity; 
	a23 = 4 * n1*n2 * n3*n3 * deltaViscosity; 
	a24 = 2 * n1*n3 * (2*n3*n3-1) * deltaViscosity;
	a25 = 2 * n2*n3 * (2*n3*n3-1) * deltaViscosity;
	 
	a33 = (4 * n1*n1 * n2*n2 - n1*n1 - n2*n2) * deltaViscosity; 
	a34 = (4 * n1*n1 * n2*n3 - n2*n3) * deltaViscosity; 
	a35 = (4 * n1*n2 * n2*n3 - n1*n3) * deltaViscosity;
	
	a44 = (4 * n1*n1 * n3*n3 - n1*n1 -n3*n3) * deltaViscosity; 
	a45 = (4 * n1*n2 * n3*n3 - n1*n2) * deltaViscosity;

	a55 = (4 * n3*n3 * n2*n2 - n3*n3 - n2*n2) * deltaViscosity;
		
	/* D_{anisotropic} to D */
	D[0][0] += a00 ; D[0][1] += a01 ; D[0][2] += a02 ; D[0][3] += a03 ; D[0][4] += a04 ; D[0][5] += a05 ;
	D[1][0] += a01 ; D[1][1] += a11 ; D[1][2] += a12 ; D[1][3] += a13 ; D[1][4] += a14 ; D[1][5] += a15 ;
	D[2][0] += a02 ; D[2][1] += a12 ; D[2][2] += a22 ; D[2][3] += a23 ; D[2][4] += a24 ; D[2][5] += a25 ;
	D[3][0] += a03 ; D[3][1] += a13 ; D[3][2] += a23 ; D[3][3] += a33 ; D[3][4] += a34 ; D[3][5] += a35 ;
	D[4][0] += a04 ; D[4][1] += a14 ; D[4][2] += a24 ; D[4][3] += a34 ; D[4][4] += a44 ; D[4][5] += a45 ;
	D[5][0] += a05 ; D[5][1] += a15 ; D[5][2] += a25 ; D[5][3] += a35 ; D[5][4] += a45 ; D[5][5] += a55 ;

	self->isDiagonal = False;
}


void _ConstitutiveMatCartesian_Refactored_Assemble_D_B( void* constitutiveMatrix, double** GNx, Node_Index node_I, double** D_B ) {
	ConstitutiveMat_Refactored* self   = (ConstitutiveMat_Refactored*) constitutiveMatrix;

	self->_assemble_D_B = ( self->dim == 2 ? 
			_ConstitutiveMatCartesian_Refactored2D_Assemble_D_B :
			_ConstitutiveMatCartesian_Refactored3D_Assemble_D_B );

	ConstitutiveMat_Refactored_Assemble_D_B( self, GNx, node_I, D_B );
}

/*
[B] = [ d/dx,     0  ]
      [    0,  d/dy  ]
      [ d/dy,  d/dx  ]  */
void _ConstitutiveMatCartesian_Refactored2D_Assemble_D_B( void* constitutiveMatrix, double** GNx, Node_Index node_I, double** D_B ){
	ConstitutiveMat_Refactored* self = (ConstitutiveMat_Refactored*) constitutiveMatrix;
	double**            D    = self->matrixData;
	double              d_dx = GNx[ I_AXIS ][ node_I ];
	double              d_dy = GNx[ J_AXIS ][ node_I ];

	if (self->isDiagonal) {
		D_B[0][0] = D[0][0] * d_dx;
		D_B[0][1] = 0.0;
				
		D_B[1][0] = 0.0;
		D_B[1][1] = D[1][1] * d_dy;
				
		D_B[2][0] = D[2][2] * d_dy;
		D_B[2][1] = D[2][2] * d_dx;		
	}
	else {
		D_B[0][0] = D[0][0] * d_dx + D[0][2] * d_dy;
		D_B[0][1] = D[0][1] * d_dy + D[0][2] * d_dx;
				
		D_B[1][0] = D[1][0] * d_dx + D[1][2] * d_dy;
		D_B[1][1] = D[1][1] * d_dy + D[1][2] * d_dx;
				
		D_B[2][0] = D[2][0] * d_dx + D[2][2] * d_dy;
		D_B[2][1] = D[2][1] * d_dy + D[2][2] * d_dx;
	}
}


/*
[B] = [ d/dx,     0,      0  ]
      [    0,  d/dy,      0  ]
      [    0,     0,   d/dx  ]
      [ d/dy,  d/dx,      0  ] 
      [ d/dz,     0,   d/dx  ]
      [    0,  d/dz,   d/dy  ] */
void _ConstitutiveMatCartesian_Refactored3D_Assemble_D_B( void* constitutiveMatrix, double** GNx, Node_Index node_I, double** D_B ){
	ConstitutiveMat_Refactored* self = (ConstitutiveMat_Refactored*) constitutiveMatrix;
	double**            D    = self->matrixData;
	double              d_dx = GNx[ I_AXIS ][ node_I ];
	double              d_dy = GNx[ J_AXIS ][ node_I ];
	double              d_dz = GNx[ K_AXIS ][ node_I ];

	if (self->isDiagonal) {
		D_B[0][0] = D[0][0] * d_dx;
		D_B[0][1] = 0.0;
		D_B[0][2] = 0.0;
		
		D_B[1][0] = 0.0;
		D_B[1][1] = D[1][1] * d_dy;
		D_B[1][2] = 0.0;

		D_B[2][0] = 0.0;
		D_B[2][1] = 0.0;
		D_B[2][2] = D[2][2] * d_dz;

		D_B[3][0] = D[3][3] * d_dy;
		D_B[3][1] = D[3][3] * d_dx;
		D_B[3][2] = 0.0;

		D_B[4][0] = D[4][4] * d_dz;
		D_B[4][1] = 0.0;
		D_B[4][2] = D[4][4] * d_dx;
			
		D_B[5][0] = 0.0;
		D_B[5][1] = D[5][5] * d_dz;
		D_B[5][2] = D[5][5] * d_dy;
	}
	else {
		D_B[0][0] = D[0][0] * d_dx + D[0][3] * d_dy + D[0][4] * d_dz;
		D_B[0][1] = D[0][1] * d_dy + D[0][3] * d_dx + D[0][5] * d_dz;
		D_B[0][2] = D[0][2] * d_dz + D[0][4] * d_dx + D[0][5] * d_dy;
		
		D_B[1][0] = D[1][0] * d_dx + D[1][3] * d_dy + D[1][4] * d_dz;
		D_B[1][1] = D[1][1] * d_dy + D[1][3] * d_dx + D[1][5] * d_dz;
		D_B[1][2] = D[1][2] * d_dz + D[1][4] * d_dx + D[1][5] * d_dy;

		D_B[2][0] = D[2][0] * d_dx + D[2][3] * d_dy + D[2][4] * d_dz;
		D_B[2][1] = D[2][1] * d_dy + D[2][3] * d_dx + D[2][5] * d_dz;
		D_B[2][2] = D[2][2] * d_dz + D[2][4] * d_dx + D[2][5] * d_dy;

		D_B[3][0] = D[3][0] * d_dx + D[3][3] * d_dy + D[3][4] * d_dz;
		D_B[3][1] = D[3][1] * d_dy + D[3][3] * d_dx + D[3][5] * d_dz;
		D_B[3][2] = D[3][2] * d_dz + D[3][4] * d_dx + D[3][5] * d_dy;

		D_B[4][0] = D[4][0] * d_dx + D[4][3] * d_dy + D[4][4] * d_dz;
		D_B[4][1] = D[4][1] * d_dy + D[4][3] * d_dx + D[4][5] * d_dz;
		D_B[4][2] = D[4][2] * d_dz + D[4][4] * d_dx + D[4][5] * d_dy;
			
		D_B[5][0] = D[5][0] * d_dx + D[5][3] * d_dy + D[5][4] * d_dz;
		D_B[5][1] = D[5][1] * d_dy + D[5][3] * d_dx + D[5][5] * d_dz;
		D_B[5][2] = D[5][2] * d_dz + D[5][4] * d_dx + D[5][5] * d_dy;
	}
}


void _ConstitutiveMatCartesian_Refactored_CalculateStress( void* constitutiveMatrix, SymmetricTensor strainRate, SymmetricTensor stress ) {
	ConstitutiveMat_Refactored* self   = (ConstitutiveMat_Refactored*) constitutiveMatrix;

	self->_calculateStress = ( self->dim == 2 ? 
			_ConstitutiveMatCartesian_Refactored2D_CalculateStress :
			_ConstitutiveMatCartesian_Refactored3D_CalculateStress );

	ConstitutiveMat_Refactored_CalculateStress( self, strainRate, stress );
}

void _ConstitutiveMatCartesian_Refactored2D_CalculateStress( void* constitutiveMatrix, SymmetricTensor strainRate, SymmetricTensor stress ) {
	ConstitutiveMat_Refactored* self = (ConstitutiveMat_Refactored*) constitutiveMatrix;
	double**            D    = self->matrixData;

	if (self->isDiagonal) {
		stress[0] = D[0][0] * strainRate[0];
		stress[1] = D[1][1] * strainRate[1];
		stress[2] = D[2][2] * 2.0 * strainRate[2];
	}
	else {
		stress[0] = D[0][0] * strainRate[0] + D[0][1] * strainRate[1] + D[0][2] * 2.0 * strainRate[2];
		stress[1] = D[1][0] * strainRate[0] + D[1][1] * strainRate[1] + D[1][2] * 2.0 * strainRate[2];
		stress[2] = D[2][0] * strainRate[0] + D[2][1] * strainRate[1] + D[2][2] * 2.0 * strainRate[2];
	}
}



void _ConstitutiveMatCartesian_Refactored3D_CalculateStress( void* constitutiveMatrix, SymmetricTensor strainRate, SymmetricTensor stress ) {
	ConstitutiveMat_Refactored* self = (ConstitutiveMat_Refactored*) constitutiveMatrix;
	double**            D    = self->matrixData;
	          	
	if (self->isDiagonal) {
		stress[0] = D[0][0] * strainRate[0];
		stress[1] = D[1][1] * strainRate[1];
		stress[2] = D[2][2] * strainRate[2];
		stress[3] = D[3][3] * 2.0 * strainRate[3];
		stress[4] = D[4][4] * 2.0 * strainRate[4];
		stress[5] = D[5][5] * 2.0 * strainRate[5];
	}
	else {
		stress[0] = D[0][0] * strainRate[0] + D[0][1] * strainRate[1] + D[0][2] * strainRate[2] 
			+ 2.0 * (D[0][3] * strainRate[3] + D[0][4] * strainRate[4] + D[0][5] * strainRate[5]);

		stress[1] = D[1][0] * strainRate[0] + D[1][1] * strainRate[1] + D[1][2] * strainRate[2] 
			+ 2.0 * (D[1][3] * strainRate[3] + D[1][4] * strainRate[4] + D[1][5] * strainRate[5]);

		stress[2] = D[2][0] * strainRate[0] + D[2][1] * strainRate[1] + D[2][2] * strainRate[2] 
			+ 2.0 * (D[2][3] * strainRate[3] + D[2][4] * strainRate[4] + D[2][5] * strainRate[5]);

		stress[3] = D[3][0] * strainRate[0] + D[3][1] * strainRate[1] + D[3][2] * strainRate[2] 
			+ 2.0 * (D[3][3] * strainRate[3] + D[3][4] * strainRate[4] + D[3][5] * strainRate[5]);

		stress[4] = D[4][0] * strainRate[0] + D[4][1] * strainRate[1] + D[4][2] * strainRate[2] 
			+ 2.0 * (D[4][3] * strainRate[3] + D[4][4] * strainRate[4] + D[4][5] * strainRate[5]);

		stress[5] = D[5][0] * strainRate[0] + D[5][1] * strainRate[1] + D[5][2] * strainRate[2] 
			+ 2.0 * (D[5][3] * strainRate[3] + D[5][4] * strainRate[4] + D[5][5] * strainRate[5]);
	}
}

