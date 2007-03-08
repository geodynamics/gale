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
** $Id: LateralViscosityAnalytic.c 358 2006-10-18 06:17:30Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Rheology/Rheology.h>

#include <assert.h>

#define SMALL 1.0e-5
#define IS_ODD(A)   ( (A) % 2 == 1 )

const Type LateralViscosityAnalytic_Type = "LateralViscosityAnalytic";

typedef struct {
	__AnalyticSolution
	double                  beta;
	int                     wavenumberX;
	int                     wavenumberY;
} LateralViscosityAnalytic;

/** Analytic Solution taken from 
 *  Shijie Zhong. Analytic solutions for Stokes' flow with lateral variations in viscosity. Geophys. J. Int., 124:18-28, 1996.
 *  All equations refer to this paper */

void LateralViscosityAnalytic_TemperatureIC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DiscretisationContext*  context            = (DiscretisationContext*)_context;
	FeVariable*             temperatureField   = (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	FiniteElement_Mesh*     mesh               = temperatureField->feMesh;
	BlockGeometry*          geometry           = (BlockGeometry*) mesh->layout->elementLayout->geometry;
	double*                 result             = (double*) _result;
	Dictionary*             dictionary         = context->dictionary;
	double*                 coord;
	double                  x; 
	double                  y;
	double                  kx;
	double                  ky;
	int                     wavenumberX;
	int                     wavenumberY;
	double                  L;
	
	/* Find coordinate of node */
	coord = Mesh_CoordAt( mesh, node_lI );

	/* Make sure that the box has right dimensions */
	assert( ( geometry->max[ J_AXIS ] - geometry->min[ J_AXIS ] - 1.0 ) < SMALL );
	L = geometry->max[ I_AXIS ] - geometry->min[ I_AXIS ];

	x = coord[ I_AXIS ] - geometry->min[ I_AXIS ];
	y = coord[ J_AXIS ] - geometry->min[ J_AXIS ];

	wavenumberX = Dictionary_GetInt_WithDefault( dictionary, "wavenumberX", 1 );
	wavenumberY = Dictionary_GetInt_WithDefault( dictionary, "wavenumberY", 1 );

	kx = wavenumberX * M_PI/ L;
	ky = wavenumberY * M_PI;

	*result = sin( ky * y ) * cos( kx * x );
}

void _LateralViscosityAnalytic_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ) {
	LateralViscosityAnalytic*         self               = (LateralViscosityAnalytic*)analyticSolution;
	BlockGeometry*          geometry           = (BlockGeometry*) analyticFeVariable->feMesh->layout->elementLayout->geometry;
	Cmplx                   N                  = {0,0};
	double                  chi1;
	double                  chi2;
	Cmplx                   sumX;
	Cmplx                   sumY;
	Cmplx                   m[5];                       /* CRAPPY FORTRAN ARRAYS */
	Cmplx                   mSquared[5];
	Cmplx                   mExp[5];
	Cmplx                   N1;
	Cmplx                   N2;
	Cmplx                   N3;
	Cmplx                   N4;
	double                  beta2_4kn;
	Cmplx                   lambda1;
	Cmplx                   lambda2;
	Cmplx                   n1;
	Cmplx                   n2;
	Cmplx                   n1_minus_beta;
	Cmplx                   n2_minus_beta;
	Cmplx                   M;
	Cmplx                   oneOnM;
	Cmplx                   tmp;
	double                  beta;
	double                  kn;
	Index                   index;
	Index                   i;
	Index                   j;
	Index                   h;
	Index                   k;
	double                  L;
	double                  kx;
	double                  ky;
	Cmplx                   A1;
	Cmplx                   A2;
	Cmplx                   A3;
	Cmplx                   A4;
	Cmplx                   exp_n1x;
	Cmplx                   exp_n2x;
	Cmplx                   a1End;
	Cmplx                   a2End;
	double                  x;
	double                  y;
	double                  factorA;
	double                  factorB;
	int                     wavenumberX;
	int                     wavenumberY;
	double                  sqRootValue;

	x = coord[ I_AXIS ] - geometry->min[ I_AXIS ];
	y = coord[ J_AXIS ] - geometry->min[ J_AXIS ];

	beta = self->beta;
	wavenumberX = self->wavenumberX;
	wavenumberY = self->wavenumberY;
	
	L    = geometry->max[ I_AXIS ] - geometry->min[ I_AXIS ];
	kx   = wavenumberX * M_PI/ L;
	ky   = wavenumberY * M_PI;
	kn   = ky; 

	beta2_4kn = beta * beta + 4.0 * kn * kn;
	sqRootValue = beta2_4kn * beta2_4kn + 16.0 * beta * beta * kn * kn;
	sqRootValue = sqrt( sqRootValue );
	chi1 = sqrt( 0.5 * (  beta2_4kn + sqRootValue ) );
	chi2 = sqrt( 0.5 * ( -beta2_4kn + sqRootValue ) );

	n1[ REAL_PART ] = 0.5 * (beta + chi1);
	n1[ IMAG_PART ] = 0.5 * (chi2);
	
	n2[ REAL_PART ] = 0.5 * (beta - chi1);
	n2[ IMAG_PART ] = 0.5 * (-chi2);

	Cmplx_Copy( n1, m[1] );
	Cmplx_Conjugate( n1, m[2] );
	Cmplx_Copy( n2, m[3] );
	Cmplx_Conjugate( n2, m[4] );
			
	for ( index = 1 ; index <= 4 ; index++ ) {
		Cmplx_Multiply( m[index], m[index], mSquared[index] );
		Cmplx_RealMultiply( m[ index ], -L, tmp );
		Cmplx_Exp( tmp, mExp[ index ] );
	}

	/* Calculate M - Equation A7 */
	Cmplx_Subtract( m[2], m[1], oneOnM );
	Cmplx_Subtract( m[3], m[1], tmp );
	Cmplx_Multiply( oneOnM, tmp, oneOnM );
	Cmplx_Subtract( m[4], m[1], tmp );
	Cmplx_Multiply( oneOnM, tmp, oneOnM );
	Cmplx_RealDivideByCmplx( oneOnM, 1.0, M );

	/* Calculate N - Equation A8 */
	for ( j = 1 ; j <= 3 ; j++ ) {
		for ( i = j+1 ; i <= 4 ; i++ ) {
			Cmplx product = {0,0};
			
			/* Get h */
			for ( h = 1 ; h <= 4 ; h++ ) {
				if ( h != i && h != j )
					break; /* Found h */
			}

			/* Get k */
			for ( k = h+1 ; k <= 4 ; k++ ) {
				if ( k != i && k != j )
					break; /* Found k */
			}			

			product[ REAL_PART ] = IS_ODD( i + j ) ? 1.0 : -1.0;
			
			Cmplx_Multiply( product, mExp[ i ], product );
			Cmplx_Multiply( product, mExp[ j ], product );

			Cmplx_Subtract( mSquared[ k ], mSquared[ h ], tmp );
			Cmplx_Multiply( product, tmp, product );

			Cmplx_Subtract( mSquared[ i ], mSquared[ j ], tmp );
			Cmplx_Multiply( product, tmp, product );

			Cmplx_Add( N, product, N );
		}
	}

	/* Build N1 - Equation A9 */
	Cmplx_Subtract( mSquared[4], mSquared[3], tmp );
	Cmplx_Multiply( tmp, mSquared[2], tmp );
	Cmplx_Multiply( tmp, mExp[2], N1 ); 

	Cmplx_Subtract( mSquared[4], mSquared[2], tmp );
	Cmplx_Multiply( tmp, mSquared[3], tmp );
	Cmplx_Multiply( tmp, mExp[3], tmp );
	Cmplx_Subtract( N1, tmp, N1 );
	
	Cmplx_Subtract( mSquared[3], mSquared[2], tmp );
	Cmplx_Multiply( tmp, mSquared[4], tmp );
	Cmplx_Multiply( tmp, mExp[4], tmp );
	Cmplx_Add( N1, tmp, N1 );
	
	/* Build N2 - Equation A10 */
	Cmplx_Subtract( mSquared[4], mSquared[2], tmp );
	Cmplx_Multiply( tmp, mSquared[1], tmp );
	Cmplx_Multiply( tmp, mExp[1], N2 );

	Cmplx_Subtract( mSquared[4], mSquared[1], tmp );
	Cmplx_Multiply( tmp, mSquared[2], tmp );
	Cmplx_Multiply( tmp, mExp[2], tmp );
	Cmplx_Subtract( N2, tmp, N2 );

	Cmplx_Subtract( mSquared[2], mSquared[1], tmp );
	Cmplx_Multiply( tmp, mSquared[4], tmp );
	Cmplx_Multiply( tmp, mExp[4], tmp );
	Cmplx_Add( N2, tmp, N2 );

	/* Build N3 - Equation A11 */
	Cmplx_Subtract( mSquared[3], mSquared[4], tmp );
	Cmplx_Multiply( tmp, mExp[2], N3 );

	Cmplx_Subtract( mSquared[4], mSquared[2], tmp );
	Cmplx_Multiply( tmp, mExp[3], tmp );
	Cmplx_Add( N3, tmp, N3 );

	Cmplx_Subtract( mSquared[2], mSquared[3], tmp );
	Cmplx_Multiply( tmp, mExp[4], tmp );
	Cmplx_Add( N3, tmp, N3 );

	/* Build N4 - Equation A12 */
	Cmplx_Subtract( mSquared[2], mSquared[4], tmp );
	Cmplx_Multiply( tmp, mExp[1], N4 );

	Cmplx_Subtract( mSquared[4], mSquared[1], tmp );
	Cmplx_Multiply( tmp, mExp[2], tmp );
	Cmplx_Add( N4, tmp, N4 );

	Cmplx_Subtract( mSquared[1], mSquared[2], tmp );
	Cmplx_Multiply( tmp, mExp[4], tmp );
	Cmplx_Add( N4, tmp, N4 );


    /* Build lambda1 - Equation B8 */
	Cmplx_AddReal( n1, - beta, n1_minus_beta );
    Cmplx_RealMultiply( n1_minus_beta,L, tmp );
	Cmplx_Exp( tmp, tmp );
    Cmplx_RealMultiply(tmp, cos( kx * L ), tmp );
	Cmplx_RealMinusCmplx( tmp, 1.0, tmp );
    Cmplx_RealMultiply(tmp, -kx, tmp );
	Cmplx_Multiply( M, tmp, lambda1 );
	
	Cmplx_Multiply( n1_minus_beta, n1_minus_beta, tmp );
	Cmplx_AddReal( tmp, kx * kx, tmp );

	Cmplx_Division( lambda1, tmp, lambda1 );
	
	/* Build lambda2 - Equation B9 */
	Cmplx_AddReal( n2, - beta, n2_minus_beta );
	Cmplx_RealMultiply( n2_minus_beta,L, tmp );
	Cmplx_Exp( tmp, tmp );
	Cmplx_RealMultiply(tmp, cos( kx * L ), tmp );
	Cmplx_RealMinusCmplx( tmp, 1.0, tmp );
	Cmplx_RealMultiply(tmp, kx, tmp );
	Cmplx_Multiply( M, tmp, lambda2 );
	
	Cmplx_Multiply( n2_minus_beta, n2_minus_beta, tmp );
	Cmplx_AddReal( tmp, kx * kx, tmp );

	Cmplx_Division( lambda2, tmp, lambda2 );



	/* Build A1 - Equation B4 */
	Cmplx_Multiply( mExp[1], lambda1, tmp );
	factorA = tmp[ REAL_PART ];
	Cmplx_Multiply( mExp[3], lambda2, tmp );
	factorA += tmp[ REAL_PART ];
	Cmplx_Division( N1, N, tmp );
	Cmplx_RealMultiply( tmp, 2.0 * factorA, A1 );

	Cmplx_Multiply( mExp[1], lambda1, tmp );
	Cmplx_Multiply( tmp, mSquared[1], tmp );
	factorB = tmp[ REAL_PART ];
	Cmplx_Multiply( mExp[3], lambda2, tmp );
	Cmplx_Multiply( tmp, mSquared[3], tmp );
	factorB += tmp[ REAL_PART ];
	Cmplx_Division( N3, N, tmp );
	Cmplx_RealMultiply( tmp, 2.0 * factorB, tmp );
	Cmplx_Add( A1, tmp, A1 );

	Cmplx_Multiply( n1_minus_beta, n1_minus_beta, a1End );
	Cmplx_AddReal( a1End, kx*kx, a1End );
	Cmplx_RealDivideByCmplx( a1End, kx, a1End );
	Cmplx_Multiply( a1End, M, a1End );
	Cmplx_Add( A1, a1End, A1 );

	/* Build A2 - Equation B5 */
	Cmplx_Division( N2, N, tmp );
	Cmplx_RealMultiply( tmp, 2.0 * factorA, A2 );

	Cmplx_Division( N4, N, tmp );
	Cmplx_RealMultiply( tmp, 2.0 * factorB, tmp );
	Cmplx_Add( A2, tmp, A2 );

	Cmplx_Multiply( n2_minus_beta, n2_minus_beta, a2End );
	Cmplx_AddReal( a2End, kx*kx, a2End );
	Cmplx_RealDivideByCmplx( a2End, kx, a2End );
	Cmplx_Multiply( a2End, M, a2End );
	Cmplx_Subtract( A2, a2End, A2 );	

	/* Build A3 - Equation B6 */
	Cmplx_Subtract( a2End, a1End, A3 );

	/* Build A4 - Equation B7 */
	Cmplx_Multiply( n1_minus_beta, n1_minus_beta, tmp );
	Cmplx_AddReal( tmp, kx*kx, tmp );
	Cmplx_Division( n1_minus_beta, tmp, tmp );
	Cmplx_Multiply( tmp, M, A4 );

	Cmplx_Multiply( n2_minus_beta, n2_minus_beta, tmp );
	Cmplx_AddReal( tmp, kx*kx, tmp );
	Cmplx_Division( n2_minus_beta, tmp, tmp );
	Cmplx_Multiply( tmp, M, tmp );
	Cmplx_Subtract( A4, tmp, A4 );
	
	/*************************************** Get Velocities - Eqn B10 - B11 *****************************************/
	
	/* Get First Term in Sums */
	Cmplx_RealMultiply( n1, -x, tmp );
	Cmplx_Exp( tmp, exp_n1x );
	Cmplx_Multiply( A1, exp_n1x, sumX );

	Cmplx_Copy( sumX, sumY );
	Cmplx_Multiply( sumY, n1, sumY );
	Cmplx_RealMultiply( sumY, -1.0, sumY );

	/* Add second term in sum */
	Cmplx_RealMultiply( n2, -x, tmp );
	Cmplx_Exp( tmp, exp_n2x );
	Cmplx_Multiply( A2, exp_n2x, tmp );
	Cmplx_Add( sumX, tmp, sumX );

	Cmplx_Multiply( tmp, n2, tmp );
	Cmplx_RealMultiply( tmp, -1.0, tmp );
	Cmplx_Add( sumY, tmp, sumY );
	
	/* Add third term in sum */
	sumX[ REAL_PART ] += (A3[ REAL_PART ] * cos( kx*x ) + A4[ REAL_PART ] * sin( kx * x ) ) * exp( -beta*x );
	sumY[ REAL_PART ] += exp( -beta * x ) * (
		(kx * A4[REAL_PART] - beta * A3[ REAL_PART ]) * cos( kx * x ) - 
		(kx * A3[REAL_PART] + beta * A4[ REAL_PART ]) * sin( kx * x ) );

	velocity[ I_AXIS ] = -2.0 * kx * ky * cos( ky * y ) * sumX[ REAL_PART ]; 
	velocity[ J_AXIS ] = 2.0 * kx * sin( ky * y ) * sumY[ REAL_PART ]; 

}

void _LateralViscosityAnalytic_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	LateralViscosityAnalytic*         self           = (LateralViscosityAnalytic*)analyticSolution;
	FeVariable*             velocityField;
	AbstractContext*        context;
	ConditionFunction*      condFunc;
	
	/* Construct Parent */
	_AnalyticSolution_Construct( self, cf, data );

	context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data ); 
	
	/* Add temperature initial condition */
	condFunc = ConditionFunction_New( LateralViscosityAnalytic_TemperatureIC, "LateralViscosityAnalytic_TemperatureIC" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );
	
	/* Create Analytic Fields */
	velocityField = Stg_ComponentFactory_ConstructByName( cf, "VelocityField", FeVariable, True, data ); 
	AnalyticSolution_CreateAnalyticField( self, velocityField, _LateralViscosityAnalytic_VelocityFunction );

	self->beta = Stg_ComponentFactory_GetRootDictDouble( cf, "beta", 0.0 );
	self->wavenumberX = Stg_ComponentFactory_GetRootDictInt( cf, "wavenumberX", 1 );
	self->wavenumberY = Stg_ComponentFactory_GetRootDictInt( cf, "wavenumberY", 1 );
}


void* _LateralViscosityAnalytic_DefaultNew( Name name ) {
	return _AnalyticSolution_New(
			sizeof(LateralViscosityAnalytic),
			LateralViscosityAnalytic_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_LateralViscosityAnalytic_DefaultNew,
			_LateralViscosityAnalytic_Construct,
			_AnalyticSolution_Build,
			_AnalyticSolution_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );
}

Index Underworld_LateralViscosityAnalytic_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, LateralViscosityAnalytic_Type, "0", _LateralViscosityAnalytic_DefaultNew );
}
