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
*+		Mirko Velic
*+		Julian Giordani
*+		Robert Turnbull
*+		Vincent Lemiale
*+		Louis Moresi
*+		David May
*+		David Stegman
*+		Patrick Sunter
** $Id: $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#ifndef __Underworld_Velic_solKz_h__
#define __Underworld_Velic_solKz_h__

	extern const Type Velic_solKz_Type;

	typedef struct {
		__FieldTest
		double sigma;
		double km;
		int n;
		double B;
	} Velic_solKz;

	Index Underworld_Velic_solKz_Register( PluginsManager* pluginsManager );
	void* _Velic_solKz_DefaultNew( Name name );
	void _Velic_solKz_AssignFromXML( void* analyticSolution, Stg_ComponentFactory* cf, void* data );
	void _Velic_solKz_Init( Velic_solKz* self, double sigma, double km, double B, int n );
	void Velic_solKz_PressureFunction( void* analyticSolution, const double *coord, double* pressure );
	void Velic_solKz_VelocityFunction( void* analyticSolution, const double *coord, double* velocity );
	void Velic_solKz_StressFunction( void* analyticSolution, const double *coord, double* stress );
	void Velic_solKz_StrainRateFunction( void* analyticSolution, const double *coord, double* strainRate );

	void _Velic_solKz( 
		const double pos[],
		double _sigma, /* density */
		double _km, int _n, /* wavelength in z, wavenumber in x */
		double _B, /* viscosity parameter */
		double vel[], double* presssure, 
		double total_stress[], double strain_rate[] );
#endif
