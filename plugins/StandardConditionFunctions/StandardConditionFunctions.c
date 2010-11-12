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
** $Id: StandardConditionFunctions.c 1196 2008-08-04 16:29:30Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <string.h>
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <assert.h>
#include "StandardConditionFunctions.h"

const Type StgFEM_StandardConditionFunctions_Type = "StgFEM_StandardConditionFunctions";

void _StgFEM_StandardConditionFunctions_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
	Codelet*		self		= (Codelet*)component;
	AbstractContext*        context;
	ConditionFunction*      condFunc;
	Dictionary*		pluginDict	= Codelet_GetPluginDictionary( component, cf->rootDict );

	context = (AbstractContext*)Stg_ComponentFactory_ConstructByName( cf, Dictionary_GetString( pluginDict, (Dictionary_Entry_Key)"Context"  ), AbstractContext, True, data );
	self->context = context;
	
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SolidBodyRotation, (Name)"Velocity_SolidBodyRotation"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_PartialRotationX, (Name)"Velocity_PartialRotationX"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
		
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_PartialRotationY, (Name)"Velocity_PartialRotationY"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_TaperedRotationX, (Name)"TaperedRotationX" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
		
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_TaperedRotationY, (Name)"TaperedRotationY" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SimpleShear, (Name)"Velocity_SimpleShear" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
        condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SimpleShearInverted, (Name)"Velocity_SimpleShearInverted" );
        ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_ShearZ, (Name)"ShearZ" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_Extension, (Name)"Velocity_Extension" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_PartialLid_TopLayer, (Name)"Velocity_PartialLid_TopLayer"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_Trigonometry, (Name)"Temperature_Trigonometry"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_LinearInterpolationLid, (Name)"Velocity_LinearInterpolationLid"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_Lid_RampWithCentralMax, (Name)"Velocity_Lid_RampWithCentralMax"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_LinearVelocityLeftWall, (Name)"LinearVelocityLeftWall"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_LinearVelocityRightWall, (Name)"LinearVelocityRightWall"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SinusoidalLid, (Name)"Velocity_SinusoidalLid"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_CornerOnly, (Name)"Velocity_Lid_CornerOnly"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_TemperatureCosineHill, (Name)"Temperature_CosineHill"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_ConvectionBenchmark, (Name)"Temperature_ConvectionBenchmark"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_LinearWithSinusoidalPerturbation, (Name)"LinearWithSinusoidalPerturbation"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_EdgeDriveConvectionIC, (Name)"EdgeDriveConvectionIC"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_ThermalEdgeDriveConvectionIC, (Name)"ThermalEdgeDriveConvectionIC"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_AnalyticalTemperatureIC, (Name)"AnalyticalTemperatureIC"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( Stg_FEM_VelicTemperatureIC, (Name)"VelicTemperatureIC" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( Stg_FEM_VelicTemperatureIC_SolB, (Name)"VelicTemperatureIC_SolB" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SinusoidalExtension, (Name)"SinusoidalExtension" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_StepFunction, (Name)"StepFunction" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StG_FEM_StandardConditionFunctions_StepFunctionProduct1, (Name)"StepFunctionProduct1");
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StG_FEM_StandardConditionFunctions_StepFunctionProduct2, (Name)"StepFunctionProduct2");
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StG_FEM_StandardConditionFunctions_StepFunctionProduct3, (Name)"StepFunctionProduct3");
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StG_FEM_StandardConditionFunctions_StepFunctionProduct4, (Name)"StepFunctionProduct4");
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_TemperatureProfile, (Name)"TemperatureProfile");
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StG_FEM_StandardConditionFunctions_Gaussian, (Name)"Gaussian");
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_ERF,
                                         (Name)"ERF");
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_ERFC,
                                         (Name)"ERFC");
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_RubberSheet,
                                         (Name)"RubberSheet");
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_MovingStepFunction, (Name)"MovingStepFunction");
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SpecRidge3D, (Name)"SpecRidge3D" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SpectralBCX, (Name)"SpectralBCX" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SpectralBCY, (Name)"SpectralBCY" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SpectralBCZ, (Name)"SpectralBCZ" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SpectralPressureBCX, (Name)"SpectralPressureBCX" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SpectralPressureBCY, (Name)"SpectralPressureBCY" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_ErrorFunc, (Name)"ErrorFunc" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_ConstantVector, (Name)"ConstantVector" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_GaussianDistribution, (Name)"GaussianDistribution" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_1DGaussianDistribution, (Name)"1DGaussianDistribution" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_HalfContainer, (Name)"HalfContainer" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_ConstantValue, (Name)"ConstantValue" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_DiagonalLine, (Name)"DiagonalLine" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_DeltaFunction, (Name)"DeltaFunction" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_InflowBottom, (Name)"InflowBottom" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

        condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_GaussianTube, (Name)"GaussianTube" );
        ConditionFunction_Register_Add( condFunc_Register, condFunc );

        condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_GravitationalPotential, (Name)"GravitationalPotential" );
        ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_WarsTemperature,
                                         (Name)"WarsTemperature");
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_Quadratic,
                                         (Name)"Quadratic");
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

	condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File1,
                                         (Name)"File1");
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File2,
                                         (Name)"File2");
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File3,
                                         (Name)"File3");
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File4,
                                         (Name)"File4");
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File5,
                                         (Name)"File5");
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File6,
                                         (Name)"File6");
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File7,
                                         (Name)"File7");
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File8,
                                         (Name)"File8");
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File9,
                                         (Name)"File9");
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File10,
                                         (Name)"File10");
	ConditionFunction_Register_Add( condFunc_Register, condFunc );

}

void _StgFEM_StandardConditionFunctions_Destroy( void* _self, void* data ) {
   /* This function will totally clean the condFunc_Register
    *
    * This could be trouble some if other code adds or deletes condition functions on this register
    */

   unsigned *refCount = &(condFunc_Register->count);

   /* first check if there are things still on the condFunc_Register, if so .... */
   if( *refCount != 0 ) {
      while( *refCount != 0 ) {

         _ConditionFunction_Delete( condFunc_Register->_cf[ *refCount-1 ] );
         condFunc_Register->_cf[ *refCount-1 ] = NULL;

         *refCount = *refCount - 1;
      }
   }
   _Codelet_Destroy( _self, data );
}
void* _StgFEM_StandardConditionFunctions_DefaultNew( Name name ) {
	return Codelet_New(
		StgFEM_StandardConditionFunctions_Type,
		_StgFEM_StandardConditionFunctions_DefaultNew,
		_StgFEM_StandardConditionFunctions_AssignFromXML,
		_Codelet_Build,
		_Codelet_Initialise,
		_Codelet_Execute,
		_StgFEM_StandardConditionFunctions_Destroy,
		name );
}

Index StgFEM_StandardConditionFunctions_Register( PluginsManager* pluginsManager ) {
	Journal_DPrintf( StgFEM_Debug, "In: %s( void* )\n", __func__ );

	return PluginsManager_Submit( pluginsManager, StgFEM_StandardConditionFunctions_Type, (Name)"0", _StgFEM_StandardConditionFunctions_DefaultNew  );
}

Bool StgFEM_StandardConditionFunctions_Init( int* argc, char** argv[] ) {
  Stg_ComponentRegister* componentsRegister = Stg_ComponentRegister_Get_ComponentRegister();
  Stg_ComponentRegister_Add(componentsRegister,
                            StgFEM_StandardConditionFunctions_Type, (Name)"0",
                            _StgFEM_StandardConditionFunctions_DefaultNew );
  RegisterParent( StgFEM_StandardConditionFunctions_Type, Stg_Component_Type );
}

#ifdef NO_ERF

/* Copied from the OpenBSD iplementation of erf.c
   (src/lib/libm/src/erf.c and src/lib/libm/src/math_private.h).
   Modified to only work on 32 bit little endian machines.
   This is just a hack for Windows machines. */

/* @(#)s_erf.c 5.1 93/09/24 */
/*
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunPro, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice 
 * is preserved.
 * ====================================================
 */

/* double erf(double x)
 * double erfc(double x)
 *                           x
 *                    2      |\
 *     erf(x)  =  ---------  | exp(-t*t)dt
 *                  sqrt(pi) \| 
 *                           0
 *
 *     erfc(x) =  1-erf(x)
 *  Note that 
 *              erf(-x) = -erf(x)
 *              erfc(-x) = 2 - erfc(x)
 *
 * Method:
 *      1. For |x| in [0, 0.84375]
 *          erf(x)  = x + x*R(x^2)
 *          erfc(x) = 1 - erf(x)           if x in [-.84375,0.25]
 *                  = 0.5 + ((0.5-x)-x*R)  if x in [0.25,0.84375]
 *         where R = P/Q where P is an odd poly of degree 8 and
 *         Q is an odd poly of degree 10.
 *                                               -57.90
 *                      | R - (erf(x)-x)/x | <= 2
 *      
 *
 *         Remark. The formula is derived by noting
 *          erf(x) = (2/sqrt(pi))*(x - x^3/3 + x^5/10 - x^7/42 + ....)
 *         and that
 *          2/sqrt(pi) = 1.128379167095512573896158903121545171688
 *         is close to one. The interval is chosen because the fix
 *         point of erf(x) is near 0.6174 (i.e., erf(x)=x when x is
 *         near 0.6174), and by some experiment, 0.84375 is chosen to
 *         guarantee the error is less than one ulp for erf.
 *
 *      2. For |x| in [0.84375,1.25], let s = |x| - 1, and
 *         c = 0.84506291151 rounded to single (24 bits)
 *              erf(x)  = sign(x) * (c  + P1(s)/Q1(s))
 *              erfc(x) = (1-c)  - P1(s)/Q1(s) if x > 0
 *                        1+(c+P1(s)/Q1(s))    if x < 0
 *              |P1/Q1 - (erf(|x|)-c)| <= 2**-59.06
 *         Remark: here we use the taylor series expansion at x=1.
 *              erf(1+s) = erf(1) + s*Poly(s)
 *                       = 0.845.. + P1(s)/Q1(s)
 *         That is, we use rational approximation to approximate
 *                      erf(1+s) - (c = (single)0.84506291151)
 *         Note that |P1/Q1|< 0.078 for x in [0.84375,1.25]
 *         where 
 *              P1(s) = degree 6 poly in s
 *              Q1(s) = degree 6 poly in s
 *
 *      3. For x in [1.25,1/0.35(~2.857143)], 
 *              erfc(x) = (1/x)*exp(-x*x-0.5625+R1/S1)
 *              erf(x)  = 1 - erfc(x)
 *         where 
 *              R1(z) = degree 7 poly in z, (z=1/x^2)
 *              S1(z) = degree 8 poly in z
 *
 *      4. For x in [1/0.35,28]
 *              erfc(x) = (1/x)*exp(-x*x-0.5625+R2/S2) if x > 0
 *                      = 2.0 - (1/x)*exp(-x*x-0.5625+R2/S2) if -6<x<0
 *                      = 2.0 - tiny                (if x <= -6)
 *              erf(x)  = sign(x)*(1.0 - erfc(x)) if x < 6, else
 *              erf(x)  = sign(x)*(1.0 - tiny)
 *         where
 *              R2(z) = degree 6 poly in z, (z=1/x^2)
 *              S2(z) = degree 7 poly in z
 *
 *      Note1:
 *         To compute exp(-x*x-0.5625+R/S), let s be a single
 *         precision number and s := x; then
 *              -x*x = -s*s + (s-x)*(s+x)
 *              exp(-x*x-0.5626+R/S) = 
 *                      exp(-s*s-0.5625)*exp((s-x)*(s+x)+R/S);
 *      Note2:
 *         Here 4 and 5 make use of the asymptotic series
 *                        exp(-x*x)
 *              erfc(x) ~ ---------- * ( 1 + Poly(1/x^2) )
 *                        x*sqrt(pi)
 *         We use rational approximation to approximate
 *              g(s)=f(1/x^2) = log(erfc(x)*x) - x*x + 0.5625
 *         Here is the error bound for R1/S1 and R2/S2
 *              |R1/S1 - f(x)|  < 2**(-62.57)
 *              |R2/S2 - f(x)|  < 2**(-61.52)
 *
 *      5. For inf > x >= 28
 *              erf(x)  = sign(x) *(1 - tiny)  (raise inexact)
 *              erfc(x) = tiny*tiny (raise underflow) if x > 0
 *                      = 2 - tiny if x<0
 *
 *      7. Special case:
 *              erf(0)  = 0, erf(inf)  = 1, erf(-inf) = -1,
 *              erfc(0) = 1, erfc(inf) = 0, erfc(-inf) = 2, 
 *                 erfc/erf(NaN) is NaN
 */

/*  Assume little endian, 32 bit machines  */

typedef int int32_t;
typedef unsigned int u_int32_t;

typedef union
{
  double value;
  struct
  {
    u_int32_t lsw;
    u_int32_t msw;
  } parts;
} ieee_double_shape_type;

/* Get the more significant 32 bit int from a double.  */

#define GET_HIGH_WORD(i,d)                                      \
do {                                                            \
  ieee_double_shape_type gh_u;                                  \
  gh_u.value = (d);                                             \
  (i) = gh_u.parts.msw;                                         \
} while (0)

/* Set the less significant 32 bits of a double from an int.  */

#define SET_LOW_WORD(d,v)                                       \
do {                                                            \
  ieee_double_shape_type sl_u;                                  \
  sl_u.value = (d);                                             \
  sl_u.parts.lsw = (v);                                         \
  (d) = sl_u.value;                                             \
} while (0)


static const double
tiny        = 1e-300,
half=  5.00000000000000000000e-01, /* 0x3FE00000, 0x00000000 */
one =  1.00000000000000000000e+00, /* 0x3FF00000, 0x00000000 */
two =  2.00000000000000000000e+00, /* 0x40000000, 0x00000000 */
        /* c = (float)0.84506291151 */
erx =  8.45062911510467529297e-01, /* 0x3FEB0AC1, 0x60000000 */
/*
 * Coefficients for approximation to  erf on [0,0.84375]
 */
efx =  1.28379167095512586316e-01, /* 0x3FC06EBA, 0x8214DB69 */
efx8=  1.02703333676410069053e+00, /* 0x3FF06EBA, 0x8214DB69 */
pp0  =  1.28379167095512558561e-01, /* 0x3FC06EBA, 0x8214DB68 */
pp1  = -3.25042107247001499370e-01, /* 0xBFD4CD7D, 0x691CB913 */
pp2  = -2.84817495755985104766e-02, /* 0xBF9D2A51, 0xDBD7194F */
pp3  = -5.77027029648944159157e-03, /* 0xBF77A291, 0x236668E4 */
pp4  = -2.37630166566501626084e-05, /* 0xBEF8EAD6, 0x120016AC */
qq1  =  3.97917223959155352819e-01, /* 0x3FD97779, 0xCDDADC09 */
qq2  =  6.50222499887672944485e-02, /* 0x3FB0A54C, 0x5536CEBA */
qq3  =  5.08130628187576562776e-03, /* 0x3F74D022, 0xC4D36B0F */
qq4  =  1.32494738004321644526e-04, /* 0x3F215DC9, 0x221C1A10 */
qq5  = -3.96022827877536812320e-06, /* 0xBED09C43, 0x42A26120 */
/*
 * Coefficients for approximation to  erf  in [0.84375,1.25] 
 */
pa0  = -2.36211856075265944077e-03, /* 0xBF6359B8, 0xBEF77538 */
pa1  =  4.14856118683748331666e-01, /* 0x3FDA8D00, 0xAD92B34D */
pa2  = -3.72207876035701323847e-01, /* 0xBFD7D240, 0xFBB8C3F1 */
pa3  =  3.18346619901161753674e-01, /* 0x3FD45FCA, 0x805120E4 */
pa4  = -1.10894694282396677476e-01, /* 0xBFBC6398, 0x3D3E28EC */
pa5  =  3.54783043256182359371e-02, /* 0x3FA22A36, 0x599795EB */
pa6  = -2.16637559486879084300e-03, /* 0xBF61BF38, 0x0A96073F */
qa1  =  1.06420880400844228286e-01, /* 0x3FBB3E66, 0x18EEE323 */
qa2  =  5.40397917702171048937e-01, /* 0x3FE14AF0, 0x92EB6F33 */
qa3  =  7.18286544141962662868e-02, /* 0x3FB2635C, 0xD99FE9A7 */
qa4  =  1.26171219808761642112e-01, /* 0x3FC02660, 0xE763351F */
qa5  =  1.36370839120290507362e-02, /* 0x3F8BEDC2, 0x6B51DD1C */
qa6  =  1.19844998467991074170e-02, /* 0x3F888B54, 0x5735151D */
/*
 * Coefficients for approximation to  erfc in [1.25,1/0.35]
 */
ra0  = -9.86494403484714822705e-03, /* 0xBF843412, 0x600D6435 */
ra1  = -6.93858572707181764372e-01, /* 0xBFE63416, 0xE4BA7360 */
ra2  = -1.05586262253232909814e+01, /* 0xC0251E04, 0x41B0E726 */
ra3  = -6.23753324503260060396e+01, /* 0xC04F300A, 0xE4CBA38D */
ra4  = -1.62396669462573470355e+02, /* 0xC0644CB1, 0x84282266 */
ra5  = -1.84605092906711035994e+02, /* 0xC067135C, 0xEBCCABB2 */
ra6  = -8.12874355063065934246e+01, /* 0xC0545265, 0x57E4D2F2 */
ra7  = -9.81432934416914548592e+00, /* 0xC023A0EF, 0xC69AC25C */
sa1  =  1.96512716674392571292e+01, /* 0x4033A6B9, 0xBD707687 */
sa2  =  1.37657754143519042600e+02, /* 0x4061350C, 0x526AE721 */
sa3  =  4.34565877475229228821e+02, /* 0x407B290D, 0xD58A1A71 */
sa4  =  6.45387271733267880336e+02, /* 0x40842B19, 0x21EC2868 */
sa5  =  4.29008140027567833386e+02, /* 0x407AD021, 0x57700314 */
sa6  =  1.08635005541779435134e+02, /* 0x405B28A3, 0xEE48AE2C */
sa7  =  6.57024977031928170135e+00, /* 0x401A47EF, 0x8E484A93 */
sa8  = -6.04244152148580987438e-02, /* 0xBFAEEFF2, 0xEE749A62 */
/*
 * Coefficients for approximation to  erfc in [1/.35,28]
 */
rb0  = -9.86494292470009928597e-03, /* 0xBF843412, 0x39E86F4A */
rb1  = -7.99283237680523006574e-01, /* 0xBFE993BA, 0x70C285DE */
rb2  = -1.77579549177547519889e+01, /* 0xC031C209, 0x555F995A */
rb3  = -1.60636384855821916062e+02, /* 0xC064145D, 0x43C5ED98 */
rb4  = -6.37566443368389627722e+02, /* 0xC083EC88, 0x1375F228 */
rb5  = -1.02509513161107724954e+03, /* 0xC0900461, 0x6A2E5992 */
rb6  = -4.83519191608651397019e+02, /* 0xC07E384E, 0x9BDC383F */
sb1  =  3.03380607434824582924e+01, /* 0x403E568B, 0x261D5190 */
sb2  =  3.25792512996573918826e+02, /* 0x40745CAE, 0x221B9F0A */
sb3  =  1.53672958608443695994e+03, /* 0x409802EB, 0x189D5118 */
sb4  =  3.19985821950859553908e+03, /* 0x40A8FFB7, 0x688C246A */
sb5  =  2.55305040643316442583e+03, /* 0x40A3F219, 0xCEDF3BE6 */
sb6  =  4.74528541206955367215e+02, /* 0x407DA874, 0xE79FE763 */
sb7  = -2.24409524465858183362e+01; /* 0xC03670E2, 0x42712D62 */

double
erf(double x) 
{
        int32_t hx,ix,i;
        double R,S,P,Q,s,y,z,r;
        GET_HIGH_WORD(hx,x);
        ix = hx&0x7fffffff;
        if(ix>=0x7ff00000) {                /* erf(nan)=nan */
            i = ((u_int32_t)hx>>31)<<1;
            return (double)(1-i)+one/x;        /* erf(+-inf)=+-1 */
        }

        if(ix < 0x3feb0000) {                /* |x|<0.84375 */
            if(ix < 0x3e300000) {         /* |x|<2**-28 */
                if (ix < 0x00800000) 
                    return 0.125*(8.0*x+efx8*x);  /*avoid underflow */
                return x + efx*x;
            }
            z = x*x;
            r = pp0+z*(pp1+z*(pp2+z*(pp3+z*pp4)));
            s = one+z*(qq1+z*(qq2+z*(qq3+z*(qq4+z*qq5))));
            y = r/s;
            return x + x*y;
        }
        if(ix < 0x3ff40000) {                /* 0.84375 <= |x| < 1.25 */
            s = fabs(x)-one;
            P = pa0+s*(pa1+s*(pa2+s*(pa3+s*(pa4+s*(pa5+s*pa6)))));
            Q = one+s*(qa1+s*(qa2+s*(qa3+s*(qa4+s*(qa5+s*qa6)))));
            if(hx>=0) return erx + P/Q; else return -erx - P/Q;
        }
        if (ix >= 0x40180000) {                /* inf>|x|>=6 */
            if(hx>=0) return one-tiny; else return tiny-one;
        }
        x = fabs(x);
        s = one/(x*x);
        if(ix< 0x4006DB6E) {        /* |x| < 1/0.35 */
            R=ra0+s*(ra1+s*(ra2+s*(ra3+s*(ra4+s*(
                                ra5+s*(ra6+s*ra7))))));
            S=one+s*(sa1+s*(sa2+s*(sa3+s*(sa4+s*(
                                sa5+s*(sa6+s*(sa7+s*sa8)))))));
        } else {        /* |x| >= 1/0.35 */
            R=rb0+s*(rb1+s*(rb2+s*(rb3+s*(rb4+s*(
                                rb5+s*rb6)))));
            S=one+s*(sb1+s*(sb2+s*(sb3+s*(sb4+s*(
                                sb5+s*(sb6+s*sb7))))));
        }
        z  = x;  
        SET_LOW_WORD(z,0);
        r  =  exp(-z*z-0.5625)*exp((z-x)*(z+x)+R/S);
        if(hx>=0) return one-r/x; else return  r/x-one;
}

double
erfc(double x) 
{
        int32_t hx,ix;
        double R,S,P,Q,s,y,z,r;
        GET_HIGH_WORD(hx,x);
        ix = hx&0x7fffffff;
        if(ix>=0x7ff00000) {                        /* erfc(nan)=nan */
                                                /* erfc(+-inf)=0,2 */
            return (double)(((u_int32_t)hx>>31)<<1)+one/x;
        }

        if(ix < 0x3feb0000) {                /* |x|<0.84375 */
            if(ix < 0x3c700000)          /* |x|<2**-56 */
                return one-x;
            z = x*x;
            r = pp0+z*(pp1+z*(pp2+z*(pp3+z*pp4)));
            s = one+z*(qq1+z*(qq2+z*(qq3+z*(qq4+z*qq5))));
            y = r/s;
            if(hx < 0x3fd00000) {          /* x<1/4 */
                return one-(x+x*y);
            } else {
                r = x*y;
                r += (x-half);
                return half - r ;
            }
        }
        if(ix < 0x3ff40000) {                /* 0.84375 <= |x| < 1.25 */
            s = fabs(x)-one;
            P = pa0+s*(pa1+s*(pa2+s*(pa3+s*(pa4+s*(pa5+s*pa6)))));
            Q = one+s*(qa1+s*(qa2+s*(qa3+s*(qa4+s*(qa5+s*qa6)))));
            if(hx>=0) {
                z  = one-erx; return z - P/Q; 
            } else {
                z = erx+P/Q; return one+z;
            }
        }
        if (ix < 0x403c0000) {                /* |x|<28 */
            x = fabs(x);
            s = one/(x*x);
            if(ix< 0x4006DB6D) {        /* |x| < 1/.35 ~ 2.857143*/
                R=ra0+s*(ra1+s*(ra2+s*(ra3+s*(ra4+s*(
                                ra5+s*(ra6+s*ra7))))));
                S=one+s*(sa1+s*(sa2+s*(sa3+s*(sa4+s*(
                                sa5+s*(sa6+s*(sa7+s*sa8)))))));
            } else {                        /* |x| >= 1/.35 ~ 2.857143 */
                if(hx<0&&ix>=0x40180000) return two-tiny;/* x < -6 */
                R=rb0+s*(rb1+s*(rb2+s*(rb3+s*(rb4+s*(
                                rb5+s*rb6)))));
                S=one+s*(sb1+s*(sb2+s*(sb3+s*(sb4+s*(
                                sb5+s*(sb6+s*sb7))))));
            }
            z  = x;
            SET_LOW_WORD(z,0);
            r  =  exp(-z*z-0.5625)*
                        exp((z-x)*(z+x)+R/S);
            if(hx>0) return r/x; else return two-r/x;
        } else {
            if(hx>0) return tiny*tiny; else return two-tiny;
        }
}

#endif


void StgFEM_StandardConditionFunctions_SolidBodyRotation( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	FeVariable*             feVariable         = NULL;
	FeMesh*			mesh               = NULL;
	double*                 result             = (double*) _result;
	double*                 coord;
	Coord                   centre;
	Coord                   vector;
	double                  omega;
	
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;

	/* Find Centre of Solid Body Rotation */
	centre[ I_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreX", 0.0  );
	centre[ J_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreY", 0.0  );
	centre[ K_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreZ", 0.0  );
	omega            = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationOmega", 1.0  );

	/* Find coordinate of node */
	coord = Mesh_GetVertex( mesh, node_lI );

	/* Find vector from centre to node */
	StGermain_VectorSubtraction( vector, coord, centre, 2 );

	result[ I_AXIS ] = -omega * vector[ J_AXIS ];
	result[ J_AXIS ] =  omega * vector[ I_AXIS ];
}


void StgFEM_StandardConditionFunctions_PartialRotationX( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	FeVariable*             feVariable         = NULL;
	FeMesh*			mesh               = NULL;
	double*                 result             = (double*) _result;
	double*                 coord;
	Coord                   centre;
	Coord                   vector;
	double                  omega;
	double			size;

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;

	/* Find Centre of Solid Body Rotation */
	centre[ I_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreX", 0.0  );
	centre[ J_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreY", 0.0  );
	centre[ K_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreZ", 0.0  );
	size             = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"RadiusCylinder", 0.0  );
	omega            = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationOmega", 1.0  );

	/* Find coordinate of node */
	coord = Mesh_GetVertex( mesh, node_lI );

	/* Find vector from centre to node */
	StGermain_VectorSubtraction( vector, coord, centre, 2 );

	/*if (context->currentTime > 1.33e-6)
	  omega=0.0;*/
	
	if ((vector[ I_AXIS ]*vector[ I_AXIS ]+vector[ J_AXIS ]*vector[ J_AXIS ])<=size*size)
		*result = -omega * vector[ J_AXIS ];
	else
		*result = 0.0;
}

void StgFEM_StandardConditionFunctions_PartialRotationY( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	FeVariable*             feVariable         = NULL;
	FeMesh*			mesh               = NULL;
	double*                 result             = (double*) _result;
	double*                 coord;
	Coord                   centre;
	Coord                   vector;
	double                  omega;
	double			size;

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;

	/* Find Centre of Solid Body Rotation */
	centre[ I_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreX", 0.0  );
	centre[ J_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreY", 0.0  );
	centre[ K_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreZ", 0.0  );
	size             = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"RadiusCylinder", 0.0  );
	omega            = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationOmega", 1.0  );

	/* Find coordinate of node */
	coord = Mesh_GetVertex( mesh, node_lI );

	/* Find vector from centre to node */
	StGermain_VectorSubtraction( vector, coord, centre, 2 );
	
	if ((vector[ I_AXIS ]*vector[ I_AXIS ]+vector[ J_AXIS ]*vector[ J_AXIS ])<=size*size)
		*result =  omega * vector[ I_AXIS ];
	else 
		*result = 0.0;
}


void StgFEM_StandardConditionFunctions_TaperedRotationX( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	FeVariable*             feVariable         = NULL;
	FeMesh*			mesh               = NULL;
	double*                 result             = (double*) _result;
	double*                 coord;
	Coord                   centre;
	Coord                   vector;
	double                  omega;
	double			size, r, taper;

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;

	/* Find Centre of Solid Body Rotation */
	centre[ I_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationCentreX", 0.0 );
	centre[ J_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationCentreY", 0.0 );
	centre[ K_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationCentreZ", 0.0 );
	size             = Dictionary_GetDouble_WithDefault( dictionary, "RadiusCylinder", 0.0 );
	omega            = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationOmega",   1.0 );

	taper            = Dictionary_GetDouble_WithDefault( dictionary, "TaperedRadius",   0.0 );

	/* Find coordinate of node */
	coord = Mesh_GetVertex( mesh, node_lI );

	/* Find vector from centre to node */
	StGermain_VectorSubtraction( vector, coord, centre, 2 );

        r=sqrt(vector[ I_AXIS ]*vector[ I_AXIS ]
               +vector[ J_AXIS ]*vector[ J_AXIS ]);
	if (r<=size)
          *result = -omega * vector[ J_AXIS ];
	else if(r<=taper)
          *result = -omega * vector[ J_AXIS ]*(taper-r)/(taper-size);
        else
          *result = 0;
}

void StgFEM_StandardConditionFunctions_TaperedRotationY( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	FeVariable*             feVariable         = NULL;
	FeMesh*			mesh               = NULL;
	double*                 result             = (double*) _result;
	double*                 coord;
	Coord                   centre;
	Coord                   vector;
	double                  omega;
	double			size, r, taper;

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;

	/* Find Centre of Solid Body Rotation */
	centre[ I_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationCentreX", 0.0 );
	centre[ J_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationCentreY", 0.0 );
	centre[ K_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationCentreZ", 0.0 );
	size             = Dictionary_GetDouble_WithDefault( dictionary, "RadiusCylinder", 0.0 );
	omega            = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationOmega",   1.0 );

	taper            = Dictionary_GetDouble_WithDefault( dictionary, "TaperedRadius",   0.0 );

	/* Find coordinate of node */
	coord = Mesh_GetVertex( mesh, node_lI );

	/* Find vector from centre to node */
	StGermain_VectorSubtraction( vector, coord, centre, 2 );


        r=sqrt(vector[ I_AXIS ]*vector[ I_AXIS ]
               +vector[ J_AXIS ]*vector[ J_AXIS ]);
	if (r<=size)
          *result = omega * vector[ I_AXIS ];
	else if(r<=taper)
          *result = omega * vector[ I_AXIS ]*(taper-r)/(taper-size);
        else
          *result = 0;
}




void StgFEM_StandardConditionFunctions_SimpleShear( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	FeVariable*             feVariable         = NULL;
	FeMesh*			mesh               = NULL;
	double*                 result             = (double*) _result;
	double*                 coord;
	double                  centre;
	double                  factor;
	
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;

	/* Find Centre of Solid Body Rotation */
	centre = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SimpleShearCentreY", 0.0  );
	factor = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SimpleShearFactor", 1.0  );

	/* Find coordinate of node */
	coord = Mesh_GetVertex( mesh, node_lI );

	*result = factor * (coord[ J_AXIS ] - centre);
}

void StgFEM_StandardConditionFunctions_ShearZ( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	FeVariable*             feVariable         = NULL;
	FeMesh*			mesh               = NULL;
	double*                 result             = (double*) _result;
	double*                 coord;
	double                  centre;
	double                  factor;
	
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;

	/* Find Centre of Solid Body Rotation */
	centre = Dictionary_GetDouble_WithDefault( dictionary, "ShearZCentre", 0.0 );
	factor = Dictionary_GetDouble_WithDefault( dictionary, "ShearZFactor", 1.0 );

	/* Find coordinate of node */
	coord = Mesh_GetVertex( mesh, node_lI );

	*result = factor * (coord[ K_AXIS ] - centre);
}

void StgFEM_StandardConditionFunctions_SimpleShearInverted( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
        DomainContext*  context            = (DomainContext*)_context;
        Dictionary*             dictionary         = context->dictionary;
        FeVariable*             feVariable         = NULL;
        FeMesh*                 mesh               = NULL;
        double*                 result             = (double*) _result;
        double*                 coord;
        double                  centre;
        double                  factor;
        double                  yAxisInvert;

        feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
        mesh       = feVariable->feMesh;

        /* Find Centre of Solid Body Rotation */
        centre = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SimpleShearCentreY", 0.0  );
        factor = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SimpleShearFactor", 1.0  );

        /* Find coordinate of node */
        coord = Mesh_GetVertex( mesh, node_lI );

        yAxisInvert = coord[ J_AXIS ] * -1.0 - 1.0;

        *result = factor * ( 1.0 - coord[ J_AXIS ] ) ;
}


void StgFEM_StandardConditionFunctions_Extension( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	FeVariable*             feVariable         = NULL;
	FeMesh*			mesh               = NULL;
	double*                 result             = (double*) _result;
	double*                 coord;
	double                  centre;
	double                  factor;
	
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;

	/* Find Centre of Solid Body Rotation */
	centre = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"ExtensionCentreX", 0.0  );
	factor = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"ExtensionFactor", 1.0  );

	/* Find coordinate of node */
	coord = Mesh_GetVertex( mesh, node_lI );

	*result = factor * (coord[ I_AXIS ] - centre);
}


void StgFEM_StandardConditionFunctions_PartialLid_TopLayer( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* result ) {
	DomainContext*	context = (DomainContext*)_context;
	FeVariable*             velVar = NULL;
	FeMesh*			mesh = NULL;
	double*			velResult = (double*)result;
	double                  margin = 0;
	double			min[3], max[3];
	
	velVar = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh = velVar->feMesh;

	Mesh_GetMinimumSeparation( mesh, &margin, NULL );
	Mesh_GetGlobalCoordRange( mesh, min, max );
	margin *= 1.1;
	if( (Mesh_GetVertex( mesh, node_lI )[I_AXIS] < (max[I_AXIS] - margin )) && 
	    (Mesh_GetVertex( mesh, node_lI )[I_AXIS] > (min[I_AXIS] + margin )))
	{
		(*velResult) = 1;
	}
	else {
		(*velResult) = 0;
	}
}

void StgFEM_StandardConditionFunctions_LinearInterpolationLid( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* result ) {
	DomainContext*	context = (DomainContext*)_context;
	FeVariable*             velVar = NULL;
	FeMesh*			mesh = NULL;
	double*			velResult = (double*)result;
	double			boxLength = 0;
	double			leftHandSideValue = 0;
	double			rightHandSideValue = 0;
	double			gradient = 0;
	double			min[3], max[3];
	
	velVar = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh = velVar->feMesh;

	Mesh_GetGlobalCoordRange( mesh, min, max );
	boxLength = max[I_AXIS] - min[I_AXIS];
	leftHandSideValue = Dictionary_GetDouble_WithDefault( context->dictionary, (Dictionary_Entry_Key)"bcLeftHandSideValue", 0.0  );
	rightHandSideValue = Dictionary_GetDouble_WithDefault( context->dictionary, (Dictionary_Entry_Key)"bcRightHandSideValue", 1.0 );
	gradient = (rightHandSideValue - leftHandSideValue) / boxLength;
	(*velResult ) = leftHandSideValue + gradient * (Mesh_GetVertex( mesh, node_lI )[I_AXIS] - min[I_AXIS] );
}


void StgFEM_StandardConditionFunctions_Lid_RampWithCentralMax( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* result ) {
	DomainContext*	context = (DomainContext*)_context;
	FeVariable*             velVar = NULL;
	FeMesh*			mesh = NULL;
	double*			velResult = (double*)result;
	double			boxLength = 0;
	double			xPosRelativeToTopLeft = 0;
	double			min[3], max[3];
	
	velVar = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh = velVar->feMesh;

	Mesh_GetGlobalCoordRange( mesh, min, max );
	xPosRelativeToTopLeft = Mesh_GetVertex( mesh, node_lI )[I_AXIS] - min[I_AXIS];
	boxLength = max[I_AXIS] - min[I_AXIS];
	if ( xPosRelativeToTopLeft < boxLength / 2 ) {
		(*velResult) =  2 * xPosRelativeToTopLeft / boxLength;
	}
	else {
		(*velResult) = 1 - 2 * ( xPosRelativeToTopLeft - (boxLength/2) );
	}
}
void StgFEM_StandardConditionFunctions_LinearVelocityLeftWall( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* result ) {
	DomainContext*	context = (DomainContext*)_context;
	FeVariable*             velVar = NULL;
	FeMesh*			mesh = NULL;
	double*			velResult = (double*)result;
	Dictionary*             dictionary         = context->dictionary;
	double			min[3], max[3];
	double			gradient, maxvel;
	velVar = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh = velVar->feMesh;

	maxvel = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"MaximumVelocity_Left", 0.0  );
	Mesh_GetGlobalCoordRange( mesh, min, max );
	gradient = maxvel/(min[1] - max[1]);
	
	(*velResult)   = gradient*Mesh_GetVertex( mesh, node_lI )[J_AXIS];
	 //printf("Left velResult is %g\n",(*velResult));
	 
}
void StgFEM_StandardConditionFunctions_LinearVelocityRightWall( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* result ) {
	DomainContext*	context = (DomainContext*)_context;
	FeVariable*             velVar = NULL;
	FeMesh*			mesh = NULL;
	double*			velResult = (double*)result;
	Dictionary*             dictionary         = context->dictionary;
	double			min[3], max[3];
	double			gradient, maxvel;
	velVar = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh = velVar->feMesh;

	maxvel = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"MaximumVelocity_Right", 0.0  );
	Mesh_GetGlobalCoordRange( mesh, min, max );
	gradient = maxvel/(max[1] - min[1]);
	 
	(*velResult)   = maxvel - gradient*Mesh_GetVertex( mesh, node_lI )[J_AXIS];
	//printf("Right velResult is %g\n",(*velResult));
}


void StgFEM_StandardConditionFunctions_SinusoidalLid( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* result ) {
	DomainContext*	context = (DomainContext*)_context;
	FeVariable*             velVar = NULL;
	FeMesh*			mesh = NULL;
	double*			velResult = (double*)result;
	double			boxLength = 0;
	double			linearInterp = 0;
	double          	wavenumber;
	double			min[3], max[3];

	wavenumber = Dictionary_GetDouble_WithDefault( context->dictionary, (Dictionary_Entry_Key)"sinusoidalLidWavenumber", 1 );
	
	velVar = (FeVariable* )FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh = velVar->feMesh;

	Mesh_GetGlobalCoordRange( mesh, min, max );
	boxLength = max[I_AXIS] - min[I_AXIS];
	linearInterp = (Mesh_GetVertex( mesh, node_lI )[I_AXIS] - min[I_AXIS] ) / boxLength;
	(*velResult) = sin( linearInterp * M_PI * wavenumber );
}


void StgFEM_StandardConditionFunctions_CornerOnly( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* result ) {
	DomainContext*	context = (DomainContext*)_context;
	FeVariable*             velVar = NULL;
	FeMesh*			feMesh = NULL;
	double*			velResult = (double*)result;
	Node_GlobalIndex	node_gI = 0;
	unsigned		inds[3];
	Grid*			elGrid;

	velVar = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	feMesh = velVar->feMesh;
	elGrid = *(Grid**)ExtensionManager_Get( feMesh->info, feMesh, 
						ExtensionManager_GetHandle( feMesh->info, (Name)"elGrid" )  );

	node_gI = Mesh_DomainToGlobal( feMesh, MT_VERTEX, node_lI );
	RegularMeshUtils_Node_1DTo3D( feMesh, node_gI, inds );
	
	if ( inds[0] == elGrid->sizes[I_AXIS] ) {
		(*velResult) = 1;
	}
	else {
		(*velResult) = 0;
	}
}

double StGermain_CosineHillValue( double* centre, double* position, double height, double diameterAtBase, Dimension_Index dim ) {
	double distanceFromCentre = StGermain_DistanceBetweenPoints( centre, position, dim );
	
	if (distanceFromCentre < diameterAtBase * 0.5 ) 
		return height * (0.5 + 0.5 * cos( 2.0 * M_PI/diameterAtBase * distanceFromCentre ) );
	else
		return 0.0;
}

void StgFEM_StandardConditionFunctions_TemperatureCosineHill( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	FeVariable*             feVariable         = NULL;
	FeMesh*			feMesh               = NULL;
	double*                 result             = (double*) _result;
	Coord                   centre;
	Coord                   rotationCentre;
	double                  omega;
	double                  hillHeight;
	double                  hillDiameter;
	
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	feMesh       = feVariable->feMesh;

	/* Read values from dictionary */
	hillHeight       = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"CosineHillHeight"  , 1.0  );
	hillDiameter     = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"CosineHillDiameter", 1.0  );
	centre[ I_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"CosineHillCentreX" , 0.0  );
	centre[ J_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"CosineHillCentreY" , 0.0  );
	centre[ K_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"CosineHillCentreZ" , 0.0  );

	if ( Dictionary_GetBool( dictionary, "RotateCosineHill" ) ) {
		/* Assume solid body rotation */
		rotationCentre[ I_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreX", 0.0  );
		rotationCentre[ J_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreY", 0.0  );
		rotationCentre[ K_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreZ", 0.0  );
		omega                    = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationOmega", 1.0  );

		StGermain_VectorSubtraction( centre, rotationCentre, centre, context->dim );
		StGermain_RotateCoordinateAxis( centre, centre, K_AXIS, omega * context->currentTime );
		StGermain_VectorAddition( centre, centre, rotationCentre, context->dim );
	}

	*result = StGermain_CosineHillValue( centre, Mesh_GetVertex( feMesh, node_lI ), hillHeight, hillDiameter, context->dim );
}


void StgFEM_StandardConditionFunctions_LinearWithSinusoidalPerturbation( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	FeVariable*             feVariable = NULL;
	FeMesh*			feMesh = NULL;
	unsigned		nDims;
	double*                 result = (double*) _result;
	double                  topLayerBC;
	double                  bottomLayerBC;
	double                  perturbationAmplitude;
	double                  horizontalWaveNumber;
	double                  verticalWaveNumber;
	double                  scaleFactor;
	double*                 coord;
	Coord                   relScaledCoord; 
	double			min[3], max[3], topLayerCoord, bottomLayerCoord;
	
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	feMesh       = feVariable->feMesh;

	nDims = Mesh_GetDimSize( feMesh );
	Mesh_GetGlobalCoordRange( feMesh, min, max );

	topLayerCoord = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_TopLayerCoord", max[J_AXIS]  );
	bottomLayerCoord = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_BottomLayerCoord", min[J_AXIS]  );

	topLayerBC = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_TopLayerBC", 0.0  );
	bottomLayerBC = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_BottomLayerBC", 1.0  );
	scaleFactor = bottomLayerBC - topLayerBC;
	perturbationAmplitude = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_PerturbationAmplitude", 0.1  );
	/* Note: these are both multiplied by pi, so wavenumber = 1 means the perturbation goes from 0 to pi, which is
	 * half a full sin or cos cycle. Wavenumber = 3 means the range is 0 -> 3pi, or 1 and a half full cycles. */
	horizontalWaveNumber = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_HorizontalWaveNumber", 1.0  );
	verticalWaveNumber = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_VerticalWaveNumber", 1.0  );

	coord = Mesh_GetVertex( feMesh, node_lI );

	/* if node is outside IC shape set to 0 temperature */
	if( coord[J_AXIS] > topLayerCoord || coord[J_AXIS] < bottomLayerCoord ) {
		*result = 0; return ;
	}

	/* make coord relative to box bottom left corner, then scale from 0 to 1 between box min & max */
	relScaledCoord[I_AXIS] = (coord[0] - min[0]) / (max[0] - min[0]);
	relScaledCoord[J_AXIS] = (coord[1] - bottomLayerCoord) / (topLayerCoord - bottomLayerCoord);


	/* Note: ok to use the 1.0 below since we've already scaled the coord to somewhere between 0 to 1 */
	*result = topLayerBC + scaleFactor * ( 1.0 - relScaledCoord[ J_AXIS ] )
		+ perturbationAmplitude * ( cos( horizontalWaveNumber * M_PI * coord[ I_AXIS ] )
					    * sin( verticalWaveNumber * M_PI * relScaledCoord[ J_AXIS ] ) );
}

void StgFEM_StandardConditionFunctions_Trigonometry( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*			feMesh               = NULL;
	double*                 result             = (double*) _result;
	double*                 coord;
	double                  height, width;
	double			min[3], max[3];

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	feMesh       = feVariable->feMesh;

	Mesh_GetGlobalCoordRange( feMesh, min, max );
	coord = Mesh_GetVertex( feMesh, node_lI );

	/* Get Aspect Ratio */
	height = max[ J_AXIS ] - min[ J_AXIS ];
	width  = max[ I_AXIS ] - min[ I_AXIS ];
	
	*result = 1.0 - 0.5 * M_PI * coord[ J_AXIS ] * sin( M_PI * coord[ I_AXIS ]/width );
}

#define SMALL 1.0e-5
void Stg_FEM_VelicTemperatureIC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*  context            = (DomainContext*)_context;
	FeVariable*             temperatureField   = (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	FeMesh*			feMesh               = temperatureField->feMesh;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double*                 coord;
	double                  x; 
	double                  y;
	double                  kx;
	double                  ky;
	int                     wavenumberX;
	double                  wavenumberY;
	double                  sigma;
	double                  Lx;
	double			min[3], max[3];
	
	/* Find coordinate of node */
	coord = Mesh_GetVertex( feMesh, node_lI );
	Mesh_GetGlobalCoordRange( feMesh, min, max );

	/* Make sure that the box has right dimensions */
	assert( ( max[ J_AXIS ] - min[ J_AXIS ] - 1.0 ) < SMALL );
	Lx = max[ I_AXIS ] - min[ I_AXIS ];

	x = coord[ I_AXIS ] - min[ I_AXIS ];
	y = coord[ J_AXIS ] - min[ J_AXIS ];

	wavenumberX = Dictionary_GetInt_WithDefault( dictionary, (Dictionary_Entry_Key)"wavenumberX", 1  );
	wavenumberY = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"wavenumberY", 1.0  );
	sigma = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"sigma", 1.0 );
	
	assert( sigma > 0.0 );
	assert( wavenumberY > 0.0 );
	assert( wavenumberX > 0.0 );
	
	kx = (double)wavenumberX * M_PI / Lx;
	ky = (double)wavenumberY * M_PI;

	*result = sigma * sin( ky * y ) * cos( kx * x  );
}

/* IC from Mirko Velic. This is the IC temperature for his solB, from his Analytic Suite. Added 22-May-2006 */
void Stg_FEM_VelicTemperatureIC_SolB( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*  context            = (DomainContext*)_context;
	FeVariable*             temperatureField   = (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	FeMesh*			feMesh               = temperatureField->feMesh;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double*                 coord;
	double                  x; 
	double                  y;
	double                  km; /*  for y-direction */
	double                  kn; /*  for x-direction */
	double                  wavenumberX;
	double                  wavenumberY;
	double                  L;
	double                  sigma;
	double			min[3], max[3];
	
	/* Find coordinate of node */
	coord = Mesh_GetVertex( feMesh, node_lI );
	Mesh_GetGlobalCoordRange( feMesh, min, max );

	/* Make sure that the box has right dimensions */
	assert( (max[ J_AXIS ] - min[ J_AXIS ] - 1.0 ) < SMALL );
	L = max[ I_AXIS ] - min[ I_AXIS ];

	x = coord[ I_AXIS ] - min[ I_AXIS ];
	y = coord[ J_AXIS ] - min[ J_AXIS ];

	wavenumberX = Dictionary_GetInt_WithDefault( dictionary, (Dictionary_Entry_Key)"wavenumberX", 1  );
	wavenumberY = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"wavenumberY", 2.0 );
	assert( wavenumberX != wavenumberY  );
	sigma = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"sigma", 1.0 );

	kn = wavenumberX * M_PI / L;
/* 	 TODO: Re-write Mirko's code and/or Documentation so the input parameters for these ICs are less confusing */
	km = wavenumberY / L;

	*result = sigma * sinh( km * y ) * cos( kn * x  );
}


/* Initial Condition derived from Boundary Layer theory -
   taken from P. E. van Keken, S. D. King, U. R. Schmeling, U. R. Christensen, D. Neumeister, and M.-P. Doin. A comparison of methods for the modeling of thermochemical convection. Journal of Geophysical Research, 102(B10):22477-22496, october 1997. */
void StgFEM_StandardConditionFunctions_AnalyticalTemperatureIC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*			feMesh               = NULL;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double*                 coord;
	double                  u0, v0, Q;
	double                  x, y;
	double                  RaT;
	double                  lambda, height, width;
	double                  Tu, Tl, Tr, Ts;
	double			min[3], max[3];

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	feMesh       = feVariable->feMesh;

	coord      = Mesh_GetVertex( feMesh, node_lI );
	Mesh_GetGlobalCoordRange( feMesh, min, max );

	/* Get Aspect Ratio */
	height = max[ J_AXIS ] - min[ J_AXIS ];
	width  = max[ I_AXIS ] - min[ I_AXIS ];
	lambda = width/height;
	
	x = coord[ I_AXIS ] - min[ I_AXIS ];
	y = coord[ J_AXIS ] - min[ J_AXIS ];
	
	/* Get thermal Rayleigh Number from Dictionary */
	RaT = Dictionary_GetDouble( dictionary, "RaT" );
	
	/* Horizontal fluid velocity at upper boundary & lower boundary - Equation A3 */
	u0 = pow( lambda , 7.0/3.0 )/ pow(1 + lambda*lambda*lambda*lambda, 2.0/3.0) * pow(0.5*RaT/sqrt(M_PI) , 2.0/3.0);

	/* Vertical velocity of the upwelling and downwellings - Modified from Van Keken to match Turcotte and Shubert */
	v0 = u0; /*lambda; */
	
	/* Total rate of heat flow out of the top of the cell per unit distance along the axis of the roll - Equation A3 */
	Q = 2.0 * sqrt(M_1_PI * lambda/u0);
	Tu = 0.5 * erf( 0.5 * ( 1 - y ) * sqrt(u0/x) );                                                      /* Equation A2a */
	Tl = 1.0 - 0.5 * erf(0.5 * y * sqrt(u0/(lambda-x)));                                                 /* Equation A2b */
	Tr = 0.5 + 0.5*Q/sqrt(M_PI) * sqrt(v0/(y+1)) * exp( -x*x*v0/(4*y+4) );                               /* Equation A2c */
	Ts = 0.5 - 0.5*Q/sqrt(M_PI) * sqrt(v0/(2-y)) * exp( -(lambda - x) * (lambda - x) * v0 / (8 - 4*y) ); /* Equation A2d */

	/* Equation A1 */
	*result = Tu + Tl + Tr + Ts - 1.5;

	/* Crop result */
	if ( *result > 1.0 ) 
		*result = 1.0;
	else if ( *result < 0.0 ) 
		*result = 0.0;
	
}

void StgFEM_StandardConditionFunctions_EdgeDriveConvectionIC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) 
{        
	DomainContext*  context = (DomainContext*)_context;        
	Dictionary*             dictionary         = context->dictionary;        
	FeVariable*             feVariable = NULL;        
	FeMesh*			mesh = NULL;        
	double*                 result = (double*) _result;        
	double                  perturbationAmplitude;        
	double                  thermalAnomalyOffset;        
	double*                 coord;        
	
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );        
	mesh       = feVariable->feMesh;        
	perturbationAmplitude = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_PerturbationAmplitude", 0.1  );        
	thermalAnomalyOffset = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"thermalAnomalyOffset", 0.0  );        
	coord = Mesh_GetVertex( mesh, node_lI );
	
	/* eqn 1 from S.D.King & D.L. Anderson, "Edge-drive convection", EPSL 160 (1998) 289-296 */        
	
	*result = 1.0 + perturbationAmplitude * sin( M_PI * coord[ J_AXIS ] ) * cos( 0.5 * M_PI * ( coord[ I_AXIS ] + thermalAnomalyOffset ) );
}

void StgFEM_StandardConditionFunctions_ThermalEdgeDriveConvectionIC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result )
{
        DomainContext*  context = (DomainContext*)_context;
        Dictionary*             dictionary         = context->dictionary;
        FeVariable*             feVariable = NULL;
        FeMesh*                 mesh = NULL;
        double*                 result = (double*) _result;
        double*                 coord;
        int                     dim;
        double                  contStartX, contEndX;
        double                  contStartY, contEndY;
        double                  contStartZ, contEndZ;
        double                  minY, maxY, interiorTemp;

        feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
        mesh       = feVariable->feMesh;
        coord = Mesh_GetVertex( mesh, node_lI );
        
	dim = Dictionary_GetInt_WithDefault( dictionary, (Dictionary_Entry_Key)"dim", 0.0  );
        contStartX = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"contStartX", 0.0  );
        contEndX = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"contEndX", 0.0  );
        contStartY = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"contStartY", 0.0  );
        contEndY = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"contEndY", 0.0  );
        minY = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"minY", 0.0  );
        maxY = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"maxY", 0.0  );
	interiorTemp = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"interiorTemp", 1.0 );
        if ( dim == 3  ) {
                contStartZ = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"contStartZ", 0.0  );
                contEndZ = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"contEndZ", 0.0 );
        }

        if(( coord[I_AXIS] >= contStartX && coord[ I_AXIS ] <= contEndX ) && ( coord[J_AXIS] >= contStartY && coord[ J_AXIS ] <= contEndY )) {
                if ( dim == 3 ) {
                        if ( coord[K_AXIS] >= contStartZ && coord[ K_AXIS ] <= contEndZ  )
                                        *result = 0.0;
                        else
                                        *result = interiorTemp;
                }
        }
        else
                        *result = interiorTemp;
}

void StgFEM_StandardConditionFunctions_SinusoidalExtension( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double                  frequency;
	double                  vel0;
	double                  amplitude;
	double                  phaseShift;

	frequency  = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalExtensionFrequency", 1.0  );
	vel0       = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalExtensionVelocity", 0.0  );
	amplitude  = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalExtensionAmplitude", 0.0  );
	phaseShift = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalExtensionPhaseShift", 0.0 );


	*result = vel0 + amplitude * cos( 2.0 * M_PI * frequency * (context->currentTime + context->dt - phaseShift )  );
}


void StgFEM_StandardConditionFunctions_StepFunction( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*			feMesh               = NULL;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
        double*                 coord;
	double                  lower_offset, upper_offset;
	double                  value, lower_value, upper_value;
	unsigned		dim;

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	feMesh       = feVariable->feMesh;
	coord      = Mesh_GetVertex( feMesh, node_lI );

	lower_offset = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionLowerOffset", 0.0 );
	upper_offset = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionUpperOffset", lower_offset );
	value = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionValue", 0.0 );
	dim = Dictionary_GetUnsignedInt_WithDefault( dictionary, "StepFunctionDim", 0 );

        lower_value = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionLowerValue", 0.0 );
        upper_value = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionUpperValue", value );

        if(dim==3)
          {
            dim=0;
            coord=&(context->currentTime);
          }

        if(coord[dim] < lower_offset) {
          *result=lower_value;
        } else if(coord[dim] < upper_offset) {
          *result=lower_value + 
            (upper_value-lower_value)
            *(coord[dim] - lower_offset)/(upper_offset-lower_offset);
        } else {
          *result=upper_value;
        }
}


void StG_FEM_StandardConditionFunctions_StepFunctionProduct1( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*     mesh               = NULL;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
        double*                 coord;
	double                  start, end;
	double                  value;
	unsigned		dim;

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;
	coord      = Mesh_GetVertex( mesh, node_lI );

	start = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct1Start", 0.0 );
	end = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct1End", 0.0 );
	value = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct1Value", 0.0 );
	dim = Dictionary_GetUnsignedInt_WithDefault( dictionary, "StepFunctionProduct1Dim", 0 );

        if( coord[dim] > start && coord[dim] < end ) {
          *result = value;
        }
        else {
          *result = 0;
        }
}

void StG_FEM_StandardConditionFunctions_StepFunctionProduct2( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*     mesh               = NULL;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
        double*                 coord;
	double                  start, end;
	double                  value;
	unsigned		dim;

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;
	coord      = Mesh_GetVertex( mesh, node_lI );

	start = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct2Start", 0.0 );
	end = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct2End", 0.0 );
	value = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct2Value", 0.0 );
	dim = Dictionary_GetUnsignedInt_WithDefault( dictionary, "StepFunctionProduct2Dim", 0 );

        if( coord[dim] > start && coord[dim] < end ) {
          *result = value;
        }
        else {
          *result = 0;
        }
}


void StG_FEM_StandardConditionFunctions_StepFunctionProduct3( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*     mesh               = NULL;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
        double*                 coord;
	double                  start, end;
	double                  value;
	unsigned		dim;

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;
	coord      = Mesh_GetVertex( mesh, node_lI );

	start = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct3Start", 0.0 );
	end = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct3End", 0.0 );
	value = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct3Value", 0.0 );
	dim = Dictionary_GetUnsignedInt_WithDefault( dictionary, "StepFunctionProduct3Dim", 1 );

        if( coord[dim] > start && coord[dim] < end ) {
          *result = value;
        }
        else {
          *result = 0;
        }
}

void StG_FEM_StandardConditionFunctions_StepFunctionProduct4( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*     mesh               = NULL;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
        double*                 coord;
	double                  start, end;
	double                  value;
	unsigned		dim;

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;
	coord      = Mesh_GetVertex( mesh, node_lI );

	start = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct4Start", 0.0 );
	end = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct4End", 0.0 );
	value = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct4Value", 0.0 );
	dim = Dictionary_GetUnsignedInt_WithDefault( dictionary, "StepFunctionProduct4Dim", 1 );

        if( coord[dim] > start && coord[dim] < end ) {
          *result = value;
        }
        else {
          *result = 0;
        }
}

/* A Gaussian GaussianHeight*exp(-((GaussianCenter-x)/GaussianWidth)^2) */

void StG_FEM_StandardConditionFunctions_Gaussian
( Node_LocalIndex node_lI, Variable_Index var_I, void* _context,
  void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*     mesh               = NULL;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
        double*                 coord;
	double                  center, width, height;
	unsigned		dim;

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;
	coord      = Mesh_GetVertex( mesh, node_lI );

        center = Dictionary_GetDouble_WithDefault( dictionary,
                                                   "GaussianCenter", 0.0 );
	width = Dictionary_GetDouble_WithDefault( dictionary,
                                                  "GaussianWidth", 1.0 );
	height = Dictionary_GetDouble_WithDefault( dictionary,
                                                   "GaussianHeight", 1.0 );
	dim = Dictionary_GetUnsignedInt_WithDefault( dictionary,
                                                     "GaussianDim", 0 );

        *result=height*exp(-(center-coord[dim])*(center-coord[dim])
                           /(width*width));
}

void StgFEM_StandardConditionFunctions_MovingStepFunction( Node_LocalIndex nodeInd, Variable_Index varInd, void* _ctx, void* _result ) {
   FiniteElementContext* ctx = (FiniteElementContext*)_ctx;
   FeVariable* velField;
   FeMesh* mesh;
   Dictionary* dict = ctx->dictionary;
   double* result = (double*)_result;
   double* coord, offsetLower, offsetUpper, left, right;
   double *wallCrd, pos;
   int dim, wallDepth;
   unsigned ijk[3];
   char* movingWall;
   Grid* grid;

   /*
   ** Get the velocity field. */
   velField = (FeVariable*)FieldVariable_Register_GetByName(
      ctx->fieldVariable_Register, "VelocityField" );

   /*
   ** Get the mesh and the coordinate of the node. */
   mesh = velField->feMesh;
   coord = Mesh_GetVertex( mesh, nodeInd );

   /*
   ** Extract all the parameters we need from the dictionary. */
   offsetLower = Dictionary_GetDouble_WithDefault( dict, (Dictionary_Entry_Key)"MovingStepFunctionOffsetLower", 0.0  );
   offsetUpper = Dictionary_GetDouble_WithDefault( dict, (Dictionary_Entry_Key)"MovingStepFunctionOffsetUpper", 0.0  );
   dim = Dictionary_GetUnsignedInt_WithDefault( dict, "MovingStepFunctionDim", 0 );
   left = Dictionary_GetDouble_WithDefault( dict, (Dictionary_Entry_Key)"MovingStepFunctionLeftSide", 0.0  );
   right = Dictionary_GetDouble_WithDefault( dict, (Dictionary_Entry_Key)"MovingStepFunctionRightSide", 0.0  );
   movingWall = Dictionary_GetString_WithDefault( dict, "MovingStepFunctionMovingWall", "lower" );
   wallDepth = Dictionary_GetInt_WithDefault( dict, (Dictionary_Entry_Key)"MovingStepFunctionWallDepth", 0  );

   /*
   ** Because we're dealing with a moving step function, we need to calculate
   ** from where the offset should be applied. */
   grid = *(Grid**)Mesh_GetExtension( mesh, Grid**, "vertexGrid" );
   assert( grid );
   memset( ijk, 0, 3 * sizeof(unsigned) );
   if( !strcmp( movingWall, "lower" ) ) {
      ijk[dim] = wallDepth;
      wallCrd = Mesh_GetVertex( mesh, Grid_Project( grid, ijk ) );
      offsetLower += wallCrd[dim];
      offsetUpper += wallCrd[dim];
   }
   else {
      ijk[dim] = grid->sizes[dim] - wallDepth - 1;
      wallCrd = Mesh_GetVertex( mesh, Grid_Project( grid, ijk ) );
      offsetLower += wallCrd[dim];
      offsetUpper += wallCrd[dim];
   }

   /*
   ** Apply the set of parameters to this node. */
   pos = coord[dim];
   if( pos <= offsetLower )
      *result = left;
   else if( pos >= offsetUpper )
      *result = right;
   else {
      *result = left + ((pos - offsetLower) / (offsetUpper - offsetLower)) * (right - left);
   }
}

void StgFEM_StandardConditionFunctions_ConvectionBenchmark( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	/* This IC is for the 2D ConvectionBenchmark defined in
	 * http://www.mcc.monash.edu.au/twiki/view/Research/ConvectionBenchmarks
	 */
	
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*			mesh;
	double*                 result             = (double*) _result;
	double			min[3], max[3];
        double*                 coord;
	double                  x,y;
	double                  Lx, Ly;

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	mesh       = (FeMesh*)feVariable->feMesh;

	Mesh_GetGlobalCoordRange( mesh, min, max );
	
	Lx = max[ I_AXIS ] - min[ I_AXIS ];
	Ly = max[ J_AXIS ] - min[ J_AXIS ];
	
	coord      = Mesh_GetVertex( mesh, node_lI );

	x = ( coord[0] - min[ I_AXIS ] ) / Lx;
	y = ( coord[1] - min[ J_AXIS ] ) / Ly;


	*result = ( 1 - y ) + ( cos( M_PI * x ) * sin( M_PI * y ) ) / 100 ;
}

void StgFEM_StandardConditionFunctions_ConstantVector( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*		context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	
	result[0] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"ConstantValueX", 0.0  );
	result[1] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"ConstantValueY", 0.0 );
  if (context->dim == 3  ) 
    result[2] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"ConstantValueZ", 0.0 );
}

/* 3D spec ridge top BC (for milestone 1 of magma project ) 
 * to be applied to the top x-z plane of the domain */
void StgFEM_StandardConditionFunctions_SpecRidge3D( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*			feMesh             = NULL;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
        double*                 coord;

	double			leftVal;
	double			rightVal;
	double			xOffset1;
	double			xOffset2;
	double			yOffset1, yOffset2;
	double			xBegin, xEnd;
	double			zBegin, zEnd;


	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	feMesh     = feVariable->feMesh;
	coord      = Mesh_GetVertex( feMesh, node_lI );

	leftVal = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SpecRidge3DLeftSide", 0.0  );
	rightVal = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SpecRidge3DRightSide", 0.0  );
	xOffset1 = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SpecRidge3DXOffset1", 0.0  );
	xOffset2 = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SpecRidge3DXOffset2", 0.0  );
	yOffset1 = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SpecRidge3DZOffset1", 0.0  );
	yOffset2 = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SpecRidge3DZOffset2", 0.0  );
	xBegin = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SpecRidge3DXBegin", 0.0  );
	xEnd = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SpecRidge3DXEnd", 0.0  );
	zBegin = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SpecRidge3DZBegin", 0.0  );
	zEnd = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SpecRidge3DZEnd", 0.0 );

	if( coord[0] < xBegin || coord[0] > xEnd ||
	    coord[2] < zBegin || coord[2] > zEnd )
	{
		*result = 0.0;
	}
	else if( coord[0] < xOffset1 )
		*result = leftVal;
	else if( coord[0] < xOffset2 && coord[2] > yOffset1 && coord[2] < yOffset2  )
		*result = leftVal;
	else
		*result = rightVal;
}

void StgFEM_StandardConditionFunctions_TemperatureProfile( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
  FiniteElementContext *	context            = (FiniteElementContext*)_context;
  FeVariable*             feVariable         = NULL;
  FeMesh*     mesh               = NULL;
  Dictionary*             dictionary         = context->dictionary;
  double*                 result             = (double*) _result;
  double*                 coord;
  double                  T_0, H_0, dH, H, H_m, A, B, C, x_min, x_max, y_max, T_m, xc, dum;
  /* G.Ito 10/08 added variables x_min, x_max, T_m, Xc, to do variation in x
     and limit maximum T */
  
  feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
  mesh       = feVariable->feMesh;
  coord      = Mesh_GetVertex( mesh, node_lI );
  
  T_0 = Dictionary_GetDouble_WithDefault( dictionary, "TemperatureProfileTop", 0.0 );
  T_m = Dictionary_GetDouble_WithDefault( dictionary, "TemperatureProfileMax", 10000.0 );
  H_0 = Dictionary_GetDouble_WithDefault( dictionary, "TemperatureProfileH0", -1.0 );
  H_m = Dictionary_GetDouble_WithDefault( dictionary, "TemperatureProfileHm", 1.0e+8 );
  dH = Dictionary_GetDouble_WithDefault( dictionary, "TemperatureProfiledH", 0.0 );     
  A = Dictionary_GetDouble_WithDefault( dictionary, "TemperatureProfileLinearCoefficient", 0.0 );
  B = Dictionary_GetDouble_WithDefault( dictionary, "TemperatureProfileExponentialCoefficient1", 0.0 );
  C = Dictionary_GetDouble_WithDefault( dictionary, "TemperatureProfileExponentialCoefficient2", 0.0 );
  y_max = Dictionary_GetDouble_WithDefault( dictionary, "maxY", 0.0 );
  x_max = Dictionary_GetDouble_WithDefault( dictionary, "maxX", 0.0 );
  x_min = Dictionary_GetDouble_WithDefault( dictionary, "minX", 0.0 );
  xc = Dictionary_GetDouble_WithDefault( dictionary, "ExtensionCentreX", 0.0 );
  
  if (H_0<0.0)
    {
      if(coord[1]>y_max)
        {
          *result=T_0;
        }
      else
        {
          *result=T_0 + A*(y_max-coord[1]) + B*(1-exp(-C*(y_max-coord[1])));
        }
    }
  else
    {
      if(coord[1]>=y_max)
        {
          *result=T_0;
        }
      else
        {
          H=H_0 + 2*fabs(coord[0]-xc)/(x_max-x_min)*dH;
          if (H>H_m) H=H_m;
          
          dum=T_0 + ((T_m-T_0)/H)*(y_max-coord[1])
            + B*(1-exp(-C*(y_max-coord[1])));
          if (dum>T_m) dum=T_m;
          *result=dum;
        }
    }

}

void StgFEM_StandardConditionFunctions_ERF( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*			feMesh               = NULL;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
        double*                 coord;
	double                  width, scale, dilate, offset, constant;
	unsigned		dim;

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	feMesh       = feVariable->feMesh;
	coord      = Mesh_GetVertex( feMesh, node_lI );

	width = Dictionary_GetDouble_WithDefault( dictionary, "ERFWidth", 0.0 );
	offset= Dictionary_GetDouble_WithDefault(dictionary, "ERFOffset",0.0 );
	constant=Dictionary_GetDouble_WithDefault(dictionary,"ERFConstant",0.0);
        scale = Dictionary_GetDouble_WithDefault( dictionary, "ERFScale", 1.0 );
	dilate = Dictionary_GetDouble_WithDefault( dictionary,"ERFDilate",1.0 );
	dim = Dictionary_GetUnsignedInt_WithDefault( dictionary, "ERFDim", 0 );

        if(dim==3)
          {
            dim=0;
            coord=&(context->currentTime);
          }

        if(coord[dim]+offset < -width && width!=0)
          *result=constant-scale;
        else if(coord[dim]+offset > width && width!=0)
          *result=constant+scale;
        else
          *result=constant+scale*erf((coord[dim]+offset)/dilate);
}

void StgFEM_StandardConditionFunctions_ERFC(Node_LocalIndex node_lI,
                                            Variable_Index var_I,
                                            void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*			feMesh               = NULL;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
        double*                 coord;
	double                  width, scale, dilate, offset, constant;
	unsigned		dim;

	feVariable = (FeVariable*)FieldVariable_Register_GetByName
          ( context->fieldVariable_Register, "VelocityField" );
	feMesh       = feVariable->feMesh;
	coord      = Mesh_GetVertex( feMesh, node_lI );

	width = Dictionary_GetDouble_WithDefault(dictionary, "ERFCWidth", 0.0 );
	offset= Dictionary_GetDouble_WithDefault(dictionary, "ERFCOffset",0.0 );
	constant=Dictionary_GetDouble_WithDefault(dictionary,"ERFCConstant",0.0);
        scale = Dictionary_GetDouble_WithDefault(dictionary, "ERFCScale", 1.0 );
	dilate = Dictionary_GetDouble_WithDefault(dictionary,"ERFCDilate",1.0 );
	dim = Dictionary_GetUnsignedInt_WithDefault(dictionary, "ERFCDim", 0 );

        if(dim==3)
          {
            dim=0;
            coord=&(context->currentTime);
          }

        if(coord[dim]+offset < -width && width!=0)
          *result=constant-scale;
        else if(coord[dim]+offset > width && width!=0)
          *result=constant+scale;
        else
          *result=constant+scale*erfc((coord[dim]+offset)/dilate);
}

void StgFEM_StandardConditionFunctions_RubberSheet( Node_LocalIndex node_lI,
                                                    Variable_Index var_I,
                                                    void* _context,
                                                    void* _result )
{
  FiniteElementContext *	context            = (FiniteElementContext*)_context;
  FeVariable*             feVariable         = NULL;
  FeMesh*			feMesh               = NULL;
  Dictionary*             dictionary         = context->dictionary;
  double*                 result             = (double*) _result;
  double*                 coord;
  double                  lower_offset, upper_offset;
  double                  lower_value, upper_value, time;
  unsigned		dim;

  feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
  feMesh       = feVariable->feMesh;
  coord      = Mesh_GetVertex( feMesh, node_lI );

  lower_offset = Dictionary_GetDouble_WithDefault( dictionary,
                                                   "RubberSheetLowerOffset",
                                                   0.0 );
  upper_offset = Dictionary_GetDouble_WithDefault( dictionary,
                                                   "RubberSheetUpperOffset",
                                                   lower_offset );
  dim = Dictionary_GetUnsignedInt_WithDefault( dictionary,
                                               "RubberSheetDim", 0 );

  lower_value = Dictionary_GetDouble_WithDefault( dictionary,
                                                  "RubberSheetLowerValue",
                                                  0.0 );
  upper_value = Dictionary_GetDouble_WithDefault( dictionary,
                                                  "RubberSheetUpperValue",
                                                  0.0 );

  time=context->currentTime;

  if(coord[dim] < lower_offset + lower_value*time)
    {
      *result=lower_value;
    }
  else if(coord[dim] < upper_offset + upper_value*time)
    {
      double min[3], max[3];
      Mesh_GetGlobalCoordRange( feMesh, min, max );
      *result=lower_value + 
        (upper_value-lower_value)
        *(coord[dim] - min[dim])/(max[dim]-min[dim]);
    }
  else
    {
      *result=upper_value;
    }
}

/* get the BC's from the analytic solution as stored on the relevant FeVariable */
void StgFEM_StandardConditionFunctions_SpectralBCX( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             analyticFeVarX     = NULL;
	FeVariable*             numericFeVar       = NULL;
	double*                 result             = (double*) _result;
	/*FeMesh*			feMesh             = NULL;
        double*                 coord;
	Node_LocalIndex		analyticNodeI;
	Element_DomainIndex	analyticElement_I;
	double			analyticLocalElementCoord[3];
	FeMesh*			analyticFeMesh;
	*/
	analyticFeVarX = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "SpectralVelocityXField" );
	numericFeVar   = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	//feMesh         = numericFeVar->feMesh;
	//coord          = Mesh_GetVertex( feMesh, node_lI );

	//analyticFeMesh = analyticFeVarX->feMesh;
	//if( Mesh_SearchElements( analyticFeMesh, coord, &analyticElement_I ) ) {
	//	FeMesh_CoordGlobalToLocal( analyticFeMesh, analyticElement_I, coord, analyticLocalElementCoord );
	//	FeVariable_InterpolateWithinElement( analyticFeVarX, analyticElement_I, analyticLocalElementCoord, result );
	//}
	//else {	/* numerical solution node outside analytic mesh - just find closest point & use that */
	//	analyticNodeI  = Mesh_NearestVertex( analyticFeMesh, coord );
	//	FeVariable_GetValueAtNode( analyticFeVarX, analyticNodeI, result );
	//}
	
	FeVariable_GetValueAtNode( analyticFeVarX, node_lI, result );
}

void StgFEM_StandardConditionFunctions_SpectralBCY( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             analyticFeVarY     = NULL;
	FeVariable*             numericFeVar       = NULL;
	double*                 result             = (double*) _result;
	/*FeMesh*			feMesh             = NULL;
        double*                 coord;
	Node_LocalIndex		analyticNodeI;
	Element_DomainIndex	analyticElement_I;
	double			analyticLocalElementCoord[3];
	FeMesh*			analyticFeMesh;
	*/
	analyticFeVarY = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "SpectralVelocityYField" );
	numericFeVar   = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	//feMesh         = numericFeVar->feMesh;
	//coord          = Mesh_GetVertex( feMesh, node_lI );

	//analyticFeMesh = analyticFeVarY->feMesh;
	//if( Mesh_SearchElements( analyticFeMesh, coord, &analyticElement_I ) ) {
	//	FeMesh_CoordGlobalToLocal( analyticFeMesh, analyticElement_I, coord, analyticLocalElementCoord );
	//	FeVariable_InterpolateWithinElement( analyticFeVarY, analyticElement_I, analyticLocalElementCoord, result );
	//}
	//else {
	//	analyticNodeI  = Mesh_NearestVertex( analyticFeMesh, coord );
	//	FeVariable_GetValueAtNode( analyticFeVarY, analyticNodeI, result );
	//}
	
	FeVariable_GetValueAtNode( analyticFeVarY, node_lI, result );
}

void StgFEM_StandardConditionFunctions_SpectralBCZ( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             analyticFeVarZ     = NULL;
	FeVariable*             numericFeVar       = NULL;
	double*                 result             = (double*) _result;
	/*
	FeMesh*			feMesh             = NULL;
        double*                 coord;
	Node_LocalIndex		analyticNodeI;
	Element_DomainIndex	analyticElement_I;
	double			analyticLocalElementCoord[3];
	FeMesh*			analyticFeMesh;
	*/
	analyticFeVarZ = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "SpectralVelocityZField" );
	numericFeVar   = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	//feMesh         = numericFeVar->feMesh;
	//coord          = Mesh_GetVertex( feMesh, node_lI );

	//analyticFeMesh = analyticFeVarZ->feMesh;
	//if( Mesh_SearchElements( analyticFeMesh, coord, &analyticElement_I ) ) {
	//	FeMesh_CoordGlobalToLocal( analyticFeMesh, analyticElement_I, coord, analyticLocalElementCoord );
	//	FeVariable_InterpolateWithinElement( analyticFeVarZ, analyticElement_I, analyticLocalElementCoord, result );
	//}
	//else {
	//	analyticNodeI  = Mesh_NearestVertex( analyticFeMesh, coord );
	//	FeVariable_GetValueAtNode( analyticFeVarZ, analyticNodeI, result );
	//}
	
	FeVariable_GetValueAtNode( analyticFeVarZ, node_lI, result );
}

void StgFEM_StandardConditionFunctions_SpectralPressureBCX( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             analyticFeVarX     = NULL;
	FeVariable*             numericFeVar       = NULL;
	FeMesh*			feMesh             = NULL;
	double*                 result             = (double*) _result;
        double*                 coord;

	analyticFeVarX = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "SpectralPressureField" );
	numericFeVar   = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "PressureField" );
	feMesh         = numericFeVar->feMesh;
	coord          = Mesh_GetVertex( feMesh, node_lI );

	FeVariable_GetValueAtNode( analyticFeVarX, node_lI, result );
}

void StgFEM_StandardConditionFunctions_SpectralPressureBCY( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             analyticFeVarY     = NULL;
	FeVariable*             numericFeVar       = NULL;
	FeMesh*			feMesh             = NULL;
	double*                 result             = (double*) _result;
        double*                 coord;

	analyticFeVarY = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "SpectralPressureField" );
	numericFeVar   = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "PressureField" );
	feMesh         = numericFeVar->feMesh;
	coord          = Mesh_GetVertex( feMesh, node_lI );

	FeVariable_GetValueAtNode( analyticFeVarY, node_lI, result );
}

/* error function for use in 3D spec ridge top BC */
double errorFunction( double z, int n ) {
	double		pi	= 3.1415926535;
	double 		a;
	double 		erf 	= 0.0;
	int		denom;
	int		i, j;

	a = 2.0/sqrt( pi );

	for( i=0 ; i<n ; i++ ) {
		denom = 1;
		for( j=1 ; j<=2*i+1 ; j+=2 ) 
			denom *= j; 
		
		erf += pow( 2, i )*pow( z, 2*i+1 )/denom;
	}

	return erf *= a*exp( -1.0*z*z );
}



void StgFEM_StandardConditionFunctions_ErrorFunc( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*			feMesh             = NULL;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
        double*                 coord;
	double			dilate;
	double			width; 

	feVariable  = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	feMesh      = feVariable->feMesh;
	coord       = Mesh_GetVertex( feMesh, node_lI );

	dilate      = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"ErrorFuncDilate", 0.0  );
	width       = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"ErrorFuncWidth", 0.0 );

	if( coord[0] < -1.0*width ) {
		*result = -1.0;
	}
	else if( coord[0] > width  ) {
		*result = 1.0;
	}
	else {
		*result     = errorFunction( coord[0]/dilate, 5 );	
	}
}

void StgFEM_StandardConditionFunctions_GaussianDistribution( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             feVariable         = NULL;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	Name			variableName;
	double*			coord;
	unsigned		nDims              = context->dim;
	unsigned		dim_I;
	double			orig[3];
	double			sigma              = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"sigma", 1.0  );
	double			gaussianScale      = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"GaussianScale", 1.0  );
	double			background         = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"backgroundValue", 1.0  );
	double			distsq             = 0.0;

	variableName = Dictionary_GetString_WithDefault( dictionary, "FieldVariable", "" );
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, variableName );
	coord = Mesh_GetVertex( feVariable->feMesh, node_lI );

	orig[0] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"x0", 0.0  );
	orig[1] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"y0", 0.0  );
	orig[2] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"z0", 0.0 );

	for( dim_I = 0; dim_I < nDims; dim_I++ )
		distsq += ( coord[dim_I] - orig[dim_I] ) * ( coord[dim_I] - orig[dim_I] );

	*result = gaussianScale * exp( -distsq / ( 2.0 * sigma * sigma )  ) + background;
}

void StgFEM_StandardConditionFunctions_GravitationalPotential( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             feVariable         = NULL;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	Name			variableName;
	double*			coord;

	variableName = Dictionary_GetString_WithDefault( dictionary, "FieldVariable", "" );
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, variableName );
	coord = Mesh_GetVertex( feVariable->feMesh, node_lI );

	*result = -1.0 * coord[J_AXIS];
}

void StgFEM_StandardConditionFunctions_1DGaussianDistribution( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             feVariable         = NULL;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	Name			variableName;
	double*			coord;
	double			orig[3];
	double			sigma              = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"sigma", 1.0  );
	double			gaussianScale      = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"GaussianScale", 1.0  );
	double			background         = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"backgroundValue", 1.0  );
	double			distsq             = 0.0;

	variableName = Dictionary_GetString_WithDefault( dictionary, "FieldVariable", "" );
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, variableName );
	coord = Mesh_GetVertex( feVariable->feMesh, node_lI );

	orig[0] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"x0", 0.0  );
	orig[1] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"y0", 0.0  );
	orig[2] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"z0", 0.0 );

	distsq = ( coord[J_AXIS] - orig[J_AXIS] ) * ( coord[J_AXIS] - orig[J_AXIS] );

	*result = gaussianScale * exp( -distsq / ( 2.0 * sigma * sigma )  ) + background;
}

void StgFEM_StandardConditionFunctions_HalfContainer( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             feVariable         = NULL;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	Name			variableName;
	double*			coord;
	double			halfPoint          = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"halfPoint", 0.0  );

	variableName = Dictionary_GetString_WithDefault( dictionary, "FieldVariable", "" );
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, variableName );
	coord = Mesh_GetVertex( feVariable->feMesh, node_lI );

	if( coord[1] < halfPoint )
		*result = 1;
	else
		*result = 0;	
}

void StgFEM_StandardConditionFunctions_ConstantValue( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double			value              = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"constantValue", 1.0  );

	*result = value;
}

void StgFEM_StandardConditionFunctions_DiagonalLine( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double			width              = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"lineWidth", 1.0  );
	double*			coord;
	Name			variableName;
	FeVariable*             feVariable         = NULL;

	variableName = Dictionary_GetString_WithDefault( dictionary, "FieldVariable", "" );
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, variableName );
	coord = Mesh_GetVertex( feVariable->feMesh, node_lI );

	if( fabs( coord[0] - coord[1] ) < width )
		*result = 1.0;
	else
		*result = 0.0;
}

void StgFEM_StandardConditionFunctions_DeltaFunction( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double			epsilon		   = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"deltaFunctionEpsilon", 0.001  );
	unsigned		dim		   = Dictionary_GetUnsignedInt_WithDefault( dictionary, "deltaFunctionDim", 0 );
	double			centre		   = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"deltaFunctionCentre", 0.5  );
	double			value		   = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"deltaFunctionValue", 1.0  );
	double*			coord;
	Name			variableName	   = Dictionary_GetString_WithDefault( dictionary, "DeltaFunctionFeVariable", "" );
	FeVariable*		feVariable	   = (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, variableName );

	coord = Mesh_GetVertex( feVariable->feMesh, node_lI );
	
	*result = (fabs( coord[dim] - centre ) < epsilon) ? value : 0.0;
}

void StgFEM_StandardConditionFunctions_InflowBottom( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	FeVariable* feVariable;
	Dictionary*             dictionary         = context->dictionary;
	FeMesh*			mesh               = NULL;
	double*                 result             = (double*) _result;
	double sideLength, wallLength, sideV;
	double min[3], max[3];
	
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;
	Mesh_GetGlobalCoordRange( mesh, min, max );
	sideLength = max[1] - min[1];
	wallLength = max[0] - min[0];
	sideV = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"InflowSideVelocity", 1.0  );

	*result = 2.0 * sideV * sideLength / wallLength;
}

void StgFEM_StandardConditionFunctions_GaussianTube( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) 
{
        DomainContext*  context = (DomainContext*)_context;
        Dictionary*             dictionary         = context->dictionary;
        FeVariable*             feVariable = NULL;
        FeMesh*                 feMesh = NULL;
        unsigned                nDims;
        double*                 result = (double*) _result;
        double                  a1,b1,c1, a2,b2,c2, x,y,z,r_y,r_yz;
        double*                 coord;
        double                  min[3], max[3];
	double                  y_shift, z_shift;

        feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
        feMesh       = feVariable->feMesh;

        nDims = Mesh_GetDimSize( feMesh );
        Mesh_GetGlobalCoordRange( feMesh, min, max );

 	a1 = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"GaussianTube_a1", 1.0  ); /* Scales the magnitude of the perturbation. */
	c1 = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"GaussianTube_c1", 0.1  ); /* Controls the smoothing length. Smaller values produce less smoothing. */

        a2 = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"GaussianTube_a2", 0.05  ); /* Controls ampltude of oscillations */
        b2 = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"GaussianTube_b2", 6.28318530718  ); /* Controls frequency of oscillations */
	c2 = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"GaussianTube_c2", 1.570796326795  ); /* Shifts oscillations */

	y_shift = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"GaussianTube_y_origin", 0.0  );
	z_shift = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"GaussianTube_z_origin", 0.0  );

        coord = Mesh_GetVertex( feMesh, node_lI );

	x = coord[ I_AXIS ];
	y = coord[ J_AXIS ];

	y = y - y_shift;
	if (nDims==2) {
		b1 = a2 * sin( b2*x - c2 );
		r_y = sqrt( (y-b1)*(y-b1) );
		*result = a1 * exp( -(r_y * r_y) / (2.0*c1*c1) );
	}
	if (nDims==3) {
		z = coord[ K_AXIS ];
		z = z - z_shift;

		b1 = a2 * sin( b2*x - c2 );
		r_yz = sqrt( (y-b1)*(y-b1) + z*z );
		*result = a1 * exp( -(r_yz * r_yz) / (2.0*c1*c1) );
	}


}


void StgFEM_StandardConditionFunctions_WarsTemperature( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) 
{
  FiniteElementContext *	context            = (FiniteElementContext*)_context;
  FeVariable*             feVariable         = NULL;
  FeMesh*     mesh               = NULL;
  Dictionary*             dictionary         = context->dictionary;
  double*                 result             = (double*) _result;
  double*                 coord;
  double                  EAEnd, WarsStart, WarsHeight, WarsTTop,
    WarsTBottom, h, maxY;
  
  feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
  mesh       = feVariable->feMesh;
  coord      = Mesh_GetVertex( mesh, node_lI );
  
  EAEnd = Dictionary_GetDouble( dictionary, "EAEnd");
  WarsStart = Dictionary_GetDouble( dictionary, "WarsStart");
  WarsHeight = Dictionary_GetDouble( dictionary, "WarsHeight");
  WarsTTop = Dictionary_GetDouble( dictionary, "WarsTTop");
  WarsTBottom = Dictionary_GetDouble( dictionary, "WarsTBottom");     
  maxY=Dictionary_GetDouble( dictionary, "maxY");


  h=WarsHeight*(coord[0]-EAEnd)/(WarsStart-EAEnd);
  if(coord[0]<EAEnd)
    h=0;
  if(coord[0]>WarsStart)
    h=WarsHeight;
  *result=WarsTBottom + ((coord[1]-h)/(maxY-h))*(WarsTTop-WarsTBottom);
}

void StgFEM_StandardConditionFunctions_Quadratic( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) 
{
  FiniteElementContext *	context            = (FiniteElementContext*)_context;
  FeVariable*             feVariable         = NULL;
  FeMesh*     mesh               = NULL;
  Dictionary*             dictionary         = context->dictionary;
  double*                 result             = (double*) _result;
  double*                 coord;
  int                     dim;
  double                  a, b, c;
  
  feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
  mesh       = feVariable->feMesh;
  coord      = Mesh_GetVertex( mesh, node_lI );
  
  dim = Dictionary_GetInt( dictionary, "Quadratic_Dim");
  a = Dictionary_GetDouble( dictionary, "Quadratic_Constant");
  b = Dictionary_GetDouble( dictionary, "Quadratic_Linear");
  c = Dictionary_GetDouble( dictionary, "Quadratic_Quadratic");

  *result= a + coord[dim]*(b + c*coord[dim]);
}

int Binary_Search(double *data, int s, int e, double value);

void StgFEM_StandardConditionFunctions_File1( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) 
{
  StgFEM_StandardConditionFunctions_FileN(node_lI,var_I,_context,_result,1);
}

void StgFEM_StandardConditionFunctions_File2( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) 
{
  StgFEM_StandardConditionFunctions_FileN(node_lI,var_I,_context,_result,2);
}

void StgFEM_StandardConditionFunctions_File3( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) 
{
  StgFEM_StandardConditionFunctions_FileN(node_lI,var_I,_context,_result,3);
}

void StgFEM_StandardConditionFunctions_File4( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) 
{
  StgFEM_StandardConditionFunctions_FileN(node_lI,var_I,_context,_result,4);
}

void StgFEM_StandardConditionFunctions_File5( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) 
{
  StgFEM_StandardConditionFunctions_FileN(node_lI,var_I,_context,_result,5);
}

void StgFEM_StandardConditionFunctions_File6( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) 
{
  StgFEM_StandardConditionFunctions_FileN(node_lI,var_I,_context,_result,6);
}

void StgFEM_StandardConditionFunctions_File7( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) 
{
  StgFEM_StandardConditionFunctions_FileN(node_lI,var_I,_context,_result,7);
}

void StgFEM_StandardConditionFunctions_File8( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) 
{
  StgFEM_StandardConditionFunctions_FileN(node_lI,var_I,_context,_result,8);
}

void StgFEM_StandardConditionFunctions_File9( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) 
{
  StgFEM_StandardConditionFunctions_FileN(node_lI,var_I,_context,_result,9);
}

void StgFEM_StandardConditionFunctions_File10( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) 
{
  StgFEM_StandardConditionFunctions_FileN(node_lI,var_I,_context,_result,10);
}

void StgFEM_StandardConditionFunctions_FileN( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result, int file_num )
{
  FiniteElementContext *	context            = (FiniteElementContext*)_context;
  FeVariable*             feVariable         = NULL;
  FeMesh*     mesh               = NULL;
  Dictionary*             dictionary         = context->dictionary;
  double*                 result             = (double*) _result;
  double*                 coord;
  int                     dim, i;
  char *filename;
  int N;
  int result_index;
  double factor;
  static double *coords=NULL;
  static double *data=NULL;
  feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
  mesh       = feVariable->feMesh;
  coord      = Mesh_GetVertex( mesh, node_lI );
  
  char fileN_number[10], fileN_dim[15], fileN_name[15], fileN_N[15];
  sprintf(fileN_number,"File%d",file_num);
  sprintf(fileN_dim,"File%d_Dim",file_num);
  sprintf(fileN_name,"File%d_Name",file_num);
  sprintf(fileN_N,"File%d_N",file_num);

  dim = Dictionary_GetInt( dictionary, fileN_dim);
  filename = Dictionary_GetString( dictionary, fileN_name);
  N = Dictionary_GetInt( dictionary, fileN_N);

  Journal_Firewall(dim>=0 && dim<3,
                   Journal_Register( Error_Type,"StgFEM_StandardConditionFunctions_FileN"),
                   "%s must be either 0, 1, or 2, but was set to %d\n",
                   fileN_dim,dim);
  Journal_Firewall(N>0,
                   Journal_Register( Error_Type,"StgFEM_StandardConditionFunctions_FileN"),
                   "%s must be greater than zero, but was set to %d.\n",
                   fileN_N,N);
  if(data==NULL)
    {
      
      FILE *fp=fopen(filename,"r");
      Journal_Firewall(fp!=NULL,
                       Journal_Register( Error_Type,"StgFEM_StandardConditionFunctions_FileN"),
                       "Bad filename for %s.  Could not open %s\n",
                       fileN_name,filename);
      data=(double *)malloc(N*sizeof(double));
      coords=(double *)malloc(N*sizeof(double));

      Journal_Firewall(data!=NULL && coords!=NULL,
                       Journal_Register( Error_Type,"StgFEM_StandardConditionFunctions_FileN"),
                       "Could not allocate enough memory for %s\n",file_num);
      for(i=0;i<N;++i)
        fscanf(fp,"%lf %lf",coords+i,data+i);
    }

  Journal_Firewall(!(coord[dim]<coords[0] || coord[dim]>coords[N-1]),
                   Journal_Register( Error_Type,"StgFEM_StandardConditionFunctions_FileN"),
                   "The range in the file '%s' does not cover this value %g\nIt only covers %g to %g.\n",
                   filename,coord[dim],coords[0],coords[1]);

  result_index=Binary_Search(coords,0,N-1,coord[dim]);
  factor=(coords[result_index+1]-coord[dim])
    / (coords[result_index+1]-coords[result_index]);
  
  *result=data[result_index]*factor + data[result_index+1]*(1-factor);
}


int Binary_Search(double *data, int s, int e, const double value)
{
  int start, end, midpoint;

  start=s;
  end=e;
  midpoint=e;
  
  midpoint=(end-start)/2 + start;
  while(start!=midpoint)
    {
      if(data[midpoint]>=value)
        end=midpoint;
      else
        start=midpoint;
      midpoint=(end-start)/2 + start;
    }
  return start;
}
