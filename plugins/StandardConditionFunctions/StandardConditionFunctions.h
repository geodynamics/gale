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
*/
/** \file
** Role:
**	Defines any header info, such as new structures, needed by this plugin
**
** Assumptions:
**
** Comments:
**
** $Id $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_StandardConditionFunctions_h__
#define __StgFEM_StandardConditionFunctions_h__

extern const Type StgFEM_StandardConditionFunctions_Type;

typedef struct {
	__Codelet
} StgFEM_StandardConditionFunctions;

Index StgFEM_StandardConditionFunctions_Register( PluginsManager* pluginsManager );

void StgFEM_StandardConditionFunctions_SolidBodyRotation(const double *coord, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_PartialRotationX(const double *coord, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_PartialRotationY(const double *coord, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_TaperedRotationX(const double *coord, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_TaperedRotationY(const double *coord, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_SimpleShear(const double *coord, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_SimpleShearInverted(const double *coord, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_ShearZ(const double *coord, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_Extension(const double *coord, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_PartialLid_TopLayer(const double *coord, void* _context, void* result ) ;
void StgFEM_StandardConditionFunctions_LinearInterpolationLid(const double *coord, void* _context, void* result ) ;
void StgFEM_StandardConditionFunctions_Lid_RampWithCentralMax(const double *coord, void* _context, void* result ) ;


void StgFEM_StandardConditionFunctions_SinusoidalLid(const double *coord, void* _context, void* result ) ;
void StgFEM_StandardConditionFunctions_LinearWithSinusoidalPerturbation(const double *coord, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_Trigonometry(const double *coord, void* _context, void* _result ) ;
void Stg_FEM_VelicTemperatureIC(const double *coord, void* _context, void* _result ) ;
void Stg_FEM_VelicTemperatureIC_SolB(const double *coord, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_AnalyticalTemperatureIC(const double *coord, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_TemperatureCosineHill(const double *coord, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_SinusoidalExtension(const double *coord, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_EdgeDriveConvectionIC(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_ThermalEdgeDriveConvectionIC(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_StepFunction(const double *coord, void* _context, void* _result ) ;
void StG_FEM_StandardConditionFunctions_StepFunctionProduct1(const double *coord, void* _context, void* _result ) ;
void StG_FEM_StandardConditionFunctions_StepFunctionProduct2(const double *coord, void* _context, void* _result ) ;
void StG_FEM_StandardConditionFunctions_StepFunctionProduct3(const double *coord, void* _context, void* _result ) ;
void StG_FEM_StandardConditionFunctions_StepFunctionProduct4(const double *coord, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_TemperatureProfile(const double *coord, void* _context, void* _result ) ;
void StG_FEM_StandardConditionFunctions_Gaussian(const double *coord, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_MovingStepFunction(const double *coord, void* _ctx, void* _result );
void StgFEM_StandardConditionFunctions_SpecRidge3D(const double *coord, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_ErrorFunc(const double *coord, void* _context, void* _result ) ;

void StgFEM_StandardConditionFunctions_LinearVelocityLeftWall(const double *coord, void* _context, void* result );
void StgFEM_StandardConditionFunctions_LinearVelocityRightWall(const double *coord, void* _context, void* result );

void StgFEM_StandardConditionFunctions_ConvectionBenchmark(const double *coord, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_ConstantVector(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_ConstantVelocity(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_ERF(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_ERFC(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_RubberSheet(const double *coord, void* _context, void* _result );

void StgFEM_StandardConditionFunctions_GaussianDistribution(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_GravitationalPotential(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_1DGaussianDistribution(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_HalfContainer(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_ConstantValue(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_DiagonalLine(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_DeltaFunction(const double *coord, void* _context, void* _result );

void StgFEM_StandardConditionFunctions_InflowBottom(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_GaussianTube(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_WarsTemperature(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_Quadratic(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_File1(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_File2(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_File3(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_File4(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_File5(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_File6(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_File7(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_File8(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_File9(const double *coord, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_File10(const double *coord, void* _context, void* _result );

#endif	
