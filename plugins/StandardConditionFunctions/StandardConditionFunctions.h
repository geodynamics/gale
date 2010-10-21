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

void StgFEM_StandardConditionFunctions_SolidBodyRotation( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_PartialRotationX( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_PartialRotationY( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_TaperedRotationX( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_TaperedRotationY( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_SimpleShear( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_SimpleShearInverted( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_ShearZ( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_Extension( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_PartialLid_TopLayer( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* result ) ;
void StgFEM_StandardConditionFunctions_LinearInterpolationLid( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* result ) ;
void StgFEM_StandardConditionFunctions_Lid_RampWithCentralMax( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* result ) ;


void StgFEM_StandardConditionFunctions_SinusoidalLid( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* result ) ;
void StgFEM_StandardConditionFunctions_CornerOnly( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* result ) ;
void StgFEM_StandardConditionFunctions_LinearWithSinusoidalPerturbation( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_Trigonometry( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void Stg_FEM_VelicTemperatureIC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void Stg_FEM_VelicTemperatureIC_SolB( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_AnalyticalTemperatureIC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_TemperatureCosineHill( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_SinusoidalExtension( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_EdgeDriveConvectionIC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_ThermalEdgeDriveConvectionIC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_StepFunction( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StG_FEM_StandardConditionFunctions_StepFunctionProduct1( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StG_FEM_StandardConditionFunctions_StepFunctionProduct2( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StG_FEM_StandardConditionFunctions_StepFunctionProduct3( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StG_FEM_StandardConditionFunctions_StepFunctionProduct4( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_TemperatureProfile( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StG_FEM_StandardConditionFunctions_Gaussian( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_MovingStepFunction( Node_LocalIndex nodeInd, Variable_Index varInd, void* _ctx, void* _result );
void StgFEM_StandardConditionFunctions_SpecRidge3D( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_SpectralBCX( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_SpectralBCY( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_SpectralBCZ( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_SpectralPressureBCX( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_SpectralPressureBCY( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_ErrorFunc( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;

void StgFEM_StandardConditionFunctions_LinearVelocityLeftWall( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* result );
void StgFEM_StandardConditionFunctions_LinearVelocityRightWall( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* result );

void StgFEM_StandardConditionFunctions_ConvectionBenchmark( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
void StgFEM_StandardConditionFunctions_ConstantVector( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_ConstantVelocity( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_ERF( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_ERFC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_RubberSheet( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );

void StgFEM_StandardConditionFunctions_GaussianDistribution( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_GravitationalPotential( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_1DGaussianDistribution( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_HalfContainer( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_ConstantValue( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_DiagonalLine( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_DeltaFunction( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );

void StgFEM_StandardConditionFunctions_InflowBottom( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_GaussianTube( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_WarsTemperature( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_Quadratic( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_File1( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_File2( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_File3( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_File4( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_File5( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_File6( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_File7( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_File8( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_File9( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );
void StgFEM_StandardConditionFunctions_File10( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result );


#endif	
