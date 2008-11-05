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
** $Id: $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __UnderworldRheology_Plateness_h__
#define __UnderworldRheology_Plateness_h__

	extern const Type Underworld_Plateness_Type;
	typedef struct {
		__Codelet
		OperatorFeVariable* reducedStrainRateFieldInvariantRoot;
		double      plateness;
	} Underworld_Plateness;

	Index Underworld_Plateness_Register( PluginsManager* pluginsManager );
	void Underworld_Plateness_Setup( UnderworldContext* context ) ;
	void Underworld_Plateness_Output( UnderworldContext* context ) ;
	void Underworld_Plateness_SymmetricTensor_LowerDimension_InvariantRoot_2d( void* operatorObject, double* operand0, double* result );
	void Underworld_Plateness_SymmetricTensor_LowerDimension_InvariantRoot_3d( void* operatorObject, double* operand0, double* result );

	double Plateness_IntegrateLayer_AxisIndependent( 
		void* feVariable, void* _swarm,
		Axis layerAxis, Index layerIndex, Dimension_Index dim, 
		Axis axis0, Axis axis1, Axis axis2, double** value_weight_Matrix); //, double minValue, double maxValue );
	double Plateness_IntegratePlane( void* feVariable, Axis planeAxis, double planeHeight, double** value_weight_Matrix); //, double minValue, double maxValue);
	double Plateness_IntegrateElement_AxisIndependent( 
			void* feVariable, void* _swarm,
			Element_DomainIndex dElement_I, Dimension_Index dim, 
			Axis axis0, Axis axis1, Axis axis2, double** value_weight_Matrix, double minValue, double maxValue) ;

#endif
