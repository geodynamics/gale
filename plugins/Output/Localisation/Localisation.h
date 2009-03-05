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

/** this plugin calculates the localisation factor of the deformation.  this is defined as the 
    ratio of the area over which 80% (by default) of deformation occurs, to the total area.
    it is calculated across the top layer of the domain.   the localisation value may give
    an indication of plate-like behavior, which typically exhibits a large proportion of 
    deformation across small areas (and so small localisation values).
    
    this value is calculated by first constructing a field consisting of the required strain-rate
    components for the given surface.  this field is integrated (note that integration functions 
    utilised here have been adapted from those found in FeVariable.c), with the contributing
    integrand evaluations stored, then qsorted (according to the evaluation value).  we then sum
    this ordered list (and also corresponding integration point weights), starting from largest value, 
    and continuing until some factor (0.8 by default) of the total integral is reached.  at this 
    point, the ratio of the total area summed (using the integration point weights) to the total
    area is the localisation factor.
*/

#ifndef __UnderworldRheology_Localisation_h__
#define __UnderworldRheology_Localisation_h__

	extern const Type Underworld_Localisation_Type;

	typedef struct {
		double integrand;
		double weight;
	} IntegrandWeightStruct;


	typedef struct {
		__Codelet
		OperatorFeVariable*    reducedStrainRateFieldInvariantRoot;
		double                 Localisation;
		double                 deformationFactor;
		IntegrandWeightStruct* integrandWeight;
	} Underworld_Localisation;

	Index Underworld_Localisation_Register( PluginsManager* pluginsManager );
	void Underworld_Localisation_Setup( UnderworldContext* context ) ;
	void Underworld_Localisation_Output( UnderworldContext* context ) ;
	void Underworld_Localisation_SymmetricTensor_LowerDimension_InvariantRoot_2d( void* operatorObject, double* operand0, double* result );
	void Underworld_Localisation_SymmetricTensor_LowerDimension_InvariantRoot_3d( void* operatorObject, double* operand0, double* result );

	double Localisation_IntegrateLayer_AxisIndependent( 
		void* feVariable, void* _swarm,
		Axis layerAxis, Index layerIndex, Dimension_Index dim, 
		Axis axis0, Axis axis1, Axis axis2, void* Localisation );
	double Localisation_IntegratePlane( void* feVariable, Axis planeAxis, double planeHeight, void* Localisation);
	double Localisation_IntegrateElement_AxisIndependent( 
			void* feVariable, void* _swarm,
			Element_DomainIndex dElement_I, Dimension_Index dim, 
			Axis axis0, Axis axis1, Axis axis2, void* Localisation ) ;
	unsigned Localisation_FindTotElements( void* feVariable, Axis layerAxis, double planeHeight  );
	int _QsortIntegrand( const void* _a, const void* _b );
	void Localisation_Create_MPI_Datatype();	
#endif
