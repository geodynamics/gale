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
** $Id: MohrCoulomb.c 11997 2008-05-22 18:16:17Z walter $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "RheologyClass.h"
#include "StrainWeakening.h"
#include "YieldRheology.h"
#include "MohrCoulomb.h"
#include "ConstitutiveMatrix.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type MohrCoulomb_Type = "MohrCoulomb";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
MohrCoulomb* _MohrCoulomb_New(  MOHRCOULOMB_DEFARGS  ) 
{
	MohrCoulomb*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(MohrCoulomb) );
	self = (MohrCoulomb*) _YieldRheology_New(  YIELDRHEOLOGY_PASSARGS  );
	
	/* Function pointers for this class that are not on the parent class should be set here */
	
	return self;
}

void _MohrCoulomb_Init(    
		MohrCoulomb*          self,
		FeVariable*           pressureField,
		FeVariable*           strainRateField,
		SwarmVariable*        swarmStrainRate,
		MaterialPointsSwarm*  materialPointsSwarm,
		double                cohesion,
		double                cohesionAfterSoftening,
		double                frictionCoefficient,
		double                frictionCoefficientAfterSoftening,
		double                minimumYieldStress)
{                                                       
	self->materialPointsSwarm     = materialPointsSwarm;
	self->pressureField           = pressureField;
	self->strainRateField  = strainRateField;
	
	self->cohesion = cohesion;
	self->frictionCoefficient = frictionCoefficient;
	
	/* Strain softening of Cohesion - (linear weakening is assumed) */
	/* needs a softening factor between +0 and 1 and a reference strain > 0 */
	self->cohesionAfterSoftening = cohesionAfterSoftening;
	
	/* Strain softening of Friction - (linear weakening is assumed) */
	/* needs a softening factor between +0 and 1 and a reference strain > 0 */
	self->frictionCoefficientAfterSoftening = frictionCoefficientAfterSoftening;

	self->minimumYieldStress = minimumYieldStress;
}

void* _MohrCoulomb_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(MohrCoulomb);
	Type                                                             type = MohrCoulomb_Type;
	Stg_Class_DeleteFunction*                                     _delete = _YieldRheology_Delete;
	Stg_Class_PrintFunction*                                       _print = _YieldRheology_Print;
	Stg_Class_CopyFunction*                                         _copy = _YieldRheology_Copy;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _MohrCoulomb_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _MohrCoulomb_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _MohrCoulomb_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _MohrCoulomb_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _YieldRheology_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _MohrCoulomb_Destroy;
	Rheology_ModifyConstitutiveMatrixFunction*  _modifyConstitutiveMatrix = _MohrCoulomb_ModifyConstitutiveMatrix;
	YieldRheology_GetYieldCriterionFunction*           _getYieldCriterion = _MohrCoulomb_GetYieldCriterion;
	YieldRheology_GetYieldIndicatorFunction*           _getYieldIndicator = _MohrCoulomb_GetYieldIndicator;
	YieldRheology_HasYieldedFunction*                         _hasYielded = _MohrCoulomb_HasYielded;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _MohrCoulomb_New(  MOHRCOULOMB_PASSARGS  );
}

void _MohrCoulomb_AssignFromXML( void* rheology, Stg_ComponentFactory* cf,
                             void *data ){
	MohrCoulomb*           self = (MohrCoulomb*)rheology;
	FeVariable*            pressureField;
	MaterialPointsSwarm*   materialPointsSwarm;
	FeVariable*            strainRateField;
	SwarmVariable*         swarmStrainRate;
	
	/* Construct Parent */
	_YieldRheology_AssignFromXML( self, cf, data );

	/* Make sure that there is strain weakening */
	Journal_Firewall(
		self->strainWeakening != NULL,
		Journal_Register( Error_Type, (Name)self->type  ),
		"Error in func '%s' for %s '%s': MohrCoulomb rheology needs strain weakening.\n", 
		__func__, self->type, self->name );
	
	materialPointsSwarm    = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"MaterialPointsSwarm", MaterialPointsSwarm, True, data  );
	pressureField          = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"PressureField", FeVariable, True, data  );
	strainRateField = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"StrainRateField", FeVariable, True, data  );

	swarmStrainRate = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"swarmStrainRate", SwarmVariable, False, data  );
	
	_MohrCoulomb_Init( 
			self,
			pressureField,
			strainRateField,
			swarmStrainRate,
			materialPointsSwarm, 
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"cohesion", 0.0  ),
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"cohesionAfterSoftening", 0.0  ),
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"frictionCoefficient", 0.0  ),
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"frictionCoefficientAfterSoftening", 0.0  ),
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"minimumYieldStress", 0.0 )  );
}

void _MohrCoulomb_Build( void* rheology, void* data ) {
	MohrCoulomb*  self               = (MohrCoulomb*) rheology;

	/* Build parent */
	_YieldRheology_Build( self, data );
	
	Stg_Component_Build(	self->materialPointsSwarm, data, False);
	Stg_Component_Build(	self->pressureField, data, False);
	Stg_Component_Build(	self->strainRateField, data, False);

}

void _MohrCoulomb_Initialise( void* rheology, void* data ) {
	MohrCoulomb*     self                  = (MohrCoulomb*) rheology;
	Particle_Index                  lParticle_I;
	Particle_Index                  particleLocalCount;

	_YieldRheology_Initialise( self, data );
	
	/* Initialise these components now just in case this function is called before their own _Initialise function */
	Stg_Component_Initialise( self->materialPointsSwarm, data, False );
	Stg_Component_Build(	self->pressureField, data, False);	
	Stg_Component_Initialise( self->strainWeakening, data, False );

	/* We don't need to Initialise hasYieldedVariable because it's a parent variable and _YieldRheology_Initialise
	 * has already been called */
	particleLocalCount = self->hasYieldedVariable->variable->arraySize;

	/* If restarting from checkpoint, don't change the parameters on the particles */
	if ( self->context->loadFromCheckPoint == False ) {
		for ( lParticle_I = 0 ; lParticle_I < particleLocalCount ; lParticle_I++ ) { 
			Variable_SetValueChar( self->hasYieldedVariable->variable, lParticle_I, False );
                }
	}
}

void _MohrCoulomb_Destroy( void* rheology, void* data ) {
	MohrCoulomb*  self               = (MohrCoulomb*) rheology;

	Stg_Component_Destroy(	self->materialPointsSwarm, data, False);
	Stg_Component_Destroy(	self->pressureField, data, False);
	Stg_Component_Destroy(	self->strainRateField, data, False);

	/* Destroy parent */
	_YieldRheology_Destroy( self, data );

}
	
void _MohrCoulomb_ModifyConstitutiveMatrix( 
		void*                                              rheology, 
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               materialPointsSwarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi )
{
	MohrCoulomb*     self                  = (MohrCoulomb*) rheology;

	/* Don't want to yield on the first ever solve */
	if ( constitutiveMatrix->previousSolutionExists == False ) {
		return;
	}

	/* Calculate and store values of stress, pressure, velocity gradients, eigenvectors so they are only calculated once */
	_MohrCoulomb_StoreCurrentParameters( self, constitutiveMatrix, materialPointsSwarm, lElement_I, materialPoint, xi );

	_YieldRheology_ModifyConstitutiveMatrix( self, constitutiveMatrix, materialPointsSwarm, lElement_I, materialPoint, xi );
}

double _MohrCoulomb_GetYieldIndicator( 
		void*                                              rheology,
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               materialPointsSwarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi )
{
	MohrCoulomb*           self             = (MohrCoulomb*) rheology;
	Eigenvector*                          eigenvectorList  = self->currentEigenvectorList;
	double                                sigma_ns;
	Dimension_Index                       dim              = constitutiveMatrix->dim;
	double                                stressMin;
	double                                stressMax;

        /* The yield indicator is simply the difference in principal
           stresses.  In 2D, this is the 2nd invariant of the stress
           tensor. This is not the actual shear stress.  Rather, it is
           the terms with deviatoric stress, because that scales
           linearly with viscosity.  So if we scale the viscosity, we
           only affect the yield indicator. */

        stressMin = eigenvectorList[0].eigenvalue;
        stressMax = (dim == 2 ? eigenvectorList[1].eigenvalue : eigenvectorList[2].eigenvalue);
        sigma_ns = 0.5 * ( stressMax - stressMin );

        if(dim==3)
          {
            double alpha=_MohrCoulomb_EffectiveFrictionCoefficient( self, materialPoint );
            double theta = 0.5 * atan( 1.0/ alpha);
            double factor=1/(sin(2*theta) - alpha*cos(2*theta));

            sigma_ns-=(stressMax+stressMin)*0.5*alpha*factor;
          }
	
	return fabs(sigma_ns);
}

double _MohrCoulomb_GetYieldCriterion( 
		void*                                              rheology,
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               materialPointsSwarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi )
{
	MohrCoulomb*          self             = (MohrCoulomb*) rheology;
	double                               minimumYieldStress;
	double                               effectiveCohesion;
	double                               effectiveFrictionCoefficient;
	double                               frictionalStrength;

	FeVariable*                       pressureField  = self->pressureField;  
	double                            pressure;
        double theta;
        double factor;

	FeVariable_InterpolateWithinElement( pressureField, lElement_I, xi, &pressure );

	/* Calculate frictional strength.  We modify the friction and
           cohesion because we have grouped terms from the normal
           stresses and moved it to the yield indicator. */
	effectiveFrictionCoefficient = _MohrCoulomb_EffectiveFrictionCoefficient( self, materialPoint );
	effectiveCohesion            = _MohrCoulomb_EffectiveCohesion( self, materialPoint );

        theta = 0.5 * atan( 1.0/ effectiveFrictionCoefficient);

        factor=sin(2*theta);
        
        frictionalStrength = effectiveFrictionCoefficient * pressure*factor + effectiveCohesion*factor;

        /* If the minimumYieldStress is not set, then use the
           effective cohesion.  Maybe it should be the modified
           effective cohesion, though that probably should not matter
           much. */
	minimumYieldStress = self->minimumYieldStress;
        if(minimumYieldStress==0.0)
          minimumYieldStress=effectiveCohesion;
	
	/* Make sure frictionalStrength is above the minimum */
	if ( frictionalStrength < minimumYieldStress*factor) 
		frictionalStrength = minimumYieldStress*factor;
	return frictionalStrength;
}

/* This scales the viscosity so that the stress is the yield stress.
   We can simply scale the viscosity by yieldCriterion/yieldIndicator
   because yieldCriterion does not depend on the deviatoric stress
   (and thus the viscosity), while the yieldIndicator depends linearly
   on the deviatoric stresses. */

void _MohrCoulomb_HasYielded( 
		void*                                              rheology,
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               materialPointsSwarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		double                                             yieldCriterion,
		double                                             yieldIndicator )
{
  double                           viscosity = ConstitutiveMatrix_GetIsotropicViscosity( constitutiveMatrix );

  double beta=1.0-yieldCriterion/(yieldIndicator);
  double corr = -viscosity * beta;

  if( (viscosity + corr) < ((MohrCoulomb*)rheology)->minVisc )
    corr = ((MohrCoulomb*)rheology)->minVisc - viscosity;

  ConstitutiveMatrix_IsotropicCorrection( constitutiveMatrix, corr );
}

double _MohrCoulomb_EffectiveCohesion( void* rheology, void* materialPoint ) {
	MohrCoulomb*    self                          = (MohrCoulomb*) rheology;
	double                         effectiveCohesion;
	
        double strainWeakeningRatio = StrainWeakening_CalcRatio( self->strainWeakening, materialPoint );
	
        effectiveCohesion =  self->cohesion * (1.0 - strainWeakeningRatio) + 
          self->cohesionAfterSoftening * strainWeakeningRatio;

	return effectiveCohesion;
}

double _MohrCoulomb_EffectiveFrictionCoefficient( void* rheology, void* materialPoint ) {
	MohrCoulomb*    self                          = (MohrCoulomb*) rheology;
	static double                  effectiveFrictionCoefficient  = 0.0;

        double strainWeakeningRatio = StrainWeakening_CalcRatio( self->strainWeakening, materialPoint );
	
        effectiveFrictionCoefficient = self->frictionCoefficient * (1.0 - strainWeakeningRatio) +
          self->frictionCoefficientAfterSoftening * strainWeakeningRatio;	

	return effectiveFrictionCoefficient;
}

void _MohrCoulomb_StoreCurrentParameters( 
		void*                                              rheology,
		ConstitutiveMatrix*                                constitutiveMatrix, 
		MaterialPointsSwarm*                               materialPointsSwarm,
		Element_LocalIndex                                 lElement_I, 
		MaterialPoint*                                     materialPoint,
		Coord                                              xi ) 
{
	MohrCoulomb*		self = (MohrCoulomb*) rheology;
	Dimension_Index	dim = constitutiveMatrix->dim;
	Eigenvector			evectors[3];
	double				trace;
	int					i;
	
	FeVariable_InterpolateWithinElement( self->pressureField, lElement_I, xi, &self->currentPressure );
	if( !self->swarmStrainRate ) {
		FeVariable_InterpolateWithinElement( self->strainRateField, lElement_I, xi, self->currentStrainRate );
	}
	else {
	   SwarmVariable_ValueAt( self->swarmStrainRate, constitutiveMatrix->currentParticleIndex, self->currentStrainRate );
	}

	SymmetricTensor_GetTrace(self->currentStrainRate, dim, &trace);

	/* Subtract the trace (which should be zero anyway).  We can
		use TensorMapST3D even for 2D, because it is the same for
		the xx and yy components */
	for(i=0;i<dim;++i)
		self->currentStrainRate[TensorMapST3D[i][i]]-=trace/dim;

	ConstitutiveMatrix_CalculateStress( constitutiveMatrix, self->currentStrainRate, self->currentStress );
	
	SymmetricTensor_CalcAllEigenvectors( self->currentStress, dim, self->currentEigenvectorList );

	SymmetricTensor_CalcAllEigenvectors( self->currentStrainRate, dim, evectors);
}


