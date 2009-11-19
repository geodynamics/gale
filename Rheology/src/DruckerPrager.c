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
** $Id: DruckerPrager.c 743 2008-06-23 01:49:43Z JulianGiordani $
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
#include "VonMises.h"
#include "DruckerPrager.h"
#include "ConstitutiveMatrix.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type DruckerPrager_Type = "DruckerPrager";

/* Public Constructor */
DruckerPrager* DruckerPrager_New(
      Name                  name,
      AbstractContext*      context,
      StrainWeakening*      strainWeakening, 
      MaterialPointsSwarm*  materialPointsSwarm, 
      double                minVisc, 
      FeVariable*           strainRateField,
      SwarmVariable*        swarmStrainRate,
      double                cohesion,
      double                cohesionAfterSoftening,
      Bool                  strainRateSoftening,
		FeVariable*           pressureField,
		SwarmVariable*        swarmPressure,
		double                minimumYieldStress,
		double                frictionCoefficient,
		double                frictionCoefficientAfterSoftening )
{
   DruckerPrager* self = (DruckerPrager*) _DruckerPrager_DefaultNew( name );

   _Rheology_Init( self, context );
   _YieldRheology_Init( self, strainWeakening, materialPointsSwarm, minVisc ); 
   _VonMises_Init( self, strainRateField, swarmStrainRate, cohesion, cohesionAfterSoftening, strainRateSoftening );
   _DruckerPrager_Init( self, pressureField, swarmPressure, materialPointsSwarm, minimumYieldStress, frictionCoefficient, frictionCoefficientAfterSoftening );
   self->isConstructed = True;
   return self;
}

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
DruckerPrager* _DruckerPrager_New( 
		SizeT                                              sizeOfSelf,
		Type                                               type,
		Stg_Class_DeleteFunction*                          _delete,
		Stg_Class_PrintFunction*                           _print,
		Stg_Class_CopyFunction*                            _copy, 
		Stg_Component_DefaultConstructorFunction*          _defaultConstructor,
		Stg_Component_ConstructFunction*                   _construct,
		Stg_Component_BuildFunction*                       _build,
		Stg_Component_InitialiseFunction*                  _initialise,
		Stg_Component_ExecuteFunction*                     _execute,
		Stg_Component_DestroyFunction*                     _destroy,
		Rheology_ModifyConstitutiveMatrixFunction*         _modifyConstitutiveMatrix,
		YieldRheology_GetYieldCriterionFunction*           _getYieldCriterion,
		YieldRheology_GetYieldIndicatorFunction*           _getYieldIndicator,
		YieldRheology_HasYieldedFunction*                  _hasYielded,
		Name                                               name ) 
{
	DruckerPrager*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(DruckerPrager) );
	self = (DruckerPrager*) _VonMises_New( 
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
			_modifyConstitutiveMatrix,
			_getYieldCriterion,
			_getYieldIndicator,
			_hasYielded,
			name );
	
	/* Function pointers for this class that are not on the parent class should be set here */
	
	return self;
}

void _DruckerPrager_Init(
		DruckerPrager*                                     self,
		FeVariable*                                        pressureField,
		SwarmVariable*                                     swarmPressure,
		MaterialPointsSwarm*                               materialPointsSwarm,
		double                                             minimumYieldStress,
		double                                             frictionCoefficient,
		double                                             frictionCoefficientAfterSoftening )
{
	DruckerPrager_Particle*   particleExt;
	StandardParticle          materialPoint;
	
	self->particleExtHandle = ExtensionManager_Add( materialPointsSwarm->particleExtensionMgr,
		DruckerPrager_Type, sizeof(DruckerPrager_Particle) );
		
	/* Assign Pointers */
	self->pressureField       = pressureField;
	self->frictionCoefficient = frictionCoefficient;
	self->minimumYieldStress  = minimumYieldStress;
	self->swarmPressure       = swarmPressure;
	
	/* Strain softening of Friction - (linear weakening is assumed) */
	/* needs a softening factor between +0 and 1 and a reference strain > 0 */
	self->frictionCoefficientAfterSoftening = frictionCoefficientAfterSoftening;

	/* Update Drawing Parameters */
	EP_PrependClassHook( Context_GetEntryPoint( self->context, AbstractContext_EP_DumpClass ),
								_DruckerPrager_UpdateDrawParameters, self );
	
	particleExt = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, &materialPoint, self->particleExtHandle );
	
	/* Setup Variables for Visualisation */
	self->brightness = Swarm_NewScalarVariable(
		materialPointsSwarm,
		"DruckerPragerBrightness",
		(ArithPointer) &particleExt->brightness - (ArithPointer) &materialPoint,
		Variable_DataType_Float );
	
	self->opacity = Swarm_NewScalarVariable(
		materialPointsSwarm,
		"DruckerPragerOpacity",
		(ArithPointer) &particleExt->opacity - (ArithPointer) &materialPoint,
		Variable_DataType_Float );
	
	self->diameter = Swarm_NewScalarVariable(
		materialPointsSwarm,
		"DruckerPragerDiameter",
		(ArithPointer) &particleExt->diameter - (ArithPointer) &materialPoint,
		Variable_DataType_Float );

	/* The tensileFailure variable allows to check whether a materialPoint has failed in tensile mode or not */
	self->tensileFailure = Swarm_NewScalarVariable(
		materialPointsSwarm,
		"DruckerPragerTensileFailure",
		(ArithPointer) &particleExt->tensileFailure - (ArithPointer) &materialPoint,
		Variable_DataType_Char );

        self->curFrictionCoef = 0.0;

}

void* _DruckerPrager_DefaultNew( Name name ) {
	return (void*) _DruckerPrager_New(
			sizeof(DruckerPrager),
			DruckerPrager_Type,
			_YieldRheology_Delete,
			_YieldRheology_Print,
			_YieldRheology_Copy,
			_DruckerPrager_DefaultNew,
			_DruckerPrager_AssignFromXML,
			_DruckerPrager_Build,
			_DruckerPrager_Initialise,
			_YieldRheology_Execute,
			_DruckerPrager_Destroy,
			_YieldRheology_ModifyConstitutiveMatrix,
			_DruckerPrager_GetYieldCriterion,
			_VonMises_GetYieldIndicator,
			_DruckerPrager_HasYielded,
			name );
}

void _DruckerPrager_AssignFromXML( void* druckerPrager, Stg_ComponentFactory* cf, void* data ){
	DruckerPrager*          self           = (DruckerPrager*)druckerPrager;
	FeVariable*             pressureField = NULL;
	SwarmVariable*          swarmPressure = NULL;
	MaterialPointsSwarm*    materialPointsSwarm = NULL;

	/* Construct Parent */
	_VonMises_AssignFromXML( self, cf, data );
	
	pressureField      = (FeVariable *) 
            Stg_ComponentFactory_ConstructByKey( cf, self->name, "PressureField", FeVariable, False, data );
   swarmPressure = Stg_ComponentFactory_ConstructByKey(
            cf, self->name, "swarmPressure", SwarmVariable, False, data );
   Journal_Firewall( 
			( pressureField || swarmPressure), 
			Journal_Register( Error_Type, self->type ), 
			"\n Error in component type %s, name '%s'.\n Must specify a PressureField OR a swarmPressure, but not both. \n", self->type, self->name ); 
			
	materialPointsSwarm     = (MaterialPointsSwarm*)
			Stg_ComponentFactory_ConstructByKey( cf, self->name, "MaterialPointsSwarm", MaterialPointsSwarm, True, data );
		
	_DruckerPrager_Init( self, 
			pressureField,
			swarmPressure,
			materialPointsSwarm, 
			Stg_ComponentFactory_GetDouble( cf, self->name, "minimumYieldStress", 0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "frictionCoefficient", 0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "frictionCoefficientAfterSoftening", 0.0 ) );
}

void _DruckerPrager_Build( void* rheology, void* data ) {
	DruckerPrager*          self               = (DruckerPrager*) rheology;

	/* Build parent */
	_VonMises_Build( self, data );

	if(self->pressureField) Stg_Component_Build( self->pressureField, data, False );
	if(self->swarmPressure) Stg_Component_Build( self->swarmPressure, data, False );
	Stg_Component_Build( self->brightness, data, False );
	Stg_Component_Build( self->opacity, data, False );
	Stg_Component_Build( self->diameter, data, False );
	Stg_Component_Build( self->tensileFailure, data, False );

}

void _DruckerPrager_Destroy( void* rheology, void* data ) {
	DruckerPrager*          self               = (DruckerPrager*) rheology;

	if(self->pressureField) Stg_Component_Destroy( self->pressureField, data, False );
	if(self->swarmPressure) Stg_Component_Destroy( self->swarmPressure, data, False );
	Stg_Component_Destroy( self->brightness, data, False );
	Stg_Component_Destroy( self->opacity, data, False );
	Stg_Component_Destroy( self->diameter, data, False );
	Stg_Component_Destroy( self->tensileFailure, data, False );

	/* Destroy parent */
	_VonMises_Destroy( self, data );

}


void _DruckerPrager_Initialise( void* rheology, void* data ) {
	DruckerPrager*                  self                  = (DruckerPrager*) rheology;
	Particle_Index                  lParticle_I;
	Particle_Index                  particleLocalCount;
	AbstractContext*                context = (AbstractContext*)data;

	_VonMises_Initialise( self, data );

	/* Initialise variables that I've created - (mainly just SwarmVariables)
	 * This will run a Variable_Update for us */
	if(self->pressureField) Stg_Component_Initialise( self->pressureField, data, False );
	if(self->swarmPressure) Stg_Component_Initialise( self->swarmPressure, data, False );
	Stg_Component_Initialise( self->brightness, data, False );
	Stg_Component_Initialise( self->opacity, data, False );
	Stg_Component_Initialise( self->diameter, data, False );
	Stg_Component_Initialise( self->tensileFailure, data, False );

	/* We should only set initial conditions if in regular non-restart mode. If in restart mode, then
	the particle-based variables will be set correcty when we re-load the Swarm. */
	if ( !(context && (True == context->loadFromCheckPoint)) ) {
		/* We don't need to Initialise hasYieldedVariable because it's a parent variable and _YieldRheology_Initialise
		 * has already been called */
		particleLocalCount = self->hasYieldedVariable->variable->arraySize;

		for ( lParticle_I = 0 ; lParticle_I < particleLocalCount ; lParticle_I++ ) { 
		
			Variable_SetValueFloat( self->brightness->variable, lParticle_I, 0.0 );
			Variable_SetValueFloat( self->opacity->variable,    lParticle_I, 0.0 );
			Variable_SetValueFloat( self->diameter->variable,   lParticle_I, 0.0 );
		
			Variable_SetValueChar( self->hasYieldedVariable->variable, lParticle_I, False );
			Variable_SetValueChar( self->tensileFailure->variable,lParticle_I, False );
		}
	}	
}

double _DruckerPrager_GetYieldCriterion( 
			void*                            druckerPrager,
			ConstitutiveMatrix*              constitutiveMatrix,
			MaterialPointsSwarm*             materialPointsSwarm,
			Element_LocalIndex               lElement_I,
			MaterialPoint*                   materialPoint,
			Coord                            xi )
{
	DruckerPrager*                    self             = (DruckerPrager*) druckerPrager;
	double                            cohesion;
	double                            cohesionAfterSoftening;
	double                            frictionCoefficient;
	double                            frictionCoefficientAfterSoftening;
	double                            minimumYieldStress;
	double                            strainWeakeningRatio;
	double                            effectiveCohesion;
	double                            effectiveFrictionCoefficient;
	double                            phi;
	double                            sinphi;
	double                            oneOverDenominator;
	double                            dpFrictionCoefficient;
	double                            dpCohesion;
	double                            frictionalStrength;
	FeVariable*                       pressureField  = self->pressureField;  
	double                            pressure;
	DruckerPrager_Particle*           particleExt;
	
	/* Get Parameters From Rheology */
	cohesion                           = self->cohesion;
	cohesionAfterSoftening             = self->cohesionAfterSoftening;
	frictionCoefficient                = self->frictionCoefficient;
	frictionCoefficientAfterSoftening  = self->frictionCoefficientAfterSoftening;
	minimumYieldStress                 = self->minimumYieldStress;
	
	particleExt = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, materialPoint, self->particleExtHandle );

        if( self->pressureField )
          FeVariable_InterpolateWithinElement( self->pressureField, lElement_I, xi, &pressure );
        else {
          SwarmVariable_ValueAt( self->swarmPressure,
                                 constitutiveMatrix->currentParticleIndex,
                                 &pressure );
        }
	
	/* Strain softening of yield stress - if the strain weakening is not defined then this function returns 
	 * a 0 value */
	strainWeakeningRatio = StrainWeakening_CalcRatio( self->strainWeakening, materialPoint );
		
	effectiveCohesion =  cohesion * (1.0 - strainWeakeningRatio) + 
			cohesionAfterSoftening * strainWeakeningRatio;
		
	effectiveFrictionCoefficient = frictionCoefficient * (1.0 - strainWeakeningRatio) +
			frictionCoefficientAfterSoftening * strainWeakeningRatio;
	
		
	phi = atan(effectiveFrictionCoefficient);
	sinphi = sin(phi);
	oneOverDenominator = 1.0 / (sqrt(3.0) * (3.0 - sinphi));
		
	dpFrictionCoefficient = 6.0 * sinphi * oneOverDenominator;
	dpCohesion = 6.0 * effectiveCohesion * cos(phi) * oneOverDenominator;

	/* plane strain */
	
	phi = atan(effectiveFrictionCoefficient);
	dpCohesion = effectiveCohesion * cos(phi);
	dpFrictionCoefficient = sin(phi);
	
	
	frictionalStrength = dpFrictionCoefficient * pressure + dpCohesion ;

	particleExt->tensileFailure = (frictionalStrength <= 0.0);

	if ( frictionalStrength < minimumYieldStress) 
		frictionalStrength = minimumYieldStress;

        self->yieldCriterion = frictionalStrength;
        self->curFrictionCoef = dpFrictionCoefficient;

	return frictionalStrength;
}

void _DruckerPrager_HasYielded( 
		void*                            rheology,
		ConstitutiveMatrix*              constitutiveMatrix,
		MaterialPointsSwarm*             materialPointsSwarm,
		Element_LocalIndex               lElement_I,
		MaterialPoint*                   materialPoint,
		double                           yieldCriterion,
		double                           yieldIndicator )
{
   DruckerPrager* self = (DruckerPrager*)rheology;
   double viscosity = ConstitutiveMatrix_GetIsotropicViscosity( constitutiveMatrix );

   _VonMises_HasYielded(
      self, constitutiveMatrix, materialPointsSwarm, lElement_I, materialPoint,
      yieldCriterion, yieldIndicator );

   if( constitutiveMatrix->sle && constitutiveMatrix->sle->nlFormJacobian ) {

      constitutiveMatrix->derivs[8] += viscosity * self->curFrictionCoef / yieldIndicator;

   }

}

void _DruckerPrager_UpdateDrawParameters( void* rheology ) {
	DruckerPrager*                   self               = (DruckerPrager*) rheology;
	Particle_Index                   lParticle_I;
	Particle_Index                   particleLocalCount;
	StrainWeakening*                 strainWeakening    = self->strainWeakening;
	MaterialPoint*                   materialPoint;
	
	double                           length;
	double                           brightness;
	double                           opacity;
	double                           strainWeakeningRatio;
	double                           localMaxStrainIncrement;
	double                           localMeanStrainIncrement;
	Particle_Index                   localFailed;
	
	double                           globalMaxStrainIncrement;
	double                           globalMeanStrainIncrement;
	Particle_Index                   globalFailed;
	
	double                           averagedGlobalMaxStrainIncrement = 0.0;

	double                           oneOverGlobalMaxStrainIncrement;
	double                           postFailureWeakeningIncrement;

	/* Note : this function defines some drawing parameters (brightness, opacity, diameter) as
	 * functions of the strain weakening - this needs to be improved since most of the parameters
	 * that define this dependency are hard coded here. We need to have a more flexible way
	 * to construct the viz parameters as functions of material parameters */

	/* We should only update the drawing parameters if the strain weakening is defined */ 
	if (strainWeakening==NULL)
		return;
	
	localMaxStrainIncrement = 0.0;
	localMeanStrainIncrement = 0.0;
	localFailed = 0;

	/* Update all variables */
	Variable_Update( self->hasYieldedVariable->variable );
	Variable_Update( self->brightness->variable );
	Variable_Update( self->opacity->variable );
	Variable_Update( self->diameter->variable );
	Variable_Update( strainWeakening->postFailureWeakeningIncrement->variable );

	particleLocalCount = self->hasYieldedVariable->variable->arraySize;
	
	for ( lParticle_I = 0 ; lParticle_I < particleLocalCount ; lParticle_I++ ) { 
		if ( Variable_GetValueChar( self->hasYieldedVariable->variable, lParticle_I )) {
			localFailed++;

			postFailureWeakeningIncrement = 
					Variable_GetValueDouble( strainWeakening->postFailureWeakeningIncrement->variable, lParticle_I );
			
			localMeanStrainIncrement += postFailureWeakeningIncrement;
		
			if(localMaxStrainIncrement < postFailureWeakeningIncrement)
				localMaxStrainIncrement = postFailureWeakeningIncrement;
		}
	}
	
	MPI_Allreduce( &localMaxStrainIncrement,  &globalMaxStrainIncrement,  1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
	MPI_Allreduce( &localMeanStrainIncrement, &globalMeanStrainIncrement, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &localFailed,              &globalFailed,              1, MPI_INT,    MPI_SUM, MPI_COMM_WORLD );
	
	if(globalFailed == 0) 
		return;
				
	globalMeanStrainIncrement /= (double) globalFailed;
	
	averagedGlobalMaxStrainIncrement = 
			0.5 * averagedGlobalMaxStrainIncrement + 
			0.25 * globalMeanStrainIncrement +
			0.25 * globalMaxStrainIncrement;
	
	/* Let's simply assume that twice the mean is a good place to truncate these values */
	oneOverGlobalMaxStrainIncrement = 1.0 / averagedGlobalMaxStrainIncrement;
	
	for ( lParticle_I = 0 ; lParticle_I < particleLocalCount ; lParticle_I++ ) { 
		materialPoint = (MaterialPoint*)Swarm_ParticleAt( strainWeakening->swarm, lParticle_I );

		if ( Variable_GetValueChar( self->hasYieldedVariable->variable, lParticle_I ) == False ||
					StrainWeakening_GetPostFailureWeakening( strainWeakening, materialPoint ) < 0.0 ) 
		{
			Variable_SetValueFloat( self->brightness->variable, lParticle_I, 0.0 );
			Variable_SetValueFloat( self->opacity->variable, lParticle_I, 0.0 );
			Variable_SetValueFloat( self->diameter->variable, lParticle_I, 0.0 );
			continue;
		}  

		postFailureWeakeningIncrement = 
				Variable_GetValueDouble( strainWeakening->postFailureWeakeningIncrement->variable, lParticle_I );
			
		strainWeakeningRatio = StrainWeakening_CalcRatio( strainWeakening, materialPoint );
		
		length     = 0.001 + 0.003 * strainWeakeningRatio;
		brightness = strainWeakeningRatio * postFailureWeakeningIncrement * oneOverGlobalMaxStrainIncrement;
		
		opacity = 0.5 * brightness; 
		
		if( brightness > 1.0 )
			brightness = 1.0;
		
		if( opacity > 0.5 )
			opacity = 0.5;
		
		if( opacity < 0.1 )
			opacity = 0.0;
		
		Variable_SetValueFloat( self->brightness->variable, lParticle_I, brightness );
		Variable_SetValueFloat( self->opacity->variable,    lParticle_I, opacity );
		Variable_SetValueFloat( self->diameter->variable,   lParticle_I, (float) length );
	}
	
}
