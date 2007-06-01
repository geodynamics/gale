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
** $Id: FaultingMoresiMuhlhaus2006.c 483 2007-06-01 00:56:27Z StuartClark $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <stdlib.h>
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "RheologyClass.h"
#include "StrainWeakening.h"
#include "YieldRheology.h"
#include "FaultingMoresiMuhlhaus2006.h"
#include "ConstitutiveMatrix.h"
#include "Director.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type FaultingMoresiMuhlhaus2006_Type = "FaultingMoresiMuhlhaus2006";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
FaultingMoresiMuhlhaus2006* _FaultingMoresiMuhlhaus2006_New( 
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
	FaultingMoresiMuhlhaus2006*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(FaultingMoresiMuhlhaus2006) );
	self = (FaultingMoresiMuhlhaus2006*) _YieldRheology_New( 
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

void _FaultingMoresiMuhlhaus2006_Init(
		FaultingMoresiMuhlhaus2006*                        self,
		FeVariable*                                        pressureField,
		FeVariable*                                        velocityGradientsField,
		MaterialPointsSwarm*                               materialPointsSwarm,
		FiniteElementContext*                              context,
		Director*                                          director,
		double                                             cohesion,
		double                                             cohesionAfterSoftening,
		double                                             frictionCoefficient,
		double                                             frictionCoefficientAfterSoftening,
		double                                             minimumYieldStress,
		Bool                                               ignoreOldOrientation )
{
	FaultingMoresiMuhlhaus2006_Particle* particleExt;
	StandardParticle                    materialPoint;
	Dimension_Index                     dim   = materialPointsSwarm->dim;
	
	self->materialPointsSwarm     = materialPointsSwarm;
	self->pressureField           = pressureField;
	self->velocityGradientsField  = velocityGradientsField;
	
	self->director                = director;
	
	self->particleExtHandle       = ExtensionManager_Add( materialPointsSwarm->particleExtensionMgr,
			FaultingMoresiMuhlhaus2006_Type, sizeof(FaultingMoresiMuhlhaus2006_Particle) );	
	
	self->cohesion = cohesion;
	self->frictionCoefficient = frictionCoefficient;
	
	/* Strain softening of Cohesion - (linear weakening is assumed) */
	/* needs a softening factor between +0 and 1 and a reference strain > 0 */
	self->cohesionAfterSoftening = cohesionAfterSoftening;
	
	/* Strain softening of Friction - (linear weakening is assumed) */
	/* needs a softening factor between +0 and 1 and a reference strain > 0 */
	self->frictionCoefficientAfterSoftening = frictionCoefficientAfterSoftening;

	self->minimumYieldStress = minimumYieldStress;
	
	/* Should orientation from previous timestep be tested ? */
	self->ignoreOldOrientation = ignoreOldOrientation;

	/* Update Drawing Parameters */
	EP_PrependClassHook( Context_GetEntryPoint( context, AbstractContext_EP_DumpClass ),
			_FaultingMoresiMuhlhaus2006_UpdateDrawParameters, self );

	particleExt = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, &materialPoint, self->particleExtHandle );

	/* Add variables for viz purpose */

	self->slipRate = Swarm_NewScalarVariable(
			materialPointsSwarm,
			"SlipRate",
			(ArithPointer) &particleExt->slipRate - (ArithPointer) &materialPoint,
			Variable_DataType_Double );

	/* slip vector gives the orientation of the failure plane */ 
	self->slip = Swarm_NewVectorVariable(
			materialPointsSwarm,
			"SlipVector",
			(ArithPointer) &particleExt->slip - (ArithPointer) &materialPoint, 
			Variable_DataType_Double, 
			dim, 
			"SlipVectorX",
			"SlipVectorY",
			"SlipVectorZ" );	
	
	/* Some visualisation parameters (brightness, opacity, length, thickness) */
	self->brightness = Swarm_NewScalarVariable(
			materialPointsSwarm,
			"FaultingMoresiMuhlhaus2006Brightness",
			(ArithPointer) &particleExt->brightness - (ArithPointer) &materialPoint,
			Variable_DataType_Float );
	
	self->opacity = Swarm_NewScalarVariable(
			materialPointsSwarm,
			"FaultingMoresiMuhlhaus2006Opacity",
			(ArithPointer) &particleExt->opacity - (ArithPointer) &materialPoint,
			Variable_DataType_Float );
	
	self->length = Swarm_NewScalarVariable(
			materialPointsSwarm,
			"FaultingMoresiMuhlhaus2006Length",
			(ArithPointer) &particleExt->length - (ArithPointer) &materialPoint,
			Variable_DataType_Float );

	self->thickness = Swarm_NewScalarVariable(
			materialPointsSwarm,
			"FaultingMoresiMuhlhaus2006Thickness",
			(ArithPointer) &particleExt->thickness - (ArithPointer) &materialPoint,
			Variable_DataType_Float );
	
	/* The tensileFailure variable allows to check whether a materialPoint has failed in tensile mode or not */
	self->tensileFailure = Swarm_NewScalarVariable(
			materialPointsSwarm,
			"FaultingMoresiMuhlhaus2006TensileFailure",
			(ArithPointer) &particleExt->tensileFailure - (ArithPointer) &materialPoint,
			Variable_DataType_Char );
}

void* _FaultingMoresiMuhlhaus2006_DefaultNew( Name name ) {
	return (void*) _FaultingMoresiMuhlhaus2006_New(
		sizeof(FaultingMoresiMuhlhaus2006),
		FaultingMoresiMuhlhaus2006_Type,
		_YieldRheology_Delete,
		_YieldRheology_Print,
		_YieldRheology_Copy,
		_FaultingMoresiMuhlhaus2006_DefaultNew,
		_FaultingMoresiMuhlhaus2006_Construct,
		_FaultingMoresiMuhlhaus2006_Build,
		_FaultingMoresiMuhlhaus2006_Initialise,
		_YieldRheology_Execute,
		_YieldRheology_Destroy,
		_FaultingMoresiMuhlhaus2006_ModifyConstitutiveMatrix,
		_FaultingMoresiMuhlhaus2006_GetYieldCriterion,
		_FaultingMoresiMuhlhaus2006_GetYieldIndicator,
		_FaultingMoresiMuhlhaus2006_HasYielded,
		name );
}

void _FaultingMoresiMuhlhaus2006_Construct( void* rheology, Stg_ComponentFactory* cf, void* data ){
	FaultingMoresiMuhlhaus2006*   self           = (FaultingMoresiMuhlhaus2006*)rheology;
	FeVariable*                   pressureField;
	MaterialPointsSwarm*          materialPointsSwarm;
	FeVariable*                   velocityGradientsField;
	FiniteElementContext*         context;
	Director*                     director;
	
	/* Construct Parent */
	_YieldRheology_Construct( self, cf, data );

	/* Make sure that there is strain weakening */
	Journal_Firewall(
		self->strainWeakening != NULL,
		Journal_Register( Error_Type, self->type ),
		"Error in func '%s' for %s '%s': FaultingMoresiMuhlhaus2006 rheology needs strain weakening.\n", 
		__func__, self->type, self->name );
	
	context                = Stg_ComponentFactory_ConstructByName( cf, 
			"context", FiniteElementContext, True, data );
	materialPointsSwarm    = Stg_ComponentFactory_ConstructByKey( cf, self->name, 
			"MaterialPointsSwarm", MaterialPointsSwarm, True, data );
	pressureField          = Stg_ComponentFactory_ConstructByKey( cf, self->name,
			"PressureField", FeVariable, True, data );
	velocityGradientsField = Stg_ComponentFactory_ConstructByKey( cf, self->name,
			"VelocityGradientsField", FeVariable, True, data );
	director               =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "Director", Director, True, data );
	
	_FaultingMoresiMuhlhaus2006_Init( 
			self,
			pressureField,
			velocityGradientsField,
			materialPointsSwarm, 
			context,
			director,
			Stg_ComponentFactory_GetDouble( cf, self->name, "cohesion", 0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "cohesionAfterSoftening", 0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "frictionCoefficient", 0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "frictionCoefficientAfterSoftening", 0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "minimumYieldStress", 0.0 ),
			Stg_ComponentFactory_GetBool(   cf, self->name, "ignoreOldOrientation", False )	);
}

void _FaultingMoresiMuhlhaus2006_Build( void* rheology, void* data ) {
	FaultingMoresiMuhlhaus2006*  self               = (FaultingMoresiMuhlhaus2006*) rheology;

	/* Build parent */
	_YieldRheology_Build( self, data );

	Stg_Component_Build( self->slipRate, data, False );
	Stg_Component_Build( self->slip, data, False );
	Stg_Component_Build( self->brightness, data, False );
	Stg_Component_Build( self->opacity, data, False );
	Stg_Component_Build( self->length, data, False );
	Stg_Component_Build( self->thickness, data, False );
	Stg_Component_Build( self->tensileFailure, data, False );

}

void _FaultingMoresiMuhlhaus2006_Initialise( void* rheology, void* data ) {
	FaultingMoresiMuhlhaus2006*     self                  = (FaultingMoresiMuhlhaus2006*) rheology;
	XYZ                             normal                = { 0.0,1.0,0.0 };
	XYZ                             slip                  = { 1.0,0.0,0.0 };
	double*                         ptr;
	Particle_Index                  lParticle_I;
	Particle_Index                  particleLocalCount;
	double                          normalLength2;
	double                          invNormalLength;
	double                          initialDamageFraction;
	double*                         normalDirector;
	MaterialPoint*                  materialPoint;
	Index                           dof_I;
	Dimension_Index                 dim                   = self->materialPointsSwarm->dim;
	AbstractContext*                context = (AbstractContext*)data;

	_YieldRheology_Initialise( self, data );
	
	/* Initialise these components now just in case this function is called before their own _Initialise function */
	Stg_Component_Initialise( self->materialPointsSwarm, data, False );
	Stg_Component_Initialise( self->director, data, False );
	Stg_Component_Initialise( self->strainWeakening, data, False );

	/* Initialise variables that I've created - (mainly just SwarmVariables)
	 * This will run a Variable_Update for us */
	Stg_Component_Initialise( self->slipRate, data, False );
	Stg_Component_Initialise( self->slip, data, False );
	Stg_Component_Initialise( self->brightness, data, False );
	Stg_Component_Initialise( self->opacity, data, False );
	Stg_Component_Initialise( self->length, data, False );
	Stg_Component_Initialise( self->thickness, data, False );
	Stg_Component_Initialise( self->tensileFailure, data, False );

	/* We don't need to Initialise hasYieldedVariable because it's a parent variable and _YieldRheology_Initialise
	 * has already been called */
	particleLocalCount = self->hasYieldedVariable->variable->arraySize;

	/* If restarting from checkpoint, don't change the parameters on the particles */
	if ( !(context && (True == context->loadFromCheckPoint) ) ) {
		for ( lParticle_I = 0 ; lParticle_I < particleLocalCount ; lParticle_I++ ) { 
			Variable_SetValueChar( self->hasYieldedVariable->variable, lParticle_I, False );
			Variable_SetValueDouble( self->slipRate->variable, lParticle_I, 0.0 );
			
			ptr = Variable_GetPtrDouble( self->slip->variable, lParticle_I );
			
			materialPoint = (MaterialPoint*)Swarm_ParticleAt( self->materialPointsSwarm, lParticle_I );
			normalDirector = Director_GetNormalPtr( self->director, materialPoint);
			initialDamageFraction = StrainWeakening_GetInitialDamageFraction( self->strainWeakening, materialPoint );
			
			/*  Is this out of sync with the other damage model ??
				Surely the orientation should be randomized for all particles where the initial damage is non-zero ??
				LM  */
		
		
			if (1 || rand() < RAND_MAX*initialDamageFraction) {
				normalLength2 = 0.0;
				
				for( dof_I=0; dof_I < dim ; dof_I++) {
					normal[dof_I] = 1.0 - (2.0 * rand())/RAND_MAX; 
					normalLength2 += normal[dof_I] * normal[dof_I];
				}
				
				invNormalLength = 1.0/sqrt(normalLength2);
				
				for( dof_I=0; dof_I < dim ; dof_I++){
					normal[dof_I] *= invNormalLength; 
				}
				
				/*TODO : improve this initialisation (is it really needed ?)
					Dear Dr TODO, If you mean "is the slip really needed" then
					I think the answer is "no" .... LM 
				*/
				slip[0] = - normal[1];
				slip[1] =   normal[0];
				slip[2] =   0.0;

			}
			
			memcpy( normalDirector, normal, sizeof(Coord) );
			memcpy( ptr, slip, sizeof(Coord) );
			
			Variable_SetValueDouble( self->slipRate->variable,     lParticle_I, 0.0 );

			Variable_SetValueFloat( self->opacity->variable,      lParticle_I, 0.0 );
			Variable_SetValueFloat( self->length->variable,       lParticle_I, 0.0 );
			Variable_SetValueFloat( self->thickness->variable,    lParticle_I, 0.0 );
			Variable_SetValueChar( self->tensileFailure->variable,lParticle_I, False );
		}
	}
}
	
void _FaultingMoresiMuhlhaus2006_ModifyConstitutiveMatrix( 
		void*                                              rheology, 
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               materialPointsSwarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi )
{
	FaultingMoresiMuhlhaus2006*     self                  = (FaultingMoresiMuhlhaus2006*) rheology;
	double                          postFailureWeakening;
	double                          yieldCriterion;
	double                          yieldIndicator;  /* A materialPoint will yield if yieldCriterion < yieldIndicator */

	/* Don't want to yield on the first ever solve */
	if ( constitutiveMatrix->previousSolutionExists == False ) {
		return;
	}

	/* Calculate and store values of stress, pressure, velocity gradients, eigenvectors so they are only calculated once */
	_FaultingMoresiMuhlhaus2006_StoreCurrentParameters( self, constitutiveMatrix, materialPointsSwarm, lElement_I, materialPoint, xi );

	postFailureWeakening = StrainWeakening_GetPostFailureWeakening( self->strainWeakening, materialPoint );
	
	/* First part : treatment of the existing weakened directions */
	if ( !self->ignoreOldOrientation && postFailureWeakening > 0.0 ) {
	
		/* tryingOldOrientation is a flag used to know where we are. This is necessary because some parts are treated differently
		 * if we are dealing with old orientation or with pristine material (see for example _FaultingMoresiMuhlhaus2006_GetYieldIndicator) */
		self->tryingOldOrientation = True;

		if ( _FaultingMoresiMuhlhaus2006_OldOrientationStillSoftening( self, materialPointsSwarm, materialPoint, constitutiveMatrix->dim ) ) {
			yieldCriterion = _FaultingMoresiMuhlhaus2006_GetYieldCriterion( self, constitutiveMatrix, materialPointsSwarm,
				lElement_I, materialPoint, xi );
			yieldIndicator = _FaultingMoresiMuhlhaus2006_GetYieldIndicator( self, constitutiveMatrix, materialPointsSwarm,
				lElement_I, materialPoint, xi );

			/* Set a bool to TRUE or FLAG depending on whether a materialPoint has failed or not */
			YieldRheology_SetParticleFlag( self, materialPointsSwarm, materialPoint, (yieldCriterion < yieldIndicator) );
			
			if ( yieldCriterion < yieldIndicator ) {
				StrainWeakening_AssignIncrement( self->strainWeakening, constitutiveMatrix,
						materialPointsSwarm, lElement_I, materialPoint, yieldCriterion, yieldIndicator );
				_FaultingMoresiMuhlhaus2006_HasYielded( self, constitutiveMatrix, materialPointsSwarm, lElement_I, materialPoint,
						yieldCriterion, yieldIndicator );
				return;
			}
		}
	}

	/* Second part : treatment of pristine material with no direction favoured */
	self->tryingOldOrientation = False;
	_YieldRheology_ModifyConstitutiveMatrix( self, constitutiveMatrix, materialPointsSwarm, lElement_I, materialPoint, xi );
}

double _FaultingMoresiMuhlhaus2006_GetYieldIndicator( 
		void*                                              rheology,
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               materialPointsSwarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi )
{
	FaultingMoresiMuhlhaus2006*           self             = (FaultingMoresiMuhlhaus2006*) rheology;
	double*                               stress           = self->currentStress;
	Eigenvector*                          eigenvectorList  = self->currentEigenvectorList;
	double                                sigma_ns;
	Dimension_Index                       dim              = constitutiveMatrix->dim;
	FaultingMoresiMuhlhaus2006_Particle*  particleExt;
	double                                stressMin;
	double                                stressMax;
	double                                theta;
	XYZ                                   normal;	
	
	Director_GetNormal( self->director, materialPoint, normal );
	
	particleExt = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, materialPoint, self->particleExtHandle );

	if ( self->tryingOldOrientation ) {
		sigma_ns = SymmetricTensor_MultiplyByVectors( stress, normal, particleExt->slip, dim );
	}
	else {
		/* Failure criterion from stress field - 
		 * we can calculate whether the material has failed at the current stress from the principle stresses 
		 * (without having to compute the orientation) if it does fail, 
		 * we can then calculate the orientation for the orthotropy model and history */
		theta = 0.5 * atan( 1.0/ _FaultingMoresiMuhlhaus2006_EffectiveFrictionCoefficient( self, materialPoint ) );
		
		stressMin = eigenvectorList[0].eigenvalue;
		stressMax = (dim == 2 ? eigenvectorList[1].eigenvalue : eigenvectorList[2].eigenvalue);
		
		sigma_ns = 0.5 *  sin( 2.0 * theta ) * ( stressMax - stressMin );
	}
	
	return fabs(sigma_ns);
}

double _FaultingMoresiMuhlhaus2006_GetYieldCriterion( 
		void*                                              rheology,
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               materialPointsSwarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi )
{
	FaultingMoresiMuhlhaus2006*          self             = (FaultingMoresiMuhlhaus2006*) rheology;
	double                               minimumYieldStress;
	double                               effectiveCohesion;
	double                               effectiveFrictionCoefficient;
	double                               frictionalStrength;
	double                               sigma_nn;
	FaultingMoresiMuhlhaus2006_Particle* particleExt;
		
	particleExt = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, materialPoint, self->particleExtHandle );
		
	/* Calculate frictional strength */
	effectiveFrictionCoefficient = _FaultingMoresiMuhlhaus2006_EffectiveFrictionCoefficient( self, materialPoint );
	effectiveCohesion            = _FaultingMoresiMuhlhaus2006_EffectiveCohesion( self, materialPoint );
	sigma_nn                     = _FaultingMoresiMuhlhaus2006_Sigma_nn( self, materialPoint, constitutiveMatrix->dim );

	frictionalStrength = effectiveFrictionCoefficient * sigma_nn + effectiveCohesion;

	/* Check if it should break in tension (strictly : frictionalStrength < tensileStrength) */
	particleExt->tensileFailure = (frictionalStrength < 0.0);
	
	/* Make sure frictionalStrength is above the minimum */
	minimumYieldStress = self->minimumYieldStress;
	
	if ( frictionalStrength < minimumYieldStress) 
		frictionalStrength = minimumYieldStress;
	
	return frictionalStrength;
}

void _FaultingMoresiMuhlhaus2006_HasYielded( 
		void*                                              rheology,
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               materialPointsSwarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		double                                             yieldCriterion,
		double                                             yieldIndicator )
{
	FaultingMoresiMuhlhaus2006*          self             = (FaultingMoresiMuhlhaus2006*) rheology;
	double                               beta;
	double*                              normal;
	double                               viscosity        = ConstitutiveMatrix_GetIsotropicViscosity( constitutiveMatrix );
	FaultingMoresiMuhlhaus2006_Particle* particleExt;

	normal = _FaultingMoresiMuhlhaus2006_UpdateNormalDirection( self, materialPointsSwarm, materialPoint, constitutiveMatrix->dim );
	
	beta = 1.0 - yieldCriterion / yieldIndicator;
	
	ConstitutiveMatrix_SetSecondViscosity( constitutiveMatrix, viscosity * beta, normal ); 

	particleExt = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, materialPoint, self->particleExtHandle );
	particleExt->slipRate = self->storedSlipRateValue;
	
	/* Set a flag to tell the director that particles which have failed don't need to have their normal updated */
	Director_SetDontUpdateParticleFlag( self->director, materialPoint, True );
}

Bool _FaultingMoresiMuhlhaus2006_OldOrientationStillSoftening( void* rheology, MaterialPointsSwarm* materialPointsSwarm, void* materialPoint, Dimension_Index dim ) {
	FaultingMoresiMuhlhaus2006*           self             = (FaultingMoresiMuhlhaus2006*) rheology;
	double*                               velocityGradient = self->currentVelocityGradient;
	double*                               stress           = self->currentStress;
	FaultingMoresiMuhlhaus2006_Particle*  particleExt;
	double*                               n;
	double*                               s;
	double                                tau_nn;
	double                                traction[3];
	double                                dVparalleldXperpendicular;
	double                                dVparalleldXperpendicular1;
	XYZ                                   normal;
	
	Director_GetNormal( self->director, materialPoint, normal );

	particleExt = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, materialPoint, self->particleExtHandle );
	n = normal;
	s = particleExt->slip;
	 
	/* If an existing weakened direction exists we test this direction for further failure
	 * The slip direction should be set to the direction of maximum shear stress and stored
	 * for plotting 
	 * (it does not need to be updated as a history variable like the director ) */
	if ( dim == 2 ) {
		traction[0] = stress[0] * n[0] + stress[2] * n[1];
		traction[1] = stress[2] * n[0] + stress[1] * n[1];
		
		tau_nn =  n[0] * traction[0] + n[1] * traction[1];
				
		s[0] = traction[0] - tau_nn * n[0];
		s[1] = traction[1] - tau_nn * n[1];
	}
	else {
		traction[0] = stress[0] * n[0] + stress[3] * n[1] + stress[4] * n[2];
		traction[1] = stress[3] * n[0] + stress[1] * n[1] + stress[5] * n[2];
		traction[2] = stress[4] * n[0] + stress[5] * n[1] + stress[2] * n[2];	
				
		tau_nn =  n[0] * traction[0] + n[1] * traction[1] + n[2] * traction[2];
				
		s[0] = traction[0] - tau_nn * n[0];
		s[1] = traction[1] - tau_nn * n[1];
		s[2] = traction[2] - tau_nn * n[2];
	}
	StGermain_VectorNormalise( s, dim );
	
	/* Store tau_nn so that we can use it in calculating sigma_nn in _FaultingMoresiMuhlhaus2006_Sigma_nn */
	self->tau_nn = tau_nn;

	/* Softening direction */
	dVparalleldXperpendicular  = fabs(TensorArray_MultiplyByVectors( velocityGradient, s, n, dim )) ;

	/* Hardening Direction */
	dVparalleldXperpendicular1 = fabs(TensorArray_MultiplyByVectors( velocityGradient, n, s, dim )) ;  

	/* Store this value on class in case we need to store it on the materialPoint if it fails */
	self->storedSlipRateValue = dVparalleldXperpendicular;

	/* We should only continue testing this plane if it is also in the 
	 * softening orientation with respect to the strain-rate gradient rather than hardening */
	
	return (dVparalleldXperpendicular1 < dVparalleldXperpendicular);
}

double* _FaultingMoresiMuhlhaus2006_UpdateNormalDirection( void* rheology, MaterialPointsSwarm* materialPointsSwarm, void* materialPoint, Dimension_Index dim ) {
	FaultingMoresiMuhlhaus2006*          self             = (FaultingMoresiMuhlhaus2006*) rheology;
	FaultingMoresiMuhlhaus2006_Particle* particleExt;
	double*                              velocityGradient = self->currentVelocityGradient;
	Eigenvector*                         eigenvectorList  = self->currentEigenvectorList;
	double                               strainRate_ns[2];
	XYZ                                  normal[2];
	XYZ                                  slip[2];
	double                               theta;
	int                                  favourablePlane;
	double*                              normalDirector;
	double                               tanPhi;
	
	normalDirector = Director_GetNormalPtr( self->director, materialPoint);
	
	particleExt = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, materialPoint, self->particleExtHandle );
	
	/* This function is specific for pristine materials -- 
	 * We don't need to calculate it for failure in old orientations */
	if ( self->tryingOldOrientation )
		return normalDirector;
	
	if ((fabs(tanPhi = _FaultingMoresiMuhlhaus2006_EffectiveFrictionCoefficient( self, materialPoint )))<0.000001)
		theta = M_PI / 4.0;
	else
		theta = 0.5 * atan( 1.0/ tanPhi );
	
	if (dim == 2){ 
		
		/* Identify potential failure directions */
		StGermain_RotateCoordinateAxis( slip[0],   eigenvectorList[0].vector, K_AXIS, +theta);
		StGermain_RotateCoordinateAxis( slip[1],   eigenvectorList[0].vector, K_AXIS, -theta);				
		StGermain_RotateCoordinateAxis( normal[0], eigenvectorList[0].vector, K_AXIS,  0.5*M_PI + theta);
		StGermain_RotateCoordinateAxis( normal[1], eigenvectorList[0].vector, K_AXIS,  0.5*M_PI - theta);
	}
	else {
			
		/* Identify potential failure directions */
		StGermain_RotateVector( slip[0],   eigenvectorList[0].vector, eigenvectorList[1].vector, + theta );
		StGermain_RotateVector( slip[1],   eigenvectorList[0].vector, eigenvectorList[1].vector, - theta );
		StGermain_RotateVector( normal[0], eigenvectorList[0].vector, eigenvectorList[1].vector,  0.5*M_PI + theta);
		StGermain_RotateVector( normal[1], eigenvectorList[0].vector, eigenvectorList[1].vector,  0.5*M_PI - theta);
	}
				
	/* Resolve shear strain-rate for the potential failure planes */
	strainRate_ns[0] = fabs(TensorArray_MultiplyByVectors( velocityGradient, slip[0], normal[0], dim ));
	strainRate_ns[1] = fabs(TensorArray_MultiplyByVectors( velocityGradient, slip[1], normal[1], dim ));
				
	/* Choose the plane which is oriented favorably for continued slip */
	favourablePlane = strainRate_ns[0] > strainRate_ns[1] ? 0 : 1;

	memcpy( normalDirector,      normal[favourablePlane], dim * sizeof(double) );
	memcpy( particleExt->slip,   slip[favourablePlane],   dim * sizeof(double) );
	
	self->storedSlipRateValue = strainRate_ns[ favourablePlane ];

	return normalDirector;
}

double _FaultingMoresiMuhlhaus2006_EffectiveCohesion( void* rheology, void* materialPoint ) {
	FaultingMoresiMuhlhaus2006*    self                          = (FaultingMoresiMuhlhaus2006*) rheology;
	double                         effectiveCohesion;
	
	if ( self->tryingOldOrientation || self->ignoreOldOrientation ) {
		/* If this is the old orientation or old orientations are ignored
		 * then the weakening should be considered here   */
		double strainWeakeningRatio = StrainWeakening_CalcRatio( self->strainWeakening, materialPoint );
	
		effectiveCohesion =  self->cohesion * (1.0 - strainWeakeningRatio) + 
				self->cohesionAfterSoftening * strainWeakeningRatio;
	}
	else {
		/*	If old orientations are given priority and this is an unbroken direction 
		 *	then consider material as completely pristine */
		effectiveCohesion =  self->cohesion;
	}

	return effectiveCohesion;
}

double _FaultingMoresiMuhlhaus2006_EffectiveFrictionCoefficient( void* rheology, void* materialPoint ) {
	FaultingMoresiMuhlhaus2006*    self                          = (FaultingMoresiMuhlhaus2006*) rheology;
	static double                  effectiveFrictionCoefficient  = 0.0;
	static int                     prevTryingOldOrientationValue = -1;
	static void*                   prevParticle                  = NULL;

	/* See if we can return cached value */
	if ( self->tryingOldOrientation == (Bool) prevTryingOldOrientationValue && materialPoint == prevParticle ) 
		return effectiveFrictionCoefficient;
	
	/* This is a different situation to the cached case - therefore recalculate */
	if ( self->tryingOldOrientation || self->ignoreOldOrientation ) {
		/* If this is the old orientation or old orientations are ignored
		 * then the weakening should be considered here   */
		double strainWeakeningRatio = StrainWeakening_CalcRatio( self->strainWeakening, materialPoint );
	
		effectiveFrictionCoefficient = self->frictionCoefficient * (1.0 - strainWeakeningRatio) +
				self->frictionCoefficientAfterSoftening * strainWeakeningRatio;	
	}
	else {
		/*	If old orientations are given priority and this is an unbroken direction 
		 *	then consider material as completely pristine */
		effectiveFrictionCoefficient =  self->frictionCoefficient;		
	}

	/* Cache values */
	prevParticle                  = materialPoint;
	prevTryingOldOrientationValue = self->tryingOldOrientation;

	return effectiveFrictionCoefficient;
}

double _FaultingMoresiMuhlhaus2006_Sigma_nn( void* rheology, void* materialPoint, Dimension_Index dim ) {
	FaultingMoresiMuhlhaus2006*      self             = (FaultingMoresiMuhlhaus2006*) rheology;
	double                           pressure         = self->currentPressure;
	Eigenvector*                     eigenvectorList  = self->currentEigenvectorList;
	double                           tau_nn;
	double                           sigma_nn;

	if ( self->tryingOldOrientation ) {
		/* tau_nn has already been calculated in _FaultingMoresiMuhlhaus2006_OldOrientationStillSoftening */
		tau_nn = self->tau_nn;
	}
	else {
		double stressMin = eigenvectorList[0].eigenvalue;
		double stressMax = (dim == 2 ? eigenvectorList[1].eigenvalue : eigenvectorList[2].eigenvalue);
		double theta = 0.5 * atan( 1.0/ _FaultingMoresiMuhlhaus2006_EffectiveFrictionCoefficient( self, materialPoint ) );

		tau_nn = 0.5 * (( stressMax + stressMin ) + cos( 2.0 * theta ) * ( stressMax - stressMin ));
	}

	sigma_nn = pressure - tau_nn;

	return sigma_nn;
}

void _FaultingMoresiMuhlhaus2006_StoreCurrentParameters( 
		void*                                              rheology,
		ConstitutiveMatrix*                                constitutiveMatrix, 
		MaterialPointsSwarm*                               materialPointsSwarm,
		Element_LocalIndex                                 lElement_I, 
		MaterialPoint*                                     materialPoint,
		Coord                                              xi ) 
{
	FaultingMoresiMuhlhaus2006*          self               = (FaultingMoresiMuhlhaus2006*) rheology;
	SymmetricTensor                      strainRate;
	Dimension_Index                      dim                = constitutiveMatrix->dim;
	FaultingMoresiMuhlhaus2006_Particle* particleExt;
	
	particleExt = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, materialPoint, self->particleExtHandle );
	particleExt->slipRate = 0.0;
	
	FeVariable_InterpolateWithinElement( self->pressureField, lElement_I, xi, &self->currentPressure );	
	FeVariable_InterpolateWithinElement( self->velocityGradientsField, lElement_I, xi, self->currentVelocityGradient );
	
	TensorArray_GetSymmetricPart( self->currentVelocityGradient, dim, strainRate );
	ConstitutiveMatrix_CalculateStress( constitutiveMatrix, strainRate, self->currentStress );
	
	SymmetricTensor_CalcAllEigenvectors( self->currentStress, dim, self->currentEigenvectorList );

	Director_SetDontUpdateParticleFlag( self->director, materialPoint, False );
}

void _FaultingMoresiMuhlhaus2006_UpdateDrawParameters( void* rheology ) {
	FaultingMoresiMuhlhaus2006*       self               = (FaultingMoresiMuhlhaus2006*) rheology;
	Particle_Index                   lParticle_I;
	Particle_Index                   particleLocalCount;
	StrainWeakening*                 strainWeakening    = self->strainWeakening;
	StandardParticle*                materialPoint;
	double                           slipRate;
	
	double 							 ratio;
	double                           length;
	double                           brightness;
	double                           opacity;
	double                           strainWeakeningRatio;
	double                           localMaxStrainIncrement;
	double                           localMinStrainIncrement;
	double                           localMaxSlipRate;
	double                           localMeanStrainIncrement;
	double                           localMeanSlipRate;
	Particle_Index                   localFailed;
	
	double                           globalMaxSlipRate;
	double                           globalMaxStrainIncrement;
	double                           globalMinStrainIncrement;
	double                           globalMeanStrainIncrement;
	double                           globalMeanSlipRate;
	Particle_Index                   globalFailed;
	
	double                           averagedGlobalMaxSlipRate = 0.0;
	double                           averagedGlobalMaxStrainIncrement = 0.0;

	double                           oneOverGlobalMaxSlipRate;
	double                           oneOverGlobalMaxStrainIncrement;
	double                           postFailureWeakeningIncrement;
	
	/* This parameter is only needed for the alternative commented out set of parameters */ 
	/* double                           globalMaxStrainWeakeningRatio    = 0.0; */

        /* Note : this function defines some drawing parameters (brightness, opacity, diameter) as
	 * functions of the strain weakening - this needs to be improved since most of the parameters
	 * that define this dependency are hard coded here. We need to have a more flexible way
	 * to construct the viz parameters as functions of material parameters */
	
	localMaxStrainIncrement = -1.0e30;
	localMinStrainIncrement =  1.0e30;
	localMeanStrainIncrement = 0.0;
	localMaxSlipRate = 0.0;
	localMeanSlipRate = 0.0;
	localFailed = 0;

	/* Update all variables */
	Variable_Update( self->hasYieldedVariable->variable );
	Variable_Update( self->slipRate->variable );
	Variable_Update( self->brightness->variable );
	Variable_Update( self->opacity->variable );
	Variable_Update( self->length->variable );
	Variable_Update( self->thickness->variable );
	Variable_Update( strainWeakening->postFailureWeakeningIncrement->variable );

	particleLocalCount = self->hasYieldedVariable->variable->arraySize;
	
	for ( lParticle_I = 0 ; lParticle_I < particleLocalCount ; lParticle_I++ ) { 
		if ( Variable_GetValueChar( self->hasYieldedVariable->variable, lParticle_I )) {
			localFailed++;
			slipRate = Variable_GetValueDouble( self->slipRate->variable, lParticle_I );

			postFailureWeakeningIncrement = 
				Variable_GetValueDouble( strainWeakening->postFailureWeakeningIncrement->variable, lParticle_I );
			
			localMeanSlipRate += slipRate; /* Always positive */
			localMeanStrainIncrement += fabs(postFailureWeakeningIncrement); /* Positive in failure, negative in healing */
			
			if(localMaxSlipRate < slipRate)
				localMaxSlipRate = slipRate;
		
			if(localMaxStrainIncrement < postFailureWeakeningIncrement)
				localMaxStrainIncrement = postFailureWeakeningIncrement;
				
			if(localMinStrainIncrement > postFailureWeakeningIncrement)
				localMinStrainIncrement = postFailureWeakeningIncrement;
		}
	}
	
	MPI_Allreduce( &localMaxSlipRate,         &globalMaxSlipRate,         1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
	MPI_Allreduce( &localMaxStrainIncrement,  &globalMaxStrainIncrement,  1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
	MPI_Allreduce( &localMinStrainIncrement,  &globalMinStrainIncrement,  1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
	MPI_Allreduce( &localMeanSlipRate,        &globalMeanSlipRate,        1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &localMeanStrainIncrement, &globalMeanStrainIncrement, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &localFailed,              &globalFailed,              1, MPI_INT,    MPI_SUM, MPI_COMM_WORLD );
	
	if(localFailed == 0 || globalFailed == 0) 
		return;
				
	globalMeanStrainIncrement /= (double) globalFailed;
	globalMeanSlipRate /= (double) globalFailed;
	
	averagedGlobalMaxStrainIncrement = 
		0.5 * averagedGlobalMaxStrainIncrement + 
		0.25 * globalMeanStrainIncrement +
		0.25 * globalMaxStrainIncrement;
	averagedGlobalMaxSlipRate = 
		0.5 * averagedGlobalMaxSlipRate + 
		0.25 * globalMaxSlipRate +
		0.25 * globalMeanSlipRate;
	
	fprintf(stderr,"globalMeanStrainIncrement = %g\n",globalMeanStrainIncrement);
	fprintf(stderr,"globalMaxStrainIncrement  = %g\n",globalMaxStrainIncrement);
	fprintf(stderr,"globalMinStrainIncrement  = %g\n",globalMinStrainIncrement);
	fprintf(stderr,"averagedGlobalMaxStrainIncrement = %g\n",averagedGlobalMaxStrainIncrement);
	
	fprintf(stderr,"globalMeanSlipRate = %g\n",globalMeanSlipRate); 
	fprintf(stderr,"averagedGlobalMaxSlipRate = %g\n",averagedGlobalMaxSlipRate);
	
	fprintf(stderr,"Number of particles failing = %d\n",localFailed);
	
	
	
	/* Let's simply assume that twice the mean is a good place to truncate these values */
	oneOverGlobalMaxSlipRate = 1.0 / averagedGlobalMaxSlipRate;
	oneOverGlobalMaxStrainIncrement = 1.0 / averagedGlobalMaxStrainIncrement;
	
	for ( lParticle_I = 0 ; lParticle_I < particleLocalCount ; lParticle_I++ ) { 
		materialPoint = Swarm_ParticleAt( strainWeakening->swarm, lParticle_I );

		if ( Variable_GetValueChar( self->hasYieldedVariable->variable, lParticle_I ) == False ||
				StrainWeakening_GetPostFailureWeakening( strainWeakening, materialPoint ) < 0.0 ) 
		{
			Variable_SetValueFloat( self->brightness->variable, lParticle_I, 0.0 );
			Variable_SetValueFloat( self->opacity->variable, lParticle_I, 0.0 );
			Variable_SetValueFloat( self->length->variable, lParticle_I, 0.0 );
			Variable_SetValueFloat( self->thickness->variable, lParticle_I, 0.0 );
			continue;
		}  
		
		slipRate = Variable_GetValueDouble( self->slipRate->variable, lParticle_I );
		strainWeakeningRatio = StrainWeakening_CalcRatio( strainWeakening, materialPoint );
		
		postFailureWeakeningIncrement = 
			Variable_GetValueDouble( strainWeakening->postFailureWeakeningIncrement->variable, lParticle_I );
		
		
		/* Distinguish between healing (e.g. dark) and softening (e.g. light) branches 
		
			1. The intensity of the colouring should be proportional to strain rate
			2. Size of square is proportional to total strain
			3. The opacity is used to mask out the "uninteresting cases"
		
		*/
		
				length = 0.025 * strainWeakeningRatio;
				
				ratio = postFailureWeakeningIncrement / globalMeanStrainIncrement; // healing v. not 
				if( ratio < -0.0) {  /* Healing - don't bother with the ones which are fading fast */
					brightness = 0.0;	
					opacity = 1.0 + ratio;
					if(opacity < 0.0)
						opacity = 0.0;
				}
				else {  /* Growing  */
					brightness = slipRate * oneOverGlobalMaxSlipRate;
					opacity = 0.25 * (ratio < 1.0 ? ratio : 1.0);
					
				}
				
					
					
		/*
		length     = 0.005 + 0.005 * strainWeakeningRatio;
		length *= 10;
		brightness = strainWeakeningRatio * slipRate * oneOverGlobalMaxSlipRate;
		if( brightness > 1.0 )
			brightness = 1.0;
		opacity    = 0.1 + 0.05 * pow(brightness,3.0); 
		*/

		Variable_SetValueFloat( self->brightness->variable, lParticle_I, brightness );
		Variable_SetValueFloat( self->opacity->variable,    lParticle_I, opacity );
		Variable_SetValueFloat( self->length->variable,     lParticle_I, (float) length );
		Variable_SetValueFloat( self->thickness->variable,  lParticle_I,  0.5 + 10.0 * (float)length );
	
	/* This parameter is only needed for the alternative commented out set of parameters */ 
	/*	if (strainWeakeningRatio>globalMaxStrainWeakeningRatio)
			globalMaxStrainWeakeningRatio = strainWeakeningRatio;
	*/
	}
	
	/* This part is another set of parameters that has been used with shear problems.
	 * it's kept here for now, but it shows that in principle we could define as many set of parameters
	 * as we want (depending on the problem ?)*/
	/*
	for ( lParticle_I = 0 ; lParticle_I < particleLocalCount ; lParticle_I++ ) { 
		materialPoint = Swarm_ParticleAt( strainWeakening->materialPointsSwarm, lParticle_I );

		if ( Variable_GetValueChar( self->hasYieldedVariable->variable, lParticle_I ) == False ||
				StrainWeakening_GetPostFailureWeakening( strainWeakening, materialPoint ) < 0.0 ) 
		{
			Variable_SetValueFloat( self->brightness->variable, lParticle_I, 0.0 );
			Variable_SetValueFloat( self->opacity->variable, lParticle_I, 0.0 );
			Variable_SetValueFloat( self->length->variable, lParticle_I, 0.0 );
			Variable_SetValueFloat( self->thickness->variable, lParticle_I, 0.0 );
			continue;
		}  
		
		slipRate = Variable_GetValueDouble( self->slipRate->variable, lParticle_I );
		strainWeakeningRatio = StrainWeakening_CalcRatio( strainWeakening, materialPoint );
		
		length     = 0.03 + 0.03 * (strainWeakeningRatio/globalMaxStrainWeakeningRatio);

		brightness = 1.0-(strainWeakeningRatio/globalMaxStrainWeakeningRatio);
		
		if (brightness > 1.0)
			brightness = 1.0;
			
		opacity    = (slipRate/globalMaxSlipRate);
		
		if (opacity > 0.90)
			opacity = 1.0;/* this condition is to make sure we have enough planes that will be clearly seen. */ /*

		Variable_SetValueFloat( self->brightness->variable, lParticle_I, brightness );
		Variable_SetValueFloat( self->opacity->variable,    lParticle_I, opacity );
		Variable_SetValueFloat( self->length->variable,     lParticle_I, (float) length );
		Variable_SetValueFloat( self->thickness->variable,  lParticle_I,  0.5 + 10.0 * (float)length );
	}
	*/
}
