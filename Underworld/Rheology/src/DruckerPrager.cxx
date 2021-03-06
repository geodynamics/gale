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
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type DruckerPrager_Type = "DruckerPrager";

/* Public Constructor */
DruckerPrager* DruckerPrager_New(
	Name                  name,
	StrainWeakening*      strainWeakening, 
	MaterialPointsSwarm*  materialPointsSwarm, 
	double                minVisc, 
	FeVariable*           strainRateField,
	double                cohesion,
	double                cohesionAfterSoftening,
	Bool                  strainRateSoftening,
	FeVariable*           pressureField,
	double                minimumYieldStress,
	double                minimumViscosity,
	double                maxStrainRate,
	double                frictionCoefficient,
	double                frictionCoefficientAfterSoftening,
        double                boundaryCohesion,
        double                boundaryCohesionAfterSoftening,
        double                boundaryFrictionCoefficient,
        double                boundaryFrictionCoefficientAfterSoftening,
        Bool                  boundaryTop,
        Bool                  boundaryBottom,
        Bool                  boundaryLeft,
        Bool                  boundaryRight,
        Bool                  boundaryFront,
        Bool                  boundaryBack,
        HydrostaticTerm*      hydrostaticTerm)
{
   DruckerPrager* self = (DruckerPrager*) _DruckerPrager_DefaultNew( name );

   _Rheology_Init( self, (PICelleratorContext*)self->context );
   _YieldRheology_Init( (YieldRheology*)self, strainWeakening, materialPointsSwarm, minVisc ); 
   _VonMises_Init( (VonMises*)self, strainRateField, cohesion,
                   cohesionAfterSoftening, strainRateSoftening );
   _DruckerPrager_Init( self, pressureField,
                        materialPointsSwarm,
                        frictionCoefficient,
                        frictionCoefficientAfterSoftening,
                        boundaryCohesion,
                        boundaryCohesionAfterSoftening,
                        boundaryFrictionCoefficient,
                        boundaryFrictionCoefficientAfterSoftening,
                        boundaryBottom,
                        boundaryTop,
                        boundaryLeft,
                        boundaryRight,
                        boundaryFront,
                        boundaryBack,
                        minimumYieldStress,
                        minimumViscosity,
                        maxStrainRate,
                        hydrostaticTerm );
   self->isConstructed = True;
   return self;
}

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
DruckerPrager* _DruckerPrager_New(  DRUCKERPRAGER_DEFARGS  ) 
{
	DruckerPrager*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(DruckerPrager) );
	self = (DruckerPrager*) _VonMises_New(  VONMISES_PASSARGS  );
	
	/* Function pointers for this class that are not on the parent class should be set here */
	
	return self;
}

void _DruckerPrager_Init(
		DruckerPrager*                                     self,
		FeVariable*                                        pressureField,
		MaterialPointsSwarm*                               materialPointsSwarm,
		double                                             frictionCoefficient,
		double                                             frictionCoefficientAfterSoftening,
                double                    boundaryCohesion,
                double                    boundaryCohesionAfterSoftening,
                double                    boundaryFrictionCoefficient,
                double                    boundaryFrictionCoefficientAfterSoftening,
                Bool                      boundaryBottom,
                Bool                      boundaryTop,
                Bool                      boundaryLeft,
                Bool                      boundaryRight,
                Bool                      boundaryFront,
                Bool                      boundaryBack,
                double                    minimumYieldStress,
                double                    minimumViscosity,
                double                    maxStrainRate,
                HydrostaticTerm*          hydrostaticTerm )
{
	/* Assign Pointers */
	self->pressureField       = pressureField;
	self->frictionCoefficient = frictionCoefficient;
	self->minimumYieldStress  = minimumYieldStress;
	self->minimumViscosity    = minimumViscosity;
	self->maxStrainRate       = maxStrainRate;
	
	/* Strain softening of Cohesion and friction - (linear
           weakening is assumed) needs a softening factor between +0
           and 1 and a reference strain > 0 */
	self->frictionCoefficientAfterSoftening = frictionCoefficientAfterSoftening;
	self->boundaryCohesion = boundaryCohesion;
	self->boundaryFrictionCoefficient = boundaryFrictionCoefficient;
	self->boundaryCohesionAfterSoftening = boundaryCohesionAfterSoftening;
	self->boundaryFrictionCoefficientAfterSoftening = boundaryFrictionCoefficientAfterSoftening;
        self->boundaryBottom=boundaryBottom;
        self->boundaryTop=boundaryTop;
        self->boundaryLeft=boundaryLeft;
        self->boundaryRight=boundaryRight;
        self->boundaryFront=boundaryFront;
        self->boundaryBack=boundaryBack;

	self->minimumYieldStress = minimumYieldStress;
	self->minimumViscosity = minimumViscosity;
	self->maxStrainRate = maxStrainRate;
        self->hydrostaticTerm=hydrostaticTerm;

	/* Update Drawing Parameters */
	EP_PrependClassHook( Context_GetEntryPoint( self->context, AbstractContext_EP_DumpClass ),
								_DruckerPrager_UpdateDrawParameters, self );
	
        self->curFrictionCoef = 0.0;
}

void* _DruckerPrager_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(DruckerPrager);
	Type                                                             type = DruckerPrager_Type;
	Stg_Class_DeleteFunction*                                     _delete = _YieldRheology_Delete;
	Stg_Class_PrintFunction*                                       _print = _YieldRheology_Print;
	Stg_Class_CopyFunction*                                         _copy = _YieldRheology_Copy;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _DruckerPrager_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _DruckerPrager_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _DruckerPrager_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _DruckerPrager_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _YieldRheology_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _DruckerPrager_Destroy;
	Rheology_ModifyConstitutiveMatrixFunction*  _modifyConstitutiveMatrix = _YieldRheology_ModifyConstitutiveMatrix;
	YieldRheology_GetYieldCriterionFunction*           _getYieldCriterion = _DruckerPrager_GetYieldCriterion;
	YieldRheology_GetYieldIndicatorFunction*           _getYieldIndicator = _VonMises_GetYieldIndicator;
	YieldRheology_HasYieldedFunction*                         _hasYielded = _DruckerPrager_HasYielded;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _DruckerPrager_New(  DRUCKERPRAGER_PASSARGS   );
}

void _DruckerPrager_AssignFromXML( void* druckerPrager, Stg_ComponentFactory* cf, void* data ){
	DruckerPrager*          self           = (DruckerPrager*)druckerPrager;
	FeVariable*             pressureField = NULL;
	MaterialPointsSwarm*    materialPointsSwarm = NULL;

	/* Construct Parent */
	_VonMises_AssignFromXML( self, cf, data );
	
	pressureField      = (FeVariable *) 
            Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"PressureField", FeVariable, True, data  );
			
	materialPointsSwarm     = (MaterialPointsSwarm*)
			Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"MaterialPointsSwarm", MaterialPointsSwarm, True, data  );
		
	_DruckerPrager_Init( self, 
			pressureField,
			materialPointsSwarm, 
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"frictionCoefficient", 0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"frictionCoefficientAfterSoftening", 0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"boundaryCohesion", 0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"boundaryCohesionAfterSoftening", 0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"boundaryFrictionCoefficient", 0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"boundaryFrictionCoefficientAfterSoftening", 0.0 ),
                        Stg_ComponentFactory_GetBool(  cf,  self->name, (Dictionary_Entry_Key)"boundaryBottom", False ),
                        Stg_ComponentFactory_GetBool(  cf,  self->name, (Dictionary_Entry_Key)"boundaryTop", False ),
                        Stg_ComponentFactory_GetBool(  cf,  self->name, (Dictionary_Entry_Key)"boundaryLeft", False ),
                        Stg_ComponentFactory_GetBool(  cf,  self->name, (Dictionary_Entry_Key)"boundaryRight", False ),
                        Stg_ComponentFactory_GetBool(  cf,  self->name, (Dictionary_Entry_Key)"boundaryFront", False ),
                        Stg_ComponentFactory_GetBool(  cf,  self->name, (Dictionary_Entry_Key)"boundaryBack", False ),
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"minimumYieldStress", 0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"minimumViscosity", 0.0),
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"maxStrainRate", 0.0),
                        Stg_ComponentFactory_ConstructByKey( cf, self->name, 
                                                             (Dictionary_Entry_Key)"HydrostaticTerm", HydrostaticTerm, False, data ) );
}

void _DruckerPrager_Build( void* rheology, void* data ) {
	DruckerPrager*          self               = (DruckerPrager*) rheology;

	/* Build parent */
	_VonMises_Build( self, data );

        Stg_Component_Build( self->pressureField, data, False );
}

void _DruckerPrager_Destroy( void* rheology, void* data ) {
	DruckerPrager*          self               = (DruckerPrager*) rheology;

        Stg_Component_Destroy( self->pressureField, data, False );

	/* Destroy parent */
	_VonMises_Destroy( self, data );
}


void _DruckerPrager_Initialise( void* rheology, void* data ) {
	DruckerPrager*                  self                  = (DruckerPrager*) rheology;
	Particle_Index                  lParticle_I;
	Particle_Index                  particleLocalCount;

	_VonMises_Initialise( self, data );

	/* Initialise variables that I've created - (mainly just SwarmVariables)
	 * This will run a Variable_Update for us */
        Stg_Component_Initialise( self->pressureField, data, False );

	/* We should only set initial conditions if in regular non-restart mode. If in restart mode, then
	the particle-based variables will be set correcty when we re-load the Swarm. */
	if ( self->context->loadFromCheckPoint == False ) {

		/* We don't need to Initialise hasYieldedVariable because it's a parent variable and _YieldRheology_Initialise
		 * has already been called */
		particleLocalCount = self->hasYieldedVariable->variable->arraySize;

		for ( lParticle_I = 0 ; lParticle_I < particleLocalCount ; lParticle_I++ ) { 
		
			Variable_SetValueChar( self->hasYieldedVariable->variable, lParticle_I, False );
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
        Dimension_Index                   dim = constitutiveMatrix->dim;
	double                            minimumYieldStress;
	double                            effectiveCohesion;
	double                            effectiveFrictionCoefficient;
	double                            frictionalStrength;
	double                            pressure;
        Cell_Index                        cell_I;
        Coord                             coord;
        Element_GlobalIndex	          element_gI = 0;
        unsigned		          inds[3];
        Grid*			          elGrid;
        Bool                              inside_boundary;
        double                            factor;
	
	/* Get Parameters From Rheology */
	minimumYieldStress                 = self->minimumYieldStress;
	
        FeVariable_InterpolateWithinElement( self->pressureField, lElement_I, xi, &pressure );

        cell_I=CellLayout_MapElementIdToCellId(materialPointsSwarm->cellLayout,
                                               lElement_I );
        FeMesh_CoordLocalToGlobal(self->pressureField->feMesh, cell_I, xi, coord);
        if(self->hydrostaticTerm)
          {
            pressure+=HydrostaticTerm_Pressure(self->hydrostaticTerm,coord);
          }

        /* Normally we add the average of the trace of the stress.
           With compressible material, you have to do it.  But with
           stabilized linear pressure elements, the non-zero trace is
           a numerical artifact.  So we do not add it. */

/*   pressure+=self->trace/dim; */

        /* Calculate frictional strength.  We modify the friction and
           cohesion because we have grouped terms from the normal
           stresses and moved it to the yield indicator. */
	

        /* Big song and dance to see if we are at a boundary that we care about */
        elGrid = *(Grid**)ExtensionManager_Get(self->pressureField->feMesh->info,
                                               self->pressureField->feMesh,
                                               ExtensionManager_GetHandle(self->pressureField->feMesh->info, "elementGrid" ) );
  
        element_gI = FeMesh_ElementDomainToGlobal( self->pressureField->feMesh, lElement_I );
        RegularMeshUtils_Element_1DTo3D( self->pressureField->feMesh, element_gI, inds );
  
        inside_boundary=
          ((self->boundaryBottom && inds[1]==0)
           || (self->boundaryTop && inds[1]==elGrid->sizes[1]-1)
           || (self->boundaryLeft && inds[0]==0)
           || (self->boundaryRight && inds[0]==elGrid->sizes[0]-1)
           || (dim==3 && ((self->boundaryBack && inds[2]==0)
                          || (self->boundaryFront && inds[2]==elGrid->sizes[2]-1))))
          ? True : False;

        effectiveFrictionCoefficient =
          _DruckerPrager_EffectiveFrictionCoefficient( self, materialPoint,
                                                       inside_boundary );
        effectiveCohesion = _DruckerPrager_EffectiveCohesion(self,materialPoint,
                                                             inside_boundary);

  if(dim==2)
    {
      /* effectiveFrictionCoefficient=tan(phi).  If
         factor=sin(atan(1/tan(phi))) =>
         factor=cos(phi)=1/sqrt(1+tan(phi)**2) */
      factor=1/sqrt(1 + effectiveFrictionCoefficient*effectiveFrictionCoefficient);
      frictionalStrength = effectiveFrictionCoefficient*pressure*factor
        + effectiveCohesion*factor;
    }
  else
    {
      double cos_phi, sin_phi;
      /* cos(phi)=1/sqrt(1+tan(phi)**2) */
      cos_phi=
        1/sqrt(1 + effectiveFrictionCoefficient*effectiveFrictionCoefficient);
      sin_phi=effectiveFrictionCoefficient*cos_phi;
      factor=2*cos_phi/(sqrt(3.0)*(3-sin_phi));

      /* The full expression is

         sqrt(J2)=p*2*sin(phi)/(sqrt(3)*(3-sin(phi)))
                  + C*6*cos(phi)/(sqrt(3)*(3-sin(phi)))

         Note the extra factor of 3 for cohesion */

      frictionalStrength = effectiveFrictionCoefficient*factor*pressure
        + effectiveCohesion*3*factor;
    }
  
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

  self->yieldCriterion = frictionalStrength;
  self->curFrictionCoef = effectiveFrictionCoefficient*factor;

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
   double strainRate=self->strainRateSecondInvariant;
   double old_viscosity = ConstitutiveMatrix_GetIsotropicViscosity( constitutiveMatrix );
   double viscosity;

   if(self->maxStrainRate>0 && strainRate>self->maxStrainRate)
     {
       viscosity=yieldCriterion/(2*self->maxStrainRate);
     }
   else
     {
       viscosity = yieldCriterion/(2*strainRate);
     }

   if(viscosity<self->minimumViscosity)
     viscosity=self->minimumViscosity;

   ConstitutiveMatrix_SetIsotropicViscosity( constitutiveMatrix, viscosity );

   if( constitutiveMatrix->sle && constitutiveMatrix->sle->nlFormJacobian ) {
      constitutiveMatrix->derivs[8] += old_viscosity * self->curFrictionCoef / yieldIndicator;
   }

}

double _DruckerPrager_EffectiveCohesion( void* rheology, void* materialPoint,
                                       Bool inside_boundary ) {
	DruckerPrager*    self                          = (DruckerPrager*) rheology;
	double                         effectiveCohesion;
	
        double strainWeakeningRatio = StrainWeakening_CalcRatio( self->strainWeakening, materialPoint );

        if(inside_boundary)
          {
            effectiveCohesion=self->boundaryCohesion*(1.0-strainWeakeningRatio)+
              self->boundaryCohesionAfterSoftening*strainWeakeningRatio;
          }
        else
          {
            effectiveCohesion = self->cohesion * (1.0 - strainWeakeningRatio) + 
              self->cohesionAfterSoftening * strainWeakeningRatio;
          }

	return effectiveCohesion;
}

double _DruckerPrager_EffectiveFrictionCoefficient( void* rheology, void* materialPoint, Bool inside_boundary ) {
	DruckerPrager*    self                          = (DruckerPrager*) rheology;
	double                  effectiveFrictionCoefficient  = 0.0;

        double strainWeakeningRatio = StrainWeakening_CalcRatio( self->strainWeakening, materialPoint );
	
        if(inside_boundary)
          {
            effectiveFrictionCoefficient =
              self->boundaryFrictionCoefficient * (1.0 - strainWeakeningRatio) +
              self->boundaryFrictionCoefficientAfterSoftening
              *strainWeakeningRatio;	
          }
        else
          {
            effectiveFrictionCoefficient =
              self->frictionCoefficient * (1.0 - strainWeakeningRatio) +
              self->frictionCoefficientAfterSoftening * strainWeakeningRatio;	
          }

	return effectiveFrictionCoefficient;
}

void _DruckerPrager_UpdateDrawParameters( void* rheology ) {
  DruckerPrager* self=(DruckerPrager*) rheology;
  StrainWeakening* strainWeakening=self->strainWeakening;
  /* We should only update the drawing parameters if the strain
     weakening is defined */ 
  if (strainWeakening==NULL)
    return;

  /* Update all variables */
  Variable_Update( self->hasYieldedVariable->variable );
  Variable_Update(strainWeakening->postFailureWeakeningIncrement->variable);
}


