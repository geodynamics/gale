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
** $Id: NonNewtonian.c 610 2007-10-11 08:09:29Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "RheologyClass.h"
#include "NonNewtonian.h"
#include "ConstitutiveMatrix.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type NonNewtonian_Type = "NonNewtonian";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
NonNewtonian* _NonNewtonian_New( 
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
		Name                                               name ) 
{
	NonNewtonian*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(NonNewtonian) );
	self = (NonNewtonian*) _Rheology_New( 
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
			name );
	
	return self;
}

void _NonNewtonian_Init( NonNewtonian* self, FeVariable* strainRateInvField,
                         FeVariable* temperatureField,
                         double n,
                         double T_0,
                         double A,
                         double refStrainRate,
                         double minViscosity,
                         double maxViscosity) {

	self->strainRateInvField = strainRateInvField;
	self->temperatureField = temperatureField;

	self->n = n;
	self->T_0 = T_0;
	self->A = A;
	self->refStrainRate = refStrainRate;
	self->minViscosity = minViscosity;
	self->maxViscosity = maxViscosity;

	Rheology_SetToNonLinear( self );
}

void* _NonNewtonian_DefaultNew( Name name ) {
	return (void*) _NonNewtonian_New(
		sizeof(NonNewtonian),
		NonNewtonian_Type,
		_Rheology_Delete,
		_Rheology_Print,
		_Rheology_Copy,
		_NonNewtonian_DefaultNew,
		_NonNewtonian_Construct,
		_Rheology_Build,
		_Rheology_Initialise,
		_Rheology_Execute,
		_Rheology_Destroy,
		_NonNewtonian_ModifyConstitutiveMatrix,
		name );
}

void _NonNewtonian_Construct( void* rheology, Stg_ComponentFactory* cf, void* data ){
	NonNewtonian*  self = (NonNewtonian*)rheology;
	FeVariable    *strainRateInvField, *temperatureField;
        double refStrainRate;
	Stream* errorStream = Journal_Register( Error_Type,
                                                NonNewtonian_Type );

	/* Construct Parent */
	_Rheology_Construct( self, cf, data );
	
	/* TODO: 'Keyfallback' soon to be deprecated/updated */
	strainRateInvField = Stg_ComponentFactory_ConstructByNameWithKeyFallback( 
		cf, 
		self->name,
                "StrainRateInvariantField", 
		"StrainRateInvariantField", 
		FeVariable, 
		True,
		data );
	temperatureField = Stg_ComponentFactory_ConstructByKey(
		cf, 
		self->name,
                "TemperatureField", 
		FeVariable, 
		True,
		data );

        refStrainRate=Stg_ComponentFactory_GetDouble( cf, self->name,
                                                      "refStrainRate", 0.0 );
        if(refStrainRate==0.0)
          Journal_Firewall(0,errorStream,
                           "refStrainRate must be set to a non-zero number in a NonNewtonian rheology");

	_NonNewtonian_Init( 
			self,
			strainRateInvField,
                        temperatureField,
			Stg_ComponentFactory_GetDouble( cf, self->name, "n", 1.0 ), 
			Stg_ComponentFactory_GetDouble( cf, self->name, "T_0", 1.0 ), 
			Stg_ComponentFactory_GetDouble( cf, self->name, "A", 1.0 ), 
                        refStrainRate,
			Stg_ComponentFactory_GetDouble( cf, self->name, "minViscosity", 0.0 ) ,
			Stg_ComponentFactory_GetDouble( cf, self->name, "maxViscosity", 0.0 ) );
}

void _NonNewtonian_ModifyConstitutiveMatrix( 
		void*                                              rheology, 
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               swarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi )
{
	NonNewtonian*	              self              = (NonNewtonian*) rheology;
	double                        strainRateInv;
	double                        minVis;
	double                        maxVis;
	double                        viscosity;
	double                        n, T_0, A, T;
	Stream* errorStream = Journal_Register( Error_Type,
                                                NonNewtonian_Type );

	n = self->n;
        T_0=self->T_0;
        A=self->A;
	minVis = self->minViscosity;
	maxVis = self->maxViscosity;

        FeVariable_InterpolateWithinElement(self->temperatureField,
                                            lElement_I,xi,&T);
        Journal_Firewall(n!=0 && T!=0 && A!=0,errorStream,
                         "Error in NonNewtonian: T, n, and A must all be non-zero:\n\tT=%g\n\tn=%g\n\tA=%g\n",T,n,A);

        viscosity=exp(T_0/(n*T))/(2*A);

	if ( fabs( n - 1.0 ) > 1.0e-10 )
          {
            if ( !constitutiveMatrix->previousSolutionExists )
              /* on the first ever solve, use a ref strainrate */
              {
                strainRateInv = self->refStrainRate;
              }
            else
              {
                /* Calculate Parameters */
                FeVariable_InterpolateWithinElement(self->strainRateInvField,
                                                    lElement_I,xi,&strainRateInv);
              }
            /* Calculate New Viscosity */
            viscosity*= pow(strainRateInv,1.0/n - 1.0);
          }
        if(maxVis!=0 && viscosity>maxVis)
          {
            viscosity=maxVis;
          }
        if(minVis!=0 && viscosity<minVis)
          {
            viscosity=minVis;
          }
	ConstitutiveMatrix_SetIsotropicViscosity(constitutiveMatrix,viscosity);
}
