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
** $Id:  $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include <string.h>

/* Each component type needs a unique identifier (as a string for us to read and as an integer for realtime comparisions) */
const Type Underworld_Ra_Scaling_Type = "Underworld_Ra_Scaling_Type";

void _Ra_CheckScalings_Func( void* context, void* ptrToContext ) {
  UnderworldContext*      self         = (UnderworldContext*) context;
  ForceTerm* bfTerm         = (ForceTerm*)LiveComponentRegister_Get( self->CF->LCRegister, "buoyancyForceTerm" );

  /* check the Rayleigh Number Scaling:
   * first check the RHS for the stokes SLE has a force term which uses Ra*/
  if( bfTerm != NULL && (bfTerm->type == ThermalBuoyancyForceTerm_Type ) ) {
    RheologyMaterial* material;
    Rheology* rheology;
    Materials_Register*     materials_Register = self->materials_Register;
    char* errorMesg = "";
    double Ra, Ra_0, eta0, diffusivity, gravity, thermalExp=1;
    int isValid = 1; /* is this scaling check valid */

    /* check if there's only one material */
    if( Materials_Register_GetCount( materials_Register ) > 1 ) {
      isValid=0; 
			Journal_RPrintf( self->info, "\n\n** Warning - The Ra_ScalingCheck plugin is on but you're using more than two materials, thus the Ra_ScalingCheck is invalid because a Rayleigh number doesn't make sense when multiply materials are used **\n" );
    }

		/* get the material rheology */
    material = (RheologyMaterial*)Materials_Register_GetByIndex( materials_Register, 0 );

    /* check if there's only one rheology on that material */
		if (material->rheology_Register->objects->count > 1 ) isValid = 0;
    
    Ra        = ((ThermalBuoyancyForceTerm*)bfTerm)->rayleighNumber;
    rheology  = Rheology_Register_GetByIndex( material->rheology_Register, 0 ); /* assume one rheology */

    /* check if this is a valid test */
    if( rheology->type == MaterialViscosity_Type && isValid ) { 
      eta0 = ((MaterialViscosity*)rheology)->eta0;

      diffusivity   = Stg_ComponentFactory_GetDouble( self->CF, "defaultResidualForceTerm", "defaultDiffusivity", 1.0 );
      gravity       = Stg_ComponentFactory_GetRootDictDouble( self->CF, "gravity", 1.0 );

      Ra_0 = (gravity * thermalExp)/(diffusivity * eta0 );
			/* check if the Ra matches the viscosity, gravity, thermal expansivity and diffusivity of the problem */
      if( abs(Ra - Ra_0) > 1e-3 ) {
        Stg_asprintf( &errorMesg, "* Error - Your combination of diffusivity, gravity and eta0 (rheology) don't agree with your Ra:"
            "input:\ndiffusivity = %g\ngravity = %g\neta0 = %g from rheology %s\nthermal expansivity force = %g\nRa = %g\n", diffusivity, gravity, eta0, rheology->name, thermalExp, Ra );
        Journal_RPrintf( self->info, "\n\n** Scaling issues detected **\n%s\nIf you believe the problem(s) above are ok in your model and you want to get rid of this error message use the command line argument --Ra_ScalingCheck=false\n\n", errorMesg );
        exit(0);
      }
    }
  }
}

void _Underworld_Ra_Scaling_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
   UnderworldContext* context = Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data ); 

   Bool checkScaling = Stg_ComponentFactory_GetRootDictBool( cf, "Ra_ScalingCheck", True ); 

   if ( checkScaling ) {
      EntryPoint_Append( Context_GetEntryPoint( context, AbstractContext_EP_Build ),
      "Underworld CheckScalings",
      _Ra_CheckScalings_Func,
      Underworld_Ra_Scaling_Type );
   }
}

/* This function will provide StGermain the abilty to instantiate (create) this codelet on demand. */
void* _Underworld_Ra_Scaling_DefaultNew( Name name ) {
	return Codelet_New(
			Underworld_Ra_Scaling_Type,
			_Underworld_Ra_Scaling_DefaultNew,
			_Underworld_Ra_Scaling_AssignFromXML, /* SQ NOTE: Used to be a construct extensions. */
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}

/* This function is automatically run by StGermain when this plugin is loaded. The name must be "<plugin-name>_Register". */
Index Underworld_Ra_Scaling_Register( PluginsManager* pluginsManager ) {
	/* A plugin is only properly registered once it returns the handle provided when submitting a codelet to StGermain. */
	return PluginsManager_Submit( pluginsManager, Underworld_Ra_Scaling_Type, "0", _Underworld_Ra_Scaling_DefaultNew );
}
