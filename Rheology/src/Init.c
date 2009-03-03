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
** $Id: Init.c 750 2008-07-07 02:26:33Z LukeHodkinson $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "Rheology.h"

Bool Underworld_Rheology_Init( int* argc, char** argv[] ) {
	Stg_ComponentRegister* componentRegister = Stg_ComponentRegister_Get_ComponentRegister();

	Journal_Printf( Journal_Register( DebugStream_Type, "Context" ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */
	
	Stg_ComponentRegister_Add( componentRegister, ConstitutiveMatrixCartesian_Type, "0", _ConstitutiveMatrixCartesian_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, MaterialViscosity_Type,       "0", _MaterialViscosity_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, RheologyMaterial_Type,        "0", _RheologyMaterial_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, MultiRheologyMaterial_Type,   "0", _MultiRheologyMaterial_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, Compressible_Type,            "0", _Compressible_DefaultNew );

	Stg_ComponentRegister_Add( componentRegister, Arrhenius_Type,               "0", _Arrhenius_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, FrankKamenetskii_Type,        "0", _FrankKamenetskii_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, NonNewtonian_Type,            "0", _NonNewtonian_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, LinearViscosity_Type,         "0", _LinearViscosity_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, DepthDependentViscosity_Type, "0", _DepthDependentViscosity_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, Anisotropic_Type,             "0", _Anisotropic_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, OrthotropicAligned_Type,      "0", _OrthotropicAligned_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, Orthotropic_Type,             "0", _Orthotropic_DefaultNew );
	
	Stg_ComponentRegister_Add( componentRegister, VonMises_Type,                "0", _VonMises_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, Byerlee_Type,                 "0", _Byerlee_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, DruckerPrager_Type,           "0", _DruckerPrager_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, FaultingMoresiMuhlhaus2006_Type,             "0", _FaultingMoresiMuhlhaus2006_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, MohrCoulomb_Type,             "0", _MohrCoulomb_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, Pouliquen_etal_Type,          "0", _Pouliquen_etal_DefaultNew );

	Stg_ComponentRegister_Add( componentRegister, StrainWeakening_Type,         "0", _StrainWeakening_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, BuiterStrainWeakening_Type,   "0", _BuiterStrainWeakening_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, Director_Type,                "0", _Director_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, AlignmentSwarmVariable_Type,  "0", _AlignmentSwarmVariable_DefaultNew );
	
	Stg_ComponentRegister_Add( componentRegister, ViscosityFieldRheology_Type,  "0", _ViscosityFieldRheology_DefaultNew );

	Stg_ComponentRegister_Add( componentRegister, StoreStress_Type,             "0", _StoreStress_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, StoreVisc_Type,               "0", _StoreVisc_DefaultNew );
	
	Stg_ComponentRegister_Add( componentRegister, ConstitutiveMatCartesian_Refactored_Type, "0", _ConstitutiveMatCartesian_Refactored_DefaultNew );

	/* Register Parents for type checking */
	RegisterParent( Rheology_Type,                Stg_Component_Type );
	RegisterParent( Arrhenius_Type,               Rheology_Type );
	RegisterParent( FrankKamenetskii_Type,        Rheology_Type );
	RegisterParent( MaterialViscosity_Type,       Rheology_Type );
	RegisterParent( NonNewtonian_Type,            Rheology_Type );
	RegisterParent( DepthDependentViscosity_Type, Rheology_Type );
	RegisterParent( LinearViscosity_Type,         Rheology_Type );
	RegisterParent( Anisotropic_Type,             Rheology_Type );
	RegisterParent( OrthotropicAligned_Type,      Rheology_Type );
	RegisterParent( Orthotropic_Type,             Rheology_Type );

	RegisterParent( StoreStress_Type,             Rheology_Type );
	RegisterParent( StoreVisc_Type,               Rheology_Type );
	
	RegisterParent( ViscosityFieldRheology_Type,     Rheology_Type );
	RegisterParent( YieldRheology_Type,              Rheology_Type );
	RegisterParent( VonMises_Type,                   YieldRheology_Type );
	RegisterParent( Byerlee_Type,                    VonMises_Type );
	RegisterParent( DruckerPrager_Type,              VonMises_Type );
	RegisterParent( FaultingMoresiMuhlhaus2006_Type, YieldRheology_Type );
	RegisterParent( MohrCoulomb_Type,                YieldRheology_Type );
	RegisterParent( Pouliquen_etal_Type,             VonMises_Type );
	
	RegisterParent( StrainWeakening_Type,         TimeIntegratee_Type );
	RegisterParent( BuiterStrainWeakening_Type,   StrainWeakening_Type );
	RegisterParent( Director_Type,                TimeIntegratee_Type );
	RegisterParent( AlignmentSwarmVariable_Type,  SwarmVariable_Type );

	RegisterParent( ConstitutiveMatrix_Type,          StiffnessMatrixTerm_Type );
	RegisterParent( ConstitutiveMatrixCartesian_Type, ConstitutiveMatrix_Type );
	RegisterParent( RheologyMaterial_Type,            Material_Type );
	RegisterParent( MultiRheologyMaterial_Type,       RheologyMaterial_Type );
	RegisterParent( Compressible_Type,                StiffnessMatrixTerm_Type );

	RegisterParent( ConstitutiveMat_Refactored_Type,          Stg_Component_Type );
	RegisterParent( ConstitutiveMatCartesian_Refactored_Type, ConstitutiveMat_Refactored_Type );

	return True;
}
