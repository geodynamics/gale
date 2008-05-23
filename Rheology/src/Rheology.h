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
** $Id: Rheology.h 736 2008-05-23 06:05:21Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __Underworld_Rheology_h__
#define __Underworld_Rheology_h__

	#include "types.h"
	#include "RheologyClass.h"
	#include "YieldRheology.h"
	#include "Rheology_Register.h"
	#include "RheologyMaterial.h"
	#include "MultiRheologyMaterial.h"
	#include "Compressible.h"
	#include "ConstitutiveMatrix.h"
	#include "ConstitutiveMatrixCartesian.h"
	#include "MaterialViscosity.h"
	#include "Arrhenius.h"
	#include "FrankKamenetskii.h"
	#include "NonNewtonian.h"
	#include "DepthDependentViscosity.h"
	#include "LinearViscosity.h"
	#include "Anisotropic.h"
        #include "OrthotropicAligned.h"
        #include "Orthotropic.h"
	#include "VonMises.h"
	#include "ViscosityFieldRheology.h"
	#include "Byerlee.h"
	#include "DruckerPrager.h"
	#include "FaultingMoresiMuhlhaus2006.h"
	#include "StrainWeakening.h"
	#include "BuiterStrainWeakening.h"
	#include "Director.h"
	#include "AlignmentSwarmVariable.h"
	#include "StoreStress.h"
	#include "StoreViscosity.h"
	#include "Init.h"
	#include "Finalise.h"
	
#endif
