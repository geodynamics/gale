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
** $Id: types.h 765 2008-07-24 05:40:02Z LukeHodkinson $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __Underworld_Utils_types_h__
#define __Underworld_Utils_types_h__
	
	typedef struct UnderworldContext		UnderworldContext;
	typedef struct PressureTemperatureOutput     PressureTemperatureOutput;
	typedef struct Underworld_SwarmOutput              Underworld_SwarmOutput;
	typedef struct RadiogenicHeatingTerm         RadiogenicHeatingTerm;
	typedef struct StressField                   StressField;
	typedef struct NodalPressureField NodalPressureField;
	typedef struct SmoothVelGradField SmoothVelGradField;
	typedef struct ViscosityField                ViscosityField;
	typedef struct DensityField                  DensityField;
   typedef struct BaseRecoveryFeVar     BaseRecoveryFeVar;
   typedef struct SPR_StrainRate SPR_StrainRate;
   typedef struct REP_Algorithm         REP_Algorithm;	
   typedef struct RecoveredFeVariable   RecoveredFeVariable;

#define REP_MAXNODESPERPATCH 30
	typedef struct {
		int nodeIDList[REP_MAXNODESPERPATCH]; /* TODO: Assume */
		int numberOfNodes;
	} LmStruct;

	/* This is a little something I need for ordering the construction of the Ax=b problem
	 * AND keeping a list of nodes which a certain patch can be applied to */
	typedef struct {
		Bool   onMeshBoundary;
		int    numOfPatches2use;
		int    patchNodes[REP_MAXNODESPERPATCH];
	} BoundaryNodesInfo;

#endif 
