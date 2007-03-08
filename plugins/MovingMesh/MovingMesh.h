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
** $Id: MovingMesh.h 354 2006-10-12 08:19:27Z SteveQuenette $
** 
** Comments:
**   Not really sure why this is a TimeIntegrator of all node coords based on the velocity field,
**	since as Vincent & I discussed it really should be based only on the BCs. Would like to
**	refactor at some stage. Perhaps its due to the coupling with particle advection.
**	-- Main. PatrickSunter
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __Underworld_MovingMesh_h__
#define __Underworld_MovingMesh_h__

	typedef enum MinOrMaxFlag { MIN_COORDS, MAX_COORDS } MinOrMaxFlag;

	/** The MeshExtender extends upon the Codelet class */
	#define __MeshExtender \
		__Codelet \
		\
		FeVariable*  velocityField; \
		AbstractContext*  context; /** Needed to get dt and calculate updated width of box */\
		Bool         remeshAccordingToAxis[3];

	typedef struct MeshExtender { __MeshExtender } MeshExtender;

	Index Underworld_MovingMesh_Register( PluginsManager* pluginsManager );
	void* _Underworld_MovingMesh_DefaultNew( Name name );
void Underworld_MovingMesh_Construct( void* component, Stg_ComponentFactory* _cf, void* data ) ;

	void Underworld_MovingMesh_Remesh( TimeIntegrator* timeIntegrator, MeshExtender* self ) ;
	void Underworld_MovingMesh_RemeshAccordingToSidewalls( MeshExtender* self, FeVariable* velocityField ) ;
	void Underworld_MovingMesh_RemeshAccordingToSidewall_SingleAxis(
		MeshExtender*    self,
		FeVariable*      velocityField,
		Dimension_Index  remeshAxis ) ;
	
	void Underworld_MovingMesh_CalculateMinOrMaxCoordsOnSidewall(
		MeshExtender*    self,
		FeVariable*      velocityField,
		Dimension_Index  remeshAxis,
		MinOrMaxFlag     minOrMaxFlag,
		double*          newWallCoordsInRemeshAxisGlobal );

#endif
