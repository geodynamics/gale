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
*%		Cecile Duboz - Cecile.Duboz@sci.monash.edu.au
*%
** Contributors:
*+		Cecile Duboz
*+		Robert Turnbull
*+		Alan Lo
*+		Louis Moresi
*+		David Stegman
*+		David May
*+		Stevan Quenette
*+		Patrick Sunter
*+		Greg Watson
*+
** $Id: Init.c 740 2007-10-11 08:05:31Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <glucifer/Base/Base.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>

#include "DrawingObjects.h"

Bool lucDrawingObjects_Init() {
	Stg_ComponentRegister* componentRegister = Stg_ComponentRegister_Get_ComponentRegister();

	Journal_Printf( Journal_Register( DebugStream_Type, (Name)"Context"  ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */
	

	Stg_ComponentRegister_Add( componentRegister, lucColourBar_Type, (Name)"0", _lucColourBar_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, lucFieldVariableBorder_Type, (Name)"0", _lucFieldVariableBorder_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, lucIsosurface_Type, (Name)"0", _lucIsosurface_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, lucCrossSection_Type, (Name)"0", _lucCrossSection_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, lucScalarFieldCrossSection_Type, (Name)"0", _lucScalarFieldCrossSection_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, lucScalarField_Type, (Name)"0", _lucScalarField_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, lucVectorArrowCrossSection_Type, (Name)"0", _lucVectorArrowCrossSection_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, lucVectorArrows_Type, (Name)"0", _lucVectorArrows_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, lucTextureMap_Type, (Name)"0", _lucTextureMap_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, lucContour_Type, (Name)"0", _lucContour_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, lucFeVariableSurface_Type, (Name)"0", _lucFeVariableSurface_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, lucSwarmViewer_Type, (Name)"0", _lucSwarmViewer_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, lucSwarmVectors_Type, (Name)"0", _lucSwarmVectors_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, lucSwarmSquares_Type, (Name)"0", _lucSwarmSquares_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, lucHistoricalSwarmTrajectory_Type, (Name)"0", _lucHistoricalSwarmTrajectory_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, lucEigenvectorsCrossSection_Type, (Name)"0", _lucEigenvectorsCrossSection_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, lucEigenvectors_Type, (Name)"0", _lucEigenvectors_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, lucSwarmRGBColourViewer_Type, (Name)"0", _lucSwarmRGBColourViewer_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, lucMeshViewer_Type, (Name)"0", _lucMeshViewer_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, lucTitle_Type, (Name)"0", _lucTitle_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, lucAxis_Type, (Name)"0", _lucAxis_DefaultNew  );
   Stg_ComponentRegister_Add( componentRegister, lucTimeStep_Type, (Name)"0", _lucTimeStep_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, lucScalarFieldOnMeshCrossSection_Type, (Name)"0", _lucScalarFieldOnMeshCrossSection_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, lucScalarFieldOnMesh_Type, (Name)"0", _lucScalarFieldOnMesh_DefaultNew );
	

	/* Register Parents for type checking */
	RegisterParent( lucOpenGLDrawingObject_Type,             lucDrawingObject_Type );
	RegisterParent( lucCrossSection_Type,                    lucOpenGLDrawingObject_Type );
	RegisterParent( lucScalarFieldCrossSection_Type,         lucCrossSection_Type );
	RegisterParent( lucScalarFieldOnMeshCrossSection_Type,   lucCrossSection_Type );
	RegisterParent( lucVectorArrowCrossSection_Type,         lucCrossSection_Type );
	RegisterParent( lucEigenvectorsCrossSection_Type,        lucCrossSection_Type );
	RegisterParent( lucScalarField_Type,                     lucScalarFieldCrossSection_Type );	
	RegisterParent( lucScalarFieldOnMesh_Type,  	 	         lucScalarFieldOnMeshCrossSection_Type );
	RegisterParent( lucVectorArrows_Type,                    lucVectorArrowCrossSection_Type );
	RegisterParent( lucEigenvectors_Type,                    lucEigenvectorsCrossSection_Type );
	
	RegisterParent( lucColourBar_Type,                       lucDrawingObject_Type );
	RegisterParent( lucFieldVariableBorder_Type,             lucOpenGLDrawingObject_Type );
	RegisterParent( lucIsosurface_Type,                      lucOpenGLDrawingObject_Type );
	RegisterParent( lucTextureMap_Type,                      lucOpenGLDrawingObject_Type );
	RegisterParent( lucContour_Type,                         lucOpenGLDrawingObject_Type );
	RegisterParent( lucFeVariableSurface_Type,               lucOpenGLDrawingObject_Type );
	RegisterParent( lucSwarmViewerBase_Type,           lucOpenGLDrawingObject_Type );
	RegisterParent( lucSwarmViewer_Type,               lucSwarmViewerBase_Type );
	RegisterParent( lucSwarmVectors_Type,              lucSwarmViewerBase_Type );
	RegisterParent( lucSwarmSquares_Type,              lucSwarmViewerBase_Type );
	RegisterParent( lucHistoricalSwarmTrajectory_Type, lucOpenGLDrawingObject_Type );
	RegisterParent( lucSwarmRGBColourViewer_Type,      lucSwarmViewerBase_Type );
	RegisterParent( lucMeshViewer_Type,                lucOpenGLDrawingObject_Type );
	RegisterParent( lucTitle_Type,                     lucOpenGLDrawingObject_Type );
	RegisterParent( lucAxis_Type,                      lucOpenGLDrawingObject_Type );
	RegisterParent( lucTimeStep_Type,                  lucOpenGLDrawingObject_Type );




	return True;
}


