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
** $Id: types.h 598 2006-08-01 03:11:05Z CecileDuboz $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __lucDrawingObjects_types_h__
#define __lucDrawingObjects_types_h__
	
	typedef struct lucOpenGLDrawingObject            lucOpenGLDrawingObject;

	typedef struct lucFieldVariableBorder            lucFieldVariableBorder;
	typedef struct lucScalarFieldCrossSection        lucScalarFieldCrossSection;
	typedef struct lucScalarField                    lucScalarField;
	typedef struct lucColourBar                      lucColourBar;
	typedef struct lucIsosurface                     lucIsosurface;
	typedef struct lucVectorArrows                   lucVectorArrows;
	typedef struct lucVectorArrowCrossSection        lucVectorArrowCrossSection;
	typedef struct lucEigenvectorsCrossSection       lucEigenvectorsCrossSection;
	typedef struct lucEigenvectors                   lucEigenvectors;
	typedef struct lucTextureMap                     lucTextureMap;
	typedef struct lucContour                        lucContour;
	typedef struct lucFeVariableSurface              lucFeVariableSurface;
	typedef struct lucSwarmViewerBase                lucSwarmViewerBase;
	typedef struct lucSwarmViewer                    lucSwarmViewer;
	typedef struct lucSwarmVectors                   lucSwarmVectors;
	typedef struct lucSwarmSquares                   lucSwarmSquares;
	typedef struct lucHistoricalSwarmTrajectory      lucHistoricalSwarmTrajectory;
	typedef struct lucSwarmRGBColourViewer           lucSwarmRGBColourViewer;
	typedef struct lucMeshViewer                     lucMeshViewer;
        typedef struct lucTitle                          lucTitle;
	typedef struct lucAxis                           lucAxis;
	typedef struct lucTimeStep                       lucTimeStep;
	typedef struct lucScalarFieldOnMeshCrossSection  lucScalarFieldOnMeshCrossSection;
	typedef struct lucScalarFieldOnMesh              lucScalarFieldOnMesh;


	

        

#endif
