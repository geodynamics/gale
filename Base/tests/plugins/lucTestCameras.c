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
** $Id: lucTestCameras.c 740 2007-10-11 08:05:31Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>

void lucTestLeftDirection( lucCamera* camera, Stream* stream ) {
	XYZ   leftDirection;
	XYZ   focusDirection;

	lucCamera_GetLeftDirection( camera, leftDirection );
	lucCamera_GetFocusDirection( camera, focusDirection );

	Journal_Firewall( fabs( StGermain_AngleBetweenVectors( leftDirection, camera->upDirection, 3 ) - M_PI/2.0 ) < 0.01 , 
			stream, "Failed %s - left direction = %f %f %f, focusDirection = %f %f %f\n", __func__,
			leftDirection[0], leftDirection[1], leftDirection[2], 
			focusDirection[0], focusDirection[1], focusDirection[2] ); 
	Journal_Firewall( fabs( StGermain_AngleBetweenVectors( leftDirection, focusDirection, 3 ) - M_PI/2.0 ) < 0.01 , 
			stream, "Failed %s - left direction = %f %f %f, focusDirection = %f %f %f\n", __func__,
			leftDirection[0], leftDirection[1], leftDirection[2], 
			focusDirection[0], focusDirection[1], focusDirection[2] ); 

}

void lucTestBroadcast( lucCamera* camera, Stream* stream, AbstractContext* context ) {
	int rootRank = 0;
	double rank = (double) context->rank;

	camera->coord[ J_AXIS ] = rank;

	camera->eyeSeparation = rank;
	
	lucCamera_Broadcast( camera, rootRank, context->communicator );

	Journal_Firewall( fabs( camera->coord[ J_AXIS ] - (double) rootRank ) < 0.01, 
			stream, "Failed %s\n", __func__ );

	Journal_Firewall( fabs( camera->eyeSeparation - (double) rootRank ) < 0.01, 
			stream, "Failed %s\n", __func__ );
}

void lucTestEye( lucCamera* camera, Stream* stream ) {
	Coord leftEye;
	Coord rightEye;
	XYZ   vector;
	XYZ   leftDirection;

	camera->buffer = lucRight;
	lucCamera_SwapStereoBuffer( camera );
	lucCamera_CurrentEyePosition( camera, leftEye );
	lucCamera_SwapStereoBuffer( camera );
	lucCamera_CurrentEyePosition( camera, rightEye );

	StGermain_VectorSubtraction( vector, leftEye, rightEye, 3 );
	if ( camera->stereoType == lucMono ) {
		Journal_Firewall( fabs(StGermain_VectorMagnitude( vector, 3 ) ) < 0.01, 
				stream, "Failed %s\n", __func__ );
		StGermain_VectorSubtraction( vector, leftEye, camera->coord, 3 );
		Journal_Firewall( fabs(StGermain_VectorMagnitude( vector, 3 ) ) < 0.01, 
				stream, "Failed %s\n", __func__ );
		return;
	}

	Journal_Firewall( fabs(StGermain_VectorMagnitude( vector, 3 ) - camera->eyeSeparation ) < 0.01*camera->eyeSeparation, 
			stream, "Failed %s\n", __func__ );

	lucCamera_GetLeftDirection( camera, leftDirection );
	Journal_Firewall( fabs( StGermain_AngleBetweenVectors( vector, leftDirection, 3 ) ) < 0.01 , 
			stream, "Failed %s\n", __func__ );
}

void testCamera( lucCamera* camera, Stream* stream, AbstractContext* context ) {
	lucTestEye( camera, stream );
	lucTestLeftDirection( camera, stream );
}


void lucTestThisCamera( AbstractContext* context, lucCamera* camera, Stream* stream ) {
	Index                         angleCount = 8;                         
	Index                         angle_I;                         
	XYZ                           leftDirection;
	XYZ                           focusDirection;
	XYZ                           originalFocusDirection;
	
	Journal_Printf( stream, "************** Checking Camera '%s' **************\n", camera->name );
	Stream_Indent( stream );
	lucCamera_GetFocusDirection( camera, originalFocusDirection );

	Journal_Printf( stream, "************** Testing lucCamera_Zoom **************\n");
	lucCamera_Zoom( camera, 2.0 ) ;
	Stg_Class_Print( camera, stream );
	testCamera( camera, stream, context );
	lucCamera_Reset( camera );

	lucCamera_Zoom( camera, 0.25 ) ;
	Stg_Class_Print( camera, stream );
	testCamera( camera, stream, context );
	lucCamera_Reset( camera );

	lucTestBroadcast( camera, stream, context );
	lucCamera_Reset( camera );

	Journal_Printf( stream, "************** Testing lucCamera_RotateAroundUpDirection **************\n");
	Stream_Indent( stream );
	for ( angle_I = 0 ; angle_I < angleCount ; angle_I++ ) {
		lucCamera_RotateAroundUpDirection( camera, 2.0 * M_PI / (double)angleCount ) ;
		lucCamera_GetFocusDirection( camera, focusDirection );
		Journal_Printf(stream, "Has rotated %.5g degrees.\n", 
				StGermain_AngleBetweenVectors( focusDirection, originalFocusDirection, 3) * 180.0 / M_PI );
		lucCamera_GetLeftDirection( camera, leftDirection );
	}
	Stream_UnIndent( stream );
	lucCamera_Reset( camera );
	
	Journal_Printf( stream, "************** Testing lucCamera_RotateTowardsUpDirection **************\n");
	Stream_Indent( stream );
	for ( angle_I = 0 ; angle_I < angleCount ; angle_I++ ) {
		lucCamera_RotateTowardsUpDirection( camera,  0.5 * M_PI / (double)angleCount + M_PI/180.0 ) ;
		lucCamera_GetFocusDirection( camera, focusDirection );
		Journal_Printf(stream, "Has rotated %.5g degrees.\n", 
				StGermain_AngleBetweenVectors( focusDirection, originalFocusDirection, 3) * 180.0 / M_PI );
		testCamera( camera, stream, context );
	}
	lucCamera_Reset( camera );
	for ( angle_I = 0 ; angle_I < angleCount ; angle_I++ ) {
		lucCamera_RotateTowardsUpDirection( camera,  -0.5 * M_PI / (double)angleCount - M_PI/180.0 ) ;
		lucCamera_GetFocusDirection( camera, focusDirection );
		Journal_Printf(stream, "Has rotated %.5g degrees.\n", 
				StGermain_AngleBetweenVectors( focusDirection, originalFocusDirection, 3) * 180.0 / M_PI );

		/* Print left direction */
		testCamera( camera, stream, context );
	}
	Stream_UnIndent( stream );
	lucCamera_Reset( camera );
	Stream_UnIndent( stream );
}

void lucTestAllCameras( AbstractContext* context ) {
	lucCamera*                    stereo;
	lucCamera*                    frontOn;
	Stream*                       stream = Journal_Register( Info_Type, CURR_MODULE_NAME );

	if ( context->rank == 0 )
		Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, "camera.txt" );
	Stream_SetPrintingRank( stream, 0 );

	/* Note: this func gets added to construct extensions... may it should be changed to take in the tummy too? */
	frontOn = Stg_ComponentFactory_ConstructByName( context->CF, "frontOn", lucCamera, True, 0 /* dummy */ ); 
	stereo = Stg_ComponentFactory_ConstructByName( context->CF, "stereo", lucCamera, True, 0 /* dummy */ ); 
	Stg_Class_Print( frontOn, stream );
	Stg_Class_Print( stereo, stream );
	lucTestThisCamera( context, frontOn, stream );
	lucTestThisCamera( context, stereo, stream );
}

const Type lucTestCameras_Type = "lucTestCameras";
typedef struct {
	__Codelet
} lucTestCameras;

void _lucTestCameras_Construct( void* components, Stg_ComponentFactory* cf, void* data ) {
	AbstractContext* context;
	context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data ); 
	ContextEP_Append( context, AbstractContext_EP_ConstructExtensions, lucTestAllCameras );
}

void* _lucTestCameras_DefaultNew( Name name ) {
	return Codelet_New(
		lucTestCameras_Type,
		_lucTestCameras_DefaultNew,
		_lucTestCameras_Construct,
		_Codelet_Build,
		_Codelet_Initialise,
		_Codelet_Execute,
		_Codelet_Destroy,
		name );
}

Index lucTestCameras_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, lucTestCameras_Type, "0", _lucTestCameras_DefaultNew );
}

