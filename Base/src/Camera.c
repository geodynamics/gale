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
** $Id: Camera.c 628 2006-10-12 08:23:07Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include "types.h"
#include "Camera.h"
#include "Init.h"

#include <string.h>
#include <assert.h>

const Type lucCamera_Type = "lucCamera";

MPI_Datatype lucCamera_MPI_Datatype;

lucCamera* lucCamera_New( 
		Name                                               name,
		Coord                                              coord, 
		Coord                                              focalPoint,
		XYZ                                                upDirection,
		double                                             focalLength,
		double                                             aperture,
		double                                             eyeSeparation,
		lucStereoType                                      stereoType,
		FieldVariable*                                     centreFieldVariable )
{
	lucCamera* self = (lucCamera*) _lucCamera_DefaultNew( name );

	lucCamera_InitAll( self, coord, focalPoint, upDirection, focalLength, aperture, eyeSeparation, stereoType, centreFieldVariable );

	return self;
}

lucCamera* _lucCamera_New(
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
		Name                                               name )
{
	lucCamera*    self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucCamera) );
	self = (lucCamera*) _Stg_Component_New( 
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
			name, 
			NON_GLOBAL );
	
	return self;
}

void lucCamera_Init(		
		lucCamera*                                         self,
		Coord                                              coord, 
		Coord                                              focalPoint,
		XYZ                                                upDirection,
		double                                             focalLength,
		double                                             aperture,
		double                                             eyeSeparation,
		lucStereoType                                      stereoType,
		FieldVariable*                                     centreFieldVariable )
{
	memcpy( self->coord, coord, sizeof(Coord) );
	memcpy( self->focalPoint, focalPoint, sizeof(Coord) );
	memcpy( self->upDirection, upDirection, sizeof(XYZ) );

	self->focalLength         = focalLength;
	self->aperture            = aperture;
	self->eyeSeparation       = eyeSeparation;
	self->stereoType          = stereoType;
	self->centreFieldVariable = centreFieldVariable;

	/* Store Original Values */
	self->needsToDraw   = True;
	self->buffer        = lucLeft;
	self->originalCamera = lucCamera_Copy( self );
}

void lucCamera_InitAll( 
		void*                                              camera,
		Coord                                              coord, 
		Coord                                              focalPoint,
		XYZ                                                upDirection,
		double                                             focalLength,
		double                                             aperture,
		double                                             eyeSeparation,
		lucStereoType                                      stereoType,
		FieldVariable*                                     centreFieldVariable )
{
	lucCamera* self        = camera;

	/* TODO Init parent */
	lucCamera_Init( self, coord, focalPoint, upDirection, focalLength, aperture, eyeSeparation, stereoType, centreFieldVariable );
}

void _lucCamera_Delete( void* camera ) {
	lucCamera* self        = camera;
	
	if ( self->originalCamera != NULL );
		Stg_Class_Delete( self->originalCamera );
	
	_Stg_Component_Delete( self );
}

void _lucCamera_Print( void* camera, Stream* stream ) {
	lucCamera* self        = camera;
	
	Journal_Printf( stream, "lucCamera: %s\n", self->name );

	Stream_Indent( stream );

	/* Print Parent */
	_Stg_Component_Print( self, stream );

	Journal_PrintArray( stream, self->coord, 3 );
	Journal_PrintArray( stream, self->focalPoint, 3 );
	Journal_PrintArray( stream, self->upDirection, 3 );
	
	Journal_PrintValue( stream, self->aperture );

	Journal_Printf( stream, "self->stereoType = ");
	switch ( self->stereoType ) {
		case lucMono:
			Journal_Printf( stream, "lucMono\n" ); break;
		case lucStereoToeIn:
			Journal_Printf( stream, "lucStereoToeIn\n" ); break;
		case lucStereoAsymmetric:
			Journal_Printf( stream, "lucStereoAsymmetric\n" ); break;
	}
	
	Journal_Printf( stream, "self->buffer = ");
	switch ( self->buffer ) {
		case lucLeft:
			Journal_Printf( stream, "lucLeft\n" ); break;
		case lucRight:
			Journal_Printf( stream, "lucRight\n" ); break;
	}
	Journal_PrintValue( stream, self->eyeSeparation );

	Stream_UnIndent( stream );
}

void* _lucCamera_Copy( void* camera, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	lucCamera* self        = camera;
	lucCamera* newCamera;

	newCamera = _Stg_Component_Copy( self, dest, deep, nameExt, ptrMap );

	memcpy( newCamera->coord,       self->coord,       sizeof(Coord) );
	memcpy( newCamera->focalPoint,  self->focalPoint,  sizeof(Coord) );
	memcpy( newCamera->upDirection, self->upDirection, sizeof(XYZ) );
	
	newCamera->aperture      = self->aperture;
	newCamera->stereoType    = self->stereoType;
	newCamera->eyeSeparation = self->eyeSeparation;
	newCamera->buffer        = self->buffer;
	
	return (void*) newCamera;
}

void* _lucCamera_DefaultNew( Name name ) {
	return _lucCamera_New( 
			sizeof( lucCamera ),
			lucCamera_Type,
			_lucCamera_Delete,
			_lucCamera_Print,
			_lucCamera_Copy,
			_lucCamera_DefaultNew,
			_lucCamera_Construct,
			_lucCamera_Build,
			_lucCamera_Initialise,
			_lucCamera_Execute,
			_lucCamera_Destroy,
			name );
}

void _lucCamera_Construct( void* camera, Stg_ComponentFactory* cf, void* data ) {
	lucCamera*             self               = (lucCamera*) camera;
	Coord                  coord;
	Coord                  focalPoint;
	XYZ                    upDirection;
	FieldVariable*         centreFieldVariable;
	double                 focalLength;
	double                 aperture;
	double                 eyeSeparation;
	lucStereoType          stereoType;
	Name                   stereoTypeName;

  #define DTOR 0.0174532925
	focalPoint[I_AXIS]  = Stg_ComponentFactory_GetDouble( cf, self->name, "focalPointX", 0.0 );
	focalPoint[J_AXIS]  = Stg_ComponentFactory_GetDouble( cf, self->name, "focalPointY", 0.0 );
	focalPoint[K_AXIS]  = Stg_ComponentFactory_GetDouble( cf, self->name, "focalPointZ", 0.0 );

	upDirection[I_AXIS] = Stg_ComponentFactory_GetDouble( cf, self->name, "upDirectionX", 0.0 );
	upDirection[J_AXIS] = Stg_ComponentFactory_GetDouble( cf, self->name, "upDirectionY", 1.0 );
	upDirection[K_AXIS] = Stg_ComponentFactory_GetDouble( cf, self->name, "upDirectionZ", 0.0 );
	
	focalLength = Stg_ComponentFactory_GetDouble( cf, self->name, "focalLength", 0.0 );

	coord[I_AXIS]  = Stg_ComponentFactory_GetDouble( cf, self->name, "coordX", 0.0 );
	coord[J_AXIS]  = Stg_ComponentFactory_GetDouble( cf, self->name, "coordY", 0.0 );
	coord[K_AXIS]  = Stg_ComponentFactory_GetDouble( cf, self->name, "coordZ", 1.0 );

	aperture = Stg_ComponentFactory_GetDouble( cf, self->name, "aperture", 45.0 );

	/* Get Stereo Type */
	stereoTypeName = Stg_ComponentFactory_GetString( cf, self->name, "stereoType", "lucMono" );
	if ( strcasecmp( stereoTypeName, "ToeIn" ) == 0 ) 
		stereoType = lucStereoToeIn;
	else if ( strcasecmp( stereoTypeName, "Asymmetric" ) == 0 ) 
		stereoType = lucStereoAsymmetric;
	else 
		stereoType = lucMono;

	eyeSeparation  = Stg_ComponentFactory_GetDouble( cf, self->name, "eyeSeparation", 0.2 );

	centreFieldVariable = Stg_ComponentFactory_ConstructByKey( cf, self->name, "CentreFieldVariable", FieldVariable,False,data);

	lucCamera_Init( self, coord, focalPoint, upDirection, focalLength, aperture, eyeSeparation, stereoType, centreFieldVariable );
}

void _lucCamera_Build( void* camera, void* data ) { }
void _lucCamera_Initialise( void* camera, void* data ) { }
void _lucCamera_Execute( void* camera, void* data ) { }
void _lucCamera_Destroy( void* camera, void* data ) { }

void lucCamera_Zoom( void* camera, double zoomFactor ) {
	lucCamera* self                = camera;
	XYZ        vectorFocusToCamera;
	
	StGermain_VectorSubtraction( vectorFocusToCamera, self->coord, self->focalPoint, 3 );
	vectorFocusToCamera[ I_AXIS ] *= zoomFactor;
	vectorFocusToCamera[ J_AXIS ] *= zoomFactor;
	vectorFocusToCamera[ K_AXIS ] *= zoomFactor;
	StGermain_VectorAddition( self->coord, vectorFocusToCamera, self->focalPoint, 3 );
	
	self->needsToDraw = True;
}

void lucCamera_RotateAroundUpDirection( void* camera, double deltaTheta ) {
	lucCamera* self                = camera;
	XYZ        vectorFocusToCamera;
	XYZ        rotatedVectorFocusToCamera;

	StGermain_VectorSubtraction( vectorFocusToCamera, self->coord, self->focalPoint, 3 );
	StGermain_RotateVector( rotatedVectorFocusToCamera, vectorFocusToCamera, self->upDirection, deltaTheta );
	StGermain_VectorAddition( self->coord, rotatedVectorFocusToCamera, self->focalPoint, 3 );
	
	self->needsToDraw = True;
}

void lucCamera_RotateTowardsUpDirection( void* camera, double deltaTheta ) {
	lucCamera*   self                = camera;
	XYZ          rotationAxis;
	XYZ          vectorFocusToCamera;
	XYZ          rotatedVectorFocusToCamera;
	double       angleBetween;
	double       predictedAngle;
	const double minAngle            = M_PI/180.0 * 0.1; /* A tenth of a degree */
	const double maxAngle            = M_PI - minAngle;

	StGermain_VectorSubtraction( vectorFocusToCamera, self->coord, self->focalPoint, 3 );
	StGermain_VectorCrossProduct( rotationAxis, vectorFocusToCamera, self->upDirection );
	StGermain_VectorNormalise( rotationAxis, 3 );

	/* Check if angle is too large */
	angleBetween = StGermain_AngleBetweenVectors( vectorFocusToCamera, self->upDirection, 3 );
	predictedAngle = angleBetween - deltaTheta;
	if ( predictedAngle < minAngle ) 
		deltaTheta = angleBetween - minAngle;
	else if ( predictedAngle > maxAngle ) 
		deltaTheta = angleBetween - maxAngle;
		
	StGermain_RotateVector( rotatedVectorFocusToCamera, vectorFocusToCamera, rotationAxis, deltaTheta );
	
	StGermain_VectorAddition( self->coord, rotatedVectorFocusToCamera, self->focalPoint, 3 );
	self->needsToDraw = True;
}

void lucCamera_Reset( void* camera ) {
	lucCamera* self                = camera;
	lucCamera* originalCamera      = self->originalCamera;
	
	Stg_Class_Copy( originalCamera, self, False, NULL, NULL );
	self->originalCamera = originalCamera;
	self->needsToDraw = True;
}

void lucCamera_SetOriginal( void* camera ) {
	lucCamera* self                = camera;
	
	Stg_Class_Copy( self, self->originalCamera, False, NULL, NULL );
	self->originalCamera->originalCamera = NULL;
}

void lucCamera_Broadcast( void* camera, int rootRank, MPI_Comm comm ) {
	lucCamera* self                = camera;

	Journal_DPrintfL( lucDebug, 2, "In %s for %s '%s'.\n", __func__, self->type, self->name );

	MPI_Bcast( self, 1, lucCamera_MPI_Datatype, rootRank, comm );
}

void lucCamera_GetFocusDirection( void* camera, XYZ focusDirection ) {
	lucCamera*   self                = camera;
	
	StGermain_VectorSubtraction( focusDirection, self->focalPoint, self->coord, 3 );
	StGermain_VectorNormalise( focusDirection, 3 );
}

void lucCamera_GetLeftDirection( void* camera, XYZ leftDirection ) {
	lucCamera*   self                = camera;
	XYZ          focusDirection;
	
	lucCamera_GetFocusDirection( self, focusDirection );
	StGermain_VectorCrossProduct( leftDirection, self->upDirection, focusDirection );
	StGermain_VectorNormalise( leftDirection, 3 );
}

void lucCamera_SwapStereoBuffer( void* camera ) {
	lucCamera* self                = camera;

	if ( self->buffer == lucLeft )
		self->buffer = lucRight;
	else 
		self->buffer = lucLeft;
	self->needsToDraw = True;
}

void lucCamera_CurrentEyePosition( void* camera, Coord currEyePos ) {
	lucCamera* self                = camera;
	XYZ        leftDirection;
	double     factor;

	memcpy( currEyePos, self->coord, sizeof(Coord) );
	if ( self->stereoType == lucMono )
		return;

	lucCamera_GetLeftDirection( camera, leftDirection );

	factor = 0.5 * self->eyeSeparation;
	if ( self->buffer == lucRight )
		factor *= -1.0; 

	currEyePos[ I_AXIS ] += factor * leftDirection[ I_AXIS ];
	currEyePos[ J_AXIS ] += factor * leftDirection[ J_AXIS ];
	currEyePos[ K_AXIS ] += factor * leftDirection[ K_AXIS ];
}


void lucCamera_CentreFromFieldVariable( void* camera ) {
	lucCamera*     self = (lucCamera*) camera;
	Coord          min, max;
	Coord          oldFocalPoint;
	XYZ            difference = {0.0,0.0,0.0};
	FieldVariable* fieldVariable = self->centreFieldVariable;

	if ( !fieldVariable )
		return;

	FieldVariable_GetMinAndMaxGlobalCoords( fieldVariable, min, max );
  
		
	/* Set up focal point */
	memcpy( oldFocalPoint, self->focalPoint, sizeof(Coord) );
	self->focalPoint[I_AXIS] = 0.5 * (min[ I_AXIS ] + max[ I_AXIS ]);
	self->focalPoint[J_AXIS] = 0.5 * (min[ J_AXIS ] + max[ J_AXIS ]);
	if (fieldVariable->dim == 2) 
		self->focalPoint[K_AXIS] = 0.0;
	else 
		self->focalPoint[K_AXIS] = 0.5 * (min[ K_AXIS ] + max[ K_AXIS ]);
	
	/* Set up camera */
	StGermain_VectorSubtraction( difference, self->focalPoint, oldFocalPoint, fieldVariable->dim );
	self->coord[ I_AXIS ] += difference[ I_AXIS ] ;
	self->coord[ J_AXIS ] += difference[ J_AXIS ] ;
	self->coord[ K_AXIS ] += difference[ K_AXIS ] ;

	/* Adjust Original Camera */
	self->originalCamera->coord[ I_AXIS ] += difference[ I_AXIS ] ;
	self->originalCamera->coord[ J_AXIS ] += difference[ J_AXIS ] ;
	self->originalCamera->coord[ K_AXIS ] += difference[ K_AXIS ] ;
	memcpy( self->originalCamera->focalPoint, self->focalPoint, sizeof( Coord ) );
}

void lucCamera_ChangeFocalPoint( void* camera, Pixel_Index startx, Pixel_Index starty, Pixel_Index xpos, Pixel_Index ypos ){
	lucCamera*     self = (lucCamera*) camera;
	XYZ             leftDirection;
	Dimension_Index dim_I;

	lucCamera_GetLeftDirection( camera, leftDirection );
	for ( dim_I = 0 ; dim_I < 3 ; dim_I++ ) {
		self->focalPoint[ dim_I ] -= 0.01 * ((double)xpos - (double)startx) * leftDirection[ dim_I ];
		self->focalPoint[ dim_I ] -= 0.01 * ((double)ypos - (double)starty) * self->upDirection[ dim_I ];
	}
	memcpy( self->originalCamera->focalPoint, self->focalPoint, sizeof( Coord ) );

	self->needsToDraw = True;
}

void lucCamera_Pickle( void* camera, Stream* stream ) {
	lucCamera* self                = camera;
	       
	Journal_Printf( stream, "<struct name=\"%s\">\n", self->name);
	Stream_Indent( stream );

	Journal_Printf( stream, "<param name=\"Type\">%s</param>\n", self->type );
	Journal_Printf( stream, "<param name=\"coordX\">%.5g</param>\n", self->coord[ I_AXIS ] );
	Journal_Printf( stream, "<param name=\"coordY\">%.5g</param>\n", self->coord[ J_AXIS ] );
	Journal_Printf( stream, "<param name=\"coordZ\">%.5g</param>\n", self->coord[ K_AXIS ] );

	Journal_Printf( stream, "<param name=\"focalPointX\">%.5g</param>\n", self->focalPoint[ I_AXIS ] );
	Journal_Printf( stream, "<param name=\"focalPointY\">%.5g</param>\n", self->focalPoint[ J_AXIS ] );
	Journal_Printf( stream, "<param name=\"focalPointZ\">%.5g</param>\n", self->focalPoint[ K_AXIS ] );

	Journal_Printf( stream, "<param name=\"upDirectionX\">%.5g</param>\n", self->upDirection[ I_AXIS ] );
	Journal_Printf( stream, "<param name=\"upDirectionY\">%.5g</param>\n", self->upDirection[ J_AXIS ] );
	Journal_Printf( stream, "<param name=\"upDirectionZ\">%.5g</param>\n", self->upDirection[ K_AXIS ] );

	Journal_Printf( stream, "<param name=\"stereoType\">%s</param>\n", 
			self->stereoType == lucStereoToeIn      ? "ToeIn" :
			self->stereoType == lucStereoAsymmetric ? "Asymmetric" : "Mono" );
	
	Journal_Printf( stream, "<param name=\"eyeSeparation\">%.5g</param>\n", self->eyeSeparation );
	Journal_Printf( stream, "<param name=\"focalLength\">%.5g</param>\n", self->focalLength );
	Stream_UnIndent( stream );
	Journal_Printf( stream, "</struct>\n");
}
void lucCamera_SetNeedsToDraw( void * camera ){
	lucCamera* self = (lucCamera*) camera;
	self->needsToDraw = True;
}

#define lucCamera_TypesCount 9
void lucCamera_Create_MPI_Datatype() {
	MPI_Datatype        typeList[lucCamera_TypesCount]     = 
		{ MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT };
	int                 blocklen[lucCamera_TypesCount]     = {3, 3, 3, 1, 1, 1, 1, 1, 1};
	MPI_Aint            displacement[lucCamera_TypesCount];
	lucCamera           camera;

	displacement[0] = GetOffsetOfMember( camera, coord );
	displacement[1] = GetOffsetOfMember( camera, focalPoint );
	displacement[2] = GetOffsetOfMember( camera, upDirection );
	displacement[3] = GetOffsetOfMember( camera, focalLength );
	displacement[4] = GetOffsetOfMember( camera, aperture );
	displacement[5] = GetOffsetOfMember( camera, eyeSeparation );
	displacement[6] = GetOffsetOfMember( camera, buffer );
	displacement[7] = GetOffsetOfMember( camera, stereoType );
	displacement[8] = GetOffsetOfMember( camera, needsToDraw );
	
	MPI_Type_struct( lucCamera_TypesCount, blocklen, displacement, typeList, &lucCamera_MPI_Datatype );
	MPI_Type_commit( & lucCamera_MPI_Datatype );
}
