/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**	Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
** Contributors:
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Lesser General Public
**  License as published by the Free Software Foundation; either
**  version 2.1 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License along with this library; if not, write to the Free Software
**  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**
** $Id: StandardConditionFunctions.c 971 2007-10-30 05:24:15Z DavidLee $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <assert.h>
#include "StandardConditionFunctions.h"

const Type StgFEM_StandardConditionFunctions_Type = "StgFEM_StandardConditionFunctions";

void _StgFEM_StandardConditionFunctions_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	AbstractContext*        context;
	ConditionFunction*      condFunc;

	context = (AbstractContext*)Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data ); 
	
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SolidBodyRotation, "Velocity_SolidBodyRotation" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );
	
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_PartialRotationX, "Velocity_PartialRotationX" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );
		
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_PartialRotationY, "Velocity_PartialRotationY" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );
	
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SimpleShear, "Velocity_SimpleShear" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_Extension, "Velocity_Extension" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_PartialLid_TopLayer, "Velocity_PartialLid_TopLayer" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_Trigonometry, "Temperature_Trigonometry" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_LinearInterpolationLid, "Velocity_LinearInterpolationLid" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_Lid_RampWithCentralMax, "Velocity_Lid_RampWithCentralMax" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SinusoidalLid, "Velocity_SinusoidalLid" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_CornerOnly, "Velocity_Lid_CornerOnly" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );
	
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_ConvectionBenchmark, "Temperature_ConvectionBenchmark" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );
	
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_LinearWithSinusoidalPerturbation, "LinearWithSinusoidalPerturbation" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_EdgeDriveConvectionIC, "EdgeDriveConvectionIC" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );
	
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_AnalyticalTemperatureIC, "AnalyticalTemperatureIC" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( Stg_FEM_VelicTemperatureIC, "VelicTemperatureIC");
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( Stg_FEM_VelicTemperatureIC_SolB, "VelicTemperatureIC_SolB");
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );
	
	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SinusoidalExtension, "SinusoidalExtension");
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_StepFunction, "StepFunction");
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SpecRidge3D, "SpecRidge3D");
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SpectralBCX, "SpectralBCX");
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SpectralBCY, "SpectralBCY");
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SpectralBCZ, "SpectralBCZ");
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SpectralPressureBCX, "SpectralPressureBCX");
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SpectralPressureBCY, "SpectralPressureBCY");
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_ErrorFunc, "ErrorFunc");
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );

	condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_ConstantVelocity, "ConstantVelocity");
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );
}

void* _StgFEM_StandardConditionFunctions_DefaultNew( Name name ) {
	return Codelet_New(
		StgFEM_StandardConditionFunctions_Type,
		_StgFEM_StandardConditionFunctions_DefaultNew,
		_StgFEM_StandardConditionFunctions_Construct,
		_Codelet_Build,
		_Codelet_Initialise,
		_Codelet_Execute,
		_Codelet_Destroy,
		name );
}

Index StgFEM_StandardConditionFunctions_Register( PluginsManager* pluginsManager ) {
	Journal_DPrintf( StgFEM_Debug, "In: %s( void* )\n", __func__ );

	return PluginsManager_Submit( pluginsManager, StgFEM_StandardConditionFunctions_Type, "0", _StgFEM_StandardConditionFunctions_DefaultNew );
}

void StgFEM_StandardConditionFunctions_SolidBodyRotation( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	FeVariable*             feVariable         = NULL;
	FeMesh*			mesh               = NULL;
	double*                 result             = (double*) _result;
	double*                 coord;
	Coord                   centre;
	Coord                   vector;
	double                  omega;
	
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;

	/* Find Centre of Solid Body Rotation */
	centre[ I_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationCentreX", 0.0 );
	centre[ J_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationCentreY", 0.0 );
	centre[ K_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationCentreZ", 0.0 );
	omega            = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationOmega",   1.0 );

	/* Find coordinate of node */
	coord = Mesh_GetVertex( mesh, node_lI );

	/* Find vector from centre to node */
	StGermain_VectorSubtraction( vector, coord, centre, 2 );

	result[ I_AXIS ] = -omega * vector[ J_AXIS ];
	result[ J_AXIS ] =  omega * vector[ I_AXIS ];
}


void StgFEM_StandardConditionFunctions_PartialRotationX( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	FeVariable*             feVariable         = NULL;
	FeMesh*			mesh               = NULL;
	double*                 result             = (double*) _result;
	double*                 coord;
	Coord                   centre;
	Coord                   vector;
	double                  omega;
	double			size;

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;

	/* Find Centre of Solid Body Rotation */
	centre[ I_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationCentreX", 0.0 );
	centre[ J_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationCentreY", 0.0 );
	centre[ K_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationCentreZ", 0.0 );
	size             = Dictionary_GetDouble_WithDefault( dictionary, "RadiusCylinder", 0.0 );
	size += 0.1;
	omega            = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationOmega",   1.0 );

	/* Find coordinate of node */
	coord = Mesh_GetVertex( mesh, node_lI );

	/* Find vector from centre to node */
	StGermain_VectorSubtraction( vector, coord, centre, 2 );

	/*if (context->currentTime > 1.33e-6)
	  omega=0.0;*/
	
	if ((vector[ I_AXIS ]*vector[ I_AXIS ]+vector[ J_AXIS ]*vector[ J_AXIS ])<=size*size)
		*result = -omega * vector[ J_AXIS ];
	else
		*result = 0.0;
}

void StgFEM_StandardConditionFunctions_PartialRotationY( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	FeVariable*             feVariable         = NULL;
	FeMesh*			mesh               = NULL;
	double*                 result             = (double*) _result;
	double*                 coord;
	Coord                   centre;
	Coord                   vector;
	double                  omega;
	double			size;

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;

	/* Find Centre of Solid Body Rotation */
	centre[ I_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationCentreX", 0.0 );
	centre[ J_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationCentreY", 0.0 );
	centre[ K_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationCentreZ", 0.0 );
	size             = Dictionary_GetDouble_WithDefault( dictionary, "RadiusCylinder", 0.0 );
	size += 0.1;
	omega            = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationOmega",   1.0 );

	/* Find coordinate of node */
	coord = Mesh_GetVertex( mesh, node_lI );

	/* Find vector from centre to node */
	StGermain_VectorSubtraction( vector, coord, centre, 2 );

	/*if (context->currentTime > 1.33e-6)
	  omega=0.0;*/
	
	if ((vector[ I_AXIS ]*vector[ I_AXIS ]+vector[ J_AXIS ]*vector[ J_AXIS ])<=size*size)
		*result =  omega * vector[ I_AXIS ];
	else 
		*result = 0.0;
}


void StgFEM_StandardConditionFunctions_SimpleShear( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	FeVariable*             feVariable         = NULL;
	FeMesh*			mesh               = NULL;
	double*                 result             = (double*) _result;
	double*                 coord;
	double                  centre;
	double                  factor;
	
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;

	/* Find Centre of Solid Body Rotation */
	centre = Dictionary_GetDouble_WithDefault( dictionary, "SimpleShearCentreY", 0.0 );
	factor = Dictionary_GetDouble_WithDefault( dictionary, "SimpleShearFactor", 1.0 );

	/* Find coordinate of node */
	coord = Mesh_GetVertex( mesh, node_lI );

	*result = factor * (coord[ J_AXIS ] - centre);
}

void StgFEM_StandardConditionFunctions_Extension( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	FeVariable*             feVariable         = NULL;
	FeMesh*			mesh               = NULL;
	double*                 result             = (double*) _result;
	double*                 coord;
	double                  centre;
	double                  factor;
	
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;

	/* Find Centre of Solid Body Rotation */
	centre = Dictionary_GetDouble_WithDefault( dictionary, "ExtensionCentreX", 0.0 );
	factor = Dictionary_GetDouble_WithDefault( dictionary, "ExtensionFactor", 1.0 );

	/* Find coordinate of node */
	coord = Mesh_GetVertex( mesh, node_lI );

	*result = factor * (coord[ I_AXIS ] - centre);
}


void StgFEM_StandardConditionFunctions_PartialLid_TopLayer( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* result ) {
	DomainContext*	context = (DomainContext*)_context;
	FeVariable*             velVar = NULL;
	FeMesh*			mesh = NULL;
	double*			velResult = (double*)result;
	double                  margin = 0;
	double			min[3], max[3];
	
	velVar = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh = velVar->feMesh;

	Mesh_GetMinimumSeparation( mesh, &margin, NULL );
	Mesh_GetGlobalCoordRange( mesh, min, max );
	margin *= 1.1;
	if( (Mesh_GetVertex( mesh, node_lI )[I_AXIS] < (max[I_AXIS] - margin )) && 
	    (Mesh_GetVertex( mesh, node_lI )[I_AXIS] > (min[I_AXIS] + margin )))
	{
		(*velResult) = 1;
	}
	else {
		(*velResult) = 0;
	}
}

void StgFEM_StandardConditionFunctions_LinearInterpolationLid( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* result ) {
	DomainContext*	context = (DomainContext*)_context;
	FeVariable*             velVar = NULL;
	FeMesh*			mesh = NULL;
	double*			velResult = (double*)result;
	double			boxLength = 0;
	double			leftHandSideValue = 0;
	double			rightHandSideValue = 0;
	double			gradient = 0;
	double			min[3], max[3];
	
	velVar = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh = velVar->feMesh;

	Mesh_GetGlobalCoordRange( mesh, min, max );
	boxLength = max[I_AXIS] - min[I_AXIS];
	leftHandSideValue = Dictionary_GetDouble_WithDefault( context->dictionary, "bcLeftHandSideValue", 0.0 );
	rightHandSideValue = Dictionary_GetDouble_WithDefault( context->dictionary, "bcRightHandSideValue", 1.0 );
	gradient = (rightHandSideValue - leftHandSideValue) / boxLength;
	(*velResult) = leftHandSideValue + gradient * (Mesh_GetVertex( mesh, node_lI )[I_AXIS] - min[I_AXIS] );
}


void StgFEM_StandardConditionFunctions_Lid_RampWithCentralMax( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* result ) {
	DomainContext*	context = (DomainContext*)_context;
	FeVariable*             velVar = NULL;
	FeMesh*			mesh = NULL;
	double*			velResult = (double*)result;
	double			boxLength = 0;
	double			xPosRelativeToTopLeft = 0;
	double			min[3], max[3];
	
	velVar = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh = velVar->feMesh;

	Mesh_GetGlobalCoordRange( mesh, min, max );
	xPosRelativeToTopLeft = Mesh_GetVertex( mesh, node_lI )[I_AXIS] - min[I_AXIS];
	boxLength = max[I_AXIS] - min[I_AXIS];
	if ( xPosRelativeToTopLeft < boxLength / 2 ) {
		(*velResult) =  2 * xPosRelativeToTopLeft / boxLength;
	}
	else {
		(*velResult) = 1 - 2 * ( xPosRelativeToTopLeft - (boxLength/2) );
	}
}


void StgFEM_StandardConditionFunctions_SinusoidalLid( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* result ) {
	DomainContext*	context = (DomainContext*)_context;
	FeVariable*             velVar = NULL;
	FeMesh*			mesh = NULL;
	double*			velResult = (double*)result;
	double			boxLength = 0;
	double			linearInterp = 0;
	double          	wavenumber;
	double			min[3], max[3];

	wavenumber = Dictionary_GetDouble_WithDefault( context->dictionary, "sinusoidalLidWavenumber", 1 );
	
	velVar = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh = velVar->feMesh;

	Mesh_GetGlobalCoordRange( mesh, min, max );
	boxLength = max[I_AXIS] - min[I_AXIS];
	linearInterp = (Mesh_GetVertex( mesh, node_lI )[I_AXIS] - min[I_AXIS] ) / boxLength;
	(*velResult) = sin( linearInterp * M_PI * wavenumber );
}


void StgFEM_StandardConditionFunctions_CornerOnly( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* result ) {
	DomainContext*	context = (DomainContext*)_context;
	FeVariable*             velVar = NULL;
	FeMesh*			feMesh = NULL;
	double*			velResult = (double*)result;
	Node_GlobalIndex	node_gI = 0;
	unsigned		inds[3];
	Grid*			elGrid;

	velVar = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	feMesh = velVar->feMesh;
	elGrid = *(Grid**)ExtensionManager_Get( feMesh->info, feMesh, 
						ExtensionManager_GetHandle( feMesh->info, "elGrid" ) );

	node_gI = Mesh_DomainToGlobal( feMesh, MT_VERTEX, node_lI );
	RegularMeshUtils_Node_1DTo3D( feMesh, node_gI, inds );
	
	if ( inds[0] == elGrid->sizes[I_AXIS] ) {
		(*velResult) = 1;
	}
	else {
		(*velResult) = 0;
	}
}

void StgFEM_StandardConditionFunctions_LinearWithSinusoidalPerturbation( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	FeVariable*             feVariable = NULL;
	FeMesh*			feMesh = NULL;
	unsigned		nDims;
	double*                 result = (double*) _result;
	double                  topLayerBC;
	double                  bottomLayerBC;
	double                  perturbationAmplitude;
	double                  horizontalWaveNumber;
	double                  verticalWaveNumber;
	double                  scaleFactor;
	double*                 coord;
	Dimension_Index         dim_I=0;
	Coord                   relScaledCoord;
	double			min[3], max[3];
	
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	feMesh       = feVariable->feMesh;

	nDims = Mesh_GetDimSize( feMesh );
	Mesh_GetGlobalCoordRange( feMesh, min, max );

	topLayerBC = Dictionary_GetDouble_WithDefault( dictionary, "SinusoidalTempIC_TopLayerBC", 0.0 );
	bottomLayerBC = Dictionary_GetDouble_WithDefault( dictionary, "SinusoidalTempIC_BottomLayerBC", 1.0 );
	scaleFactor = bottomLayerBC - topLayerBC;
	perturbationAmplitude = Dictionary_GetDouble_WithDefault( dictionary, "SinusoidalTempIC_PerturbationAmplitude", 0.1 );
	/* Note: these are both multiplied by pi, so wavenumber = 1 means the perturbation goes from 0 to pi, which is
	 * half a full sin or cos cycle. Wavenumber = 3 means the range is 0 -> 3pi, or 1 and a half full cycles. */
	horizontalWaveNumber = Dictionary_GetDouble_WithDefault( dictionary, "SinusoidalTempIC_HorizontalWaveNumber", 1.0 );
	verticalWaveNumber = Dictionary_GetDouble_WithDefault( dictionary, "SinusoidalTempIC_VerticalWaveNumber", 1.0 );

	coord = Mesh_GetVertex( feMesh, node_lI );
	/* make coord relative to box bottom left corner, then scale from 0 to 1 between box min & max */
	for( dim_I = 0; dim_I < nDims; dim_I++ ) {
		relScaledCoord[dim_I] = (coord[dim_I] - min[dim_I]) / (max[dim_I] - min[dim_I]);
	}

	/* Note: ok to use the 1.0 below since we've already scaled the coord to somewhere between 0 to 1 */
	*result = topLayerBC + scaleFactor * ( 1.0 - relScaledCoord[ J_AXIS ] )
		+ perturbationAmplitude * ( cos( horizontalWaveNumber * M_PI * coord[ I_AXIS ] )
					    * sin( verticalWaveNumber * M_PI * coord[ J_AXIS ] ) );
}

void StgFEM_StandardConditionFunctions_Trigonometry( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*			feMesh               = NULL;
	double*                 result             = (double*) _result;
	double*                 coord;
	double                  height, width;
	double			min[3], max[3];

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	feMesh       = feVariable->feMesh;

	Mesh_GetGlobalCoordRange( feMesh, min, max );
	coord = Mesh_GetVertex( feMesh, node_lI );

	/* Get Aspect Ratio */
	height = max[ J_AXIS ] - min[ J_AXIS ];
	width  = max[ I_AXIS ] - min[ I_AXIS ];
	
	*result = 1.0 - 0.5 * M_PI * coord[ J_AXIS ] * sin( M_PI * coord[ I_AXIS ]/width );
}

#define SMALL 1.0e-5
void Stg_FEM_VelicTemperatureIC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*  context            = (DomainContext*)_context;
	FeVariable*             temperatureField   = (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	FeMesh*			feMesh               = temperatureField->feMesh;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double*                 coord;
	double                  x; 
	double                  y;
	double                  kx;
	double                  ky;
	int                     wavenumberX;
	double                  wavenumberY;
	double                  sigma;
	double                  Lx;
	double			min[3], max[3];
	
	/* Find coordinate of node */
	coord = Mesh_GetVertex( feMesh, node_lI );
	Mesh_GetGlobalCoordRange( feMesh, min, max );

	/* Make sure that the box has right dimensions */
	assert( ( max[ J_AXIS ] - min[ J_AXIS ] - 1.0 ) < SMALL );
	Lx = max[ I_AXIS ] - min[ I_AXIS ];

	x = coord[ I_AXIS ] - min[ I_AXIS ];
	y = coord[ J_AXIS ] - min[ J_AXIS ];

	wavenumberX = Dictionary_GetInt_WithDefault( dictionary, "wavenumberX", 1 );
	wavenumberY = Dictionary_GetDouble_WithDefault( dictionary, "wavenumberY", 1.0 );
	sigma = Dictionary_GetDouble_WithDefault( dictionary, "sigma", 1.0 );
	
	assert( sigma > 0.0 );
	assert( wavenumberY > 0.0 );
	assert( wavenumberX > 0.0 );
	
	kx = (double)wavenumberX * M_PI / Lx;
	ky = (double)wavenumberY * M_PI;

	*result = sigma * sin( ky * y ) * cos( kx * x );
}

/* IC from Mirko Velic. This is the IC temperature for his solB, from his Analytic Suite. Added 22-May-2006 */
void Stg_FEM_VelicTemperatureIC_SolB( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*  context            = (DomainContext*)_context;
	FeVariable*             temperatureField   = (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	FeMesh*			feMesh               = temperatureField->feMesh;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double*                 coord;
	double                  x; 
	double                  y;
	double                  km; /*  for y-direction */
	double                  kn; /*  for x-direction */
	double                  wavenumberX;
	double                  wavenumberY;
	double                  L;
	double                  sigma;
	double			min[3], max[3];
	
	/* Find coordinate of node */
	coord = Mesh_GetVertex( feMesh, node_lI );
	Mesh_GetGlobalCoordRange( feMesh, min, max );

	/* Make sure that the box has right dimensions */
	assert( (max[ J_AXIS ] - min[ J_AXIS ] - 1.0 ) < SMALL );
	L = max[ I_AXIS ] - min[ I_AXIS ];

	x = coord[ I_AXIS ] - min[ I_AXIS ];
	y = coord[ J_AXIS ] - min[ J_AXIS ];

	wavenumberX = Dictionary_GetInt_WithDefault( dictionary, "wavenumberX", 1 );
	wavenumberY = Dictionary_GetDouble_WithDefault( dictionary, "wavenumberY", 2.0 );
	assert( wavenumberX != wavenumberY );
	sigma = Dictionary_GetDouble_WithDefault( dictionary, "sigma", 1.0 );

	kn = wavenumberX * M_PI / L;
/* 	 TODO: Re-write Mirko's code and/or Documentation so the input parameters for these ICs are less confusing */
	km = wavenumberY / L;

	*result = sigma * sinh( km * y ) * cos( kn * x );
}


/* Initial Condition derived from Boundary Layer theory -
   taken from P. E. van Keken, S. D. King, U. R. Schmeling, U. R. Christensen, D. Neumeister, and M.-P. Doin. A comparison of methods for the modeling of thermochemical convection. Journal of Geophysical Research, 102(B10):22477-22496, october 1997. */
void StgFEM_StandardConditionFunctions_AnalyticalTemperatureIC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*			feMesh               = NULL;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double*                 coord;
	double                  u0, v0, Q;
	double                  x, y;
	double                  RaT;
	double                  lambda, height, width;
	double                  Tu, Tl, Tr, Ts;
	double			min[3], max[3];

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	feMesh       = feVariable->feMesh;

	coord      = Mesh_GetVertex( feMesh, node_lI );
	Mesh_GetGlobalCoordRange( feMesh, min, max );

	/* Get Aspect Ratio */
	height = max[ J_AXIS ] - min[ J_AXIS ];
	width  = max[ I_AXIS ] - min[ I_AXIS ];
	lambda = width/height;
	
	x = coord[ I_AXIS ] - min[ I_AXIS ];
	y = coord[ J_AXIS ] - min[ J_AXIS ];
	
	/* Get thermal Rayleigh Number from Dictionary */
	RaT = Dictionary_GetDouble( dictionary, "RaT" );
	
	/* Horizontal fluid velocity at upper boundary & lower boundary - Equation A3 */
	u0 = pow( lambda , 7.0/3.0 )/ pow(1 + lambda*lambda*lambda*lambda, 2.0/3.0) * pow(0.5*RaT/sqrt(M_PI) , 2.0/3.0);

	/* Vertical velocity of the upwelling and downwellings - Modified from Van Keken to match Turcotte and Shubert */
	v0 = u0; //lambda;
	
	/* Total rate of heat flow out of the top of the cell per unit distance along the axis of the roll - Equation A3 */
	Q = 2.0 * sqrt(M_1_PI * lambda/u0);
	Tu = 0.5 * erf( 0.5 * ( 1 - y ) * sqrt(u0/x) );                                                      /* Equation A2a */
	Tl = 1.0 - 0.5 * erf(0.5 * y * sqrt(u0/(lambda-x)));                                                 /* Equation A2b */
	Tr = 0.5 + 0.5*Q/sqrt(M_PI) * sqrt(v0/(y+1)) * exp( -x*x*v0/(4*y+4) );                               /* Equation A2c */
	Ts = 0.5 - 0.5*Q/sqrt(M_PI) * sqrt(v0/(2-y)) * exp( -(lambda - x) * (lambda - x) * v0 / (8 - 4*y) ); /* Equation A2d */

	/* Equation A1 */
	*result = Tu + Tl + Tr + Ts - 1.5;

	/* Crop result */
	if ( *result > 1.0 ) 
		*result = 1.0;
	else if ( *result < 0.0 ) 
		*result = 0.0;
	
}

void StgFEM_StandardConditionFunctions_EdgeDriveConvectionIC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) 
{        
	DomainContext*  context = (DomainContext*)_context;        
	Dictionary*             dictionary         = context->dictionary;        
	FeVariable*             feVariable = NULL;        
	FeMesh*			mesh = NULL;        
	double*                 result = (double*) _result;        
	double                  perturbationAmplitude;        
	double                  thermalAnomalyOffset;        
	double*                 coord;        
	
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );        
	mesh       = feVariable->feMesh;        
	perturbationAmplitude = Dictionary_GetDouble_WithDefault( dictionary, "SinusoidalTempIC_PerturbationAmplitude", 0.1 );        
	thermalAnomalyOffset = Dictionary_GetDouble_WithDefault( dictionary, "thermalAnomalyOffset", 0.0 );        
	coord = Mesh_GetVertex( mesh, node_lI );
	
	/* eqn 1 from S.D.King & D.L. Anderson, "Edge-drive convection", EPSL 160 (1998) 289-296 */        
	
	*result = 1.0 + perturbationAmplitude * sin( M_PI * coord[ J_AXIS ] ) * cos( 0.5 * M_PI * ( coord[ I_AXIS ] + thermalAnomalyOffset ) );
}

void StgFEM_StandardConditionFunctions_SinusoidalExtension( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double                  frequency;
	double                  vel0;
	double                  amplitude;
	double                  phaseShift;

	frequency  = Dictionary_GetDouble_WithDefault( dictionary, "SinusoidalExtensionFrequency", 1.0 );
	vel0       = Dictionary_GetDouble_WithDefault( dictionary, "SinusoidalExtensionVelocity", 0.0 );
	amplitude  = Dictionary_GetDouble_WithDefault( dictionary, "SinusoidalExtensionAmplitude", 0.0 );
	phaseShift = Dictionary_GetDouble_WithDefault( dictionary, "SinusoidalExtensionPhaseShift", 0.0 );


	*result = vel0 + amplitude * cos( 2.0 * M_PI * frequency * (context->currentTime + context->dt - phaseShift ) );
}


void StgFEM_StandardConditionFunctions_StepFunction( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*			feMesh               = NULL;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
        double*                 coord;
	double                  offset;
	/* double                  value; */
	unsigned		dim;
	/* Bool			less; */
	double			left;
	double			right;
	double			grad;
	Bool			useGrad;
	double			leftBegin, rightEnd;

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	feMesh       = feVariable->feMesh;
	coord      = Mesh_GetVertex( feMesh, node_lI );

	offset = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionOffset", 0.0 );
	/* value = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionValue", 0.0 ); */
	dim = Dictionary_GetUnsignedInt_WithDefault( dictionary, "StepFunctionDim", 0 );
	/* less = Dictionary_GetBool_WithDefault( dictionary, "StepFunctionLessThan", True ); */
	left = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionLeftSide", 0.0 );
	right = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionRightSide", 0.0 );
	leftBegin = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionLeftBegin", 0.0 );
	rightEnd = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionRightEnd", 0.0 );
	grad = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionGradient", 0.0 );
	useGrad = Dictionary_GetBool_WithDefault( dictionary, "StepFunctionUseGradient", False );

	/*if( less ) {
		if( coord[dim] < offset ) {
			if( useGrad )
				*result = (offset - coord[dim])*grad;
			else
				*result = value;
		}
		else
			*result = 0;
	}
	else {
		if( coord[dim] > offset ) {
			if( useGrad )
				*result = (coord[dim] - offset)*grad;
			else
				*result = value;
		}
		else
			*result = 0;
	}*/
	/*
	 * have changed the step function such that it takes two different input values, one for each side
	 * of the offset - dave, 24.05.07
	 */
	if( coord[dim] < leftBegin || coord[dim] > rightEnd ) {
		*result = 0.0;
	}
	else if( coord[dim] < offset ) {
		if( useGrad )
			*result = (offset - coord[dim])*grad;
		else
			*result = left;
	}
	else if( coord[dim] >= offset ) {
		if( useGrad )
			*result = (coord[dim] - offset)*grad;
		else
			*result = right;
	}

/*         if(coord[0] < 5.0) */
/*           { */
/*             *result=0; */
/*           } */
/*         else if(coord[0] < 15.0) */
/*           { */
/*             *result=value*((coord[0]-5.0)/10); */
/*           } */
/*         else */
/*           { */
/*             *result=value; */
/*           } */
}

void StgFEM_StandardConditionFunctions_ConvectionBenchmark( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	/* This IC is for the 2D ConvectionBenchmark defined in
	 * http://www.mcc.monash.edu.au/twiki/view/Research/ConvectionBenchmarks
	 */
	
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*			mesh;
	double*                 result             = (double*) _result;
	double			min[3], max[3];
        double*                 coord;
	double                  x,y;
	double                  Lx, Ly;

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	mesh       = (FeMesh*)feVariable->feMesh;

	Mesh_GetGlobalCoordRange( mesh, min, max );
	
	Lx = max[ I_AXIS ] - min[ I_AXIS ];
	Ly = max[ J_AXIS ] - min[ J_AXIS ];
	
	coord      = Mesh_GetVertex( mesh, node_lI );

	x = ( coord[0] - min[ I_AXIS ] ) / Lx;
	y = ( coord[1] - min[ J_AXIS ] ) / Ly;


	*result = ( 1 - y ) + ( cos( M_PI * x ) * sin( M_PI * y ) ) / 100 ;
}

void StgFEM_StandardConditionFunctions_ConstantVelocity( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	FeVariable*             feVariable         = NULL;
	FeMesh*			mesh               = NULL;
	double*                 result             = (double*) _result;
	double                  velocity[3];
	
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;

	velocity[ I_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "ConstantVelocity_Vx", 0.0 );
	velocity[ J_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "ConstantVelocity_Vy", 1.0 );
	velocity[ K_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "ConstantVelocity_Vz", 0.0 );

	result[ I_AXIS ] = velocity[ I_AXIS ];
	result[ J_AXIS ] = velocity[ J_AXIS ];
	if( feVariable->dim == 3 )
		result[ K_AXIS ] = velocity[ K_AXIS ];
}


/* 3D spec ridge top BC (for milestone 1 of magma project) 
 * to be applied to the top x-z plane of the domain */
void StgFEM_StandardConditionFunctions_SpecRidge3D( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*			feMesh             = NULL;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
        double*                 coord;

	double			leftVal;
	double			rightVal;
	double			xOffset1;
	double			xOffset2;
	double			yOffset;


	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	feMesh     = feVariable->feMesh;
	coord      = Mesh_GetVertex( feMesh, node_lI );

	leftVal = Dictionary_GetDouble_WithDefault( dictionary, "SpecRidge3DLeftSide", 0.0 );
	rightVal = Dictionary_GetDouble_WithDefault( dictionary, "SpecRidge3DRightSide", 0.0 );
	xOffset1 = Dictionary_GetDouble_WithDefault( dictionary, "SpecRidge3DXOffset1", 0.0 );
	xOffset2 = Dictionary_GetDouble_WithDefault( dictionary, "SpecRidge3DXOffset2", 0.0 );
	yOffset = Dictionary_GetDouble_WithDefault( dictionary, "SpecRidge3DZOffset", 0.0 );

	if( coord[0] < xOffset1 )
		*result = leftVal;
	else if( coord[0] < xOffset2 && coord[2] > yOffset )
		*result = leftVal;
	else
		*result = rightVal;
}

/* get the BC's from the analytic solution as stored on the relevant FeVariable */
void StgFEM_StandardConditionFunctions_SpectralBCX( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             analyticFeVarX     = NULL;
	FeVariable*             numericFeVar       = NULL;
	FeMesh*			feMesh             = NULL;
	double*                 result             = (double*) _result;
        double*                 coord;
	Node_LocalIndex		analyticNodeI;
	Element_DomainIndex	analyticElement_I;
	double			analyticLocalElementCoord[3];
	FeMesh*			analyticFeMesh;
	
	analyticFeVarX = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "SpectralVelocityXField" );
	numericFeVar   = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	//feMesh         = numericFeVar->feMesh;
	//coord          = Mesh_GetVertex( feMesh, node_lI );

	//analyticFeMesh = analyticFeVarX->feMesh;
	//if( Mesh_SearchElements( analyticFeMesh, coord, &analyticElement_I ) ) {
	//	FeMesh_CoordGlobalToLocal( analyticFeMesh, analyticElement_I, coord, analyticLocalElementCoord );
	//	FeVariable_InterpolateWithinElement( analyticFeVarX, analyticElement_I, analyticLocalElementCoord, result );
	//}
	//else {	/* numerical solution node outside analytic mesh - just find closest point & use that */
	//	analyticNodeI  = Mesh_NearestVertex( analyticFeMesh, coord );
	//	FeVariable_GetValueAtNode( analyticFeVarX, analyticNodeI, result );
	//}

	/* same mesh used for spectral and numeric velocities - just get the value at the node */
	FeVariable_GetValueAtNode( analyticFeVarX, node_lI, result );

}

void StgFEM_StandardConditionFunctions_SpectralBCY( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             analyticFeVarY     = NULL;
	FeVariable*             numericFeVar       = NULL;
	FeMesh*			feMesh             = NULL;
	double*                 result             = (double*) _result;
        double*                 coord;
	Node_LocalIndex		analyticNodeI;
	Element_DomainIndex	analyticElement_I;
	double			analyticLocalElementCoord[3];
	FeMesh*			analyticFeMesh;
	
	analyticFeVarY = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "SpectralVelocityYField" );
	numericFeVar   = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	//feMesh         = numericFeVar->feMesh;
	//coord          = Mesh_GetVertex( feMesh, node_lI );

	//analyticFeMesh = analyticFeVarY->feMesh;
	//if( Mesh_SearchElements( analyticFeMesh, coord, &analyticElement_I ) ) {
	//	FeMesh_CoordGlobalToLocal( analyticFeMesh, analyticElement_I, coord, analyticLocalElementCoord );
	//	FeVariable_InterpolateWithinElement( analyticFeVarY, analyticElement_I, analyticLocalElementCoord, result );
	//}
	//else {
	//	analyticNodeI  = Mesh_NearestVertex( analyticFeMesh, coord );
	//	FeVariable_GetValueAtNode( analyticFeVarY, analyticNodeI, result );
	//}

	/* same mesh used for spectral and numeric velocities - just get the value at the node */
	FeVariable_GetValueAtNode( analyticFeVarY, node_lI, result );
}

void StgFEM_StandardConditionFunctions_SpectralBCZ( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             analyticFeVarZ     = NULL;
	FeVariable*             numericFeVar       = NULL;
	FeMesh*			feMesh             = NULL;
	double*                 result             = (double*) _result;
        double*                 coord;
	Node_LocalIndex		analyticNodeI;
	Element_DomainIndex	analyticElement_I;
	double			analyticLocalElementCoord[3];
	FeMesh*			analyticFeMesh;

	analyticFeVarZ = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "SpectralVelocityZField" );
	numericFeVar   = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	//feMesh         = numericFeVar->feMesh;
	//coord          = Mesh_GetVertex( feMesh, node_lI );

	//analyticFeMesh = analyticFeVarZ->feMesh;
	//if( Mesh_SearchElements( analyticFeMesh, coord, &analyticElement_I ) ) {
	//	FeMesh_CoordGlobalToLocal( analyticFeMesh, analyticElement_I, coord, analyticLocalElementCoord );
	//	FeVariable_InterpolateWithinElement( analyticFeVarZ, analyticElement_I, analyticLocalElementCoord, result );
	//}
	//else {
	//	analyticNodeI  = Mesh_NearestVertex( analyticFeMesh, coord );
	//	FeVariable_GetValueAtNode( analyticFeVarZ, analyticNodeI, result );
	//}

	/* same mesh used for spectral and numeric velocities - just get the value at the node */
	FeVariable_GetValueAtNode( analyticFeVarZ, node_lI, result );
}

void StgFEM_StandardConditionFunctions_SpectralPressureBCX( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             analyticFeVarX     = NULL;
	FeVariable*             numericFeVar       = NULL;
	FeMesh*			feMesh             = NULL;
	double*                 result             = (double*) _result;
        double*                 coord;

	analyticFeVarX = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "SpectralPressureField" );
	numericFeVar   = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "PressureField" );
	feMesh         = numericFeVar->feMesh;
	coord          = Mesh_GetVertex( feMesh, node_lI );

	FeVariable_GetValueAtNode( analyticFeVarX, node_lI, result );
}

void StgFEM_StandardConditionFunctions_SpectralPressureBCY( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             analyticFeVarY     = NULL;
	FeVariable*             numericFeVar       = NULL;
	FeMesh*			feMesh             = NULL;
	double*                 result             = (double*) _result;
        double*                 coord;

	analyticFeVarY = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "SpectralPressureField" );
	numericFeVar   = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "PressureField" );
	feMesh         = numericFeVar->feMesh;
	coord          = Mesh_GetVertex( feMesh, node_lI );

	FeVariable_GetValueAtNode( analyticFeVarY, node_lI, result );
}

/* error function for use in 3D spec ridge top BC */
double errorFunction( double z, int n ) {
	double		pi	= 3.1415926535;
	double 		a;
	double 		erf 	= 0.0;
	int		denom;
	int		i, j;

	a = 2.0/sqrt( pi );

	for( i=0 ; i<n ; i++ ) {
		denom = 1;
		for( j=1 ; j<=2*i+1 ; j+=2 ) 
			denom *= j; 
		
		erf += pow( 2, i )*pow( z, 2*i+1 )/denom;
	}

	return erf *= a*exp( -1.0*z*z );
}



void StgFEM_StandardConditionFunctions_ErrorFunc( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*			feMesh             = NULL;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
        double*                 coord;
	double			dilate;
	double			width; 

	feVariable  = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	feMesh      = feVariable->feMesh;
	coord       = Mesh_GetVertex( feMesh, node_lI );

	dilate      = Dictionary_GetDouble_WithDefault( dictionary, "ErrorFuncDilate", 0.0 );
	width       = Dictionary_GetDouble_WithDefault( dictionary, "ErrorFuncWidth", 0.0 );

	if( coord[0] < -1.0*width ) {
		*result = -1.0;
	}
	else if( coord[0] > width ) {
		*result = 1.0;
	}
	else {
		*result     = errorFunction( coord[0]/dilate, 5 );	
	}

}


