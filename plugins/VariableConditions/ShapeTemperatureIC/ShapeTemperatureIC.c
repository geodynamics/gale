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
** $Id: ShapeTemperatureIC.c 610 2007-10-11 08:09:29Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include <assert.h>

const Type Underworld_ShapeTemperatureIC_Type = "Underworld_ShapeTemperatureIC";
typedef struct {
	__Codelet
} Underworld_ShapeTemperatureIC;

void Underworld_ShapeTemperatureICFunction( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	UnderworldContext*      context            = (UnderworldContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	FeMesh*			mesh               = NULL;
	double*                 result             = (double*) _result;
	Stg_Shape*              shape;
	Name                    shapeName;
	double*                 coord;
	
	mesh       = context->temperatureField->feMesh;

	shapeName = Dictionary_GetString( dictionary, "temperatureICShape" );
	shape = (Stg_Shape*) LiveComponentRegister_Get( context->CF->LCRegister, shapeName );
	assert( shape );

	/* Find coordinate of node */
	coord = Mesh_GetVertex( mesh, node_lI );

	if ( Stg_Shape_IsCoordInside( shape, coord ) ) 
		*result = 1.0;
	else 
		*result = 0.0;
}

void _Underworld_ShapeTemperatureIC_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	ConditionFunction*      condFunc;
	UnderworldContext*      context;

	context = (UnderworldContext*)Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data ); 
	
	condFunc = ConditionFunction_New( Underworld_ShapeTemperatureICFunction, "ShapeTemperatureIC" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );
}

void* _Underworld_ShapeTemperatureIC_DefaultNew( Name name ) {
	return Codelet_New(
		Underworld_ShapeTemperatureIC_Type,
		_Underworld_ShapeTemperatureIC_DefaultNew,
		_Underworld_ShapeTemperatureIC_Construct,
		_Codelet_Build,
		_Codelet_Initialise,
		_Codelet_Execute,
		_Codelet_Destroy,
		name );
}

Index Underworld_ShapeTemperatureIC_Register( PluginsManager* pluginsManager ) {
	Journal_DPrintf( Underworld_Debug, "In: %s( void* )\n", __func__ );

	return PluginsManager_Submit( pluginsManager, Underworld_ShapeTemperatureIC_Type, "0", _Underworld_ShapeTemperatureIC_DefaultNew );
}

