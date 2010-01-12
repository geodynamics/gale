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
** $Id: HomogeneousEssentialBCs.c 968 2007-10-23 07:53:39Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

const Type HomogeneousEssentialBCs_Type = "HomogeneousEssentialBCs";

typedef struct { 
	__AnalyticSolution
	double angle;
	FeVariable* temperatureField;
} HomogeneousEssentialBCs;


void HomogeneousEssentialBCs_Velocity_SkewToMesh( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context = (DomainContext*)_context;
	HomogeneousEssentialBCs*  self    = Stg_ComponentFactory_ConstructByName( 
		context->CF, 
		HomogeneousEssentialBCs_Type, 
		HomogeneousEssentialBCs, 
		True,
		0 );
	double*                 result  = (double*) _result;
	
	result[ I_AXIS ] =  cos( self->angle );
	result[ J_AXIS ] =  sin( self->angle );
}


void HomogeneousEssentialBCs_TemperatureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* temperature ) {
	HomogeneousEssentialBCs *self = (HomogeneousEssentialBCs*)analyticSolution;

	if ( coord[ J_AXIS ] < tan( self->angle ) * coord[ I_AXIS ] + 0.25 )
		*temperature = 1.0;
	else 
		*temperature = 0.0;

	/* Perform Essential Boundary Conditions */
	if ( (coord[ I_AXIS ] > 0.99 && coord[ J_AXIS ] > 0.01) || coord[ J_AXIS ] > 0.99 ) {
		*temperature = 0.0;
	}
}
	
void HomogeneousEssentialBCs_TemperatureBC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context    = (DomainContext*)_context;
	HomogeneousEssentialBCs*  self       = Stg_ComponentFactory_ConstructByName( 
		context->CF, 
		HomogeneousEssentialBCs_Type, 
		HomogeneousEssentialBCs, 
		True,
		0 );
	FeVariable*             feVariable = NULL;
	FeMesh*			mesh       = NULL;
	double*                 result     = (double*) _result;
	double*                 coord;
	
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	mesh       = feVariable->feMesh;
	coord = Mesh_GetVertex( mesh, node_lI );

	HomogeneousEssentialBCs_TemperatureFunction( self, feVariable, coord, result );
}

void _HomogeneousEssentialBCs_AssignFromXML( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	HomogeneousEssentialBCs* self = (HomogeneousEssentialBCs*)analyticSolution;
	AbstractContext*       context;
	ConditionFunction*     condFunc;

	_AnalyticSolution_AssignFromXML( self, cf, data );

	self->temperatureField = Stg_ComponentFactory_ConstructByName( cf, "TemperatureField", FeVariable, True, data ); 
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, self->temperatureField, HomogeneousEssentialBCs_TemperatureFunction );

	self->angle = StGermain_DegreeToRadian (Stg_ComponentFactory_GetRootDictDouble( cf, "VelocitySkewAngle", 45.0 ) );

	/* Create Condition Functions */
	context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data ); 
	condFunc = ConditionFunction_New( HomogeneousEssentialBCs_Velocity_SkewToMesh, "Velocity_SkewToMesh" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( HomogeneousEssentialBCs_TemperatureBC, "Temperature_StepFunction" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );
}

void _HomogeneousEssentialBCs_Build( void* analyticSolution, void* data ) {
	HomogeneousEssentialBCs* self = (HomogeneousEssentialBCs*)analyticSolution;	

	_AnalyticSolution_Build( self, data );
}

void* _HomogeneousEssentialBCs_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(HomogeneousEssentialBCs);
	Type                                                      type = HomogeneousEssentialBCs_Type;
	Stg_Class_DeleteFunction*                              _delete = _AnalyticSolution_Delete;
	Stg_Class_PrintFunction*                                _print = _AnalyticSolution_Print;
	Stg_Class_CopyFunction*                                  _copy = _AnalyticSolution_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _HomogeneousEssentialBCs_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _HomogeneousEssentialBCs_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _HomogeneousEssentialBCs_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _AnalyticSolution_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _AnalyticSolution_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _AnalyticSolution_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _AnalyticSolution_New(  ANALYTICSOLUTION_PASSARGS  );
}

/* This function is automatically run by StGermain when this plugin is loaded. The name must be "<plugin-name>_Register". */
Index StgFEM_HomogeneousEssentialBCs_Register( PluginsManager* pluginsManager ) {
	/* A plugin is only properly registered once it returns the handle provided when submitting a codelet to StGermain. */
	return PluginsManager_Submit( pluginsManager, HomogeneousEssentialBCs_Type, "0", _HomogeneousEssentialBCs_DefaultNew );
}



