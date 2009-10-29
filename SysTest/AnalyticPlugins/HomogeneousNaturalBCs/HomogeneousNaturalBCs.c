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
** $Id: HomogeneousNaturalBCs.c 967 2007-10-23 05:27:09Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

const Type HomogeneousNaturalBCs_Type = "HomogeneousNaturalBCs";

typedef struct { 
	__FieldTest
	double angle;
	FeVariable* temperatureField;
} HomogeneousNaturalBCs;


void HomogeneousNaturalBCs_Velocity_SkewToMesh( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context = (DomainContext*)_context;
	HomogeneousNaturalBCs*  self    = Stg_ComponentFactory_ConstructByName( 
		context->CF, 
		HomogeneousNaturalBCs_Type, 
		HomogeneousNaturalBCs, 
		True,
		0 );
	double*                 result  = (double*) _result;
	
	result[ I_AXIS ] =  cos( self->angle );
	result[ J_AXIS ] =  sin( self->angle );
}


void HomogeneousNaturalBCs_TemperatureFunction( void* analyticSolution, double* coord, double* temperature ) {
	HomogeneousNaturalBCs *self = (HomogeneousNaturalBCs*)analyticSolution;

	if ( coord[ J_AXIS ] < tan( self->angle ) * coord[ I_AXIS ] + 0.25 )
		*temperature = 1.0;
	else 
		*temperature = 0.0;
}
	
void HomogeneousNaturalBCs_TemperatureBC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context    = (DomainContext*)_context;
	HomogeneousNaturalBCs*  self       = Stg_ComponentFactory_ConstructByName( 
		context->CF, 
		HomogeneousNaturalBCs_Type, 
		HomogeneousNaturalBCs, 
		True,
		0 );
	FeVariable*             feVariable = NULL;
	FeMesh*			mesh       = NULL;
	double*                 result     = (double*) _result;
	double*                 coord;
	
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	mesh       = feVariable->feMesh;
	coord = Mesh_GetVertex( mesh, node_lI );

	HomogeneousNaturalBCs_TemperatureFunction( self, coord, result );
}

void _HomogeneousNaturalBCs_AssignFromXML( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	HomogeneousNaturalBCs* self = (HomogeneousNaturalBCs*)analyticSolution;
	AbstractContext*       context;
	ConditionFunction*     condFunc;

	_FieldTest_AssignFromXML( self, cf, data );

	self->angle = StGermain_DegreeToRadian (Stg_ComponentFactory_GetRootDictDouble( cf, "VelocitySkewAngle", 45.0 ) );

	/* Create Condition Functions */
/*
	context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data ); 
*/
	condFunc = ConditionFunction_New( HomogeneousNaturalBCs_Velocity_SkewToMesh, "Velocity_SkewToMesh" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( HomogeneousNaturalBCs_TemperatureBC, "Temperature_StepFunction" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
}

void _HomogeneousNaturalBCs_Build( void* analyticSolution, void* data ) {
	HomogeneousNaturalBCs* self = (HomogeneousNaturalBCs*)analyticSolution;	

	_FieldTest_Build( self, data );

	/* here we assign the memory and the func ptr for analytic sols */
	self->_analyticSolutionList = Memory_Alloc_Array_Unnamed( FieldTest_AnalyticSolutionFunc*, 1 );
	/* this order MUST be consistent with the xml file definition */
	self->_analyticSolutionList[0] = HomogeneousNaturalBCs_TemperatureFunction;
}

void* _HomogeneousNaturalBCs_DefaultNew( Name name ) {
	return (void*) _FieldTest_New( 
			sizeof(HomogeneousNaturalBCs),
			HomogeneousNaturalBCs_Type,
			_FieldTest_Delete,
			_FieldTest_Print,
			_FieldTest_Copy,
			_HomogeneousNaturalBCs_DefaultNew,
			_HomogeneousNaturalBCs_AssignFromXML,
			_HomogeneousNaturalBCs_Build,
			_FieldTest_Initialise,
			_FieldTest_Execute,
			_FieldTest_Destroy,
			name );
}

/* This function is automatically run by StGermain when this plugin is loaded. The name must be "<plugin-name>_Register". */
Index StgFEM_HomogeneousNaturalBCs_Register( PluginsManager* pluginsManager ) {
	/* A plugin is only properly registered once it returns the handle provided when submitting a codelet to StGermain. */
	return PluginsManager_Submit( pluginsManager, HomogeneousNaturalBCs_Type, "0", _HomogeneousNaturalBCs_DefaultNew );
}

