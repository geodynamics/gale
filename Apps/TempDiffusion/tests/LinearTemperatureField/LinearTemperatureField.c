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
** $Id: LinearTemperatureField.c 830 2007-05-10 06:45:16Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <string.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>

const Type LinearTemperatureField_Type = "LinearTemperatureField";

typedef struct { __AnalyticSolution } LinearTemperatureField;

void LinearTemperatureField_TemperatureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* value ) {
	*value = 1.0 - coord[ J_AXIS ];
}

void _LinearTemperatureField_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	LinearTemperatureField *self = (LinearTemperatureField*)analyticSolution;
	FeVariable* temperatureField;

	_AnalyticSolution_Construct( self, cf, data );

	temperatureField = Stg_ComponentFactory_ConstructByName( cf, "TemperatureField", FeVariable, True, data ); 
	Stg_ObjectList_Append( self->feVariableList, temperatureField );
}

void _LinearTemperatureField_Build( void* analyticSolution, void* data ) {
	LinearTemperatureField *self = (LinearTemperatureField*)analyticSolution;
	FeVariable* feVariable       = (FeVariable*) Stg_ObjectList_At( self->feVariableList, 0 );
	char* feName                 = feVariable->name;

	Stg_Component_Build( feVariable, data, False );
	if( !strcmp( feName, "TemperatureField" ) )
		AnalyticSolution_CreateAnalyticVectorField( self, feVariable, LinearTemperatureField_TemperatureFunction );
	else
		assert(0);

	_AnalyticSolution_Build( self, data );
}

void* _LinearTemperatureField_DefaultNew( Name name ) {
	return (void*) _AnalyticSolution_New( 
			sizeof(LinearTemperatureField),
			LinearTemperatureField_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_LinearTemperatureField_DefaultNew,
			_LinearTemperatureField_Construct,
			_LinearTemperatureField_Build,
			_AnalyticSolution_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );
}

/* This function is automatically run by StGermain when this plugin is loaded. The name must be "<plugin-name>_Register". */
Index StgFEM_LinearTemperatureField_Register( PluginsManager* pluginsManager ) {
	/* A plugin is only properly registered once it returns the handle provided when submitting a codelet to StGermain. */
	return PluginsManager_Submit( pluginsManager, LinearTemperatureField_Type, "0", _LinearTemperatureField_DefaultNew );
}
