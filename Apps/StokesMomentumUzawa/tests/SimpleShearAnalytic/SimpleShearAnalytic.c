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
** $Id: SimpleShearAnalytic.c 964 2007-10-11 08:03:06Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

const Type SimpleShearAnalytic_Type = "SimpleShearAnalytic";

typedef struct { 
	__AnalyticSolution 
	double centreY;
	double factor;
} SimpleShearAnalytic;


void SimpleShearAnalytic_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ) {
	SimpleShearAnalytic *self = (SimpleShearAnalytic*)analyticSolution;
	
	velocity[ I_AXIS ] = self->factor * (coord[ J_AXIS ] - self->centreY);
	velocity[ J_AXIS ] = 0.0;
}


void SimpleShearAnalytic_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure ) {
	*pressure = 0.0;
}

void _SimpleShearAnalytic_Build( void* analyticSolution, void* data ) {
	SimpleShearAnalytic *self = (SimpleShearAnalytic*)analyticSolution;

	_AnalyticSolution_Build( self, data );
}

void _SimpleShearAnalytic_AssignFromXML( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	SimpleShearAnalytic *self = (SimpleShearAnalytic*)analyticSolution;
	FeVariable*       velocityField;
	FeVariable*       pressureField;

	_AnalyticSolution_AssignFromXML( self, cf, data );

	velocityField = Stg_ComponentFactory_ConstructByName( cf, (Name)"VelocityField", FeVariable, True, data  ); 
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, velocityField, SimpleShearAnalytic_VelocityFunction );
	
	pressureField = Stg_ComponentFactory_ConstructByName( cf, (Name)"PressureField", FeVariable, True, data  ); 
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, pressureField, SimpleShearAnalytic_PressureFunction );
	
	/* Set constants */
	self->centreY = Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"simpleShearCentreY", 0.0  );
	self->factor  = Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"SimpleShearFactor", 1.0 );
}

void* _SimpleShearAnalytic_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(SimpleShearAnalytic);
	Type                                                      type = SimpleShearAnalytic_Type;
	Stg_Class_DeleteFunction*                              _delete = _AnalyticSolution_Delete;
	Stg_Class_PrintFunction*                                _print = _AnalyticSolution_Print;
	Stg_Class_CopyFunction*                                  _copy = _AnalyticSolution_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _SimpleShearAnalytic_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _SimpleShearAnalytic_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _SimpleShearAnalytic_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _AnalyticSolution_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _AnalyticSolution_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _AnalyticSolution_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _AnalyticSolution_New(  ANALYTICSOLUTION_PASSARGS  );
}

/* This function is automatically run by StGermain when this plugin is loaded. The name must be "<plugin-name>_Register". */
Index StgFEM_SimpleShearAnalytic_Register( PluginsManager* pluginsManager  ) {
	/* A plugin is only properly registered once it returns the handle provided when submitting a codelet to StGermain. */
	return PluginsManager_Submit( pluginsManager, SimpleShearAnalytic_Type, (Name)"0", _SimpleShearAnalytic_DefaultNew  );
}


