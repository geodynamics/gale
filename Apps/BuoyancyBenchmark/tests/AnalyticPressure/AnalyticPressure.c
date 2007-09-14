/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
** Copyright (c) 2005-2006, Monash Cluster Computing, Building 28, Monash University Clayton Campus,
**	Victoria, 3800, Australia
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**
** Contributors:
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Patrick D. Sunter, Software Engineer, VPAC. (patrick@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	David Stegman, Postdoctoral Fellow, Monash University. (david.stegman@sci.monash.edu.au)
**	Wendy Sharples, PhD Student, Monash University (wendy.sharples@sci.monash.edu.au)
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
** $Id: AnalyticPressure.c 509 2007-09-14 06:00:17Z RobertTurnbull $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include <assert.h>

const Type AnalyticPressure_Type = "AnalyticPressure";

typedef struct {
	__AnalyticSolution
	FeVariable*		pressureField;
	double                  density;
	double                  gravity;
	double                  maxY;
	double                  minY;
} AnalyticPressure;

void _AnalyticPressure_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure ) {
	AnalyticPressure*  self    = (AnalyticPressure*)analyticSolution;
	double             density = self->density;
	double             gravity = self->gravity;
        double             y;

	/* Find coordinate of node */
	y = -coord[ J_AXIS ] + (self->minY+self->maxY)/2.0;

	*pressure = y * density * gravity;
}

	
void _AnalyticPressure_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	AnalyticPressure*        self           = (AnalyticPressure*)analyticSolution;

	/* Construct Parent */
	_AnalyticSolution_Construct( self, cf, data );

	/* Create Analytic Fields */
	self->pressureField = Stg_ComponentFactory_ConstructByName( cf, "PressureField", FeVariable, True, data ); 
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, self->pressureField, _AnalyticPressure_PressureFunction );
	
	
	self->density  = Stg_ComponentFactory_GetDouble( cf, "layer", "density", 0.0 );
	self->gravity  = Stg_ComponentFactory_GetRootDictDouble( cf, "gravity", 0 );
	self->maxY     = Stg_ComponentFactory_GetRootDictDouble( cf, "maxY", 0 );
	self->minY     = Stg_ComponentFactory_GetRootDictDouble( cf, "minY", 0 );
}

void _AnalyticPressure_Build( void* analyticSolution, void* data ) {
	AnalyticPressure*	self = (AnalyticPressure*)analyticSolution;

	assert( self && Stg_CheckType( self, AnalyticPressure ) );

	_AnalyticSolution_Build( self, data );
}

void* _AnalyticPressure_DefaultNew( Name name ) {
	return _AnalyticSolution_New(
			sizeof(AnalyticPressure),
			AnalyticPressure_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_AnalyticPressure_DefaultNew,
			_AnalyticPressure_Construct,
			_AnalyticPressure_Build,
			_AnalyticSolution_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );
}

Index _PICellerator_AnalyticPressure_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, AnalyticPressure_Type, "0", _AnalyticPressure_DefaultNew );
}
