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
** $Id: IncompressibleExtensionAnalytic.c 650 2008-02-19 01:27:56Z RobertTurnbull $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include "../../IncompressibleExtensionBC.h"

const Type IncompressibleExtensionAnalytic_Type = "IncompressibleExtensionAnalytic";

typedef struct { 
	__AnalyticSolution 
	FeVariable* velocityField;
	FeVariable* strainRateInvField;
	double constantHeight;
} IncompressibleExtensionAnalytic;


void IncompressibleExtensionAnalytic_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ) {
	IncompressibleExtensionAnalytic *self = (IncompressibleExtensionAnalytic*)analyticSolution;

	GetVelocity( self->velocityField, self->constantHeight, coord, velocity );
}

void IncompressibleExtensionAnalytic_StrainRateInvariantFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* inv ) {
	IncompressibleExtensionAnalytic *self = (IncompressibleExtensionAnalytic*)analyticSolution;
	FeVariable*         velocityField = self->velocityField;
	double              V_b = GetLeftWallVelocity( velocityField );
	double              V_a = GetRightWallVelocity( velocityField );
	double              width;
	XYZ                 min, max;
	
	/* Calculate Width and Height */
	FieldVariable_GetMinAndMaxGlobalCoords( velocityField, min, max );
	width  = (max[ I_AXIS ] - min[ I_AXIS ]);

	*inv = (V_a - V_b)/width;
}

void _IncompressibleExtensionAnalytic_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	IncompressibleExtensionAnalytic *self = (IncompressibleExtensionAnalytic*)analyticSolution;
	UnderworldContext*  context  = Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data );

	_AnalyticSolution_Construct( self, cf, data );

	self->velocityField = Stg_ComponentFactory_ConstructByName( cf, "VelocityField", FeVariable, True, data ); 
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, self->velocityField, IncompressibleExtensionAnalytic_VelocityFunction );
	
	self->strainRateInvField = Stg_ComponentFactory_ConstructByName( cf, "StrainRateInvariantField", FeVariable, True, data ); 
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, self->strainRateInvField, IncompressibleExtensionAnalytic_StrainRateInvariantFunction );

	/* Set constants */
	self->constantHeight = Stg_ComponentFactory_GetRootDictDouble( cf, "constantHeight", 0.0 );

	/* Rearrange EP */
	EP_Remove( Context_GetEntryPoint( context, AbstractContext_EP_UpdateClass ), self->type );
	EntryPoint_PrependClassHook( Context_GetEntryPoint( context, AbstractContext_EP_UpdateClass ), 
			self->type, AnalyticSolution_Update, self->type, self );
//	EntryPoint_PrependClassHook( Context_GetEntryPoint( context, AbstractContext_EP_UpdateClass ), 
//			self->type, Stg_Component_Initialise, self->type, self );
}

void _IncompressibleExtensionAnalytic_Build( void* analyticSolution, void* data ) {
	IncompressibleExtensionAnalytic *self = (IncompressibleExtensionAnalytic*)analyticSolution;

	_AnalyticSolution_Build( self, data );
}

void* _IncompressibleExtensionAnalytic_DefaultNew( Name name ) {
	return (void*) _AnalyticSolution_New( 
			sizeof(IncompressibleExtensionAnalytic),
			IncompressibleExtensionAnalytic_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_IncompressibleExtensionAnalytic_DefaultNew,
			_IncompressibleExtensionAnalytic_Construct,
			_IncompressibleExtensionAnalytic_Build, 
			_AnalyticSolution_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );
}

/* This function is automatically run by StGermain when this plugin is loaded. The name must be "<plugin-name>_Register". */
Index Underworld_IncompressibleExtensionAnalytic_Register( PluginsManager* pluginsManager ) {
	/* A plugin is only properly registered once it returns the handle provided when submitting a codelet to StGermain. */
	return PluginsManager_Submit( pluginsManager, IncompressibleExtensionAnalytic_Type, "0", _IncompressibleExtensionAnalytic_DefaultNew );
}
