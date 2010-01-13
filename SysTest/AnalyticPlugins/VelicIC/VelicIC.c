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
** $Id: VelicIC.c 610 2007-10-11 08:09:29Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <assert.h>
#include "VelicIC.h"

const Type Underworld_VelicIC_Type = "Underworld_VelicIC";

#define SMALL 1.0e-5

/* Works with SolA */
void Underworld_VelicIC_Sinusoidal( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*  context            = (DomainContext*)_context;
	FeVariable*             temperatureField   = (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	FeMesh*			feMesh             = temperatureField->feMesh;
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

	wavenumberX = Dictionary_GetInt_WithDefault( dictionary, (Dictionary_Entry_Key)"wavenumberX", 1  );
	wavenumberY = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"wavenumberY", 1.0  );
	sigma = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"sigma", 1.0 );
	
	assert( sigma > 0.0 );
	assert( wavenumberY > 0.0 );
	assert( wavenumberX > 0.0 );
	
	kx = (double)wavenumberX * M_PI / Lx;
	ky = (double)wavenumberY * M_PI;

	*result = sigma * sin( ky * y ) * cos( kx * x  );
}

/* Works with SolB */
void Underworld_VelicIC_Hyperbolic( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*  context            = (DomainContext*)_context;
	FeVariable*             temperatureField   = (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	FeMesh*			feMesh             = temperatureField->feMesh;
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

	wavenumberX = Dictionary_GetInt_WithDefault( dictionary, (Dictionary_Entry_Key)"wavenumberX", 1  );
	wavenumberY = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"wavenumberY", 2.0 );
	assert( wavenumberX != wavenumberY  );
	sigma = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"sigma", 1.0 );

	kn = wavenumberX * M_PI / L;
/* 	 TODO: Re-write Mirko's code and/or Documentation so the input parameters for these ICs are less confusing */
	km = wavenumberY / L;

	*result = sigma * sinh( km * y ) * cos( kn * x  );
}



void _Underworld_VelicIC_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
	AbstractContext*        context;
	ConditionFunction*      condFunc;

	context = (AbstractContext*)Stg_ComponentFactory_ConstructByName( cf, (Name)"context", AbstractContext, True, data  ); 
	
	condFunc = ConditionFunction_New( Underworld_VelicIC_Sinusoidal, (Name)"VelicIC_Sinusoidal" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( Underworld_VelicIC_Hyperbolic, (Name)"VelicIC_Hyperbolic" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
}	

void* _Underworld_VelicIC_DefaultNew( Name name ) {
	return Codelet_New(
		Underworld_VelicIC_Type,
		_Underworld_VelicIC_DefaultNew,
		_Underworld_VelicIC_AssignFromXML,
		_Codelet_Build,
		_Codelet_Initialise,
		_Codelet_Execute,
		_Codelet_Destroy,
		name );
}

Index Underworld_VelicIC_Register( PluginsManager* pluginsManager ) {
	Journal_DPrintf( StgFEM_Debug, "In: %s( void* )\n", __func__ );

	return PluginsManager_Submit( pluginsManager, Underworld_VelicIC_Type, (Name)"0", _Underworld_VelicIC_DefaultNew  );
}


