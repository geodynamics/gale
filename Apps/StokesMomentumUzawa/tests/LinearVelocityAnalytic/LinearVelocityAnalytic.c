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
** $Id: LinearVelocityAnalytic.c 1099 2008-04-08 05:51:09Z RobertTurnbull $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <string.h>

const Type LinearVelocityAnalytic_Type = "LinearVelocityAnalytic";

typedef struct { 
	__AnalyticSolution
	FeVariable* velocityField;
	double nodeVelocity[8][3];
	int    cornerNodeCount;
} LinearVelocityAnalytic;

Index Grid_ProjectIJK( Grid* grid, Index i, Index j, Index k ) {
	IJK ijk = {0,0,0};
	
	ijk[0] = i;
	ijk[1] = j;
	ijk[2] = k;

	return Grid_Project( grid, ijk );
}
Index Grid_ProjectIJK_MinMax( Grid* grid, Bool iIsMax, Bool jIsMax, Bool kIsMax ) {
	IJK ijk = {0,0,0};
	
	if ( iIsMax )
		ijk[0] = grid->sizes[0] - 1;
	if ( jIsMax )
		ijk[1] = grid->sizes[1] - 1;
	if ( kIsMax )
		ijk[2] = grid->sizes[2] - 1;

	return Grid_Project( grid, ijk );
}

void LinearVelocityAnalytic_GetCornerNodeVelocities(void* analyticSolution) {
	LinearVelocityAnalytic *self = (LinearVelocityAnalytic*)analyticSolution;
	Grid*                   vertGrid;
	Node_GlobalIndex        nodeMapper[8];
	FeVariable*             velocityField = self->velocityField;
	FeMesh*                 mesh = velocityField->feMesh;
	Dimension_Index         dim = velocityField->dim;
	Node_Index              globalNode_I;
	Node_Index              ii;
	
	vertGrid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

	/* Find global indicies of nodes */
	self->cornerNodeCount = 4;
	nodeMapper[0] = Grid_ProjectIJK_MinMax( vertGrid, 0, 0, 0 );
	nodeMapper[1] = Grid_ProjectIJK_MinMax( vertGrid, 1, 0, 0 );
	nodeMapper[2] = Grid_ProjectIJK_MinMax( vertGrid, 0, 1, 0 );
	nodeMapper[3] = Grid_ProjectIJK_MinMax( vertGrid, 1, 1, 0 );
	if ( dim == 3 ) {
		self->cornerNodeCount = 8;
		nodeMapper[4] = Grid_ProjectIJK_MinMax( vertGrid, 0, 0, 1 );
		nodeMapper[5] = Grid_ProjectIJK_MinMax( vertGrid, 1, 0, 1 );
		nodeMapper[6] = Grid_ProjectIJK_MinMax( vertGrid, 0, 1, 1 );
		nodeMapper[7] = Grid_ProjectIJK_MinMax( vertGrid, 1, 1, 1 );
	}

	/* Loop over corner nodes */
	for ( ii = 0 ; ii < self->cornerNodeCount ; ii++ ) {
		globalNode_I = nodeMapper[ ii ];
		FeVariable_GetValueAtNodeGlobal( velocityField, globalNode_I, self->nodeVelocity[ii] );
	}
}

/* Do a normal linear interpolation as if the box were a FEM element */
void LinearVelocityAnalytic_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ) {
	LinearVelocityAnalytic *self = (LinearVelocityAnalytic*)analyticSolution;
	FeVariable*             velocityField = self->velocityField;
	FeMesh*                 mesh = velocityField->feMesh;
	Dimension_Index         dim = velocityField->dim;
	Dimension_Index         dim_I;
	XYZ                     xi;
	XYZ                     min;
	XYZ                     max;
	double                  Ni[8];
	Node_Index              ii;
	ElementType*            elementType;

	/* Transform the coordinate into a master coordinate system */
	_FeVariable_GetMinAndMaxGlobalCoords( velocityField, min, max );
	for ( dim_I = 0 ; dim_I < dim ; dim_I++ ) {
		xi[ dim_I ] = 2.0 * (coord[ dim_I ] - min[ dim_I ])/(max[dim_I] - min[dim_I]) - 1;
	}	

	/* Get Shape Functions */
	elementType = FeMesh_GetElementType( mesh, 0 );
	ElementType_EvaluateShapeFunctionsAt( elementType, xi, Ni );

	/* Do interpolation */
	/* Loop over corner nodes */
	memset( velocity, 0, dim*sizeof(double) );
	for ( ii = 0 ; ii < self->cornerNodeCount ; ii++ ) {
		for ( dim_I = 0 ; dim_I < dim ; dim_I++ ) {
			velocity[ dim_I ] += Ni[ ii ] * self->nodeVelocity[ii][ dim_I ];
		}
	}
}
void LinearVelocityAnalytic_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure ) {
	*pressure = 0.0;
}

void _LinearVelocityAnalytic_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	LinearVelocityAnalytic *self = (LinearVelocityAnalytic*)analyticSolution;
	FeVariable*       pressureField;

	_AnalyticSolution_Construct( self, cf, data );

	self->velocityField = Stg_ComponentFactory_ConstructByName( cf, "VelocityField", FeVariable, True, data ); 
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, self->velocityField, LinearVelocityAnalytic_VelocityFunction );
	
	pressureField = Stg_ComponentFactory_ConstructByName( cf, "PressureField", FeVariable, True, data ); 
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, pressureField, LinearVelocityAnalytic_PressureFunction );
	
}

void _LinearVelocityAnalytic_Initialise( void* analyticSolution, void* data ) {
	LinearVelocityAnalytic *self = (LinearVelocityAnalytic*)analyticSolution;
	
	Stg_Component_Initialise( self->velocityField, data,  False );
	LinearVelocityAnalytic_GetCornerNodeVelocities( self );

	_AnalyticSolution_Initialise( self, data );
}

void* _LinearVelocityAnalytic_DefaultNew( Name name ) {
	return (void*) _AnalyticSolution_New( 
			sizeof(LinearVelocityAnalytic),
			LinearVelocityAnalytic_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_LinearVelocityAnalytic_DefaultNew,
			_LinearVelocityAnalytic_Construct,
			_AnalyticSolution_Build,
			_LinearVelocityAnalytic_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );
}

/* This function is automatically run by StGermain when this plugin is loaded. The name must be "<plugin-name>_Register". */
Index StgFEM_LinearVelocityAnalytic_Register( PluginsManager* pluginsManager ) {
	/* A plugin is only properly registered once it returns the handle provided when submitting a codelet to StGermain. */
	return PluginsManager_Submit( pluginsManager, LinearVelocityAnalytic_Type, "0", _LinearVelocityAnalytic_DefaultNew );
}
