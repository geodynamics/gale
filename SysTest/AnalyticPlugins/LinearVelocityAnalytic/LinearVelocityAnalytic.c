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
** $Id: LinearVelocityAnalytic.c 1111 2008-04-23 04:12:36Z RobertTurnbull $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <string.h>

const Type LinearVelocityAnalytic_Type = "LinearVelocityAnalytic";

typedef struct { 
	__FieldTest
	FeVariable* velocityField;
	double  nodeVelocity[8][3];
	double  nodeCoords[8][3];
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
	
	vertGrid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, ExtensionManager_GetHandle( mesh->info, (Name)"vertexGrid" )  );

	/* Find global indicies of nodes */
	self->cornerNodeCount = 4;
	nodeMapper[0] = Grid_ProjectIJK_MinMax( vertGrid, False, False, False );
	nodeMapper[1] = Grid_ProjectIJK_MinMax( vertGrid, True, False, False );
	nodeMapper[2] = Grid_ProjectIJK_MinMax( vertGrid, False, True, False );
	nodeMapper[3] = Grid_ProjectIJK_MinMax( vertGrid, True, True, False );
	if ( dim == 3 ) {
		self->cornerNodeCount = 8;
		nodeMapper[4] = Grid_ProjectIJK_MinMax( vertGrid, False, False, True );
		nodeMapper[5] = Grid_ProjectIJK_MinMax( vertGrid, True, False, True );
		nodeMapper[6] = Grid_ProjectIJK_MinMax( vertGrid, False, True, True );
		nodeMapper[7] = Grid_ProjectIJK_MinMax( vertGrid, True, True, True );
	}

	/* Loop over corner nodes */
	for ( ii = 0 ; ii < self->cornerNodeCount ; ii++ ) {
		globalNode_I = nodeMapper[ ii ];
		FeVariable_GetValueAtNodeGlobal( velocityField, globalNode_I, self->nodeVelocity[ii] );
		FeVariable_GetCoordAtNodeGlobal( velocityField, globalNode_I, self->nodeCoords[ii] );
	}
}

void GetLocalCoords( LinearVelocityAnalytic* self, double* coord, double* xi ) {
	FeVariable*             velocityField = self->velocityField;
	XYZ                     min;
	XYZ                     max;
	Dimension_Index         dim = velocityField->dim;
	Dimension_Index         dim_I;

	_FeVariable_GetMinAndMaxGlobalCoords( velocityField, min, max );
	for ( dim_I = 0 ; dim_I < dim ; dim_I++ ) {
		xi[ dim_I ] = 2.0 * (coord[ dim_I ] - min[ dim_I ])/(max[dim_I] - min[dim_I]) - 1;
	}	
}	

/* Do a normal linear interpolation as if the box were a FEM element */
void LinearVelocityAnalytic_VelocityFunction( void* analyticSolution, double* coord, double* velocity ) {
	LinearVelocityAnalytic *self = (LinearVelocityAnalytic*)analyticSolution;
	FeVariable*             velocityField = self->velocityField;
	FeMesh*                 mesh = velocityField->feMesh;
	Dimension_Index         dim = velocityField->dim;
	Dimension_Index         dim_I;
	XYZ                     xi;
	double                  Ni[8];
	Node_Index              ii;
	ElementType*            elementType;

	/* Transform the coordinate into a master coordinate system */
	GetLocalCoords( self, coord, xi );

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
void LinearVelocityAnalytic_PressureFunction( void* analyticSolution, double* coord, double* pressure ) {
	*pressure = 0.0;
}

void LinearVelocityAnalytic_VelocityGradientsFunction( void* analyticSolution, double* coord, double* velocityGradients ) {
	LinearVelocityAnalytic *self = (LinearVelocityAnalytic*)analyticSolution;
	FeVariable*             velocityField = self->velocityField;
	FeMesh*                 mesh = velocityField->feMesh;
	ElementType*            elementType;
	XYZ                     xi;
	double                  jac[3][3];
	double                  cof[3][3];	/* cofactors */
	double                  detJac;
	Node_Index              node_I;
	double**                GNi; 
	double**                GNx; 
	double*                 nodeCoord;
	double                  nodeValue;
	Dimension_Index         dim = velocityField->dim;
	Dimension_Index         i, j, dx, dxi;
	Dimension_Index         dim_I;
	
	/* Transform the coordinate into a master coordinate system */
	GetLocalCoords( self, coord, xi );

	GNi = Memory_Alloc_2DArray( double, dim, self->cornerNodeCount, (Name)"GNi"  );
	GNx = Memory_Alloc_2DArray( double, dim, self->cornerNodeCount, (Name)"GNx"  );

	/* Get Shape Functions */
	elementType = FeMesh_GetElementType( mesh, 0 );
	elementType->_evaluateShapeFunctionLocalDerivsAt( elementType, xi, GNi );

	/* build the jacobian matrix */
	/*
	jac = 	\sum_i d/d\xi( N_i ) x_i 		\sum_i d/d\xi( N_i ) y_i
			\sum_i d/d\eta( N_i ) x_i 		\sum_i d/d\eta( N_i ) y_i
	*/
	if( dim == 2 ) {
		jac[0][0] = jac[0][1] = jac[1][0] = jac[1][1] = 0.0;
		for(  node_I =0;  node_I < self->cornerNodeCount ;  node_I ++){	
			nodeCoord = self->nodeCoords[ node_I ];
			jac[0][0] = jac[0][0] + GNi[0][node_I] * nodeCoord[0];
			jac[0][1] = jac[0][1] + GNi[0][node_I] * nodeCoord[1];
			
			jac[1][0] = jac[1][0] + GNi[1][node_I] * nodeCoord[0];
			jac[1][1] = jac[1][1] + GNi[1][node_I] * nodeCoord[1];
		}
	}
	
	if( dim == 3 ) {
		jac[0][0] = jac[0][1] = jac[0][2] = 0.0;
		jac[1][0] = jac[1][1] = jac[1][2] = 0.0;
		jac[2][0] = jac[2][1] = jac[2][2] = 0.0;
		for(  node_I =0;  node_I < self->cornerNodeCount ;  node_I ++){	
			nodeCoord = self->nodeCoords[ node_I ];
			jac[0][0] = jac[0][0] + GNi[0][node_I] * nodeCoord[0];
			jac[0][1] = jac[0][1] + GNi[0][node_I] * nodeCoord[1];
			jac[0][2] = jac[0][2] + GNi[0][node_I] * nodeCoord[2];
			
			jac[1][0] = jac[1][0] + GNi[1][node_I] * nodeCoord[0];
			jac[1][1] = jac[1][1] + GNi[1][node_I] * nodeCoord[1];
			jac[1][2] = jac[1][2] + GNi[1][node_I] * nodeCoord[2];
			
			jac[2][0] = jac[2][0] + GNi[2][node_I] * nodeCoord[0];
			jac[2][1] = jac[2][1] + GNi[2][node_I] * nodeCoord[1];
			jac[2][2] = jac[2][2] + GNi[2][node_I] * nodeCoord[2];
		}
	}
	
	/* get determinant of the jacobian matrix */
	if( dim == 2 ) {
		detJac = jac[0][0]*jac[1][1] - jac[0][1]*jac[1][0]; 
	}		
	if( dim == 3 ) {
		detJac = jac[0][0]*( jac[1][1]*jac[2][2] - jac[1][2]*jac[2][1] ) 
				  - jac[0][1]*( jac[1][0]*jac[2][2] - jac[1][2]*jac[2][0] ) 
				  + jac[0][2]*( jac[1][0]*jac[2][1] - jac[1][1]*jac[2][0] );
	}
	
	/* invert the jacobian matrix A^-1 = adj(A)/det(A) */
	if( dim == 2 ) {
		double tmp = jac[0][0];
		jac[0][0] = jac[1][1]/detJac;
		jac[1][1] = tmp/detJac;
		jac[0][1] = -jac[0][1]/detJac;
		jac[1][0] = -jac[1][0]/detJac;		
	}
	if( dim == 3 ) {
		/*
		00 01 02
		10 11 12
		20 21 22		
		*/		
		cof[0][0] = jac[1][1]*jac[2][2] - jac[1][2]*jac[2][1];
		cof[1][0] = -(jac[1][0]*jac[2][2] - jac[1][2]*jac[2][0]);
		cof[2][0] = jac[1][0]*jac[2][1] - jac[1][1]*jac[2][0];
		
		cof[0][1] = -(jac[0][1]*jac[2][2] - jac[0][2]*jac[2][1]);
		cof[1][1] = jac[0][0]*jac[2][2] - jac[0][2]*jac[2][0];
		cof[2][1] = -(jac[0][0]*jac[2][1] - jac[0][1]*jac[2][0]);
		
		cof[0][2] = jac[0][1]*jac[1][2] - jac[0][2]*jac[1][1];
		cof[1][2] = -(jac[0][0]*jac[1][2] - jac[0][2]*jac[1][0]);
		cof[2][2] = jac[0][0]*jac[1][1] - jac[0][1]*jac[1][0];
		
		for( i=0; i<dim; i++ ) {
			for( j=0; j<dim; j++ ) {
				jac[i][j] = cof[i][j]/detJac;
			}
		}
		
		
	}
	
	/* get global derivs Ni_x, Ni_y and Ni_z if dim == 3 */
	for( dx=0; dx<dim; dx++ ) {
		for(  node_I =0;  node_I < self->cornerNodeCount ;  node_I ++){	
			double globalSF_DerivVal = 0.0;
			for(dxi=0; dxi<dim; dxi++) {
				globalSF_DerivVal = globalSF_DerivVal + GNi[dxi][node_I] * jac[dx][dxi];
			}
			
			GNx[dx][node_I] = globalSF_DerivVal;
		}
	}
	
	/* Initialise velocity gradients */
	memset( velocityGradients, 0, sizeof( double ) * dim * dim );
	
	for ( dim_I = 0 ; dim_I < dim ; dim_I++ ) {
		/* Interpolate derivative from nodes */
		for( node_I =0;  node_I < self->cornerNodeCount ;  node_I ++){	
			nodeValue    = self->nodeVelocity[ node_I ][ dim_I ];
			
			velocityGradients[dim_I*dim + 0] += GNx[0][node_I] * nodeValue;
			velocityGradients[dim_I*dim + 1] += GNx[1][node_I] * nodeValue;
			if( dim == 3 ) 
				velocityGradients[dim_I*dim + 2] += GNx[2][node_I] * nodeValue;	
		}
	}
	Memory_Free( GNi );
	Memory_Free( GNx );
}

void LinearVelocityAnalytic_StrainRateFunction( void* analyticSolution, double* coord, double* strainRate ) {
	LinearVelocityAnalytic *self = (LinearVelocityAnalytic*)analyticSolution;
	Dimension_Index         dim = self->velocityField->dim;
	TensorArray             velocityGradients;

	/* Get Velocity Gradients */
	LinearVelocityAnalytic_VelocityGradientsFunction( self, coord, velocityGradients );
	
	/* Get Strain Rate */
	TensorArray_GetSymmetricPart( velocityGradients, dim, strainRate );
}

void LinearVelocityAnalytic_StrainRateInvFunction( void* analyticSolution, double* coord, double* strainRateInv ) {
	LinearVelocityAnalytic *self = (LinearVelocityAnalytic*)analyticSolution;
	Dimension_Index         dim = self->velocityField->dim;
	SymmetricTensor         strainRate;

	/* Get Strain Rate */
	LinearVelocityAnalytic_StrainRateFunction( self, coord, strainRate );
	
	/* Get Invariant */
	*strainRateInv = SymmetricTensor_2ndInvariant( strainRate, dim );
}

void _LinearVelocityAnalytic_AssignFromXML( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	LinearVelocityAnalytic *self = (LinearVelocityAnalytic*)analyticSolution;

	_FieldTest_AssignFromXML( self, cf, data );

	self->velocityField = Stg_ComponentFactory_ConstructByName( cf, (Name)"VelocityField", FeVariable, True, data  ); 
}

void _LinearVelocityAnalytic_Build( void* analyticSolution, void* data ) {
	LinearVelocityAnalytic *self = (LinearVelocityAnalytic*)analyticSolution;
	
	_FieldTest_Build( self, data );

	/* here we assign the memory and the func ptr for analytic sols */
	self->_analyticSolutionList = Memory_Alloc_Array_Unnamed( FieldTest_AnalyticSolutionFunc*, 4 );
	/* this order MUST be consistent with the xml file definition */
	self->_analyticSolutionList[0] = LinearVelocityAnalytic_VelocityFunction;
	self->_analyticSolutionList[1] = LinearVelocityAnalytic_PressureFunction;
	self->_analyticSolutionList[2] = LinearVelocityAnalytic_StrainRateFunction;
	self->_analyticSolutionList[3] = LinearVelocityAnalytic_StrainRateInvFunction;
}

void _LinearVelocityAnalytic_Initialise( void* analyticSolution, void* data ) {
	LinearVelocityAnalytic *self = (LinearVelocityAnalytic*)analyticSolution;
	
	Stg_Component_Initialise( self->velocityField, data,  False );
	LinearVelocityAnalytic_GetCornerNodeVelocities( self );

	_FieldTest_Initialise( self, data );
}

void* _LinearVelocityAnalytic_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(LinearVelocityAnalytic);
	Type                                                      type = LinearVelocityAnalytic_Type;
	Stg_Class_DeleteFunction*                              _delete = _FieldTest_Delete;
	Stg_Class_PrintFunction*                                _print = _FieldTest_Print;
	Stg_Class_CopyFunction*                                  _copy = _FieldTest_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _LinearVelocityAnalytic_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _LinearVelocityAnalytic_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _LinearVelocityAnalytic_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _LinearVelocityAnalytic_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _FieldTest_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _FieldTest_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _FieldTest_New(  FIELDTEST_PASSARGS  );
}

/* This function is automatically run by StGermain when this plugin is loaded. The name must be "<plugin-name>_Register". */
Index StgFEM_LinearVelocityAnalytic_Register( PluginsManager* pluginsManager ) {
	/* A plugin is only properly registered once it returns the handle provided when submitting a codelet to StGermain. */
	return PluginsManager_Submit( pluginsManager, LinearVelocityAnalytic_Type, (Name)"0", _LinearVelocityAnalytic_DefaultNew  );
}


