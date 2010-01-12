/* Suite of functions for evaluating diagnostic properties from StG.
 * components (ie: Swarms, FeVariables, OperatorFeVariables)
 *
 * */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "types.h"
#include "FeMesh.h"
#include "FeVariable.h"
#include "FunctionSuite.h"

/* root mean squared velocity, takes in a velocity squared
 * OperatorFeVariable and an integration Swarm */
double StgFEM_Vrms( FeVariable* velsq, Swarm* swarm ) {
	Mesh*		mesh		= (Mesh*) velsq->feMesh;
	double		max[3], min[3];
	double		volume;
	double		integral;
	double		vrms;
	unsigned	dim		= Mesh_GetDimSize( mesh );

	Mesh_GetGlobalCoordRange( mesh, min, max );

	integral = FeVariable_Integrate( velsq, swarm );

	volume = ( max[I_AXIS] - min[I_AXIS] ) * ( max[J_AXIS] - min[J_AXIS] );
	if( dim == 3 )
		volume *= ( max[K_AXIS] - min[K_AXIS] );

	vrms = sqrt( integral / volume );

	return vrms;
}

/* optimised version of the feVariable interpolation function. assumes constant number of dofs for each node */
void StgFEM_InterpolateValue_WithNi( void* _feVariable, Element_LocalIndex lElement_I, double* Ni, double* value ) {
	FeVariable*             self        = (FeVariable*) _feVariable;
	Node_ElementLocalIndex  elLocalNode_I;
	Node_LocalIndex         lNode_I;
	Dof_Index               dof_I;
	Dof_Index               dofCount;
	Variable*               dofVariable;
	double                  nodeValue;
	unsigned		nInc, *inc;

	/* Gets number of degrees of freedom - assuming it is the same throughout the mesh */
	dofCount = self->dofLayout->dofCounts[0];

	/* Initialise */
	memset( value, 0, sizeof( double )*dofCount );

	FeMesh_GetElementNodes( self->feMesh, lElement_I, self->inc );
	nInc = IArray_GetSize( self->inc );
	inc = IArray_GetPtr( self->inc );

	for ( dof_I = 0; dof_I < dofCount; dof_I++ ) {
		dofVariable = DofLayout_GetVariable( self->dofLayout, 0, dof_I );
		for ( elLocalNode_I = 0 ; elLocalNode_I < nInc ; elLocalNode_I++) {
			lNode_I      = inc[ elLocalNode_I ];
			nodeValue    = Variable_GetValueDouble( dofVariable, lNode_I );
			value[dof_I] += Ni[elLocalNode_I] * nodeValue;
		}
	}
}

/* optimised version of the feVariable interpolate derivatives function. assumes constant number of dofs for each node */
void StgFEM_InterpolateDerivatives_WithGNx( void* _feVariable, Element_LocalIndex lElement_I, double** GNx, double* value ) {
	FeVariable*             self        = (FeVariable*) _feVariable;
	Node_ElementLocalIndex  elLocalNode_I;
	Node_LocalIndex         lNode_I;
	Dof_Index               dof_I, dim_I;
	Dof_Index               dofCount;
	Variable*               dofVariable;
	double                  nodeValue;
	unsigned		nInc, *inc;
	Dimension_Index         dim         = self->dim;

	/* Gets number of degrees of freedom - assuming it is the same throughout the mesh */
	dofCount = self->dofLayout->dofCounts[0];

	/* Initialise */
	memset( value, 0, sizeof( double )*dofCount*dim );

	FeMesh_GetElementNodes( self->feMesh, lElement_I, self->inc );
	nInc = IArray_GetSize( self->inc );
	inc = IArray_GetPtr( self->inc );

	for ( dof_I = 0 ; dof_I < dofCount ; dof_I++ ) {
		dofVariable  = DofLayout_GetVariable( self->dofLayout, 0, dof_I );
		/* Interpolate derivative from nodes */
		for ( elLocalNode_I = 0; elLocalNode_I < nInc; elLocalNode_I++) {
			lNode_I      = inc[ elLocalNode_I ];
			nodeValue    = Variable_GetValueDouble( dofVariable, lNode_I );
			
			for( dim_I = 0; dim_I < dim; dim_I++ )
				value[dof_I*dim + dim_I] += GNx[dim_I][elLocalNode_I] * nodeValue;
		}
	}
}



