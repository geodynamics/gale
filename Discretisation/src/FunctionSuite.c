/* Suite of functions for evaluating diagnostic properties from StG.
 * components (ie: Swarms, FeVariables, OperatorFeVariables)
 *
 * */

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
