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
** $Id: Timestep.c 964 2007-10-11 08:03:06Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include "mpi.h"
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/LinearAlgebra/LinearAlgebra.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"
#include "types.h"
#include "Timestep.h"

#include "AdvectionDiffusionSLE.h"
#include "Residual.h"
#include <math.h>
#include <assert.h>
#include <string.h>


double AdvectionDiffusionSLE_CalculateDt( void* advectionDiffusionSLE, FiniteElementContext* context ) {
	AdvectionDiffusionSLE*     self   = (AdvectionDiffusionSLE*) advectionDiffusionSLE;
	double                     advectionTimestep;
	double                     diffusionTimestep;
	double                     advectionTimestepGlobal;
	double                     diffusionTimestepGlobal;
	double                     timestep;
	
	Journal_DPrintf( self->debug, "In func: %s\n", __func__ );
	
	/* Calculate Courant Number */
	advectionTimestep = AdvectionDiffusionSLE_AdvectiveTimestep( self );
	diffusionTimestep = AdvectionDiffusionSLE_DiffusiveTimestep( self );

	MPI_Allreduce( &advectionTimestep, &advectionTimestepGlobal, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
	MPI_Allreduce( &diffusionTimestep, &diffusionTimestepGlobal, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );

	Journal_Printf( self->debug, "%s Dominating. - Advective Timestep = %g - Diffusive Timestep = %g\n", 
			advectionTimestepGlobal < diffusionTimestepGlobal ? "Advection" : "Diffusion",
			advectionTimestepGlobal, diffusionTimestepGlobal);
	
	/* Calculate Time Step */
	timestep = MIN( advectionTimestepGlobal, diffusionTimestepGlobal );

	return timestep;
}


double AdvectionDiffusionSLE_DiffusiveTimestep( void* advectionDiffusionSLE ) {
	AdvectionDiffusionSLE*    self              = (AdvectionDiffusionSLE*) advectionDiffusionSLE;
	double                    minSeparation;
	double                    minSeparationEachDim[3];

	Journal_DPrintf( self->debug, "In func: %s\n", __func__ );
	
	FeVariable_GetMinimumSeparation( self->phiField, &minSeparation, minSeparationEachDim );

	return self->courantFactor * minSeparation * minSeparation / self->maxDiffusivity;
}


double AdvectionDiffusionSLE_AdvectiveTimestep( void* advectionDiffusionSLE ) {
	AdvectionDiffusionSLE*    self              = (AdvectionDiffusionSLE*) advectionDiffusionSLE;
	AdvDiffResidualForceTerm* residualForceTerm = self->advDiffResidualForceTerm;
	FeVariable*               velocityField     = residualForceTerm->velocityField;
	Node_LocalIndex           nodeLocalCount    = FeMesh_GetNodeLocalSize( self->phiField->feMesh );
	Node_LocalIndex           node_I;
	Dimension_Index           dim               = self->dim;
	Dimension_Index           dim_I;
	double                    timestep          = HUGE_VAL;
	XYZ                       velocity;
	double                    minSeparation;
	double                    minSeparationEachDim[3];
	double*                   meshCoord;
	
	Journal_DPrintf( self->debug, "In func: %s\n", __func__ );

	FeVariable_GetMinimumSeparation( self->phiField, &minSeparation, minSeparationEachDim );

	for( node_I = 0 ; node_I < nodeLocalCount ; node_I++ ){
		meshCoord = Mesh_GetVertex( self->phiField->feMesh, node_I );
		FieldVariable_InterpolateValueAt( velocityField, meshCoord, velocity );
		
		for ( dim_I = 0 ; dim_I < dim ; dim_I++ ) {
			if( velocity[ dim_I ] == 0.0 ) 
				continue;
			timestep = MIN( timestep, fabs( minSeparationEachDim[ dim_I ]/velocity[ dim_I ] ) );
		}
	}

	return self->courantFactor * timestep;
}
