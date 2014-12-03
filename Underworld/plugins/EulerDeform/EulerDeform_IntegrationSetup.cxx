/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	David May, PhD Student Monash University, VPAC. (david.may@sci.maths.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
** $Id: LevelSetPlg.c 200 2005-07-08 08:24:41Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include "types.h"
#include "Context.h"
#include "EulerDeform.h"

void EulerDeform_IntegrationSetup( void* _timeIntegrator, void* context ) {
	EulerDeform_Context*	edCtx = (EulerDeform_Context*)context;
	unsigned					sys_i;

	FeVariable_SyncShadowValues( edCtx->systems[0].velField );

	/* 
	** We'll need to store side values that we require to be static here, for later
	** return to the mesh.
	*/

	for( sys_i = 0; sys_i < edCtx->nSystems; sys_i++ ) {
		EulerDeform_System* sys = edCtx->systems + sys_i;

		if( sys->staticSides ) {
			IndexSet	*tmpIndSet;
			unsigned	nInds, *inds;
			unsigned	nDims;
			unsigned	ind_i;

			/* Collect indices of all the sides. */

                        tmpIndSet = EulerDeform_CreateStaticSet(sys);
                        IndexSet_GetMembers( tmpIndSet, &nInds, &inds );

			/* Copy coords to temporary array. */
			nDims = Mesh_GetDimSize( sys->mesh );
			sys->sideCoords = AllocArray2D( double, nInds, nDims );
			for( ind_i = 0; ind_i < nInds; ind_i++ )
				memcpy( sys->sideCoords[ind_i], sys->mesh->verts[inds[ind_i]], nDims * sizeof(double) );
                        FreeObject( tmpIndSet );
                        FreeArray( inds );
		}
	}

	/* Update advection arrays. */
	for( sys_i = 0; sys_i < edCtx->nSystems; sys_i++ ) {
		EulerDeform_System*	sys = edCtx->systems + sys_i;
		unsigned					nDims;
		unsigned					nLocalNodes;
		unsigned					n_i;

		nDims = Mesh_GetDimSize( sys->mesh );
		nLocalNodes = Mesh_GetLocalSize( sys->mesh, MT_VERTEX );

		for( n_i = 0; n_i < nLocalNodes; n_i++ )
			memcpy( sys->verts + n_i * nDims, sys->mesh->verts[n_i], nDims * sizeof(double) );
	}
}


