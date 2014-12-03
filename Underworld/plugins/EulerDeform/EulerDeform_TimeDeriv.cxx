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

Bool EulerDeform_TimeDeriv( void* crdAdvector, Index arrayInd, double* timeDeriv ) {
	TimeIntegrand*		self = (TimeIntegrand*)crdAdvector;
	FeVariable*				velocityField = (FeVariable*)self->data[0];
	InterpolationResult	result = LOCAL;

	/* check if the node information is on the local proc */
	if (arrayInd >= Mesh_GetDomainSize(velocityField->feMesh,MT_VERTEX) )
		result = OTHER_PROC;

	FeVariable_GetValueAtNode( velocityField, arrayInd, timeDeriv );

	/* Check if periodic */
	if ( Stg_Class_IsInstance( velocityField->feMesh->generator, CartesianGenerator_Type ) ) {
		CartesianGenerator* cartesianGenerator = (CartesianGenerator*) velocityField->feMesh->generator;
		if ( cartesianGenerator->periodic[ I_AXIS ] )
			timeDeriv[I_AXIS] = 0.0;
		if ( cartesianGenerator->periodic[ J_AXIS ] )
			timeDeriv[J_AXIS] = 0.0;
		if ( cartesianGenerator->periodic[ K_AXIS ] )
			timeDeriv[K_AXIS] = 0.0;
	}
		
	if ( result == OTHER_PROC || result == OUTSIDE_GLOBAL || isinf(timeDeriv[0]) || isinf(timeDeriv[1]) || 
	     ( velocityField->dim == 3 && isinf(timeDeriv[2]) ) ) 
	{
#if 0
		Journal_Printf( Journal_Register( Error_Type, (Name)self->type  ),
				"Error in func '%s' for particle with index %u.\n\tPosition (%g, %g, %g)\n\tVelocity here is (%g, %g, %g)."
				"\n\tInterpolation result is %s.\n",
				__func__, array_I, coord[0], coord[1], coord[2], 

				InterpolationResultToStringMap[result]  );
#endif	
		return False;	
	}

	return True;
}
