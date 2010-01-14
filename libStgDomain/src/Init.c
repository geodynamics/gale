/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
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
** $Id: Init.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>

#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>
#include <StgDomain/Utils/Utils.h>
#include <StgDomain/Swarm/Swarm.h>

#include "Init.h"

#include <stdio.h>

Bool StgDomain_Init( int* argc, char** argv[] ) {
	/* This init function tells StGermain of all the component types, etc this module contributes. Because it can be linked at compile
	   time or linked in by a toolbox at runtime, we need to make sure it isn't run twice (compiled in and loaded through a toolbox.*/
	if( !ToolboxesManager_IsInitialised( stgToolboxesManager, "StgDomain" ) ) {
		int tmp;
	
		Journal_Printf( Journal_Register( DebugStream_Type, (Name)"Context"  ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */
		tmp = Stream_GetPrintingRank( Journal_Register( InfoStream_Type, (Name)"Context" )  );
		Stream_SetPrintingRank( Journal_Register( InfoStream_Type, (Name)"Context"  ), 0 );
		Journal_Printf( /* DO NOT CHANGE OR REMOVE */
			Journal_Register( InfoStream_Type, (Name)"Context"  ), 
			"StGermain Domain Library revision %s. Copyright (C) 2003-2007 VPAC.\n", VERSION );
		Stream_Flush( Journal_Register( InfoStream_Type, (Name)"Context" )  );
		Stream_SetPrintingRank( Journal_Register( InfoStream_Type, (Name)"Context"  ), tmp );
	
		StgDomainGeometry_Init( argc, argv );
		StgDomainShape_Init( argc, argv );
		StgDomainMesh_Init( argc, argv );
		StgDomainUtils_Init( argc, argv );
		StgDomainSwarm_Init( argc, argv );

		ToolboxesManager_SetInitialised( stgToolboxesManager, "StgDomain" );
		return True;
	}
	return False;
}


