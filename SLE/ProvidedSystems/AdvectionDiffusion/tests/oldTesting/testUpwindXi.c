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
** $Id: testUpwindXi.c 964 2007-10-11 08:03:06Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"
#include "StgFEM/SLE/ProvidedSystems/AdvectionDiffusion/AdvectionDiffusion.h"

int main( int argc, char* argv[] ) {
	Stream* dataStream;
	double  minPercletNumber = -15.0;
	double  maxPercletNumber = 15.0;
	double  dPerceltNumber = 0.5;
	double  perceltNumber;
	char*   dataFileName = "./output/output.dat";

	MPI_Init( &argc, &argv );
	StGermain_Init( &argc, &argv );

	/* Create Data File */
	dataStream    = Journal_Register( Info_Type, "DataStream" );
	Stream_RedirectFile( dataStream, dataFileName );

	Journal_Printf( dataStream, "# File to compare code with Brooks, Hughes 1982 - Fig 3.3\n");
	Journal_Printf( dataStream, "# Integration rule for the optimal upwind scheme, doubly asymptotic approximation and critical approximation.\n");
	Journal_Printf( dataStream, "# Plot using line:\n");
	Journal_Printf( dataStream, "# plot \"%s\" using 1:2 title \"Exact\" with line, ", dataFileName );
	Journal_Printf( dataStream, "\"%s\" using 1:3 title \"DoublyAsymptoticAssumption\" with line, ", dataFileName );
	Journal_Printf( dataStream, "\"%s\" using 1:4 title \"CriticalAssumption\" with line\n", dataFileName );

	Journal_Printf( dataStream, "# Perclet Number \t Exact \t DoublyAsymptoticAssumption \t CriticalAssumption\n" );
	for ( perceltNumber = minPercletNumber ; perceltNumber < maxPercletNumber ; perceltNumber += dPerceltNumber )
		Journal_Printf( dataStream, "%0.3g \t\t %0.3g \t\t %0.3g \t\t %0.3g\n", 
				perceltNumber, 
				AdvDiffResidualForceTerm_UpwindXiExact( NULL, perceltNumber),
				AdvDiffResidualForceTerm_UpwindXiDoublyAsymptoticAssumption( NULL, perceltNumber),
				AdvDiffResidualForceTerm_UpwindXiCriticalAssumption( NULL, perceltNumber) );

	StGermain_Finalise();
	MPI_Finalize();

	return EXIT_SUCCESS;
}
