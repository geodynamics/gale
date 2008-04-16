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
** $Id: CPUTime.c 1107 2008-04-16 01:54:15Z BelindaMay $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/FrequentOutput/FrequentOutput.h>

#include "CPUTime.h"
  
const Type StgFEM_CPUTime_Type = "StgFEM_CPUTime";

void StgFEM_CPUTime_PrintTimeInfo( AbstractContext* context ) {
	StgFEM_CPUTime* self = (StgFEM_CPUTime*)LiveComponentRegister_Get( context->CF->LCRegister, StgFEM_CPUTime_Type );

	/* Print Current Time Taken */
	StgFEM_FrequentOutput_PrintValue( context, MPI_Wtime() - self->initialTime );
}

void _StgFEM_CPUTime_Construct( void* componment, Stg_ComponentFactory* cf, void* data ) {
	StgFEM_CPUTime* self = (StgFEM_CPUTime*)componment;
	AbstractContext* context;

	context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data ); 
	
	/* Initialise Timer */
	self->initialTime = MPI_Wtime();

	/* Print Header to file */
	StgFEM_FrequentOutput_PrintString( context, "CPU_Time" );
	
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput ,StgFEM_CPUTime_PrintTimeInfo );
}

void* _StgFEM_CPUTime_DefaultNew( Name name ) {
	return _Codelet_New(
			sizeof( StgFEM_CPUTime ),
			StgFEM_CPUTime_Type,
			_Codelet_Delete,
			_Codelet_Print,
			_Codelet_Copy,
			_StgFEM_CPUTime_DefaultNew,
			_StgFEM_CPUTime_Construct,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}
   
Index StgFEM_CPUTime_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, StgFEM_CPUTime_Type, "0", _StgFEM_CPUTime_DefaultNew );
}


