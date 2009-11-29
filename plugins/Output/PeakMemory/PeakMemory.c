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
** $Id:  $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include <StgFEM/FrequentOutput/FrequentOutput.h>

#include "PeakMemory.h"
  
const Type StgFEM_PeakMemory_Type = "StgFEM_PeakMemory";

void StgFEM_PeakMemory_PrintMemoryInfo( AbstractContext* context ) {
	PetscLogDouble stgPeak, totalMem, petscMem, ave;
	
	stgPeak = (PetscLogDouble)stgMemory->stgPeakMemory;
	MPI_Allreduce( &stgPeak, &ave, 1, MPI_DOUBLE, MPI_SUM, context->communicator );
	ave /= 1024 * 1024;
	StgFEM_FrequentOutput_PrintValue( context, ave );

	PetscMallocGetMaximumUsage( &petscMem );
	MPI_Allreduce( &petscMem, &ave, 1, MPI_DOUBLE, MPI_SUM, context->communicator );
	ave /= 1024 * 1024;
	StgFEM_FrequentOutput_PrintValue( context, ave );

	PetscMemoryGetMaximumUsage( &totalMem );
	MPI_Allreduce( &totalMem, &ave, 1, MPI_DOUBLE, MPI_SUM, context->communicator );
	ave /= 1024 * 1024;
	StgFEM_FrequentOutput_PrintValue( context, ave );
}

void _StgFEM_PeakMemory_AssignFromXML( void* componment, Stg_ComponentFactory* cf, void* data ) {
	AbstractContext* context;

	/* Turn on the magical petsc logging */
	PetscMemorySetGetMaximumUsage();

	context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data ); 
	StgFEM_FrequentOutput* self = (StgFEM_FrequentOutput*)LiveComponentRegister_Get(
                  context->CF->LCRegister,
                  StgFEM_FrequentOutput_Type );

	/* set the stupid stream column width so I don't get "..." behaviour */
	self->columnWidth = 15;
	
	/* Print Header to file */
	StgFEM_FrequentOutput_PrintString( context, "StgPeakMem(Mb)" );
	StgFEM_FrequentOutput_PrintString( context, "PetscMem(Mb)" );
	StgFEM_FrequentOutput_PrintString( context, "ProgMem(Mb)" );
	
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput ,StgFEM_PeakMemory_PrintMemoryInfo );
}

void* _StgFEM_PeakMemory_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof( Codelet );
	Type                                                      type = StgFEM_PeakMemory_Type;
	Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
	Stg_Class_PrintFunction*                                _print = _Codelet_Print;
	Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _StgFEM_PeakMemory_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _StgFEM_PeakMemory_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Codelet_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _Codelet_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _Codelet_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	return _Codelet_New(  CODELET_PASSARGS  );
}
   
Index StgFEM_PeakMemory_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, StgFEM_PeakMemory_Type, "0", _StgFEM_PeakMemory_DefaultNew );
}



