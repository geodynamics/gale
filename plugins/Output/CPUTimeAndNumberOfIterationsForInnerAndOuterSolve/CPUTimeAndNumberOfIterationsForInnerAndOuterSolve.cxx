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
** $Id: CPUTimeAndNumberOfIterationsForInnerAndOuterSolve.c 1107 2008-04-16 01:54:15Z BelindaMay $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/FrequentOutput/FrequentOutput.h>
#include "CPUTimeAndNumberOfIterationsForInnerAndOuterSolve.h"
#include "StgFEM/Discretisation/Discretisation.h"
#include "../../../SLE/SystemSetup/src/SystemSetup.h"
#include "../../../SLE/ProvidedSystems/StokesFlow/src/StokesFlow.h"
//#include "../../../SLE/ProvidedSystems/StokesFlow/src/Stokes_SLE_UzawaSolve.h"

  
const Type StgFEM_CPUTimeAndNumberOfIterationsForInnerAndOuterSolve_Type = "StgFEM_CPUTimeAndNumberOfIterationsForInnerAndOuterSolve";

void StgFEM_CPUTimeAndNumberOfIterationsForInnerAndOuterSolve_PrintTimeInfo( AbstractContext* context ) {
	Stokes_SLE*	sle  = (Stokes_SLE*) LiveComponentRegister_Get( context->CF->LCRegister, (Name)"stokesEqn");
	SLE_Solver*	solver = (SLE_Solver* ) LiveComponentRegister_Get( context->CF->LCRegister, (Name)"uzawa");

	/* Print Current Average InnerIteration Time Taken */
	StgFEM_FrequentOutput_PrintValue( context, solver->avgtimeinnerits);
	/* Print Current Average OuterIteration Time Taken */
	StgFEM_FrequentOutput_PrintValue( context, solver->avgtimeouterits);
	/* Print Current Average NonLinearIteration Time Taken */
	if(sle->isNonLinear == True){
		StgFEM_FrequentOutput_PrintValue( context, solver->avgtimenonlinearits);
	}	
	/* Print Average Number of Inner Iterations*/
	StgFEM_FrequentOutput_PrintValue( context, solver->avgnuminnerits);
	/* Print Number of Outer Iterations */
	StgFEM_FrequentOutput_PrintValue( context, solver->avgnumouterits);
	/* Print Number of NonLinear Iterations*/
	if(sle->isNonLinear == True){
		StgFEM_FrequentOutput_PrintValue( context, solver->totalnumnonlinearits);
	}
}

void _StgFEM_CPUTimeAndNumberOfIterationsForInnerAndOuterSolve_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
	StgFEM_CPUTimeAndNumberOfIterationsForInnerAndOuterSolve* self = (StgFEM_CPUTimeAndNumberOfIterationsForInnerAndOuterSolve*)component;
	AbstractContext* context;
	context = Stg_ComponentFactory_ConstructByName( cf, (Name)"context", AbstractContext, True, data  ); 
	self->context = context;
		
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput ,StgFEM_CPUTimeAndNumberOfIterationsForInnerAndOuterSolve_PrintTimeInfo );
}

void _StgFEM_CPUTimeAndNumberOfIterationsForInnerAndOuterSolve_Initialise( void* component, void* data ) {
	StgFEM_CPUTimeAndNumberOfIterationsForInnerAndOuterSolve* self = (StgFEM_CPUTimeAndNumberOfIterationsForInnerAndOuterSolve*)component;
	Stokes_SLE*  sle  = (Stokes_SLE*) LiveComponentRegister_Get( self->context->CF->LCRegister, (Name)"stokesEqn" );
		
	/*this isn't set to true before the initialise phase*/
	
	/* Print Header to file */
	StgFEM_FrequentOutput_PrintString( self->context, "AvgCPUInner" );
	StgFEM_FrequentOutput_PrintString( self->context, "AvgCPUOuter" );
	if(sle->isNonLinear == True){
		StgFEM_FrequentOutput_PrintString( self->context, "AvgCPUNonLin" );
	}
	StgFEM_FrequentOutput_PrintString( self->context, "AvgInIts" );
	StgFEM_FrequentOutput_PrintString( self->context, "AvgOutIts" );
	if(sle->isNonLinear == True){
		StgFEM_FrequentOutput_PrintString( self->context, "NonLinIts" );
	}
}

void* _StgFEM_CPUTimeAndNumberOfIterationsForInnerAndOuterSolve_DefaultNew( Name name ) {
	SizeT                                              _sizeOfSelf = sizeof( StgFEM_CPUTimeAndNumberOfIterationsForInnerAndOuterSolve );
	Type                                                      type = StgFEM_CPUTimeAndNumberOfIterationsForInnerAndOuterSolve_Type;
	Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
	Stg_Class_PrintFunction*                                _print = _Codelet_Print;
	Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _StgFEM_CPUTimeAndNumberOfIterationsForInnerAndOuterSolve_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _StgFEM_CPUTimeAndNumberOfIterationsForInnerAndOuterSolve_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Codelet_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _StgFEM_CPUTimeAndNumberOfIterationsForInnerAndOuterSolve_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _Codelet_Destroy;
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return _Codelet_New( CODELET_PASSARGS );
}
   
Index StgFEM_CPUTimeAndNumberOfIterationsForInnerAndOuterSolve_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, StgFEM_CPUTimeAndNumberOfIterationsForInnerAndOuterSolve_Type, (Name)"0", _StgFEM_CPUTimeAndNumberOfIterationsForInnerAndOuterSolve_DefaultNew  );
}


