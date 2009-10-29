/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
** Copyright (c) 2005-2006, Monash Cluster Computing, Building 28, Monash University Clayton Campus,
**	Victoria, 3800, Australia
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**
** Contributors:
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Patrick D. Sunter, Software Engineer, VPAC. (patrick@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	David Stegman, Postdoctoral Fellow, Monash University. (david.stegman@sci.monash.edu.au)
**	Wendy Sharples, PhD Student, Monash University (wendy.sharples@sci.monash.edu.au)
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
** $Id: Context.c 518 2007-10-11 08:07:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include "types.h"
#include "Materials_Register.h"
#include "Material.h"
#include "PICelleratorContext.h"

/* Textual name of this class */
const Type PICelleratorContext_Type = "PICelleratorContext";

/* Constructors ------------------------------------------------------------------------------------------------*/
void* _PICelleratorContext_DefaultNew( Name name ) {
	return (void*) _PICelleratorContext_New(
		sizeof(PICelleratorContext),
		PICelleratorContext_Type,
		_PICelleratorContext_Delete,
		_PICelleratorContext_Print,
		NULL, 
		_PICelleratorContext_DefaultNew,
		_PICelleratorContext_AssignFromXML,
		_AbstractContext_Build,
		_AbstractContext_Initialise,
		_AbstractContext_Execute,
		_AbstractContext_Destroy,
		_PICelleratorContext_SetDt,
		name,
		False,
		0,
		0,
		MPI_COMM_WORLD,
		NULL );
}

PICelleratorContext* PICelleratorContext_New( 
		Name                                            name,
		double                                          start,
		double                                          stop,
		MPI_Comm                                        communicator,
		Dictionary*                                     dictionary )
{
	return _PICelleratorContext_New(
		sizeof(PICelleratorContext),
		PICelleratorContext_Type,
		_PICelleratorContext_Delete,
		_PICelleratorContext_Print,
		NULL, 
		_PICelleratorContext_DefaultNew,
		_PICelleratorContext_AssignFromXML,
		_AbstractContext_Build,
		_AbstractContext_Initialise,
		_AbstractContext_Execute,
		_AbstractContext_Destroy,
		_PICelleratorContext_SetDt,
		name,
		True,
		start,
		stop,
		communicator,
		dictionary );
}	


PICelleratorContext* _PICelleratorContext_New( 
		SizeT                                           sizeOfSelf,
		Type                                            type,
		Stg_Class_DeleteFunction*                       _delete,
		Stg_Class_PrintFunction*                        _print,
		Stg_Class_CopyFunction*                         _copy,
		Stg_Component_DefaultConstructorFunction*       _defaultConstructor,
		Stg_Component_ConstructFunction*                _construct,
		Stg_Component_BuildFunction*                    _build,
		Stg_Component_InitialiseFunction*               _initialise,
		Stg_Component_ExecuteFunction*                  _execute,
		Stg_Component_DestroyFunction*                  _destroy,
		AbstractContext_SetDt*                          _setDt,
		Name                                            name,
		Bool                                            initFlag,
		double                                          start,
		double                                          stop,
		MPI_Comm                                        communicator,
		Dictionary*                                     dictionary )
{
	PICelleratorContext* self;
	
	/* Allocate memory */
	self = (PICelleratorContext*)_FiniteElementContext_New( 
		sizeOfSelf, 
		type, 
		_delete, 
		_print, 
		_copy,
		_defaultConstructor,
		_construct,
		_build,
		_initialise,
		_execute,
		_destroy,
		name,
		initFlag,
		_setDt, 
		start, 
		stop, 
		communicator, 
		dictionary );
	
	/* General info */
	
	/* Virtual info */
	
	if( initFlag ){
		_PICelleratorContext_Init( self );
	}
	
	return self;
}


void _PICelleratorContext_Init( void* context ) {
	PICelleratorContext* self = (PICelleratorContext*)context;
	self->isConstructed = True;

	self->materials_Register = Materials_Register_New();

	ContextEP_Prepend( self, AbstractContext_EP_AssignFromXMLExtensions, PICelleratorContext_CreateDefaultMaterial );

/* 	 TODO want to append an EP to the end of the time integration that makes sure that after all integration  */
/* 	 swarms have been updated, each element has at least one integration point from one swarm in it, or else */
/* 	 we will fail the next timestep */
}


/* Virtual Functions -------------------------------------------------------------------------------------------------------------*/

void _PICelleratorContext_Delete( void* context ) {
	PICelleratorContext* self = (PICelleratorContext*)context;
	
	Journal_DPrintf( self->debug, "In: %s()\n", __func__ );

	Stg_Class_Delete( self->materials_Register );

	/* Stg_Class_Delete parent */
	_FiniteElementContext_Delete( self );
}

void _PICelleratorContext_Destroy( void* context ) {
	PICelleratorContext* self = (PICelleratorContext*)context;
	
	_FiniteElementContext_Destroy( self );
}

void _PICelleratorContext_Print( void* context, Stream* stream ) {
	PICelleratorContext* self = (PICelleratorContext*)context;
	
	/* General info */
	Journal_Printf( (void*) stream, "PICelleratorContext (ptr): %p\n", self );
	
	/* Print parent */
	_FiniteElementContext_Print( self, stream );

	Journal_PrintPointer( stream, self->materials_Register );

}

void _PICelleratorContext_SetDt( void* context, double dt ) {
	PICelleratorContext* self = (PICelleratorContext*)context;
	
	self->dt = dt;
}


/* Public Functions ----------------------------------------------------------------------------------------------------*/



/* EntryPoint Hooks ----------------------------------------------------------------------------------------------------*/
void PICelleratorContext_CreateDefaultMaterial( void* context ) {
	PICelleratorContext*      self         = (PICelleratorContext*) context;
	
	if ( Materials_Register_GetCount( self->materials_Register ) == 0 ) {
		Stg_Shape* everywhereShape = (Stg_Shape*) Everywhere_New( "defaultShape", self->dim );

		Material_New( "backgroundMaterial", everywhereShape, self->dictionary, self->materials_Register );
	}
}

void _PICelleratorContext_AssignFromXML( void* context, Stg_ComponentFactory *cf, void* data ){
	PICelleratorContext* self = (PICelleratorContext*) context;
	
	_FiniteElementContext_AssignFromXML( context, cf, data );
	_PICelleratorContext_Init( self );
}
