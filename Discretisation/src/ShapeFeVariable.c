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
** $Id: ShapeFeVariable.c 532 2006-04-04 00:21:59Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <string.h>
#include <mpi.h>
#include <StGermain/StGermain.h>

#include "types.h"
#include "Mesh.h"
#include "FeVariable.h"
#include "ShapeFeVariable.h"

#include <assert.h>

const Type ShapeFeVariable_Type = "ShapeFeVariable";

void* ShapeFeVariable_DefaultNew( Name name ) {
	return (ShapeFeVariable*)
		_ShapeFeVariable_New(
			sizeof(ShapeFeVariable), 
			ShapeFeVariable_Type, 
			_ShapeFeVariable_Delete,
			_ShapeFeVariable_Print,
			_ShapeFeVariable_Copy,
			(Stg_Component_DefaultConstructorFunction*)ShapeFeVariable_DefaultNew,
			_ShapeFeVariable_Construct,
			_ShapeFeVariable_Build,
			_ShapeFeVariable_Initialise,
			_ShapeFeVariable_Execute,
			_ShapeFeVariable_Destroy,
			name,
			False,
			_FeVariable_InterpolateValueAt,
			_FeVariable_GetMinGlobalFieldMagnitude,
			_FeVariable_GetMaxGlobalFieldMagnitude, 
			_FeVariable_GetMinAndMaxLocalCoords,
			_FeVariable_GetMinAndMaxGlobalCoords,
			_FeVariable_InterpolateNodeValuesToElLocalCoord,
			_FeVariable_GetValueAtNode,
			NULL,
			NULL,
			NULL,
			0,	/* dim */
			False,
			0,	/* communicator */
			NULL	/* fv_Register */
			);
}			

ShapeFeVariable* _ShapeFeVariable_New(
 		SizeT                                             _sizeOfSelf,
		Type                                              type,
		Stg_Class_DeleteFunction*                         _delete,
		Stg_Class_PrintFunction*                          _print,
		Stg_Class_CopyFunction*                           _copy, 
		Stg_Component_DefaultConstructorFunction*         _defaultConstructor,
		Stg_Component_ConstructFunction*                  _construct,
		Stg_Component_BuildFunction*                      _build,
		Stg_Component_InitialiseFunction*                 _initialise,
		Stg_Component_ExecuteFunction*                    _execute,
		Stg_Component_DestroyFunction*                    _destroy,
		Name                                              name,
		Bool                                              initFlag,
		FieldVariable_InterpolateValueAtFunction*         _interpolateValueAt,
		FieldVariable_GetValueFunction*	                  _getMinGlobalFeMagnitude,
		FieldVariable_GetValueFunction*                   _getMaxGlobalFeMagnitude,
		FieldVariable_GetCoordFunction*                   _getMinAndMaxLocalCoords,
		FieldVariable_GetCoordFunction*                   _getMinAndMaxGlobalCoords,		
		FeVariable_InterpolateWithinElementFunction*      _interpolateWithinElement,	
		FeVariable_GetValueAtNodeFunction*                _getValueAtNode,
		void*                                              feMesh,
		void*                                              geometryMesh,
		void*                                              dofLayout,
		Dimension_Index                                    dim,
		Bool                                               isCheckpointedAndReloaded,
		MPI_Comm                                           communicator,
		FieldVariable_Register*                            fV_Register
	       	)
{
	ShapeFeVariable*		self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(ShapeFeVariable) );
	self = (ShapeFeVariable*)
		_FeVariable_New(
			_sizeOfSelf, 
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
			False,
			_interpolateValueAt,
			_getMinGlobalFeMagnitude, 
			_getMaxGlobalFeMagnitude,
			_getMinAndMaxLocalCoords, 
			_getMinAndMaxGlobalCoords,
			_interpolateWithinElement,
			_getValueAtNode,
			feMesh,
			geometryMesh,
			dofLayout, /* dofLayout */
			NULL,   /* bcs */
			NULL,   /* ics */
			NULL,   /* linkedDofInfo */
			NULL,   /* templateFeVariable */
			1,      /* fieldComponentCount */ 
			dim,	/* dim */
			isCheckpointedAndReloaded, /* Checkpointing boolean */
			StgFEM_Native_ImportExportType,	/* import format type */
			StgFEM_Native_ImportExportType,	/* export format type */
			communicator,	/* communicator */
			fV_Register	/* fv_Register */
			);

	return self;
}

void _ShapeFeVariable_Init( void* shapeFeVariable, Stg_Shape* shape ) {
	ShapeFeVariable*         self              = (ShapeFeVariable*) shapeFeVariable;

	self->shape = shape;
	
	//EP_AppendClassHook( Context_GetEntryPoint( context, AbstractContext_EP_UpdateClass ),	ParticleFeVariable_Update, self );
}

void _ShapeFeVariable_Delete( void* _shapeFeVariable ) {
	ShapeFeVariable* self = (ShapeFeVariable*) _shapeFeVariable;

	_FeVariable_Delete( self );
}

void _ShapeFeVariable_Print( void* _shapeFeVariable, Stream* stream ) {
	ShapeFeVariable* self = (ShapeFeVariable*) _shapeFeVariable;

	/* Print parent */
	_FeVariable_Print( self, stream );

	Journal_PrintPointer( stream, self->shape );
}


void* _ShapeFeVariable_Copy( void* shapeFeVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	//ShapeFeVariable*	self = (ShapeFeVariable*)shapeFeVariable;
	ShapeFeVariable*	newShapeFeVariable;
	
	assert(0);
	return (void*)newShapeFeVariable;
}

void _ShapeFeVariable_Construct( void* shapeFeVariable, Stg_ComponentFactory* cf, void* data ) {
	ShapeFeVariable*        self       = (ShapeFeVariable*) shapeFeVariable;

	_FeVariable_Construct( self, cf, data );

	_ShapeFeVariable_Init( 
		self, 
		Stg_ComponentFactory_ConstructByKey(  cf,  self->name,  "Shape", Stg_Shape,  True, data )  ) ;
}

void _ShapeFeVariable_Build( void* shapeFeVariable, void* data ) {
	ShapeFeVariable* self = (ShapeFeVariable*) shapeFeVariable;

	Stg_Component_Build( self->shape, data, False );
	_FeVariable_Build( self, data );
}

void _ShapeFeVariable_Initialise( void* shapeFeVariable, void* data ) {
	ShapeFeVariable* self = (ShapeFeVariable*) shapeFeVariable;
	Node_DomainIndex node_dI = 0;

	Stg_Component_Initialise( self->shape, data, False );
	_FeVariable_Initialise( self, data );

	/* Set up the basic "level set" describiing if nodes are inside the shape or not */
	for ( node_dI = 0; node_dI < self->feMesh->nodeDomainCount; node_dI++ ) {
		if ( True == Stg_Shape_IsCoordInside( self->shape, self->feMesh->nodeCoord[node_dI] ) ) {
			//set value = 1
			FeVariable_SetComponentAtNode( self, node_dI, 0, 1 );
		}		
		else {
			//set value = 0
			FeVariable_SetComponentAtNode( self, node_dI, 0, 0 );
		}
	}
}

void _ShapeFeVariable_Execute( void* shapeFeVariable, void* data ) {
	ShapeFeVariable* self = (ShapeFeVariable*) shapeFeVariable;
	
	_FeVariable_Execute( self, data );
}

void _ShapeFeVariable_Destroy( void* shapeFeVariable, void* data ) {
	ShapeFeVariable* self = (ShapeFeVariable*) shapeFeVariable;
	
	_FeVariable_Destroy( self, data );
}
