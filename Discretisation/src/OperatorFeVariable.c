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
** $Id: OperatorFeVariable.c 920 2007-07-20 06:19:34Z RobertTurnbull $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <string.h>
#include <mpi.h>
#include <StGermain/StGermain.h>

#include "types.h"
#include "FeVariable.h"
#include "OperatorFeVariable.h"
#include "FeMesh.h"

#include <assert.h>

const Type OperatorFeVariable_Type = "OperatorFeVariable";

OperatorFeVariable* OperatorFeVariable_NewUnary( 
		Name                                               name ,
		void*                                              _feVariable,
		Name                                               operatorName )
{
	FeVariable* feVariable  = (FeVariable*) _feVariable;
	Stream*     errorStream = Journal_Register( Error_Type, OperatorFeVariable_Type );

	Journal_Firewall( feVariable != NULL, errorStream, "In func %s: Trying to operate on NULL field.\n", __func__ );
       	
	return OperatorFeVariable_New( 
			name,
			OperatorFeVariable_UnaryInterpolationFunc, 
			OperatorFeVariable_UnaryValueAtNodeFunc,
			operatorName,
			1,
			&feVariable, 
			feVariable->dim,
			feVariable->isCheckpointedAndReloaded,
			feVariable->communicator,
			feVariable->fieldVariable_Register );
}

OperatorFeVariable* OperatorFeVariable_NewBinary( 
		Name                                               name ,
		void*                                              _feVariable1,
		void*                                              _feVariable2,
		Name                                               operatorName )
{
	FeVariable* feVariableList[2];
	Stream*     errorStream = Journal_Register( Error_Type, OperatorFeVariable_Type );
	
	Journal_Firewall( _feVariable1 != NULL, errorStream, "In func %s: First field to operate on is NULL.\n", __func__ );
	Journal_Firewall( _feVariable2 != NULL, errorStream, "In func %s: Second field to operate on is NULL.\n", __func__ );
       
	feVariableList[0] = (FeVariable*) _feVariable1;
	feVariableList[1] = (FeVariable*) _feVariable2;
	
	return OperatorFeVariable_New( 
			name,
			OperatorFeVariable_BinaryInterpolationFunc, 
			OperatorFeVariable_BinaryValueAtNodeFunc,
			operatorName,
			2, 
			feVariableList, 
			feVariableList[0]->dim,
			feVariableList[0]->isCheckpointedAndReloaded,
			feVariableList[0]->communicator,
			feVariableList[0]->fieldVariable_Register );
}

void* OperatorFeVariable_DefaultNew( Name name ) {
		return _OperatorFeVariable_New( 
			sizeof(OperatorFeVariable), 
			OperatorFeVariable_Type, 
			_OperatorFeVariable_Delete, 
			_OperatorFeVariable_Print,
			_OperatorFeVariable_Copy, 
			(Stg_Component_DefaultConstructorFunction*)OperatorFeVariable_DefaultNew,
			_OperatorFeVariable_Construct,
			_OperatorFeVariable_Build, 
			_OperatorFeVariable_Initialise, 
			_OperatorFeVariable_Execute,
			_OperatorFeVariable_Destroy,
			name,
			False,
			_OperatorFeVariable_InterpolateValueAt,
			_FeVariable_GetMinGlobalFieldMagnitude,
			_FeVariable_GetMaxGlobalFieldMagnitude, 
			_FeVariable_GetMinAndMaxLocalCoords,
			_FeVariable_GetMinAndMaxGlobalCoords,
			_OperatorFeVariable_InterpolateWithinElement,
			_OperatorFeVariable_GetValueAtNode,
			NULL,
			0,
			NULL,
			NULL,
			NULL,
			0,
			False,
			MPI_COMM_WORLD,
			NULL );
}

OperatorFeVariable* OperatorFeVariable_New( 
		Name                                               name ,
		FeVariable_InterpolateWithinElementFunction*       interpolateWithinElement,
		FeVariable_GetValueAtNodeFunction*                 getValueAtNode,
		Name                                               operatorName,
		Index                                              feVariableCount,
		FeVariable**                                       feVariableList,
		Dimension_Index                                    dim,
		Bool                                               isCheckpointedAndReloaded,
		MPI_Comm                                           communicator,
		FieldVariable_Register*                            fV_Register )
{
	FeVariable* feVariable = feVariableList[0];

		return _OperatorFeVariable_New( 
			sizeof(OperatorFeVariable), 
			OperatorFeVariable_Type, 
			_OperatorFeVariable_Delete, 
			_OperatorFeVariable_Print,
			_OperatorFeVariable_Copy, 
			(Stg_Component_DefaultConstructorFunction*)OperatorFeVariable_DefaultNew,
			_OperatorFeVariable_Construct,
			_OperatorFeVariable_Build, 
			_OperatorFeVariable_Initialise, 
			_OperatorFeVariable_Execute,
			_OperatorFeVariable_Destroy,
			name,
			True,
			_OperatorFeVariable_InterpolateValueAt,
			_FeVariable_GetMinGlobalFieldMagnitude,
			_FeVariable_GetMaxGlobalFieldMagnitude, 
			_FeVariable_GetMinAndMaxLocalCoords,
			_FeVariable_GetMinAndMaxGlobalCoords,
			interpolateWithinElement,
			getValueAtNode,
			operatorName,
			feVariableCount,
			feVariableList,
			feVariable->feMesh,
			feVariable->geometryMesh,
			dim,
			isCheckpointedAndReloaded,
			communicator,
			fV_Register );
}

OperatorFeVariable* _OperatorFeVariable_New(
 		SizeT                                              _sizeOfSelf, 
		Type                                               type,
		Stg_Class_DeleteFunction*                          _delete,
		Stg_Class_PrintFunction*                           _print, 
		Stg_Class_CopyFunction*                            _copy, 
		Stg_Component_DefaultConstructorFunction*          _defaultConstructor,
		Stg_Component_ConstructFunction*                   _construct,
		Stg_Component_BuildFunction*                       _build,
		Stg_Component_InitialiseFunction*                  _initialise,
		Stg_Component_ExecuteFunction*                     _execute,
		Stg_Component_DestroyFunction*                     _destroy,
		Name                                               name,
		Bool										       initFlag,
		FieldVariable_InterpolateValueAtFunction*          _interpolateValueAt,
		FieldVariable_GetValueFunction*			   _getMinGlobalFieldMagnitude,
		FieldVariable_GetValueFunction*                    _getMaxGlobalFieldMagnitude,
		FieldVariable_GetCoordFunction*                    _getMinAndMaxLocalCoords,
		FieldVariable_GetCoordFunction*                    _getMinAndMaxGlobalCoords,
		FeVariable_InterpolateWithinElementFunction*       interpolateWithinElement,
		FeVariable_GetValueAtNodeFunction*                 getValueAtNode,
		Name                                               operatorName,
		Index                                              feVariableCount,
		FeVariable**                                       feVariableList,
		void*                                              feMesh,
		void*                                              geometryMesh,
		Dimension_Index                                    dim,
		Bool                                               isCheckpointedAndReloaded,
		MPI_Comm                                           communicator,
		FieldVariable_Register*                            fV_Register )
{
	OperatorFeVariable*		    self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(OperatorFeVariable) );
	self = (OperatorFeVariable*) _FeVariable_New( 
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
			initFlag,
			_interpolateValueAt,
			_getMinGlobalFieldMagnitude,
			_getMaxGlobalFieldMagnitude,
			_getMinAndMaxLocalCoords,
			_getMinAndMaxGlobalCoords,
			interpolateWithinElement,
			getValueAtNode,
			_OperatorFeVariable_SyncShadowValues, 
			feMesh,
			geometryMesh,
			NULL, /* dofLayout */
			NULL, /* BCs */
			NULL, /* ICs */
			NULL, /* Linked dof info */
			NULL, /* Template FeVariable */
			0 /* fieldComponentCount - this will be reset later */,
			dim,
			False, /* isCheckpointedAndReloaded */
			/* TODO: hack as always StgFEM_Native for now - PatrickSunter - 7 July 2006 */
			StgFEM_Native_ImportExportType,
			StgFEM_Native_ImportExportType,
			NULL,
			/* TODO: as above, potentially might want to save to non-standard paths... fix
			 * later though.- PatrickSunter - 9 July 2007 */
			NULL,
			False, /* Not a reference variable, so this line and following False */
			False, 
			communicator,
			fV_Register );

	if( initFlag ){
		_OperatorFeVariable_Init( self, operatorName, feVariableCount, feVariableList );
	}

	return self;
}

void _OperatorFeVariable_Init( void* oFeVar, Name operatorName, Index feVariableCount, FeVariable** feVariableList ) {
	OperatorFeVariable*         self              = (OperatorFeVariable*) oFeVar;
	FeVariable*                 feVariable;
	Index                       feVariable_I;
	Stream*                     errorStream       = Journal_Register( Error_Type, self->type );

	/* Assign values to object */
	self->feVariableCount     = feVariableCount;
	self->operatorName = operatorName;

	/* Copy field variable list */
	self->feVariableList      = Memory_Alloc_Array( FeVariable*, feVariableCount, "Array of Field Variables" );
	memcpy( self->feVariableList, feVariableList, feVariableCount * sizeof( FeVariable* ) );

	for ( feVariable_I = 0 ; feVariable_I < feVariableCount ; feVariable_I++ ) {
		feVariable = feVariableList[ feVariable_I ];
		Journal_Firewall( feVariable != NULL,
				errorStream, "In func %s: Field Variable %u in list is NULL.\n", __func__, feVariable_I );
		Journal_Firewall( feVariable->fieldComponentCount <= MAX_FIELD_COMPONENTS,
				errorStream, "In func %s: Field Variable '%s' has too many components.\n", __func__, feVariable->name );
	}
}

void _OperatorFeVariable_Delete( void* _feVariable ) {
	OperatorFeVariable* self = (OperatorFeVariable*) _feVariable;

	Memory_Free( self->feVariableList );

	/* TODO - HACK - Should be FeVariable_Delete */
	_FeVariable_Delete( self );
}

void _OperatorFeVariable_Print( void* _feVariable, Stream* stream ) {
	OperatorFeVariable* self = (OperatorFeVariable*) _feVariable;
	Index               feVariable_I;

	_FeVariable_Print( self, stream );

	Journal_PrintValue( stream, self->feVariableCount );
	for ( feVariable_I = 0 ; feVariable_I < self->feVariableCount ; feVariable_I++ ) 
		Journal_Printf( stream, "\tFeVariable %u - '%s'\n", feVariable_I, self->feVariableList[ feVariable_I ]->name );

}


void* _OperatorFeVariable_Copy( void* feVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	OperatorFeVariable*	self = (OperatorFeVariable*)feVariable;
	OperatorFeVariable*	newOperatorFeVariable;
	
	newOperatorFeVariable = _FeVariable_Copy( self, dest, deep, nameExt, ptrMap );
	
	newOperatorFeVariable->_operator              = self->_operator;
	newOperatorFeVariable->feVariableCount     = self->feVariableCount;
	
	if (deep) {
		newOperatorFeVariable->feVariableList = Memory_Alloc_Array( FeVariable*, self->feVariableCount, 
				"Array of Field Variables" );
		memcpy( newOperatorFeVariable->feVariableList, self->feVariableList, 
				self->feVariableCount * sizeof( FeVariable* ) );
	}
	else 
		newOperatorFeVariable->feVariableList = self->feVariableList;
	
	return (void*)newOperatorFeVariable;
}

void _OperatorFeVariable_Construct( void* feVariable, Stg_ComponentFactory* cf, void* data ) {
	OperatorFeVariable*     self       = (OperatorFeVariable*) feVariable;
	Dictionary*             dictionary = Dictionary_GetDictionary( cf->componentDict, self->name );
	Dictionary_Entry_Value* list;
	Index                   feVariableCount = 0;
	Index                   feVariable_I;
	FieldVariable_Register* fV_Register     = Stg_ObjectList_Get( cf->registerRegister, "FieldVariable_Register" );
	Name                    feVariableName;
	Name                    operatorName;
	FeVariable**            feVariableList;
	
	/* Construct Parent */
	_FieldVariable_Construct( self, cf, data );

	operatorName = Stg_ComponentFactory_GetString( cf, self->name, "Operator", "" );

	list = Dictionary_Get( dictionary, "FeVariables" );
	
	feVariableCount = ( list ? Dictionary_Entry_Value_GetCount(list) : 1 );
	feVariableList = Memory_Alloc_Array( FeVariable*, feVariableCount, "FeVars" );

	for ( feVariable_I = 0 ; feVariable_I < feVariableCount ; feVariable_I++ ) {
		feVariableName = (list ? 
				Dictionary_Entry_Value_AsString( Dictionary_Entry_Value_GetElement( list, feVariable_I ) ) :
				Dictionary_GetString( dictionary, "FeVariable" ) );

		/* Check in fV_Register first before assuming in LiveComponentRegister */
		Journal_PrintfL( cf->infoStream, 2, "Looking for FeVariable '%s' in fieldVariable_Register.\n",
				feVariableName );
		feVariableList[feVariable_I] = (FeVariable*) FieldVariable_Register_GetByName( fV_Register, feVariableName );
		
		if ( !feVariableList[feVariable_I] )
			feVariableList[feVariable_I] = 
	Stg_ComponentFactory_ConstructByName( cf, feVariableName, FeVariable, True, data ); 
		
	}

	_FeVariable_Init( (FeVariable*) self, feVariableList[0]->feMesh, feVariableList[0]->geometryMesh,
		feVariableList[0]->dofLayout, NULL, NULL, NULL, NULL,
		/* TODO: hack as always StgFEM native for now - PatrickSunter 7/7/2006 */
		StgFEM_Native_ImportExportType, StgFEM_Native_ImportExportType,
		/* TODO: hack as always the default path - PatrickSunter 9/7/2007 */
		NULL, NULL,
		False, False );
	_OperatorFeVariable_Init( self, operatorName, feVariableCount, feVariableList );

	Memory_Free( feVariableList );
}

void _OperatorFeVariable_Build( void* feVariable, void* data ) {
	OperatorFeVariable* self = (OperatorFeVariable*) feVariable;
	Index                  feVariable_I;
	Stream*                     errorStream       = Journal_Register( Error_Type, self->type );

	for ( feVariable_I = 0 ; feVariable_I < self->feVariableCount ; feVariable_I++ ) 
		Stg_Component_Build( self->feVariableList[ feVariable_I ] , data, False );

	/* Check if we are using a gradient operator */
	if ( strcasecmp( self->operatorName, "gradient" ) == 0 ) {
		self->useGradient = True;
		assert( self->feVariableCount == 1 );
		self->fieldComponentCount = self->feVariableList[0]->fieldComponentCount * self->dim;
	}
	else {
		/* just use normal operator */
		self->useGradient = False;
		/* Added 5 May 2006 by P Sunter: in the case of VectorScale, the fieldComponentCount should be based
		on the 2nd operator. Also make sure the 2nd operator has at least as may dofs per node as the first. */
		if ( self->feVariableCount == 2 ) {
			Journal_Firewall( self->feVariableList[1]->fieldComponentCount >= self->feVariableList[0]->fieldComponentCount,
				errorStream, "Error - in %s: tried to create a %s operator from feVariables \"%s\" "
				"and \"%s\" - who have fieldCounts %d and %d - unsupported since operations "
				"such as VectorScale require the 2nd feVariable to have >= the first's field count.\n",
				__func__, self->operatorName, self->feVariableList[0]->name, self->feVariableList[1]->name,
				self->feVariableList[0]->fieldComponentCount, self->feVariableList[1]->fieldComponentCount );
			self->_operator  = Operator_NewFromName( self->operatorName, self->feVariableList[1]->fieldComponentCount,
				self->dim );
		}
		else {	
			self->_operator  = Operator_NewFromName( self->operatorName, self->feVariableList[0]->fieldComponentCount,
				self->dim );
		}

		self->fieldComponentCount = self->_operator->resultDofs; /* reset this value from that which is from operator */
	}
	_OperatorFeVariable_SetFunctions( self );
}

void _OperatorFeVariable_Initialise( void* feVariable, void* data ) {
	OperatorFeVariable* self = (OperatorFeVariable*) feVariable;
	Index                  feVariable_I;

	for ( feVariable_I = 0 ; feVariable_I < self->feVariableCount ; feVariable_I++ ) 
		Stg_Component_Initialise( self->feVariableList[ feVariable_I ] , data, False );
}

void _OperatorFeVariable_Execute( void* feVariable, void* data ) {}
void _OperatorFeVariable_Destroy( void* feVariable, void* data ) {}

void _OperatorFeVariable_SetFunctions( void* feVariable ) {
	OperatorFeVariable* self            = (OperatorFeVariable*) feVariable;
	Stream*             error           = Journal_Register( Error_Type, self->type );

	if ( self->useGradient ) {
		Journal_Firewall( self->feVariableCount == 1, error, "Cannot use gradient operators for multiple variables.\n" );
		self->_interpolateWithinElement = OperatorFeVariable_GradientInterpolationFunc;
		self->_getValueAtNode = OperatorFeVariable_GradientValueAtNodeFunc;
	}
	else {
		switch ( self->feVariableCount ) {
			case 1:
				self->_interpolateWithinElement = OperatorFeVariable_UnaryInterpolationFunc; 
				self->_getValueAtNode = OperatorFeVariable_UnaryValueAtNodeFunc; break;
			case 2:
				self->_interpolateWithinElement = OperatorFeVariable_BinaryInterpolationFunc; 
				self->_getValueAtNode = OperatorFeVariable_BinaryValueAtNodeFunc; break;
			default: {
				Journal_Firewall( False, error,
						"Cannot use '%s' with feVariableCount = %u.\n", __func__, self->feVariableCount );
			}
		}
	}
}
void _OperatorFeVariable_InterpolateWithinElement( void* feVariable, Element_DomainIndex dElement_I, Coord coord, double* value ) {
	OperatorFeVariable* self            = (OperatorFeVariable*) feVariable;

	_OperatorFeVariable_SetFunctions( self );
	FeVariable_InterpolateWithinElement( self, dElement_I, coord, value );
}

void _OperatorFeVariable_GetValueAtNode( void* feVariable, Node_DomainIndex dNode_I, double* value ) {
	OperatorFeVariable* self            = (OperatorFeVariable*) feVariable;

	_OperatorFeVariable_SetFunctions( self );
	FeVariable_GetValueAtNode( self, dNode_I, value );
}

void _OperatorFeVariable_SyncShadowValues( void* feVariable ) {
	OperatorFeVariable* self            = (OperatorFeVariable*) feVariable;
	int v_i;

	assert( self );
	for( v_i = 0; v_i < self->feVariableCount; v_i++ )
		FeVariable_SyncShadowValues( self->feVariableList[v_i] );
	self->shadowValuesSynchronised = True;
}


/** Private Functions */
Bool _OperatorFeVariable_CheckIfValidToInterpolateInShadowSpace( OperatorFeVariable* self ){
	int                     feVar_I;
	FeVariable*             currFeVar;
	Bool                    parentIsValid;

	for ( feVar_I=0; feVar_I < self->feVariableCount; feVar_I++ ) {
		currFeVar = self->feVariableList[feVar_I];
		if( Stg_Class_IsInstance( currFeVar, FeVariable_Type ) ) {
			if ( False == currFeVar->shadowValuesSynchronised ) {
        	        	return False;
        	}
			else {
				return True;
			}
		}
		else {
			parentIsValid = _OperatorFeVariable_CheckIfValidToInterpolateInShadowSpace(
				(OperatorFeVariable*)currFeVar );
			if ( False == parentIsValid ) {
				return False;
			}
		} 
	}

	return True;
}


/** Needed to over-ride the standard _FeVariable_InterpolateValueAt() since for operator fe variables, as long as the
base fe variable has been shadowed, this operatorFeVar can be calculated in shadow space */
InterpolationResult _OperatorFeVariable_InterpolateValueAt( void* variable, Coord globalCoord, double* value ) {
	OperatorFeVariable*     self = (OperatorFeVariable*)variable;
	Element_DomainIndex     elementCoordIn = (unsigned)-1;
	Coord                   elLocalCoord={0,0,0};
	InterpolationResult     retValue;
	Bool                    validToInterpolateInShadowSpace;


	retValue = FeVariable_GetElementLocalCoordAtGlobalCoord( self, globalCoord, elLocalCoord, &elementCoordIn );

	if ( retValue == LOCAL ) {
		/* Now interpolate the value at that coordinate, using shape functions */
		self->_interpolateWithinElement( self, elementCoordIn, elLocalCoord, value );
	}
	else if ( retValue == SHADOW ) {
		validToInterpolateInShadowSpace = _OperatorFeVariable_CheckIfValidToInterpolateInShadowSpace(self);	
			
		if ( False == validToInterpolateInShadowSpace ) {
			Stream* warningStr = Journal_Register( Error_Type, self->type );
			Journal_Printf( warningStr, "Warning - in %s, for OperatorFeVar \"%s\": user asking to "
				"interpolate a value at "
				"coord (%g,%g,%g), which is in shadow space, but "
				"FeVariable_SyncShadowValues() hasn't been called on FeVariables this one is "
				"derived from yet.\n",
				__func__, self->name, globalCoord[0], globalCoord[1], globalCoord[2] );
			return retValue;
		}
		/* Now interpolate the value at that coordinate, using shape functions */
		self->_interpolateWithinElement( self, elementCoordIn, elLocalCoord, value );
	}

	return retValue;
}



void OperatorFeVariable_UnaryInterpolationFunc( void* feVariable, Element_DomainIndex dElement_I, Coord coord, double* value ) {
	OperatorFeVariable* self            = (OperatorFeVariable*) feVariable;
	FeVariable*         field0          = self->feVariableList[0];
	double              fieldValue0[ MAX_FIELD_COMPONENTS ]; 

	field0->_interpolateWithinElement( field0, dElement_I, coord, fieldValue0 );
	Operator_CarryOutUnaryOperation( self->_operator, fieldValue0, value );
}

void OperatorFeVariable_BinaryInterpolationFunc( void* feVariable, Element_DomainIndex dElement_I, Coord localCoord, double* value ) {
	OperatorFeVariable* self            = (OperatorFeVariable*) feVariable;
	FeVariable*         field0          = self->feVariableList[0];
	FeVariable*         field1          = self->feVariableList[1];
	FeMesh*             mesh            = self->feMesh;
	double              fieldValue0[ MAX_FIELD_COMPONENTS ]; 
	double              fieldValue1[ MAX_FIELD_COMPONENTS ]; 
	Coord               globalCoord;
	
	/* Get value for field0 */
	if ( field0->feMesh == mesh )
		field0->_interpolateWithinElement( field0, dElement_I, localCoord, fieldValue0 );
	else {
		/* Get Global Coordinates */
		FeMesh_CoordLocalToGlobal( mesh, dElement_I, localCoord, globalCoord );
		field0->_interpolateValueAt( field0, globalCoord, fieldValue0 );
	}

	/* Get value for field1 */
	if ( field1->feMesh == mesh )
		field1->_interpolateWithinElement( field1, dElement_I, localCoord, fieldValue1 );
	else {
		/* Get Global Coordinates */
		FeMesh_CoordLocalToGlobal( mesh, dElement_I, localCoord, globalCoord );
		field1->_interpolateValueAt( field1, globalCoord, fieldValue1 );
	}
	
	Operator_CarryOutBinaryOperation( self->_operator, fieldValue0, fieldValue1, value );
}

void OperatorFeVariable_GradientInterpolationFunc( void* feVariable, Element_DomainIndex dElement_I, Coord coord, double* value ) {
	OperatorFeVariable* self            = (OperatorFeVariable*) feVariable;
	FeVariable*         field0          = self->feVariableList[0];
	
	FeVariable_InterpolateDerivativesToElLocalCoord( field0, dElement_I, coord, value );
}

void OperatorFeVariable_UnaryValueAtNodeFunc( void* feVariable, Node_DomainIndex dNode_I, double* value ) {
	OperatorFeVariable* self            = (OperatorFeVariable*) feVariable;
	FeVariable*         field0          = self->feVariableList[0];
	double              fieldValue0[ MAX_FIELD_COMPONENTS ]; 
	
	FeVariable_GetValueAtNode( field0, dNode_I, fieldValue0 );
	Operator_CarryOutUnaryOperation( self->_operator, fieldValue0, value );
}

void OperatorFeVariable_BinaryValueAtNodeFunc( void* feVariable, Node_DomainIndex dNode_I, double* value ) {
	OperatorFeVariable* self            = (OperatorFeVariable*) feVariable;
	FeVariable*         field0          = self->feVariableList[0];
	FeVariable*         field1          = self->feVariableList[1];
	double              fieldValue0[ MAX_FIELD_COMPONENTS ]; 
	double              fieldValue1[ MAX_FIELD_COMPONENTS ]; 
	
	FeVariable_GetValueAtNode( field0, dNode_I, fieldValue0 );
	FeVariable_GetValueAtNode( field1, dNode_I, fieldValue1 );
	Operator_CarryOutBinaryOperation( self->_operator, fieldValue0, fieldValue1, value );
}

void OperatorFeVariable_GradientValueAtNodeFunc( void* feVariable, Node_DomainIndex dNode_I, double* value ) {
	OperatorFeVariable* self            = (OperatorFeVariable*) feVariable;
	Mesh*               mesh            = (Mesh*) self->feMesh;
	double*             coord;

	coord = Mesh_GetVertex( mesh, dNode_I );

	memset( value, 0, self->fieldComponentCount * sizeof(double) );
	FeVariable_InterpolateDerivativesAt( self->feVariableList[0], coord, value );
}
