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
** $Id: OperatorFeVariable.c 1137 2008-05-23 05:57:48Z RobertTurnbull $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <string.h>
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "types.h"
#include "FeVariable.h"
#include "OperatorFeVariable.h"
#include "FeMesh.h"

#include <assert.h>

const Type OperatorFeVariable_Type = "OperatorFeVariable";

OperatorFeVariable* OperatorFeVariable_NewUnary(
	Name				name,
	DomainContext*	context,
	void*				_feVariable,
	Name				operatorName )
{
	FeVariable*	feVariable = (FeVariable*) _feVariable;
	Stream*		errorStream = Journal_Register( Error_Type, OperatorFeVariable_Type );

	Journal_Firewall( feVariable != NULL, errorStream, "In func %s: Trying to operate on NULL field.\n", __func__ );
       	
	return OperatorFeVariable_New( 
		name,
		context,
		feVariable->feMesh,
  		feVariable->geometryMesh,
  		feVariable->dofLayout,                                                                                 
  		feVariable->bcs,
  		feVariable->ics,
  		feVariable->linkedDofInfo,
  		feVariable->templateFeVariable,    
		OperatorFeVariable_UnaryInterpolationFunc, 
		OperatorFeVariable_UnaryValueAtNodeFunc,
		operatorName,
		NULL,
		1,
		&feVariable, 
		feVariable->dim,
		feVariable->isCheckpointedAndReloaded,
		feVariable->communicator,
		feVariable->fieldVariable_Register );
}

OperatorFeVariable* OperatorFeVariable_NewUnary_OwnOperator(
	Name				name,
	DomainContext*	context,
	void*				_feVariable,
	Operator*		ownOperator )
{
	FeVariable*	feVariable = (FeVariable*) _feVariable;
	Stream*		errorStream = Journal_Register( Error_Type, OperatorFeVariable_Type );

	Journal_Firewall( feVariable != NULL, errorStream, "In func %s: Trying to operate on NULL field.\n", __func__ );
       	
	return OperatorFeVariable_New( 
		name,
		context,
		feVariable->feMesh,
  		feVariable->geometryMesh,
  		feVariable->dofLayout,                                                                                 
  		feVariable->bcs,
  		feVariable->ics,
  		feVariable->linkedDofInfo,
  		feVariable->templateFeVariable,    
		OperatorFeVariable_UnaryInterpolationFunc, 
		OperatorFeVariable_UnaryValueAtNodeFunc,
		ownOperator->name,
		ownOperator,
		1,
		&feVariable, 
		feVariable->dim,
		feVariable->isCheckpointedAndReloaded,
		feVariable->communicator,
		feVariable->fieldVariable_Register );
}

OperatorFeVariable* OperatorFeVariable_NewBinary(
	Name				name,
	DomainContext*	context,
	void*				_feVariable1,
	void*				_feVariable2,
	Name				operatorName )
{
	FeVariable*	feVariableList[2];
	Stream*		errorStream = Journal_Register( Error_Type, OperatorFeVariable_Type );
	
	Journal_Firewall( _feVariable1 != NULL, errorStream, "In func %s: First field to operate on is NULL.\n", __func__ );
	Journal_Firewall( _feVariable2 != NULL, errorStream, "In func %s: Second field to operate on is NULL.\n", __func__ );
       
	feVariableList[0] = (FeVariable*) _feVariable1;
	feVariableList[1] = (FeVariable*) _feVariable2;
	
	return OperatorFeVariable_New( 
		name,
		context,
		feVariableList[0]->feMesh,
  		feVariableList[0]->geometryMesh,
  		feVariableList[0]->dofLayout,                                                                                 
  		feVariableList[0]->bcs,
  		feVariableList[0]->ics,
  		feVariableList[0]->linkedDofInfo,
  		feVariableList[0]->templateFeVariable,    
		OperatorFeVariable_BinaryInterpolationFunc, 
		OperatorFeVariable_BinaryValueAtNodeFunc,
		operatorName,
		NULL,
		2, 
		feVariableList, 
		feVariableList[0]->dim,
		feVariableList[0]->isCheckpointedAndReloaded,
		feVariableList[0]->communicator,
		feVariableList[0]->fieldVariable_Register );
}

OperatorFeVariable* OperatorFeVariable_New( 
	Name														name,
	DomainContext*											context,
	void*														feMesh,
  	void*														geometryMesh,
  	DofLayout*												dofLayout,                                                                                 
  	void*														bcs,
  	void*														ics,
  	void*														linkedDofInfo,
  	void*														templateFeVariable,    
	FeVariable_InterpolateWithinElementFunction*	interpolateWithinElement,
	FeVariable_GetValueAtNodeFunction*				getValueAtNode,
	Name														operatorName,
	Operator*                                    ownOperator,
	Index														feVariableCount,
	FeVariable**											feVariableList,
	Dimension_Index										dim,
	Bool														isCheckpointedAndReloaded,
	MPI_Comm													communicator,
	FieldVariable_Register*								fieldVariable_Register )
{
	OperatorFeVariable* self = _OperatorFeVariable_DefaultNew( name );

	self->isConstructed = True;
	_FieldVariable_Init( (FieldVariable*)self, context, feVariableCount, dim, isCheckpointedAndReloaded, communicator, fieldVariable_Register );                                                                                                          
   _FeVariable_Init( (FeVariable*)self, feMesh, geometryMesh, dofLayout, bcs, ics, linkedDofInfo, templateFeVariable, False, False );
	_OperatorFeVariable_Init( self, operatorName, feVariableCount, feVariableList, ownOperator );

	return self;
}

void* _OperatorFeVariable_DefaultNew( Name name ) {
	return _OperatorFeVariable_New( 
		sizeof(OperatorFeVariable), 
		OperatorFeVariable_Type, 
		_OperatorFeVariable_Delete, 
		_OperatorFeVariable_Print,
		_OperatorFeVariable_Copy, 
		(Stg_Component_DefaultConstructorFunction*)_OperatorFeVariable_DefaultNew,
		_OperatorFeVariable_AssignFromXML,
		_OperatorFeVariable_Build, 
		_OperatorFeVariable_Initialise, 
		_OperatorFeVariable_Execute,
		_OperatorFeVariable_Destroy,
		name,
		NON_GLOBAL,
		_OperatorFeVariable_InterpolateValueAt, /* _interpolateValueAt */
		_FeVariable_GetMinGlobalFieldMagnitude, /* _getMinGlobalFieldMagnitude */
		_FeVariable_GetMaxGlobalFieldMagnitude, /* _getMaxGlobalFieldMagnitude */
		_FeVariable_GetMinAndMaxLocalCoords, /* _getMinAndMaxLocalCoords */
		_FeVariable_GetMinAndMaxGlobalCoords, /* _getMinAndMaxGlobalCoords */
		0, /* fieldComponentCount */
		0, /* dim */
		False, /* isCheckpointedAndReloaded */
		MPI_COMM_WORLD, /* communicator */
		NULL, /* fieldVariable_Register */
		_OperatorFeVariable_InterpolateWithinElement, /*_interpolateWithinElement */
		_OperatorFeVariable_GetValueAtNode, /* _getValueAtNode */
		_OperatorFeVariable_SyncShadowValues, /* _syncShadowValues */
		NULL, /* feMesh */
		NULL, /* geometryMesh */
		NULL, /* bcs */
		NULL, /* ics */
		NULL, /* linkedDofInfo */
		NULL, /* templateFeVariable */
		NULL, /* dofLayout */
		False, /* referenceSoulution */
		False, /* loadReferenceEachTimestep */
		NULL, /* operatorName */
		NULL, /* ownOperator */
		0, /* feVariableCount */
		NULL ); /* feVariableList */
}

OperatorFeVariable* _OperatorFeVariable_New( OPERATORFEVARIABLE_DEFARGS ) {
	OperatorFeVariable* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(OperatorFeVariable) );
	self = (OperatorFeVariable*) _FeVariable_New( FEVARIABLE_PASSARGS );

	return self;
}

void _OperatorFeVariable_Init( void* oFeVar, Name operatorName, Index feVariableCount, FeVariable** feVariableList, Operator* ownOperator ) {
	OperatorFeVariable*	self = (OperatorFeVariable*) oFeVar;
	FeVariable*				feVariable;
	Index						feVariable_I;
	Stream*					errorStream = Journal_Register( Error_Type, self->type );

	/* Assign values to object */
	self->feVariableCount = feVariableCount;
	self->operatorName =operatorName;
	self->_operator = ownOperator;

	/* Copy field variable list */
	self->feVariableList = Memory_Alloc_Array( FeVariable*, feVariableCount, "Array of Field Variables" );
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

void _OperatorFeVariable_AssignFromXML( void* feVariable, Stg_ComponentFactory* cf, void* data ) {
	OperatorFeVariable*     self       = (OperatorFeVariable*) feVariable;
	Dictionary*             dictionary = Dictionary_GetDictionary( cf->componentDict, self->name );
	Dictionary_Entry_Value* list;
	Index                   feVariableCount = 0;
	Index                   feVariable_I;
	FieldVariable_Register* fV_Register;     
	Name                    feVariableName;
	Name                    operatorName;
	FeVariable**            feVariableList;
	
	/* Construct Parent */
	_FieldVariable_AssignFromXML( self, cf, data );

	fV_Register = self->context->fieldVariable_Register;

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
		feVariableList[0]->dofLayout, NULL, NULL, NULL, NULL, False, False );
	_OperatorFeVariable_Init( self, operatorName, feVariableCount, feVariableList, NULL );

	Memory_Free( feVariableList );
}

void _OperatorFeVariable_Build( void* feVariable, void* data ) {
	OperatorFeVariable* self = (OperatorFeVariable*) feVariable;
	Index                  feVariable_I;
	Stream*                     errorStream       = Journal_Register( Error_Type, self->type );

        {
          int dim;
          int numNodes;

          Stg_Component_Build(self->feMesh, NULL, False);

          dim = Mesh_GetDimSize(self->feMesh);
          numNodes = FeMesh_GetElementNodeSize(self->feMesh, 0);
          self->GNx = Memory_Alloc_2DArray( double, dim, numNodes, "Global Shape Function Derivatives" );
        }

	for ( feVariable_I = 0 ; feVariable_I < self->feVariableCount ; feVariable_I++ ) 
		Stg_Component_Build( self->feVariableList[ feVariable_I ] , data, False );

	if ( !self->_operator){
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
	} else {
			self->useGradient = False;
			self->fieldComponentCount = self->_operator->resultDofs; /* reset this value from that which is from operator */
	}		
		
	_OperatorFeVariable_SetFunctions( self );
}

void _OperatorFeVariable_Initialise( void* feVariable, void* data ) {
	OperatorFeVariable*     self = (OperatorFeVariable*) feVariable;
	DomainContext*          context = self->context;
	Index                   feVariable_I;
	Dictionary_Entry_Value* feVarsList = NULL;

	for ( feVariable_I = 0 ; feVariable_I < self->feVariableCount ; feVariable_I++ ) 
		Stg_Component_Initialise( self->feVariableList[ feVariable_I ] , data, False );

   /* also include check to see if this fevariable should be checkpointed, just incase it didn't go through the fieldvariable construct phase */ 
   feVarsList = Dictionary_Get( context->dictionary, "fieldVariablesToCheckpoint" );
   if ( NULL == feVarsList ) {
      feVarsList = Dictionary_Get( context->dictionary, "FieldVariablesToCheckpoint" );
   }
   if (feVarsList != NULL ) {
      Index                    listLength = Dictionary_Entry_Value_GetCount( feVarsList );
      Index                    var_I = 0;
      Dictionary_Entry_Value*  feVarDictValue = NULL;
      char*                    fieldVariableName;
   
      for ( var_I = 0; var_I < listLength; var_I++ ) {
         feVarDictValue = Dictionary_Entry_Value_GetElement( feVarsList, var_I );
         fieldVariableName = Dictionary_Entry_Value_AsString( feVarDictValue ); 
         if ( 0 == strcmp( self->name, fieldVariableName ) ) {
            self->isCheckpointedAndReloaded = True;
            break;
         }
      }
   }

   feVarsList = NULL;
   /** also include check to see if this fevariable should be saved for analysis purposes */ 
   feVarsList = Dictionary_Get( context->dictionary, "fieldVariablesToSave" );
   if ( NULL == feVarsList ) {
      feVarsList = Dictionary_Get( context->dictionary, "FieldVariablesToSave" );
   }
   if (feVarsList != NULL ) {
      Index                    listLength = Dictionary_Entry_Value_GetCount( feVarsList );
      Index                    var_I = 0;
      Dictionary_Entry_Value*  feVarDictValue = NULL;
      char*                    fieldVariableName;
   
      for ( var_I = 0; var_I < listLength; var_I++ ) {
         feVarDictValue = Dictionary_Entry_Value_GetElement( feVarsList, var_I );
         fieldVariableName = Dictionary_Entry_Value_AsString( feVarDictValue ); 
         if ( 0 == strcmp( self->name, fieldVariableName ) ) {
            self->isSavedData = True;
            break;
         }
      }
   }
}

void _OperatorFeVariable_Execute( void* feVariable, void* data ) {}

void _OperatorFeVariable_Destroy( void* feVariable, void* data ) {
	OperatorFeVariable* self = (OperatorFeVariable*) feVariable;

	Memory_Free( self->feVariableList );

	_FeVariable_Destroy( self, data );
}

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



void OperatorFeVariable_UnaryInterpolationFunc( void* feVariable, Element_DomainIndex dElement_I, Coord localCoord, double* value ) {
	OperatorFeVariable* self            = (OperatorFeVariable*) feVariable;
	FeVariable*         field0          = self->feVariableList[0];
	double              fieldValue0[ MAX_FIELD_COMPONENTS ]; 

	field0->_interpolateWithinElement( field0, dElement_I, localCoord, fieldValue0 );
	Operator_CarryOutUnaryOperation( self->_operator, fieldValue0, value );

#ifdef DEBUG
	if ( Stream_IsPrintableLevel( self->debug, 3 ) ) {
		Dimension_Index dim_I;
		Dof_Index       dof_I;
		Journal_DPrintf( self->debug, "%s Unary Op: El %d, xi = ( ", self->name, dElement_I );
		for ( dim_I = 0 ; dim_I < self->dim ; dim_I++ ) 
			Journal_DPrintf( self->debug, "%g ", localCoord[ dim_I ] );

		/* Field 0 */
		Journal_DPrintf( self->debug, ") - %s = [ ", field0->name );
		for ( dof_I = 0 ; dof_I < field0->fieldComponentCount ; dof_I++ ) 
			Journal_DPrintf( self->debug, "%g ", fieldValue0[ dof_I ] );

		/* Result */
		Journal_DPrintf( self->debug, "] - Result = [ " );
		for ( dof_I = 0 ; dof_I < self->fieldComponentCount ; dof_I++ ) 
			Journal_DPrintf( self->debug, "%g ", value[ dof_I ] );
		Journal_DPrintf( self->debug, "]\n" );
	}
#endif
	
}

void OperatorFeVariable_BinaryInterpolationFunc( void* feVariable, Element_DomainIndex dElement_I, Coord localCoord, double* value ) {
	OperatorFeVariable* self            = (OperatorFeVariable*) feVariable;
	FeVariable*         field0          = self->feVariableList[0];
	FeVariable*         field1          = self->feVariableList[1];
	FeMesh*             mesh            = self->feMesh;
	double              fieldValue0[ MAX_FIELD_COMPONENTS ]; 
	double              fieldValue1[ MAX_FIELD_COMPONENTS ]; 
	
	FeVariable_InterpolateFromMeshLocalCoord( field0, mesh, dElement_I, localCoord, fieldValue0 );
	FeVariable_InterpolateFromMeshLocalCoord( field1, mesh, dElement_I, localCoord, fieldValue1 );
	
	Operator_CarryOutBinaryOperation( self->_operator, fieldValue0, fieldValue1, value );

#ifdef DEBUG
	if ( Stream_IsPrintableLevel( self->debug, 3 ) ) {
		Dimension_Index dim_I;
		Dof_Index       dof_I;
		Journal_DPrintf( self->debug, "%s Binary Op: El %d, xi = ( ", self->name, dElement_I );
		for ( dim_I = 0 ; dim_I < self->dim ; dim_I++ ) 
			Journal_DPrintf( self->debug, "%g ", localCoord[ dim_I ] );

		/* Field 0 */
		Journal_DPrintf( self->debug, ") - %s = [ ", field0->name );
		for ( dof_I = 0 ; dof_I < field0->fieldComponentCount ; dof_I++ ) 
			Journal_DPrintf( self->debug, "%g ", fieldValue0[ dof_I ] );

		/* Field 1 */
		Journal_DPrintf( self->debug, "] - %s = [ ", field1->name );
		for ( dof_I = 0 ; dof_I < field1->fieldComponentCount ; dof_I++ ) 
			Journal_DPrintf( self->debug, "%g ", fieldValue1[ dof_I ] );

		/* Result */
		Journal_DPrintf( self->debug, "] - Result = [ " );
		for ( dof_I = 0 ; dof_I < self->fieldComponentCount ; dof_I++ ) 
			Journal_DPrintf( self->debug, "%g ", value[ dof_I ] );
		Journal_DPrintf( self->debug, "]\n" );
	}
#endif
}

void OperatorFeVariable_GradientInterpolationFunc( void* feVariable, Element_DomainIndex dElement_I, Coord localCoord, double* value ) {
	OperatorFeVariable* self            = (OperatorFeVariable*) feVariable;
	FeVariable*         field0          = self->feVariableList[0];
	
	FeVariable_InterpolateDerivativesToElLocalCoord( field0, dElement_I, localCoord, value );
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
	double*		    	coord;
	Element_DomainIndex	field0Element;
	Element_DomainIndex	field1Element;
	Coord			field0LocalCoord;
	Coord			field1LocalCoord;
	Node_LocalIndex		field0NearestNode;
	Node_LocalIndex		field1NearestNode;

	if( field0->feMesh == self->feMesh && field1->feMesh == self->feMesh ) {
		FeVariable_GetValueAtNode( field0, dNode_I, fieldValue0 );
		FeVariable_GetValueAtNode( field1, dNode_I, fieldValue1 );
	}
	/* field variables are using different meshes, so interpolate the values of each based on the global coord */
	else {  
		coord = Mesh_GetVertex( self->feMesh, dNode_I );
		assert( Mesh_SearchElements( field0->feMesh, coord, &field0Element ) );
		FeMesh_CoordGlobalToLocal( field0->feMesh, field0Element, coord, field0LocalCoord );
		FeVariable_InterpolateWithinElement( field0, field0Element, field0LocalCoord, fieldValue0 );
		
		assert( Mesh_SearchElements( field1->feMesh, coord, &field1Element ) ); 
		FeMesh_CoordGlobalToLocal( field1->feMesh, field1Element, coord, field1LocalCoord );
		FeVariable_InterpolateWithinElement( field1, field1Element, field1LocalCoord, fieldValue1 );
	}
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
