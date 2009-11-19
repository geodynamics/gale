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
** $Id: AnalyticSolution.c 1137 2008-05-23 05:57:48Z RobertTurnbull $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "types.h"
#include "AnalyticSolution.h"
#include "FeVariable.h"
#include "OperatorFeVariable.h"
#include "FeMesh.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

const Type AnalyticSolution_Type = "AnalyticSolution";
/* Singleton for analytical solution */
AnalyticSolution* mySingleton = NULL;

void* _AnalyticSolution_DefaultNew( Name name ) {
	return _AnalyticSolution_New(
		sizeof(AnalyticSolution),
		AnalyticSolution_Type,
		_AnalyticSolution_Delete, 
		_AnalyticSolution_Print,
		_AnalyticSolution_Copy,
		_AnalyticSolution_DefaultNew,
		_AnalyticSolution_AssignFromXML,
		_AnalyticSolution_Build,
		_AnalyticSolution_Initialise,
		_AnalyticSolution_Execute, 
		_AnalyticSolution_Destroy,
		name );
}

AnalyticSolution* _AnalyticSolution_New( 
		SizeT                                       _sizeOfSelf,
		Type                                        type,
		Stg_Class_DeleteFunction*                   _delete,
		Stg_Class_PrintFunction*                    _print,
		Stg_Class_CopyFunction*                     _copy, 
		Stg_Component_DefaultConstructorFunction*   _defaultConstructor,
		Stg_Component_ConstructFunction*            _construct,
		Stg_Component_BuildFunction*                _build,
		Stg_Component_InitialiseFunction*           _initialise,
		Stg_Component_ExecuteFunction*              _execute,
		Stg_Component_DestroyFunction*              _destroy,
		Name                                        name )
{
	AnalyticSolution*			self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(AnalyticSolution) );
	/* Construct using parent */
	self = (AnalyticSolution*)_Stg_Component_New( 
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
			NON_GLOBAL );

	/* Assign singleton ptr */
	mySingleton = self;
	return self;
}

void _AnalyticSolution_Init( AnalyticSolution* self, Swarm* integrationSwarm, LiveComponentRegister* LC_Register, AbstractContext* context, Bool verboseMode ) {
	self->LC_Register = LC_Register;
	self->integrationSwarm = integrationSwarm;
	self->context = (DomainContext*)context;

	/* Initialise AnalyticFunctions */
	self->_getAnalyticVelocity = NULL;
	self->_getAnalyticPressure = NULL;
	self->_getAnalyticTotalStress = NULL;
	self->_getAnalyticStrainRate = NULL;

	/* Create Lists */
	self->feVariableList             = Stg_ObjectList_New();
	self->analyticFeVariableList     = Stg_ObjectList_New();
	self->analyticFeVariableFuncList = Stg_ObjectList_New();
	self->errorMagnitudeFieldList    = Stg_ObjectList_New();
	self->relativeErrorMagnitudeFieldList    = Stg_ObjectList_New();
	self->streamList                 = Stg_ObjectList_New();

	/* Add functions to entry points */
	EntryPoint_AppendClassHook( Context_GetEntryPoint( context, AbstractContext_EP_UpdateClass ), 
			self->type, AnalyticSolution_Update, self->type, self );
	EP_AppendClassHook( Context_GetEntryPoint( context, AbstractContext_EP_DumpClass ), 
			AnalyticSolution_TestAll, self );

	if ( verboseMode ) {
		Stream* infoStream = Journal_MyStream( Info_Type, self );
		Stream_SetLevel( infoStream, 2 );
	}
}

void _AnalyticSolution_Delete( void* analyticSolution ) {
	AnalyticSolution* self = (AnalyticSolution*)analyticSolution;
	
	/* Stg_Class_Delete parent*/
	_Stg_Component_Delete( self );
}

void _AnalyticSolution_Print( void* analyticSolution, Stream* stream ) {
	AnalyticSolution* self = (AnalyticSolution*)analyticSolution;
	
	/* Print parent */
	_Stg_Component_Print( self, stream );
	
}

void* _AnalyticSolution_Copy( void* analyticSolution, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	abort();
	return NULL;
}

void _AnalyticSolution_AssignFromXML( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) 
{
	AnalyticSolution* self = (AnalyticSolution*)analyticSolution;
	AbstractContext*  context;
	Swarm*            integrationSwarm;
	Bool              verboseMode;

	self->context = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Context", DomainContext, False, data );
	if( !self->context )
		self->context = Stg_ComponentFactory_ConstructByName( cf, "context", DomainContext, True, data );

	integrationSwarm = Stg_ComponentFactory_ConstructByName( cf, "gaussSwarm", Swarm, True, data ); 
	verboseMode = Stg_ComponentFactory_GetRootDictBool( cf, "analyticSolutionVerbose", False );

	_AnalyticSolution_Init( self, integrationSwarm, cf->LCRegister, context, verboseMode );
}

void _AnalyticSolution_Build( void* analyticSolution, void* data ) {
	AnalyticSolution* self = (AnalyticSolution*) analyticSolution;
	Index             analyticFeVariableCount;
	Index             analyticFeVariable_I;

	/* Build all the analytic fields registered to this AnalyticSolution */
	AnalyticSolution_BuildAllAnalyticFields( self, data );
	
	analyticFeVariableCount = Stg_ObjectList_Count( self->analyticFeVariableList );
	assert( analyticFeVariableCount == Stg_ObjectList_Count( self->analyticFeVariableFuncList ) );

	for ( analyticFeVariable_I = 0 ; analyticFeVariable_I < analyticFeVariableCount ; analyticFeVariable_I++ ) {
		Stg_Component_Build( Stg_ObjectList_At( self->analyticFeVariableList, analyticFeVariable_I ), data, False ) ;
		Stg_Component_Build( Stg_ObjectList_At( self->errorMagnitudeFieldList, analyticFeVariable_I ), data, False ) ;
		Stg_Component_Build( Stg_ObjectList_At( self->relativeErrorMagnitudeFieldList, analyticFeVariable_I ), data, False ) ;
	}
}

void _AnalyticSolution_Initialise( void* analyticSolution, void* data ) {
	AnalyticSolution* self = (AnalyticSolution*) analyticSolution;
	Index             analyticFeVariableCount = Stg_ObjectList_Count( self->analyticFeVariableList );
	Index             analyticFeVariable_I;

	assert( analyticFeVariableCount == Stg_ObjectList_Count( self->analyticFeVariableFuncList ) );

	/* Assign values to all analytic fields */
	for ( analyticFeVariable_I = 0 ; analyticFeVariable_I < analyticFeVariableCount ; analyticFeVariable_I++ ) {
		Stg_Component_Initialise( Stg_ObjectList_At( self->feVariableList, analyticFeVariable_I ), data, False ) ;
		Stg_Component_Initialise( Stg_ObjectList_At( self->analyticFeVariableList, analyticFeVariable_I ), data, False ) ;
		Stg_Component_Initialise( Stg_ObjectList_At( self->errorMagnitudeFieldList, analyticFeVariable_I ), data, False ) ;
		Stg_Component_Initialise( Stg_ObjectList_At( self->relativeErrorMagnitudeFieldList, analyticFeVariable_I ), data, False ) ;

		AnalyticSolution_PutAnalyticSolutionOntoNodes( self, analyticFeVariable_I );
	}
}

/* This function is called when the 'Update' phase happens */
void AnalyticSolution_Update( void* analyticSolution ) {
	AnalyticSolution* self = (AnalyticSolution*) analyticSolution;

	self->_initialise( self, NULL );
}

void _AnalyticSolution_Execute( void* analyticSolution, void* data ) {
}

void _AnalyticSolution_Destroy( void* analyticSolution, void* data ) {
	AnalyticSolution* self = (AnalyticSolution*)analyticSolution;

	Stg_Class_Delete( self->feVariableList );
	Stg_Class_Delete( self->analyticFeVariableList );
	Stg_Class_Delete( self->analyticFeVariableFuncList );
	Stg_Class_Delete( self->errorMagnitudeFieldList );
	Stg_Class_Delete( self->relativeErrorMagnitudeFieldList );
	Stg_Class_Delete( self->streamList );

	if ( self->toleranceList )
		Memory_Free( self->toleranceList );

	Stg_Component_Destroy( self, data, False );
}

void AnalyticSolution_PutAnalyticSolutionOntoNodes( void* analyticSolution, Index analyticFeVariable_I ) {
	AnalyticSolution*                            self = (AnalyticSolution*) analyticSolution;
	AnalyticSolution_SolutionFunction*           solutionFunction;
	FeVariable*                                  analyticFeVariable;
	double*                                      coord;
	Dof_Index                                    dofAtEachNodeCount;
	Node_DomainIndex                             dNode_I;
	double*                                      value;
	Stream*                                      infoStream = Journal_MyStream( Info_Type, self );
	FeVariable*                                  feVariable;
	FeMesh*						mesh;

	/* Do some error checking */
	assert( Stg_ObjectList_Count( self->analyticFeVariableList ) >= analyticFeVariable_I );

	/* Grab pointers */
	analyticFeVariable = 
		Stg_CheckType( Stg_ObjectList_At( self->analyticFeVariableList, analyticFeVariable_I ), FeVariable );
	mesh = analyticFeVariable->feMesh;
	solutionFunction   = (AnalyticSolution_SolutionFunction*) 
		Stg_ObjectList_ObjectAt( self->analyticFeVariableFuncList, analyticFeVariable_I );
	feVariable = AnalyticSolution_GetFeVariableFromAnalyticFeVariable( self, analyticFeVariable );

	/* Get number of degrees of freedom at each node (assuming they are the same) */
	dofAtEachNodeCount = analyticFeVariable->fieldComponentCount;
	value = Memory_Alloc_Array( double, dofAtEachNodeCount, "value" );

	/* Loop over all the nodes - applying the analytic solution */
	for ( dNode_I = 0 ; dNode_I < Mesh_GetDomainSize( mesh, MT_VERTEX ); dNode_I++ ) {
		coord = Mesh_GetVertex( mesh, dNode_I );

		/* Calculate value at node */
		memset( value, 0, dofAtEachNodeCount * sizeof(double) );
		solutionFunction( self, analyticFeVariable, coord, value );

		/* Put value on node */
		FeVariable_SetValueAtNode( analyticFeVariable, dNode_I, value );

		/* Print verbose information */
		if ( Stream_IsPrintableLevel( infoStream, 2 ) ) {
			Dimension_Index dim = analyticFeVariable->dim;
			Dimension_Index dim_I;
			Dof_Index       dof_I;

			/* Print Coord */
			Journal_Printf( infoStream, "Node = %u Coord ( ", dNode_I );
			for ( dim_I = 0 ; dim_I < dim - 1 ; dim_I++ ) {
				Journal_Printf( infoStream, "%g, ", coord[ dim_I ] );
			}
			Journal_Printf( infoStream, "%g ) ", coord[ dim_I ] );

			/* Print Analytic Values */
			Journal_Printf( infoStream, "%s ( ", analyticFeVariable->name );
			for ( dof_I = 0 ; dof_I < dofAtEachNodeCount - 1 ; dof_I++ ) {
				Journal_Printf( infoStream, "%g, ", value[ dof_I ] );
			}
			Journal_Printf( infoStream, "%g ) ", value[ dof_I ] );

			/* Print FeVariable Values */
			memset( value, 0, dofAtEachNodeCount * sizeof(double) );
			FeVariable_GetValueAtNode( feVariable, dNode_I, value );
			Journal_Printf( infoStream, "%s ( ", feVariable->name );
			for ( dof_I = 0 ; dof_I < dofAtEachNodeCount - 1 ; dof_I++ ) {
				Journal_Printf( infoStream, "%g, ", value[ dof_I ] );
			}
			Journal_Printf( infoStream, "%g )\n", value[ dof_I ] );			
		}
	}

	Memory_Free( value );
}

void AnalyticSolution_Test( void* analyticSolution, Index analyticFeVariable_I ) {
	AnalyticSolution*                            self = (AnalyticSolution*) analyticSolution;
	FeVariable*                                  errorMagnitudeField;
	double                                       result;
	Stream*                                      stream;
	double                                       tolerance;
	Stream*                                      infoStream = Journal_MyStream( Info_Type, self );

	/* Do some error checking */
	assert( Stg_ObjectList_Count( self->analyticFeVariableList ) >= analyticFeVariable_I );

	/* Grab pointers */
	errorMagnitudeField = 
		Stg_CheckType( Stg_ObjectList_At( self->errorMagnitudeFieldList, analyticFeVariable_I ), FeVariable );

	stream = (Stream*) Stg_ObjectList_ObjectAt( self->streamList, analyticFeVariable_I );
	tolerance = self->toleranceList[ analyticFeVariable_I ];

	result = FeVariable_Integrate( errorMagnitudeField, self->integrationSwarm );

	Journal_Printf( stream, "Timestep %u: Total integrated value of '%s' is %s a tolerance %.5g.\n", 
				self->context->timeStep, errorMagnitudeField->name, result <= tolerance ? "within" : "outside", tolerance );
	Stream_Flush( stream );

	Journal_Printf( infoStream, "Timestep %u: Total integrated value of '%s' is %g\n", 
			self->context->timeStep, errorMagnitudeField->name, result );
}
	

void AnalyticSolution_TestAll( void* analyticSolution, void* data ) {
	AnalyticSolution* self = (AnalyticSolution*) analyticSolution;
	Index             analyticFeVariableCount = Stg_ObjectList_Count( self->analyticFeVariableList );
	Index             analyticFeVariable_I;

	assert( analyticFeVariableCount == Stg_ObjectList_Count( self->analyticFeVariableFuncList ) );

	/* Assign values to all analytic fields */
	for ( analyticFeVariable_I = 0 ; analyticFeVariable_I < analyticFeVariableCount ; analyticFeVariable_I++ ) {
		AnalyticSolution_Test( self, analyticFeVariable_I );
	}
}

void AnalyticSolution_RegisterFeVariableWithAnalyticFunction( void* analyticSolution, FeVariable* feVariable, AnalyticSolution_SolutionFunction* solutionFunction) {
	AnalyticSolution* self    = (AnalyticSolution*) analyticSolution;
	char*             tmpName = Stg_Object_AppendSuffix( feVariable, "Analytic" );

	/* Add feVariable to list */
	Stg_ObjectList_Append( self->feVariableList, feVariable );
	/* Add function to list */
	Stg_ObjectList_GlobalPointerAppend( self->analyticFeVariableFuncList, solutionFunction, tmpName );
	Memory_Free( tmpName );	
}


FeVariable* AnalyticSolution_RegisterFeVariableFromCF( void* analyticSolution, char* fieldName, AnalyticSolution_SolutionFunction* solutionFunction, Stg_ComponentFactory* cf, Bool isEssential, void* data ) {
	AnalyticSolution* self    = (AnalyticSolution*) analyticSolution;
	FeVariable*       field;

	field = Stg_ComponentFactory_ConstructByName( cf, fieldName, FeVariable, isEssential, data ); 
	if ( field )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, field, solutionFunction );

	return field;
}


void AnalyticSolution_BuildAllAnalyticFields( void* analyticSolution, void* data ) {
	AnalyticSolution* self       = (AnalyticSolution*) analyticSolution;
	FeVariable*       feVariable = NULL;
	Stream*           errStream  = Journal_Register( Error_Type, "AnalyticSolution" );
	unsigned          feVar_I, feVarCount;

	feVarCount = Stg_ObjectList_Count( self->feVariableList );

	for( feVar_I = 0 ; feVar_I < feVarCount ; feVar_I++ ) {
		feVariable = (FeVariable*)Stg_ObjectList_At( self->feVariableList, feVar_I );
		
		/* Check to see whether this field has already been added */
		

		/* Build the FeVariable here ensuring the analytic "copy" of it has valid values */
		Stg_Component_Build( feVariable, data, False ) ;
		
		if( feVariable->fieldComponentCount == 1 ) {
			if( NULL == AnalyticSolution_CreateAnalyticField( self, feVariable ) ) {

				Journal_Firewall( 0 , errStream,
				"Error in function %s: Error in building analyticSolution for the feVariable %s,\n",
				__func__ , feVariable );
			}
		} else {
			if( NULL == AnalyticSolution_CreateAnalyticSymmetricTensorField( self, feVariable ) ) {

				Journal_Firewall( 0 , errStream,
				"Error in function %s: Error in building analyticSolution for the feVariable %s,\n",
				__func__ , feVariable );
			}
		}

	}


}

FeVariable* AnalyticSolution_CreateAnalyticField( void* analyticSolution, FeVariable* feVariable ) {
	AnalyticSolution*                            self = (AnalyticSolution*) analyticSolution;
	Name                                         tmpName;
	FeVariable*                                  analyticFeVariable;
	Variable*                                    dataVariable;
	DofLayout*                                   dofLayout;
	Variable_Register*                           variable_Register = self->context->variable_Register;
	Dof_Index                                    componentsCount   = feVariable->fieldComponentCount;
	Name                                         variableName[9];
	Variable_Index                               variable_I;
	Node_DomainIndex                             node_I;
	Bool                                         scalar             = ( componentsCount == 1 );
	OperatorFeVariable*                          analyticMagField;
	OperatorFeVariable*                          errorField;
	OperatorFeVariable*                          errorMagnitudeField;
	OperatorFeVariable*                          relativeErrorMagnitudeField;
	Stream*                                      stream;
	Index                                        count;

	Stg_Component_Build( feVariable->feMesh, NULL, False );

	/* Create new data Variable */
	tmpName = Stg_Object_AppendSuffix( feVariable, "Analytic-DataVariable" );
	if ( scalar ) {
		Sync*	sync;

		sync = Mesh_GetSync( feVariable->feMesh, MT_VERTEX );
		dataVariable = Variable_NewScalar( 	
			tmpName,
			Variable_DataType_Double, 
			(unsigned*)&sync->nDomains, 
			NULL,
			(void**)NULL, 
			variable_Register );
	}
	else {
		Sync*		sync;
		unsigned	c_i;

		sync = Mesh_GetSync( feVariable->feMesh, MT_VERTEX );

		/* Create names of variables */
		assert( componentsCount <= 9 );
		for ( variable_I = 0 ; variable_I < componentsCount ; variable_I++ ) {
			Stg_asprintf( &variableName[ variable_I ], "%s-Analytic-ComponentVariable%d", feVariable->name, variable_I );
		}
		dataVariable = Variable_NewVector( 	
				tmpName,
				Variable_DataType_Double, 
				componentsCount,
				(unsigned*)&sync->nDomains, 
				NULL, 
				(void**)NULL, 
				variable_Register,
				variableName[0],
				variableName[1],
				variableName[2],
				variableName[3],
				variableName[4],
				variableName[5],
				variableName[6],
				variableName[7],
				variableName[8] );
		for( c_i = 0; c_i < dataVariable->dataTypeCounts[0]; c_i++ )
			dataVariable->components[c_i]->allocateSelf = True;
	}
	Memory_Free( tmpName );
	dataVariable->allocateSelf = True;
	
	/* Create new dof layout */
	tmpName = Stg_Object_AppendSuffix( feVariable, "Analytic-DofLayout" );
	dofLayout = DofLayout_New( tmpName, variable_Register, Mesh_GetDomainSize( feVariable->feMesh, MT_VERTEX ), NULL );
	if ( scalar ) {
		DofLayout_AddAllFromVariableArray( dofLayout, 1, &dataVariable );
	}
	else {
		for ( variable_I = 0 ; variable_I < componentsCount ; variable_I++ ) {

			/* We have to set the array ptr ptr for these guys manually - this should be fixed */
			Variable* variable = Variable_Register_GetByName( variable_Register, variableName[ variable_I ] );
			variable->arrayPtrPtr = &dataVariable->arrayPtr;

			/* Assign variable to each node */
			for( node_I = 0; node_I < Mesh_GetDomainSize( feVariable->feMesh, MT_VERTEX ); node_I++ ) {
				DofLayout_AddDof_ByVarName( dofLayout, variableName[variable_I], node_I );
			}
			/* Free Name */
			Memory_Free( variableName[ variable_I ] );
		}
	}
	Memory_Free( tmpName );

	/* Create new FeVariable */
	tmpName = Stg_Object_AppendSuffix( feVariable, "Analytic" );
	analyticFeVariable = FeVariable_New( tmpName, self->context, feVariable->feMesh, feVariable->geometryMesh, dofLayout,
		NULL, NULL, NULL, feVariable->dim, feVariable->isCheckpointedAndReloaded, 
		False, False,
		feVariable->fieldVariable_Register );

	/* Add new analyticFeVariable to list */
	Stg_ObjectList_Append( self->analyticFeVariableList, analyticFeVariable );

	/* Create Magnitude Field */
	tmpName = Stg_Object_AppendSuffix( analyticFeVariable, "Magnitude" );
	analyticMagField = OperatorFeVariable_NewUnary( tmpName, self->context, analyticFeVariable, "Magnitude" );
	Memory_Free( tmpName );

	/* Create Error field - The the calculated field minus the analytic field */
	tmpName = Stg_Object_AppendSuffix( feVariable, "ErrorField" );
	errorField = OperatorFeVariable_NewBinary( tmpName, self->context, feVariable, analyticFeVariable, "Subtraction" );
	Memory_Free( tmpName );
	
	/* Create Error magnitude field */
	tmpName = Stg_Object_AppendSuffix( feVariable, "ErrorMagnitudeField" );
	errorMagnitudeField = OperatorFeVariable_NewUnary( tmpName, self->context, errorField, "Magnitude" );
	Memory_Free( tmpName );
	Stg_ObjectList_Append( self->errorMagnitudeFieldList, errorMagnitudeField ); /* Add it to list */

	/* Create Relative Error magnitude field - The magnitude of relative error */
	tmpName = Stg_Object_AppendSuffix( feVariable, "RelativeErrorMagnitudeField" );
	relativeErrorMagnitudeField = OperatorFeVariable_NewBinary( tmpName, self->context, errorMagnitudeField, analyticMagField, "ScalarDivision" );
	Memory_Free( tmpName );
	Stg_ObjectList_Append( self->relativeErrorMagnitudeFieldList, relativeErrorMagnitudeField ); /* Add it to list */

	/* Create Stream for field to dump error information to */
	tmpName = Stg_Object_AppendSuffix( feVariable, "ErrorFile" );
	stream = Journal_Register( Dump_Type, tmpName );
	Stg_ObjectList_GlobalPointerAppend( self->streamList, stream, tmpName );
	Stream_Enable( stream, True );
	Stream_RedirectFile_WithPrependedPath( stream, self->context->outputPath, tmpName );
	Memory_Free( tmpName );

	/* Get tolerance from dictionary */
	count = Stg_ObjectList_Count( self->analyticFeVariableList );
	self->toleranceList = Memory_Realloc_Array( self->toleranceList, double, count );
	tmpName = Stg_Object_AppendSuffix( feVariable, "Tolerance" );
	self->toleranceList[ count - 1 ] = Dictionary_GetDouble_WithDefault( self->context->dictionary, tmpName, 0.0 );
	Memory_Free( tmpName );

	/* Add components to LiveComponentRegister so they will be visible to other components 
	 * and will be built, initialised and  deleted */
	LiveComponentRegister_Add( self->LC_Register, (Stg_Component*) dataVariable );
	LiveComponentRegister_Add( self->LC_Register, (Stg_Component*) dofLayout );
	LiveComponentRegister_Add( self->LC_Register, (Stg_Component*) analyticFeVariable );
	LiveComponentRegister_Add( self->LC_Register, (Stg_Component*) analyticMagField );
	LiveComponentRegister_Add( self->LC_Register, (Stg_Component*) errorField );
	LiveComponentRegister_Add( self->LC_Register, (Stg_Component*) errorMagnitudeField );
	LiveComponentRegister_Add( self->LC_Register, (Stg_Component*) relativeErrorMagnitudeField );

	return analyticFeVariable;
}

FeVariable* AnalyticSolution_CreateAnalyticVectorField( void* analyticSolution, FeVariable* vectorField, AnalyticSolution_SolutionFunction* solutionFunction ) {
	AnalyticSolution*	self = (AnalyticSolution*) analyticSolution;
	FeVariable*			analyticVectorField;

	analyticVectorField = AnalyticSolution_CreateAnalyticField( self, vectorField );

	return analyticVectorField;
}

FeVariable* AnalyticSolution_CreateAnalyticSymmetricTensorField( void* analyticSolution, FeVariable* vectorField ) {
	AnalyticSolution*		self = (AnalyticSolution*) analyticSolution;
	FeVariable*				analyticVectorField;
	OperatorFeVariable*	analyticVectorInvField;
	Name						tmpName, tmpName2;
	DofLayout*				dofLayout;

	analyticVectorField = AnalyticSolution_CreateAnalyticField( self, vectorField );

	/* Create new dof layout */
	tmpName = Stg_Object_AppendSuffix( analyticVectorField, "Analytic-DofLayout" );
	dofLayout = DofLayout_New( tmpName, self->context->variable_Register, Mesh_GetDomainSize( analyticVectorField->feMesh, MT_VERTEX ), NULL );

	/* Create Invariant Field */
	tmpName2 = Stg_Object_AppendSuffix( analyticVectorField, "Invariant" );
	analyticVectorInvField = OperatorFeVariable_NewUnary( tmpName2, self->context, analyticVectorField, "SymmetricTensor_Invariant" );

	Memory_Free( tmpName );
	Memory_Free( tmpName2 );

	LiveComponentRegister_Add( self->LC_Register, (Stg_Component*) analyticVectorInvField );

	return analyticVectorField;
}

FeVariable* AnalyticSolution_GetFeVariableFromAnalyticFeVariable( void* analyticSolution, FeVariable* analyticFeVariable ) {
	AnalyticSolution* self = (AnalyticSolution*) analyticSolution;
	Index             analyticFeVariableCount = Stg_ObjectList_Count( self->analyticFeVariableList );
	Index             analyticFeVariable_I;

	assert( analyticFeVariableCount == Stg_ObjectList_Count( self->feVariableList ) );

	for ( analyticFeVariable_I = 0 ; analyticFeVariable_I < analyticFeVariableCount ; analyticFeVariable_I++ ) {
		/* find the index of analytic feVariable and then return the corresponding feVariable */
		if ( analyticFeVariable == (FeVariable*) Stg_ObjectList_At( self->analyticFeVariableList, analyticFeVariable_I ) )
			return (FeVariable*) Stg_ObjectList_At( self->feVariableList, analyticFeVariable_I );
	}

	return NULL;
}
	
InterpolationResult AnalyticSolution_InterpolateValueFromNormalFeVariable( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* value ) {
	AnalyticSolution*                            self = (AnalyticSolution*) analyticSolution;
	FeVariable*       normalFeVariable;

	normalFeVariable = AnalyticSolution_GetFeVariableFromAnalyticFeVariable( self, analyticFeVariable );

	return FieldVariable_InterpolateValueAt( normalFeVariable, coord, value );
}

AnalyticSolution* AnalyticSolution_GetAnalyticSolution() {
	Journal_Firewall( mySingleton != NULL , Journal_Register( Error_Type, "AnalyticSolution" ),
			"Error in function %s: The Singleton Ptr is NULL, meaning the AnalyticSolution has not been created yet\n",
			__func__ );
	return mySingleton;
}


