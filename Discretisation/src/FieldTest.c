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
** $Id: FieldTest.c 1095 2008-04-03 06:29:29Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "types.h"
#include "FieldTest.h"
#include "FeVariable.h"
#include "ElementType.h"
#include "FeMesh.h"
#include "OperatorFeVariable.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

const Type FieldTest_Type = "FieldTest";

FieldTest* fieldTestSingleton = NULL;

void* _FieldTest_DefaultNew( Name name ) {
	return _FieldTest_New(
		sizeof(FieldTest),
		FieldTest_Type,
		_FieldTest_Delete, 
		_FieldTest_Print,
		_FieldTest_Copy,
		_FieldTest_DefaultNew,
		_FieldTest_Construct,
		_FieldTest_Build,
		_FieldTest_Initialise,
		_FieldTest_Execute, 
		_FieldTest_Destroy,
		name );
}

FieldTest* _FieldTest_New( 
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
	FieldTest*			self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(FieldTest) );
	/* Construct using parent */
	self = (FieldTest*)_Stg_Component_New( 
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

	self->normalise = False;
	self->epsilon = 0.0001;
	self->referenceSolnFromFile = False;

	/* Assign singleton ptr */
	fieldTestSingleton = self;
	return self;
}

void _FieldTest_Delete( void* fieldTest ) {
	FieldTest* 	self 	= (FieldTest*)fieldTest;
	
	if( self->integrationSwarm ) Stg_Class_Delete( self->integrationSwarm );

	Memory_Free( self->gAnalyticSq );
	Memory_Free( self->gErrorSq    );
	Memory_Free( self->gError      );
	Memory_Free( self->gErrorNorm  );

	Memory_Free( self->analyticSolnForFeVarKey );
	Memory_Free( self->_analyticSolutionList );
	
	/* Stg_Class_Delete parent*/
	_Stg_Component_Delete( self );
}

void _FieldTest_Print( void* fieldTest, Stream* stream ) {
	FieldTest* self = (FieldTest*)fieldTest;
	
	/* Print parent */
	_Stg_Component_Print( self, stream );
	
}

void* _FieldTest_Copy( void* fieldTest, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	abort();
	return NULL;
}

void _FieldTest_Construct( void* fieldTest, Stg_ComponentFactory* cf, void* data ) {
	FieldTest* 			self 			= (FieldTest*)fieldTest;
	Dictionary*			dict			= cf->rootDict;
	Dictionary_Entry_Value*		dictEntryVal		= Dictionary_Get( dict, "pluginData" );
	Dictionary_Entry_Value*		pluginDict		= Dictionary_Entry_Value_AsDictionary( dictEntryVal );
	Dictionary_Entry_Value*		fieldList;
	Dictionary_Entry_Value*		swarmVarList		= Dictionary_Get( dict, "NumericSwarmVariableNames" );
	FieldVariable_Register* 	fV_Register     	= Stg_ObjectList_Get( cf->registerRegister, "FieldVariable_Register" );
	FieldVariable_Register* 	sW_Register     	= Stg_ObjectList_Get( cf->registerRegister, "SwarmVariable_Register" );
	Index				feVariable_I;
	Index				swarmVar_I;
	Name                    	fieldName;
	Hook*				generateErrorFields;
	Hook*				physicsTestHook;

	fieldList = Dictionary_Get( pluginDict, "NumericFields" );
	self->fieldCount = fieldList ? Dictionary_Entry_Value_GetCount( fieldList ) : 0;

	if( self->fieldCount ) {
		self->numericFieldList   	= Memory_Alloc_Array( FeVariable*, self->fieldCount, "numeric fields" );
		self->referenceFieldList 	= Memory_Alloc_Array( FeVariable*, self->fieldCount, "reference fields" );
		self->errorFieldList     	= Memory_Alloc_Array( FeVariable*, self->fieldCount, "error fields" );
		self->referenceMagFieldList	= Memory_Alloc_Array( OperatorFeVariable*, self->fieldCount, "reference field magnitudes" );
		self->errorMagFieldList		= Memory_Alloc_Array( OperatorFeVariable*, self->fieldCount, "error field magnitudes" );
	
		for( feVariable_I = 0; feVariable_I < self->fieldCount; feVariable_I++ ) {
			fieldName = ( fieldList ) ? 
				    Dictionary_Entry_Value_AsString( Dictionary_Entry_Value_GetElement( fieldList, feVariable_I ) ) :
				    Dictionary_GetString( pluginDict, "FeVariable" );
			
			self->numericFieldList[feVariable_I] = (FeVariable*) FieldVariable_Register_GetByName( fV_Register, fieldName );

			if( !self->numericFieldList[feVariable_I] )
				self->numericFieldList[feVariable_I] = Stg_ComponentFactory_ConstructByName( cf, fieldName, FeVariable, True, data ); 
		}
	}

	self->integrationSwarm 	= LiveComponentRegister_Get( cf->LCRegister, Dictionary_Entry_Value_AsString( Dictionary_Get( pluginDict, "IntegrationSwarm" ) ) );
	self->constantMesh     	= LiveComponentRegister_Get( cf->LCRegister, Dictionary_Entry_Value_AsString( Dictionary_Get( pluginDict, "ConstantMesh"     ) ) );
	self->elementMesh      	= LiveComponentRegister_Get( cf->LCRegister, Dictionary_Entry_Value_AsString( Dictionary_Get( pluginDict, "ElementMesh"      ) ) );

	self->swarmVarCount = swarmVarList ? Dictionary_Entry_Value_GetCount( swarmVarList ) : 0;
	if( self->swarmVarCount ) {
		self->swarmVarNameList = Memory_Alloc_Array( Name, self->swarmVarCount, "numeric swarm variable names" );
	
		for( swarmVar_I = 0; swarmVar_I < self->swarmVarCount; swarmVar_I++ ) {
			self->swarmVarNameList[swarmVar_I] = ( swarmVarList ) ? 
					Dictionary_Entry_Value_AsString( Dictionary_Entry_Value_GetElement( swarmVarList, swarmVar_I ) ) :
					Dictionary_GetString( pluginDict, "SwarmVariable" );
		}	
	}
	
	self->referenceSolnPath		= Dictionary_Entry_Value_AsString( Dictionary_Get( pluginDict, "referenceSolutionFilePath" ) );
	self->normalise			= Dictionary_Entry_Value_AsBool( Dictionary_Get( pluginDict, "normaliseByReferenceSolution" ) );
	self->epsilon			= Dictionary_Entry_Value_AsDouble( Dictionary_Get( pluginDict, "epsilon" ) );
	self->referenceSolnFromFile	= Dictionary_Entry_Value_AsBool( Dictionary_Get( pluginDict, "useReferenceSolutionFromFile" ) );
	self->context			= Stg_ComponentFactory_ConstructByName( cf, "context", DomainContext, True, data );

	/* for the physics test */
	self->expectedFileName = Dictionary_Entry_Value_AsString( Dictionary_Get( pluginDict, "expectedFileName" ) );
	self->expectedFilePath = Dictionary_Entry_Value_AsString( Dictionary_Get( pluginDict, "expectedFilePath" ) );
	self->dumpExpectedFileName = Dictionary_Entry_Value_AsString( Dictionary_Get( pluginDict, "expectedOutputFileName" ) );
	self->expectedPass     = False;

	/* set up the entry point */
	generateErrorFields = Hook_New( "Generate error fields hook", FieldTest_GenerateErrFields, self->name );
	_EntryPoint_AppendHook( Context_GetEntryPoint( self->context, AbstractContext_EP_Solve ), generateErrorFields );

	/* entry point for the fisix test func */
	if( strlen(self->expectedFileName) > 1 ) {
		physicsTestHook = Hook_New( "Physics test hook", FieldTest_EvaluatePhysicsTest, self->name );
		_EntryPoint_AppendHook( Context_GetEntryPoint( self->context, AbstractContext_EP_Solve ), physicsTestHook );
	}

	self->LCRegister = cf->LCRegister;
}

void _FieldTest_Build( void* fieldTest, void* data ) {
	FieldTest* 		self 		= (FieldTest*) fieldTest;
	Index			field_I;
	//Index			swarm_I;

	if( self->constantMesh ) Stg_Component_Build( self->constantMesh,  data, False );

	for( field_I = 0; field_I < self->fieldCount; field_I++ ) {
		FieldTest_BuildReferenceField( self, field_I );
		FieldTest_BuildErrField( self, field_I );
	
		Stg_Component_Build( self->numericFieldList[field_I], data, False );
		Stg_Component_Build( self->errorFieldList[field_I], data, False );

		Stg_Component_Build( self->referenceFieldList[field_I], data, False );
		Stg_Component_Build( self->errorFieldList[field_I], data, False );
	}

	/*for( swarm_I = 0; swarm_I < self->swarmCount; swarm_I++ ) {
		FieldTest_BuildReferenceSwarm( self->referenceMesh, self->numericSwarmList[swarm_I], 
					       self->context, &self->referenceSwarmList[swarm_I] );
		FieldTest_BuildErrSwarm( self->constantMesh, self->numericSwarmList[swarm_I], 
					 self->context, &self->referenceSwarmList[swarm_I] );
	}*/

	if( self->fieldCount ) {
		self->gAnalyticSq = Memory_Alloc_2DArray( double, self->fieldCount, 3, "global reference solution squared" );
		self->gErrorSq 	  = Memory_Alloc_2DArray( double, self->fieldCount, 3, "global L2 error squared" );
		self->gError	  = Memory_Alloc_2DArray( double, self->fieldCount, 3, "global L2 error" );
		self->gErrorNorm  = Memory_Alloc_2DArray( double, self->fieldCount, 3, "global L2 error normalised" );

		if( !self->referenceSolnFromFile ) {
			self->analyticSolnForFeVarKey = Memory_Alloc_Array( unsigned, self->fieldCount, 
								    "analytic solution index for ith feVariable" );
			self->_analyticSolutionList = Memory_Alloc_Array_Unnamed( FieldTest_AnalyticSolutionFunc*, self->fieldCount );
		}
	}
}

void _FieldTest_Initialise( void* fieldTest, void* data ) {
	FieldTest* 		self 		= (FieldTest*) fieldTest;
	Index			field_I;
	FILE*			expected_fp;
	char*			expectedFilename;
	int			num_time_steps;
	int			dof_i, dim_i;
	int			expected_i 	= 0;

	for( field_I = 0; field_I < self->fieldCount; field_I++ ) {
		Stg_Component_Initialise( self->numericFieldList[field_I], data, False );
		Stg_Component_Initialise( self->errorFieldList[field_I], data, False );
		Stg_Component_Initialise( self->referenceFieldList[field_I], data, False );
		if( self->referenceMagFieldList[field_I] )
			Stg_Component_Initialise( self->referenceMagFieldList[field_I], data, False );
		if( self->errorMagFieldList[field_I] )
			Stg_Component_Initialise( self->errorMagFieldList[field_I], data, False );
	}

	/* load the reference solution from file if req'd */
	if( self->referenceSolnFromFile ) {
		for( field_I = 0; field_I < self->fieldCount; field_I++ ) {
			FieldTest_LoadReferenceSolutionFromFile( self->referenceFieldList[field_I],
				       				 self->referenceSolnFileList[field_I],
								 self->referenceSolnPath, self->context );
		}
	}
	/* calculate the analytic solutions */
	else {
		for( field_I = 0; field_I < self->fieldCount; field_I++ )
			FieldTest_CalculateAnalyticSolutionForField( self, field_I );
	}

	if( strlen(self->expectedFileName) > 1 ) {
		expectedFilename = Memory_Alloc_Array_Unnamed( char, strlen(self->expectedFilePath) + strlen(self->expectedFileName) + 1 );
		sprintf( expectedFilename, "%s%s", self->expectedFilePath, self->expectedFileName );

		expected_fp = fopen( expectedFilename, "r" );

		fscanf( expected_fp, "%d %d", &self->expectedDofs, &num_time_steps );
	
		self->expected  = Memory_Alloc_Array_Unnamed( Event, num_time_steps + 1 );
		self->numeric   = Memory_Alloc_Array_Unnamed( Event, self->context->maxTimeSteps + 1 );
		self->tolerance = Memory_Alloc_Array_Unnamed( Event, num_time_steps + 1 );

		while ( !feof( expected_fp ) ) {
			fscanf( expected_fp, "%lf ", &self->expected[expected_i].time );
			
			for( dim_i = 0; dim_i < self->context->dim; dim_i++ )
				fscanf( expected_fp, "%lf ", &self->expected[expected_i].place[dim_i] );

			for( dof_i = 0; dof_i < self->expectedDofs; dof_i++ )
				fscanf( expected_fp, "%lf ", &self->expected[expected_i].value[dof_i] );

			fscanf( expected_fp, "%lf ", &self->tolerance[expected_i].time );
			
			for( dim_i = 0; dim_i < self->context->dim; dim_i++ )
				fscanf( expected_fp, "%lf ", &self->tolerance[expected_i].place[dim_i] );

			for( dof_i = 0; dof_i < self->expectedDofs; dof_i++ )
				fscanf( expected_fp, "%lf ", &self->tolerance[expected_i].value[dof_i] );

			expected_i++;
		}
		fclose( expected_fp );

		Memory_Free( expectedFilename );
	}
}

void FieldTest_CalculateAnalyticSolutionForField( void* fieldTest, Index field_I ) {
	FieldTest* 			self 		= (FieldTest*) fieldTest;
	FeVariable*			analyticField	= self->referenceFieldList[field_I];
	FeMesh*				analyticMesh	= analyticField->feMesh;
	FieldTest_AnalyticSolutionFunc*	analyticSolution;
	Index				lNode_I;
	Index				lMeshSize	= Mesh_GetLocalSize( analyticMesh, MT_VERTEX );
	double*				coord;
	double*				value;
	
	analyticSolution = self->_analyticSolutionList[self->analyticSolnForFeVarKey[field_I]];
	value = Memory_Alloc_Array_Unnamed( double, analyticField->fieldComponentCount );
	memset( value, 0, analyticField->fieldComponentCount * sizeof(double) );

	for( lNode_I = 0; lNode_I < lMeshSize; lNode_I++ ) {
		coord = Mesh_GetVertex( analyticMesh, lNode_I );
		analyticSolution( self, coord, value );
		FeVariable_SetValueAtNode( analyticField, lNode_I, value );
	}

	Memory_Free( value );
}
void _FieldTest_Execute( void* fieldTest, void* data ) {
	FieldTest* 		self 		= (FieldTest*) fieldTest;
}

void _FieldTest_Destroy( void* fieldTest, void* data ) {
	FieldTest* 		self 		= (FieldTest*) fieldTest;
	
	Memory_Free( self->numericFieldList );
	Memory_Free( self->referenceFieldList );
	Memory_Free( self->errorFieldList );
	Memory_Free( self->referenceMagFieldList );
	Memory_Free( self->errorMagFieldList );

	Memory_Free( self->gAnalyticSq );
	Memory_Free( self->gErrorSq );
	Memory_Free( self->gError );
	Memory_Free( self->gErrorNorm );
	if( strlen(self->expectedFileName) > 1 ) {
		Memory_Free( self->expected );
		Memory_Free( self->numeric );
		Memory_Free( self->tolerance );
	}
	if( !self->referenceSolnFromFile ) {
		Memory_Free( self->analyticSolnForFeVarKey );
		Memory_Free( self->_analyticSolutionList );
	}
}

void FieldTest_BuildReferenceField( void* fieldTest, Index field_I ) {
	FieldTest* 		self 			= (FieldTest*) fieldTest;
	FeVariable*		numericField		= self->numericFieldList[field_I];
	FeMesh*			referenceMesh		= numericField->feMesh;
	DomainContext*		context			= self->context;
	Variable_Register*	variable_Register	= context->variable_Register;
	Sync*			sync			= Mesh_GetSync( referenceMesh, MT_VERTEX );
	Name			tmpName;
	Dof_Index		componentsCount		= numericField->fieldComponentCount;
	Name			varName[9];
	unsigned		var_I;
	unsigned		node_I;
	Variable*		variable;
	Variable*		baseVariable		= NULL;
	DofLayout*		referenceDofLayout	= NULL;

	unsigned		nDomainVerts		= Mesh_GetDomainSize( referenceMesh, MT_VERTEX );
	static double*		arrayPtr;

	tmpName = Stg_Object_AppendSuffix( numericField, "ReferenceVariable" );

	if( componentsCount == 1 ) {
		arrayPtr = Memory_Alloc_Array_Unnamed( double, nDomainVerts );
		baseVariable = Variable_NewScalar( tmpName, Variable_DataType_Double, &nDomainVerts, NULL, (void**)&arrayPtr, 
						   variable_Register );
	}
	else {
		for( var_I = 0; var_I < componentsCount; var_I++ )
			Stg_asprintf( &varName[var_I], "%s-Component-%d", tmpName, var_I );

		arrayPtr = Memory_Alloc_Array_Unnamed( double, nDomainVerts * componentsCount );
		baseVariable = Variable_NewVector( tmpName, 
						   Variable_DataType_Double, 
						   componentsCount, 
						   &nDomainVerts,
						   NULL, 
						   (void**)&arrayPtr, 
						   variable_Register,
						   varName[0], varName[1], varName[2], 
						   varName[3], varName[4], varName[5], 
						   varName[6], varName[7], varName[8] );
	}
	Memory_Free( tmpName );

	tmpName = Stg_Object_AppendSuffix( numericField, "ReferenceDofLayout" );

	referenceDofLayout = DofLayout_New( tmpName, variable_Register, Mesh_GetDomainSize( referenceMesh, MT_VERTEX ), referenceMesh );
	if( componentsCount == 1 )
		DofLayout_AddAllFromVariableArray( referenceDofLayout, 1, &baseVariable );
	else {
		for( var_I = 0; var_I < componentsCount; var_I++ ) {
			variable = Variable_Register_GetByName( variable_Register, varName[var_I] );
			variable->arrayPtrPtr = &baseVariable->arrayPtr;

			for( node_I = 0; node_I < Mesh_GetDomainSize( referenceMesh, MT_VERTEX ); node_I++ )
				DofLayout_AddDof_ByVarName( referenceDofLayout, varName[var_I], node_I );

			Memory_Free( varName[var_I] );
		}

		//referenceDofLayout->baseVariables = Memory_Alloc_Array_Unnamed( Variable*, componentsCount );
		//referenceDofLayout->baseVariables[0] = baseVariable;
	}

	Stg_Component_Build( referenceDofLayout, NULL, False );
	Stg_Component_Initialise( referenceDofLayout, NULL, False );

	Memory_Free( tmpName );

	tmpName = Stg_Object_AppendSuffix( numericField, "Reference" );

	self->referenceFieldList[field_I] = FeVariable_New( tmpName, referenceMesh, referenceMesh, referenceDofLayout, NULL, NULL, NULL, 
							Mesh_GetDimSize( referenceMesh ), False, "StgFEM_Native", "StgFEM_Native",
							"", "", False, False, context->fieldVariable_Register );
	/* so that the eqnation numbers don't get built for this guy */
	self->referenceFieldList[field_I]->buildEqNums = False;

	tmpName = Stg_Object_AppendSuffix( self->referenceFieldList[field_I], "Magnitude" );
	self->referenceMagFieldList[field_I] = OperatorFeVariable_NewUnary( tmpName, self->referenceFieldList[field_I], "Magnitude" );
	Memory_Free( tmpName );
	Stg_Component_Build( self->referenceMagFieldList[field_I], context, False );
	self->referenceMagFieldList[field_I]->_operator = Operator_NewFromName( self->referenceMagFieldList[field_I]->operatorName, 
							self->referenceFieldList[field_I]->fieldComponentCount, context->dim );
	self->referenceMagFieldList[field_I]->fieldComponentCount = self->referenceMagFieldList[field_I]->_operator->resultDofs;
	_OperatorFeVariable_SetFunctions( self->referenceMagFieldList[field_I] );

	LiveComponentRegister_Add( context->CF->LCRegister, (Stg_Component*) baseVariable );
	LiveComponentRegister_Add( context->CF->LCRegister, (Stg_Component*) referenceDofLayout );
	LiveComponentRegister_Add( context->CF->LCRegister, (Stg_Component*) self->referenceFieldList[field_I] );
	LiveComponentRegister_Add( context->CF->LCRegister, (Stg_Component*) self->referenceMagFieldList[field_I] );
}

void FieldTest_BuildErrField( void* fieldTest, Index field_I ) {
	FieldTest* 		self 			= (FieldTest*) fieldTest;
	FeMesh*			constantMesh		= self->constantMesh;
	FeVariable*		numericField		= self->numericFieldList[field_I];
	DomainContext*		context			= self->context;
	Variable_Register*	variable_Register	= context->variable_Register;
	Sync*			sync			= Mesh_GetSync( constantMesh, MT_VERTEX );
	Name			tmpName;
	Dof_Index		componentsCount		= numericField->fieldComponentCount;
	Name			varName[9];
	unsigned		var_I;
	unsigned		node_I;
	Variable*		variable;
	Variable*		baseVariable		= NULL;
	DofLayout*		errorDofLayout		= NULL;
	unsigned		nDomainVerts		= Mesh_GetDomainSize( constantMesh, MT_VERTEX );
	static void*		arrayPtr;

	tmpName = Stg_Object_AppendSuffix( numericField, "ErrorVariable" );

	if( componentsCount == 1 ) {
		arrayPtr = Memory_Alloc_Array_Unnamed( double, nDomainVerts );
		baseVariable = Variable_NewScalar( tmpName, Variable_DataType_Double, &nDomainVerts, NULL, &arrayPtr, 
						   variable_Register );
	}
	else {
		for( var_I = 0; var_I < componentsCount; var_I++ )
			Stg_asprintf( &varName[var_I], "%s-Component-%d", tmpName, var_I );

		arrayPtr = Memory_Alloc_Array_Unnamed( double, nDomainVerts * componentsCount );
		baseVariable = Variable_NewVector( tmpName, Variable_DataType_Double, componentsCount, &nDomainVerts, 
						   NULL, (void**)&arrayPtr, variable_Register,
						   varName[0], varName[1], varName[2], varName[3], varName[4],
						   varName[5], varName[6], varName[7], varName[8] );
	}
	Memory_Free( tmpName );

	tmpName = Stg_Object_AppendSuffix( numericField, "ErrorDofLayout" );

	errorDofLayout = DofLayout_New( tmpName, variable_Register, Mesh_GetDomainSize( constantMesh, MT_VERTEX ), constantMesh );
	if( componentsCount == 1 )
		DofLayout_AddAllFromVariableArray( errorDofLayout, 1, &baseVariable );
	else {
		for( var_I = 0; var_I < componentsCount; var_I++ ) {
			variable = Variable_Register_GetByName( variable_Register, varName[var_I] );
			variable->arrayPtrPtr = &baseVariable->arrayPtr;

			for( node_I = 0; node_I < Mesh_GetDomainSize( constantMesh, MT_VERTEX ); node_I++ )
				DofLayout_AddDof_ByVarName( errorDofLayout, varName[var_I], node_I );

			Memory_Free( varName[var_I] );
		}
		
		//errorDofLayout->baseVariables = Memory_Alloc_Array_Unnamed( Variable*, componentsCount );
		//errorDofLayout->baseVariables[0] = baseVariable;
	}
	
	Stg_Component_Build( errorDofLayout, NULL, False );
	Stg_Component_Initialise( errorDofLayout, NULL, False );

	Memory_Free( tmpName );

	tmpName = Stg_Object_AppendSuffix( numericField, "Error" );

	self->errorFieldList[field_I] = FeVariable_New( tmpName, constantMesh, constantMesh, errorDofLayout, NULL, NULL, NULL, 
							Mesh_GetDimSize( constantMesh ), False, "StgFEM_Native", "StgFEM_Native",
							"", "", False, False, context->fieldVariable_Register );
	/* so that the eqnation numbers don't get built for this guy */
	self->errorFieldList[field_I]->buildEqNums = False;

	tmpName = Stg_Object_AppendSuffix( self->errorFieldList[field_I], "Magnitude" );
	self->errorMagFieldList[field_I] = OperatorFeVariable_NewUnary( tmpName, self->errorFieldList[field_I], "Magnitude" );
	Memory_Free( tmpName );
	Stg_Component_Build( self->errorMagFieldList[field_I], context, False );
	self->errorMagFieldList[field_I]->_operator = Operator_NewFromName( self->errorMagFieldList[field_I]->operatorName, 
							self->errorFieldList[field_I]->fieldComponentCount, context->dim );
	self->errorMagFieldList[field_I]->fieldComponentCount = self->errorMagFieldList[field_I]->_operator->resultDofs;
	_OperatorFeVariable_SetFunctions( self->errorMagFieldList[field_I] );

	LiveComponentRegister_Add( context->CF->LCRegister, (Stg_Component*) baseVariable );
	LiveComponentRegister_Add( context->CF->LCRegister, (Stg_Component*) errorDofLayout );
	LiveComponentRegister_Add( context->CF->LCRegister, (Stg_Component*) self->errorFieldList[field_I] );
	LiveComponentRegister_Add( context->CF->LCRegister, (Stg_Component*) self->errorMagFieldList[field_I] );
}

void FieldTest_LoadReferenceSolutionFromFile( FeVariable* feVariable, Name referenceSolnName, Name referenceSolnPath, DomainContext* context ) {
	FeMesh*			feMesh			= feVariable->feMesh;
	char*			filename;
	unsigned		nx = 0, ny = 0, nz = 0, total;
	unsigned		lineNum = 0;
	double			resolution[3];
	double*			coord;
	Index			node_I, dim_I;
	unsigned		increments[3];
	double			value[3];
	unsigned		lineNum0;
	unsigned		dofAtEachNodeCount, dof_I;
	unsigned		meshSize 		= Mesh_GetLocalSize( feMesh, MT_VERTEX );
	unsigned		nDims			= Mesh_GetDimSize( feMesh );
	double			vertex0[3], coordPrime[3];
	double			Ni[27], values[27][3];
	unsigned		shapeFunc_I;
	double			*posx, *posy, *posz;
	double			**variables;
	unsigned		numShapeFuncs 		= ( nDims == 3 ) ? 27 : 9;
	double 			xi, eta, zeta;
	double 			a0, b0, c0;
	double 			a1, b1, c1;
	double 			a2, b2, c2;
	double 			m0, m1, m2, m3, m4, m5, m6;
#ifdef READ_HDF5
	hid_t			inputFile;
	hid_t 			dataSet, memSpace, dataSpace;
	hsize_t 		start[2], count[2], hSize;
#endif
	int 			sizes[3];
	double* 		data;
	int 			dataPos = 0;
	double			nodeDummy;

	Stg_Component_Initialise( feMesh,     context, False );
	Stg_Component_Initialise( feVariable, context, False );

	dofAtEachNodeCount = feVariable->fieldComponentCount;
													  /* .  h5  \0 */
	filename = Memory_Alloc_Array_Unnamed( char, strlen(referenceSolnPath) + strlen(referenceSolnName) + 1 + 2 + 1 );
	sprintf( filename, "%s%s.h5", referenceSolnPath, referenceSolnName );
#ifdef READ_HDF5
	inputFile = H5Fopen( filename, H5F_ACC_RDONLY, H5P_DEFAULT );
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
	dataSet = H5Dopen( inputFile, "/size" );
#else
	dataSet = H5Dopen( inputFile, "/size", H5P_DEFAULT );
#endif
	H5Dread( dataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sizes );
	nx = sizes[0];
	ny = sizes[1];
	total = nx * ny;
	if( nDims == 3 ) {
		nz = sizes[2];
		total *= nz;
	}
	H5Dclose( dataSet );

	posx = Memory_Alloc_Array_Unnamed( double, total );
	posy = Memory_Alloc_Array_Unnamed( double, total );
	if( nDims == 3 ) posz = Memory_Alloc_Array_Unnamed( double, total );
	variables = Memory_Alloc_2DArray_Unnamed( double, total, dofAtEachNodeCount );
	data = Memory_Alloc_Array_Unnamed( double, nDims + dofAtEachNodeCount );

	hSize = nDims + dofAtEachNodeCount;
	memSpace = H5Screate_simple( 1, &hSize, NULL );
	H5Sselect_all( memSpace );
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
	dataSet = H5Dopen( inputFile, "/data" );
#else
	dataSet = H5Dopen( inputFile, "/data", H5P_DEFAULT );
#endif
	dataSpace = H5Dget_space( dataSet );
	start[0] = 0;
	start[1] = 0;
	count[0] = 1;
	count[1] = nDims + dofAtEachNodeCount;
	H5Sselect_hyperslab( dataSpace, H5S_SELECT_SET, start, NULL, count, NULL );
	for( lineNum = 0; lineNum < total; lineNum++ ) {
		start[0] = lineNum;
		H5Sselect_hyperslab( dataSpace, H5S_SELECT_SET, start, NULL, count, NULL );
		H5Dread( dataSet, H5T_NATIVE_DOUBLE, memSpace, dataSpace, H5P_DEFAULT, data );

		nodeDummy = data[dataPos++];
		posx[lineNum] = data[dataPos++];
		posy[lineNum] = data[dataPos++];
		if( nDims == 3 ) posz[lineNum] = data[dataPos++];

		for( dof_I = 0; dof_I < dofAtEachNodeCount; dof_I++ )
			variables[lineNum][dof_I] = data[dataPos++];
	}

	H5Sclose( memSpace );
	H5Sclose( dataSpace );
	H5Dclose( dataSet );
#endif
	Memory_Free( data );

	resolution[0] = posx[1]  - posx[0];
	resolution[1] = posy[nx] - posy[0];
	if( nDims == 3 ) resolution[2] = posz[nx*ny] - posz[0];

	for( node_I = 0; node_I < meshSize; node_I++ ) {
		coord = Mesh_GetVertex( feMesh, node_I );
		increments[0] = (unsigned)( ( coord[0] - posx[0] ) / resolution[0] );
		increments[1] = (unsigned)( ( coord[1] - posy[0] ) / resolution[1] );
		if( nDims == 3 ) increments[2] = (unsigned)( ( coord[2] - posz[0] ) / resolution[2] );
			
		for( dim_I = 0; dim_I < nDims; dim_I++ ) 
			if( increments[dim_I] % 2 == 1 )
				increments[dim_I]--;
		if( increments[0] >= nx - 2 )
			increments[0] = nx - 3;
		if( increments[1] >= ny - 2 )
			increments[1] = ny - 3;
		if( nDims == 3 && increments[2] >= nz - 2 )
			increments[2] = nz - 3;

		lineNum0 = increments[0] + nx * increments[1];
		if( nDims == 3 ) lineNum0 += nx * ny * increments[2];
		if( lineNum0 >= total )
			Journal_Printf( context->info, "interpolation error: node value: %d resolution size: %d\n", lineNum0, total );
			
		vertex0[0] = posx[lineNum0];
		vertex0[1] = posy[lineNum0];
		if( nDims == 3 ) vertex0[2] = posz[lineNum0];
		
		/* for quadratic elements the resolution is twice the distance between the nodes */
		for( dim_I = 0; dim_I < nDims; dim_I++ )
			coordPrime[dim_I] = ( coord[dim_I] - vertex0[dim_I] ) / resolution[dim_I] - 1.0;

		/* assign the shape functions & interpolate quadratically */
		if( nDims == 2 ) {
			xi = coordPrime[0]; eta = coordPrime[1];
			a0 = xi - 1.0; b0 = eta - 1.0;
			a1 = 1.0 - xi * xi; b1 = 1.0 - eta * eta;
			a2 = xi + 1.0; b2 = eta + 1.0;
			m0 = 0.5 * xi; m1 = 0.5 * eta; m2 = 0.25 * xi * eta;

			Ni[0] = m2 * a0 * b0; Ni[1] = m1 * a1 * b0; Ni[2] = m2 * a2 * b0;
			Ni[3] = m0 * a0 * b1; Ni[4] =      a1 * b1; Ni[5] = m0 * a2 * b1;
			Ni[6] = m2 * a0 * b2; Ni[7] = m1 * a1 * b2; Ni[8] = m2 * a2 * b2;

			for( dof_I = 0; dof_I < dofAtEachNodeCount; dof_I++ ) {
				values[0][dof_I] = variables[lineNum0][dof_I];
				values[1][dof_I] = variables[lineNum0+1][dof_I];
				values[2][dof_I] = variables[lineNum0+2][dof_I];
				values[3][dof_I] = variables[lineNum0+nx][dof_I];
				values[4][dof_I] = variables[lineNum0+nx+1][dof_I];
				values[5][dof_I] = variables[lineNum0+nx+2][dof_I];
				values[6][dof_I] = variables[lineNum0+(2*nx)][dof_I];
				values[7][dof_I] = variables[lineNum0+(2*nx)+1][dof_I];
				values[8][dof_I] = variables[lineNum0+(2*nx)+2][dof_I];
					
				value[dof_I] = 0.0;
				for( shapeFunc_I = 0; shapeFunc_I < numShapeFuncs; shapeFunc_I++ )
					value[dof_I] += Ni[shapeFunc_I] * values[shapeFunc_I][dof_I];
			}
		}
		else {
			xi = coordPrime[0]; eta = coordPrime[1]; zeta = coordPrime[2];
			a0 = xi - 1.0; b0 = eta - 1.0; c0 = zeta - 1.0;
			a1 = 1.0 - xi * xi; b1 = 1.0 - eta * eta; c1 = 1.0 - zeta * zeta;
			a2 = xi + 1.0; b2 = eta + 1.0; c2 = zeta + 1.0;
			m0 = 0.5 * xi; m1 = 0.5 * eta; m2 = 0.5 * zeta;
			m3 = 0.25 * xi * eta; m4 = 0.25 * xi * zeta; m5 = 0.25 * eta * zeta;
			m6 = 0.125 * xi * eta * zeta;

			Ni[0]  = m6 * a0 * b0 * c0; Ni[1]  = m5 * a1 * b0 * c0; Ni[2]  = m6 * a2 * b0 * c0;
			Ni[3]  = m4 * a0 * b1 * c0; Ni[4]  = m2 * a1 * b1 * c0; Ni[5]  = m4 * a2 * b1 * c0;
			Ni[6]  = m6 * a0 * b2 * c0; Ni[7]  = m5 * a1 * b2 * c0; Ni[8]  = m6 * a2 * b2 * c0;

			Ni[9]  = m3 * a0 * b0 * c1; Ni[10] = m1 * a1 * b0 * c1; Ni[11] = m3 * a2 * b0 * c1;
			Ni[12] = m0 * a0 * b1 * c1; Ni[13] =      a1 * b1 * c1; Ni[14] = m0 * a2 * b1 * c1;
			Ni[15] = m3 * a0 * b2 * c1; Ni[16] = m1 * a1 * b2 * c1; Ni[17] = m3 * a2 * b2 * c1;

			Ni[18] = m6 * a0 * b0 * c2; Ni[19] = m5 * a1 * b0 * c2; Ni[20] = m6 * a2 * b0 * c2;
			Ni[21] = m4 * a0 * b1 * c2; Ni[22] = m2 * a1 * b1 * c2; Ni[23] = m4 * a2 * b1 * c2;
			Ni[24] = m6 * a0 * b2 * c2; Ni[25] = m5 * a1 * b2 * c2; Ni[26] = m6 * a2 * b2 * c2;
				
			for( dof_I = 0; dof_I < dofAtEachNodeCount; dof_I++ ) {
				values[0][dof_I]  = variables[lineNum0][dof_I];
				values[1][dof_I]  = variables[lineNum0+1][dof_I];
				values[2][dof_I]  = variables[lineNum0+2][dof_I];
				values[3][dof_I]  = variables[lineNum0+nx][dof_I];
				values[4][dof_I]  = variables[lineNum0+nx+1][dof_I];
				values[5][dof_I]  = variables[lineNum0+nx+2][dof_I];
				values[6][dof_I]  = variables[lineNum0+(2*nx)][dof_I];
				values[7][dof_I]  = variables[lineNum0+(2*nx)+1][dof_I];
				values[8][dof_I]  = variables[lineNum0+(2*nx)+2][dof_I];

				values[9][dof_I]  = variables[lineNum0+(nx*ny)][dof_I];
				values[10][dof_I] = variables[lineNum0+(nx*ny)+1][dof_I];
				values[11][dof_I] = variables[lineNum0+(nx*ny)+2][dof_I];
				values[12][dof_I] = variables[lineNum0+(nx*ny)+nx][dof_I];
				values[13][dof_I] = variables[lineNum0+(nx*ny)+nx+1][dof_I];
				values[14][dof_I] = variables[lineNum0+(nx*ny)+nx+2][dof_I];
				values[15][dof_I] = variables[lineNum0+(nx*ny)+(2*nx)][dof_I];
				values[16][dof_I] = variables[lineNum0+(nx*ny)+(2*nx)+1][dof_I];
				values[17][dof_I] = variables[lineNum0+(nx*ny)+(2*nx)+2][dof_I];

				values[18][dof_I] = variables[lineNum0+(2*nx*ny)][dof_I];
				values[19][dof_I] = variables[lineNum0+(2*nx*ny)+1][dof_I];
				values[20][dof_I] = variables[lineNum0+(2*nx*ny)+2][dof_I];
				values[21][dof_I] = variables[lineNum0+(2*nx*ny)+nx][dof_I];
				values[22][dof_I] = variables[lineNum0+(2*nx*ny)+nx+1][dof_I];
				values[23][dof_I] = variables[lineNum0+(2*nx*ny)+nx+2][dof_I];
				values[24][dof_I] = variables[lineNum0+(2*nx*ny)+(2*nx)][dof_I];
				values[25][dof_I] = variables[lineNum0+(2*nx*ny)+(2*nx)+1][dof_I];
				values[26][dof_I] = variables[lineNum0+(2*nx*ny)+(2*nx)+2][dof_I];

				value[dof_I] = 0.0;
				for( shapeFunc_I = 0; shapeFunc_I < numShapeFuncs; shapeFunc_I++ )
					value[dof_I] += Ni[shapeFunc_I] * values[shapeFunc_I][dof_I];
			}
		}

		FeVariable_SetValueAtNode( feVariable, node_I, value );
	}

	Memory_Free( filename );
	Memory_Free( posx );
	Memory_Free( posy );
	if( nDims == 3 ) Memory_Free( posz );
	Memory_Free( variables );

#ifdef READ_HDF5
	H5Fclose( inputFile );
#endif
}

/* by default, success of the physics test is set to false. this is reset if the test passes */
void FieldTest_EvaluatePhysicsTest( void* _context, void* data ) {
	DomainContext*			context			= (DomainContext*)_context;
	/* dodgy!!! - not sure how else to pass the self reference at an entry point */
	FieldTest*			self			= fieldTestSingleton;
	FieldTest_ExpectedResultFunc*	expectedFunc		= self->expectedFunc;
	FILE*				dumpExpectedFilePtr;
	char*				dumpExpectedFileName;
	int				dim_i, dof_i;

	if( expectedFunc( self->expectedData, context, self->expected, self->numeric, self->tolerance ) )
		self->expectedPass = True;

	if( strlen(self->dumpExpectedFileName) > 1 ) {
		dumpExpectedFileName = Memory_Alloc_Array_Unnamed( char, strlen(self->expectedFilePath) + 
									 strlen(self->dumpExpectedFileName) + 5 );
		sprintf( dumpExpectedFileName, "%s%s.out", self->expectedFilePath, self->dumpExpectedFileName );
		dumpExpectedFilePtr = fopen( dumpExpectedFileName, "a" );

		fprintf( dumpExpectedFilePtr, "%.8e ", self->numeric[context->timeStep].time );
		for( dim_i = 0; dim_i < context->dim; dim_i++ )
			fprintf( dumpExpectedFilePtr, "%.8e ", self->numeric[context->timeStep].place[dim_i] );
		for( dof_i = 0; dof_i < self->expectedDofs; dof_i++ )
			fprintf( dumpExpectedFilePtr, "%.8e ", self->numeric[context->timeStep].value[dof_i] );

		fprintf( dumpExpectedFilePtr, "\n" );

		if( context->timeStep == context->maxTimeSteps ) {
			if( self->expectedPass )
				fprintf( dumpExpectedFilePtr, "test result: PASS\n" );
			else
				fprintf( dumpExpectedFilePtr, "test result: FAIL\n" );
		}

		Memory_Free( dumpExpectedFileName );
		fclose( dumpExpectedFilePtr );
	}
}

void FieldTest_GenerateErrFields( void* _context, void* data ) {
	DomainContext*		context			= (DomainContext*)_context;
	/* this a really dodgy way to get the self ptr, as will only work if the textual name is consistent with that in 
	 * the XML - need to find a way to add an entry point which allows the self ptr to be passed as a void * */
	//FieldTest* 		self 			= LiveComponentRegister_Get( context->CF->LCRegister, "NumericFields" );
	/* this is also a dodgy way to get the self ptr, as its obtained from a global variable */
	FieldTest*		self			= fieldTestSingleton;
	FeVariable*		errorField;
	Index			lMeshSize, lElement_I;
	double			elErrorSq[3], elNormSq[3], elError[3];
	double			lAnalyticSq[3], gAnalyticSq[3];
	double			lErrorSq[3], gErrorSq[3];
	Bool			normalise		= self->normalise;
	Index			numDofs, dof_I;
	Index			field_I;
	Index			expected_I;
	Bool			pass;
	double			numericTestResult;
	double			eps			= self->epsilon;
	
	for( field_I = 0; field_I < self->fieldCount; field_I++ ) {
		/* should be using MT_VOLUME for the reference field mesh, but seems to have a bug */
		lMeshSize  = Mesh_GetLocalSize( self->constantMesh, MT_VERTEX );
		errorField = self->errorFieldList[field_I];
		numDofs	   = self->numericFieldList[field_I]->fieldComponentCount;

		assert( !strcmp( errorField->feMesh->name, "constantMesh" ) );

		for( dof_I = 0; dof_I < numDofs; dof_I++ ) {
			lAnalyticSq[dof_I] = 0.0;
			lErrorSq[dof_I]    = 0.0;
		}
	
		for( lElement_I = 0; lElement_I < lMeshSize; lElement_I++ ) {
			for( dof_I = 0; dof_I < numDofs; dof_I++ ) {
				elErrorSq[dof_I] = 0.0;
				elNormSq[dof_I]  = 0.0;
			}

			if( self->referenceSolnFromFile )
				FieldTest_ElementErrReferenceFromField( self, field_I, lElement_I, elErrorSq, elNormSq );
			else
				FieldTest_ElementErrAnalyticFromField( self, field_I, lElement_I, elErrorSq, elNormSq );

			for( dof_I = 0; dof_I < numDofs; dof_I++ ) {
				if( !self->normalise || elNormSq[dof_I] > eps ) {
					lAnalyticSq[dof_I] += elNormSq[dof_I];
					lErrorSq[dof_I]    += elErrorSq[dof_I];
					//elError[dof_I] = normalise ? sqrt( elErrorSq[dof_I] / ( elNormSq[dof_I] + eps ) ) : sqrt( elErrorSq[dof_I] );
					elError[dof_I] = normalise ? sqrt( elErrorSq[dof_I] / ( elNormSq[dof_I] ) ) : sqrt( elErrorSq[dof_I] );
				}
			}

			/* constant mesh, so node and element indices map 1:1 */
			FeVariable_SetValueAtNode( errorField, lElement_I, elError );
		}
	
		MPI_Allreduce( &lAnalyticSq, &gAnalyticSq, numDofs, MPI_DOUBLE, MPI_SUM, self->referenceFieldList[field_I]->communicator );
		MPI_Allreduce( &lErrorSq,    &gErrorSq,    numDofs, MPI_DOUBLE, MPI_SUM, self->referenceFieldList[field_I]->communicator );

		for( dof_I = 0; dof_I < numDofs; dof_I++ ) {
			self->gAnalyticSq[field_I][dof_I] = gAnalyticSq[dof_I];
			self->gErrorSq[field_I][dof_I]    = gErrorSq[dof_I];
			self->gErrorNorm[field_I][dof_I]  = sqrt( gErrorSq[dof_I] / gAnalyticSq[dof_I] );

			if( normalise )
				Journal_Printf( context->info, "%s - dof %d normalised global error: %.8e\n", 
				     	self->numericFieldList[field_I]->name, dof_I, self->gErrorNorm[field_I][dof_I] );
			else
				Journal_Printf( context->info, "%s - dof %d global error: %.8e\n",
					self->numericFieldList[field_I]->name, dof_I, sqrt( self->gErrorSq[field_I][dof_I] ) );
		}
	}
}

void FieldTest_ElementErrReferenceFromField( void* fieldTest, Index field_I, Index lElement_I, double* elErrorSq, double* elNormSq ) {
	FieldTest* 		self 			= (FieldTest*) fieldTest;
	FeVariable*		referenceField		= self->referenceFieldList[field_I];
	FeVariable*		numericField		= self->numericFieldList[field_I];
	FeMesh*			referenceMesh		= referenceField->feMesh;
	Index			constantElNode		= lElement_I;
	double*			coord			= Mesh_GetVertex( self->constantMesh, constantElNode );
	unsigned		nDims			= Mesh_GetDimSize( referenceMesh );
	Index			el_I;
	ElementType*		elType;
	Swarm*			intSwarm		= self->integrationSwarm;
	Index			cell_I;
	unsigned		cellParticleCount;
	Index			cParticle_I;
	IntegrationPoint*	cParticle;
	double			*xi, weight;
	double			globalCoord[3];
	double			detJac;
	double			reference[3], numeric[3];
	Index			numDofs			= numericField->fieldComponentCount;
	Index			dof_I;

	/* don't assume that the constant error field mesh & reference field mesh necessarily map 1:1 */
	Mesh_SearchElements( referenceMesh, coord, &el_I );
	elType = FeMesh_GetElementType( referenceMesh, el_I );

	cell_I = CellLayout_MapElementIdToCellId( intSwarm->cellLayout, el_I );
	cellParticleCount = intSwarm->cellParticleCountTbl[cell_I];

	for( cParticle_I = 0; cParticle_I < cellParticleCount; cParticle_I++ ) {
		cParticle = (IntegrationPoint*) Swarm_ParticleInCellAt( intSwarm, cell_I, cParticle_I );
		xi 	  = cParticle->xi;
		weight	  = cParticle->weight;

		FeMesh_CoordLocalToGlobal( referenceMesh, el_I, xi, globalCoord );
		FieldVariable_InterpolateValueAt( referenceField, globalCoord, reference );
		FieldVariable_InterpolateValueAt( numericField,   globalCoord, numeric   );

		detJac = ElementType_JacobianDeterminant( elType, referenceMesh, el_I, xi, nDims );

		for( dof_I = 0; dof_I < numDofs; dof_I++ ) {
			elErrorSq[dof_I] += ( numeric[dof_I] - reference[dof_I] ) * ( numeric[dof_I] - reference[dof_I] ) 
					    * weight * detJac;
			elNormSq[dof_I]  += reference[dof_I] * reference[dof_I] * weight * detJac;
		}
	}
}

void FieldTest_ElementErrAnalyticFromField( void* fieldTest, Index field_I, Index lElement_I, double* elErrorSq, double* elNormSq ) {
	FieldTest* 		self 			= (FieldTest*) fieldTest;
	FeVariable*		referenceField		= self->referenceFieldList[field_I];
	FeVariable*		numericField		= self->numericFieldList[field_I];
	FeMesh*			elementMesh		= self->elementMesh;
	Index			constantElNode		= lElement_I;
	double*			coord			= Mesh_GetVertex( self->constantMesh, constantElNode );
	unsigned		nDims			= Mesh_GetDimSize( elementMesh );
	Index			el_I;
	ElementType*		elType;
	Swarm*			intSwarm		= self->integrationSwarm;
	Index			cell_I;
	unsigned		cellParticleCount;
	Index			cParticle_I;
	IntegrationPoint*	cParticle;
	double			*xi, weight;
	double			globalCoord[3];
	double			detJac;
	double			analytic[3], numeric[3];
	Index			numDofs			= numericField->fieldComponentCount;
	Index			dof_I;
	/* corresponding analytic solution function for this feVariable, as assigned in the plugin */
	FieldTest_AnalyticSolutionFunc*	analyticSolution = self->_analyticSolutionList[self->analyticSolnForFeVarKey[field_I]];

	/* don't assume that the constant error field mesh & reference field mesh necessarily map 1:1 */
	Mesh_SearchElements( elementMesh, coord, &el_I );
	elType = FeMesh_GetElementType( elementMesh, el_I );

	cell_I = CellLayout_MapElementIdToCellId( intSwarm->cellLayout, el_I );
	cellParticleCount = intSwarm->cellParticleCountTbl[cell_I];

	for( cParticle_I = 0; cParticle_I < cellParticleCount; cParticle_I++ ) {
		cParticle = (IntegrationPoint*) Swarm_ParticleInCellAt( intSwarm, cell_I, cParticle_I );
		xi	  = cParticle->xi;
		weight    = cParticle->weight;

		FeMesh_CoordLocalToGlobal( elementMesh, el_I, xi, globalCoord );
		analyticSolution( self, globalCoord, analytic );
		FieldVariable_InterpolateValueAt( numericField, globalCoord, numeric );

		detJac = ElementType_JacobianDeterminant( elType, elementMesh, el_I, xi, nDims );

		for( dof_I = 0; dof_I < numDofs; dof_I++ ) {
			elErrorSq[dof_I] += ( numeric[dof_I] - analytic[dof_I] ) * ( numeric[dof_I] - analytic[dof_I] ) 
					    * weight * detJac;
			elNormSq[dof_I]  += analytic[dof_I] * analytic[dof_I] * weight * detJac;
		}
	}
}

void FieldTest_ElementErrAnalyticFromSwarm( void* fieldTest, Index var_I, Index lElement_I, double* elErrorSq, double* elNormSq ) {
	FieldTest* 		self 		= (FieldTest*) fieldTest;

	
}

void FieldTest_ElementErrReferenceFromSwarm( void* fieldTest, Index var_I, Index lElement_I, double* elErrorSq, double* elNormSq ) {
	FieldTest* 		self 		= (FieldTest*) fieldTest;


}

/* the first index 'func_I' denotes the index of the function in the analytic solution list to be applied to calculate
 * the analytic field at index 'field_I' in the analytic field list.
 *
 * the analytic fields are in the same order as their numeric counterparts are read in from the XML */
void FieldTest_AddAnalyticSolutionFuncToListAtIndex( void* fieldTest, Index func_I, FieldTest_AnalyticSolutionFunc* func, Index field_I ) {
	FieldTest* 	self 	= (FieldTest*) fieldTest;

	self->_analyticSolutionList[func_I] = func;
	self->analyticSolnForFeVarKey[field_I] = func_I;
}


