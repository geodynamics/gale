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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

const Type FieldTest_Type = "FieldTest";

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

	/* Assign singleton ptr */
	//mySingleton = self;
	return self;
}

void _FieldTest_Delete( void* fieldTest ) {
	FieldTest* self = (FieldTest*)fieldTest;
	
	if( self->numericField     ) Stg_Class_Delete( self->numericField );
	if( self->referenceField   ) Stg_Class_Delete( self->referenceField );
	if( self->numericSwarm     ) Stg_Class_Delete( self->numericSwarm );
	if( self->constantMesh     ) Stg_Class_Delete( self->constantMesh );
	if( self->integrationSwarm ) Stg_Class_Delete( self->integrationSwarm );

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
	FieldTest* 	self 			= (FieldTest*)fieldTest;
	Hook*		generateErrorFieldHook	= Hook_New( "Generate error field", FieldTest_GenerateErrorField, self->name );

	self->numericField	= Stg_ComponentFactory_ConstructByKey( cf, self->name, "NumericVariable",     FeVariable, False, data );
	self->numericSwarm	= Stg_ComponentFactory_ConstructByKey( cf, self->name, "NumericSwarm",        Swarm,      False, data );
	self->integrationSwarm 	= Stg_ComponentFactory_ConstructByKey( cf, self->name, "IntegrationSwarm",    Swarm,      True,  data );
	self->constantMesh	= Stg_ComponentFactory_ConstructByKey( cf, self->name, "ConstantElementMesh", FeMesh,     True,  data );
	self->referenceMesh	= Stg_ComponentFactory_ConstructByKey( cf, self->name, "ReferenceMesh",       FeMesh,     True,  data );

	self->normalise			= Stg_ComponentFactory_GetBool( cf, self->name, "normaliseByReferenceSoln", True );
	self->swarmVariableName 	= Stg_ComponentFactory_GetString( cf, self->name, "numericVariableName",       "" );
	self->referenceSolnFileName 	= Stg_ComponentFactory_GetString( cf, self->name, "referenceSolutionFileName", "" );
	self->referenceSolnPath		= Stg_ComponentFactory_GetString( cf, self->name, "referenceSolutionFilePath", "" );

	self->context		= Stg_ComponentFactory_ConstructByName( cf, "context", DomainContext, True, data );

	/* set up the entry point */
	_EntryPoint_AppendHook( Context_GetEntryPoint( self->context, AbstractContext_EP_Solve ), generateErrorFieldHook );
}

void _FieldTest_Build( void* fieldTest, void* data ) {
	FieldTest* 		self 		= (FieldTest*) fieldTest;

	/* create the analytic feVariable & its dofLayout */
	FieldTest_BuildReferenceField( fieldTest, data );

	FieldTest_BuildErrorField( fieldTest, data );	
}

void _FieldTest_Initialise( void* fieldTest, void* data ) {
	FieldTest* self = (FieldTest*) fieldTest;

	/* load the reference solution from file if req'd */
	if( self->referenceSolnFileName )
		FieldTest_LoadReferenceSolutionFromFile( self );
}

void _FieldTest_Execute( void* fieldTest, void* data ) {
	FieldTest* 		self 		= (FieldTest*) fieldTest;
}

void _FieldTest_Destroy( void* fieldTest, void* data ) {
	FieldTest* 		self 		= (FieldTest*) fieldTest;
}

void FieldTest_BuildReferenceField( void* fieldTest, void* data ) {
	FieldTest* 		self 			= (FieldTest*) fieldTest;
	FeMesh*			referenceMesh		= self->referenceMesh;
	FeVariable*		numericField		= self->numericField;
	Swarm*			numericSwarm		= self->numericSwarm;
	Variable_Register*	variable_Register	= self->context->variable_Register;
	Sync*			sync			= Mesh_GetSync( referenceMesh, MT_VERTEX );
	Name			tmpName;
	Dof_Index		componentsCount;
	SwarmVariable*		numericSwarmVariable	= NULL;
	unsigned		swarmVar_I;
	Name			varName[9];
	unsigned		var_I;
	unsigned		node_I;
	Variable*		variable;
	Variable*		baseVariable		= NULL;
	DofLayout*		referenceDofLayout	= NULL;
	FeVariable*		referenceField		= NULL;

	if( numericField )
		componentsCount = self->numericField->fieldComponentCount;

	if( numericSwarm ) {
		for( swarmVar_I = 0; swarmVar_I < numericSwarm->nSwarmVars; swarmVar_I++ )
			if( !strcmp( numericSwarm->swarmVars[swarmVar_I]->name, self->swarmVariableName ) )
				numericSwarmVariable = numericSwarm->swarmVars[swarmVar_I];

		if( numericSwarmVariable )
			componentsCount = numericSwarmVariable->dofCount;
	}

	Stg_Component_Build( referenceMesh, data, False );
	
	if( numericField )
		tmpName = Stg_Object_AppendSuffix( numericField, "-ReferenceVariable" );
	if( numericSwarm )
		tmpName = Stg_Object_AppendSuffix( numericSwarm, "-ReferenceVariable" );

	if( componentsCount == 1 ) {
		baseVariable = Variable_NewScalar( tmpName, 
						   Variable_DataType_Double, 
						   (unsigned*)&sync->nDomains,
						   NULL,
						   (void**)NULL,
						   variable_Register );
	}
	else {
		for( var_I = 0; var_I < componentsCount; var_I++ )
			Stg_asprintf( &varName[var_I], "%s-Component-%d", self->name, var_I );

		baseVariable = Variable_NewVector( tmpName, 
						   Variable_DataType_Double, 
						   componentsCount, 
						   (unsigned*)&sync->nDomains, 
						   NULL,
						   (void**)NULL, 
						   variable_Register,
						   varName[0],
						   varName[1],
						   varName[2],
						   varName[3],
						   varName[4],
						   varName[5],
						   varName[6],
						   varName[7],
						   varName[8] );
	}
	Memory_Free( tmpName );

	if( numericField )
		tmpName = Stg_Object_AppendSuffix( numericField, "-ReferenceDofLayout" );
	if( numericSwarm )
		tmpName = Stg_Object_AppendSuffix( numericSwarm, "-ReferenceDofLayout" );

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
	}
	Memory_Free( tmpName );

	if( numericField )
		tmpName = Stg_Object_AppendSuffix( numericField, "-ReferenceFieldVariable" );
	if( numericSwarm )
		tmpName = Stg_Object_AppendSuffix( numericSwarm, "-ReferenceFieldVariable" );

	referenceField = FeVariable_New( tmpName, referenceMesh, referenceMesh, referenceDofLayout, NULL, NULL, NULL, 
					 Mesh_GetDimSize( referenceMesh ), False, "StgFEM_Native", "StgFEM_Native",
					 "", "", False, False, self->context->fieldVariable_Register );

	LiveComponentRegister_Add( self->context->CF->LCRegister, (Stg_Component*) referenceDofLayout );
	LiveComponentRegister_Add( self->context->CF->LCRegister, (Stg_Component*) referenceField );

	self->referenceField = referenceField;
}

void FieldTest_BuildErrorField( void* fieldTest, void* data ) {
	FieldTest* 		self 			= (FieldTest*) fieldTest;
	FeMesh*			constantMesh		= self->constantMesh;
	FeVariable*		numericField		= self->numericField;
	Swarm*			numericSwarm		= self->numericSwarm;
	Variable_Register*	variable_Register	= self->context->variable_Register;
	Sync*			sync			= Mesh_GetSync( constantMesh, MT_VERTEX );
	Name			tmpName;
	Dof_Index		componentsCount;
	SwarmVariable*		numericSwarmVariable	= NULL;
	unsigned		swarmVar_I;
	Name			varName[9];
	unsigned		var_I;
	unsigned		node_I;
	Variable*		variable;
	Variable*		baseVariable		= NULL;
	DofLayout*		errorDofLayout		= NULL;
	FeVariable*		errorField		= NULL;

	if( numericField )
		componentsCount = self->numericField->fieldComponentCount;

	if( numericSwarm ) {
		for( swarmVar_I = 0; swarmVar_I < numericSwarm->nSwarmVars; swarmVar_I++ )
			if( !strcmp( numericSwarm->swarmVars[swarmVar_I]->name, self->swarmVariableName ) )
				numericSwarmVariable = numericSwarm->swarmVars[swarmVar_I];

		if( numericSwarmVariable )
			componentsCount = numericSwarmVariable->dofCount;
	}

	Stg_Component_Build( constantMesh, data, False );
	
	if( numericField )
		tmpName = Stg_Object_AppendSuffix( numericField, "-ErrorVariable" );
	if( numericSwarm )
		tmpName = Stg_Object_AppendSuffix( numericSwarm, "-ErrorVariable" );

	if( componentsCount == 1 ) {
		baseVariable = Variable_NewScalar( tmpName, 
						   Variable_DataType_Double, 
						   (unsigned*)&sync->nDomains,
						   NULL,
						   (void**)NULL,
						   variable_Register );
	}
	else {
		for( var_I = 0; var_I < componentsCount; var_I++ )
			Stg_asprintf( &varName[var_I], "%s-Component-%d", self->name, var_I );

		baseVariable = Variable_NewVector( tmpName, 
						   Variable_DataType_Double, 
						   componentsCount, 
						   (unsigned*)&sync->nDomains, 
						   NULL,
						   (void**)NULL, 
						   variable_Register,
						   varName[0],
						   varName[1],
						   varName[2],
						   varName[3],
						   varName[4],
						   varName[5],
						   varName[6],
						   varName[7],
						   varName[8] );
	}
	Memory_Free( tmpName );

	if( numericField )
		tmpName = Stg_Object_AppendSuffix( numericField, "-ErrorDofLayout" );
	if( numericSwarm )
		tmpName = Stg_Object_AppendSuffix( numericSwarm, "-ErrorDofLayout" );

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
	}
	Memory_Free( tmpName );

	if( numericField )
		tmpName = Stg_Object_AppendSuffix( numericField, "-ReferenceFieldVariable" );
	if( numericSwarm )
		tmpName = Stg_Object_AppendSuffix( numericSwarm, "-ReferenceFieldVariable" );

	errorField = FeVariable_New( tmpName, constantMesh, constantMesh, errorDofLayout, NULL, NULL, NULL, 
				     Mesh_GetDimSize( constantMesh ), False, "StgFEM_Native", "StgFEM_Native",
				     "", "", False, False, self->context->fieldVariable_Register );

	LiveComponentRegister_Add( self->context->CF->LCRegister, (Stg_Component*) errorDofLayout );
	LiveComponentRegister_Add( self->context->CF->LCRegister, (Stg_Component*) errorField );

	self->errorField = errorField;
}

void FieldTest_LoadReferenceSolutionFromFile( void* fieldTest ) {
	FieldTest* 		self 			= (FieldTest*) fieldTest;
	FeVariable*		feVariable		= self->referenceField;
	FeMesh*			feMesh			= self->referenceField->feMesh;
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
#ifdef HAVE_HDF5
	hid_t			inputFile;
	hid_t 			dataSet, memSpace, dataSpace;
	hsize_t 		start[2], count[2], hSize;
#endif
	int 			sizes[3];
	double* 		data;
	int 			dataPos = 0;

	Stg_Component_Initialise( feMesh,     self->context, False );
	Stg_Component_Initialise( feVariable, self->context, False );

	dofAtEachNodeCount = feVariable->fieldComponentCount;
													 	      /* .  h5  \0 */
	filename = Memory_Alloc_Array_Unnamed( char, strlen(self->referenceSolnPath) + strlen(self->referenceSolnFileName) + 1 + 2 + 1 );
	sprintf( filename, "%s%s.h5", self->referenceSolnPath, self->referenceSolnFileName );
#ifdef HAVE_HDF5
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
	if( nDims == 3 ) Memory_Alloc_Array_Unnamed( double, total );
	Memory_Alloc_2DArray_Unnamed( double, total, dofAtEachNodeCount );
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
			Journal_Printf( self->context->info, "interpolation error: node value: %d resolution size: %d\n", lineNum0, total );
			
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

#ifdef HAVE_HDF5
	H5Fclose( inputFile );
#endif
}

void FieldTest_GenerateErrorField( void* fieldTest, void* data ) {
	FieldTest* 		self 		= (FieldTest*) fieldTest;
	FeVariable*		referenceField	= self->referenceField;
	FeMesh*			referenceMesh	= referenceField->feMesh;
	unsigned		meshLocalSize	= Mesh_GetLocalSize( referenceMesh, MT_VOLUME );
	unsigned		lElement_I;
	FeVariable*		numericField	= self->numericField;
	Swarm*			numericSwarm	= self->numericSwarm;
	double			elErrorSq[3];
	double			elNormSq[3];
	double			elError[3];
	double			lAnalyticSq[3], gAnalyticSq[3];
	double			lErrorSq[3], gErrorSq[3];
	Bool			normalise	= self->normalise;
	Index			numDofs		= numericField->fieldComponentCount;
	Index			dof_I;
	
	assert( !strcmp( referenceMesh->name, "constantMesh" ) );

	for( dof_I = 0; dof_I < numDofs; dof_I++ ) {
		lAnalyticSq[dof_I] = 0.0;
		lErrorSq[dof_I]    = 0.0;
	}
	
	for( lElement_I = 0; lElement_I < meshLocalSize; lElement_I++ ) {
		for( dof_I = 0; dof_I < numDofs; dof_I++ ) {
			elErrorSq[dof_I] = 0.0;
			elNormSq[dof_I]  = 0.0;
		}

		if( numericField && self->_analyticSolution )
			FieldTest_ElementErrorAnalyticFromField( self, lElement_I, elErrorSq, elNormSq );

		else if( numericField && self->referenceSolnFileName )
			FieldTest_ElementErrorReferenceFromField( self, lElement_I, elErrorSq, elNormSq );

		else if( numericSwarm && self->_analyticSolution )
			FieldTest_ElementErrorAnalyticFromSwarm( self, lElement_I, elErrorSq, elNormSq );

		else if( numericSwarm && self->referenceSolnFileName )
			FieldTest_ElementErrorReferenceFromSwarm( self, lElement_I, elErrorSq, elNormSq );

		for( dof_I = 0; dof_I < numDofs; dof_I++ ) {
			lAnalyticSq[dof_I] += elNormSq[dof_I];
			lErrorSq[dof_I]    += elErrorSq[dof_I];
			elError[dof_I] = ( normalise ) ? sqrt( elErrorSq[dof_I] / elNormSq[dof_I] ) : 
							 sqrt( elErrorSq[dof_I] );
		}

		/* constant mesh, so node and element indices map 1:1 */
		FeVariable_SetValueAtNode( referenceMesh, lElement_I, elError );
	}
	
	MPI_Allreduce( &lAnalyticSq, &gAnalyticSq, 3, MPI_DOUBLE, MPI_SUM, referenceField->communicator );
	MPI_Allreduce( &lErrorSq,    &gErrorSq,    3, MPI_DOUBLE, MPI_SUM, referenceField->communicator );

	for( dof_I = 0; dof_I < numDofs; dof_I++ ) {
		self->gAnalyticSq[dof_I] = gAnalyticSq[dof_I];
		self->gErrorSq[dof_I]    = gErrorSq[dof_I];
		self->gErrorNorm[dof_I]  = sqrt( gErrorSq[dof_I] / gAnalyticSq[dof_I] );
	}
}

void FieldTest_ElementErrorReferenceFromField( void* fieldTest, Element_LocalIndex lElement_I, double* elErrorSq, double* elNormSq ) {
	FieldTest* 		self 			= (FieldTest*) fieldTest;
	FeVariable*		referenceField		= self->referenceField;
	FeVariable*		numericField		= self->numericField;
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
	double			reference[3], numeric[3]; /* scalar fields only */
	Index			numDofs			= numericField->fieldComponentCount;
	Index			dof_I;

	/* don't assume that the constant error field mesh & reference field mesh necessarily map 1:1 */
	Mesh_SearchElements( referenceField, coord, &el_I );
	elType = FeMesh_GetElementType( referenceMesh, el_I );

	cell_I = CellLayout_MapElementIdToCellId( intSwarm, el_I );
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

void FieldTest_ElementErrorAnalyticFromField( void* fieldTest, Element_LocalIndex lElement_I, double* elErrorSq, double* elNormSq ) {
	FieldTest* 		self 			= (FieldTest*) fieldTest;
	FeVariable*		referenceField		= self->referenceField;
	FeVariable*		numericField		= self->numericField;
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
	double			analytic[3], numeric[3]; /* scalar fields only */
	Index			numDofs			= numericField->fieldComponentCount;
	Index			dof_I;
	FieldTest_AnalyticSolutionFunc*		analyticSolution = self->_analyticSolution;

	/* don't assume that the constant error field mesh & reference field mesh necessarily map 1:1 */
	Mesh_SearchElements( referenceField, coord, &el_I );
	elType = FeMesh_GetElementType( referenceMesh, el_I );

	cell_I = CellLayout_MapElementIdToCellId( intSwarm, el_I );
	cellParticleCount = intSwarm->cellParticleCountTbl[cell_I];

	for( cParticle_I = 0; cParticle_I < cellParticleCount; cParticle_I++ ) {
		cParticle = (IntegrationPoint*) Swarm_ParticleInCellAt( intSwarm, cell_I, cParticle_I );
		xi	  = cParticle->xi;
		weight    = cParticle->weight;

		FeMesh_CoordLocalToGlobal( referenceMesh, el_I, xi, globalCoord );
		analyticSolution( self, globalCoord, analytic );
		FieldVariable_InterpolateValueAt( numericField, globalCoord, numeric );

		detJac = ElementType_JacobianDeterminant( elType, referenceMesh, el_I, xi, nDims );

		for( dof_I = 0; dof_I < numDofs; dof_I++ ) {
			elErrorSq[dof_I] += ( numeric[dof_I] - analytic[dof_I] ) * ( numeric[dof_I] - analytic[dof_I] ) 
					    * weight * detJac;
			elNormSq[dof_I]  += analytic[dof_I] * analytic[dof_I] * weight * detJac;
		}
	}
}

void FieldTest_ElementErrorAnalyticFromSwarm( void* fieldTest, Element_LocalIndex lElement_I, double* elErrorSq, double* elNormSq ) {
	FieldTest* 		self 		= (FieldTest*) fieldTest;

	
}

void FieldTest_ElementErrorReferenceFromSwarm( void* fieldTest, Element_LocalIndex lElement_I, double* elErrorSq, double* elNormSq ) {
	FieldTest* 		self 		= (FieldTest*) fieldTest;


}


